import os
import gzip
import json
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
import tempfile
import re

def _parse_cigar(cigar: str):
    return [(op, int(length))
            for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]

from Bio.Align import PairwiseAligner

def align_to_four_bases(
    read_window: str,
    ref_seq: str,
    bowtie_index_prefix: str,
    motif: str = "CCAGGTG",
):
    """
    Align `read_window` to `ref_seq` using Biopython PairwiseAligner (local mode),
    then return the bases in the READ that align to `motif` in the REFERENCE.

    - Gaps in the *reference* are penalized extremely heavily so the motif
      should not be broken by gaps in the aligned reference.
    - If a reference position within the motif aligns to a gap in the read,
      we return '-' at that position.
    - If the motif is not found in `ref_seq`, returns None.

    NOTE: bowtie_index_prefix is ignored; it's only kept for API compatibility
    with the previous Bowtie2-based implementation.
    """
    read_window = str(read_window).upper()
    ref_seq = str(ref_seq).upper()
    motif = str(motif).upper()

    # Find motif in the ungapped reference sequence
    motif_start = ref_seq.find(motif)
    if motif_start == -1:
        # Motif not present in this reference
        return None

    # Local alignment: reference = target, read = query
    aligner = PairwiseAligner()
    aligner.mode = "local"

    # Allow gaps in the read with moderate penalty
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    aligner.query_open_gap_score = -5
    aligner.query_extend_gap_score = -1

    # Make gaps in the REFERENCE essentially forbidden
    aligner.target_open_gap_score = -1e6
    aligner.target_extend_gap_score = -1e6

    # Perform alignment
    best = aligner.align(ref_seq, read_window)[0]

    # Build mapping: ref_index -> read_index from aligned blocks
    ref_blocks, read_blocks = best.aligned
    ref_to_read = {}

    for (ref_start, ref_end), (read_start, read_end) in zip(ref_blocks, read_blocks):
        length = ref_end - ref_start
        for i in range(length):
            ref_pos = ref_start + i
            read_pos = read_start + i
            ref_to_read[ref_pos] = read_pos

    read_seq = str(best.query)

    # Collect bases in the read that align to each motif position in the reference
    bases_across = []
    for ref_pos in range(motif_start, motif_start + len(motif)):
        q_pos = ref_to_read.get(ref_pos)
        if q_pos is None or q_pos < 0 or q_pos >= len(read_seq):
            # This reference position is aligned to a gap or not covered
            bases_across.append("-")
        else:
            bases_across.append(read_seq[q_pos])

    return "".join(bases_across)


# =======================
# Alignment setup
# =======================

def hamming_distance(s1: str, s2: str) -> int:
    if len(s1) != len(s2):
        raise ValueError("Strings must be of equal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

aligner = PairwiseAligner()
aligner.mode = "local"

# Gap penalties set very high to discourage gaps
aligner.open_gap_score = -1e6
aligner.extend_gap_score = -1e6
aligner.query_open_gap_score = -1e6
aligner.query_extend_gap_score = -1e6
aligner.target_open_gap_score = -1e6
aligner.target_extend_gap_score = -1e6

coding_reference = "GGATGCTTCCAGCCAGATGAATCTCTTCGTGCTCACAGGTCAGGAGATCAGGCTGTCTTGCCACGGAGCCTTCGGCGGAAATACCTGATGATGAAGCTGGCCCAGGTGCTGCAGAAGGGATTCCTGCAGGTGTACAGCTA"
coding_reference_edited = "GGATGCTTCCAGCCAGATGAATCTCTTCGTGCTCACAGGTCAGGAGATCAGGCTGTCTTGCCACGGAGCCTTCGGCGGAAATACCTGATGATGAAGCTGGCATTCGTGCTGCAGAAGGGATTCCTGCAGGTGTACAGCTA"
template_reference = "ACGAGCAGAGCGCCGTGTGGCTGGATGCTTCCAGCCAGATGAATCTCTTCGTGCTCACAGGTCAGGAGATCAGGCTGTCTTGCCACGGAGCCTTCGGCGGAAATACCTGATGATGAAGCTGGCCCAGGTGCTGCAGAAGGGATTCCTGCAGGTGTACAGCTAGCT"
template_reference_edited = "ACGAGCAGAGCGCCGTGTGGCTGGATGCTTCCAGCCAGATGAATCTCTTCGTGCTCACAGGTCAGGAGATCAGGCTGTCTTGCCACGGAGCCTTCGGCGGAAATACCTGATGATGAAGCTGGCATTCGTGCTGCAGAAGGGATTCCTGCAGGTGTACAGCTAGCT"
quantification_window_unedited = "ATGATGAAGCTGGCCCAGGTGCTGCAGAAGGGATTCCTGC"
quantification_window_edited = "ATGATGAAGCTGGCATTCGTGCTGCAGAAGGGATTCCTGC"

FASTQ_DIR = "2025-11-20_analyze_reads"
RESULTS_DIR = "2025-11-24_calculate_efficiency"
os.makedirs(RESULTS_DIR, exist_ok=True)

# How many reads per chunk to send to worker processes
CHUNK_SIZE = 50000


# =======================
# Helpers
# =======================

def get_gapped_strings(best):
    lines = best.format().splitlines()

    target_candidates = [l for l in lines if l.startswith("target")]
    query_candidates  = [l for l in lines if l.startswith("query")]

    if target_candidates and query_candidates:
        target_line = target_candidates[0]
        query_line  = query_candidates[0]
        gapped_target = target_line.split()[2]
        gapped_query  = query_line.split()[2]
    else:
        nonempty = [l for l in lines if l.strip()]
        gapped_target = nonempty[0].replace(" ", "")
        gapped_query  = nonempty[2].replace(" ", "")

    return gapped_target, gapped_query


def alignment_identity(best):
    gT, gQ = get_gapped_strings(best)
    matches = 0
    aligned = 0
    for t, q in zip(gT, gQ):
        if t == "-" or q == "-":
            continue
        aligned += 1
        if t == q:
            matches += 1
    if aligned == 0:
        return 0.0
    return matches / aligned


def safe_load_json(path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path, "r") as f:
            return json.load(f)
    return {}


def atomic_save_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=2)
    os.replace(tmp, path)


def get_aligned_strings(aln):
    """
    Return the aligned portions of target (reference) and query (read)
    as two same-length strings (no gaps).
    """
    ref_chunks = []
    read_chunks = []

    # aln.aligned is a tuple: (coords_for_target, coords_for_query)
    ref_blocks, read_blocks = aln.aligned

    for (ref_start, ref_end), (read_start, read_end) in zip(ref_blocks, read_blocks):
        ref_chunks.append(str(aln.target[ref_start:ref_end]))
        read_chunks.append(str(aln.query[read_start:read_end]))

    ref_aln = "".join(ref_chunks)
    read_aln = "".join(read_chunks)
    return ref_aln, read_aln


# ====================================
# Worker processes chunk
# ====================================

def process_read_chunk(chunk, quant_start):
    """
    Process a chunk of reads for one FASTQ file.

    chunk: list of (seq_str, umi) tuples, seq_str already uppercased.
    quant_start: index in the read where the quantification window starts
                 (reference.find(quantification_window_unedited)).

    Returns:
        total, hits, partial_hits, misses, partial_misses, umi_results_chunk
    """
    quant_window_un = quantification_window_unedited.upper()
    quant_window_ed = quantification_window_edited.upper()
    window_len = len(quant_window_un)

    file_total = 0
    file_hits = 0
    file_partial_hits = 0
    file_misses = 0
    file_partial_misses = 0

    umi_results_chunk = {}

    for seq_str, umi in chunk:
        file_total += 1

        if umi not in umi_results_chunk:
            umi_results_chunk[umi] = {
                "total": 0,
                "perfect_hits": 0,
                "partial_hits": 0,
                "perfect_misses": 0,
                "partial_misses": 0,
            }

        umi_results_chunk[umi]["total"] += 1

        quant_window_read = seq_str[quant_start: quant_start + window_len]
        # if len(quant_window_read) != window_len:
        #     print(quant_window_read)
        #     continue

        hamming_un = hamming_distance(quant_window_read, quant_window_un)
        hamming_ed = hamming_distance(quant_window_read, quant_window_ed)

        # Perfect edited
        if hamming_ed == 0:
            file_hits += 1
            umi_results_chunk[umi]["perfect_hits"] += 1

        # Perfect unedited
        if hamming_un == 0:
            file_misses += 1
            umi_results_chunk[umi]["perfect_misses"] += 1

        # Partial edited vs unedited logic
        if (
            hamming_ed != 0 and
            hamming_un != 0 and
            hamming_un < 20 and
            hamming_ed < 20
        ):
            bowtie_index = "window_ref_index"
            bases_in_align = align_to_four_bases(
                seq_str,
                quant_window_un,
                bowtie_index_prefix=bowtie_index,
                motif="CCAGGTG",
            )

            if bases_in_align[:4] == "ATTC":
                file_partial_hits += 1
                umi_results_chunk[umi]["partial_hits"] += 1
            else:
                file_partial_misses += 1
                umi_results_chunk[umi]["partial_misses"] += 1

    return (
        file_total,
        file_hits,
        file_partial_hits,
        file_misses,
        file_partial_misses,
        umi_results_chunk,
    )


# ====================================
# Worker: process ONE FASTQ file
# ====================================

def process_fastq_file(filename, reference, kind, reference_edited,
                       executor, max_in_flight):
    """
    Process one FASTQ (coding or template) in chunks and write:
      - <RESULTS_DIR>/<safe_name>_files_results.json
      - <RESULTS_DIR>/<safe_name>_umi_results.json

    Parallelism happens over chunks of reads using the provided executor.
    """
    fastq_path = os.path.join(FASTQ_DIR, filename)
    print(f"[{kind}] Processing {filename}", flush=True)

    safe_name = filename.replace("/", "_").replace(":", "_")
    if safe_name.endswith(".gz"):
        safe_name = os.path.splitext(safe_name)[0]  # strip .gz

    results_path_files = os.path.join(RESULTS_DIR, f"{safe_name}_files_results.json")
    results_path_umi   = os.path.join(RESULTS_DIR, f"{safe_name}_umi_results.json")

    files_results = safe_load_json(results_path_files)
    umi_results   = safe_load_json(results_path_umi)

    # Keep same structure: { filename: {total, hits, ...} }
    if filename not in files_results:
        files_results[filename] = {
            "total": 0,
            "hits": 0,
            "partial_hits": 0,
            "misses": 0,
            "partial_misses": 0,
        }

    file_total = files_results[filename]["total"]
    file_hits  = files_results[filename]["hits"]
    file_partial_hits = files_results[filename]["partial_hits"]
    file_misses  = files_results[filename]["misses"]
    file_partial_misses = files_results[filename]["partial_misses"]

    # Precompute where the quantification window starts in this reference
    quant_start = reference.find(quantification_window_unedited)
    if quant_start == -1:
        raise ValueError(
            "Quantification window not found in reference for file {} (kind={})"
            .format(filename, kind)
        )

    opener = gzip.open if fastq_path.endswith(".gz") else open

    futures = []
    chunk = []

    with opener(fastq_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            # --- extract UMI from header ---
            header_fields = record.description.split()
            if len(header_fields) >= 3:
                umi = header_fields[2]  # "ATGTATCGAAAG"
            else:
                umi = "UNKNOWN"

            seq_str = str(record.seq).upper()
            chunk.append((seq_str, umi))

            if len(chunk) >= CHUNK_SIZE:
                # Submit chunk to worker processes
                futures.append(executor.submit(process_read_chunk, chunk, quant_start))
                chunk = []

                # Keep at most max_in_flight futures in flight
                if len(futures) >= max_in_flight:
                    done_future = next(as_completed(futures))
                    futures.remove(done_future)

                    (
                        c_total,
                        c_hits,
                        c_partial_hits,
                        c_misses,
                        c_partial_misses,
                        c_umi_results,
                    ) = done_future.result()

                    # Merge file-level counts
                    file_total += c_total
                    file_hits += c_hits
                    file_partial_hits += c_partial_hits
                    file_misses += c_misses
                    file_partial_misses += c_partial_misses

                    files_results[filename]["total"] = file_total
                    files_results[filename]["hits"] = file_hits
                    files_results[filename]["partial_hits"] = file_partial_hits
                    files_results[filename]["misses"] = file_misses
                    files_results[filename]["partial_misses"] = file_partial_misses

                    # Merge UMI-level counts
                    for umi_key, ucounts in c_umi_results.items():
                        if umi_key not in umi_results:
                            umi_results[umi_key] = {
                                "total": 0,
                                "perfect_hits": 0,
                                "partial_hits": 0,
                                "perfect_misses": 0,
                                "partial_misses": 0,
                            }
                        umi_results[umi_key]["total"] += ucounts["total"]
                        umi_results[umi_key]["perfect_hits"] += ucounts["perfect_hits"]
                        umi_results[umi_key]["partial_hits"] += ucounts["partial_hits"]
                        umi_results[umi_key]["perfect_misses"] += ucounts["perfect_misses"]
                        umi_results[umi_key]["partial_misses"] += ucounts["partial_misses"]

                    # incremental save after merging each completed chunk
                    atomic_save_json(results_path_files, files_results)
                    atomic_save_json(results_path_umi, umi_results)

    # Submit any remaining partial chunk
    if chunk:
        futures.append(executor.submit(process_read_chunk, chunk, quant_start))

    # Drain remaining futures
    for fut in as_completed(futures):
        (
            c_total,
            c_hits,
            c_partial_hits,
            c_misses,
            c_partial_misses,
            c_umi_results,
        ) = fut.result()

        file_total += c_total
        file_hits += c_hits
        file_partial_hits += c_partial_hits
        file_misses += c_misses
        file_partial_misses += c_partial_misses

        files_results[filename]["total"] = file_total
        files_results[filename]["hits"] = file_hits
        files_results[filename]["partial_hits"] = file_partial_hits
        files_results[filename]["misses"] = file_misses
        files_results[filename]["partial_misses"] = file_partial_misses

        for umi_key, ucounts in c_umi_results.items():
            if umi_key not in umi_results:
                umi_results[umi_key] = {
                    "total": 0,
                    "perfect_hits": 0,
                    "partial_hits": 0,
                    "perfect_misses": 0,
                    "partial_misses": 0,
                }
            umi_results[umi_key]["total"] += ucounts["total"]
            umi_results[umi_key]["perfect_hits"] += ucounts["perfect_hits"]
            umi_results[umi_key]["partial_hits"] += ucounts["partial_hits"]
            umi_results[umi_key]["perfect_misses"] += ucounts["perfect_misses"]
            umi_results[umi_key]["partial_misses"] += ucounts["partial_misses"]

        # incremental save after each remaining chunk
        atomic_save_json(results_path_files, files_results)
        atomic_save_json(results_path_umi, umi_results)

    # Final save to be extra-safe
    atomic_save_json(results_path_files, files_results)
    atomic_save_json(results_path_umi, umi_results)

    print(f"[{kind}] Done {filename}: hits={file_hits}, total={file_total}", flush=True)
    return filename, kind, file_total, file_hits



# ====================================
# Main: discover files & parallelize over chunks
# ====================================

if __name__ == "__main__":
    demuxed_files = os.listdir(FASTQ_DIR)

    coding_files = [
        f for f in demuxed_files
        if "demux_Cod" in f and "R1" in f and "2sort" in f
    ]

    template_files = [
        f for f in demuxed_files
        if "demux_Temp" in f and "R4" in f and "2sort" in f
    ]

    tasks = []
    for f in template_files:
        tasks.append(("template", f, template_reference, template_reference_edited))
    for f in coding_files:
        tasks.append(("coding", f, coding_reference, coding_reference_edited))

    if not tasks:
        print("No FASTQ files found to process.")
    else:
        max_workers = os.cpu_count() or 1
        max_in_flight = max_workers * 2  # how many chunks in flight at once

        print(f"Found {len(coding_files)} coding FASTQs and {len(template_files)} template FASTQs.")
        print(f"Running with up to {max_workers} worker processes.\n", flush=True)

        # One global process pool
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            for kind, fname, ref, ref_edited in tasks:
                fname_out, kind_out, total, hits = process_fastq_file(
                    fname, ref, kind, ref_edited, executor, max_in_flight
                )
                # summary line after each file finishes
                print(f"[summary] {kind_out} {fname_out}: hits={hits}, total={total}", flush=True)
