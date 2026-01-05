import pandas as pd
import numpy as np
import gzip
import os
import time
import math
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

from Bio import pairwise2
from Bio.pairwise2 import format_alignment

import json
import sys 

print('HELLO WORLD!', flush=True)

rows = ['A','B','C','D','E','F','G','H']
cols = [str(x) for x in np.arange(1,13,1)]

# for NovaSeq created by JS on 11.20.25
i5s = [
    "TAGATCGC", "CTCTCTAT", "TATCCTCT", "AGAGTAGA", "GTAAGGAG", "ACTGCATA", "AAGGAGTA",
    "CTAAGCCT", "CGTCTAAT", "TCTCTCCG", "TCGACTAG", "TTCTAGCT"
]

match_dic = {}
match_dic[('A', 'A')] = 1
match_dic[('A', 'C')] = 0
match_dic[('A', 'T')] = 0
match_dic[('A', 'G')] = 0
match_dic[('C', 'A')] = 0
match_dic[('C', 'C')] = 1
match_dic[('C', 'T')] = 0
match_dic[('C', 'G')] = 0
match_dic[('G', 'A')] = 0
match_dic[('G', 'C')] = 0
match_dic[('G', 'T')] = 0
match_dic[('G', 'G')] = 1
match_dic[('T', 'A')] = 0
match_dic[('T', 'C')] = 0
match_dic[('T', 'T')] = 1
match_dic[('T', 'G')] = 0
match_dic[('N', 'N')] = 1
match_dic[('N', 'A')] = 1
match_dic[('N', 'T')] = 1
match_dic[('N', 'G')] = 1
match_dic[('N', 'C')] = 1
match_dic[('A', 'N')] = 1
match_dic[('T', 'N')] = 1
match_dic[('G', 'N')] = 1
match_dic[('C', 'N')] = 1

def reverse_complement(seq):
    rev_comp = ""
    for c, char in enumerate(seq[::-1]):
        if char.upper() == "N":
            rev_comp += "N"
        elif char.upper() == 'A':
            rev_comp += "T"
        elif char.upper() == 'T':
            rev_comp += "A"
        elif char.upper() == 'G':
            rev_comp += "C"
        elif char.upper() == 'C':
            rev_comp += "G"
        else:
            rev_comp += "N"
    return rev_comp

print(reverse_complement('AAKAAAATTTTCCCCCGGGGN'))

def check_quality_and_convert(read):
    # read is [header, seq, plus, qual]
    phred_33_scoring = list('!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJ')
    updated_seq = ""
    
    for q, base_qual in enumerate(read[3]):
        if phred_33_scoring.index(base_qual) < 20:
            updated_seq += "N"
        else:
            updated_seq += read[1][q].upper()
       
    read[1] = updated_seq
    return read

def _count_reads_fastq_gz(fname):
    n_lines = 0
    with gzip.open(fname, 'rt') as f:
        for n_lines, _ in enumerate(f, start=1):
            pass
    return n_lines // 4

def fastq_iter(fname, label="", report_every=1_000_000):
    """Stream a FASTQ.gz file, yielding one [h, s, p, q] read at a time."""
    with gzip.open(fname, 'rt') as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline()
            p = f.readline()
            q = f.readline()
            if not q:
                break
            read = [h.rstrip(), s.rstrip(), p.rstrip(), q.rstrip()]
            yield read

def hamming_distance(seq1, seq2):
    hdist = 0
    for c, char in enumerate(seq1):
        if char.lower() != seq2[c].lower() or char.lower() == 'n' or seq2[c].lower() == 'n':
            hdist += 1
    return hdist

def find_run_fastqs(filename_string):
    R1_fname = R2_fname = R3_fname = R4_fname = None
    for file in os.listdir('.'):
        if filename_string in file and file.endswith(".fastq.gz"):
            if 'Read_1_' in file:   # R1
                R1_fname = file
            elif 'Read_2_' in file: # R2 (UMI / i7)
                R2_fname = file
            elif 'Read_3_' in file: # i5
                R3_fname = file
            elif 'Read_4_' in file: # R2
                R4_fname = file
    if not all([R1_fname, R2_fname, R3_fname, R4_fname]):
        raise RuntimeError("Could not find all four FASTQ.gz files for run " + filename_string)
    return R1_fname, R2_fname, R3_fname, R4_fname


BASES = ['A', 'C', 'G', 'T']

def neighbors_1mm(seq):
    seq = seq.upper()
    out = set()
    out.add(seq)
    for i, b in enumerate(seq):
        for nb in BASES:
            if nb != b:
                out.add(seq[:i] + nb + seq[i+1:])
    return out

def build_barcode_map(sample_lines):

    barcode_map = {}
    lib_barcodes = []

    for line in sample_lines:
        if not line:
            continue
        lib = line[0]
        i5_ref = line[1].strip().upper()
        lib_barcodes.append((lib, i5_ref))

        for nm in neighbors_1mm(i5_ref):
            if nm in barcode_map:
                continue
            barcode_map[nm] = (lib, i5_ref)
    return barcode_map, lib_barcodes

def demux_and_count(R1_fname, R2_fname, R3_fname, R4_fname, sample_lines, outdir):

    os.makedirs(outdir, exist_ok=True)

    barcode_map, lib_barcodes = build_barcode_map(sample_lines)
    per_barcode_counts = {i5: 0 for (_, i5) in lib_barcodes}
    barcode_counts = Counter()

    out_R1 = {}
    out_R4 = {}

    def get_handles(lib, i5_ref):
        key = (lib, i5_ref)
        if key not in out_R1:
            prefix = os.path.join(outdir, f"custom_demux_{lib}_{i5_ref}")
            out_R1[key] = open(prefix + "_R1.fastq", "w")
            out_R4[key] = open(prefix + "_R4.fastq", "w")
        return out_R1[key], out_R4[key]

    unmatched_R1 = open(os.path.join(outdir, "unmatched_read_1.fastq"), "w")
    unmatched_R4 = open(os.path.join(outdir, "unmatched_read_4.fastq"), "w")

    R1_iter = fastq_iter(R1_fname, "R1")
    R2_iter = fastq_iter(R2_fname, "R2/UMI")   # <-- add back
    R3_iter = fastq_iter(R3_fname, "R3/i5")
    R4_iter = fastq_iter(R4_fname, "R4")

    total_reads = 0
    total_assigned = 0
    t0 = time.time()

    # zip R1, R2, R3, R4
    for r, (r1, r2, r3, r4) in enumerate(zip(R1_iter, R2_iter, R3_iter, R4_iter), start=1):
        if r % 1_000_000 == 0:
            dt = time.time() - t0
            print(f"{r:,} reads processed ({dt:.1f} s)", flush=True)
            t0 = time.time()

        total_reads += 1

        # i5 handling as before
        i5_raw = r3[1].strip()
        i5_canonical = reverse_complement(i5_raw).upper()
        barcode_counts[i5_canonical] += 1

        hit = barcode_map.get(i5_canonical)
        if hit is None:
            for line in r1:
                unmatched_R1.write(line + "\n")
            for line in r4:
                unmatched_R4.write(line + "\n")
            continue

        lib, i5_ref = hit
        per_barcode_counts[i5_ref] += 1
        total_assigned += 1

        R1_out, R4_out = get_handles(lib, i5_ref)

        # UMI/i7 sequence from read 2
        umi_seq = r2[1].strip().upper()   # keep raw; no RC unless you know it's needed

        # add BOTH umi_seq and i5_canonical to headers
        r1_header = f"{r1[0]} {umi_seq} {i5_canonical}"
        r4_header = f"{r4[0]} {umi_seq} {i5_canonical}"

        for line in [r1_header, r1[1], r1[2], r1[3]]:
            R1_out.write(line + "\n")
        for line in [r4_header, r4[1], r4[2], r4[3]]:
            R4_out.write(line + "\n")

    unmatched_R1.close()
    unmatched_R4.close()
    for f in out_R1.values():
        f.close()
    for f in out_R4.values():
        f.close()

    return barcode_counts, per_barcode_counts, total_reads, total_assigned


# MAIN

if __name__ == "__main__":
    outdir = '2025-11-20_analyze_reads'
    os.makedirs(outdir, exist_ok=True)

    run_id = '22FL77LT1'
    R1_fname, R2_fname, R3_fname, R4_fname = find_run_fastqs(run_id)
    print("Using FASTQs:")
    print("  R1:", R1_fname)
    print("  R2 (UMI/i7):", R2_fname)
    print("  R3 (i5):", R3_fname)
    print("  R4:", R4_fname)

    with open('eWY043_sample_sheet.txt', 'r') as f:
        data = [x.rstrip().split() for x in f.readlines() if x.strip()]

    print("Sample sheet entries:")
    for line in data:
        print(line)

    print(f"Starting single-pass demux on {len(data)} samples...")
    start = time.time()
    barcode_counts, per_barcode_counts, total_reads, total_assigned = demux_and_count(
        R1_fname, R2_fname, R3_fname, R4_fname, data, outdir
    )
    elapsed = time.time() - start


    print("\nDemux summary:")
    print("  Total reads:      ", total_reads)
    print("  Total assigned:   ", total_assigned,
          f"({round(total_assigned / total_reads * 100, 1)}% assigned)")
    print("  Total demux time: ", round(elapsed, 1), "s")

    print("\nPer-barcode assigned counts:")
    for i5_ref, cnt in per_barcode_counts.items():
        print(i5_ref, cnt)

    # Barcode histogram
    plt.figure()
    plt.hist(list(barcode_counts.values()), bins=100)
    plt.yscale('log')
    plt.title("Raw reads (canonical i5 barcodes)")
    plt.xlabel('# Reads')
    plt.ylabel('Unique barcode count')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "barcode_hist.png"))

    # Valid barcodes
    valid_barcodes = []
    count = 0
    print("\nBarcodes with > 1000 reads:")
    for k, v in barcode_counts.items():
        if v > 1000:
            count += 1
            valid_barcodes.append(k)
            print(count, k, v)
    print("Total valid barcodes:", len(valid_barcodes))
