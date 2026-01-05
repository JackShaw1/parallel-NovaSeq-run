import os
import json
import matplotlib.pyplot as plt

files = [
    "2025-11-24_calculate_efficiency/custom_demux_Cod_DMSO_2sort_G1_AAGGAGTA_R1.fastq_umi_results.json",
    "2025-11-24_calculate_efficiency/custom_demux_Cod_NoDMSO_2sort_EdU_CTAAGCCT_R1.fastq_umi_results.json",
    "2025-11-24_calculate_efficiency/custom_demux_Temp_DMSO_2sort_G1_TAGATCGC_R4.fastq_umi_results.json",
    "2025-11-24_calculate_efficiency/custom_demux_Temp_NoDMSO_2sort_EdU_CTCTCTAT_R4.fastq_umi_results.json"
]


for path in files:
    total_counter = 0
    hit_counter = 0
    all_efficiencies = []
    stringer = path.split("/")[1].split("custom_demux_")[1].split(".")[0]
    print(path)
    with open(path, "r") as f:
        data = json.load(f)
    print(len(data))
    for umi, stats in data.items():
        total = int(stats.get("perfect_hits", 0)) + int(stats.get("partial_hits", 0)) + int(stats.get("perfect_misses", 0)) + int(stats.get("partial_misses", 0))
        # Minimum number of reads for UMI to be considered
        # CHANGE THIS 
        if total > 100:
            hits = int(stats.get("perfect_hits", 0)) + int(stats.get("partial_hits", 0))
            eff = hits / total
            all_efficiencies.append(eff)
            # Editing efficiency thresholds (for UMI 'hit'): .4 and .6
            if eff > .6:
                hit_counter += 1
                total_counter += 1
            elif eff < .4:
                total_counter += 1

    print(f"Number of UMI entries: {len(all_efficiencies)}")
    print(hit_counter / total_counter)
    print(total_counter)

    plt.figure(figsize=(7, 5))
    plt.hist(all_efficiencies, bins=20, range=(0, 1))
    plt.xlabel("Editing efficiency (hits / total)")
    plt.ylabel("Number of UMIs")
    plt.title(f"Distribution of UMI editing efficiencies ({stringer})")
    plt.tight_layout()
    plt.savefig(f"2025-11-24_calculate_efficiency/UMI_hists/{stringer}_umi_editing_efficiencies_histogram.png", dpi=300)

