import csv
import json
import matplotlib.pyplot as plt
import statistics
import numpy as np

files = [
    "2025-11-24_calculate_efficiency/custom_demux_Cod_DMSO_2sort_G1_AAGGAGTA_R1.fastq_umi_results.json",
    "2025-11-24_calculate_efficiency/custom_demux_Cod_NoDMSO_2sort_EdU_CTAAGCCT_R1.fastq_umi_results.json",
    "2025-11-24_calculate_efficiency/custom_demux_Temp_DMSO_2sort_G1_TAGATCGC_R4.fastq_umi_results.json",
    "2025-11-24_calculate_efficiency/custom_demux_Temp_NoDMSO_2sort_EdU_CTCTCTAT_R4.fastq_umi_results.json"
]

filestringers = [
    "Cod_DMSO_2sort_G1",
    "Cod_NoDMSO_2sort_EdU",
    "Temp_DMSO_2sort_G1",
    "Temp_NoDMSO_2sort_EdU"
]


counter = 0 

for file in files:

    print(file)

    counter2 = 0
    total_counter = 0

    totals = []

    with open(file, 'r') as infile:
        data = json.load(infile)

    for key, value in data.items():
        if key != "GGGGGGGGGGGG":
            if value["perfect_hits"] + value["partial_hits"] + value["perfect_misses"] + value["partial_misses"] > 0:
                totals.append(value["perfect_hits"] + value["partial_hits"] + value["perfect_misses"] + value["partial_misses"])
            else:
                counter2 += value["perfect_hits"] + value["partial_hits"] + value["perfect_misses"] + value["partial_misses"]
            total_counter += value["perfect_hits"] + value["partial_hits"] + value["perfect_misses"] + value["partial_misses"]
        else:
            pass

    # print(len(totals))

    # print(sum(totals))

    # print(statistics.mean(totals))

    binwidth = 500  # choose your width
    bins = np.arange(min(totals), max(totals) + binwidth, binwidth)

    plt.xlim(0, 15000)
    plt.yscale("log", base=10)
    plt.hist(totals, bins=bins, density=True, edgecolor="black")
    plt.xlabel("Total reads per UMI")
    plt.ylabel("Density")
    plt.title(f"Distribution of UMI totals – {filestringers[counter]}")

    plt.savefig(f"2025-11-24_calculate_efficiency/hist_UMI_counts/{filestringers[counter]}.png")
    plt.close()

    counter += 1


    print(1 - (counter2 / total_counter))