To demultiplex, run the demultiplex.py script (or submit is as a SLURM job using demultiplex.sh). After this, you can run the align_and_quantify.py script (or use the align_and_quantify.sh script) to create a .json for each unique UMI across each of the samples. In these .json files, the following information will be stored:

- Total number of reads mapped to UMI
- Total number of perfect hits to edited reference (hamming distance = 0, edited)
- Total number of perfect misses to unedited reference (hamming distance = 0, unedited)
- Total number partial hits (strong alignment score to edited reference)
- Total number partial misses (strong alignment score to unedited reference)

To visualize the editing effiencies across different treshold requirements for the number of reads mapped to a UMI, use the hist_UMI_editing_efficiencies.py script. In addition, you can visualize the distribution of number of reads mapped to UMIs using the number_in_each_UMI_hist.py script.
