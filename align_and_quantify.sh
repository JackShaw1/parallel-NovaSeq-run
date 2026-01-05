#!/bin/bash
#SBATCH --job-name=results_NovaSeq_WYJS
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G

source env/bin/activate

cd 2025-11-24_calculate_efficiency

# conda activate bowtie2_env

cd ..

python3 2025-11-24_calculate_efficiency/align_and_quantify.py
