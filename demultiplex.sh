#!/bin/bash
#SBATCH --job-name=demux
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G

source env/bin/activate

python3 -u 2025-11-20_analyze_reads/demultiplex.py
