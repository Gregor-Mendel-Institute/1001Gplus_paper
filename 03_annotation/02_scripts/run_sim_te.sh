#!/bin/bash
#SBATCH --job-name=sim_te
#SBATCH --output=echo/output_%x_output.txt
#SBATCH --error=echo/error_%x_error.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=100G
#SBATCH --time=08:00:00

# ------------------------------------
# Modules
module load anaconda3/2019.03

cores=30

conda activate pannagram

parallel -j ${cores} simsearch -in_seq {} -on_seq /groups/nordborg/projects/the1001genomesplus/01_data/09_tair10/tair10_tes.fasta -out out_te_{/.} ::: mrnas_*.fasta
