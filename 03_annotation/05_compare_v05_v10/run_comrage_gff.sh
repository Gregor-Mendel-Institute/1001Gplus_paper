#!/bin/bash
#SBATCH --job-name=compare_gffs
#SBATCH --output=compare_gffs.log
#SBATCH --error=compare_gffs.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G

# ------------------------------------
# Modules
module load anaconda3/2019.03

conda activate pannagram
Rscript compare_gffs_duplicated.R

