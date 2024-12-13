#!/bin/bash
#SBATCH --job-name=sim_genes
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

PATH_ANN="../../../01_data/02_alignment/pannagram_v10_4/intermediate/annotation/fasta/"

conda activate pannagram
parallel -j ${cores} simsearch -in_seq {} -on_seq ${PATH_ANN}genes.fasta -out ${PATH_ANN}out_{/.} ::: ${PATH_ANN}genes_*.fasta

