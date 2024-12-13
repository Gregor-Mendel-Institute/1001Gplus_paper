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


PATH_ANN="../../../01_data/02_alignment/pannagram_v10_4/intermediate/annotation/fasta/"
PATH_SIM="../../../01_data/02_alignment/pannagram_v10_4/intermediate/annotation/simsearch/"

FILE_DB="/groups/nordborg/projects/the1001genomesplus/01_data/09_tair10/tair10_tes.fasta"


makeblastdb -in ${FILE_DB} -dbtype nucl

parallel -j ${cores} simsearch -in_seq {} -on_seq ${FILE_DB} -out ${PATH_SIM}out_{/.} ::: ${PATH_ANN}mrnas_*.fasta
