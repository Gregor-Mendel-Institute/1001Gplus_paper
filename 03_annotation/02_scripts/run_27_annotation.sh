#!/bin/bash
#SBATCH --job-name=run_27_analys
#SBATCH --output=echo/analys_output.txt
#SBATCH --error=echo/analys_error.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=80G
#SBATCH --time=08:00:00

# ------------------------------------
# Modules
module load anaconda3/2019.03


# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

# Exit immediately if any command returns a non-zero status
set -e

# Keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# Define a trap for the EXIT signal
trap 'catch $?' EXIT

# Function to handle the exit signal
catch() {
    # Check if the exit code is non-zero
    if [ $1 -ne 0 ]; then
        echo "\"${last_command}\" command failed with exit code $1."
    fi
}


# ----------------------------------------------------------------------------
#            MAIN
# ----------------------------------------------------------------------------

# ------------------------------------
# VARIABLES
source activate pannagram
PATH_WORK="../../../"
PATH_ALN="${PATH_WORK}01_data/02_alignment/"
PATH_OUT="${PATH_ALN}pannagram_v10_4/"

PATH_ANNOTATION="${PATH_WORK}01_data/03_annotation/01_raw_max/"  #Annotation of genes

# ------------------------------------
# Get annotation groups
analys -path_msa ${PATH_OUT}/intermediate/annotation/ -annogroup ${PATH_ANNOTATION} -aln_type extra2_


# ------------------------------------
# Get fixed annotation
# Rscript run_gff_final.R








