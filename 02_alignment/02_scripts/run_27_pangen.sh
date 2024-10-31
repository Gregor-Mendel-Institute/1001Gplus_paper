#!/bin/bash
#SBATCH --job-name=run_27
#SBATCH --output=echo/output_%x_output.txt
#SBATCH --error=echo/error_%x_error.txt
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


# ------------------------------------
# Main
source activate pannagram

PATH_WORK="../../../"
PATH_PAN="${PATH_WORK}03_tools/pannagram/"
PATH_ALN="${PATH_WORK}01_data/02_alignment/"
PATH_ASSEMBLY="${PATH_WORK}01_data/01_assembly/01_fasta/"

${PATH_PAN}inst/pannagram.sh \
	-path_in  ${PATH_ASSEMBLY} \
    -path_out "${PATH_ALN}pannagram_v10/" \
    -nchr 5 \
    -cores 30 \
    -one2one \
    -purge_reps \
    -log 3 \
    -max_len_gap 25000 


