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
#                   MAIN
# ----------------------------------------------------------------------------

path_base='../../../01_data/02_alignment/pannagram_v10_4/intermediate/'
path_version="${path_base}version_2024_12_18/"

mkdir -p ${path_version}

# ------------------------------------------------------------------------
# Copy alignments
echo "Copy alignments..."

path_aln="${path_version}aln/"
mkdir -p ${path_aln}

cp ${path_base}/consensus/extra2_* ${path_aln}

# ------------------------------------------------------------------------
# Copy gff annotations
echo "Copy annotations..."

path_gff="${path_version}gff/"
mkdir -p ${path_gff}

cp ${path_base}/annotation/gff_pan_merged_renamed.gff ${path_gff}
cp -r ${path_base}/annotation/own ${path_gff}

# ------------------------------------------------------------------------
# Copy gff features
echo "Copy gff features..."

path_features="${path_version}gff_features/"
mkdir -p ${path_features}

cp ${path_base}/annotation/features/* ${path_features}

# ------------------------------------------------------------------------
# Copy consensus sequences
echo "Copy consensus sequences..."

path_seq="${path_version}seq/"
mkdir -p ${path_seq}

cp ${path_base}/consensus/seq/*fasta ${path_seq}


# ------------------------------------------------------------------------
# Copy SVs
echo "Copy SVs..."

path_sv="${path_version}sv/"
mkdir -p ${path_sv}

cp -r ${path_base}/consensus/sv/gff ${path_sv}
cp ${path_base}/consensus/sv/seq_sv_big.fasta ${path_sv}
cp ${path_base}/consensus/sv/sv_pangen_beg.rds ${path_sv}
cp ${path_base}/consensus/sv/sv_pangen_end.rds ${path_sv}
cp ${path_base}/consensus/sv/sv_pangen_pos.rds ${path_sv}

# ------------------------------------------------------------------------
# Archive
echo "Archive..."
folder_name=$(basename "${path_version}")
tar -czf "${path_base}${folder_name}.tar.gz" -C "${path_version}" . > /dev/null 2>&1


echo "Done!"

