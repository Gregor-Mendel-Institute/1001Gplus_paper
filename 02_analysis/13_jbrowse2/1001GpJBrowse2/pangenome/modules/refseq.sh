#
# refseq.sh - Helper functions for reference sequence.  
# 
# Authors:
#   Joffrey Fitz (joffrey.fitz@tuebingen.mpg.de)
# 
# Copyright (c) 2024 Max Planck Institute for Biology, 
#   Tübingen, Germany, https://weigelworld.org
#  
# This file is part of 1001GpJBrowse2.
#   
# 1001GpJBrowse2 is free software: you can redistribute it 
# and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of 
# the License, or (at your option) any later version.
# 
# 1001GpJBrowse2 is distributed in the hope that it will be 
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with 1001GpJBrowse2.  
# If not, see <https://www.gnu.org/licenses/>.
# 


# process_refseq()
#
# Copy the pangenome sequence to $pangenome_seq_work_dir defined in
# config.sh, compress and index

function process_refseq() {
    # Copy to work dir
    echo -n "Copying pangenome sequence to ${pangenome_seq_work_dir} ... "
    cp ${pangenome_seq_dir}/${pangenome_seq_filename} ${pangenome_seq_work_dir}/
    echo "done."

    # Compress and index sequence
    echo -n "Compressing pangenome sequence ... "
    bgzip ${pangenome_seq_work_dir}/${pangenome_seq_filename}
    echo "done."

    echo -n "Indexing pangenome sequence ... "
    samtools faidx ${pangenome_seq_work_dir}/${pangenome_seq_filename}.gz
    echo "done."
}


# load_refseq()
#
# Load refseq to jbrowse with the "jbrowse add-track" command.
# Track data is a remote file defined in config.sh

function load_refseq() {
    # Load JBrowse
    url="${server_url}/${server_refseq_dir}/${pangenome_seq_filename}.gz"
    echo "Generating JBrowse2 config for $url ... "

    jbrowse add-assembly \
    $url \
    --displayName "Pangenome maximum consensus sequence" \
    --name pangenome \
    --out $jbrowse_dir

    echo "done."
}


# deploy_refseq_data()
#
# Rsync pangenome sequence and indices to server

function deploy_refseq() {
    # Rsync to public server
    dest="${user}@${hostname}:${server_fs_path}/${server_refseq_dir}/"
    echo "Deploying to $dest ..."
    rsync -ravz $pangenome_seq_work_dir/pangenome.fasta.gz* "$dest" 
}