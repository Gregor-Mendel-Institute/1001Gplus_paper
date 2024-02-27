# JBrowse directory
jbrowse_dir=/tmp/1001Gp_pangenome

# Server URLs and paths
server_url=https://1001genomes.org/data/1001G+/releases/current/pangenome
user="root"
hostname="arnaud.tuebingen.mpg.de"
server_fs_path=/data/backup/public/pangenome
server_annotations_dir=annotations
server_refseq_dir=reference_sequence
server_rnaseq_dir=rnaseq
server_methylation_dir=methylation
server_jbrowse_dir=/var/www/html/jbrowse2/pangenome_test

# Assemblies table
assemblies_table=./assemblies_1001Gp.tsv

# pangenome assembly
pangenome_seq_filename=pangenome.fasta
pangenome_seq_work_dir=/tmp/pangenome_browser/reference_sequence

# Annotations
ann_work_dir=/tmp/pangenome_browser/gff
# Annotations are named "genes_v05_pangen_<ass_id>.gff.gz"
annotations_prefix=genes_v05_pangen_

# RNASeq
rnaseq_url_prefix=https://1001genomes.org/data/1001G+/releases/current/rnaseq/jbrowse
rnaseq_file_suffix=.FR.bw
rnaseq_work_dir=/tmp/pangenome_browser/rnaseq

# Methylation
methylation_url_prefix=https://1001genomes.org/data/1001G+/releases/current/methylation/
methylation_work_dir=/tmp/pangenome/methylation


# Define assemblies with methylation data
meth_ass=(
    "10002"
    "10015"
    "6024"
    "6909"
    "6966"
    "8236"
    "9075"
    "9537"
    "9543"
    "9728"
    "9888"
    "9905"
)
meth_acc=(
    "TueWal-2"
    "Sha"
    "Fly2-2"
    "Col-0"
    "Sq-1"
    "HSm"
    "Lerik1-4"
    "IP-Cum-1"
    "IP-Gra-0"
    "Stiav-1"
    "IP-Pva-1"
    "Ven-0"
)
