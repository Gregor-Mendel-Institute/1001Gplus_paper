# 1001Gplus_paper


This repository contains all the scripts used to produce the data for the paper.

## Structure of the repository

```
.
├── 01_general
├── 02_alignment
├── 03_annotation
├── 04_sv
├── 05_synteny
├── 06_snps
├── 07_saturation
├── 08_organellas
├── 09_expression
├── 10_epigenetics
├── 11_mod_genome
├── 12_graph
├── 13_jbrowse2
├── 14_repeats
├── 15_phylogeny
```

## Structure of the DATA folder

It is assumed that the initial data for the analysis is located in the data folder `../01_data/`, and the scripts, such as Pannagram, are in the folder `../03_tools/`.
The structure of the data folder is the following:
```
.
├── 01_assembly
│   ├── 00_raw_assembly
│   └── 01_fasta       	<- Genome .fasta files for all A.thaliana accessions
├── 02_alignment
│   └── pannagram_v10   <- Multiple genome alignment for A.thaliana
├── 03_graph
├── 04_annotation
│   ├── 01_raw_max
│   ├── 02_pannagram    <- Annotation groups annotation
│   ├── 03_edta
│   └── 04_new_genes
├── 05_expression
├── 06_methylation
├── 07_syri
├── 08_lyrata
│   ├── genomes
│   └── pannagram_v10  <- Genome alignments for A.thaliana on A.lyrata
└── 09_tair10
```



## Important comments
* Reconstructed 22001 is called 220011
* ..



## Create folders
If you want to create an analysis folder, you can use the script utils/create_directories.sh:

```
./utils/create_directories.sh 02_analysis/xx_yyy
```

The structure of the folder will be automatically created.