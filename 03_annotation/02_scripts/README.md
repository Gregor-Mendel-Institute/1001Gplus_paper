# Scripts for Annotation and Analysis

This repository contains scripts for performing preliminary annotation, similarity analysis, and feature extraction. Below are detailed instructions for running the scripts.

## Preliminary Annotation (define annotation groups)
To perform the preliminary annotation, use the following command:
```bash
sbatch run_27_annotation.sh
```

## Similarity Analysis
To perform similarity analysis for different datasets, run the respective scripts:

- **Gene Similarity Analysis:**
  ```bash
  sbatch run_sim_genes.sh
  ```

- **mRNA Similarity Analysis:**
  ```bash
  sbatch run_sim_mrna.sh
  ```

- **Transposable Elements (TE) Similarity Analysis:**
  ```bash
  sbatch run_sim_te.sh
  ```

## Feature Extraction
To extract features, execute the following R script:
```bash
Rscript run_get_features.R
```

## Rename accessions
Rename:  
220011 -> 22001_mod
10015 -> 6962
```bash
Rscript fix_renamed_acc.R
```

## Combine output in the version_XX archive
```bash
./combine_output.sh 
```
