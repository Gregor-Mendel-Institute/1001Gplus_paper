## Modifying accession 22001
- Check the Snakemake pipeline for further details, everything is commented there
- Old workflow can be found in the workflow_modgenome.md 

## Comments
- Reason for change: "Better" (easier) alignment of the two genomes
- Yes, the chromosome arm switch is real
- Genome can be found here: data/final_genome/22001.scaffolds_corrected.v2.1.mod.fasta
- Split coordinates are here: data/final_bed/merged_mm.bed
- To run the pipeline, install [snakemake](https://snakemake.readthedocs.io/en/stable/), then run the pipeline: 
    - "cd scripts"
    - "snakemake"
- Pipeline based on minimap2, fpa and seqkit