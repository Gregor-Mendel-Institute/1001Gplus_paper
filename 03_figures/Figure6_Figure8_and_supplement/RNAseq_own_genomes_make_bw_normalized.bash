#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=08:00:00
#SBATCH --array=1-125
#SBATCH --output=/groups/nordborg/user/aleksandra.kornienko/analyses/logs/array_job_slurm_%A_%a.out



export numberofline=$SLURM_ARRAY_TASK_ID 

export working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus
export eracapsaccessions=/groups/nordborg/projects/cegs/alexandra/1001Gplus/annotations/accessions.txt
export TAIR10=/groups/nordborg/projects/cegs/alexandra/GENOMES/1001.TAIR10.genome/TAIR10_all.fa
export genomefolder=/groups/nordborg/projects/cegs/alexandra/1001Gplus/genomes

export stargenome_folder=/groups/nordborg/user/aleksandra.kornienko/analyses/STAR/STAR2.7.1a/1001Gplus_genomes
export list_RNAseq_samples=/groups/nordborg/user/aleksandra.kornienko/analyses/ERA-CAPs/eracaps_all_samples_separately.txt
export samples_vs_fastq=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/eracaps_all_samples_separately_plus_corresponding_fastq_20011.txt

export fastq_folder=/groups/nordborg/projects/nordborg_common/1001g_plus/RNAseq/fastq
samples=/groups/nordborg/projects/cegs/alexandra/1001Gplus/BAM_on_own_genome/samples.txt 

#accession=`sed -n $numberofline,"$numberofline"p $samples_vs_fastq | awk '{print $1}'`
namesample=`sed -n $numberofline,"$numberofline"p $samples | awk '{print $1}'`

echo $accession
echo $namesample

#####
wc -l /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/eracaps_all_samples_separately_plus_corresponding_fastq_20011.txt
#125 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/eracaps_all_samples_separately_plus_corresponding_fastq_20011.txt


#check BAM quality 
bam=$working_folder/BAM_on_own_genome/$namesample/$namesample.Aligned.sortedByCoord.out.bam  
 

#make bigwig 
module load deeptools/3.1.2-foss-2018b-python-2.7.15
#bamCoverage -b $bam --normalizeUsing  BPM  --filterRNAstrand forward --effectiveGenomeSize 119481543 --ignoreForNormalization ChrX ChrM -o $working_folder/BAM_on_own_genome/$namesample/$namesample.normalized.F.bw --binSize=10
bamCoverage -b $bam --normalizeUsing  BPM  --filterRNAstrand reverse --effectiveGenomeSize 119481543 --ignoreForNormalization ChrX ChrM -o $working_folder/BAM_on_own_genome/$namesample/$namesample.normalized.R.bw --binSize=10








