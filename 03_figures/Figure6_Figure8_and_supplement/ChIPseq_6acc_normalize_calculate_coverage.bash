#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=2:00:00
#SBATCH --array=1-66
#SBATCH --output=/groups/nordborg/pub/Mirjam/logs/array_job_slurm_%A_%a.out


export numberofline=$SLURM_ARRAY_TASK_ID

ml deeptools/3.1.2-foss-2018b-python-2.7.15
ml bedtools/2.27.1-foss-2018b


#list of samples
list=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/ChIP_samples_6acc_eracaps_names_only_antibodies.txt
# wc -l /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/ChIP_samples_6acc_eracaps_names_only_antibodies.txt
#66 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/ChIP_samples_6acc_eracaps_names_only_antibodies.txt

#extract the sample name

#get the sample name from the list
export sample=`sed -n $numberofline,"$numberofline"p $list | awk '{print $1}'`
export accession=`sed -n $numberofline,"$numberofline"p $list | awk '{print $2}'`
echo $sample $accession 
#determine which input sample to normalize against
export input=`awk -v a="$numberofline" '(NR==a) {split($1,b,"."); print b[1]"."b[2]".INPUT"}' $list`
echo $input
#set working directory

export working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus/Chipseq 

cd $working_folder/BAM_on_own_genome/$sample

bamCompare -b1 $sample.Aligned.sortedByCoord.out.bam -b2 $working_folder/BAM_on_own_genome/$input/$input.Aligned.sortedByCoord.out.bam   --outFileName  $sample.log2.input_norm.bedGraph --operation log2 --effectiveGenomeSize 119481543 -p 4 --ignoreDuplicates --outFileFormat bedgraph
#bamCompare -b1 $sample.Aligned.sortedByCoord.out.bam -b2 $working_folder/BAM_on_own_genome/$input/$input.Aligned.sortedByCoord.out.bam --outFileName  $sample.subtr.input_norm.bedGraph --operation subtract --effectiveGenomeSize 119481543 -p 4 --ignoreDuplicates --outFileFormat bedgraph


#bamCompare -b1 $sample.Aligned.sortedByCoord.out.bam -b2  $working_folder/BAM_on_own_genome/$input/$input.Aligned.sortedByCoord.out.bam --outFileName   $sample.subtr.input_norm.bw --operation subtract --effectiveGenomeSize 119481543 -p 4 --ignoreDuplicates
bamCompare -b1 $sample.Aligned.sortedByCoord.out.bam -b2  $working_folder/BAM_on_own_genome/$input/$input.Aligned.sortedByCoord.out.bam --outFileName   $sample.log2.input_norm.bw --operation log2 --effectiveGenomeSize 119481543 -p 4 --ignoreDuplicates


export mRNAs=/groups/nordborg/projects/cegs/alexandra/1001Gplus/annotations/mRNAs/mRNA.$accession.bed
export genes=/groups/nordborg/projects/cegs/alexandra/1001Gplus/annotations/genes/loci.$accession.bed
export SV=/groups/nordborg/projects/cegs/alexandra/1001Gplus/annotations/SVs/svs_v01.$accession.bed

# calculate chipseq coverage 

bedtools sort -i $mRNAs | bedtools map -a stdin -b $working_folder/BAM_on_own_genome/$sample/$sample.log2.input_norm.bedGraph -c 4 -o mean| awk -v OFS="\t"  '{print $1,$2,$3,$4,$5,$6,$7}'  > $working_folder/BAM_on_own_genome/coverage_on_own_genome/$sample.mRNAs.own.log2.mean_cov.bed

bedtools sort -i $genes | bedtools map -a stdin -b $working_folder/BAM_on_own_genome/$sample/$sample.log2.input_norm.bedGraph -c 4 -o mean | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7}' >$working_folder/BAM_on_own_genome/coverage_on_own_genome/$sample.loci.log2.mean_cov.bed

bedtools sort -i $SV | bedtools map -a stdin -b $working_folder/BAM_on_own_genome/$sample/$sample.log2.input_norm.bedGraph -c 4 -o mean| awk -v OFS="\t"  '{print $1,$2,$3,$4,$5,$6,$7}'  > $working_folder/BAM_on_own_genome/coverage_on_own_genome/$sample.SVs.log2.mean_cov.bed







