#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=8:00:00
#SBATCH --array=1-14
#SBATCH --output=/groups/nordborg/pub/Mirjam/logs/array_job_slurm_%A_%a.out


export numberofline=$SLURM_ARRAY_TASK_ID

#module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/

ml star/2.7.1a-foss-2018b
ml cutadapt/1.18-foss-2018b-python-3.6.6
ml bedtools/2.27.1-foss-2018b
ml samtools/1.9-foss-2018b
 
 

export working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus/miRNA
export accessions_and_unmappedBam=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/RNA-seq_data/miRNA/14_miRNA_samples_flowers.txt
export accessions=/groups/nordborg/projects/cegs/alexandra/miRNAseq/accessions.miRNA.txt

export stargenome_folder=/groups/nordborg/user/aleksandra.kornienko/analyses/STAR/STAR2.7.1a/1001Gplus_genomes
export name_accession=`sed -n $numberofline,"$numberofline"p $accessions | awk '{print $1}'`

echo $name_accession

#export unmapped_bamfile=$unmappedbam_folder/$unmapped_BAM1
#echo $unmapped_bamfile 

#java -jar ${EBROOTPICARD}/picard.jar SamToFastq I=$unmapped_bamfile      FASTQ=$working_folder/fastq/$name_accession.fastq  


 #while read name_accession
 #do
 #echo $name_accession
 #cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.bed | awk '{print $1}' | uniq > $stargenome_folder/$name_accession.STARgenome_noSJDB/chromosomenames.txt
 #paste  $stargenome_folder/$name_accession.STARgenome_noSJDB/chrName.txt $stargenome_folder/$name_accession.STARgenome_noSJDB/chrLength.txt > $stargenome_folder/$name_accession.STARgenome_noSJDB/chrLength_2cols.txt
 #done< /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/RNA-seq_data/miRNA/accessions_miRNAseq.txt
 


#cutadapt -a AACTGTAGGCACCATCAAT --minimum-length 18 -o $fastq.trimmed.fastq $fastq.fastq

cd $working_folder


mkdir $working_folder/BAM_on_own_genome/$name_accession
export stargenomedir=$stargenome_folder/$name_accession.STARgenome_noSJDB
export out=$working_folder/BAM_on_own_genome/$name_accession/$name_accession.
fastq=/groups/nordborg/projects/cegs/alexandra/miRNAseq/fastq/$name_accession


STAR --runMode alignReads --readFilesCommand zcat \
--runThreadN 4 \
--runRNGseed 12345 \
--genomeDir $stargenomedir \
--readFilesIn $fastq.trimmed.fastq.gz  \
--outFileNamePrefix $out \
--limitBAMsortRAM 30000000000 \
--alignEndsType Extend5pOfRead1 \
--alignIntronMax 5000 \
--alignSJDBoverhangMin 1 \
--outReadsUnmapped Fastx \
--outSAMtype BAM Unsorted \
--outSAMmultNmax 100 \
--outSAMprimaryFlag AllBestScore \
--outSAMattributes NH HI AS nM NM MD jM jI XS \
--outFilterMultimapNmax 10 \
--outFilterMatchNmin 16 \
--outFilterMatchNminOverLread 0.66 \
--outFilterMismatchNmax 2 \
--outFilterMismatchNoverReadLmax 0.05 \
--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
--twopassMode None

  
samtools view $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Aligned.out.bam |  awk -v OFS="\t" '{print $3,$4,$5,$6,$10}' | awk -v OFS="\t" '($4=="21M"||$4=="22M"){print $1,$2,$2+$4,$5,$4}' |sortBed -i stdin  > $working_folder/BAM_on_own_genome/$name_accession/$name_accession.21-22nt.bed

samtools view $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Aligned.out.bam |  awk -v OFS="\t" '{print $3,$4,$5,$6,$10}' |  awk -v OFS="\t" '($4=="24M"){print $1,$2,$2+$4,$5,$4}' | sortBed -i stdin  > $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.bed


genomeCoverageBed -i $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.bed -g $stargenome_folder/$name_accession.STARgenome_noSJDB/chrLength_2cols.txt  -bg > $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.bedgraph

/groups/nordborg/projects/cegs/alexandra/software/bedGraphToBigWig $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.bedgraph $stargenome_folder/$name_accession.STARgenome_noSJDB/chrLength_2cols.txt $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.bw

genomeCoverageBed -i $working_folder/BAM_on_own_genome/$name_accession/$name_accession.21-22nt.bed -g $stargenome_folder/$name_accession.STARgenome_noSJDB/chrLength_2cols.txt  -bg > $working_folder/BAM_on_own_genome/$name_accession/$name_accession.21-22nt.bedgraph
/groups/nordborg/projects/cegs/alexandra/software/bedGraphToBigWig $working_folder/BAM_on_own_genome/$name_accession/$name_accession.21-22nt.bedgraph $stargenome_folder/$name_accession.STARgenome_noSJDB/chrLength_2cols.txt $working_folder/BAM_on_own_genome/$name_accession/$name_accession.21-22nt.bw


# read number for normalisation 
cd $working_folder/BAM_on_own_genome/
readN=`cat $name_accession/$name_accession.Log.final.out | grep "Uniquely mapped reads number"| awk '{print $6}'` 
readNtotal=`cat $name_accession/$name_accession.Log.final.out | grep "Number of input reads"| awk '{print $6}'` 
readNmulti=`cat $name_accession/$name_accession.Log.final.out | grep "Number of reads mapped to multiple loci"| awk '{print $9}'` 


genomeCoverageBed -i $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.bed -g $stargenome_folder/$name_accession.STARgenome_noSJDB/chrLength_2cols.txt  -d | awk -v OFS="\t" -v readN="$readN" '{print $1,$2-1,$2,$3*1000000/readN}' > $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.perbase.bed


genomeCoverageBed -i $working_folder/BAM_on_own_genome/$name_accession/$name_accession.21-22nt.bed -g $stargenome_folder/$name_accession.STARgenome_noSJDB/chrLength_2cols.txt  -d | awk -v OFS="\t" -v readN="$readN" '{print $1,$2-1,$2,$3*1000000/readN}' > $working_folder/BAM_on_own_genome/$name_accession/$name_accession.21-22nt.perbase.bed





export mRNAs=/groups/nordborg/projects/cegs/alexandra/1001Gplus/annotations/mRNAs/mRNA.$name_accession.bed
export genes=/groups/nordborg/projects/cegs/alexandra/1001Gplus/annotations/genes/loci.$name_accession.bed
export SV=/groups/nordborg/projects/cegs/alexandra/1001Gplus/annotations/SVs/svs_v01.$name_accession.bed
#export copies=/groups/nordborg/projects/cegs/alexandra/1001Gplus/annotations/SVs/svs_v01.$name_accession.bed


export working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus/miRNA

bedtools sort -i $genes | bedtools map -a stdin -b $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.perbase.bed -c 4 -o mean | sort  -k4,4 > $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.loci.coverage.perbase_calc.bed

bedtools sort -i $mRNAs | bedtools map -a stdin -b $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.perbase.bed -c 4 -o mean | sort  -k4,4 > $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.mRNAs.coverage.perbase_calc.bed


bedtools sort -i $SV | bedtools map -a stdin -b $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.perbase.bed -c 4 -o mean | sort  -k4,4 > $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.SVs.coverage.perbase_calc.bed





