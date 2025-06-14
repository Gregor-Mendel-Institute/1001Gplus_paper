

samples=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/eracaps_125_RNAseq_samples.txt

export working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus

cd $working_folder

rm table 
touch table
echo "sample" "own_g_total_readN" "own_g_%_uniq" "own_g_%_nonuniq" "own_g_%_tooshort" "sample" "tair10_total_readN" "tair10_%_uniq" "tair10_%_nonuniq" "tair10_%_tooshort" > table
while read accession
do
cd $working_folder/BAM_on_own_genome/$accession
echo $accession
#number of unique and nonunique reads
nonunique=`cat $accession.Log.final.out| grep "% of reads mapped to multiple loci" | awk -v OFS="\t" '{print $9}' `
unique=`cat $accession.Log.final.out| grep "Uniquely mapped reads %" |grep -v non | awk -v OFS="\t" '{print $6}' `
total=`cat $accession.Log.final.out| grep "Number of input reads"| awk -v OFS="\t" '{print $6}' `	   
tooshort=`cat $accession.Log.final.out| grep "% of reads unmapped: too short"| awk -v OFS="\t" '{print $8}' `	

cd $working_folder/BAM_on_TAIR10/$accession
echo $accession
#number of unique and nonunique reads
nonunique1=`cat $accession.Log.final.out| grep "% of reads mapped to multiple loci" | awk -v OFS="\t" '{print $9}' `
unique1=`cat $accession.Log.final.out| grep "Uniquely mapped reads %" |grep -v non | awk -v OFS="\t" '{print $6}' `
total1=`cat $accession.Log.final.out| grep "Number of input reads"| awk -v OFS="\t" '{print $6}' `	   
tooshort1=`cat $accession.Log.final.out| grep "% of reads unmapped: too short"| awk -v OFS="\t" '{print $8}' `	
   
echo $accession $total $unique $nonunique $tooshort $accession $total1 $unique1 $nonunique1 $tooshort1
echo $accession $total $unique $nonunique $tooshort $accession $total1 $unique1 $nonunique1 $tooshort1 >> $working_folder/table
done < $samples


cd $working_folder

cat $working_folder/table > $working_folder/RNAseq_readN_on_own_genomes_and_TAIR10.bed

cp $working_folder/RNAseq_readN_on_own_genomes_and_TAIR10.bed /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/ 

