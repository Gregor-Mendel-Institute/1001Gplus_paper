
bedtools/2.27.1-foss-2018b
cd /groups/nordborg/projects/cegs/alexandra/1001Gplus/annotations
mkdir SVs
mkdir genes
mkdir SAFs
mkdir mRNAs


# annotations from Anna 
/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/version_2024_12_18.tar.gz
tar -xvzf version_2024_12_18.tar.gz
tar -xvzf svs_2024_12_23.tar.gz

#annotations are here 

export gfffolder=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/gff/own
export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations
export accessions=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/accessions.txt
mkdir $annotationfolder/SVs
mkdir $annotationfolder/genes
mkdir $annotationfolder/SAFs
mkdir $annotationfolder/mRNAs
mkdir $annotationfolder/copies


# make bed files with loci (gene - from pangenome consenseus locus onto every genome) annotations 
cd $gfffolder
while read accession
do
echo $accession
gff=gff_
gff+=$accession.gff
cat  $gff | grep gene | awk -v OFS="\t" '{split($9,a,";");print $1, $4, $5, a[1],0,$7}'| sed 's/ID=//g' | grep -v "e+" > $annotationfolder/genes/loci.$accession.bed
#make saf 
cat  $gff | grep gene | awk -v OFS="\t" '{split($9,a,";");print $1, $4, $5, a[1],0,$7}'| sed 's/ID=//g'| awk -v OFS="\t" '{split($4,a,".");print a[1], $1, $2, $3, $6}'> $annotationfolder/genes/loci.$accession.saf
done < $accessions


# make bed files with mRNA loci annotations 
cd $gfffolder
while read accession
do
echo $accession
gff=gff_
gff+=$accession.gff
#make bed6 of mRNA annotation
cat  $gff | grep mRNA | awk -v OFS="\t" '{split($9,a,";");print $1, $4, $5, a[1],0,$7}'| sed 's/ID=//g' | awk -v OFS="\t" '{split($4,a,".");print $1, $2, $3,a[1], $5, $6}'| grep -v "e+" >  $annotationfolder/mRNAs/mRNA.$accession.bed
#make saf files for featureCounts
cat  $gff | grep exon | awk -v OFS="\t" '{split($9,a,";");print $1, $4, $5, a[1],0,$7}'| sed 's/ID=//g'| awk -v OFS="\t" '{split($4,a,".");print a[1], $1, $2, $3, $6}'> $annotationfolder/mRNAs/mRNA_exons.$accession.saf
done < $accessions


# make annotations of SVs
cd /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/sv/gff

while read accession
do
echo $accession

SV=svs_acc_
SV+=$accession
SV+=_v06
#convert to bed (fix negative coordinate problem)
cat $SV.gff | awk -v OFS="\t" '{split($9,a,".");print $1,$4,$5,a[1],$3,$9}' | awk -v OFS="\t" '{split($4,a,"=");print $1,$2,$3,a[2],$5,$6}'| awk -v OFS="\t" '{sub("-","",$3);sub("-","",$2);print $1,$2,$3,$4,"type="$5";"$6,"."  }'| awk -v OFS="\t" '($3>0 && $2>0){if ($3>$2) print $1,$2,$3,$4,$5,$6;if ($3<$2) print $1,$3,$2,$4,$5,$6;  }'| grep -v "e+" >  $annotationfolder/SVs/svs_v06.$accession.bed
#convert to saf 
cat  $annotationfolder/SVs/svs_v06.$accession.bed |awk -v OFS="\t" '{print $4, $1, $2, $3,"."}' > $annotationfolder/SVs/svs_v06.$accession.saf
done < $accessions

# make info table about SVs 
#cat svs_v03.gff | grep -v Pan| awk -v OFS="\t" '{split($9,a,".");print $1,$4,$5,$3,a[1],a[2]}' | awk -v OFS="\t" '{split($5,a,"=");split($6,b,"=");print $1,$2,$3,$4,a[2],b[2]}'> sv_info.bed 

#fix 220011 
accession=220011
cat $annotationfolder/SVs/svs_v06.$accession.bed |sed 's/22001_mod/220011/g'  > tmp 
cat tmp > $annotationfolder/SVs/svs_v06.$accession.bed

cat $annotationfolder/SVs/svs_v06.$accession.saf |sed 's/22001_mod/220011/g'  > tmp 
cat tmp > $annotationfolder/SVs/svs_v06.$accession.saf

#convert to saf 
cat  $annotationfolder/SVs/svs_v06.$accession.bed |awk -v OFS="\t" '{print $4, $1, $2, $3,"."}' > $annotationfolder/SVs/svs_v06.$accession.saf
done < $accessions


# combine genes (loci) and SVs for TPM calculation 


cd /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/sv/gff

while read accession
do
echo $accession
cat $annotationfolder/genes/loci.$accession.saf $annotationfolder/SVs/svs_v06.$accession.saf > $annotationfolder/loci_and_SVs.$accession.saf
done < $accessions


#check how annotations look like 
cd $annotationfolder/SVs


head svs_v06.22006.saf
SVgr_1_id_000002        22006_Chr1      7837    7838    .
SVgr_1_id_000003        22006_Chr1      7847    7856    .
SVgr_1_id_000005        22006_Chr1      8452    8475    .
SVgr_1_id_000007        22006_Chr1      8719    8724    .
SVgr_1_id_000008        22006_Chr1      8909    8941    .
SVgr_1_id_000012        22006_Chr1      9557    9558    .
SVgr_1_id_000019        22006_Chr1      10865   10887   .
SVgr_1_id_000020        22006_Chr1      10895   10910   .
SVgr_1_id_000022        22006_Chr1      10992   11027   .
SVgr_1_id_000024        22006_Chr1      11096   11128   .


head svs_v06.22006.bed  
22006_Chr1      7837    7838    SVgr_1_id_000002        type=deletion;ID=SVgr_1_id_000002.22006;len_init=2;len_acc=2       .
22006_Chr1      7847    7856    SVgr_1_id_000003        type=deletion;ID=SVgr_1_id_000003.22006;len_init=10;len_acc=10     .
22006_Chr1      8452    8475    SVgr_1_id_000005        type=deletion;ID=SVgr_1_id_000005.22006;len_init=24;len_acc=24     .
22006_Chr1      8719    8724    SVgr_1_id_000007        type=deletion;ID=SVgr_1_id_000007.22006;len_init=6;len_acc=6       .
22006_Chr1      8909    8941    SVgr_1_id_000008        type=indel;ID=SVgr_1_id_000008.22006;len_init=33;len_acc=33.

cd $annotationfolder/genes

head loci.6244.saf
AT1SG10000001   6244_Chr1       3723    5579    +
AT1SG20000001   6244_Chr1       6180    8612    -
AT1SG20000002   6244_Chr1       11874   12950   -
AT1SG10000003   6244_Chr1       23563   31121   +
AT1SG20000003   6244_Chr1       31424   32712   -
AT1SG20000004   6244_Chr1       34034   37677   -
AT1SG20000005   6244_Chr1       38922   40901   -
AT1G01073       6244_Chr1       44701   44811   +
AT1SG20000006   6244_Chr1       45292   46578   -
AT1SG20000007   6244_Chr1       47494   48955   -
head loci.6244.bed
6244_Chr1       3723    5579    AT1SG10000001   0       +
6244_Chr1       6180    8612    AT1SG20000001   0       -
6244_Chr1       11874   12950   AT1SG20000002   0       -
6244_Chr1       23563   31121   AT1SG10000003   0       +
6244_Chr1       31424   32712   AT1SG20000003   0       -
6244_Chr1       34034   37677   AT1SG20000004   0       -
6244_Chr1       38922   40901   AT1SG20000005   0       -
6244_Chr1       44701   44811   AT1G01073       0       +
6244_Chr1       45292   46578   AT1SG20000006   0       -



cd $annotationfolder/mRNAs
head mRNA.22006.bed
22006_Chr1      5830    7793    AT1SG10000001   0       +
22006_Chr1      8592    10862   AT1SG20000001   0       -
22006_Chr1      14082   15158   AT1SG20000002   0       -
22006_Chr1      23987   31546   AT1SG10000003   0       +
22006_Chr1      31849   33137   AT1SG20000003   0       -
22006_Chr1      34467   37536   AT1SG20000004   0       -
22006_Chr1      39368   41341   AT1SG20000005   0       -
22006_Chr1      45971   47257   AT1SG20000006   0       -
22006_Chr1      48172   49637   AT1SG20000007   0       -
22006_Chr1      50772   51451   AT1SG20000008   0       -

head mRNA_exons.22006.saf
AT1SG10000001   22006_Chr1      5830    5983    +
AT1SG10000001   22006_Chr1      6066    6346    +
AT1SG10000001   22006_Chr1      6649    6768    +
AT1SG10000001   22006_Chr1      6869    7258    +
AT1SG10000001   22006_Chr1      7337    7489    +



## make annotations of copies 
 ml bedtools/2.27.1-foss-2018b
 wc -l  /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/version_X_eracaps/genes_remove_nestgr.txt
#2001 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/version_X_eracaps/genes_remove_nestgr.txt

export badgenes=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/version_X_eracaps/genes_remove_nestgr.txt
cd /groups/nordborg/projects/the1001genomesplus/01_data/02_alignment/pannagram_v10_4/intermediate/annotation/cnv
export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations
export accessions=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/accessions.txt


while read accession
do 
echo $accession
index1=$accession
index1+=_chr1_85_85.gff
index2=$accession
index2+=_chr2_85_85.gff
index3=$accession
index3+=_chr3_85_85.gff
index4=$accession
index4+=_chr4_85_85.gff
index5=$accession
index5+=_chr5_85_85.gff
cat */simsearch.$index1 */simsearch.$index2 */simsearch.$index3 */simsearch.$index4 */simsearch.$index5   | grep -v "e+" | awk  -v OFS="\t" -v acc="$accession" '($9 ~ acc){print $0}' |awk  -v OFS="\t" '{split($9,a,"|");split(a[2],b,".");print $1, $4-1, $5, b[1], ".", $7}'  | sortBed -i stdin >   $annotationfolder/copies/$accession.copies.bed 
done < $accessions

accession=22001_mod
index1=$accession
index1+=_chr1_85_85.gff
index2=$accession
index2+=_chr2_85_85.gff
index3=$accession
index3+=_chr3_85_85.gff
index4=$accession
index4+=_chr4_85_85.gff
index5=$accession
index5+=_chr5_85_85.gff
cat */simsearch.$index1 */simsearch.$index2 */simsearch.$index3 */simsearch.$index4 */simsearch.$index5   | grep -v "e+" | awk  -v OFS="\t" -v acc="220011" '($9 ~ acc){print $0}' |awk  -v OFS="\t" '{split($9,a,"|");split(a[2],b,".");print $1, $4-1, $5, b[1], ".", $7}'  | sortBed -i stdin >   $annotationfolder/copies/220011.copies.bed 





##################################################
############## find tandem copies 
#########################################################
ml bedtools/2.27.1-foss-2018b

mkdir  $annotationfolder/copies/tandem 
while read accession
do
echo $accession
closestBed -d -io -a $annotationfolder/copies/$accession.copies.bed  -b $annotationfolder/copies/$accession.copies.bed   | awk -v OFS="\t" '($13<10000){split ($4,a,".");split ($10,b,"."); print $0, a[1],b[1]}' | awk -v OFS="\t" '($14==$15){print $14}' | sort -u > $annotationfolder/copies/tandem/genes_with_tandem.DUP.$accession.txt 

closestBed -d -io -a $annotationfolder/copies/$accession.copies.bed  -b $annotationfolder/copies/$accession.copies.bed   | awk -v OFS="\t" '($13<10000){split ($4,a,".");split ($10,b,"."); print $0, a[1],b[1]}' | awk -v OFS="\t" '($14==$15){print $14}' | sort -u | grep -v -w -f $badgenes > $annotationfolder/copies/tandem/genes_with_tandem.DUP.$accession.exclude_badgenes.txt 
done < $accessions


wc -l $annotationfolder/copies/tandem/genes_with_tandem.DUP.*.txt

 427 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.10002.txt
   312 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.10015.txt
   320 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.10024.txt
   340 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.1741.txt
   324 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.220011.txt
   308 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.22002.txt
   330 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.22003.txt
   313 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.22004.txt
   290 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.22005.txt
   298 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.22006.txt
   321 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.22007.txt
   336 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.6024.txt
   344 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.6069.txt
   321 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.6124.txt
   333 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.6244.txt
   335 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.6909.txt
   318 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.6966.txt
   320 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.8236.txt
   326 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9075.txt
   315 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9537.txt
   318 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9543.txt
   331 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9638.txt
   318 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9728.txt
   288 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9764.txt
   318 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9888.txt
   308 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9905.txt
   325 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9981.txt


wc -l $annotationfolder/copies/tandem/genes_with_tandem.DUP.*.txt | grep exclude
   382 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.10002.exclude_badgenes.txt
   298 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.10015.exclude_badgenes.txt
   303 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.10024.exclude_badgenes.txt
   321 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.1741.exclude_badgenes.txt
   303 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.220011.exclude_badgenes.txt
   284 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.22002.exclude_badgenes.txt
   310 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.22003.exclude_badgenes.txt
   293 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.22004.exclude_badgenes.txt
   281 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.22005.exclude_badgenes.txt
   284 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.22006.exclude_badgenes.txt
   303 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.22007.exclude_badgenes.txt
   314 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.6024.exclude_badgenes.txt
   319 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.6069.exclude_badgenes.txt
   301 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.6124.exclude_badgenes.txt
   313 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.6244.exclude_badgenes.txt
   312 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.6909.exclude_badgenes.txt
   303 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.6966.exclude_badgenes.txt
   307 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.8236.exclude_badgenes.txt
   311 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9075.exclude_badgenes.txt
   298 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9537.exclude_badgenes.txt
   303 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9543.exclude_badgenes.txt
   310 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9638.exclude_badgenes.txt
   305 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9728.exclude_badgenes.txt
   272 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9764.exclude_badgenes.txt
   298 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9888.exclude_badgenes.txt
   291 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9905.exclude_badgenes.txt
   304 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.9981.exclude_badgenes.txt











###############################################################################
###### calculate TPM on own genomes
###############################################################################
ml subread/2.0.1-gcc-7.3.0-2.30
 ml bedtools/2.27.1-foss-2018b

working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus
cd $working_folder

while read accession
do
echo $accession
ls $working_folder/BAM_on_own_genome/ | grep -v table| grep $accession > $working_folder/$accession.RNAseq_samplelist.txt

while read namesample
do
echo $namesample
export bam=$working_folder/BAM_on_own_genome/$namesample/$namesample.Aligned.sortedByCoord.out.bam  
export featureCountsfolder=$working_folder/expression

export saf=$annotationfolder/mRNAs/mRNA_exons.$accession.saf
export out=$featureCountsfolder/$namesample.mRNA_exons
featureCounts -p -C  -T 8  -F SAF -O -s 2  -t exon -g gene_id -a $saf -o $out.txt $bam
# the program does not count multimapping reads! 
#convert to simpler columns
sum=`cat $out.txt| awk -v OFS="\t" '{ sum += $7} END {print sum}'`
cat $out.txt| grep -v feature | grep -v Strand| awk -v OFS="\t" -v sum="$sum" ' {print $1,$6,$7,$7*1000000*1000/($6*sum)}' > $out.counts_tpm.bed

export saf=$annotationfolder/genes/loci.$accession.saf
export out=$featureCountsfolder/$namesample.consensus_locus
featureCounts -p -T 8  -F SAF -O -s 2  -t exon -g gene_id -a $saf -o $out.txt $bam
# the program does not count multimapping reads! 
#convert to simpler columns
sum=`cat $out.txt| awk -v OFS="\t" '{ sum += $7} END {print sum}'`
cat $out.txt| grep -v feature | grep -v Strand| awk -v OFS="\t" -v sum="$sum" ' {print $1,$6,$7,$7*1000000*1000/($6*sum)}' > $out.counts_tpm.bed

export saf=$annotationfolder/SVs/svs_v06.$accession.saf
export out=$featureCountsfolder/$namesample.SVs
featureCounts -p -T 8  -F SAF -O -s 0  -t exon -g gene_id -a $saf -o $out.txt $bam
# the program does not count multimapping reads! 
#convert to simpler columns
sum=`cat $out.txt| awk -v OFS="\t" '{ sum += $7} END {print sum}'`
cat $out.txt| grep -v feature | grep -v Strand| awk -v OFS="\t" -v sum="$sum" ' {print $1,$6,$7,$7*1000000*1000/($6*sum)}' > $out.counts_tpm.bed

done < $working_folder/$accession.RNAseq_samplelist.txt
done < $accessions

##########################################
#genes and SVs together + redo SVs
##########################################
ml subread/2.0.1-gcc-7.3.0-2.30
 ml bedtools/2.27.1-foss-2018b

working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus
cd $working_folder
export gfffolder=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/gff/own
export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations
export accessions=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/accessions.txt

while read accession
do
echo $accession
while read namesample
do
echo $namesample
export bam=$working_folder/BAM_on_own_genome/$namesample/$namesample.Aligned.sortedByCoord.out.bam  
export featureCountsfolder=$working_folder/expression
cd $featureCountsfolder
export saf=$annotationfolder/loci_and_SVs.$accession.saf
export out=$featureCountsfolder/$namesample.genes_SVs
featureCounts -p -T 8  -F SAF -O -s 0  -t exon -g gene_id -a $saf -o $out.txt $bam
# the program does not count multimapping reads! 
#convert to simpler columns
sum=`cat $out.txt| awk -v OFS="\t" '{ sum += $7} END {print sum}'`
cat $out.txt| grep -v feature | grep -v Strand| awk -v OFS="\t" -v sum="$sum" ' {print $1,$6,$7,$7*1000000*1000/($6*sum)}' > $out.counts_tpm.bed


export saf=$annotationfolder/SVs/svs_v06.$accession.saf
export out=$featureCountsfolder/$namesample.SVs
featureCounts -p -T 8  -F SAF -O -s 0  -t exon -g gene_id -a $saf -o $out.txt $bam
# the program does not count multimapping reads! 
#convert to simpler columns
sum=`cat $out.txt| awk -v OFS="\t" '{ sum += $7} END {print sum}'`
cat $out.txt| grep -v feature | grep -v Strand| awk -v OFS="\t" -v sum="$sum" ' {print $1,$6,$7,$7*1000000*1000/($6*sum)}' > $out.counts_tpm.bed

done < $working_folder/$accession.RNAseq_samplelist.txt
done < $accessions




###############################################################################
#combine TPM  tables for different accessions - on OWN GENOMES
###############################################################################

working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus
cd $working_folder
export featureCountsfolder=$working_folder/expression

cd $featureCountsfolder
while read accession
do
echo $accession
echo "gene" > tpm_table_locus
cat R.$accession.consensus_locus.counts_tpm.bed  |  awk -v OFS="\t" '{print $1}' >> tpm_table_locus
echo "gene" > tpm_table_exons
cat R.$accession.mRNA_exons.counts_tpm.bed  |  awk -v OFS="\t" '{print $1}' >> tpm_table_exons
echo "gene" > tpm_table_SV
cat R.$accession.genes_SVs.counts_tpm.bed  |  awk -v OFS="\t" '{print $1}' >> tpm_table_SV

while read namesample
do
echo $namesample
echo $namesample > tpm
cat $namesample.consensus_locus.counts_tpm.bed | awk -v OFS="\t" '{print $4}' >> tpm
paste tpm_table_locus tpm > inter
cat inter > tpm_table_locus

echo $namesample > tpm
cat $namesample.mRNA_exons.counts_tpm.bed | awk -v OFS="\t" '{print $4}' >> tpm
paste tpm_table_exons tpm > inter
cat inter > tpm_table_exons

echo $namesample > tpm
cat $namesample.genes_SVs.counts_tpm.bed | awk -v OFS="\t" '{print $4}' >> tpm
paste tpm_table_SV tpm > inter
cat inter > tpm_table_SV
done < $working_folder/$accession.RNAseq_samplelist.txt

cat tpm_table_exons > $accession.TPMs.by_exons.bed
cat tpm_table_locus > $accession.TPMs.by_full_loci.bed
cat tpm_table_SV > $accession.TPMs.by_full_loci_and_SVs.bed

done < $accessions


cp *TPMs.by_* /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs


#combine count  tables for different accessions

working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus
cd $working_folder
export featureCountsfolder=$working_folder/expression

cd $featureCountsfolder
while read accession
do
echo $accession
echo "gene" > tpm_table_locus
cat R.$accession.consensus_locus.counts_tpm.bed  |  awk -v OFS="\t" '{print $1}' >> tpm_table_locus
echo "gene" > tpm_table_exons
cat R.$accession.mRNA_exons.counts_tpm.bed  |  awk -v OFS="\t" '{print $1}' >> tpm_table_exons
echo "gene" > tpm_table_SV
cat R.$accession.genes_SVs.counts_tpm.bed  |  awk -v OFS="\t" '{print $1}' >> tpm_table_SV

while read namesample
do
echo $namesample
echo $namesample > tpm
cat $namesample.consensus_locus.counts_tpm.bed | awk -v OFS="\t" '{print $3}' >> tpm
paste tpm_table_locus tpm > inter
cat inter > tpm_table_locus
echo $namesample > tpm
cat $namesample.mRNA_exons.counts_tpm.bed | awk -v OFS="\t" '{print $3}' >> tpm
paste tpm_table_exons tpm > inter
cat inter > tpm_table_exons

echo $namesample > tpm
cat $namesample.genes_SVs.counts_tpm.bed | awk -v OFS="\t" '{print $3}' >> tpm
paste tpm_table_SV tpm > inter
cat inter > tpm_table_SV
done < $working_folder/$accession.RNAseq_samplelist.txt

cat tpm_table_exons > $accession.counts.by_exons.bed
cat tpm_table_locus > $accession.counts.by_full_loci.bed
cat tpm_table_SV > $accession.counts.by_full_loci_and_SVs.bed

done < $accessions


cp $featureCountsfolder/*counts.by_* /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs










###############################################################################
###### calculate Araport genes TPM on TAIR10
###############################################################################
ml subread/2.0.1-gcc-7.3.0-2.30
ml bedtools/2.27.1-foss-2018b

#convert bed6 to saf
#cat /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.bed | awk -v OFS="\t" '{print $4, $1, $2, $3,$6}'| grep -v ChrM | grep -v ChrC > /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.loci.saf

Araport_PC_genes_loci=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.loci.saf

Araport_PC_genes_mRNAs=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.nochrM_C.saf


working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus
export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations
export accessions=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/accessions.txt

working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus
cd $working_folder

while read accession
do
echo $accession
ls $working_folder/BAM_on_TAIR10/ | grep -v table| grep $accession > $working_folder/$accession.RNAseq_samplelist.txt

while read namesample
do
echo $namesample
export bam=$working_folder/BAM_on_TAIR10/$namesample/$namesample.Aligned.sortedByCoord.out.bam  
export featureCountsfolder=$working_folder/expression/onTAIR10

export saf=$Araport_PC_genes_mRNAs
export out=$featureCountsfolder/$namesample.Araport11.mRNA_exons.on_TAIR10
featureCounts -p -C  -T 8  -F SAF -O -s 2  -t exon -g gene_id -a $saf -o $out.txt $bam
# the program does not count multimapping reads! 
#convert to simpler columns
sum=`cat $out.txt| awk -v OFS="\t" '{ sum += $7} END {print sum}'`
cat $out.txt| grep -v feature | grep -v Strand| awk -v OFS="\t" -v sum="$sum" ' {print $1,$6,$7,$7*1000000*1000/($6*sum)}' > $out.counts_tpm.bed

export saf=$Araport_PC_genes_loci
export out=$featureCountsfolder/$namesample.Araport11.locus.on_TAIR10
featureCounts -p -T 8  -F SAF -O -s 2  -t exon -g gene_id -a $saf -o $out.txt $bam
# the program does not count multimapping reads! 
#convert to simpler columns
sum=`cat $out.txt| awk -v OFS="\t" '{ sum += $7} END {print sum}'`
cat $out.txt| grep -v feature | grep -v Strand| awk -v OFS="\t" -v sum="$sum" ' {print $1,$6,$7,$7*1000000*1000/($6*sum)}' > $out.counts_tpm.bed

done < $working_folder/$accession.RNAseq_samplelist.txt
done < $accessions



###############################################################################
#combine TPM  tables for different accessions - calculation on TAIR10
###############################################################################

cat *.RNAseq_samplelist.ontair10.txt > RNAseq_samplelist.ontair10.txt

working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus
cd $working_folder
cat *.RNAseq_samplelist.ontair10.txt > RNAseq_samplelist.ontair10.txt
export featureCountsfolder=$working_folder/expression/onTAIR10

cd $featureCountsfolder

echo "gene" > tpm_table_locus
cat F.22002.batch1.Araport11.locus.on_TAIR10.counts_tpm.bed  |  awk -v OFS="\t" '{print $1}' >> tpm_table_locus
echo "gene" > tpm_table_exons
cat F.22002.batch1.Araport11.mRNA_exons.on_TAIR10.counts_tpm.bed   |  awk -v OFS="\t" '{print $1}' >> tpm_table_exons

while read namesample
do
echo $namesample

echo $namesample > tpm
cat $namesample.Araport11.locus.on_TAIR10.counts_tpm.bed| awk -v OFS="\t" '{print $4}' >> tpm
paste tpm_table_locus tpm > inter
cat inter > tpm_table_locus

echo $namesample > tpm
cat $namesample.Araport11.mRNA_exons.on_TAIR10.counts_tpm.bed| awk -v OFS="\t" '{print $4}' >> tpm
paste tpm_table_exons tpm > inter
cat inter > tpm_table_exons

done < $working_folder/RNAseq_samplelist.ontair10.txt

cat tpm_table_exons > allsamples.TPMs.by_exons.on_TAIR10.bed
cat tpm_table_locus > allsamples.TPMs.by_full_loci.on_TAIR10.bed



cp $featureCountsfolder/allsamples* /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/onTAIR10


##############################################################################
#combine count  tables for different accessions - on TAIR10
###############################################################################

working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus
cd $working_folder
export featureCountsfolder=$working_folder/expression/onTAIR10

cd $featureCountsfolder
echo "gene" > tpm_table_locus
cat F.22002.batch1.Araport11.locus.on_TAIR10.counts_tpm.bed  |  awk -v OFS="\t" '{print $1}' >> tpm_table_locus
echo "gene" > tpm_table_exons
cat F.22002.batch1.Araport11.mRNA_exons.on_TAIR10.counts_tpm.bed   |  awk -v OFS="\t" '{print $1}' >> tpm_table_exons

while read namesample
do
echo $namesample

echo $namesample > tpm
cat $namesample.Araport11.locus.on_TAIR10.counts_tpm.bed| awk -v OFS="\t" '{print $3}' >> tpm
paste tpm_table_locus tpm > inter
cat inter > tpm_table_locus

echo $namesample > tpm
cat $namesample.Araport11.mRNA_exons.on_TAIR10.counts_tpm.bed| awk -v OFS="\t" '{print $3}' >> tpm
paste tpm_table_exons tpm > inter
cat inter > tpm_table_exons

done < $working_folder/RNAseq_samplelist.ontair10.txt

cat tpm_table_exons > allsamples.counts.by_exons.on_TAIR10.bed
cat tpm_table_locus > allsamples.counts.by_full_loci.on_TAIR10.bed

cp $featureCountsfolder/allsamples* /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/onTAIR10



######################################


##################################################################
##################################################################
# calculate METHYLATION 
##################################################################
##################################################################

ml bedtools/2.27.1-foss-2018b
accessions=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/accessions_methylation1001.txt
wd=/groups/nordborg/projects/cegs/alexandra/1001Gplus/methylation
#methylation calls and alignment by Haijun 

#cp /groups/nordborg/projects/the1001genomesplus/003_methylation_calls/001.bsseq_on_pacbio/002.col0-Based/*.ref.new.tsv.gz $wd

#cp /groups/nordborg/projects/the1001genomesplus/003_methylation_calls/001.bsseq_on_pacbio/003.pacbiov2.1_Based/*.own.new.tsv.gz $wd


cd  $wd
mkdir CG
mkdir CHG
mkdir CHH

#gunzip *.gz

while read name_accession
do
echo $name_accession
index=allc_
index+=$name_accession
##############################
echo $index
# all 5 chromosomes, take only CG methylation info, reformat the file into bed6-like file 
cat  $index.own.new.tsv | awk -v OFS="\t"   '{if ($4=="CGG" || $4=="CGA" || $4=="CGC"|| $4=="CGT") print $1,$2,$2, $4, 0, $3, $5,$6,$7 }'  > CG/$index.own.CG.bed
cat  $index.ref.new.tsv | awk -v OFS="\t"   '{if ($4=="CGG" || $4=="CGA" || $4=="CGC"|| $4=="CGT")  print $1,$2,$2, $4, 0, $3, $5,$6,$7 }' |sed 's/TAIR10_//g'   > CG/$index.TAIR10.CG.bed
# all 5 chromosomes, take only CHG methylation info, reformat the file into bed6-like file 
cat  $index.ref.new.tsv	| awk -v OFS="\t"   '{if ($4=="CAG" || $4=="CCG" || $4=="CTG") print $1,$2,$2, $4, 0, $3, $5,$6,$7 }'  |sed 's/TAIR10_//g'> CHG/$index.TAIR10.CHG.bed
cat  $index.own.new.tsv	| awk -v OFS="\t"   '{if ($4=="CAG" || $4=="CCG" || $4=="CTG") print $1,$2,$2, $4, 0, $3, $5,$6,$7 }'  > CHG/$index.own.CHG.bed

cat  $index.ref.new.tsv	| awk -v OFS="\t"   '{if ($4=="CAA" || $4=="CAT" || $4=="CAC"|| $4=="CTA"|| $4=="CTT"|| $4=="CTC"|| $4=="CCA"|| $4=="CCT"|| $4=="CCC") print $1,$2,$2, $4, 0, $3, $5,$6,$7 }'|sed 's/TAIR10_//g'  > CHH/$index.TAIR10.CHH.bed
cat  $index.own.new.tsv	| awk -v OFS="\t"   '{if ($4=="CAA" || $4=="CAT" || $4=="CAC"|| $4=="CTA"|| $4=="CTT"|| $4=="CTC"|| $4=="CCA"|| $4=="CCT"|| $4=="CCC") print $1,$2,$2, $4, 0, $3, $5,$6,$7 }'  > CHH/$index.own.CHH.bed
done < $accessions

# count methylation on genes 


# LOCI - in own genomes

export gfffolder=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/gff/own
export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations
#export accessions=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/accessions.txt


ml bedtools/2.27.1-foss-2018b
accessions=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/accessions_methylation1001.txt
wd=/groups/nordborg/projects/cegs/alexandra/1001Gplus/methylation

cd $wd

# genes from annotations 
while read name_accession
do
echo $name_accession
index=allc_
index+=$name_accession
##############################
echo $index
#genes 
genes=$annotationfolder/genes/loci.$name_accession.bed
sortBed -i $genes |bedtools map -a stdin -b CG/$index.own.CG.bed  -c 7,8 -o sum -null "NA" |awk -v OFS="\t" '{if($7!="NA"&&$8!="NA") print $1,$2,$3,$4,$5,$6,$7/$8;else print $1,$2,$3,$4,$5,$6,"NA"}' > CG/$name_accession.CG.own.loci.bed

sortBed -i $genes |bedtools map -a stdin -b CHG/$index.own.CHG.bed  -c 7,8 -o sum -null "NA" |awk -v OFS="\t" '{if($7!="NA"&&$8!="NA") print $1,$2,$3,$4,$5,$6,$7/$8;else print $1,$2,$3,$4,$5,$6,"NA"}' > CHG/$name_accession.CHG.own.loci.bed

sortBed -i $genes |bedtools map -a stdin -b CHH/$index.own.CHH.bed  -c 7,8 -o sum -null "NA" |awk -v OFS="\t" '{if($7!="NA"&&$8!="NA") print $1,$2,$3,$4,$5,$6,$7/$8;else print $1,$2,$3,$4,$5,$6,"NA"}' > CHH/$name_accession.CHH.own.loci.bed

done < $accessions

# LOCI on TAIR10
 
ml bedtools/2.27.1-foss-2018b
accessions=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/accessions_methylation1001.txt
wd=/groups/nordborg/projects/cegs/alexandra/1001Gplus/methylation

wc -l /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.bed
#27129 /groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.bed

cd $wd

while read name_accession
do
echo $name_accession
index=allc_
index+=$name_accession
##############################
echo $index
#genes 
genes=/groups/nordborg/user/aleksandra.kornienko/analyses/Annotation/Araport11_protein_coding.201606.genes.bed
sortBed -i $genes |bedtools map -a stdin -b CG/$index.TAIR10.CG.bed  -c 7,8 -o sum -null "NA" |awk -v OFS="\t" '{if($7!="NA"&&$8!="NA") print $1,$2,$3,$4,$5,$6,$7/$8;else print $1,$2,$3,$4,$5,$6,"NA"}' > CG/$name_accession.CG.TAIR10.loci.bed

sortBed -i $genes |bedtools map -a stdin -b CHG/$index.TAIR10.CHG.bed  -c 7,8 -o sum -null "NA" |awk -v OFS="\t" '{if($7!="NA"&&$8!="NA") print $1,$2,$3,$4,$5,$6,$7/$8;else print $1,$2,$3,$4,$5,$6,"NA"}' > CHG/$name_accession.CHG.TAIR10.loci.bed

sortBed -i $genes |bedtools map -a stdin -b CHH/$index.TAIR10.CHH.bed  -c 7,8 -o sum -null "NA" |awk -v OFS="\t" '{if($7!="NA"&&$8!="NA") print $1,$2,$3,$4,$5,$6,$7/$8;else print $1,$2,$3,$4,$5,$6,"NA"}' > CHH/$name_accession.CHH.TAIR10.loci.bed

done < $accessions


#SVs on own genomes

ml bedtools/2.27.1-foss-2018b

export gfffolder=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/gff/own
export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations

accessions=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/accessions_methylation1001.txt
wd=/groups/nordborg/projects/cegs/alexandra/1001Gplus/methylation

cd $wd

while read name_accession
do
echo $name_accession
index=allc_
index+=$name_accession
##############################
echo $index
#SVs
SVs=$annotationfolder/SVs/svs_v06.$name_accession.bed
sortBed -i $SVs |bedtools map -a stdin -b CG/$index.own.CG.bed  -c 7,8 -o sum -null "NA" |awk -v OFS="\t" '{if($7!="NA"&&$8!="NA") print $1,$2,$3,$4,$5,$6,$7/$8;else print $1,$2,$3,$4,$5,$6,"NA"}' > CG/$name_accession.CG.own.SVs.bed

sortBed -i $SVs |bedtools map -a stdin -b CHG/$index.own.CHG.bed  -c 7,8 -o sum -null "NA" |awk -v OFS="\t" '{if($7!="NA"&&$8!="NA") print $1,$2,$3,$4,$5,$6,$7/$8;else print $1,$2,$3,$4,$5,$6,"NA"}'  > CHG/$name_accession.CHG.own.SVs.bed

sortBed -i $SVs |bedtools map -a stdin -b CHH/$index.own.CHH.bed  -c 7,8 -o sum -null "NA" |awk -v OFS="\t" '{if($7!="NA"&&$8!="NA") print $1,$2,$3,$4,$5,$6,$7/$8;else print $1,$2,$3,$4,$5,$6,"NA"}'  > CHH/$name_accession.CHH.own.SVs.bed

done < $accessions


cp /groups/nordborg/projects/cegs/alexandra/1001Gplus/methylation/C*/*.loci.bed  /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation


cp /groups/nordborg/projects/cegs/alexandra/1001Gplus/methylation/C*/*.SVs.bed  /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation


##################################################################
##################################################################
# calculate HISTONE modification CHIPseq 
##################################################################
##################################################################

ml bedtools/2.27.1-foss-2018b
### chipseq 

#alignment and calculation for genes and older SV annotation
#Z:\01_POSTDOC\03_Projects\ERA-CAPS\20230410_ChIPseq_6acc_normalize_calculate_coverage.bash
#combine coverage into one table 
#Z:\01_POSTDOC\03_Projects\ERA-CAPS\20230410_combine_chip_coverage_eracaps_onowngenome.bash

ml bedtools/2.27.1-foss-2018b

#list of samples
list=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/ChIP_samples_6acc_eracaps_names_only_antibodies.txt
# wc -l /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/ChIP_samples_6acc_eracaps_names_only_antibodies.txt
#66 /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/ChIP_samples_6acc_eracaps_names_only_antibodies.txt

export annotationfolder=/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations

for numberofline in {1..66}
do
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

#bamCompare -b1 $sample.Aligned.sortedByCoord.out.bam -b2 $working_folder/BAM_on_own_genome/$input/$input.Aligned.sortedByCoord.out.bam   --outFileName  $sample.log2.input_norm.bedGraph --operation log2 --effectiveGenomeSize 119481543 -p 4 --ignoreDuplicates --outFileFormat bedgraph
#bamCompare -b1 $sample.Aligned.sortedByCoord.out.bam -b2  $working_folder/BAM_on_own_genome/$input/$input.Aligned.sortedByCoord.out.bam --outFileName   $sample.log2.input_norm.bw --operation log2 --effectiveGenomeSize 119481543 -p 4 --ignoreDuplicates

export SV=$annotationfolder/SVs/svs_v06.$accession.bed 
export locus=$annotationfolder/genes/loci.$accession.bed
export mRNA=$annotationfolder/mRNAs/mRNA.$accession.bed

# calculate chipseq coverage 
bedtools sort -i $SV | bedtools map -a stdin -b $working_folder/BAM_on_own_genome/$sample/$sample.log2.input_norm.bedGraph -c 4 -o mean| awk -v OFS="\t"  '{print $1,$2,$3,$4,$5,$6,$7}'  > $working_folder/BAM_on_own_genome/coverage_on_own_genome/$sample.svs_v06.log2.mean_cov.bed
bedtools sort -i $locus | bedtools map -a stdin -b $working_folder/BAM_on_own_genome/$sample/$sample.log2.input_norm.bedGraph -c 4 -o mean| awk -v OFS="\t"  '{print $1,$2,$3,$4,$5,$6,$7}'  > $working_folder/BAM_on_own_genome/coverage_on_own_genome/$sample.locus.log2.mean_cov.bed
bedtools sort -i $mRNA | bedtools map -a stdin -b $working_folder/BAM_on_own_genome/$sample/$sample.log2.input_norm.bedGraph -c 4 -o mean| awk -v OFS="\t"  '{print $1,$2,$3,$4,$5,$6,$7}'  > $working_folder/BAM_on_own_genome/coverage_on_own_genome/$sample.mRNA.log2.mean_cov.bed
done

#################

# combine chipseq for each accession
cd /groups/nordborg/projects/cegs/alexandra/1001Gplus/Chipseq/BAM_on_own_genome/coverage_on_own_genome

while read accession
do 
echo $accession

cat /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/ChIP_samples_6acc_eracaps_names_only_antibodies.txt | grep $accession |awk -v OFS="\t" '{print $1}'  >tmp 

sample1=`head -1 tmp `

#log2
echo -e "chr\t" "start\t" "end\t" "SV_name\t" "SV_info" > cov_SV_table
cat $sample1.svs_v06.log2.mean_cov.bed|  grep Chr |  awk -v OFS="\t" '{print $1,$2,$3,$4,$5}'  >> cov_SV_table

echo -e "chr\t" "start\t" "end\t" "group\t" "strand" > cov_loci_table
cat $sample1.locus.log2.mean_cov.bed|  grep Chr |  awk -v OFS="\t" '{print $1,$2,$3,$4,$6}'  >> cov_loci_table

echo -e "chr\t" "start\t" "end\t" "group\t" "strand"  > cov_mrna_table
cat $sample1.mRNA.log2.mean_cov.bed|  grep Chr |  awk -v OFS="\t" '{print $1,$2,$3,$4,$6}'  >> cov_mrna_table

while read name 
do 
echo $name> cov_mrna
echo $name> cov_loci
echo $name> cov_SV

cat $name.svs_v06.log2.mean_cov.bed |  grep Chr |  awk -v OFS="\t" '{print $7}'  >>cov_SV
paste cov_SV_table cov_SV > inter 
cat inter > cov_SV_table

cat $name.locus.log2.mean_cov.bed |  grep Chr |  awk -v OFS="\t" '{print $7}'  >>cov_loci
paste cov_loci_table cov_loci > inter 
cat inter > cov_loci_table

cat $name.mRNA.log2.mean_cov.bed |  grep Chr |  awk -v OFS="\t" '{print $7}'  >>cov_mrna
paste cov_mrna_table cov_mrna > inter 
cat inter > cov_mrna_table

done < tmp

cat cov_SV_table > Chipseq_coverage.owngenome.$accession.svs_v06.log2.bed
cat cov_loci_table > Chipseq_coverage.owngenome.$accession.loci.log2.bed
cat cov_mrna_table > Chipseq_coverage.owngenome.$accession.mRNAs.log2.bed

done < /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/ChIP_samples_6acc_eracaps.txt


 cp Chipseq_cov* /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/


##### miRNA seq 
#calculate coverage


export working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus/miRNA
cd $working_folder


while read name_accession
do
echo $name_accession

export SV=$annotationfolder/SVs/svs_v06.$name_accession.bed 
export locus=$annotationfolder/genes/loci.$name_accession.bed
export mRNA=$annotationfolder/mRNAs/mRNA.$name_accession.bed

bedtools sort -i $locus | bedtools map -a stdin -b $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.perbase.bed -c 4 -o mean | sort  -k4,4 > $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.loci.coverage.perbase_calc.bed

bedtools sort -i $mRNA | bedtools map -a stdin -b $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.perbase.bed -c 4 -o mean | sort  -k4,4 > $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.mRNAs.coverage.perbase_calc.bed

bedtools sort -i $SV | bedtools map -a stdin -b $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.perbase.bed -c 4 -o mean | sort  -k4,4 > $working_folder/BAM_on_own_genome/$name_accession/$name_accession.24nt.SVs.coverage.perbase_calc.bed

done < /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/accessions_miRNAseq.txt


cp $working_folder/BAM_on_own_genome/*/*.24nt.*.coverage.perbase_calc.bed /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/mirna



##########################
# combine read coverage tables 

# RNAseq 


#miRNAseq 


export working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus/miRNA
cd $working_folder

echo "accession" "tissue" "total_readN" "uniquely_mapped_N" "uniquely_mapped_%" "multimapped_%" "unmapped_too_many_loci_%" "unmapped_too_many_mismatches_%" "unmapped_too_short_%" "unmapped_other_%"  > /groups/nordborg/projects/cegs/alexandra/1001Gplus/miRNA/readN.miRNA.txt
while read name_accession
do
echo $name_accession

total_readN=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "Number of input reads" | awk '{print $6}' `
uniq_readN=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "Uniquely mapped reads number" | awk '{print $6}' `
uniq_perc=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "Uniquely mapped reads %" | awk '{print $6}' `
multimap=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "% of reads mapped to multiple loci" | awk '{print $9}' `
toomany=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "% of reads mapped to too many loci" | awk '{print $10}' `
unmap_mism=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "% of reads unmapped: too many mismatches" | awk '{print $9}' `
unmap_short=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "% of reads unmapped: too short" | awk '{print $8}' `
unmap_other=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "% of reads unmapped: other" | awk '{print $7}' `


echo $name_accession "flowers" $total_readN $uniq_readN $uniq_perc  $multimap $toomany $unmap_mism $unmap_short $unmap_other >> /groups/nordborg/projects/cegs/alexandra/1001Gplus/miRNA/readN.miRNA.txt

done < /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/accessions_miRNAseq.txt


cp /groups/nordborg/projects/cegs/alexandra/1001Gplus/miRNA/readN.miRNA.txt /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/




##### chip seq 

export working_folder=/groups/nordborg/projects/cegs/alexandra/1001Gplus/Chipseq 

cd $working_folder/BAM_on_own_genome/$sample


echo "accession" "tissue" "total_readN" "uniquely_mapped_N" "uniquely_mapped_%" "multimapped_%" "unmapped_too_many_loci_%" "unmapped_too_many_mismatches_%" "unmapped_too_short_%" "unmapped_other_%"  > $working_folder/readN.chipseq.txt
while read name_accession
do
echo $name_accession

#cat 6909.Log.final.out | grep "Number of input reads" | awk '{print $6}'


total_readN=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "Number of input reads" | awk '{print $6}' `
uniq_readN=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "Uniquely mapped reads number" | awk '{print $6}' `
uniq_perc=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "Uniquely mapped reads %" | awk '{print $6}' `
multimap=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "% of reads mapped to multiple loci" | awk '{print $9}' `
toomany=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "% of reads mapped to too many loci" | awk '{print $10}' `
unmap_mism=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "% of reads unmapped: too many mismatches" | awk '{print $9}' `
unmap_short=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "% of reads unmapped: too short" | awk '{print $8}' `
unmap_other=`cat $working_folder/BAM_on_own_genome/$name_accession/$name_accession.Log.final.out  | grep "% of reads unmapped: other" | awk '{print $7}' `


echo $name_accession $total_readN $uniq_readN $uniq_perc  $multimap $toomany $unmap_mism $unmap_short $unmap_other >> $working_folder/readN.chipseq.txt

echo $name_accession  $total_readN $uniq_readN $uniq_perc  $multimap $toomany $unmap_mism $unmap_short $unmap_other
done < /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ChIP-seq_data/ChIP_samplenames_66.txt


cp $working_folder/readN.chipseq.txt /groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/


