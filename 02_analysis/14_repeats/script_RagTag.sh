# USAGE:  # ./script_RagTag.sh <WORKDIR> <CONTIGS> <ACC> <CORES> <CUTOFF> <Q> <LENGTH> <GROUPING>
# EXAMPLE # ./script_RagTag.sh /ebio/abt6_projects8/1001g_plus_scaffolding/data/RagTag /ebio/abt6_projects8/1001g_plus_data_transfer/data/1001g_plus_assemblies/freeze/v1/Elh-2/Elh-2.polished_contigs.v1.fasta 22004 8 150000 60 10000 0.6

# Parse parameteres
WORKDIR=$1
CONTIGS=$2
ACC=$3
CORES=$4
CUTOFF=$5
Q=$6
LENGTH=$7
GROUPING=$8

# Declare other variables
WORKTMP=/tmp/global2/frabanal/${WORKDIR/\/ebio\/abt6_projects8\/}
cutoff=$(( CUTOFF / 1000 ))kb

# In this scrip, I'm fixing the REFERENCE to not having to specify it as an input parameter
MASKED=yes
if [ $MASKED == 'no' ];then
	REF=/ebio/abt6_projects7/small_projects/frabanal/reference_genomes/fasta_files/at.fa
	NAME=unmasked_"$cutoff"_q"$Q"_f"$LENGTH"_i"$GROUPING"
elif [ $MASKED == 'yes' ];then
	REF=/tmp/global2/frabanal/1001g_plus_scaffolding/data/RagTag_testing/1_masked_TAIR10/TAIR10.hard_masked.fa
	NAME=masked_"$cutoff"_q"$Q"_f"$LENGTH"_i"$GROUPING"
fi

outputTMP=$WORKTMP/2_Ragtag_$NAME/$ACC
output=$WORKDIR/2_Ragtag_$NAME/$ACC

# Start
echo ""
echo "Start script"
date

# Remove previous results
rm -rf $outputTMP
rm -rf $output

# Create necessary directories
mkdir -p $output
mkdir -p $outputTMP

# Activate virtual environment
source activate /ebio/abt6_projects8/alopecurus_genome/bin/anaconda2/envs/RagTag_v1.1.1
#source activate /ebio/abt6_projects8/alopecurus_genome/bin/anaconda2/envs/RagTag_v2.0.1
SAMTOOLS=/ebio/abt6_projects9/abt6_software/bin/samtools-1.9/bin/samtools

# Scaffolds with RagTag
echo "Soft link CONTIGS..."
ln -s $CONTIGS $outputTMP/$ACC.contigs.fa

echo "Filtering CONTIGS by size..."
perl $PWD/keep_SMALL_contigs.pl $CUTOFF $outputTMP/$ACC.contigs.fa > $outputTMP/$ACC.contigs.SMALL.$cutoff.fa 
perl $PWD/keep_LARGE_contigs.pl $CUTOFF $outputTMP/$ACC.contigs.fa > $outputTMP/$ACC.contigs.LARGE.$cutoff.fa 

echo "Scaffolding MASKED..."
ragtag.py scaffold $REF $outputTMP/$ACC.contigs.LARGE.$cutoff.fa \
	-t $CORES \
	-o $outputTMP \
	-q $Q -f $LENGTH -i $GROUPING \
	--remove-small \
	--debug

mv $outputTMP/ragtag.scaffolds.fasta $outputTMP/$ACC.ragtag_scaffolds.$cutoff.fa
mv $outputTMP/ragtag.scaffolds.stats $outputTMP/$ACC.ragtag_scaffolds.$cutoff.stats

echo "Creating joint file for analysis..."
grep -e "RagTag" $outputTMP/ragtag.scaffolds.agp | grep -e "W" > $outputTMP/TMP.1.txt
awk "NR>1" $outputTMP/ragtag.confidence.txt > $outputTMP/TMP.2.txt
paste $outputTMP/TMP.1.txt $outputTMP/TMP.2.txt > $outputTMP/$ACC.ragtag_scaffolds.$cutoff.txt

rm $outputTMP/TMP.*.txt

echo "Indexing RagTag REF..."
sed -i "s/_RagTag//g" $outputTMP/$ACC.ragtag_scaffolds.$cutoff.fa
sed -i 's/Chr/'$ACC'_Chr/g' $outputTMP/$ACC.ragtag_scaffolds.$cutoff.fa
$SAMTOOLS faidx $outputTMP/$ACC.ragtag_scaffolds.$cutoff.fa

echo "Extracting only scaffolded chromosomes..."
$SAMTOOLS faidx $outputTMP/$ACC.ragtag_scaffolds.$cutoff.fa "$ACC"_Chr1 "$ACC"_Chr2 "$ACC"_Chr3 "$ACC"_Chr4 "$ACC"_Chr5 > $outputTMP/$ACC.ragtag_scaffolds.Chr.fa
$SAMTOOLS faidx $outputTMP/$ACC.ragtag_scaffolds.Chr.fa

echo "Concatenating and indexing scaffolds plus SMALL contigs..."
cat $outputTMP/$ACC.ragtag_scaffolds.$cutoff.fa $outputTMP/$ACC.contigs.SMALL.$cutoff.fa > $outputTMP/$ACC.ragtag_scaffolds.fa
$SAMTOOLS faidx $outputTMP/$ACC.ragtag_scaffolds.fa

conda deactivate
# Transfer relevant data
rsync -av $outputTMP/$ACC.ragtag_scaffolds.$cutoff.* $output/
rsync -av $outputTMP/$ACC.*.fa* $output/

##############################################
################### MiniTV ################### 

# Declare other variables
#ChrM=/ebio/abt6_projects7/small_projects/frabanal/reference_genomes/fasta_files/at_mitochondria.fa
#ChrC=/ebio/abt6_projects7/small_projects/frabanal/reference_genomes/fasta_files/at_chloroplast.fa
TAIR10=/ebio/abt6_projects7/small_projects/frabanal/reference_genomes/fasta_files/at.fa
TAIR_Repeats=/ebio/abt6_projects8/1001g_plus_scaffolding/data/annotation_repeats/contigs_v1.1/TAIR10/TAIR10.Repeats.gff

acc_bionano=(6069 6966 8236 9638 9764 9994 10024 22001 22006)
hybrid_dir=/tmp/global2/frabanal/1001g_plus_scaffolding/data/optical_maps/1_repeats

fake_gff=/ebio/abt6_projects8/1001g_plus_data_transfer/data/QC_contigs/fake.gff
CONTIG_gff=/ebio/abt6_projects8/1001g_plus_data_transfer/data/QC_scaffolding_v1.1/output_repeats_contigs/$ACC/$ACC.Repeats.gff

# Activate virtual environment
source activate /ebio/abt6_projects8/alopecurus_genome/bin/anaconda2/envs/MiniTV

# NO hybrid
echo "Aligning SCAFFOLDS for MiniTV..."

minitv --aligner_args "-x asm5 -g 150 -t $CORES" --no_tree \
	--annotation-types gene,N_stretch,rRNA,tRNA,centromere,telomere,organelle \
	-a $TAIR_Repeats \
	-a $CONTIG_gff \
	-a $fake_gff \
	$TAIR10 \
	$outputTMP/$ACC.contigs.fa \
	$outputTMP/$ACC.ragtag_scaffolds.$cutoff.fa \
	> $outputTMP/$ACC.ragtag_contig.alitv.json

# Colour of chromosomes
sed -i 's/1d91c0/d9d9d9/g' $outputTMP/$ACC.ragtag_contig.alitv.json
# Modify "Distance between ticks (bp)"
sed -i -e 's/tickDistance": 10000/tickDistance": 100000/g' $outputTMP/$ACC.ragtag_contig.alitv.json
# Modify "Spacer between chromosomes (px)"
sed -i -e 's/karyoDistance": 2000/karyoDistance": 500000/g' $outputTMP/$ACC.ragtag_contig.alitv.json
# Modify "Width" of the canvas
sed -i -e 's/canvasWidth": 900/canvasWidth": 2200/g' $outputTMP/$ACC.ragtag_contig.alitv.json
# Modify "Height" of the canvas
sed -i -e 's/canvasHeight": 900/canvasHeight": 900/g' $outputTMP/$ACC.ragtag_contig.alitv.json
# Modify "color" of 'gene' type
sed -i -e '0,/E2EDFF/ s/E2EDFF/2b8cbe/' $outputTMP/$ACC.ragtag_contig.alitv.json
# Modify "color" of 'N_stretch' type
sed -i -e '0,/E2EDFF/ s/E2EDFF/fed976/' $outputTMP/$ACC.ragtag_contig.alitv.json
# Modify "color" of 'rRNA' type
sed -i -e '0,/E2EDFF/ s/E2EDFF/b2182b/' $outputTMP/$ACC.ragtag_contig.alitv.json
# Modify "color" of 'tRNA' type
sed -i -e '0,/E2EDFF/ s/E2EDFF/dd1c77/' $outputTMP/$ACC.ragtag_contig.alitv.json
# Modify "color" of 'centromere' type
sed -i -e '0,/E2EDFF/ s/E2EDFF/252525/' $outputTMP/$ACC.ragtag_contig.alitv.json
# Modify "color" of 'telomere' type
sed -i -e '0,/E2EDFF/ s/E2EDFF/8856a7/' $outputTMP/$ACC.ragtag_contig.alitv.json
# Modify "color" of 'organelle' type
sed -i -e '0,/E2EDFF/ s/E2EDFF/35978f/' $outputTMP/$ACC.ragtag_contig.alitv.json

# With HYBRID

if [[ " ${acc_bionano[@]} " =~ " ${ACC} " ]]; then
	
	echo "Aligning SCAFFOLDS (when HYBRID present)  for MiniTV..."
	
	minitv --aligner_args "-x asm5 -g 150 -t $CORES" --no_tree \
		--annotation-types gene,N_stretch,rRNA,tRNA,centromere,telomere,organelle \
		-a $TAIR_Repeats \
		-a $CONTIG_gff \
		-a $fake_gff \
		-a $hybrid_dir/$ACC/$ACC.Repeats.HYBRID.gff \
		$TAIR10 \
		$outputTMP/$ACC.contigs.fa \
		$outputTMP/$ACC.ragtag_scaffolds.$cutoff.fa \
		$hybrid_dir/$ACC/$ACC.fa \
		> $outputTMP/$ACC.ragtag_hybrid_contig.alitv.json

	echo "Editing json file for Mitochondria..."

	# Colour of chromosomes
	sed -i 's/1d91c0/d9d9d9/g' $outputTMP/$ACC.ragtag_hybrid_contig.alitv.json
	# Modify "Distance between ticks (bp)"
	sed -i -e 's/tickDistance": 10000/tickDistance": 100000/g' $outputTMP/$ACC.ragtag_hybrid_contig.alitv.json
	# Modify "Spacer between chromosomes (px)"
	sed -i -e 's/karyoDistance": 2000/karyoDistance": 500000/g' $outputTMP/$ACC.ragtag_hybrid_contig.alitv.json
	# Modify "Width" of the canvas
	sed -i -e 's/canvasWidth": 900/canvasWidth": 2200/g' $outputTMP/$ACC.ragtag_hybrid_contig.alitv.json
	# Modify "Height" of the canvas
	sed -i -e 's/canvasHeight": 900/canvasHeight": 900/g' $outputTMP/$ACC.ragtag_hybrid_contig.alitv.json
	# Modify "color" of 'gene' type
	sed -i -e '0,/E2EDFF/ s/E2EDFF/2b8cbe/' $outputTMP/$ACC.ragtag_hybrid_contig.alitv.json
	# Modify "color" of 'N_stretch' type
	sed -i -e '0,/E2EDFF/ s/E2EDFF/fed976/' $outputTMP/$ACC.ragtag_hybrid_contig.alitv.json
	# Modify "color" of 'rRNA' type
	sed -i -e '0,/E2EDFF/ s/E2EDFF/b2182b/' $outputTMP/$ACC.ragtag_hybrid_contig.alitv.json
	# Modify "color" of 'tRNA' type
	sed -i -e '0,/E2EDFF/ s/E2EDFF/dd1c77/' $outputTMP/$ACC.ragtag_hybrid_contig.alitv.json
	# Modify "color" of 'centromere' type
	sed -i -e '0,/E2EDFF/ s/E2EDFF/252525/' $outputTMP/$ACC.ragtag_hybrid_contig.alitv.json
	# Modify "color" of 'telomere' type
	sed -i -e '0,/E2EDFF/ s/E2EDFF/8856a7/' $outputTMP/$ACC.ragtag_hybrid_contig.alitv.json
	# Modify "color" of 'organelle' type
	sed -i -e '0,/E2EDFF/ s/E2EDFF/35978f/' $outputTMP/$ACC.ragtag_hybrid_contig.alitv.json
else
	echo "There is NO HYBRID-SCAFFOLD for accession $ACC"
fi

conda deactivate

# Transfer relevant files to project directoies
rsync -av $outputTMP/$ACC.*.alitv.json $output

date
echo "END script"
echo ""

