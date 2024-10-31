# USAGE:  # ./script_overlap_REPEATS_EDTA_contigs.sh <input_EDTA> <input_REPEATS> <ACC>
# EXAMPLE1 # ./script_overlap_REPEATS_EDTA_contigs.sh /ebio/abt6_projects8/1001g_plus_scaffolding/data/RagTag/3_EDTA_contigs_v2.1 /ebio/abt6_projects8/1001g_plus_scaffolding/data/RagTag/3_RepeatAnnotation_contigs_v2.1 22001

# Parse parameteres
#inputDIR_EDTA=$1
#inputDIR_REPEATS=$2
#ACC=$3

# Declare other variables
#WORKTMP=/tmp/global2/frabanal/${WORKDIR/\/ebio\/abt6_projects7\/}

#inputDIR=$(echo "${input_REPEATS}" | rev | cut --complement -d'/' -f1 | rev)
#inputTMP=/tmp/global2/frabanal/${inputDIR_REPEATS/\/ebio\/abt6_projects8\/}


# Required Software
BEDTOOLS=/ebio/abt6_projects9/abt6_software/bin/bedtools/bin/bedtools

# Start
echo ""
echo "START script"
date

########################################
############ Simplify TEs ##############
########################################

echo "Soft linking to original file..."
ln -s $inputDIR_EDTA/$ACC/$ACC.polished_contigs.new_names.v2.1.fasta.mod.EDTA.TEanno.split.gff3 $inputDIR_REPEATS/$ACC/

echo "Removing signatures of LTRs..."
grep -v "ID=repeat_" $inputDIR_EDTA/$ACC/$ACC.polished_contigs.new_names.v2.1.fasta.mod.EDTA.TEanno.split.gff3 > $inputDIR_REPEATS/$ACC/$ACC.EDTA.TEanno.edit.gff3
#grep -v "long_terminal_repeat" $inputDIR_REPEATS/$ACC/$ACC.EDTA.TMP1.gff3 > $inputDIR_REPEATS/$ACC/$ACC.EDTA.TMP2.gff3
#grep -v "target_site_duplication" $inputDIR_REPEATS/$ACC/$ACC.EDTA.TMP2.gff3 > $inputDIR_REPEATS/$ACC/$ACC.EDTA.TEanno.edit.gff3

echo "Classifying TEs..."
sed -i 's/Copia_LTR_retrotransposon/ClassI_LTR/' $inputDIR_REPEATS/$ACC/$ACC.EDTA.TEanno.edit.gff3
sed -i 's/Gypsy_LTR_retrotransposon/ClassI_LTR/' $inputDIR_REPEATS/$ACC/$ACC.EDTA.TEanno.edit.gff3
sed -i 's/LTR_retrotransposon/ClassI_LTR/' $inputDIR_REPEATS/$ACC/$ACC.EDTA.TEanno.edit.gff3

sed -i 's/LINE_element/ClassI_nonLTR/' $inputDIR_REPEATS/$ACC/$ACC.EDTA.TEanno.edit.gff3

sed -i 's/CACTA_TIR_transposon/ClassII_TIR/' $inputDIR_REPEATS/$ACC/$ACC.EDTA.TEanno.edit.gff3
sed -i 's/hAT_TIR_transposon/ClassII_TIR/' $inputDIR_REPEATS/$ACC/$ACC.EDTA.TEanno.edit.gff3
sed -i 's/Mutator_TIR_transposon/ClassII_TIR/' $inputDIR_REPEATS/$ACC/$ACC.EDTA.TEanno.edit.gff3
sed -i 's/PIF_Harbinger_TIR_transposon/ClassII_TIR/' $inputDIR_REPEATS/$ACC/$ACC.EDTA.TEanno.edit.gff3
sed -i 's/Tc1_Mariner_TIR_transposon/ClassII_TIR/' $inputDIR_REPEATS/$ACC/$ACC.EDTA.TEanno.edit.gff3

sed -i 's/helitron/ClassII_Helitrons/' $inputDIR_REPEATS/$ACC/$ACC.EDTA.TEanno.edit.gff3

#rm $inputDIR_REPEATS/$ACC/$ACC.EDTA.TMP*

##########################################
############ EDTA + REPEATS ##############
##########################################

input_EDTA=$inputDIR_REPEATS/$ACC/$ACC.EDTA.TEanno.edit.gff3
input_REPEATS=$inputDIR_REPEATS/$ACC/$ACC.Repeats_merged.gff

echo "# Pre-processing EDTA.gff to filter out 'TEs' overlapping rRNAs..."
cd $inputDIR_REPEATS/$ACC

# "Only report those entries in A that have no overlap in B. Restricted by -f and -r."
$BEDTOOLS intersect -v -a $input_EDTA -b $input_REPEATS > EDTA_rDNAPurged.v.gff

# "Perform a “left outer join”. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B." 
$BEDTOOLS intersect -loj -a $input_EDTA -b $input_REPEATS > EDTA_rDNAPurged.loj.gff

# Concatenate Repeats and TE annotations (only those that do not overlap) 
cat $input_REPEATS EDTA_rDNAPurged.v.gff > $ACC.Repeats_TEanno.TMP1.gff3

# Remove header and sort
grep -v "^#"  $ACC.Repeats_TEanno.TMP1.gff3 > $ACC.Repeats_TEanno.TMP2.gff3
sort -k1,1 -k4n $ACC.Repeats_TEanno.TMP2.gff3 > $ACC.Repeats_TEanno.gff3

sed -i '1i ##gff-version 3' $ACC.Repeats_TEanno.gff3
sed -i "2i ##date $(date +%Y-%m-%d)" $ACC.Repeats_TEanno.gff3

# BEFORE and AFTER TE count
echo "TEs BEFORE purging:"
grep -c -v "^#" $input_EDTA
echo "TEs AFTER purging:"
grep -c -v "^#" EDTA_rDNAPurged.v.gff

rm $ACC.Repeats_TEanno.TMP*

echo "Transferring relevant files to project directory..."
#rsync -av EDTA_rDNAPurged.* $inputDIR/
#rsync -av $ACC.Repeats_TEanno.gff3 $inputDIR/

date
echo "END script"
echo ""
