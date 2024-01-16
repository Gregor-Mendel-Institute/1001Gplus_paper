# Workflow for modification of genomes 
This is the manual way how to do the rearrangement. 


```
cd /tmp/global2/svorbrugg/1001GP/scripts/1001GP/mod_genome/data/raw/
seqkit grep -n -r -p "Chr3" 22001.scaffolds_corrected.v2.1.fasta > 22001.scaffolds_corrected.v2.1.chr3.fasta
seqkit grep -n -r -p "Chr5" 22001.scaffolds_corrected.v2.1.fasta > 22001.scaffolds_corrected.v2.1.chr5.fasta


```
Chr3 
```
cd /ebio/abt6_projects9/1001g_plus_pan/data/data/v2.1/scaffolds
for x in *; do seqkit grep -n -r -p "Chr3" $x | seqkit grep -v -n -r -p "22001" -  > ~/tmp_global2/1001GP/scripts/1001GP/mod_genome/data/other/chr3/$x.fasta; done 

cd ~/tmp_global2/1001GP/scripts/1001GP/mod_genome/data/other/chr3/
for x in *; do minimap2 -x asm5 $x ../../raw/22001.scaffolds_corrected.v2.1.chr3.fasta > ../../paf/chr3/$x.3.paf; done

cd ../../paf/chr3/
for x in *; do cat $x | fpa drop -l 50000 > ../../pafr/chr3/$x.drop.paf; done  

cd ../../
cat pafr/chr3/* > chr3.sum.txt

./../make_bed.py -p chr3.sum.txt --min > min.bed 
```


Chr5
```
cd /ebio/abt6_projects9/1001g_plus_pan/data/data/v2.1/scaffolds
for x in *; do seqkit grep -n -r -p "Chr5" $x | seqkit grep -v -n -r -p "22001" -  > ~/tmp_global2/1001GP/scripts/1001GP/mod_genome/data/other/chr5/$x.fasta; done 

cd ~/tmp_global2/1001GP/scripts/1001GP/mod_genome/data/other/chr5/
for x in *; do minimap2 -x asm5 $x ../../raw/22001.scaffolds_corrected.v2.1.chr5.fasta > ../../paf/chr5/$x.5.paf; done

cd ../../paf/chr5/
for x in *; do cat $x | fpa drop -l 50000 > ../../pafr/chr5/$x.drop.paf; done  

cd ../../
cat pafr/chr5/* > chr5.sum.txt

./../make_bed.py -p chr5.sum.txt --max > max.bed 
```

```das
cd ~/tmp_global2/1001GP/scripts/1001GP/mod_genome/data/
cat min.bed max.bed > ./confirmation/mm.bed 
cd confirmation
./../../remove_and_add.py -f /ebio/abt6_projects9/1001g_plus_pan/data/data/v2.1/scaffolds/22001.scaffolds_corrected.v2.1.fasta -u mm.bed -o 22001_mod.scaffolds_corrected.v2.1.fasta
minimap2 -x asm5 /ebio/abt6_projects9/1001g_plus_pan/data/data/v2.1/scaffolds/22002.scaffolds_corrected.v2.1.fasta 22001_mod.scaffolds_corrected.v2.1.fasta > mod.paf
cat mod.paf | fpa drop -l 50000 > mod.50k.paf
~/serverbin/miniasm/minidot mod.50k.paf > mod.22002.50k.eps


minimap2 -x asm5 /ebio/abt6_projects9/1001g_plus_pan/data/data/v2.1/scaffolds/22002.scaffolds_corrected.v2.1.fasta /ebio/abt6_projects9/1001g_plus_pan/data/data/v2.1/scaffolds/22001.scaffolds_corrected.v2.1.fasta > mod_before.paf
cat mod_before.paf | fpa drop -l 50000 > mod_before.50k.paf

~/serverbin/miniasm/minidot mod_before.50k.paf > mod.before.50k.22002.eps
```

