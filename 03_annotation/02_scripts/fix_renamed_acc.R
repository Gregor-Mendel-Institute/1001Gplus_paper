# Script to fix IDs


library(rhdf5)
library(crayon)
library(pannagram)

library(foreach)
library(doParallel)

numCores = 30


# ***********************************************************************
# ---- Variables ----

aln.type = 'extra2_'
accessions = c('0', "10002", "10015", "10024", "1741", "22001_mod", "22002", "22003", "22004", 
               "22005", "22006", "22007", "6024", "6069", "6124", "6244", "6909", 
               "6966", "8236", "9075", "9537", "9543", "9638", "9728", "9764", 
               "9888", "9905", "9981")

accessions.true = c(setdiff(accessions, '22001_mod'), '220011')

gr.accs.e <- "accs/"

# ***********************************************************************

# ---- Paths ----

path.work = "../../../"

path.lyrata = paste0(path.work, "01_data/08_lyrata/pannagram_v10/intermediate/consensus/")

path.out = paste0(path.work, '01_data/02_alignment/pannagram_v10_4/intermediate/')
path.chr = paste0(path.out,'chromosomes/')
path.msa = paste0(path.out,'consensus/')
path.annotation = paste0(path.out, 'annotation/')
path.ann.tmp = paste0(path.annotation, 'tmp/')
path.ann.own = paste0(path.annotation, 'own/')
path.fasta = paste0(path.annotation, 'fasta/')

if (!dir.exists(path.ann.tmp)) dir.create(path.ann.tmp)
if (!dir.exists(path.ann.own)) dir.create(path.ann.own)
if (!dir.exists(path.fasta)) dir.create(path.fasta)

# ---- Files ----

file.ann.tair = paste0(path.work, "01_data/09_tair10/TAIR10_GFF3_genes.gff")
file.pan.merged = paste0(path.annotation, 'gff_pan_merged.gff')
file.genes.only = paste0(path.ann.tmp, 'gff_pan_genes.rds')



# ***********************************************************************
# Accession 220011


acc.sub = c('\\.220011' = '\\.22001_mod', 
            '\\.10015' = '\\.6962')

acc.sub.p = c('220011' = '22001_mod', 
            '10015' = '6962')


# Pan file to change
file.gff = paste0(path.annotation, 'gff_pan_merged.gff')
file.gff.fix = paste0(path.annotation, 'gff_pan_merged_renamed.gff')

gff = read.table(file.gff, stringsAsFactors = F)

for(i.acc in 1:length(accessions)){
  gff$V9 = gsub(names(acc.sub)[i.acc], acc.sub[i.acc], gff$V9)
}

options(scipen = 999)
cat("##gff-version 3\n", file = file.gff.fix)
write.table(gff, file = file.gff.fix, row.names = F, col.names = F, quote = F, sep = '\t', append = T)
options(scipen = 0)

# Own file to change


for(i.acc in 1:length(accessions)){
  acc = names(acc.sub.p)[i.acc]
  
  file.gff = paste0(path.ann.own, 'gff_', acc, '.gff')
  file.gff.fix = paste0(path.ann.own, 'gff_', acc.sub.p[i.acc], '_renamed.gff')
  
  # Read
  gff = read.table(file.gff, stringsAsFactors = F)
  
  # Substitute
  gff$V9 = gsub(names(acc.sub)[i.acc], acc.sub[i.acc], gff$V9)
  gff$V1 = gsub(acc, acc.sub.p[i.acc], gff$V1)
  
  # Save
  options(scipen = 999)
  cat("##gff-version 3\n", file = file.gff.fix)
  write.table(gff, file = file.gff.fix, row.names = F, col.names = F, quote = F, sep = '\t', append = T)
  options(scipen = 0)
}






