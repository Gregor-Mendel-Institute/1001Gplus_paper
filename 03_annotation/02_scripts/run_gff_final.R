
library(rhdf5)
library(crayon)
library(pannagram)

library(foreach)
library(doParallel)

numCores = 30

mergeAnn <- function(gff.gene, gff.rest){
  
  gff.all = rbind(gff.gene, gff.rest)
  
  gff.all$gr = gff.all$V9
  
  gff.all$gr = sapply(gff.all$gr, function(s) strsplit(s, '\\.')[[1]][1])
  gff.all$gr = sapply(gff.all$gr, function(s) strsplit(s, ';')[[1]][1])
  gff.all$gr = gsub('ID=', '',gff.all$gr)
  
  gr.pos = tapply(gff.all$V4, gff.all$gr, min)
  
  gff.all$pos = gr.pos[gff.all$gr]
  
  gff.all = gff.all[order(gff.all$pos),]
  gff.all = gff.all[order(gff.all$V1),]
  
  return(gff.all)
}


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

initial.vars <- ls()

# ***********************************************************************
# ---- Tair10 to pangen coordinates ----
pokaz('* Tair10 to pangen coordinates')

# Read the data
gff.tair = read.table(file.ann.tair, stringsAsFactors = F)
gff.tair = gff.tair[gff.tair$V3 == 'gene',]
pokaz('Number of genes', nrow(gff.tair))

# TAIR10 -> pangenome coordinate
file.ann.tair.pan = paste0(path.annotation, 'gff_pan_tair.rds')
if(!file.exists(file.ann.tair.pan)){
  pokaz('gff2gff...')
  gff2 = gff2gff(path.msa, acc1 = '0', acc2 = 'pangen', gff1 = gff.tair, n.chr = 5, exact.match = F, aln.type = aln.type,s.chr = 'Chr')
  saveRDS(gff2, file.ann.tair.pan)
} else {
  pokaz('Reading the pre-calculations...')
  gff2 = readRDS(file.ann.tair.pan)
}

# Feature extraction
pokaz('Feature extraction..')

pokaz('Number of genes', nrow(gff2))

gff2$chr = as.numeric(gsub('Chr', '', gff2$V1))
gff2$gr = sapply(gff2$V9, function(s) strsplit(s, ';')[[1]][1])

gff2$gr = gsub('ID=', '', gff2$gr)

pokaz('Length of annotations IDs trair10', max(nchar(gff2$gr)), min(nchar(gff2$gr)))

gff2$note = sapply(gff2$V9, function(s) strsplit(s, ';')[[1]][2])
gff2 = gff2[gff2$note == "Note=protein_coding_gene",]


# ***********************************************************************
# ---- Read pangenome annotation ----
pokaz('* Read pangenome annotation')

# Read the pangenome annotation
gff.pan.all = read.table(paste0(path.annotation,'gff_pan.gff'), stringsAsFactors = F)

# Get genes only

idx.mrna = gff.pan.all$V3 == 'mRNA'
gff.pan = gff.pan.all[idx.mrna,]
gff.pan$gr = sapply(gff.pan$V9, function(s) strsplit(s, '\\.')[[1]][1])
gff.pan$gr = gsub('ID=', '',gff.pan$gr)

gff.pan$V9 = paste0(gff.pan$V9, ';Parent=', gff.pan$gr)
gff.pan.all$V9[idx.mrna] = gff.pan$V9

pokaz('Length of the annotation IDs', max(nchar(gff.pan$gr)))

df.pangen = data.frame(beg = tapply(gff.pan$V4, gff.pan$gr, min),
                       end = tapply(gff.pan$V5, gff.pan$gr, max),
                       strand = tapply(gff.pan$V7, gff.pan$gr, unique))
df.pangen$gr = row.names(df.pangen)
df.pangen$chr = as.numeric(substr(df.pangen$gr, 3, 3))



# ***********************************************************************
# ---- Get pangenome coordinates of all genes ----

pokaz('* Get pangenome coordinates of all genes')

tair10.used = c()

df.pangen$tair = ''
df.pangen$tair.n = 0
for(i.chr in 1:5){
  for(i.s in c('-', '+')){
    pokaz(i.chr, i.s)
    idx = which((df.pangen$strand == i.s) & (df.pangen$chr == i.chr))
    df.tmp = df.pangen[idx,]

    gff2.tmp = gff2[(gff2$V7 == i.s) & (gff2$chr == i.chr),]
    pos = rep(0, max(max(df.tmp$end), max(gff2.tmp$V5)))

    for(irow in 1:nrow(gff2.tmp)){
      pos[gff2.tmp$V4[irow]:gff2.tmp$V5[irow]] = irow
    }

    for(irow in 1:nrow(df.tmp)){
      p = unique(pos[df.tmp$beg[irow]:df.tmp$end[irow]])
      p = setdiff(p, 0)
      if(length(p) == 0) next
      tair10.used = c(tair10.used, gff2.tmp$gr[p])
      tmp.names = paste0(gff2.tmp$gr[p], collapse = ',')
      df.tmp$tair[irow] = tmp.names
      df.tmp$tair.n[irow] = length(p)
    }

    df.pangen$tair[idx] = df.tmp$tair
    df.pangen$tair.n[idx] = df.tmp$tair.n

  }
}

df.pangen = df.pangen[order(df.pangen$beg),]
df.pangen = df.pangen[order(df.pangen$chr),]


df.pangen$tair.name = paste0(';Name=', df.pangen$tair)
df.pangen$tair.name[df.pangen$tair.n == 0] = ''
df.pangen$name = paste0('ID=', df.pangen$gr,df.pangen$tair.name)

tair10.used = unique(tair10.used)

gff2.add = gff2[!(gff2$gr %in% tair10.used),]

pokaz('Additional TAIR10 genes', nrow(gff2.add))

df.genes.new = data.frame(V1 = paste0('Pannagram_Chr', df.pangen$chr),
                          V2 = 'Pannagram',
                          V3 = 'gene',
                          V4 = df.pangen$beg,
                          V5 = df.pangen$end,
                          V6 = '.',
                          V7 = df.pangen$strand,
                          V8 = '.',
                          V9 = df.pangen$name)

df.genes.tair = data.frame(V1 = paste0('Pannagram_Chr', gff2.add$chr),
                           V2 = 'Pan_tair10',
                           V3 = 'gene',
                           V4 = gff2.add$V4,
                           V5 = gff2.add$V5,
                           V6 = '.',
                           V7 = gff2.add$V7,
                           V8 = '.',
                           V9 = paste0("ID=", gff2.add$gr, ";Name=", gff2.add$gr))

gff.new = rbind(df.genes.new, df.genes.tair)

gff.new = gff.new[order(gff.new$V4),]
gff.new = gff.new[order(gff.new$V1),]


# ***********************************************************************
# ---- Merge annotation ----
pokaz('* Merge annotation')

# save(list = ls(), file = "tmp_workspace_annotation.RData")
if(!file.exists(file.genes.only)){
  gff.pan.merged = mergeAnn(gff.new, gff.pan.all)

  options(scipen = 999)
  write.table(gff.pan.merged[,1:9], file = file.pan.merged, row.names = F, col.names = F, quote = F, sep = '\t')
  options(scipen = 0)

}

if(!file.exists(file.genes.only)){
  # Save genes only
  saveRDS(gff.new, file.genes.only)
}

# ***********************************************************************

# Cleanup variables
final.vars <- ls()
new.vars <- setdiff(final.vars, initial.vars)
pokaz('Veriables to remove', new.vars)
rm(list = new.vars)
gc()

# ***********************************************************************
# ---- Genes in own coordinates ----
pokaz('* Genes in own coordinates')

# Read own annotations
gff.own = read.table(paste0(path.annotation, 'gff_own.gff'), stringsAsFactors = F)
gff.new = readRDS(file.genes.only)

cl <- makeCluster(numCores)
registerDoParallel(cl)

tmp <- foreach(acc = accessions, .packages = c('pannagram', 'crayon', 'rhdf5')) %dopar% {
  pokaz(acc)

  file.own.tmp = paste0(path.ann.tmp, 'gff_', acc, '.rds')
  if(!file.exists(file.own.tmp)){
    gff.own.genes <- gff2gff(
      path.msa,
      acc1 = 'Pannagram',
      acc2 = acc,
      gff1 = gff.new,
      n.chr = 5,
      exact.match = FALSE,
      aln.type = aln.type,
      s.chr = '_Chr'
    )
    saveRDS(gff.own.genes, file.own.tmp)
  }
  return(NULL)
}

pokaz('Merge')

tmp <- foreach(acc = accessions, .packages = c('pannagram', 'crayon', 'rhdf5')) %dopar% {
  pokaz(acc)

  file.own.tmp = paste0(path.ann.tmp, 'gff_', acc, '.rds')

  if(acc == '22001_mod'){
    file.own.merged = paste0(path.ann.own, 'gff_', '220011', '.gff')
  } else {
    file.own.merged = paste0(path.ann.own, 'gff_', acc, '.gff')
  }


  if(file.exists(file.own.merged)) return(NULL)

  gff.own.genes = readRDS(file.own.tmp)

  gff.own.genes$V1 = gsub('22001_mod', '220011', gff.own.genes$V1)

  # Merge with mRNA and exons
  if(acc == '22001_mod'){
    gff.own.rest = gff.own[gff.own$V1 %in% paste0('220011', '_Chr', 1:5),]
  } else {
    gff.own.rest = gff.own[gff.own$V1 %in% paste0(acc, '_Chr', 1:5),]
  }
  gff.own.merged = mergeAnn(gff.own.genes, gff.own.rest)

  # Save
  options(scipen = 999)
  write.table(gff.own.merged[,1:9], file = file.own.merged, row.names = F, col.names = F, quote = F, sep = '\t')
  options(scipen = 0)

  return(NULL)
}
stopCluster(cl)

# ***********************************************************************
# ---- Confusing genes ----
pokaz('* Confusing genes')

file.gene.confusing = paste0(path.ann.tmp, 'idx_confusing.rds')

if(!file.exists(file.gene.confusing)){
  gene.confusing = c()
  for(acc in accessions){

    if(acc == '22001_mod'){
      file.own.merged = paste0(path.ann.own, 'gff_', '220011', '.gff')
    } else {
      file.own.merged = paste0(path.ann.own, 'gff_', acc, '.gff')
    }

    gff.all = read.table(file.own.merged, stringsAsFactors = F)

    gff.all$gr = sapply(gff.all$V9, function(s) strsplit(s, '\\.')[[1]][1])
    gff.all$gr = sapply(gff.all$gr, function(s) strsplit(s, ';')[[1]][1])
    gff.all$gr = gsub('ID=', '',gff.all$gr)


    cnt.strand = tapply(gff.all$V7, gff.all$gr, function(s) length(unique(s)))
    idx.confusing = which(cnt.strand > 1)
    gene.confusing = unique(c(gene.confusing, names(idx.confusing)))

    pokaz(acc, 'Number of confusing genes', length(idx.confusing), 'Total number of confusing', length(gene.confusing))
  }

  # Save
  saveRDS(gene.confusing, file.gene.confusing)
} else {
  gene.confusing = readRDS(file.gene.confusing)
}


# ***********************************************************************
# ---- Remove confusing genes from all of the files ----
pokaz('* Remove confusing genes from own')
for(acc in accessions){
  pokaz('Accession', acc)

  if(acc == '22001_mod'){
    file.own.merged = paste0(path.ann.own, 'gff_', '220011', '.gff')
  } else {
    file.own.merged = paste0(path.ann.own, 'gff_', acc, '.gff')
  }

  gff.all = read.table(file.own.merged, stringsAsFactors = F)
  gff.all$gr = sapply(gff.all$V9, function(s) strsplit(s, '\\.')[[1]][1])
  gff.all$gr = sapply(gff.all$gr, function(s) strsplit(s, ';')[[1]][1])
  gff.all$gr = gsub('ID=', '',gff.all$gr)

  pokaz('Before', nrow(gff.all))
  gff.all = gff.all[!(gff.all$gr %in% gene.confusing),]
  pokaz('After', nrow(gff.all))

  # Save
  options(scipen = 999)
  write.table(gff.all[,1:9], file = file.own.merged, row.names = F, col.names = F, quote = F, sep = '\t')
  options(scipen = 0)
}


# Remove from the pangenome annotation
pokaz('* Remove confusing genes from pangenome')

gff.all = read.table(file.pan.merged, stringsAsFactors = F)
gff.all$gr = sapply(gff.all$V9, function(s) strsplit(s, '\\.')[[1]][1])
gff.all$gr = sapply(gff.all$gr, function(s) strsplit(s, ';')[[1]][1])
gff.all$gr = gsub('ID=', '',gff.all$gr)

pokaz('Before', nrow(gff.all))
gff.all = gff.all[!(gff.all$gr %in% gene.confusing),]
pokaz('After', nrow(gff.all))

options(scipen = 999)
write.table(gff.all[,1:9], file = file.pan.merged, row.names = F, col.names = F, quote = F, sep = '\t')
options(scipen = 0)


# ***********************************************************************
# ---- Get mRNA and Gene sequences ----
pokaz('* Get sequences')

cl <- makeCluster(numCores)
registerDoParallel(cl)

tmp <- foreach(acc = accessions.true, .packages = c('pannagram', 'crayon', 'rhdf5')) %dopar% {
# for(acc in accessions.true) {

  file.fasta.genes = paste0(path.fasta, 'genes_',acc,'.fasta')
  file.fasta.mrnas = paste0(path.fasta, 'mrnas_',acc,'.fasta')

  write("", file = file.fasta.genes)
  write("", file = file.fasta.mrnas)

  pokaz(acc)

  file.own.merged = paste0(path.ann.own, 'gff_', acc, '.gff')

  gff.all = read.table(file.own.merged, stringsAsFactors = F)
  gff.all$gr = sapply(gff.all$V9, function(s) strsplit(s, '\\.')[[1]][1])
  gff.all$gr = sapply(gff.all$gr, function(s) strsplit(s, ';')[[1]][1])
  gff.all$gr = gsub('ID=', '',gff.all$gr)

  gff.all$name = paste(gff.all$gr, acc, gff.all$V3, gff.all$V4, gff.all$V5, gff.all$V7, gff.all$V5 -  gff.all$V4 + 1, sep = '|')

  for(i.chr in 1:5){
    if(acc == '220011'){
      s.chr = seq2nt(readFasta(paste0(path.chr, '22001_mod', '_chr', i.chr, '.fasta')))
    } else {
      s.chr = seq2nt(readFasta(paste0(path.chr, acc, '_chr', i.chr, '.fasta')))
    }

    for(s.type in c('gene', 'mRNA')){
      pokaz(i.chr, s.type)

      # save(list = ls(), file = "tmp_workspace_annotation.RData")

      gff.tmp = gff.all[(gff.all$V3 == s.type) & (gff.all$V1 == paste0(acc, '_Chr', i.chr)),]
      if(nrow(gff.tmp) == 0) next

      seqs = c()
      for(irow in 1:nrow(gff.tmp)){
        seqs = c(seqs, nt2seq(s.chr[gff.tmp$V4[irow]:gff.tmp$V5[irow]]))
      }
      names(seqs) = gff.tmp$name

      if(s.type == 'gene'){
        writeFasta(seqs, file.fasta.genes, append = T)
      } else {
        writeFasta(seqs, file.fasta.mrnas, append = T)
      }

    }
  }

}

stopCluster(cl)

# Combine to a common file:

file.fasta.genes = paste0(path.fasta, 'genes.fasta')
file.fasta.mrnas = paste0(path.fasta, 'mrnas.fasta')

write("", file = file.fasta.genes)
write("", file = file.fasta.mrnas)

for(acc in accessions) {
  pokaz('Accession', acc)

  file.fasta.genes.acc = paste0(path.fasta, 'genes_',acc,'.fasta')
  file.fasta.mrnas.acc = paste0(path.fasta, 'mrnas_',acc,'.fasta')

  if (file.exists(file.fasta.genes.acc)) {
    cat(readLines(file.fasta.genes.acc), file = file.fasta.genes, sep = "\n", append = TRUE)
  }
  if (file.exists(file.fasta.mrnas.acc)) {
    cat(readLines(file.fasta.mrnas.acc), file = file.fasta.mrnas, sep = "\n", append = TRUE)
  }

}

# stop('Done')

# ***********************************************************************
# ---- Get in Lyrata ----
pokaz('* Get in Lyrata')

file.lyrata.res = paste0(path.annotation, 'cov_in_lyrata.txt')
write("", file = file.lyrata.res)

# Combinations
s.comb.lyrata = c('1_1', '1_2', '2_3', '2_4', '3_3', '3_5', '4_6', '4_7', '5_6', '5_7', '5_8')
s.combs <- strsplit(s.comb.lyrata, "_")
s.combs <- data.frame(do.call(rbind, s.combs))
s.combs$s = s.comb.lyrata

# cl <- makeCluster(numCores)
# registerDoParallel(cl)

# results <- foreach(acc = accessions.true, .combine = rbind, .packages = c('pannagram', 'crayon', 'rhdf5')) %dopar% {
for(acc in accessions.true){
  pokaz('Accession', acc)

  file.own.merged = paste0(path.ann.own, 'gff_', acc, '.gff')
  gff.all = read.table(file.own.merged, stringsAsFactors = F)
  gff.all = gff.all[gff.all$V3 %in% c("mRNA", "gene"), ]
  gff.all$gr = sapply(gff.all$V9, function(s) strsplit(s, '\\.')[[1]][1])
  gff.all$gr = sapply(gff.all$gr, function(s) strsplit(s, ';')[[1]][1])
  gff.all$gr = gsub('ID=', '',gff.all$gr)
  gff.all.all = gff.all

  # accession_results <- list()
  
  for(i.chr in 1:5){
    pokaz('Chromosome', i.chr)

    gff.all = gff.all.all[grep(paste0(acc, '_Chr', i.chr), gff.all.all$V1),]

    lyrata.comb = s.combs$s[s.combs$X1 == i.chr]
    for(s.c in lyrata.comb){

      pokaz('Comb', s.c)

      file = paste0(path.lyrata, 'ref_', s.c, '_MN47.h5')
      v0 = h5read(file, paste(gr.accs.e, 'MN47', sep = ''))
      v0 = abs(v0)

      if(acc == '220011'){
        v1 = h5read(file, paste(gr.accs.e, '22001_mod', sep = ''))
      }else {
        v1 = h5read(file, paste(gr.accs.e, acc, sep = ''))
      }
      v1 = abs(v1)

      v = cbind(v1, v0)

      max.chr.len = max(nrow(v), max(abs(v[!is.na(v)])))

      v = v[v[,1]!=0,]
      v = v[!is.na(v[,1]),]
      v = v[!is.na(v[,2]),]
      idx.v.neg = which(v[,1] < 0)
      if(length(idx.v.neg) > 0){
        v[idx.v.neg,] = v[idx.v.neg,] * (-1)
      }

      # Get correspondence between two accessions
      v.corr = rep(0, max.chr.len)
      v.corr[v[,1]] = v[,2]

      # save(list = ls(), file = "tmp_workspace_lyrata.RData")
      lyrata.cov = c()
      for(irow in 1:nrow(gff.all)){
        value = mean(v.corr[gff.all$V4[irow]:gff.all$V5[irow]] > 0)
        lyrata.cov = c(lyrata.cov, value)
      }

      df.lyrata = data.frame(comb = s.c, acc = acc, gr = gff.all$gr, cov = lyrata.cov, type = gff.all$V3)
  
      write.table(df.lyrata, file.lyrata.res, append = T, row.names = F, col.names = F, sep = '\t', quote = F)
      # accession_results <- append(accession_results, list(df.lyrata))
    }
  }
  # return(accession_results)
}

# stopCluster(cl)

# write.table(results, file.lyrata.res, row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)














