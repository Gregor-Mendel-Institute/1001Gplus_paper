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
path.out = paste0(path.work, '01_data/02_alignment/pannagram_v10_4/intermediate/')
path.annotation = paste0(path.out, 'annotation/')
path.ann.own = paste0(path.annotation, 'own/')

path.features = paste0(path.annotation, 'features/')
if (!dir.exists(path.features)) dir.create(path.features)


# ---- Files ----

file.pan.merged = paste0(path.annotation, 'gff_pan_merged.gff')


# ***********************************************************************
# ---- mRNA ans Gene frequences ----

cl <- makeCluster(30)
registerDoParallel(cl)

df <- foreach(acc = accessions.true, .combine = rbind, .packages = c('pannagram', 'crayon', 'rhdf5', 'utils')) %dopar% {
  file.own.merged <- paste0(path.ann.own, 'gff_', acc, '.gff')
  
  gff.all <- read.table(file.own.merged, stringsAsFactors = FALSE)
  gff.all <- gff.all[gff.all$V3 %in% c('gene', 'mRNA'), ]
  gff.all$gr <- sapply(gff.all$V9, function(s) strsplit(s, '\\.')[[1]][1])
  gff.all$gr <- sapply(gff.all$gr, function(s) strsplit(s, ';')[[1]][1])
  gff.all$gr <- gsub('ID=', '', gff.all$gr)
  
  df.acc = data.frame(acc = acc, gr = gff.all$gr, type = gff.all$V3)
  return(df.acc)
}

stopCluster(cl)

cnt.genes = table(df$gr[df$type == 'gene'], df$acc[df$type == 'gene'])
cnt.mrnas = table(df$gr[df$type == 'mRNA'], df$acc[df$type == 'mRNA'])

save(list = ls(), file = "tmp_workspace_counts.RData")

write.table(cnt.genes, paste0(path.features, 'counts_gene.txt'), sep = '\t', quote = F)
write.table(cnt.mrnas, paste0(path.features, 'counts_mrna.txt'), sep = '\t', quote = F)




write.table(cnt.genes, paste0(path.features, 'counts_total.txt'), sep = '\t', quote = F)



# ***********************************************************************
# ---- Lyrata ----


# ***********************************************************************
# ---- mRNA on TEs ----


# ***********************************************************************
# ---- mRNA on mRNA ----


# ***********************************************************************
# ---- Genes on Genes ----


# ***********************************************************************
# ---- mRNA ans Gene frequences ----













