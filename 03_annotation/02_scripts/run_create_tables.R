library(rhdf5)
library(crayon)
library(pannagram)
library(igraph)

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

path.simsearch = paste0(path.annotation, 'simsearch/')

path.fasta = paste0(path.annotation, 'fasta/')
path.features = paste0(path.annotation, 'features/')
if (!dir.exists(path.features)) dir.create(path.features)


# ---- Files ----

file.pan.merged = paste0(path.annotation, 'gff_pan_merged.gff')


# # ***********************************************************************
# # ---- mRNA and Gene frequencies ----
# pokaz('* mRNA and Gene frequencies')
# 
# cl <- makeCluster(30)
# registerDoParallel(cl)
# 
# df <- foreach(acc = accessions.true, .combine = rbind, .packages = c('pannagram', 'crayon', 'rhdf5', 'utils')) %dopar% {
#   file.own.merged <- paste0(path.ann.own, 'gff_', acc, '.gff')
#   
#   gff.all <- read.table(file.own.merged, stringsAsFactors = FALSE)
#   gff.all <- gff.all[gff.all$V3 %in% c('gene', 'mRNA'), ]
#   gff.all$gr <- sapply(gff.all$V9, function(s) strsplit(s, '\\.')[[1]][1])
#   gff.all$gr <- sapply(gff.all$gr, function(s) strsplit(s, ';')[[1]][1])
#   gff.all$gr <- gsub('ID=', '', gff.all$gr)
#   
#   df.acc = data.frame(acc = acc, gr = gff.all$gr, type = gff.all$V3)
#   return(df.acc)
# }
# 
# stopCluster(cl)
# 
# cnt.genes = table(df$gr[df$type == 'gene'], df$acc[df$type == 'gene'])
# cnt.mrnas = table(df$gr[df$type == 'mRNA'], df$acc[df$type == 'mRNA'])
# 
# # save(list = ls(), file = "tmp_workspace_counts.RData")
# 
# write.table(cnt.genes, paste0(path.features, 'counts_gene.txt'), sep = '\t', quote = F)
# write.table(cnt.mrnas, paste0(path.features, 'counts_mrna.txt'), sep = '\t', quote = F)
# 
# tot.genes = rowSums(cnt.genes)
# tot.mrnas = rowSums(cnt.mrnas)
# 
# gr.unique = unique(c(names(tot.genes), names(tot.mrnas)))
# 
# df.tot = data.frame(matrix(0, nrow = length(gr.unique), ncol = 2, dimnames = list(gr.unique, c('gene', 'mrna'))))
# df.tot[names(tot.genes),'gene'] = tot.genes
# df.tot[names(tot.mrnas),'mrna'] = tot.mrnas
# 
# write.table(df.tot, paste0(path.features, 'counts_total.txt'), sep = '\t', quote = F)
# 
# # ***********************************************************************
# # ---- Lyrata ----
# pokaz('* Lyrata')
# 
# df = read.table(paste0(path.annotation, 'cov_in_lyrata.txt'), stringsAsFactors = F)
# df = df[!is.na(df$V4),]
# 
# df$comb = paste(df$V3, df$V2, df$V5, sep  = '|')
# 
# x = tapply(df$V4, df$comb, max)
# 
# x.gene = tapply(df$V4[df$V5 == 'gene'], df$V3[df$V5 == 'gene'], max)
# x.mrna = tapply(df$V4[df$V5 == 'mRNA'], df$V3[df$V5 == 'mRNA'], max)
# 
# gr.unique = unique(c(names(x.gene), names(x.mrna)))
# 
# df.tot = data.frame(matrix(NA, nrow = length(gr.unique), ncol = 2, dimnames = list(gr.unique, c('gene', 'mrna'))))
# df.tot[names(x.gene),'gene'] = x.gene
# df.tot[names(x.mrna),'mrna'] = x.mrna
# 
# write.table(df.tot, paste0(path.features, 'lyrata_coverage.txt'), sep = '\t', quote = F)
# 
# # ***********************************************************************
# # ---- mRNA on TEs ----
# pokaz('* TEs')
# 
# cl <- makeCluster(30)
# registerDoParallel(cl)
# df <- foreach(acc = setdiff(accessions.true, '0'), .combine = rbind, .packages = c('pannagram', 'crayon', 'rhdf5', 'utils')) %dopar% {
#   
#   file.res = paste0(path.simsearch, 'out_temrnas_', acc, '/simsearch.tair10_tes.rds')
#   x = readRDS(file.res)
# 
#   x$cov = x$C1 / x$len1
#   
#   x$acc = sapply(x$V1, function(s) strsplit(s, '\\|')[[1]][1])
#   
#   y = tapply(x$cov, x$acc, max)
#   
#   df.acc = data.frame(gr = paste(names(y), acc, sep = '.'), group = names(y), acc = acc, cov = y)
#   return(df.acc)
# }
# stopCluster(cl)
# 
# # save(list = ls(), file = "tmp_workspace_counts.RData")
# 
# df = df[order(df$gr),]
# rownames(df) = NULL
# write.table(df, paste0(path.features, 'te_coverage_mrna.txt'), sep = '\t', quote = F, row.names = F)
# 
# df.stat = data.frame(mean = tapply(df$cov, df$group, mean), max = tapply(df$cov, df$group, max))
# df.stat$group = rownames(df.stat)
# df.stat$is.te = (df.stat$max >= 0.5) * 1
# 
# write.table(df.stat[,c(3,1,2, 4)], paste0(path.features, 'te_coverage_mrna_stat.txt'), sep = '\t', quote = F, row.names = F)
# 
# 
# ***********************************************************************
# ---- mRNA on mRNA ----
pokaz('* mRNA on mRNA')

cov.cutoff = 0.85

cl <- makeCluster(30)
registerDoParallel(cl)
df <- foreach(acc = setdiff(accessions.true, '0'), .combine = rbind, .packages = c('pannagram', 'crayon', 'rhdf5', 'utils')) %dopar% {

  file.res = paste0(path.simsearch, 'out_mrnas_', acc, '/simsearch.mrnas.rds')
  x = readRDS(file.res)

  x = x[(x$p1 > cov.cutoff) & (x$p8 > cov.cutoff),]

  return(x[, c('V1', 'V8')])
}
stopCluster(cl)
# save(list = ls(), file = "tmp_workspace_counts.RData")

# Create the graph and get connected components

x.graph <- igraph::make_graph(t(df), directed = F)
x.comp <- igraph::components(x.graph)

pokaz('Number of similarity groups', x.comp$no)

width <- nchar(as.character(x.comp$no))
simgr.names = paste0('SimGrM_', sprintf(paste0("%0", width+1, "d"), 1:x.comp$no))

df.sim = data.frame(group = names(x.comp$members), sim.gr = simgr.names[x.comp$members])

# Additional names
file.fasta.mrnas = paste0(path.fasta, 'mrnas.fasta')
x = readFasta(file.fasta.mrnas)
x.names = names(x)

x.add = setdiff(x.names, df.sim$group)
simgr.names.add = paste0('SimGrM_u', sprintf(paste0("%0", width, "d"), 1:length(x.add)))
df.add = data.frame(group = x.add, sim.gr = simgr.names.add)
pokaz('Number of additional groups', nrow(df.add))

df.sim = rbind(df.sim, df.add)
df.sim = df.sim[order(df.sim$group),]

write.table(df.sim, paste0(path.features, 'sim_groups_mrnas.txt'), sep = '\t', quote = F, row.names = F)

# ***********************************************************************
# ---- Genes on Genes ----
pokaz('* Genes on Genes')

cov.cutoff = 0.85

cl <- makeCluster(30)
registerDoParallel(cl)
df <- foreach(acc = setdiff(accessions.true, '0'), .combine = rbind, .packages = c('pannagram', 'crayon', 'rhdf5', 'utils')) %dopar% {
  
  file.res = paste0(path.simsearch, 'out_genes_', acc, '/simsearch.genes.rds')
  x = readRDS(file.res)
  
  x = x[(x$p1 > cov.cutoff) & (x$p8 > cov.cutoff),]
  
  return(x[, c('V1', 'V8')])
}
stopCluster(cl)

# Create the graph and get connected components

x.graph <- igraph::make_graph(t(df), directed = F)
x.comp <- igraph::components(x.graph)

pokaz('Number of similarity groups', x.comp$no)

width <- nchar(as.character(x.comp$no))
simgr.names = paste0('SimGrG_', sprintf(paste0("%0", width+1, "d"), 1:x.comp$no))

df.sim = data.frame(group = names(x.comp$members), sim.gr = simgr.names[x.comp$members])

# Additional names
file.fasta.mrnas = paste0(path.fasta, 'genes.fasta')
x = readFasta(file.fasta.mrnas)
x.names = names(x)

save(list = ls(), file = "tmp_workspace_counts.RData")

pokaz(setdiff(df.sim$group, x.names))

x.add = setdiff(x.names, df.sim$group)
simgr.names.add = paste0('SimGrG_u', sprintf(paste0("%0", width, "d"), 1:length(x.add)))
df.add = data.frame(group = x.add, sim.gr = simgr.names.add)
pokaz('Number of additional groups', nrow(df.add))

df.sim = rbind(df.sim, df.add)
df.sim = df.sim[order(df.sim$group),]

write.table(df.sim, paste0(path.features, 'sim_groups_genes.txt'), sep = '\t', quote = F, row.names = F)


# ***********************************************************************
# ---- Confusing genes and mRNAs ----
pokaz('* Confusing genes and mRNAs')

df.mrnas = read.table(paste0(path.features, 'sim_groups_mrnas.txt'), stringsAsFactors = F, header = 1)
df.genes = read.table(paste0(path.features, 'sim_groups_genes.txt'), stringsAsFactors = F, header = 1)


df.mrnas$gr = sapply(df.mrnas$group, function(s) strsplit(s, '\\|')[[1]][1])
df.genes$gr = sapply(df.genes$group, function(s) strsplit(s, '\\|')[[1]][1])

cnt.mrnas = tapply(df.mrnas$sim.gr, df.mrnas$gr, function(s) length(unique(s)))
cnt.genes = tapply(df.genes$sim.gr, df.genes$gr, function(s) length(unique(s)))


gr.unique = unique(c(names(cnt.genes), names(cnt.mrnas)))

df.tot = data.frame(matrix(NA, nrow = length(gr.unique), ncol = 2, dimnames = list(gr.unique, c('gene', 'mrna'))))
df.tot[names(cnt.genes),'gene'] = cnt.genes
df.tot[names(cnt.mrnas),'mrna'] = cnt.mrnas

df.tot[df.tot == 0] = NA
df.tot = df.tot[df.tot > 1] * 1

write.table(df.tot, paste0(path.features, 'sim_groups_confusion.txt'), sep = '\t', quote = F)







