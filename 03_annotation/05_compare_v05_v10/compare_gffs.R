library(foreach)
library(doParallel)


path.annot="/groups/nordborg/projects/the1001genomesplus/01_data/02_alignment/pannagram_v10_4/intermediate/annotation/"
path.prev="/groups/nordborg/projects/the1001genomesplus/01_data/04_annotation/02_pannagram/genes_v05/"

# Load data
gff.new.own <- read.table(paste0(path.annot, 'gff_own.gff'), stringsAsFactors = FALSE)
gff.new.pan <- read.table(paste0(path.annot, 'gff_pan.gff'), stringsAsFactors = FALSE)

# Process accessions
accessions <- unique(gff.new.own$V1)
accessions <- unique(sapply(accessions, function(s) strsplit(s, '_')[[1]][1]))
cat("Number of unique accessions:", length(accessions), "\n")

# Update gff.new.own$V9
for (i.chr in 1:5) {
  s.prev <- paste0('AT', i.chr, "Gr")
  s.new <- paste0('AT', i.chr, "SG")
  gff.new.own$V9 <- gsub(s.prev, s.new, gff.new.own$V9)
}

# Prepare parallel backend
num_cores <- 30
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Process accessions in parallel
gr.table.all <- foreach(acc = accessions, .combine = rbind, .packages = "base") %dopar% {
  cat('Processing accession:', acc, "\n")
  
  x.init <- read.table(paste0(path.prev, 'genes_v05_', acc, '.gff'))
  y.init <- gff.new.own[grep(paste0(acc, '_Chr'), gff.new.own$V1),]
  
  y <- y.init[y.init$V3 == 'mRNA',]
  x <- x.init[x.init$V3 == 'mRNA',]
  
  x$cover <- 0
  x$id <- 1:nrow(x)
  
  y$V10 <- ''
  y$id <- 1:nrow(y)
  
  # Nested loops for chromosomes and strands
  for (i.chr in 1:5) {
    for (s.strand in c('+', '-')) {
      cat('Processing chromosome:', i.chr, 'strand:', s.strand, "\n")
      
      xi <- x[x$V1 == paste0(acc, '_Chr', i.chr) & x$V7 == s.strand,]
      yi <- y[y$V1 == paste0(acc, '_Chr', i.chr) & y$V7 == s.strand,]
      
      pos.x <- integer(40000000)
      for (irow in seq_len(nrow(xi))) {
        pos.x[xi$V4[irow]:xi$V5[irow]] <- irow
      }
      
      for (irow in seq_len(nrow(yi))) {
        pp <- pos.x[yi$V4[irow]:yi$V5[irow]]
        m <- mean(pp > 0)
        if (m < 0.5) next
        
        p <- unique(pp[pp > 0])
        if (length(p) == 1) {
          y[yi$id[irow], 'V10'] <- xi$V9[p]
        } else if (length(p) > 1) {
          y[yi$id[irow], 'V10'] <- paste(xi$V9[p], collapse = '|')
        }
      }
      
      pos.y <- integer(40000000)
      for (irow in seq_len(nrow(yi))) {
        pos.y[yi$V4[irow]:yi$V5[irow]] <- irow
      }
      
      for (irow in seq_len(nrow(xi))) {
        m <- mean(pos.y[xi$V4[irow]:xi$V5[irow]])
        x$cover[xi$id[irow]] <- m
      }
    }
  }
  
  # Post-processing
  gr.table <- y[, c('V9', 'V10')]
  gr.table <- gr.table[!grepl("\\|", gr.table$V10) & gr.table$V10 != '',]
  for (icol in seq_len(ncol(gr.table))) {
    gr.table[, icol] <- sapply(gr.table[, icol], function(s) strsplit(s, ';')[[1]][1])
    gr.table[, icol] <- sapply(gr.table[, icol], function(s) strsplit(s, '\\.')[[1]][1])
    gr.table[, icol] <- sapply(gr.table[, icol], function(s) strsplit(s, '=')[[1]][2])
  }
  
  if (any(is.na(gr.table))) stop("NA detected")
  
  gr.table <- gr.table[!duplicated(gr.table$V10),]
  gr.table$acc <- acc
  gr.table
}

# Stop cluster
stopCluster(cl)

# Save or further process gr.table.all as needed
saveRDS(gr.table.all, paste0(path.annot,'gr_table_all.rds'))

