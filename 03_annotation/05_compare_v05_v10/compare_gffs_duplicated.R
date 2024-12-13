library(pannagram)
library(rhdf5)
library(ggplot2)
library(stringr)
library(crayon)
library(foreach)
library(doParallel)
library(gridExtra)


num.cores = 30
myCluster <- makeCluster(num.cores, type = "PSOCK") 
registerDoParallel(myCluster) 


gff.new.own = read.table('gff_own.gff', stringsAsFactors = F)
gff.new.pan = read.table('gff_pan.gff', stringsAsFactors = F)


gr.table.all = readRDS('gr_table_all.rds')
gr.table1 = unique(gr.table.all[,1:2])
gr.dup9 = gr.table1[duplicated(gr.table1[,2]),2]


gr.accs.e <- "accs/"
path.figures = 'figures/'

if (!dir.exists(path.figures)) dir.create(path.figures)

file.gr.dup = "gr_dup.rds"
if(file.exists(file.gr.dup)){
  pokaz("Reading...")
  gr.dup = readRDS(file.gr.dup)
} else {
  pokaz("Analysing...")
  gr.dup = data.frame(gr.prev = gr.dup9)
  gr.dup$gr.new = ''
  gr.dup$n.new = 0
  gr.dup$pos.beg = 0
  gr.dup$pos.end = 0
  gr.dup$len = 0
  
  gr.dup <- foreach(i.gr = 1:nrow(gr.dup), .combine = rbind, .packages = c('crayon', 'pannagram', 'rhdf5', 'stringr')) %dopar% {
    
    # Get the previous group for the current row
    g <- gr.dup$gr.prev[i.gr]
    
    # Extract and sort unique group identifiers
    gr2 <- sort(unique(gr.table1[gr.table1[,2] == g, 1]))
    
    # Transform group names from 'SG' to 'Gr'
    groups.target <- gsub('SG', 'Gr', gr2)
    
    # Create a new row for the current group with updated values
    gr.dup.row <- gr.dup[i.gr,]
    gr.dup.row$gr.new <- paste0(groups.target, collapse = '|') # Combine group targets into a string
    gr.dup.row$n.new <- length(groups.target) # Count the number of unique groups
    
    # Aggregate data for all groups in 'groups.target'
    tmp <- do.call(rbind, lapply(groups.target, function(g) {
      idx <- grep(g, gff.new.pan$V9) # Find matching rows in the dataset
      gff.new.pan[idx, ] # Extract corresponding rows
    }))
    
    # Define position shift (adjust as needed, default is 0 here)
    pos.shift <- 0
    pos.beg <- min(tmp$V4) - pos.shift # Compute beginning position
    pos.end <- max(tmp$V5) + pos.shift # Compute ending position
    
    # Update position and length information for the current group
    gr.dup.row$pos.beg <- pos.beg
    gr.dup.row$pos.end <- pos.end
    gr.dup.row$len <- pos.end - pos.beg + 1 # Calculate the length
    
    # Return the updated row
    gr.dup.row
  }
  
  gr.dup$chr = as.numeric(str_extract(gr.dup$gr.prev, "(?<=AT).*?(?=Gr)"))
  if(sum(is.na(gr.dup$chr)) > 0) stop("Some chromosome names are wrong")
  
  gr.dup$pattern = paste0(gr.dup$chr, '_', gr.dup$pos.beg, '_', gr.dup$pos.end)
  gr.dup$f.aln = paste0(path.figures, 'group_', gr.dup$pattern, '.fasta')
  gr.dup$f.msaplot = paste0(path.figures, 'msaplot_', gr.dup$pattern, '.png')
  gr.dup$f.dotplot = paste0(path.figures, 'dotplot_', gr.dup$pattern, '.png')
  
  saveRDS(gr.dup, file.gr.dup)
}

path.chr = '/groups/nordborg/projects/the1001genomesplus/01_data/02_alignment/pannagram_v10_4/intermediate/chromosomes/'
path.msa = '/groups/nordborg/projects/the1001genomesplus/01_data/02_alignment/pannagram_v10_4/intermediate/consensus/'

accessions.true = c('0', "10002", "10015", "10024", "1741", "220011", "22002", "22003", "22004", 
                    "22005", "22006", "22007", "6024", "6069", "6124", "6244", "6909", 
                    "6966", "8236", "9075", "9537", "9543", "9638", "9728", "9764", 
                    "9888", "9905", "9981")

pokaz('Number of groups', nrow(gr.dup))
gr.dup = gr.dup[,-c(1,2,3)]
gr.dup = unique(gr.dup)

pokaz('Number of groups', nrow(gr.dup))

for(acc in accessions.true){

  pokaz('Accession', acc)


  save(list = ls(), file = "tmp_workspace_cmp.RData")
  
  for(i.chr in 1:5){

    file.comb = paste0(path.msa, '/extra2_',i.chr, '_', i.chr, '.h5')

    if(acc == '220011'){
      v = h5read(file.comb, paste0(gr.accs.e, '22001_mod'))
      file.chr = paste0(path.chr, '22001_mod', '_chr', i.chr, '.fasta')
    } else {
      v = h5read(file.comb, paste0(gr.accs.e, acc))
      file.chr = paste0(path.chr, acc, '_chr', i.chr, '.fasta')
    }

    # Read the chromosome
    s.chr = seq2nt(readFasta(file.chr))

    idx.aln.chr = which(gr.dup$chr == i.chr)
    for(i.gr in idx.aln.chr){

      pos.beg = gr.dup$pos.beg[i.gr]
      pos.end = gr.dup$pos.end[i.gr]

      v.pos = v[pos.beg:pos.end]

      s = rep('-', length(v.pos))
      s[v.pos != 0] = s.chr[abs(v.pos[v.pos != 0])]
      if(sum(v.pos < 0) > 0){
        s[v.pos < 0] = justCompl(s[v.pos < 0])
      }

      s = nt2seq(s)
      names(s) = acc

      # save(list = ls(), file ="tmp_workspace.RData")
      # stop()

      if(file.exists(gr.dup$f.aln[i.gr])){
        writeFasta(s, gr.dup$f.aln[i.gr], append = T)
      } else {
        writeFasta(s, gr.dup$f.aln[i.gr])
      }

      pokaz(length(readFasta(gr.dup$f.aln[i.gr])))
    }

  }
}

pokaz('')

tmp <- foreach(i.gr = 1:nrow(gr.dup), .packages = c('crayon', 'pannagram', 'rhdf5', 'stringr', 'ggplot2', 'gridExtra')) %dopar% {
# for(i.gr in 1:nrow(gr.dup)){
  
  pos.beg = gr.dup$pos.beg[i.gr]
    pos.end = gr.dup$pos.end[i.gr]
  
  x = gff.new.pan[(gff.new.pan$V4 >=pos.beg) & 
                    (gff.new.pan$V5 <= pos.end) & 
                    (gff.new.pan$V3 == 'mRNA'),]
  x$chr = as.numeric(sapply(x$V1, function(s) strsplit(s, '_Chr')[[1]][2]))
  i.chr = 1
  x = x[x$chr == i.chr,]
  x = data.frame(beg = x$V4, end = x$V5)
  
  # x = unique(x)
  # orfplot(x)
  
  pokaz(gr.dup$f.aln[i.gr])
  aln = aln2mx(readFasta(gr.dup$f.aln[i.gr]))
  
  p = msaplot(aln) 
  # p
  
  p = p +  annotate("rect", xmin = x$beg - pos.beg, xmax = x$end - pos.beg,
                    # ymin = -Inf, ymax = Inf, alpha = 0.2,
                    ymin = 0, ymax = 0.5,  alpha = 0.8,
                    fill = 'blue')
  
  s = mx2cons(aln)
  
  # p.dot = dotplot(s, s, 15, 13)
  # p.dot = p.dot +  annotate("rect", xmin = x$beg - pos.beg, xmax = x$end - pos.beg,
  #                           # ymin = -Inf, ymax = Inf, alpha = 0.2,
  #                           ymin = 0, ymax = 0.5,  alpha = 0.8,
  #                           fill = 'blue')

  png(gr.dup$f.msaplot[i.gr], width = 8, height = 16, units = "in", res = 300) 
  # grid.arrange(p.dot, p, ncol = 1) 
  grid.arrange(p, ncol = 1) 
  dev.off()
  
  gc()
    
  # png(gr.dup$f.msaplot[i.gr], 
  #     width = 8, height = 6, units = "in", res = 300)
  # print(p)    
  # dev.off()
  # 
  # png(gr.dup$f.dotplot[i.gr], 
  #     width = 8, height = 8, units = "in", res = 300)
  # print(p)    
  # dev.off()
  
}

stopCluster(myCluster)
