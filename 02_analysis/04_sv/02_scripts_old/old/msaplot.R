msaplot <- function(seqs.mx, 
                    msa.cols = c("A" = "#8ACD9D", "C" = "#EE7571", "G" = "#7E9CC8", "T" = "#FFD97C", '-'='#EEEDEF')){
  
  
  # seqs.mx = as.matrix(alignment)
  df <- reshape2::melt(seqs.mx)
  df$Var1 = factor(df$Var1, levels = rev(rownames(seqs.mx)))
  df$Var2 = as.numeric(df$Var2)
  
  g.msa = ggplot(df, aes(x = Var2, y = Var1, fill = value, color = value)) + 
    geom_tile() +
    scale_fill_manual(values = msa.cols) +
    scale_color_manual(values = msa.cols) +
    theme_bw() + 
    scale_x_continuous(limits = c(0, ncol(seqs.mx)+1), expand = c(0, 0)) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank()) + ylab('') + 
    xlab(NULL) +
    theme(legend.position = "none")
  
  return(g.msa)
  
}