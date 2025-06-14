setwd("/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS")
setwd("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/")
setwd("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/")

library(scales)
colorMap=c("Core" = "#D7191C",
           "SoftCore" = "#FDAE61",
           "MedFreq" = "#92C5DE",
           "LowFreq" = "#A6D96A",
           "Private" = "#1A9641")
#eracaps_125_RNAseq_samples.txt

colorMap_new=c("pc_known"="#5b9bd5",
               "tair10_highsim" = "#00b050",
               "tair10_medsim" = "#a6d96a",
               "Brassicaceae" = "#b8e2b0",
               "partial_or_no_homology" = "lightgrey",
               "TE_protein_similarity" = "#cc91df",
               "TE_high_nucl_similarity" = "#9966ff")

colorMap_new1=c("PC_TAIR10"="#5b9bd5",
               "PC_new" = "#ffc000",
               "TE_TAIR10" = "#a049c7",
                "TE_new" = "#cc91df")
"TE_new_nuc" = "#9966ff"





###################################
# define gene classes
###################################
#PC = not TE 
#TE = TE (>50% TE content) or TE-like (protein similarity to TE proteins)

filteroutgenes<-as.vector(gene_classes$group[gene_classes$nestgr_Remove_flag==T|gene_classes$AtGrp_overlap==T])
#2183

goodgenes<- as.vector(gene_classes$group[gene_classes$nestgr_Remove_flag==F&gene_classes$AtGrp_overlap==F]) #34153

  
  
pseudogenes_all<-as.vector(gene_classes$group[gene_classes$TEgene==F&gene_classes$is.TAIR10_Pseudogene==T])
pseudogenes<-pseudogenes_all[!(pseudogenes_all %in% filteroutgenes)]

pcgenes_all<-as.vector(gene_classes$group[gene_classes$TEgene==F & gene_classes$is.TAIR10_Pseudogene==F])
pcgenes<-pcgenes_all[!(pcgenes_all %in% filteroutgenes)]

tegenes_all<-as.vector(gene_classes$group[gene_classes$TEgene==T])
tegenes<-tegenes_all[!(tegenes_all %in% filteroutgenes)]

length(pseudogenes) #341
length(pcgenes) #28138
length(tegenes)# 5674

# total number of genes 
length(c(pcgenes_all,tegenes_all,pseudogenes_all)) #[1] 36336

# number of genes after filtering out "unanalyzable genes"
length(c(pcgenes,tegenes,pseudogenes)) #[1] 34153



newgenes_all<-as.vector(gene_classes$group[gene_classes$is.new==TRUE]) #8054
newgenes<-newgenes_all[!(newgenes_all %in% filteroutgenes)] #8054

tair10genes_all<-as.vector(gene_classes$group[gene_classes$is.new==FALSE]) #28282
tair10genes<-tair10genes_all[!(tair10genes_all %in% filteroutgenes)]  #26099  

############## Figure 6A categories 
#all 

genes_PC_TAIR10_all<-as.vector(gene_classes$group[gene_classes$group %in% pcgenes_all & gene_classes$is.new==F])
genes_PC_new_all<-as.vector(gene_classes$group[gene_classes$group %in% pcgenes_all & gene_classes$is.new==T])
genes_TE_TAIR10_all<-as.vector(gene_classes$group[gene_classes$group %in% tegenes_all & gene_classes$is.new==F])
genes_TE_new_all<-as.vector(gene_classes$group[gene_classes$group %in% tegenes_all & gene_classes$is.new==T])
genes_Pseudogene_all<-as.vector(gene_classes$group[gene_classes$group %in% pseudogenes_all ])



length(genes_PC_TAIR10_all) #26856
length(genes_PC_new_all) # 2661
length(genes_TE_TAIR10_all) #1286
length(genes_TE_new_all) #5109
length(genes_Pseudogene_all) #424


#filtered - what is used in the paper and plotted 
genes_PC_TAIR10<-as.vector(gene_classes$group[gene_classes$group %in% pcgenes & gene_classes$is.new==F])
genes_PC_new<-as.vector(gene_classes$group[gene_classes$group %in% pcgenes & gene_classes$is.new==T])
genes_TE_TAIR10<-as.vector(gene_classes$group[gene_classes$group %in% tegenes & gene_classes$is.new==F])
genes_TE_new<-as.vector(gene_classes$group[gene_classes$group %in% tegenes & gene_classes$is.new==T])
genes_Pseudogene<-as.vector(gene_classes$group[gene_classes$group %in% pseudogenes ])




# types of new genes 
new_tair_high<-as.vector(newgenes_classes$AtGrp[newgenes_classes$cat=="TAIR10_High_Similarity"])
new_tairP_medium<-as.vector(newgenes_classes$AtGrp[newgenes_classes$cat=="TAIR10_Medium_Similarity"])
new_brass<-as.vector(newgenes_classes$AtGrp[newgenes_classes$cat=="Brassicaceae"])
new_other<-as.vector(newgenes_classes$AtGrp[newgenes_classes$cat=="Other_Species"])
new_nohomol<-as.vector(newgenes_classes$AtGrp[newgenes_classes$cat=="Low or No homology"])


length(new_tair_high) #744
length(new_tairP_medium) #670
length(new_brass)#817
length(new_other) #9
length(new_nohomol)#421

#types of new TE genes 

telike<-as.vector(gene_classes$group[gene_classes$is.TE_like==T & !is.na(gene_classes$is.TE_like)& gene_classes$remove_final==F])
tenuc<-as.vector(gene_classes$group[gene_classes$Pannagram.is.te==1 & !is.na(gene_classes$Pannagram.is.te)& gene_classes$remove_final==F])
###############################################

# make lists of genes used in figures for Haim
###############################################
table_with_cat<-gene_classes[,1:2]
table_with_cat$genecat_old_new<-
  
write.table(genes_TE_new,)

write.table(genes_PC_new,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/genes_PC_new.txt",quote = F, sep="\t",col.names = F,row.names = F)

write.table(genes_TE_new,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/genes_TE_new.txt",quote = F, sep="\t",col.names = F,row.names = F)

write.table(genes_TE_TAIR10,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/genes_TE_TAIR10.txt",quote = F, sep="\t",col.names = F,row.names = F)

write.table(genes_PC_TAIR10,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/genes_PC_TAIR10.txt",quote = F, sep="\t",col.names = F,row.names = F)


write.table(pc.core,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/pc.core.txt",quote = F, sep="\t",col.names = F,row.names = F)
write.table(pc.highfreq,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/pc.highfreq.txt",quote = F, sep="\t",col.names = F,row.names = F)
write.table(pc.midfreq,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/pc.midfreq.txt",quote = F, sep="\t",col.names = F,row.names = F)
write.table(pc.lowfreq,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/pc.lowfreq.txt",quote = F, sep="\t",col.names = F,row.names = F)
write.table(pc.priv,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/pc.priv.txt",quote = F, sep="\t",col.names = F,row.names = F)


write.table(te.core,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/te.core.txt",quote = F, sep="\t",col.names = F,row.names = F)
write.table(te.highfreq,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/te.highfreq.txt",quote = F, sep="\t",col.names = F,row.names = F)
write.table(te.midfreq,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/te.midfreq.txt",quote = F, sep="\t",col.names = F,row.names = F)
write.table(te.lowfreq,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/te.lowfreq.txt",quote = F, sep="\t",col.names = F,row.names = F)
write.table(te.priv,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/te.priv.txt",quote = F, sep="\t",col.names = F,row.names = F)


write.table(genes_ancestral,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/genes_ancestral.txt",quote = F, sep="\t",col.names = F,row.names = F)
write.table(genes_lyrata_seqsim,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/genes_lyrata_seqsim.txt",quote = F, sep="\t",col.names = F,row.names = F)
write.table(genes_non_ancestral,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/genes_non_ancestral.txt",quote = F, sep="\t",col.names = F,row.names = F)


write.table(filteroutgenes,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/genelists/filteroutgenes.txt",quote = F, sep="\t",col.names = F,row.names = F)


##############################################

#categories gene numbers for pie chart Fig6A
###############################################
length(genes_PC_TAIR10) #25477
length(genes_PC_new) # 2661
length(genes_TE_TAIR10) #565
length(genes_TE_new) #5109
length(genes_Pseudogene) #341

#before revision previous annotation version :
#length(genes_PC_TAIR10) #25946
#length(genes_PC_new) # 2265
#length(genes_TE_TAIR10) #1126
#length(genes_TE_new) #5076
#length(genes_Pseudogene) #406)
# numbers are very similar to what we get now

###############################################

####################################################
# categorize PC genes by frequency 
####################################################


#frequency of locus (gene) presence : PC genes - TAIR10 and new (filtered)
#total N PC genes : 28,138
pc.core<-as.vector(gene_frequency$group[gene_frequency$freq_loc==27 &gene_frequency$group %in% pcgenes ])
pc.highfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=25 & gene_frequency$freq_loc<=26 & gene_frequency$group %in% pcgenes ])
pc.midfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=4 & gene_frequency$freq_loc<=24 & gene_frequency$group %in% pcgenes ])
pc.lowfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=2 & gene_frequency$freq_loc<=3 & gene_frequency$group %in% pcgenes ])
pc.priv<-as.vector(gene_frequency$group[gene_frequency$freq_loc==1  & gene_frequency$group %in% pcgenes ])
####################################################

#gene numbers for pie chart  Fig 6B
###############################################
length(pc.core)#24435
length(pc.highfreq)#601
length(pc.midfreq)#1413
length(pc.lowfreq)#572
length(pc.priv)#1084
###############################################

#check frequency of tair10 vs new genes
####################################################
pc_new.core<-as.vector(gene_frequency$group[gene_frequency$freq_loc==27 &gene_frequency$group %in% genes_PC_new ])
pc_new.highfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=25 & gene_frequency$freq_loc<=26 & gene_frequency$group %in% genes_PC_new ])
pc_new.midfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=4 & gene_frequency$freq_loc<=24 & gene_frequency$group %in% genes_PC_new ])
pc_new.lowfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=2 & gene_frequency$freq_loc<=3 & gene_frequency$group %in% genes_PC_new ])
pc_new.priv<-as.vector(gene_frequency$group[gene_frequency$freq_loc==1  & gene_frequency$group %in% genes_PC_new ])

pc_tair10.core<-as.vector(gene_frequency$group[gene_frequency$freq_loc==27 &gene_frequency$group %in% genes_PC_TAIR10 ])
pc_tair10.highfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=25 & gene_frequency$freq_loc<=26 & gene_frequency$group %in% genes_PC_TAIR10 ])
pc_tair10.midfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=4 & gene_frequency$freq_loc<=24 & gene_frequency$group %in% genes_PC_TAIR10 ])
pc_tair10.lowfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=2 & gene_frequency$freq_loc<=3 & gene_frequency$group %in% genes_PC_TAIR10 ])
pc_tair10.priv<-as.vector(gene_frequency$group[gene_frequency$freq_loc==1  & gene_frequency$group %in% genes_PC_TAIR10 ])


length(pc_new.core)#456
length(pc_new.highfreq)#76
length(pc_new.midfreq)#672
length(pc_new.lowfreq)#489
length(pc_new.priv)#968


length(pc_tair10.core)#23979
length(pc_tair10.highfreq)#525
length(pc_tair10.midfreq)#741
length(pc_tair10.lowfreq)#83
length(pc_tair10.priv)#116


TE_new.core<-as.vector(gene_frequency$group[gene_frequency$freq_loc==27 &gene_frequency$group %in% genes_TE_new ])
TE_new.highfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=25 & gene_frequency$freq_loc<=26 & gene_frequency$group %in% genes_TE_new ])
TE_new.midfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=4 & gene_frequency$freq_loc<=24 & gene_frequency$group %in% genes_TE_new ])
TE_new.lowfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=2 & gene_frequency$freq_loc<=3 & gene_frequency$group %in% genes_TE_new ])
TE_new.priv<-as.vector(gene_frequency$group[gene_frequency$freq_loc==1  & gene_frequency$group %in% genes_TE_new ])

TE_tair10.core<-as.vector(gene_frequency$group[gene_frequency$freq_loc==27 &gene_frequency$group %in% genes_TE_TAIR10 ])
TE_tair10.highfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=25 & gene_frequency$freq_loc<=26 & gene_frequency$group %in% genes_TE_TAIR10 ])
TE_tair10.midfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=4 & gene_frequency$freq_loc<=24 & gene_frequency$group %in% genes_TE_TAIR10 ])
TE_tair10.lowfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=2 & gene_frequency$freq_loc<=3 & gene_frequency$group %in% genes_TE_TAIR10 ])
TE_tair10.priv<-as.vector(gene_frequency$group[gene_frequency$freq_loc==1  & gene_frequency$group %in% genes_TE_TAIR10 ])


length(TE_new.core)#592
length(TE_new.highfreq)#394
length(TE_new.midfreq)#1743
length(TE_new.lowfreq)#926
length(TE_new.priv)#1454


length(TE_tair10.core)#295
length(TE_tair10.highfreq)#40
length(TE_tair10.midfreq)#126
length(TE_tair10.lowfreq)#31
length(TE_tair10.priv)#73

#all TEs
te.core<-as.vector(gene_frequency$group[gene_frequency$freq_loc==27 &gene_frequency$group %in% tegenes ])
te.highfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=25 & gene_frequency$freq_loc<=26 & gene_frequency$group %in% tegenes ])
te.midfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=4 & gene_frequency$freq_loc<=24 & gene_frequency$group %in% tegenes ])
te.lowfreq<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=2 & gene_frequency$freq_loc<=3 & gene_frequency$group %in% tegenes ])
te.priv<-as.vector(gene_frequency$group[gene_frequency$freq_loc==1  & gene_frequency$group %in% tegenes ])

length(te.core)#887
length(te.highfreq)#434
length(te.midfreq)#1869
length(te.lowfreq)#957
length(te.priv)#1527
###############################################################################

#combine lowfreq and private for silencing plots 
####################################################
pc.low_and_priv<-c(pc.lowfreq,pc.priv)
length(pc.low_and_priv)#1656

te.low_and_priv<-c(te.lowfreq,te.priv)
length(te.low_and_priv)#2484

####################################################

#Extended Data figure 20B - compare frequencies of our locus/mRNA to Mercier's paper frequencies
###############################################
# use all genes (no filtering as Mercier's paper also didn't filter)
#4 frequency categories - core, softcore, dispensable, private
##softcore - (25-26) = highfreq


#frequency of mRNA (gene model) presence : PC genes - TAIR10 and new (filtered)

pc_all.core<-as.vector(gene_frequency$group[gene_frequency$freq_loc==27 &gene_frequency$group %in% pcgenes_all & !is.na(gene_frequency$freq_loc)])
pc_all.softcore<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=25 & gene_frequency$freq_loc<=26 & gene_frequency$group %in% pcgenes_all& !is.na(gene_frequency$freq_loc) ])
pc_all.dispensable<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=2 & gene_frequency$freq_loc<=24 & gene_frequency$group %in% pcgenes_all & !is.na(gene_frequency$freq_loc)])
pc_all.priv<-as.vector(gene_frequency$group[gene_frequency$freq_loc==1  & gene_frequency$group %in% pcgenes_all & !is.na(gene_frequency$freq_loc) ])

te_all.core<-as.vector(gene_frequency$group[gene_frequency$freq_loc==27 &gene_frequency$group %in% tegenes_all & !is.na(gene_frequency$freq_loc)])
te_all.softcore<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=25 & gene_frequency$freq_loc<=26 & gene_frequency$group %in% tegenes_all& !is.na(gene_frequency$freq_loc) ])
te_all.dispensable<-as.vector(gene_frequency$group[gene_frequency$freq_loc>=2 & gene_frequency$freq_loc<=24 & gene_frequency$group %in% tegenes_all & !is.na(gene_frequency$freq_loc)])
te_all.priv<-as.vector(gene_frequency$group[gene_frequency$freq_loc==1  & gene_frequency$group %in% tegenes_all & !is.na(gene_frequency$freq_loc) ])

# define frequency categories mRNA-wise


#frequency of mRNA (gene model) presence : PC genes - TAIR10 and new (filtered)

pc.mRNA.core<-as.vector(gene_frequency$group[gene_frequency$freq_mrna==27 &gene_frequency$group %in% pcgenes_all & !is.na(gene_frequency$freq_mrna)])
pc.mRNA.softcore<-as.vector(gene_frequency$group[gene_frequency$freq_mrna>=25 & gene_frequency$freq_loc<=26 & gene_frequency$group %in% pcgenes_all& !is.na(gene_frequency$freq_mrna) ])
pc.mRNA.dispensable<-as.vector(gene_frequency$group[gene_frequency$freq_mrna>=2 & gene_frequency$freq_loc<=24 & gene_frequency$group %in% pcgenes_all & !is.na(gene_frequency$freq_mrna)])
pc.mRNA.priv<-as.vector(gene_frequency$group[gene_frequency$freq_mrna==1  & gene_frequency$group %in% pcgenes_all& !is.na(gene_frequency$freq_mrna) ])


te.mRNA.core<-as.vector(gene_frequency$group[gene_frequency$freq_mrna==27 &gene_frequency$group %in% tegenes_all & !is.na(gene_frequency$freq_mrna)])
te.mRNA.softcore<-as.vector(gene_frequency$group[gene_frequency$freq_mrna>=25 & gene_frequency$freq_loc<=26 & gene_frequency$group %in% tegenes_all& !is.na(gene_frequency$freq_mrna) ])
te.mRNA.dispensable<-as.vector(gene_frequency$group[gene_frequency$freq_mrna>=2 & gene_frequency$freq_loc<=24 & gene_frequency$group %in% tegenes_all & !is.na(gene_frequency$freq_mrna)])
te.mRNA.priv<-as.vector(gene_frequency$group[gene_frequency$freq_mrna==1  & gene_frequency$group %in% tegenes_all& !is.na(gene_frequency$freq_mrna) ])
##################################

#gene numbers for bar plot  Ext data Fig "Details of gene-ome analysis" B
###############################################
length(pc_all.core)#25479
length(pc_all.softcore)#698
length(pc_all.dispensable)#2214
length(pc_all.priv)#1091

length(pc.mRNA.core)#20536
length(pc.mRNA.softcore)#292
length(pc.mRNA.dispensable)#1608
length(pc.mRNA.priv)#1626

length(te_all.core)#1193
length(te_all.softcore)#538
length(te_all.dispensable)#3136
length(te_all.priv)#1528

length(te.mRNA.core)#281
length(te.mRNA.softcore)#47
length(te.mRNA.dispensable)#2192
length(te.mRNA.priv)#2721
###############################################

### define ancestral status 
#################### 
table(gene_classes$ancestral_cat)
##ancestral ancestral_onlyseq     non_ancestral 
#22357              3418             10561 
#Ancestral status 
genes_ancestral<-as.vector(gene_classes$group[gene_classes$ancestral_cat=="ancestral"])
genes_lyrata_seqsim<-as.vector(gene_classes$group[gene_classes$ancestral_cat=="ancestral_onlyseq"])
genes_non_ancestral<-as.vector(gene_classes$group[gene_classes$ancestral_cat=="non_ancestral"])

####################################################

# Figure 6D - ancestral status (by comparison with Lyrata)
####################################################
#barplots Fig 6D

length( intersect(genes_PC_TAIR10,genes_ancestral)  )#21077
length( intersect(genes_PC_TAIR10,genes_lyrata_seqsim)  )#1759
length( intersect(genes_PC_TAIR10,genes_non_ancestral)  )#2641

length( intersect(genes_PC_new,genes_ancestral)  )#291
length( intersect(genes_PC_new,genes_lyrata_seqsim)  )#622
length( intersect(genes_PC_new,genes_non_ancestral)  )#1748

length( intersect(genes_TE_new,genes_ancestral)  )#82
length( intersect(genes_TE_new,genes_lyrata_seqsim)  )#600
length( intersect(genes_TE_new,genes_non_ancestral)  )#4427

length( intersect(genes_TE_TAIR10,genes_ancestral)  )#170
length( intersect(genes_TE_TAIR10,genes_lyrata_seqsim)  )#81
length( intersect(genes_TE_TAIR10,genes_non_ancestral)  )#314

length( intersect(genes_Pseudogene,genes_ancestral)  )#98
length( intersect(genes_Pseudogene,genes_lyrata_seqsim)  )#70
length( intersect(genes_Pseudogene,genes_non_ancestral)  )#173


## Fig6E - ancestral status for PC genes of different frequency

length( intersect(pc.core,genes_ancestral)  )#20914
length( intersect(pc.core,genes_lyrata_seqsim)  )#1404
length( intersect(pc.core,genes_non_ancestral)  )#2117

length( intersect(pc.highfreq,genes_ancestral)  )#220
length( intersect(pc.highfreq,genes_lyrata_seqsim)  )#129
length( intersect(pc.highfreq,genes_non_ancestral)  )#252

length( intersect(pc.midfreq,genes_ancestral)  )#124
length( intersect(pc.midfreq,genes_lyrata_seqsim)  )#367
length( intersect(pc.midfreq,genes_non_ancestral)  )#922

length( intersect(pc.lowfreq,genes_ancestral)  )#22
length( intersect(pc.lowfreq,genes_lyrata_seqsim)  )#135
length( intersect(pc.lowfreq,genes_non_ancestral)  )#415

length( intersect(pc.priv,genes_ancestral)  )#88
length( intersect(pc.priv,genes_lyrata_seqsim)  )#346
length( intersect(pc.priv,genes_non_ancestral)  )#650
###############################

#### distribution of protein functions - barplots Fig 6-H 
################################################################

# gene enrichments 
table(gene_classes$uniprot_class)
#DNA integration 
#95 
#F-box 
#905 
#Leucine-rich repeat 
#392 
#Membrane / Transmembrane 
#5401 
#Plant defense 
#747 
#Proteases / Peptidase 
#810 
#Reverse transcriptase 
#1101 
#Signal peptide / Secreted 
#1814 
#TE / Transposon / Transposase / DNA transposition 
#1736 
#Transcription regulation 
#1906 
#Transferase 
#1838 
#Transport 
#1553 
#Zinc-finger / Zinc knuckle 
#1428 




Uniprot_table<-as.data.frame(table(gene_classes$uniprot_class,useNA = "ifany"))
Uniprot_table$uniprot_class<-Uniprot_table$Var1
Uniprot_table$all_genes<-Uniprot_table$Freq

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% filteroutgenes],useNA = "ifany"))
a$filterout<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[!(gene_classes$group %in% filteroutgenes)],useNA = "ifany"))
a$keptgenes<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)


a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% pcgenes],useNA = "ifany"))
a$pcgenes<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% tegenes],useNA = "ifany"))
a$tegenes<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% pc.core],useNA = "ifany"))
a$pc.core<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% pc.highfreq],useNA = "ifany"))
a$pc.highfreq<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% pc.midfreq],useNA = "ifany"))
a$pc.midfreq<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% pc.low_and_priv],useNA = "ifany"))
a$pc.low_and_priv<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% genes_ancestral & gene_classes$group %in% pcgenes],useNA = "ifany"))
a$pc.ancestral<-a$Freq 
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% genes_lyrata_seqsim & gene_classes$group %in% pcgenes],useNA = "ifany"))
a$pc.ancest_seqsim<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% genes_non_ancestral & gene_classes$group %in% pcgenes],useNA = "ifany"))
a$pc.non_ancestral<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)


a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% genes_ancestral & gene_classes$group %in% tegenes],useNA = "ifany"))
a$te.ancestral<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% genes_lyrata_seqsim & gene_classes$group %in% tegenes],useNA = "ifany"))
a$te.ancest_seqsim<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% genes_non_ancestral & gene_classes$group %in% tegenes],useNA = "ifany"))
a$te.non_ancestral<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% genes_ancestral & gene_classes$group %in% tegenes],useNA = "ifany"))
a$te.ancestral<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)


a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% genes_PC_new],useNA = "ifany"))
a$pc.new<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% genes_PC_TAIR10],useNA = "ifany"))
a$pc.tair10<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)


a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% genes_TE_new],useNA = "ifany"))
a$te.new<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

a<- as.data.frame(table(gene_classes$uniprot_class[gene_classes$group %in% genes_TE_TAIR10],useNA = "ifany"))
a$te.tair10<-a$Freq
Uniprot_table<-merge(Uniprot_table,a[,c(1,3)],by="Var1",all=T)

write.table(Uniprot_table,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/Uniprot_table.txt",quote = F, sep="\t",col.names = T,row.names = F)
##############################################


#######################################
#######################################
#Fig 6 and Supplement - epigenetic and expression properties  of genes
##########################################
############################################

#boxplots and lineplots for H3K9me2 6 accessions 
################################################
for (i in 1:6) {
  acc=as.character(accessions_chip$V1[i])
  
  a1<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.core & !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  a2<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.highfreq & !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  a3<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.midfreq & !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  a4<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.low_and_priv & !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/K9.",acc,".freq_cat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
   boxplot(a1,a2,a3,a4,notch = T, names=c("27","25-26","4-24","1-3"),xlab="frequency",las=2,main=paste("H3K9me2, acc",acc),outline = F, col=colorMap,cex.main=0.8)
  #add p values 
  ###########################
 try(  rm(d))
   try(  rm(w))
   if (length(a1)>10 & length(a2)>10){
  w<-wilcox.test(a1,a2)
  d<-w$p.value
   if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=2.5,cex=0.8)
   }
   rm(d)
   rm(w)
   if (length(a3)>10 & length(a2)>10){
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=2.5,cex=0.8)
   }
   rm(d)
   rm(w)
   if (length(a3)>10 & length(a4)>10){
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=2.5,cex=0.8)
   }
   ###########################
  dev.off()
  
  a1<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% te.core & !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  a2<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% te.highfreq & !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  a3<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% te.midfreq & !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  a4<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% te.low_and_priv & !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  #a5<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% private.loci.PC,paste("K9",acc,sep=".")]
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/K9.TEgenes.",acc,".freq_cat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  boxplot(a1,a2,a3,a4,notch = T, names=c("27","25-26","4-24","1-3"),xlab="frequency",las=2,main=paste("H3K9me2 TEgenes, acc",acc),outline = F, col=colorMap,cex.main=0.8)
  #add p values 
####################################
   rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
    w<-wilcox.test(a1,a2)
  d<-w$p.value
   if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=3,cex=0.8)
  }
  rm(d)
  rm(w)
  if (length(a3)>10 & length(a2)>10){
    w<-wilcox.test(a3,a2)
  d<-w$p.value
  
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=3.5,cex=0.8)
  }
  rm(d)
  rm(w)
  if (length(a3)>10 & length(a4)>10){
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=3,cex=0.8)
  }
   ###########################
  dev.off()
  
  b1<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% intersect(genes_ancestral, pcgenes),paste("K9",acc,sep=".")]
  b2<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% intersect(genes_lyrata_seqsim, pcgenes),paste("K9",acc,sep=".")]
  b3<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% intersect(genes_non_ancestral, pcgenes),paste("K9",acc,sep=".")]
  
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/K9.",acc,".ancestral_cat.pdf",sep=""),height =3,width = 2)
  par(mar=c(6,3,3,2),mgp=c(3,1,0)) 
   boxplot(b1, b2, b3,notch = T, las=2,main=paste("H3K9me2, acc",acc),outline = F, col=c('#3274B6','#9EC2E6','#F7CBAD'),names=c("full","seq","no"),xlab="ancestry status",cex.main=0.8)
  #add p values 
  ###########################
  w<-flexible_wilcox(b1,b2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=2,cex=0.8)
  w<-flexible_wilcox(b3,b2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=2,cex=0.8)
  
  
  ###########################
dev.off()
  
  #line plots for 12 categories by freq-ancestry
######################################################
  #fixed
  anc1<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_ancestral & chipseq.loci.quantstan$gene %in% pc.core& !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  anc2<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_lyrata_seqsim & chipseq.loci.quantstan$gene %in% pc.core& !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  anc3<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_non_ancestral & chipseq.loci.quantstan$gene %in% pc.core& !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  
  #softcore
  anc4<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_ancestral & chipseq.loci.quantstan$gene %in% pc.highfreq& !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  anc5<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_lyrata_seqsim & chipseq.loci.quantstan$gene %in% pc.highfreq& !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  anc6<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_non_ancestral & chipseq.loci.quantstan$gene %in% pc.highfreq& !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  
  #midfreq
  anc7<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_ancestral & chipseq.loci.quantstan$gene %in% pc.midfreq& !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  anc8<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_lyrata_seqsim & chipseq.loci.quantstan$gene %in% pc.midfreq& !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  anc9<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_non_ancestral & chipseq.loci.quantstan$gene %in% pc.midfreq& !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  
  #lowfreq and private
  anc10<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_ancestral & chipseq.loci.quantstan$gene %in% pc.low_and_priv & !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  anc11<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_lyrata_seqsim & chipseq.loci.quantstan$gene %in% pc.low_and_priv& !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  anc12<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_non_ancestral & chipseq.loci.quantstan$gene %in% pc.low_and_priv& !is.na(chipseq.loci.quantstan[,paste("K9",acc,sep=".")]),paste("K9",acc,sep=".")]
  
 
 m1<- median(anc1,na.rm = T)
 m2<- median(anc2,na.rm = T)
 m3<- median(anc3,na.rm = T)
 m4<- median(anc4,na.rm = T)
 m5<- median(anc5,na.rm = T)
 m6<- median(anc6,na.rm = T)
 m7<- median(anc7,na.rm = T)
 m8<- median(anc8,na.rm = T)
 m9<- median(anc9,na.rm = T)
 m10<- median(anc10,na.rm = T)
 m11<- median(anc11,na.rm = T)
 m12<- median(anc12,na.rm = T)
 
 values <- na.omit(c(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12))
 
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/K9.",acc,".lineplot_medians_anc_and_freq.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  
  plot(1:4,c(m1,m4,m7,m10),pch=19, type="o", col="#3274B6", lwd=1.2,xaxt='n', xlab="frequency", ylab="", main=paste("H3K9me2, acc",acc),ylim=c(min(values)-0.05,max(values)+0.05),cex.main=0.8,xlim=c(0.5,4.5))
  axis(1,at=1:4,labels=c("27","25-26","4-24","1-3"),las=2)
  points(1:4,c(m2,m5,m7,m11), pch=19,type="o",lwd=1.2, col="#9EC2E6")
  points(1:4,c(m3,m6,m9,m12), pch=19,type="o", lwd=1.2,col='#F7CBAD')

  #add p values
###################################
  a1=anc1
  a2=anc2
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=1,y=posit,cex=0.8)}
  }
  a1=anc2
  a2=anc3
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=1,y=posit,cex=0.8)}
  }
  a1=anc4
  a2=anc5
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=2,y=posit,cex=0.8)}
}
  a1=anc5
  a2=anc6
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=2,y=posit,cex=0.8)}
  }
  a1=anc7
  a2=anc8
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=3,y=posit,cex=0.8)}
  }
  
  a1=anc8
  a2=anc9
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=3,y=posit,cex=0.8)}
}
  a1=anc10
  a2=anc11
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=4,y=posit,cex=0.8)}
}
  a1=anc11
  a2=anc12
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=4,y=posit,cex=0.8)}
  }
  
  ###########################
  dev.off()
  
  
}
################################################

#boxplots and lineplots for H3K27me3 5 accessions 
################################################
for (i in 1:5) {
  acc=c("6909","6024","9905","9888",'1741')[i]
  
  a1<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.core,paste("K27",acc,sep=".")]
  a2<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.highfreq,paste("K27",acc,sep=".")]
  a3<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.midfreq,paste("K27",acc,sep=".")]
  a4<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.low_and_priv,paste("K27",acc,sep=".")]
  #a5<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% private.loci.PC,paste("K27",acc,sep=".")]
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/K27.",acc,".freq_cat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  boxplot(a1,a2,a3,a4,notch = T, names=c("27","25-26","4-24","1-3"),xlab="frequency",las=2,main=paste("H3K27me3, acc",acc),outline = F, col=colorMap,cex.main=0.8)
  #add p values 
  ###########################
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=2.5,cex=0.8)
  w<-flexible_wilcox(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=2.5,cex=0.8)
  w<-flexible_wilcox(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=2.5,cex=0.8)
 
  ###########################
  dev.off()
  
  
  a1<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% te.core,paste("K27",acc,sep=".")]
  a2<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% te.highfreq,paste("K27",acc,sep=".")]
  a3<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% te.midfreq,paste("K27",acc,sep=".")]
  a4<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% te.low_and_priv,paste("K27",acc,sep=".")]
  #a5<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% private.loci.PC,paste("K27",acc,sep=".")]
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/K27.TEgenes.",acc,".freq_cat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  boxplot(a1,a2,a3,a4,notch = T, names=c("27","25-26","4-24","1-3"),xlab="frequency",las=2,main=paste("H3K27me3, TEgenes, acc",acc),outline = F, col=colorMap,cex.main=0.8)
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=1.5,cex=0.8)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=1.5,cex=0.8)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=1.5,cex=0.8)
   ###########################
  dev.off()
  
  b1<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% intersect(genes_ancestral, pcgenes),paste("K27",acc,sep=".")]
  b2<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% intersect(genes_lyrata_seqsim, pcgenes),paste("K27",acc,sep=".")]
  b3<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% intersect(genes_non_ancestral, pcgenes),paste("K27",acc,sep=".")]
  
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/K27.",acc,".ancestral_cat.pdf",sep=""),height =3,width = 2)
  par(mar=c(6,3,3,2),mgp=c(3,1,0)) 
  boxplot(b1, b2, b3,notch = T, las=2,main=paste("H3K27me3, acc",acc),outline = F, col=c('#3274B6','#9EC2E6','#F7CBAD'),names=c("full","seq","no"),xlab="ancestry status",cex.main=0.8)
  #add p values 
  ###########################
  w<-flexible_wilcox(b1,b2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=2,cex=0.8)
  w<-flexible_wilcox(b3,b2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=2,cex=0.8)
  ###########################
  dev.off()
  
  
  #fixed
  anc1<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_ancestral & chipseq.loci.quantstan$gene %in% pc.core & !is.na(chipseq.loci.quantstan[,paste("K27",acc,sep=".")]),paste("K27",acc,sep=".")]
  anc2<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_lyrata_seqsim & chipseq.loci.quantstan$gene %in% pc.core & !is.na(chipseq.loci.quantstan[,paste("K27",acc,sep=".")]),paste("K27",acc,sep=".")]
  anc3<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_non_ancestral & chipseq.loci.quantstan$gene %in% pc.core & !is.na(chipseq.loci.quantstan[,paste("K27",acc,sep=".")]),paste("K27",acc,sep=".")]
  
  #softcore
  anc4<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_ancestral & chipseq.loci.quantstan$gene %in% pc.highfreq & !is.na(chipseq.loci.quantstan[,paste("K27",acc,sep=".")]),paste("K27",acc,sep=".")]
  anc5<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_lyrata_seqsim & chipseq.loci.quantstan$gene %in% pc.highfreq & !is.na(chipseq.loci.quantstan[,paste("K27",acc,sep=".")]),paste("K27",acc,sep=".")]
  anc6<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_non_ancestral & chipseq.loci.quantstan$gene %in% pc.highfreq & !is.na(chipseq.loci.quantstan[,paste("K27",acc,sep=".")]),paste("K27",acc,sep=".")]
  
  #midfreq
  anc7<-na.omit(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_ancestral & chipseq.loci.quantstan$gene %in% pc.midfreq & !is.na(chipseq.loci.quantstan[,paste("K27",acc,sep=".")]),paste("K27",acc,sep=".")])
  anc8<-na.omit(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_lyrata_seqsim & chipseq.loci.quantstan$gene %in% pc.midfreq & !is.na(chipseq.loci.quantstan[,paste("K27",acc,sep=".")]),paste("K27",acc,sep=".")])
  anc9<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_non_ancestral & chipseq.loci.quantstan$gene %in% pc.midfreq & !is.na(chipseq.loci.quantstan[,paste("K27",acc,sep=".")]),paste("K27",acc,sep=".")]
  
  #lowfreq and private
  anc10<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_ancestral & chipseq.loci.quantstan$gene %in% pc.low_and_priv & !is.na(chipseq.loci.quantstan[,paste("K27",acc,sep=".")]),paste("K27",acc,sep=".")]
  anc11<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_lyrata_seqsim & chipseq.loci.quantstan$gene %in% pc.low_and_priv & !is.na(chipseq.loci.quantstan[,paste("K27",acc,sep=".")]),paste("K27",acc,sep=".")]
  anc12<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_non_ancestral & chipseq.loci.quantstan$gene %in% pc.low_and_priv & !is.na(chipseq.loci.quantstan[,paste("K27",acc,sep=".")]),paste("K27",acc,sep=".")]
  
  
  m1<- median(anc1,na.rm = T)
  m2<- median(anc2,na.rm = T)
  m3<- median(anc3,na.rm = T)
  m4<- median(anc4,na.rm = T)
  m5<- median(anc5,na.rm = T)
  m6<- median(anc6,na.rm = T)
  m7<- median(anc7,na.rm = T)
  m8<- median(anc8,na.rm = T)
  m9<- median(anc9,na.rm = T)
  m10<- median(anc10,na.rm = T)
  m11<- median(anc11,na.rm = T)
  m12<- median(anc12,na.rm = T)
  
  
  values <- na.omit(c(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12))
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/K27.",acc,".lineplot_medians_anc_and_freq.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  plot(1:4,c(m1,m4,m7,m10),pch=19, type="o", col="#3274B6", lwd=1.2,xaxt='n', xlab="frequency", ylab="", main=paste("H3K27me3, acc",acc),ylim=c(min(values)-0.05,max(values)+0.05),cex.main=0.8,xlim=c(0.5,4.5))
  axis(1,at=1:4,labels=c("27","25-26","4-24","1-3"),las=2)
  points(1:4,c(m2,m5,m7,m11), pch=19,type="o",lwd=1.2, col="#9EC2E6")
  points(1:4,c(m3,m6,m9,m12), pch=19,type="o", lwd=1.2,col='#F7CBAD')
  
  
  #add p values 
  ###########################
  a1=anc1
  a2=anc2
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=1,y=posit,cex=0.8)}
  }
  a1=anc2
  a2=anc3
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=1,y=posit,cex=0.8)}
  }
  a1=anc4
  a2=anc5
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=2,y=posit,cex=0.8)}
  }
  a1=anc5
  a2=anc6
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=2,y=posit,cex=0.8)}
  }
  a1=anc7
  a2=anc8
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=3,y=posit,cex=0.8)}
  }
  a1=anc8
  a2=anc9
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=3,y=posit,cex=0.8)}
}
  a1=anc10
  a2=anc11
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=4,y=posit,cex=0.8)}
  }
  
  
  a1=anc11
  a2=anc12
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
    w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(median(a1,na.rm = T),median(a2,na.rm = T))-(max(median(a1,na.rm = T),median(a2,na.rm = T))-min(median(a1,na.rm = T),median(a2,na.rm = T)))/2
  if (d<0.01) {text(b,x=4,y=posit,cex=0.8)}
  }
  
  ###########################
  dev.off()
    }
################################################

install.packages("vioplot")
library(vioplot)

#boxplots and lineplots for CG,CHG, CHH in 12 accessions 
################################################
for (i in 2:13) {
  acc=as.character(accessions_methylation$V1[i-1])
 
  
  a1<-CG.loci[CG.loci$gene %in% pc.core,i]
  a2<-CG.loci[CG.loci$gene %in% pc.highfreq,i]
  a3<-CG.loci[CG.loci$gene %in% pc.midfreq,i]
  a4<-CG.loci[CG.loci$gene %in% pc.low_and_priv,i]
  #a5<-CHH.loci[CHH.loci$gene %in% private.loci.PC,i]
  
  
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CG.",acc,".freq_cat.pdf",sep=""),height =3,width = 2.5)
par(mar=c(6,3,3,2),mgp=c(2,1,0)) 
vioplot(a1,a2,a3,a4,col=colorMap, names=c("27","25-26","4-24","1-3"),ylab="methylation level",xlab="frequency",las=2,main=paste("CG, leaves acc",acc) ,cex.main=0.8)

#add p values 
###########################
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
w<-flexible_wilcox(a3,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
w<-flexible_wilcox(a3,a4)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)

###########################
dev.off()

a1<-CG.loci[CG.loci$gene %in% te.core,i]
a2<-CG.loci[CG.loci$gene %in% te.highfreq,i]
a3<-CG.loci[CG.loci$gene %in% te.midfreq,i]
a4<-CG.loci[CG.loci$gene %in% te.low_and_priv,i]
#a5<-CHH.loci[CHH.loci$gene %in% private.loci.PC,i]


pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CG.TEgenes.",acc,".freq_cat.pdf",sep=""),height =3,width = 2.5)
par(mar=c(6,3,3,2),mgp=c(6,1,0)) 
library(vioplot)
vioplot(a1,a2,a3,a4,col=colorMap, names=c("27","25-26","4-24","1-3"),ylab="methylation level",xlab="frequency",las=2,main=paste("CG TEgenes, leaves acc",acc) ,cex.main=0.8)

#add p values 
###########################
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
w<-flexible_wilcox(a3,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
w<-flexible_wilcox(a3,a4)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=1)

###########################
dev.off()



a1<-CG.loci[CG.loci$gene %in% intersect(genes_ancestral,pcgenes),i]
a2<-CG.loci[CG.loci$gene %in% intersect(genes_lyrata_seqsim,pcgenes),i]
a3<-CG.loci[CG.loci$gene %in% intersect(genes_non_ancestral,pcgenes),i]

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CG.",acc,".ancest_cat.pdf",sep=""),height =3,width = 2.5)
par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
vioplot(a1,a2,a3,ylab="methylation level",las=2,main=paste("CG, leaves acc",acc) ,col=c('#3274B6','#9EC2E6','#F7CBAD'),names=c("full","seq","no"),xlab="ancestry status",cex.main=0.8)

#add p values 
###########################
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=1)
w<-flexible_wilcox(a3,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=1)
###########################
dev.off()



a1<-CHG.loci[CHG.loci$gene %in% pc.core,i]
a2<-CHG.loci[CHG.loci$gene %in% pc.highfreq,i]
a3<-CHG.loci[CHG.loci$gene %in% pc.midfreq,i]
a4<-CHG.loci[CHG.loci$gene %in% pc.low_and_priv,i]
#a5<-CHH.loci[CHH.loci$gene %in% private.loci.PC,i]

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CHG.",acc,".freq_cat.pdf",sep=""),height =3,width = 2.5)
par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
boxplot(a1,a2,a3,a4,col=colorMap, names=c("27","25-26","4-24","1-3"),ylab="methylation level",xlab="frequency",las=2,main=paste("CHG, leaves acc",acc) ,cex.main=0.8,notch = T,outline = F)
#add p values 
###########################
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.38,cex=0.8)
w<-flexible_wilcox(a3,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.38,cex=0.8)
w<-flexible_wilcox(a3,a4)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.38,cex=0.8)
###########################

dev.off()


a1<-CHH.loci[CHH.loci$gene %in% pc.core,i]
a2<-CHH.loci[CHH.loci$gene %in% pc.highfreq,i]
a3<-CHH.loci[CHH.loci$gene %in% pc.midfreq,i]
a4<-CHH.loci[CHH.loci$gene %in% pc.low_and_priv,i]
#a5<-CHH.loci[CHH.loci$gene %in% private.loci.PC,i]

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CHH.",acc,".freq_cat.pdf",sep=""),height =3,width = 2.5)
par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
boxplot(a1,a2,a3,a4,col=colorMap, names=c("27","25-26","4-24","1-3"),ylab="methylation level",xlab="frequency",las=2,main=paste("CHH, leaves acc",acc) ,cex.main=0.8,notch = T,outline = F)
#add p values 
###########################
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.08,cex=0.8)
w<-flexible_wilcox(a3,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.08,cex=0.8)
w<-flexible_wilcox(a3,a4)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.08,cex=0.8)
###########################

dev.off()




a1<-CHH.loci[CHH.loci$gene %in% te.core,i]
a2<-CHH.loci[CHH.loci$gene %in% te.highfreq,i]
a3<-CHH.loci[CHH.loci$gene %in% te.midfreq,i]
a4<-CHH.loci[CHH.loci$gene %in% te.low_and_priv,i]
#a5<-CHH.loci[CHH.loci$gene %in% private.loci.PC,i]

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CHH.TEgenes",acc,".freq_cat.pdf",sep=""),height =3,width = 2.5)
par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
library(vioplot)
boxplot(a1,a2,a3,a4,col=colorMap, names=c("27","25-26","4-24","1-3"),ylab="methylation level",xlab="frequency",las=2,main=paste("CHH,TEgenes,leaves acc",acc) ,cex.main=0.8,notch = T,outline = F)
#add p values 
###########################
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.18,cex=0.8)
w<-flexible_wilcox(a3,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.18,cex=0.8)
w<-flexible_wilcox(a3,a4)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=3.5,y=0.18,cex=0.8)
###########################

dev.off()


a1<-CHH.loci[CHH.loci$gene %in% intersect(genes_ancestral,pcgenes),i]
a2<-CHH.loci[CHH.loci$gene %in% intersect(genes_lyrata_seqsim,pcgenes),i]
a3<-CHH.loci[CHH.loci$gene %in% intersect(genes_non_ancestral,pcgenes),i]

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CHH.",acc,".ancest_cat.pdf",sep=""),height =3,width = 2.5)
par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
library(vioplot)
boxplot(a1,a2,a3,ylab="methylation level",las=2,main=paste("CHH, leaves acc",acc) ,col=c('#3274B6','#9EC2E6','#F7CBAD'),names=c("full","seq","no"),xlab="ancestry status",cex.main=0.8,outline = F,notch = T)

#add p values 
###########################
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.01)
w<-flexible_wilcox(a3,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.01)
###########################
dev.off()



a1<-CHG.loci[CHG.loci$gene %in% intersect(genes_ancestral,pcgenes),i]
a2<-CHG.loci[CHG.loci$gene %in% intersect(genes_lyrata_seqsim,pcgenes),i]
a3<-CHG.loci[CHG.loci$gene %in% intersect(genes_non_ancestral,pcgenes),i]

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CHG.",acc,".ancest_cat.pdf",sep=""),height =3,width = 2.5)
par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
library(vioplot)
boxplot(a1,a2,a3,ylab="methylation level",las=2,main=paste("CHG, leaves acc",acc) ,col=c('#3274B6','#9EC2E6','#F7CBAD'),names=c("full","seq","no"),xlab="ancestry status",cex.main=0.8,outline = F,notch = T)

#add p values 
###########################
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=1.5,y=0.1)
w<-flexible_wilcox(a3,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=0.1)
###########################
dev.off()




#line plots for 12 categories by freq-ancestry

#fixed
anc1<-na.omit(CHH.loci[CHH.loci$gene %in% genes_ancestral & CHH.loci$gene %in% pc.core,acc])
anc2<-na.omit(CHH.loci[CHH.loci$gene %in% genes_lyrata_seqsim & CHH.loci$gene %in% pc.core,acc])
anc3<-na.omit(CHH.loci[CHH.loci$gene %in% genes_non_ancestral & CHH.loci$gene %in% pc.core,acc])

#softcore
anc4<-na.omit(CHH.loci[CHH.loci$gene %in% genes_ancestral & CHH.loci$gene %in% pc.highfreq,acc])
anc5<-na.omit(CHH.loci[CHH.loci$gene %in% genes_lyrata_seqsim & CHH.loci$gene %in% pc.highfreq,acc])
anc6<-na.omit(CHH.loci[CHH.loci$gene %in% genes_non_ancestral & CHH.loci$gene %in% pc.highfreq,acc])

#midfreq
anc7<-na.omit(CHH.loci[CHH.loci$gene %in% genes_ancestral & CHH.loci$gene %in% pc.midfreq,acc])
anc8<-na.omit(CHH.loci[CHH.loci$gene %in% genes_lyrata_seqsim & CHH.loci$gene %in% pc.midfreq,acc])
anc9<-na.omit(CHH.loci[CHH.loci$gene %in% genes_non_ancestral & CHH.loci$gene %in% pc.midfreq,acc])

#lowfreq and private
anc10<-na.omit(CHH.loci[CHH.loci$gene %in% genes_ancestral & CHH.loci$gene %in% pc.low_and_priv,acc])
anc11<-na.omit(CHH.loci[CHH.loci$gene %in% genes_lyrata_seqsim & CHH.loci$gene %in% pc.low_and_priv,acc])
anc12<-na.omit(CHH.loci[CHH.loci$gene %in% genes_non_ancestral & CHH.loci$gene %in% pc.low_and_priv,acc])

#median too low, so

m1<- mean(anc1,na.rm = T)
m2<- mean(anc2,na.rm = T)
m3<- mean(anc3,na.rm = T)
m4<- mean(anc4,na.rm = T)
m5<- mean(anc5,na.rm = T)
m6<- mean(anc6,na.rm = T)
m7<- mean(anc7,na.rm = T)
m8<- mean(anc8,na.rm = T)
m9<- mean(anc9,na.rm = T)
m10<- mean(anc10,na.rm = T)
m11<- mean(anc11,na.rm = T)
m12<- mean(anc12,na.rm = T)

values<-na.omit(c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12))

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CHH.",acc,".lineplot_medians_anc_and_freq.pdf",sep=""),height =3,width = 2.5)
par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
plot(1:4,c(m1,m4,m7,m10),pch=19, type="o", col="#3274B6", lwd=1.2,xaxt='n', xlab="frequency", main=paste("CHH meth, acc",acc),ylab="mean CHHme level",cex.main=0.8,ylim=c(min(values)-0.005,max(values)+0.005),xlim=c(0.5,4.5))
axis(1,at=1:4,labels=c("27","25-26","4-24","1-3"),las=2)
points(1:4,c(m2,m5,m7,m11), pch=19,type="o",lwd=1.2, col="#9EC2E6")
points(1:4,c(m3,m6,m9,m12), pch=19,type="o", lwd=1.2,col='#F7CBAD')


#add p values 
###########################
a1=anc1
a2=anc2
rm(d)
rm(w)
if (length(a1)>10 & length(a2)>10){
m1=mean(a1,na.rm = T)
m2=mean(a2,na.rm = T)
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
if (d<0.01) {text(b,x=1,y=posit,cex=0.8)}
}
a1=anc2
a2=anc3
rm(d)
rm(w)
if (length(a1)>10 & length(a2)>10){
m1=mean(a1,na.rm = T)
m2=mean(a2,na.rm = T)
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
if (d<0.01) {text(b,x=1,y=posit,cex=0.8)}
}
a1=anc4
a2=anc5
rm(d)
rm(w)
if (length(a1)>10 & length(a2)>10){
  m1=mean(a1,na.rm = T)
m2=mean(a2,na.rm = T)
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
if (d<0.01) {text(b,x=2,y=posit,cex=0.8)}
}
a1=anc5
a2=anc6
rm(d)
rm(w)
if (length(a1)>10 & length(a2)>10){
m1=mean(a1,na.rm = T)
m2=mean(a2,na.rm = T)
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
if (d<0.01) {text(b,x=2,y=posit,cex=0.8)}
}
a1=anc7
a2=anc8
rm(d)
rm(w)
if (length(a1)>10 & length(a2)>10){
m1=mean(a1,na.rm = T)
m2=mean(a2,na.rm = T)
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
if (d<0.01) {text(b,x=3,y=posit,cex=0.8)}
}
a1=anc8
a2=anc9
rm(d)
rm(w)
if (length(a1)>10 & length(a2)>10){
m1=mean(a1,na.rm = T)
m2=mean(a2,na.rm = T)
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
if (d<0.01) {text(b,x=3,y=posit,cex=0.8)}
}
a1=anc10
a2=anc11
rm(d)
rm(w)
if (length(a1)>10 & length(a2)>10){
m1=mean(a1,na.rm = T)
m2=mean(a2,na.rm = T)
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
if (d<0.01) {text(b,x=4,y=posit,cex=0.8)}
}
a1=anc11
a2=anc12
rm(d)
rm(w)
if (length(a1)>10 & length(a2)>10){
m1=mean(a1,na.rm = T)
m2=mean(a2,na.rm = T)
w<-flexible_wilcox(a1,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
if (d<0.01) {text(b,x=4,y=posit,cex=0.8)}
}

###########################
dev.off()


}
################################################

#boxplots and lineplots for TPM 27 accessions 
################################################
for (i in 1:27) {
  acc=accessions_w_220011$V1[i]
  
  a1<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% pc.core,paste("R",acc,sep=".")]
  a2<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% pc.highfreq,paste("R",acc,sep=".")]
  a3<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% pc.midfreq,paste("R",acc,sep=".")]
  a4<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% pc.low_and_priv,paste("R",acc,sep=".")]
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/TPM.R.",acc,".freq_cat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  boxplot(a1,a2,a3,a4,notch = T, names=c("27","25-26","4-24","1-3"),xlab="frequency",las=2,main=paste("expression,",acc),outline = F, pch=1,outcol=alpha("black",alpha=0.3), col=colorMap,cex.main=0.8,ylim=c(0,35),ylab="expression, TPM")
 
  #add p values 
  ###########################
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=15,cex=0.8)
  w<-flexible_wilcox(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=15,cex=0.8)
  w<-flexible_wilcox(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=15,cex=0.8)
  ###########################
  dev.off()
  
#ancestral status   
###########################
  b1<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_ancestral&TPMs.by_full_loci$gene %in% pcgenes,paste("R",acc,sep=".")]
  b2<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_lyrata_seqsim&TPMs.by_full_loci$gene %in% pcgenes,paste("R",acc,sep=".")]
  b3<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_non_ancestral&TPMs.by_full_loci$gene %in% pcgenes,paste("R",acc,sep=".")]
  
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/TPM.R.",acc,".ancestral_cat.pdf",sep=""),height =3,width = 2)
  par(mar=c(6,3,3,2),mgp=c(3,1,0)) 
  boxplot(b1, b2, b3,notch = T, las=2,main=paste("expression, acc",acc),outline = F, col=c('#3274B6','#9EC2E6','#F7CBAD'),names=c("full","seq","no"),xlab="ancestry status",cex.main=0.8,ylab="expression, TPM")
  #add p values 
  ###########################
  w<-flexible_wilcox(b1,b2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=20,cex=0.8)
  w<-flexible_wilcox(b3,b2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=20,cex=0.8)
  
  
  ###########################
  dev.off()
  }
################################################

#boxplots and lineplots for 24nt sRNA 
################################################
for (i in 2:15) {
  acc=as.character(accessions_miRNA$V1[i-1])
  
  a1<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% pc.core,acc]
  a2<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% pc.highfreq,acc]
  a3<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% pc.midfreq,acc]
  a4<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% pc.low_and_priv,acc]
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/24nt.",acc,".freq_cat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  boxplot(a1,a2,a3,a4,notch = T, names=c("27","25-26","4-24","1-3"),xlab="frequency",las=2,main=paste("24nt sRNA targeting,",acc),outline = F, pch=1,outcol=alpha("black",alpha=0.3), col=colorMap,cex.main=0.8,ylab="24nt coverage, RPM")
  
  #add p values 
  ###########################
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=0.1,cex=0.8)
  w<-flexible_wilcox(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=0.1,cex=0.8)
  w<-flexible_wilcox(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=0.1,cex=0.8)
  ###########################
  dev.off()
  
  a1<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% te.core,acc]
  a2<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% te.highfreq,acc]
  a3<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% te.midfreq,acc]
  a4<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% te.low_and_priv,acc]
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/24nt.TEgenes.",acc,".freq_cat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  boxplot(a1,a2,a3,a4,notch = T, names=c("27","25-26","4-24","1-3"),xlab="frequency",las=2,main=paste("24nt sRNA targeting TE genes,",acc),outline = F, pch=1,outcol=alpha("black",alpha=0.3), col=colorMap,cex.main=0.8,ylab="24nt coverage, RPM")
  
  #add p values 
  ###########################
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=0.3,cex=0.8)
  w<-flexible_wilcox(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=0.3,cex=0.8)
  w<-flexible_wilcox(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=0.3,cex=0.8)
  ###########################
  dev.off()
  
  
  b1<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_ancestral &sRNA.24nt.coverage.loci$gene %in% pcgenes ,acc]
  b2<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_lyrata_seqsim&sRNA.24nt.coverage.loci$gene %in% pcgenes,acc]
  b3<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_non_ancestral&sRNA.24nt.coverage.loci$gene %in% pcgenes,acc]
  
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/24nt.",acc,".ancestral_cat.pdf",sep=""),height =3,width = 2)
  par(mar=c(6,3,3,2),mgp=c(3,1,0)) 
  boxplot(b1, b2, b3,notch = T, las=2,main=paste("24nt sRNA targeting, acc",acc),outline = F, col=c('#3274B6','#9EC2E6','#F7CBAD'),names=c("full","seq","no"),xlab="ancestry status",cex.main=0.8,ylab="24nt coverage, RPM")
  #add p values 
  ###########################
  w<-flexible_wilcox(b1,b2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=0.04,cex=0.8)
  w<-flexible_wilcox(b3,b2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=0.04,cex=0.8)
  
  
  ###########################
  dev.off()
  
  #line plots for 12 categories by freq-ancestry
  
  #fixed
  anc1<-na.omit(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_ancestral & sRNA.24nt.coverage.loci$gene %in% pc.core,acc])
  anc2<-na.omit(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_lyrata_seqsim & sRNA.24nt.coverage.loci$gene %in% pc.core,acc])
  anc3<-na.omit(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_non_ancestral & sRNA.24nt.coverage.loci$gene %in% pc.core,acc])
  
  #softcore
  anc4<-na.omit(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_ancestral & sRNA.24nt.coverage.loci$gene %in% pc.highfreq,acc])
  anc5<-na.omit(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_lyrata_seqsim & sRNA.24nt.coverage.loci$gene %in% pc.highfreq,acc])
  anc6<-na.omit(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_non_ancestral & sRNA.24nt.coverage.loci$gene %in% pc.highfreq,acc])
  
  #midfreq
  anc7<-na.omit(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_ancestral & sRNA.24nt.coverage.loci$gene %in% pc.midfreq,acc])
  anc8<-na.omit(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_lyrata_seqsim & sRNA.24nt.coverage.loci$gene %in% pc.midfreq,acc])
  anc9<-na.omit(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_non_ancestral & sRNA.24nt.coverage.loci$gene %in% pc.midfreq,acc])
  
  #lowfreq and private
  anc10<-na.omit(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_ancestral & sRNA.24nt.coverage.loci$gene %in% pc.low_and_priv,acc])
  anc11<-na.omit(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_lyrata_seqsim & sRNA.24nt.coverage.loci$gene %in% pc.low_and_priv,acc])
  anc12<-na.omit(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_non_ancestral & sRNA.24nt.coverage.loci$gene %in% pc.low_and_priv,acc])
  
  #median is always 0 so take mean
  m1<- mean(anc1,na.rm = T)
  m2<- mean(anc2,na.rm = T)
  m3<- mean(anc3,na.rm = T)
  m4<- mean(anc4,na.rm = T)
  m5<- mean(anc5,na.rm = T)
  m6<- mean(anc6,na.rm = T)
  m7<- mean(anc7,na.rm = T)
  m8<- mean(anc8,na.rm = T)
  m9<- mean(anc9,na.rm = T)
  m10<- mean(anc10,na.rm = T)
  m11<- mean(anc11,na.rm = T)
  m12<- mean(anc12,na.rm = T)
  
  values<-na.omit(c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12))
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/24nt.",acc,".lineplot_medians_anc_and_freq.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  plot(1:4,c(m1,m4,m7,m10),pch=19, type="o", col="#3274B6", lwd=1.2,xaxt='n', xlab="frequency", ylab="", main=paste("mean 24nt sRNA, acc",acc),cex.main=0.8,xlim=c(0.5,4.5),ylim=c(min(values)-0.05,max(values)+0.05))
  axis(1,at=1:4,labels=c("27","25-26","4-24","1-3"),las=2)
  points(1:4,c(m2,m5,m7,m11), pch=19,type="o",lwd=1.2, col="#9EC2E6")
  points(1:4,c(m3,m6,m9,m12), pch=19,type="o", lwd=1.2,col='#F7CBAD')
  
  
  #add p values 
  ###########################
  a1=anc1
  a2=anc2
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
    w<-wilcox.test(a1,a2)
  m1=mean(a1,na.rm = T)
  m2=mean(a2,na.rm = T)
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
  if (d<0.01) {text(b,x=1,y=posit,cex=0.8)}
  }
  a1=anc2
  a2=anc3
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
    w<-wilcox.test(a1,a2)
  m1=mean(a1,na.rm = T)
  m2=mean(a2,na.rm = T)
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
  if (d<0.01) {text(b,x=1,y=posit,cex=0.8)}
  }
  a1=anc4
  a2=anc5
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
    w<-wilcox.test(a1,a2)
  m1=mean(a1,na.rm = T)
  m2=mean(a2,na.rm = T)
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
  if (d<0.01) {text(b,x=2,y=posit,cex=0.8)}
  }
  a1=anc5
  a2=anc6
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
    w<-wilcox.test(a1,a2)
  m1=mean(a1,na.rm = T)
  m2=mean(a2,na.rm = T)
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
  if (d<0.01) {text(b,x=2,y=posit,cex=0.8)}
  }
  a1=anc7
  a2=anc8
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
    w<-wilcox.test(a1,a2)
  m1=mean(a1,na.rm = T)
  m2=mean(a2,na.rm = T)
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
  if (d<0.01) {text(b,x=3,y=posit,cex=0.8)}
  }
  a1=anc8
  a2=anc9
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
  w<-wilcox.test(a1,a2)
  m1=mean(a1,na.rm = T)
  m2=mean(a2,na.rm = T)
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
  if (d<0.01) {text(b,x=3,y=posit,cex=0.8)}
  }
  a1=anc10
  a2=anc11
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
    w<-wilcox.test(a1,a2)
  m1=mean(a1,na.rm = T)
  m2=mean(a2,na.rm = T)
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
  if (d<0.01) {text(b,x=4,y=posit,cex=0.8)}
  }
  a1=anc11
  a2=anc12
  rm(d)
  rm(w)
  if (length(a1)>10 & length(a2)>10){
    w<-wilcox.test(a1,a2)
  m1=mean(a1,na.rm = T)
  m2=mean(a2,na.rm = T)
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  posit=max(m1,m2)-(max(m1,m2)-min(m1,m2))/2
  if (d<0.01) {text(b,x=4,y=posit,cex=0.8)}
  }
  
  ###########################
  dev.off()
  
}
################################################



#filtered out genes vs those left 
#############
# frequency 
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/Suppl_filtered_genes_frequency.pdf",sep=""),height =3,width = 2.5)
a1<-gene_frequency$freq_loc[gene_frequency$group %in% pcgenes]
a2<-gene_frequency$freq_loc[gene_frequency$group %in% pcgenes_all & gene_frequency$group %in% filteroutgenes]
a3<-gene_frequency$freq_loc[gene_frequency$group %in% tegenes]
a4<-gene_frequency$freq_loc[gene_frequency$group %in% tegenes_all& gene_frequency$group %in% filteroutgenes]

a1<-gene_frequency$freq_loc[gene_frequency$group %in% c(pcgenes,tegenes)]
a2<-gene_frequency$freq_loc[gene_frequency$group %in% filteroutgenes]


par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
boxplot(a1,a2,a3,a4,notch = T, names=c("PC", "PC removed", "TE", "TE removed"),xlab="",las=2,main="gene frequency",outline = F, col=colorMap_new1,cex.main=0.8)

#add p values 
#
w<-flexible_wilcox(a3,a2)
d<-w$p.value
if (d<0.0000000001){b="***"} 
if (d>=0.0000000001 &d<0.00001 ){b="**"} 
if (d>=0.00001  &d<0.01){b="*"} 
if (d>=0.01){b="n.s."} 
text(b,x=2.5,y=2,cex=0.8)
##################################


###########################################################
#########################Extended Data figure " PC gene silencing by gene frequency and ancestral status"#####
### heatmaps with median values for all accessions 
###########################################################

#chipseq
silencing_summary_CHIP<-accessions_chip
rownames(silencing_summary_CHIP)<-silencing_summary_CHIP$V1
silencing_summary_CHIP$K9_median_core<-0
silencing_summary_CHIP$K9_median_softcore<-0
silencing_summary_CHIP$K9_median_midfreq<-0
silencing_summary_CHIP$K9_median_lowfreq<-0
#silencing_summary_CHIP$K9_median_priv<-0

silencing_summary_CHIP$K27_median_core<-0
silencing_summary_CHIP$K27_median_softcore<-0
silencing_summary_CHIP$K27_median_midfreq<-0
silencing_summary_CHIP$K27_median_lowfreq<-0
#silencing_summary_CHIP$K27_median_priv<-0


for ( i in 1:6){
  acc=as.character(accessions_chip$V1[i])
  
  silencing_summary_CHIP$K9_median_core[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.core,paste("K9",acc,sep=".")],na.rm= T)
  silencing_summary_CHIP$K9_median_softcore[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.highfreq,paste("K9",acc,sep=".")],na.rm= T)
  silencing_summary_CHIP$K9_median_midfreq[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.midfreq,paste("K9",acc,sep=".")],na.rm= T)
  silencing_summary_CHIP$K9_median_lowfreq[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.low_and_priv,paste("K9",acc,sep=".")],na.rm= T)

  
  try(silencing_summary_CHIP$K27_median_core[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.core,paste("K27",acc,sep=".")],na.rm= T))
  try(silencing_summary_CHIP$K27_median_softcore[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.highfreq,paste("K27",acc,sep=".")],na.rm= T))
      try( silencing_summary_CHIP$K27_median_midfreq[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.midfreq,paste("K27",acc,sep=".")],na.rm= T))
           try( silencing_summary_CHIP$K27_median_lowfreq[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% pc.low_and_priv,paste("K27",acc,sep=".")],na.rm= T))
                
  
}
install.packages("pheatmap")
library(pheatmap)

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.K9.6acc.pdf",sep=""),height =2.5,width = 3)
pheatmap(silencing_summary_CHIP[,2:5],cluster_rows = F,cluster_cols =F,main="median H3K9me2 level, leaves",labels_col = c("27","25-26","4-24","1-3"),scale = "row",cex.main=0.7)
dev.off()

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.K27.5acc.pdf",sep=""),height =2.5,width = 3)
pheatmap(silencing_summary_CHIP[c(1,2,3,5,6),6:9],cluster_rows = F,cluster_cols =F,main="median H3K27me3 level, leaves",labels_col = c("27","25-26","4-24","1-3"),scale = "row",cex.main=0.7)
#mtext("scaled by row",side=3, line=0, cex=0.6)
dev.off()


# ancestral status
install.packages("pheatmap")
library(pheatmap)
#chipseq
silencing_summary_CHIP<-accessions_chip
rownames(silencing_summary_CHIP)<-silencing_summary_CHIP$V1
silencing_summary_CHIP$K9_median_anc<-0
silencing_summary_CHIP$K9_median_seq<-0
silencing_summary_CHIP$K9_median_non<-0
silencing_summary_CHIP$K27_median_anc<-0
silencing_summary_CHIP$K27_median_seq<-0
silencing_summary_CHIP$K27_median_non<-0

for ( i in 1:6){
  acc=as.character(accessions_chip$V1[i])
  
  silencing_summary_CHIP$K9_median_anc[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_ancestral & chipseq.loci.quantstan$gene %in% pcgenes,paste("K9",acc,sep=".")],na.rm= T)
  silencing_summary_CHIP$K9_median_seq[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_lyrata_seqsim& chipseq.loci.quantstan$gene %in% pcgenes,paste("K9",acc,sep=".")],na.rm= T)
  silencing_summary_CHIP$K9_median_non[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_non_ancestral& chipseq.loci.quantstan$gene %in% pcgenes,paste("K9",acc,sep=".")],na.rm= T)
  
  
  try(silencing_summary_CHIP$K27_median_anc[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_ancestral& chipseq.loci.quantstan$gene %in% pcgenes,paste("K27",acc,sep=".")],na.rm= T))
  try(silencing_summary_CHIP$K27_median_seq[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_lyrata_seqsim& chipseq.loci.quantstan$gene %in% pcgenes,paste("K27",acc,sep=".")],na.rm= T))
  try( silencing_summary_CHIP$K27_median_non[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_non_ancestral& chipseq.loci.quantstan$gene %in% pcgenes,paste("K27",acc,sep=".")],na.rm= T))
 
}


pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.K9.6acc.ancestralPC.pdf",sep=""),height =2.5,width = 3)
pheatmap(silencing_summary_CHIP[,2:4],cluster_rows = F,cluster_cols =F,main="median H3K9me2 level, leaves",labels_col = c("full","seq","no"),scale = "row",cex.main=0.7)
dev.off()

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.K27.5acc.ancestralPC.pdf",sep=""),height =2.5,width = 3)
pheatmap(silencing_summary_CHIP[c(1,2,3,5,6),5:7],cluster_rows = F,cluster_cols =F,main="median H3K27me3 level, leaves",labels_col = c("full","seq","no"),scale = "row",cex.main=0.7)
#mtext("scaled by row",side=3, line=0, cex=0.6)
dev.off()




#methylation

silencing_summary_table_meth<-accessions_methylation
rownames(silencing_summary_table_meth)<-silencing_summary_table_meth$V1
silencing_summary_table_meth$CG_median_core<-0
silencing_summary_table_meth$CG_median_softcore<-0
silencing_summary_table_meth$CG_median_midfreq<-0
silencing_summary_table_meth$CG_median_lowfreq<-0
#silencing_summary_table_meth$CG_median_priv<-0

silencing_summary_table_meth$CHH_median_core<-0
silencing_summary_table_meth$CHH_median_softcore<-0
silencing_summary_table_meth$CHH_median_midfreq<-0
silencing_summary_table_meth$CHH_median_lowfreq<-0
#silencing_summary_table_meth$CHH_median_priv<-0

silencing_summary_table_meth$CHG_median_core<-0
silencing_summary_table_meth$CHG_median_softcore<-0
silencing_summary_table_meth$CHG_median_midfreq<-0
silencing_summary_table_meth$CHG_median_lowfreq<-0
#silencing_summary_table_meth$CHG_median_priv<-0

for ( i in 1:12){
  acc=as.character(accessions_methylation$V1[i])
  
  silencing_summary_table_meth$CG_median_core[i]<-median(CG.loci[CG.loci$gene %in% pc.core, acc],na.rm = T)
  silencing_summary_table_meth$CG_median_softcore[i]<-median(CG.loci[CG.loci$gene %in% pc.highfreq, acc],na.rm = T)
  silencing_summary_table_meth$CG_median_midfreq[i]<-median(CG.loci[CG.loci$gene %in% pc.midfreq, acc],na.rm = T)
  silencing_summary_table_meth$CG_median_lowfreq[i]<-median(CG.loci[CG.loci$gene %in% pc.low_and_priv, acc],na.rm = T)
#  silencing_summary_table_meth$CG_median_priv[i]<-median(CG.loci[CG.loci$gene %in% private.loci.PC, acc],na.rm = T)
  
  silencing_summary_table_meth$CHG_median_core[i]<-median(CHG.loci[CHG.loci$gene %in% pc.core, acc],na.rm = T)
  silencing_summary_table_meth$CHG_median_softcore[i]<-median(CHG.loci[CHG.loci$gene %in% pc.highfreq, acc],na.rm = T)
  silencing_summary_table_meth$CHG_median_midfreq[i]<-median(CHG.loci[CHG.loci$gene %in% pc.midfreq, acc],na.rm = T)
  silencing_summary_table_meth$CHG_median_lowfreq[i]<-median(CHG.loci[CHG.loci$gene %in% pc.low_and_priv, acc],na.rm = T)
 # silencing_summary_table_meth$CHG_median_priv[i]<-median(CHG.loci[CHG.loci$gene %in% private.loci.PC, acc],na.rm = T)
  
  silencing_summary_table_meth$CHH_median_core[i]<-median(CHH.loci[CHH.loci$gene %in% pc.core, acc],na.rm = T)
  silencing_summary_table_meth$CHH_median_softcore[i]<-median(CHH.loci[CHH.loci$gene %in% pc.highfreq, acc],na.rm = T)
  silencing_summary_table_meth$CHH_median_midfreq[i]<-median(CHH.loci[CHH.loci$gene %in% pc.midfreq, acc],na.rm = T)
  silencing_summary_table_meth$CHH_median_lowfreq[i]<-median(CHH.loci[CHH.loci$gene %in% pc.low_and_priv, acc],na.rm = T)
 # silencing_summary_table_meth$CHH_median_priv[i]<-median(CHH.loci[CHH.loci$gene %in% private.loci.PC, acc],na.rm = T)
  
}

library(pheatmap)
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.CG.12acc.pdf",sep=""),height =3,width = 3)
pheatmap(silencing_summary_table_meth[,2:5],cluster_rows = F,cluster_cols =F,main="median CG methylation",labels_col = c("27","25-26","4-24","1-3"), scale="row")
dev.off()

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.CHG.12acc.pdf",sep=""),height =3,width = 3)
pheatmap(silencing_summary_table_meth[,10:13],cluster_rows = F,cluster_cols =F,main="median CHG methylation",labels_col = c("27","25-26","4-24","1-3"), scale="row")
dev.off()
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.CHH.12acc.pdf",sep=""),height =3,width = 3)
pheatmap(silencing_summary_table_meth[,6:9],cluster_rows = F,cluster_cols =F,main="median CHH methylation",display_numbers = F,labels_col = c("27","25-26","4-24","1-3"), scale="row")
dev.off()


#ancestral 


#methylation

silencing_summary_table_meth<-accessions_methylation
rownames(silencing_summary_table_meth)<-silencing_summary_table_meth$V1
silencing_summary_table_meth$CG_median_anc<-0
silencing_summary_table_meth$CG_median_seq<-0
silencing_summary_table_meth$CG_median_non<-0

silencing_summary_table_meth$CHG_median_anc<-0
silencing_summary_table_meth$CHG_median_seq<-0
silencing_summary_table_meth$CHG_median_non<-0


silencing_summary_table_meth$CHH_median_anc<-0
silencing_summary_table_meth$CHH_median_seq<-0
silencing_summary_table_meth$CHH_median_non<-0


for ( i in 1:12){
  acc=as.character(accessions_methylation$V1[i])
  
  silencing_summary_table_meth$CG_median_anc[i]<-median(CG.loci[CG.loci$gene %in% intersect(genes_ancestral,pcgenes) , acc],na.rm = T)
  silencing_summary_table_meth$CG_median_seq[i]<-median(CG.loci[CG.loci$gene %in% intersect(genes_lyrata_seqsim,pcgenes), acc],na.rm = T)
  silencing_summary_table_meth$CG_median_non[i]<-median(CG.loci[CG.loci$gene %in% intersect(genes_non_ancestral,pcgenes), acc],na.rm = T)
  
  silencing_summary_table_meth$CHG_median_anc[i]<-median(CHG.loci[CHG.loci$gene %in% intersect(genes_ancestral,pcgenes), acc],na.rm = T)
  silencing_summary_table_meth$CHG_median_seq[i]<-median(CHG.loci[CHG.loci$gene %in% intersect(genes_lyrata_seqsim,pcgenes), acc],na.rm = T)
  silencing_summary_table_meth$CHG_median_non[i]<-median(CHG.loci[CHG.loci$gene %in% intersect(genes_non_ancestral,pcgenes), acc],na.rm = T)
  
  
  silencing_summary_table_meth$CHH_median_anc[i]<-median(CHH.loci[CHH.loci$gene %in% intersect(genes_ancestral,pcgenes), acc],na.rm = T)
  silencing_summary_table_meth$CHH_median_seq[i]<-median(CHH.loci[CHH.loci$gene %in% intersect(genes_lyrata_seqsim,pcgenes), acc],na.rm = T)
  silencing_summary_table_meth$CHH_median_non[i]<-median(CHH.loci[CHH.loci$gene %in% intersect(genes_non_ancestral,pcgenes), acc],na.rm = T)
  
}

library(pheatmap)
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.CG.12acc.PCanc.pdf",sep=""),height =3,width = 3)
pheatmap(silencing_summary_table_meth[,2:4],cluster_rows = F,cluster_cols =F,main="median CG methylation",labels_col = c("full","seq","no"), scale="row")
dev.off()

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.CHG.12acc.PCanc.pdf",sep=""),height =3,width = 3)
pheatmap(silencing_summary_table_meth[,5:7],cluster_rows = F,cluster_cols =F,main="median CHG methylation",labels_col = c("full","seq","no"),scale="row")
dev.off()
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.CHH.12acc.PCanc.pdf",sep=""),height =3,width = 3)
pheatmap(silencing_summary_table_meth[,8:10],cluster_rows = F,cluster_cols =F,main="median CHH methylation",display_numbers = F,labels_col = c("full","seq","no"), scale="row")
dev.off()













#24nt sRNA

silencing_summary_table_24nt<-accessions_miRNA
rownames(silencing_summary_table_24nt)<-silencing_summary_table_24nt$V1
silencing_summary_table_24nt$sRNA_median_core<-0
silencing_summary_table_24nt$sRNA_median_softcore<-0
silencing_summary_table_24nt$sRNA_median_midfreq<-0
silencing_summary_table_24nt$sRNA_median_lowfreq<-0
#silencing_summary_table_24nt$sRNA_median_priv<-0

for ( i in 1:14){
  acc=as.character(accessions_miRNA$V1[i])
 
  silencing_summary_table_24nt$sRNA_median_core[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% pc.core, acc],na.rm = T)
  silencing_summary_table_24nt$sRNA_median_softcore[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% pc.highfreq, acc],na.rm = T)
  silencing_summary_table_24nt$sRNA_median_midfreq[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% pc.midfreq, acc],na.rm = T)
  silencing_summary_table_24nt$sRNA_median_lowfreq[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% pc.low_and_priv, acc],na.rm = T)
#  silencing_summary_table_24nt$sRNA_median_priv[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% private.loci.PC, acc],na.rm = T)
  
}

library(pheatmap)
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.24nt.14acc.pdf",sep=""),height =3,width = 3)
pheatmap(silencing_summary_table_24nt[,2:5],cluster_rows = F,cluster_cols =F,main="24nt sRNA coverage, flowers",labels_col = c("27","25-26","4-24","1-3"),scale = "row")
#pheatmap(silencing_summary_table_24nt[,2:5],cluster_rows = F,cluster_cols =F,main="24nt sRNA coverage, flowers",scale = "row")
dev.off()




#24nt sRNA - ancestral

silencing_summary_table_24nt<-accessions_miRNA
rownames(silencing_summary_table_24nt)<-silencing_summary_table_24nt$V1
silencing_summary_table_24nt$sRNA_median_anc<-0
silencing_summary_table_24nt$sRNA_median_seq<-0
silencing_summary_table_24nt$sRNA_median_no<-0

for ( i in 1:14){
  acc=as.character(accessions_miRNA$V1[i])
  
  silencing_summary_table_24nt$sRNA_median_anc[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% intersect(genes_ancestral,pcgenes), acc],na.rm = T)
  silencing_summary_table_24nt$sRNA_median_seq[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% intersect(genes_lyrata_seqsim,pcgenes), acc],na.rm = T)
  silencing_summary_table_24nt$sRNA_median_no[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% intersect(genes_non_ancestral,pcgenes), acc],na.rm = T)
 
}

library(pheatmap)
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.24nt.14acc.PCanc.pdf",sep=""),height =3,width = 3)
pheatmap(silencing_summary_table_24nt[,2:4],cluster_rows = F,cluster_cols =F,main="24nt sRNA coverage, flowers",labels_col = c("full","seq","no"),scale = "row")
#pheatmap(silencing_summary_table_24nt[,2:5],cluster_rows = F,cluster_cols =F,main="24nt sRNA coverage, flowers",scale = "row")
dev.off()






#expression 
silencing_summary_table_TPM<-accessions_w_220011
rownames(silencing_summary_table_TPM)<-silencing_summary_table_TPM$V1
silencing_summary_table_TPM$TPM_rosette_median_core<-0
silencing_summary_table_TPM$TPM_rosette_median_softcore<-0
silencing_summary_table_TPM$TPM_rosette_median_midfreq<-0
silencing_summary_table_TPM$TPM_rosette_median_lowfreq<-0
#silencing_summary_table_TPM$TPM_rosette_median_priv<-0
for ( i in 1:27){
  acc=as.character(accessions_w_220011$V1[i])
  silencing_summary_table_TPM$TPM_rosette_median_core[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% pc.core,paste("R",acc,sep=".")],na.rm = T)
  silencing_summary_table_TPM$TPM_rosette_median_softcore[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% pc.highfreq,paste("R",acc,sep=".")],na.rm = T)
  silencing_summary_table_TPM$TPM_rosette_median_midfreq[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% pc.midfreq,paste("R",acc,sep=".")],na.rm = T)
  silencing_summary_table_TPM$TPM_rosette_median_lowfreq[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% pc.low_and_priv,paste("R",acc,sep=".")],na.rm = T)
#  silencing_summary_table_TPM$TPM_rosette_median_priv[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% private.loci.PC,paste("R",acc,sep=".")],na.rm = T)
  
}

library(pheatmap)
#pheatmap(silencing_summary_table_TPM[,2:5],cluster_rows = F,cluster_cols =F,main="expression in rosette",scale = "row")
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.TPM.27acc.pdf",sep=""),height =4,width = 3)
pheatmap(silencing_summary_table_TPM[,2:5],cluster_rows = F,cluster_cols =F,main="expression in rosette (TPM)",labels_col = c("27","25-26","4-24","1-3"))
dev.off()




#expression ancestral categories
silencing_summary_table_TPM_anc<-accessions_w_220011
rownames(silencing_summary_table_TPM_anc)<-silencing_summary_table_TPM_anc$V1
silencing_summary_table_TPM_anc$anc<-0
silencing_summary_table_TPM_anc$seq<-0
silencing_summary_table_TPM_anc$no<-0

for ( i in 1:27){
  acc=as.character(accessions_w_220011$V1[i])
  silencing_summary_table_TPM_anc$anc[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_ancestral,paste("R",acc,sep=".")],na.rm = T)
  silencing_summary_table_TPM_anc$seq[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_lyrata_seqsim,paste("R",acc,sep=".")],na.rm = T)
  silencing_summary_table_TPM_anc$no[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_non_ancestral,paste("R",acc,sep=".")],na.rm = T)
 
}

library(pheatmap)
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.TPM_PC_ancest_cat.27acc.pdf",sep=""),height =4,width = 3)
pheatmap(silencing_summary_table_TPM_anc[,2:4],cluster_rows = F,cluster_cols =F,main="expression in rosette (TPM)",labels_col = c("full","seq","no"))
dev.off()

###############################################


###########################################################
#########################FIGURE - Supplemental figure - NEW genes#####
###########################################################

##################################################
### NEW GENES
###################################################

#new gene categories
length(genes_PC_new)#2661
length(genes_TE_new)#5109
length(new_tair_high) #744
length(new_tairP_medium) #670
length(new_brass)#817
length(new_other) #9
length(new_nohomol)#421

length(genes_PC_TAIR10)#25477
length(genes_TE_TAIR10)#565


##################################################
#do new genes come from duplications? 
##################################################


# does the gene have a copy? 

cnvgenes<- as.vector(copyN$group[copyN$sd_CN>0])
duplicatedgenes<- as.vector(copyN$group[copyN$max_CN>=2])


length(intersect(cnvgenes,new_tair_high))#151
length(intersect(cnvgenes,new_tairP_medium))#137
length(intersect(cnvgenes,new_brass))#107
length(intersect(cnvgenes,new_nohomol))#24

length(intersect(cnvgenes,genes_PC_TAIR10))#1485
length(intersect(cnvgenes,genes_PC_new))#419
length(intersect(cnvgenes,genes_TE_TAIR10))#179
length(intersect(cnvgenes,genes_TE_new))#1704



length(intersect(duplicatedgenes,new_tair_high))#596
length(intersect(duplicatedgenes,new_tairP_medium))#310
length(intersect(duplicatedgenes,new_brass))#281
length(intersect(duplicatedgenes,new_nohomol))#100

length(intersect(duplicatedgenes,genes_PC_TAIR10))#2043
length(intersect(duplicatedgenes,genes_PC_new))#1290
length(intersect(duplicatedgenes,genes_TE_TAIR10))#260
length(intersect(duplicatedgenes,genes_TE_new))#3547

#### new genes - tandem dup 
length(intersect(genesnames_tandem_dup,new_tair_high))#117
length(intersect(genesnames_tandem_dup,new_tairP_medium))#47
length(intersect(genesnames_tandem_dup,new_brass))#9
length(intersect(genesnames_tandem_dup,new_nohomol))#14

length(intersect(genesnames_tandem_dup,genes_PC_TAIR10))#617
length(intersect(genesnames_tandem_dup,genes_PC_new))#189
length(intersect(genesnames_tandem_dup,genes_TE_TAIR10))#40
length(intersect(genesnames_tandem_dup,genes_TE_new))#321

##################################################################

############################
# silencing of new genes 
#######################################

colorMap_new
colorMap_new1
#boxplots for H3K9me2 6 accessions - new /old and 4 new categories 
########
for (i in 1:6) {
  acc=accessions_chip$V1[i]
  a1<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_PC_TAIR10,paste("K9",acc,sep=".")]
  a2<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_PC_new,paste("K9",acc,sep=".")]
  a3<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_TE_TAIR10,paste("K9",acc,sep=".")]
  a4<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_TE_new,paste("K9",acc,sep=".")]
  
    pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/K9.",acc,".old_new.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  boxplot(a1,a2,a3,a4,notch = T, names=names(colorMap_new1),xlab="",las=2,main=paste("H3K9me2, acc",acc),outline = F, col=colorMap_new1,cex.main=0.8)
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=3,cex=0.8)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=3,cex=0.8)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=3,cex=0.8)
       ###########################
  dev.off()
  
  
  a1<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% new_tair_high,paste("K9",acc,sep=".")]
  a2<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% new_tairP_medium,paste("K9",acc,sep=".")]
  a3<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% new_brass,paste("K9",acc,sep=".")]
  a4<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% new_nohomol,paste("K9",acc,sep=".")]
 
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/K9.",acc,".new_4Cat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  boxplot(a1,a2,a3,a4,notch = T, names=names(colorMap_new[2:5]),xlab="",las=2,main=paste("H3K9me2, acc",acc),outline = F, col=colorMap_new[2:5],cex.main=0.8,cex.axis=0.6,ylim=c(-1.5,3.2))
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=3,cex=0.8)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=3,cex=0.8)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=3,cex=0.8)
    ###########################
  dev.off()
  
  
}
########


#boxplots  H3K27me3 5 accessions - new /old and 4 new categories 
########
for (i in 1:5) {
  acc=c("6909","6024","9905","9888",'1741')[i]
  
  a1<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_PC_TAIR10,paste("K27",acc,sep=".")]
  a2<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_PC_new,paste("K27",acc,sep=".")]
  a3<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_TE_TAIR10,paste("K27",acc,sep=".")]
  a4<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_TE_new,paste("K27",acc,sep=".")]
  #a5<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% newgen_telike,paste("K27",acc,sep=".")]
  
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/K27.",acc,".new_old.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  boxplot(a1,a2,a3,a4,notch = T, names=names(colorMap_new1),xlab="",las=2,main=paste("H3K27me3, acc",acc),outline = F, col=colorMap_new1,cex.main=0.8)
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=2,cex=0.8)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=2,cex=0.8)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=2,cex=0.8)
 
  
  ###########################
  dev.off()
 
  
  
  a1<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% new_tair_high,paste("K27",acc,sep=".")]
  a2<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% new_tairP_medium,paste("K27",acc,sep=".")]
  a3<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% new_brass,paste("K27",acc,sep=".")]
  a4<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% new_nohomol,paste("K27",acc,sep=".")]
  #a5<-chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% newgen_telike,paste("K27",acc,sep=".")]
  
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/K27.",acc,".new_cat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  boxplot(a1,a2,a3,a4,notch = T, names=names(colorMap_new[2:5]),xlab="",las=2,main=paste("H3K27me3, acc",acc),outline = F, col=colorMap_new[2:5],cex.main=0.8,cex.axis=0.6)
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=2,cex=0.8)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=2,cex=0.8)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=2,cex=0.8)
  
  ###########################
  dev.off()
 }
########

#boxplots for CG,CHG, CHH in 12 accessions  - new /old and 4 new categories 
########
for (i in 2:13) {
  acc=as.character(accessions_methylation$V1[i-1])
  
  
  a1<-CG.loci[CG.loci$gene %in% genes_PC_TAIR10,i]
  a2<-CG.loci[CG.loci$gene %in% genes_PC_new,i]
  a3<-CG.loci[CG.loci$gene %in% genes_TE_TAIR10,i]
  a4<-CG.loci[CG.loci$gene %in% genes_TE_new,i]
  
 
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CG.",acc,".old_new.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  library(vioplot)
  vioplot(a1,a2,a3,a4,col=colorMap_new1, names=names(colorMap_new1),ylab="methylation level",las=2,main=paste("CG, leaves acc",acc) ,cex.main=0.8)
  
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=1)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=1)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=1)
   ###########################
  dev.off()
  
 
  
  
  # 4 categories of new PC genes
  a1<-CG.loci[CG.loci$gene %in% new_tair_high,i]
  a2<-CG.loci[CG.loci$gene %in% new_tairP_medium,i]
  a3<-CG.loci[CG.loci$gene %in% new_brass,i]
  a4<-CG.loci[CG.loci$gene %in% new_nohomol,i]
  
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CG.",acc,".newPC_4cat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  library(vioplot)
  vioplot(a1,a2,a3,a4,col=colorMap_new[2:5], names=names(colorMap_new[2:5]),ylab="methylation level",las=2,main=paste("CG, leaves acc",acc) ,cex.main=0.8)
  
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=1)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=1)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=1)
  ###########################
  dev.off()
  
  
  a1<-CHH.loci[CHH.loci$gene %in% genes_PC_TAIR10,i]
  a2<-CHH.loci[CHH.loci$gene %in% genes_PC_new,i]
  a3<-CHH.loci[CHH.loci$gene %in% genes_TE_TAIR10,i]
  a4<-CHH.loci[CHH.loci$gene %in% genes_TE_new,i]
  #a5<-CHH.loci[CHH.loci$gene %in% newgen_telike,i]
  
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CHH.",acc,".new_old.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  library(vioplot)
  boxplot(a1,a2,a3,a4,col=colorMap_new1, names=names(colorMap_new1),ylab="methylation level",las=2,main=paste("CHH, leaves acc",acc) ,cex.main=0.8,notch = T,outline = F)
  
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=0.15)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=0.15)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=0.15)
  
  ###########################
  dev.off()
  
  
  a1<-CHH.loci[CHH.loci$gene %in% new_tair_high,i]
  a2<-CHH.loci[CHH.loci$gene %in% new_tairP_medium,i]
  a3<-CHH.loci[CHH.loci$gene %in% new_brass,i]
  a4<-CHH.loci[CHH.loci$gene %in% new_nohomol,i]
  #a5<-CHH.loci[CHH.loci$gene %in% newgen_telike,i]
  
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CHH.",acc,".newPC_4cat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  library(vioplot)
  boxplot(a1,a2,a3,a4,col=colorMap_new[2:5], names=names(colorMap_new[2:5]),ylab="methylation level",las=2,main=paste("CHH, leaves acc",acc) ,cex.main=0.8,notch = T,outline = F)
  
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=0.05)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=0.05)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=0.05)
  
  ###########################
  dev.off()
  
  
  a1<-CHG.loci[CHG.loci$gene %in% genes_PC_TAIR10,i]
  a2<-CHG.loci[CHG.loci$gene %in% genes_PC_new,i]
  a3<-CHG.loci[CHG.loci$gene %in% genes_TE_TAIR10,i]
  a4<-CHG.loci[CHG.loci$gene %in% genes_TE_new,i]
  #a5<-CHH.loci[CHH.loci$gene %in% newgen_telike,i]
  
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CHG.",acc,".new_old.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  library(vioplot)
  boxplot(a1,a2,a3,a4,col=colorMap_new1, names=names(colorMap_new1),ylab="methylation level",las=2,main=paste("CHG., leaves acc",acc) ,cex.main=0.8,notch = T,outline = F)
  
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=0.15)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=0.15)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=0.15)
  
  ###########################
  dev.off()
  
  
  a1<-CHG.loci[CHG.loci$gene %in% new_tair_high,i]
  a2<-CHG.loci[CHG.loci$gene %in% new_tairP_medium,i]
  a3<-CHG.loci[CHG.loci$gene %in% new_brass,i]
  a4<-CHG.loci[CHG.loci$gene %in% new_nohomol,i]
  #a5<-CHH.loci[CHH.loci$gene %in% newgen_telike,i]
  
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/CHG.",acc,".newPC_4cat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  library(vioplot)
  boxplot(a1,a2,a3,a4,col=colorMap_new[2:5], names=names(colorMap_new[2:5]),ylab="methylation level",las=2,main=paste("CHG, leaves acc",acc) ,cex.main=0.8,notch = T,outline = F)
  
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=0.15)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=0.15)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=0.15)
  
  ###########################
  dev.off()
  }

###########


#boxplots for TPM 27 accessions 
########
for (i in 1:27) {
  acc=accessions_w_220011$V1[i]
  
  a1<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_PC_TAIR10,paste("R",acc,sep=".")]
  a2<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_PC_new,paste("R",acc,sep=".")]
  a3<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_TE_TAIR10,paste("R",acc,sep=".")]
  a4<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_TE_new,paste("R",acc,sep=".")]
#  a5<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% newgen_telike,paste("R",acc,sep=".")]
  
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/TPM.R.",acc,".old_new.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  
  boxplot(a1,a2,a3,a4,notch = T,col=colorMap_new1, names=names(colorMap_new1),las=2,main=paste("expression,",acc),outline = T, pch=1,outcol=alpha("black",alpha=0.3), cex.main=0.8,ylim=c(0,35),ylab="expression, TPM")
  
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=15,cex=0.8)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=15,cex=0.8)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=15,cex=0.8)
    ###########################
  dev.off()
  
  
 
  a1<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% new_tair_high,paste("R",acc,sep=".")]
  a2<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% new_tairP_medium,paste("R",acc,sep=".")]
  a3<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% new_brass,paste("R",acc,sep=".")]
  a4<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% new_nohomol,paste("R",acc,sep=".")]
  #  a5<-TPMs.by_full_loci[TPMs.by_full_loci$gene %in% newgen_telike,paste("R",acc,sep=".")]
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/TPM.R.",acc,".newgencat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  
  boxplot(a1,a2,a3,a4,notch = T,col=colorMap_new[2:5], names=names(colorMap_new[2:5]),las=2,main=paste("expression,",acc),outline = T, pch=1,outcol=alpha("black",alpha=0.3),cex.axis=0.6, cex.main=0.8,ylim=c(0,3),ylab="expression, TPM")
  
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=2,cex=0.8)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=2,cex=0.8)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
   #text(b,x=4.5,y=15,cex=0.8)
  ###########################
  dev.off() 
}
########


#boxplots and lineplots for 24nt sRNA 
########
for (i in 2:15) {
  acc=as.character(accessions_miRNA$V1[i-1])
  
  a1<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_PC_TAIR10,acc]
  a2<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_PC_new,acc]
  a3<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_TE_TAIR10,acc]
  a4<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_TE_new,acc]

  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/24nt.",acc,".newold.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  boxplot(a1,a2,a3,a4,notch = T, names=names(colorMap_new1),las=2,main=paste("24nt sRNA targeting,",acc),outline = F, pch=1,outcol=alpha("black",alpha=0.3), col=colorMap_new1,cex.main=0.8,ylab="24nt coverage, RPM")
  
  #add p values 
  ###########################
  w<-flexible_wilcox(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=0.3,cex=0.8)
  w<-flexible_wilcox(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=0.3,cex=0.8)
  w<-flexible_wilcox(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=0.3,cex=0.8)
   ###########################
  dev.off()
  
  
  
# new gene categories
  a1<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% new_tair_high,acc]
  a2<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% new_tairP_medium,acc]
  a3<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% new_brass,acc]
  a4<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% new_nohomol,acc]
  # a5<-sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% newgen_telike,acc]
  
  

  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/24nt.",acc,".newgen_cat.pdf",sep=""),height =3,width = 2.5)
  par(mar=c(6,3,3,2),mgp=c(4,1,0)) 
  boxplot(a1,a2,a3,a4,notch = T, col=colorMap_new[2:5], names=names(colorMap_new[2:5]),las=2,main=paste("24nt sRNA targeting,",acc),outline = F, pch=1,outcol=alpha("black",alpha=0.3), cex.axis=0.6,cex.main=0.8,ylab="24nt coverage, RPM")
  
  #add p values 
  ###########################
  w<-wilcox.test(a1,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=0.03,cex=0.8)
  w<-wilcox.test(a3,a2)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=0.03,cex=0.8)
  w<-wilcox.test(a3,a4)
  d<-w$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=0.03,cex=0.8)
  ###########################
  dev.off()
  
 }
########



Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/
############## 
### heatmaps with median values for all accessions 
##############

#chipseq
silencing_summary_CHIP<-accessions_chip
rownames(silencing_summary_CHIP)<-silencing_summary_CHIP$V1
silencing_summary_CHIP$K9_median_1<-0
silencing_summary_CHIP$K9_median_2<-0
silencing_summary_CHIP$K9_median_3<-0
silencing_summary_CHIP$K9_median_4<-0

silencing_summary_CHIP$K27_median_1<-0
silencing_summary_CHIP$K27_median_2<-0
silencing_summary_CHIP$K27_median_3<-0
silencing_summary_CHIP$K27_median_4<-0


for ( i in 1:6){
  acc=as.character(accessions_chip$V1[i])
  
  silencing_summary_CHIP$K9_median_1[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_PC_TAIR10,paste("K9",acc,sep=".")],na.rm= T)
  silencing_summary_CHIP$K9_median_2[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_PC_new,paste("K9",acc,sep=".")],na.rm= T)
  silencing_summary_CHIP$K9_median_3[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_TE_TAIR10,paste("K9",acc,sep=".")],na.rm= T)
  silencing_summary_CHIP$K9_median_4[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_TE_new,paste("K9",acc,sep=".")],na.rm= T)
  
 try( silencing_summary_CHIP$K27_median_1[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_PC_TAIR10,paste("K27",acc,sep=".")],na.rm= T))
      try( silencing_summary_CHIP$K27_median_2[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_PC_new,paste("K27",acc,sep=".")],na.rm= T))
           try( silencing_summary_CHIP$K27_median_3[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_TE_TAIR10,paste("K27",acc,sep=".")],na.rm= T))
                try(  silencing_summary_CHIP$K27_median_4[i]<-median(chipseq.loci.quantstan[chipseq.loci.quantstan$gene %in% genes_TE_new,paste("K27",acc,sep=".")],na.rm= T))
}

library(pheatmap)

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.K9.new_old.6acc.pdf",sep=""),height =2.5,width = 3)
pheatmap(silencing_summary_CHIP[,2:5],cluster_rows = F,cluster_cols =F,main="median H3K9me2 level, leaves",labels_col = names(colorMap_new1),scale="row")
dev.off()

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.K27.new_old.5acc.pdf",sep=""),height =2.5,width = 3)
pheatmap(silencing_summary_CHIP[c(1,2,3,5,6),6:9],cluster_rows = F,cluster_cols =F,main="median H3K27me3 level, leaves",labels_col = names(colorMap_new1),scale="row")
dev.off()


#methylation

silencing_summary_table_meth<-accessions_methylation
rownames(silencing_summary_table_meth)<-silencing_summary_table_meth$V1
silencing_summary_table_meth$CG_median_core<-0
silencing_summary_table_meth$CG_median_softcore<-0
silencing_summary_table_meth$CG_median_midfreq<-0
silencing_summary_table_meth$CG_median_lowfreq<-0
silencing_summary_table_meth$CG_median_priv<-0

silencing_summary_table_meth$CHH_median_core<-0
silencing_summary_table_meth$CHH_median_softcore<-0
silencing_summary_table_meth$CHH_median_midfreq<-0
silencing_summary_table_meth$CHH_median_lowfreq<-0
silencing_summary_table_meth$CHH_median_priv<-0

for ( i in 1:12){
  acc=as.character(accessions_methylation$V1[i])
  
  silencing_summary_table_meth$CG_median_core[i]<-median(CG.loci[CG.loci$gene %in% pc.core, acc],na.rm = T)
  silencing_summary_table_meth$CG_median_softcore[i]<-median(CG.loci[CG.loci$gene %in% pc.highfreq, acc],na.rm = T)
  silencing_summary_table_meth$CG_median_midfreq[i]<-median(CG.loci[CG.loci$gene %in% pc.midfreq, acc],na.rm = T)
  silencing_summary_table_meth$CG_median_lowfreq[i]<-median(CG.loci[CG.loci$gene %in% pc.low_and_priv, acc],na.rm = T)
  #  silencing_summary_table_meth$CG_median_priv[i]<-median(CG.loci[CG.loci$gene %in% private.loci.PC, acc],na.rm = T)
  
  silencing_summary_table_meth$CHG_median_core[i]<-median(CHG.loci[CHG.loci$gene %in% pc.core, acc],na.rm = T)
  silencing_summary_table_meth$CHG_median_softcore[i]<-median(CHG.loci[CHG.loci$gene %in% pc.highfreq, acc],na.rm = T)
  silencing_summary_table_meth$CHG_median_midfreq[i]<-median(CHG.loci[CHG.loci$gene %in% pc.midfreq, acc],na.rm = T)
  silencing_summary_table_meth$CHG_median_lowfreq[i]<-median(CHG.loci[CHG.loci$gene %in% pc.low_and_priv, acc],na.rm = T)
  # silencing_summary_table_meth$CHG_median_priv[i]<-median(CHG.loci[CHG.loci$gene %in% private.loci.PC, acc],na.rm = T)
  
  silencing_summary_table_meth$CHH_median_core[i]<-median(CHH.loci[CHH.loci$gene %in% pc.core, acc],na.rm = T)
  silencing_summary_table_meth$CHH_median_softcore[i]<-median(CHH.loci[CHH.loci$gene %in% pc.highfreq, acc],na.rm = T)
  silencing_summary_table_meth$CHH_median_midfreq[i]<-median(CHH.loci[CHH.loci$gene %in% pc.midfreq, acc],na.rm = T)
  silencing_summary_table_meth$CHH_median_lowfreq[i]<-median(CHH.loci[CHH.loci$gene %in% pc.low_and_priv, acc],na.rm = T)
  # silencing_summary_table_meth$CHH_median_priv[i]<-median(CHH.loci[CHH.loci$gene %in% private.loci.PC, acc],na.rm = T)
  
}

library(pheatmap)
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.CG.12acc.pdf",sep=""),height =3,width = 3)
pheatmap(silencing_summary_table_meth[,2:5],cluster_rows = F,cluster_cols =F,main="median CG methylation",scale = "row",labels_col = c("27","25-26","4-24","1-3"))
dev.off()

pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.CHG.12acc.pdf",sep=""),height =3,width = 3)
pheatmap(silencing_summary_table_meth[,10:13],cluster_rows = F,cluster_cols =F,main="median CHG methylation",scale = "row",labels_col = c("27","25-26","4-24","1-3"))
dev.off()
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.CHH.12acc.pdf",sep=""),height =3,width = 3)
pheatmap(silencing_summary_table_meth[,6:9],cluster_rows = F,cluster_cols =F,main="median CHH methylation",scale = "row",display_numbers = F,labels_col = c("27","25-26","4-24","1-3"))
dev.off()


#24nt sRNA

silencing_summary_table_24nt<-accessions_miRNA
rownames(silencing_summary_table_24nt)<-silencing_summary_table_24nt$V1
silencing_summary_table_24nt$sRNA_median_core<-0
silencing_summary_table_24nt$sRNA_median_softcore<-0
silencing_summary_table_24nt$sRNA_median_midfreq<-0
silencing_summary_table_24nt$sRNA_median_lowfreq<-0
#silencing_summary_table_24nt$sRNA_median_priv<-0

for ( i in 1:14){
  acc=as.character(accessions_miRNA$V1[i])
  
  silencing_summary_table_24nt$sRNA_median_core[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% pc.core, acc],na.rm = T)
  silencing_summary_table_24nt$sRNA_median_softcore[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% pc.highfreq, acc],na.rm = T)
  silencing_summary_table_24nt$sRNA_median_midfreq[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% pc.midfreq, acc],na.rm = T)
  silencing_summary_table_24nt$sRNA_median_lowfreq[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% pc.low_and_priv, acc],na.rm = T)
  #  silencing_summary_table_24nt$sRNA_median_priv[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% private.loci.PC, acc],na.rm = T)
  
}

library(pheatmap)
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.24nt.14acc.pdf",sep=""),height =3,width = 3)
pheatmap(silencing_summary_table_24nt[,2:5],cluster_rows = F,cluster_cols =F,main="24nt sRNA coverage, flowers",labels_col = c("27","25-26","4-24","1-3"))
#pheatmap(silencing_summary_table_24nt[,2:5],cluster_rows = F,cluster_cols =F,main="24nt sRNA coverage, flowers",scale = "row")
dev.off()




silencing_summary_table_24nt<-accessions_miRNA
rownames(silencing_summary_table_24nt)<-silencing_summary_table_24nt$V1
silencing_summary_table_24nt$sRNA_median_core<-0
silencing_summary_table_24nt$sRNA_median_softcore<-0
silencing_summary_table_24nt$sRNA_median_midfreq<-0
silencing_summary_table_24nt$sRNA_median_lowfreq<-0
#silencing_summary_table_24nt$sRNA_median_priv<-0

for ( i in 1:14){
  acc=as.character(accessions_miRNA$V1[i])
  
  silencing_summary_table_24nt$sRNA_median_core[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_PC_TAIR10, acc],na.rm = T)
  silencing_summary_table_24nt$sRNA_median_softcore[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_PC_new, acc],na.rm = T)
  silencing_summary_table_24nt$sRNA_median_midfreq[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_TE_TAIR10, acc],na.rm = T)
  silencing_summary_table_24nt$sRNA_median_lowfreq[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% genes_TE_new, acc],na.rm = T)
  #  silencing_summary_table_24nt$sRNA_median_priv[i]<-median(sRNA.24nt.coverage.loci[sRNA.24nt.coverage.loci$gene %in% private.loci.PC, acc],na.rm = T)
  
}

library(pheatmap)
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.24nt.14acc.oldnew.pdf",sep=""),height =3,width = 3)
pheatmap(silencing_summary_table_24nt[,2:5],cluster_rows = F,cluster_cols =F,main="24nt sRNA coverage, flowers",labels_col = c("PC TAIR10","PC new","TE TAIR10","TE new"),scale = "row")
#pheatmap(silencing_summary_table_24nt[,2:5],cluster_rows = F,cluster_cols =F,main="24nt sRNA coverage, flowers",scale = "row")
dev.off()

#expression 
silencing_summary_table_TPM<-accessions_w_220011
rownames(silencing_summary_table_TPM)<-silencing_summary_table_TPM$V1
silencing_summary_table_TPM$TPM_rosette_median_core<-0
silencing_summary_table_TPM$TPM_rosette_median_softcore<-0
silencing_summary_table_TPM$TPM_rosette_median_midfreq<-0
silencing_summary_table_TPM$TPM_rosette_median_lowfreq<-0
#silencing_summary_table_TPM$TPM_rosette_median_priv<-0
for ( i in 1:27){
  acc=as.character(accessions_w_220011$V1[i])
  silencing_summary_table_TPM$TPM_rosette_median_core[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% pc.core,paste("R",acc,sep=".")],na.rm = T)
  silencing_summary_table_TPM$TPM_rosette_median_softcore[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% pc.highfreq,paste("R",acc,sep=".")],na.rm = T)
  silencing_summary_table_TPM$TPM_rosette_median_midfreq[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% pc.midfreq,paste("R",acc,sep=".")],na.rm = T)
  silencing_summary_table_TPM$TPM_rosette_median_lowfreq[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% pc.low_and_priv,paste("R",acc,sep=".")],na.rm = T)
  #  silencing_summary_table_TPM$TPM_rosette_median_priv[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% private.loci.PC,paste("R",acc,sep=".")],na.rm = T)
  
}

library(pheatmap)
#pheatmap(silencing_summary_table_TPM[,2:5],cluster_rows = F,cluster_cols =F,main="expression in rosette",scale = "row")
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.TPM.27acc.pdf",sep=""),height =4,width = 3)
pheatmap(silencing_summary_table_TPM[,2:5],cluster_rows = F,cluster_cols =F,main="expression in rosette (TPM)",labels_col = c("27","25-26","4-24","1-3"))
dev.off()




#expression 
silencing_summary_table_TPM_anc<-accessions_w_220011
rownames(silencing_summary_table_TPM_anc)<-silencing_summary_table_TPM_anc$V1
silencing_summary_table_TPM_anc$anc<-0
silencing_summary_table_TPM_anc$seq<-0
silencing_summary_table_TPM_anc$no<-0

for ( i in 1:27){
  acc=as.character(accessions_w_220011$V1[i])
  silencing_summary_table_TPM_anc$anc[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_ancestral_PC,paste("R",acc,sep=".")],na.rm = T)
  silencing_summary_table_TPM_anc$seq[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_lyrata_seqsim_PC,paste("R",acc,sep=".")],na.rm = T)
  silencing_summary_table_TPM_anc$no[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_non_ancestral_PC,paste("R",acc,sep=".")],na.rm = T)
  
}

library(pheatmap)
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.TPM_PC_ancest_cat.27acc.pdf",sep=""),height =4,width = 3)
pheatmap(silencing_summary_table_TPM_anc[,2:4],cluster_rows = F,cluster_cols =F,main="expression in rosette (TPM)",labels_col = c("full","seq","no"))
dev.off()



#expression 
silencing_summary_table_TPM<-accessions_w_220011
rownames(silencing_summary_table_TPM)<-silencing_summary_table_TPM$V1
silencing_summary_table_TPM$TPM_rosette_median_pcold<-0
silencing_summary_table_TPM$TPM_rosette_median_pcnew<-0
silencing_summary_table_TPM$TPM_rosette_median_teold<-0
silencing_summary_table_TPM$TPM_rosette_median_tenew<-0
#silencing_summary_table_TPM$TPM_rosette_median_priv<-0
for ( i in 1:27){
  acc=as.character(accessions_w_220011$V1[i])
  silencing_summary_table_TPM$TPM_rosette_median_pcold[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_PC_TAIR10,paste("R",acc,sep=".")],na.rm = T)
  silencing_summary_table_TPM$TPM_rosette_median_pcnew[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_PC_new,paste("R",acc,sep=".")],na.rm = T)
  silencing_summary_table_TPM$TPM_rosette_median_teold[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_TE_TAIR10,paste("R",acc,sep=".")],na.rm = T)
  silencing_summary_table_TPM$TPM_rosette_median_tenew[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% genes_TE_new,paste("R",acc,sep=".")],na.rm = T)
  #  silencing_summary_table_TPM$TPM_rosette_median_priv[i]<-median(TPMs.by_full_loci[TPMs.by_full_loci$gene %in% private.loci.PC,paste("R",acc,sep=".")],na.rm = T)
  
}

library(pheatmap)
#pheatmap(silencing_summary_table_TPM[,2:5],cluster_rows = F,cluster_cols =F,main="expression in rosette",scale = "row")
pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/heatmap.TPM.27acc.new_old.pdf",sep=""),height =4,width = 3)
pheatmap(silencing_summary_table_TPM[,2:5],cluster_rows = F,cluster_cols =F,main="expression in rosette (TPM)",labels_col = c("PC TAIR10","PC new","TE TAIR10","TE new"),scale="row")
dev.off()



############################

###



# Supplemental Figure 10 E and F - new genes barplots  
###############################################
# make figure for expression on gene ontology categories - expression in 6909 - should not be expressed

# look at universal genes calculation in TAIR10 - we found them in tair10 - why were they not found before? 

# of the novel genes - how many genes have no expression in 6909? 
# how many of the genes (or novel genes) are not expressed in any accession/sample
# how many genes are the rest (expressed in 6909)

# tair10 only genes

#presence/absence of expression 
gene_table<-gene_classes[gene_classes$group %in% c(tegenes,pcgenes),1:2]
gene_table<-gene_table[gene_table$group %in% gene_frequency$group[gene_frequency$freq_loc>0],]
gene_table$expressed_in_6909<-NA
gene_table$expressed_nowhere<-NA
gene_table$anyreads_in_6909<-NA
gene_table$noreads_anywhere<-NA
gene_table$nolocus_in_6909<-NA
gene_table$N_acc_where_expressed.S<-0
gene_table$N_acc_where_expressed.R<-0
gene_table$N_acc_where_expressed.F<-0
gene_table$N_acc_where_expressed.P<-0

for (i in 1:length(gene_table$group)){
        gene<-as.character(gene_table$group[i])
    #    print(gene)
        gene_table$expressed_in_6909[i]<-apply(TPMs.by_full_loci[TPMs.by_full_loci$gene == gene,c("R.6909", "S.6909",    "F.6909",  "P.6909")] ,1, max, na.rm=T)>0.25
        gene_table$expressed_nowhere[i]<-apply(TPMs.by_full_loci[TPMs.by_full_loci$gene == gene,2:105] ,1, max, na.rm=T)<=0.25
        
        gene_table$anyreads_in_6909[i]<-apply(TPMs.by_full_loci[TPMs.by_full_loci$gene == gene,c("R.6909", "S.6909",    "F.6909",  "P.6909")] ,1, max, na.rm=T)>0
        gene_table$noreads_anywhere[i]<-apply(TPMs.by_full_loci[TPMs.by_full_loci$gene == gene,2:105] ,1, max, na.rm=T)<=0
        gene_table$nolocus_in_6909[i]<-is.na(TPMs.by_full_loci$R.6909[TPMs.by_full_loci$gene == gene])
        
}

for (i in 1:length(gene_table$group)){
        gene<-as.character(gene_table$group[i])
        #print(gene)
        gene_table$N_acc_where_expressed.F[i]<-TPMs.by_full_loci$Nacc_expressed.flowers[TPMs.by_full_loci$gene == gene]
        gene_table$N_acc_where_expressed.R[i]<-TPMs.by_full_loci$Nacc_expressed.rosette[TPMs.by_full_loci$gene == gene]
        gene_table$N_acc_where_expressed.P[i]<-TPMs.by_full_loci$Nacc_expressed.pollen[TPMs.by_full_loci$gene == gene]
        gene_table$N_acc_where_expressed.S[i]<-TPMs.by_full_loci$Nacc_expressed.seedl[TPMs.by_full_loci$gene == gene]
}

write.table(gene_table,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/expression_and_presence_stats.txt", sep="\t", quote = F,col.names = T,row.names = F)


# gather gene numbers for plotting 

supl10E_table<-as.data.frame(c("PC_TAIR10","PC_new","TE_TAIR10","TE_new"))
supl10E_table$genetype<-supl10E_table[,1]
supl10E_table$absent_in_Col0<-0
supl10E_table$expressed_in_Col0<-0
supl10E_table$expressed_in_nonCol0<-0
supl10E_table$silent_everywhere<-0

a<-gene_table[gene_table$group %in% genes_PC_TAIR10, ]
i=1
supl10E_table$absent_in_Col0[i]<-length(a$group[a$nolocus_in_6909==T])
supl10E_table$expressed_in_Col0[i]<-length(a$group[a$expressed_in_6909==T & a$nolocus_in_6909==F])
supl10E_table$expressed_in_nonCol0[i]<-length(a$group[a$expressed_in_6909==F & a$expressed_nowhere==F  & a$nolocus_in_6909==F])
supl10E_table$silent_everywhere[i]<-length(a$group[a$expressed_nowhere==T  & a$nolocus_in_6909==F])

a<-gene_table[gene_table$group %in% genes_PC_new, ]
i=2
supl10E_table$absent_in_Col0[i]<-length(a$group[a$nolocus_in_6909==T])
supl10E_table$expressed_in_Col0[i]<-length(a$group[a$expressed_in_6909==T & a$nolocus_in_6909==F])
supl10E_table$expressed_in_nonCol0[i]<-length(a$group[a$expressed_in_6909==F & a$expressed_nowhere==F & a$nolocus_in_6909==F])
supl10E_table$silent_everywhere[i]<-length(a$group[a$expressed_nowhere==T & a$nolocus_in_6909==F])

a<-gene_table[gene_table$group %in% genes_TE_TAIR10, ]
i=3
supl10E_table$absent_in_Col0[i]<-length(a$group[a$nolocus_in_6909==T])
supl10E_table$expressed_in_Col0[i]<-length(a$group[a$expressed_in_6909==T & a$nolocus_in_6909==F])
supl10E_table$expressed_in_nonCol0[i]<-length(a$group[a$expressed_in_6909==F & a$expressed_nowhere==F & a$nolocus_in_6909==F])
supl10E_table$silent_everywhere[i]<-length(a$group[a$expressed_nowhere==T & a$nolocus_in_6909==F])

a<-gene_table[gene_table$group %in% genes_TE_new, ]
i=4
supl10E_table$absent_in_Col0[i]<-length(a$group[a$nolocus_in_6909==T])
supl10E_table$expressed_in_Col0[i]<-length(a$group[a$expressed_in_6909==T & a$nolocus_in_6909==F])
supl10E_table$expressed_in_nonCol0[i]<-length(a$group[a$expressed_in_6909==F & a$expressed_nowhere==F & a$nolocus_in_6909==F])
supl10E_table$silent_everywhere[i]<-length(a$group[a$expressed_nowhere==T & a$nolocus_in_6909==F])

write.table(supl10E_table,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/expression_and_presence_old_new_gene_numbers.txt", sep="\t", quote = F,col.names = T,row.names = F)

##################
#  new PC genes 4 categories 
###############################################
#presence/absence of expression 
gene_table<-gene_classes[gene_classes$group %in% genes_PC_new,1:2]
gene_table<-gene_table[gene_table$group %in% gene_frequency$group[gene_frequency$freq_loc>0],]
gene_table$expressed_in_6909<-NA
gene_table$expressed_nowhere<-NA
gene_table$anyreads_in_6909<-NA
gene_table$noreads_anywhere<-NA
gene_table$nolocus_in_6909<-NA
gene_table$N_acc_where_expressed.S<-0
gene_table$N_acc_where_expressed.R<-0
gene_table$N_acc_where_expressed.F<-0
gene_table$N_acc_where_expressed.P<-0

for (i in 1:length(gene_table$group)){
  gene<-as.character(gene_table$group[i])
  #    print(gene)
  gene_table$expressed_in_6909[i]<-apply(TPMs.by_full_loci[TPMs.by_full_loci$gene == gene,c("R.6909", "S.6909",    "F.6909",  "P.6909")] ,1, max, na.rm=T)>0.25
  gene_table$expressed_nowhere[i]<-apply(TPMs.by_full_loci[TPMs.by_full_loci$gene == gene,2:105] ,1, max, na.rm=T)<=0.25
  
  gene_table$anyreads_in_6909[i]<-apply(TPMs.by_full_loci[TPMs.by_full_loci$gene == gene,c("R.6909", "S.6909",    "F.6909",  "P.6909")] ,1, max, na.rm=T)>0
  gene_table$noreads_anywhere[i]<-apply(TPMs.by_full_loci[TPMs.by_full_loci$gene == gene,2:105] ,1, max, na.rm=T)<=0
  gene_table$nolocus_in_6909[i]<-is.na(TPMs.by_full_loci$R.6909[TPMs.by_full_loci$gene == gene])
  
}

for (i in 1:length(gene_table$group)){
  gene<-as.character(gene_table$group[i])
  #print(gene)
  gene_table$N_acc_where_expressed.F[i]<-TPMs.by_full_loci$Nacc_expressed.flowers[TPMs.by_full_loci$gene == gene]
  gene_table$N_acc_where_expressed.R[i]<-TPMs.by_full_loci$Nacc_expressed.rosette[TPMs.by_full_loci$gene == gene]
  gene_table$N_acc_where_expressed.P[i]<-TPMs.by_full_loci$Nacc_expressed.pollen[TPMs.by_full_loci$gene == gene]
  gene_table$N_acc_where_expressed.S[i]<-TPMs.by_full_loci$Nacc_expressed.seedl[TPMs.by_full_loci$gene == gene]
}

write.table(gene_table,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/expression_and_presence_stats.PCnew_4categ.txt", sep="\t", quote = F,col.names = T,row.names = F)


# gather gene numbers for plotting 

supl10E_table<-as.data.frame(c("highsim","midsim","brass","nohom"))
supl10E_table$genetype<-supl10E_table[,1]
supl10E_table$absent_in_Col0<-0
supl10E_table$expressed_in_Col0<-0
supl10E_table$expressed_in_nonCol0<-0
supl10E_table$silent_everywhere<-0

a<-gene_table[gene_table$group %in% new_tair_high, ]
i=1
supl10E_table$absent_in_Col0[i]<-length(a$group[a$nolocus_in_6909==T])
supl10E_table$expressed_in_Col0[i]<-length(a$group[a$expressed_in_6909==T & a$nolocus_in_6909==F])
supl10E_table$expressed_in_nonCol0[i]<-length(a$group[a$expressed_in_6909==F & a$expressed_nowhere==F  & a$nolocus_in_6909==F])
supl10E_table$silent_everywhere[i]<-length(a$group[a$expressed_nowhere==T  & a$nolocus_in_6909==F])

a<-gene_table[gene_table$group %in% new_tairP_medium, ]
i=2
supl10E_table$absent_in_Col0[i]<-length(a$group[a$nolocus_in_6909==T])
supl10E_table$expressed_in_Col0[i]<-length(a$group[a$expressed_in_6909==T & a$nolocus_in_6909==F])
supl10E_table$expressed_in_nonCol0[i]<-length(a$group[a$expressed_in_6909==F & a$expressed_nowhere==F & a$nolocus_in_6909==F])
supl10E_table$silent_everywhere[i]<-length(a$group[a$expressed_nowhere==T & a$nolocus_in_6909==F])

a<-gene_table[gene_table$group %in% new_brass, ]
i=3
supl10E_table$absent_in_Col0[i]<-length(a$group[a$nolocus_in_6909==T])
supl10E_table$expressed_in_Col0[i]<-length(a$group[a$expressed_in_6909==T & a$nolocus_in_6909==F])
supl10E_table$expressed_in_nonCol0[i]<-length(a$group[a$expressed_in_6909==F & a$expressed_nowhere==F & a$nolocus_in_6909==F])
supl10E_table$silent_everywhere[i]<-length(a$group[a$expressed_nowhere==T & a$nolocus_in_6909==F])

a<-gene_table[gene_table$group %in% new_nohomol, ]
i=4
supl10E_table$absent_in_Col0[i]<-length(a$group[a$nolocus_in_6909==T])
supl10E_table$expressed_in_Col0[i]<-length(a$group[a$expressed_in_6909==T & a$nolocus_in_6909==F])
supl10E_table$expressed_in_nonCol0[i]<-length(a$group[a$expressed_in_6909==F & a$expressed_nowhere==F & a$nolocus_in_6909==F])
supl10E_table$silent_everywhere[i]<-length(a$group[a$expressed_nowhere==T & a$nolocus_in_6909==F])

write.table(supl10E_table,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/expression_and_presence_old_new_gene_numbers.PCnew_4categ.txt", sep="\t", quote = F,col.names = T,row.names = F)







################################################## 



################################################## 
##### Supplemental Figure 10 G - L 
##### SILENCING of new genes 
##################################################

#define colors
colorMap_new=c("pc_known"="#5b9bd5",
               "tair10_highsim" = "#00b050",
               "tair10_medsim" = "#a6d96a",
               "Brassicaceae" = "#b8e2b0",
               "partial_or_no_homology" = "lightgrey",
               "TE_protein_similarity" = "#cc91df",
               "TE_high_nucl_similarity" = "#9966ff")

colorMap_new1=c("PC_TAIR10"="#5b9bd5",
                "PC_new" = "#ffc000",
                "TE_TAIR10" = "#a049c7",
                "TE_new" = "#cc91df")
"TE_new_nuc" = "#9966ff"


### Expression in 27 accessions 



#############################################################
# Silencing of old vs new genes 
############################################################




# TPM boxplots - old-new genes
#########################################
a<-TPMs.by_exons

for (i in 1:27){
  acc=as.character(accessions_w_220011$V1[i])
  pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/expression/boxplot_TPM.",acc,".old_new.pdf",sep=""),height =3,width = 2.7)
  par(mar=c(8,5,3,2),mgp=c(3,1,0)) 
  a1<-a[a$gene %in% PCgenes_tair10,paste("K9",acc,sep=".")]
  a2<-a[a$gene %in% PCgenes_new_tair10_high_sim,paste("K9",acc,sep=".")]
  a3<-a[a$gene %in% PCgenes_new_tair10_med_sim,paste("K9",acc,sep=".")]
  a31<-a[a$gene %in% PCgenes_new_Brass,paste("K9",acc,sep=".")]
  a32<-a[a$gene %in% PCgenes_new_low_or_nohomology,paste("K9",acc,sep=".")]
  a4<-a[a$gene %in% PCgenes_new_TElike,paste("K9",acc,sep=".")]
  a5<-a[a$gene %in% TEgenes_tair10,paste("K9",acc,sep=".")]
  a6<-a[a$gene %in% TEgenes_new,paste("K9",acc,sep=".")]
  boxplot(a1,a2,a3,a31,a32,a4,a5,a6,
          outline=F, notch=T, main=paste("H3K9me2",acc),names=c("PC_TAIR10","new_tair10_high_sim","new_tair10_med_sim","new_Brass","new_low_or_nohomology","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="H3K9me2", col=c("#c4af31","#92de67","#42de67","#a686a0","#542e59","#6e4773"),cex.lab=0.7,cex.names=0.8,cex.main=0.8)
  aa<-c(a3,a31,a32)
  boxplot(a1,a2,aa,a4,a5,a6,
          outline=F, notch=T, main=paste("H3K9me2",acc),names=c("PC_TAIR10","new_tair10_high_sim","new","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="H3K9me2", col=c("#c4af31","#92de67","#42de67","#a686a0","#542e59","#6e4773"),cex.lab=0.7,cex.names=0.8,cex.main=0.8)
  
  
  vioplot(a1,a2,a3,a31,a32,a4,a5,a6,names=c("PC_TAIR10","new_tair10_high_sim","new_tair10_med_sim","new_Brass","new_low_or_nohomology","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="CG meth")
  #################
  #add p values   #
  #################
  a<-flexible_wilcox(a1,a2)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=0.2,cex=0.8)
  a<-flexible_wilcox(a3,a2)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=0.2,cex=0.8)
  a<-flexible_wilcox(a3,a4)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=0.2,cex=0.8)
  a<-flexible_wilcox(a4,a5)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=4.5,y=0.2,cex=0.8)
  a<-wilcox.test(a6,a5)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=5.5,y=0.2,cex=0.8)
  ####################################
  dev.off()
  
}
#########################################

# H3K9me2 boxplots - old-new genes
#########################################
a<-chipseq.loci.quantstan

for (i in 1:6){
        acc=as.character(accessions_chip$V1[i])
        pdf(file=paste("figures/expression/boxplot_K9.",acc,".old_new_3newcat.pdf",sep=""),height =3,width = 2.7)
        par(mar=c(8,5,3,2),mgp=c(3,1,0)) 
        a1<-a[a$gene %in% PCgenes_tair10,paste("K9",acc,sep=".")]
        a2<-a[a$gene %in% PCgenes_new_tair10_high_sim,paste("K9",acc,sep=".")]
        a3<-a[a$gene %in% PCgenes_new_tair10_med_sim,paste("K9",acc,sep=".")]
        a31<-a[a$gene %in% PCgenes_new_Brass,paste("K9",acc,sep=".")]
        a32<-a[a$gene %in% PCgenes_new_low_or_nohomology,paste("K9",acc,sep=".")]
        a4<-a[a$gene %in% PCgenes_new_TElike,paste("K9",acc,sep=".")]
        a5<-a[a$gene %in% TEgenes_tair10,paste("K9",acc,sep=".")]
        a6<-a[a$gene %in% TEgenes_new,paste("K9",acc,sep=".")]
        boxplot(a1,a2,a3,a31,a32,a4,a5,a6,
                outline=F, notch=T, main=paste("H3K9me2",acc),names=c("PC_TAIR10","new_tair10_high_sim","new_tair10_med_sim","new_Brass","new_low_or_nohomology","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="H3K9me2", col=c("#c4af31","#92de67","#42de67","#a686a0","#542e59","#6e4773"),cex.lab=0.7,cex.names=0.8,cex.main=0.8)
        aa<-c(a3,a31,a32)
        boxplot(a1,a2,aa,a4,a5,a6,
                outline=F, notch=T, main=paste("H3K9me2",acc),names=c("PC_TAIR10","new_tair10_high_sim","new","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="H3K9me2", col=c("#c4af31","#92de67","#42de67","#a686a0","#542e59","#6e4773"),cex.lab=0.7,cex.names=0.8,cex.main=0.8)
        
        
        vioplot(a1,a2,a3,a31,a32,a4,a5,a6,names=c("PC_TAIR10","new_tair10_high_sim","new_tair10_med_sim","new_Brass","new_low_or_nohomology","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="CG meth")
        #################
        #add p values   #
        #################
        a<-flexible_wilcox(a1,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=1.5,y=0.2,cex=0.8)
        a<-flexible_wilcox(a3,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=2.5,y=0.2,cex=0.8)
        a<-flexible_wilcox(a3,a4)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=3.5,y=0.2,cex=0.8)
        a<-flexible_wilcox(a4,a5)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=4.5,y=0.2,cex=0.8)
        a<-wilcox.test(a6,a5)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=5.5,y=0.2,cex=0.8)
        ####################################
        dev.off()
        
}
#######################################

#CG methylation - old-new genes
#########################################
a<-CG.loci

for (i in 1:12){
        acc=as.character(accessions_methylation$V1[i])
        pdf(file=paste("figures/expression/boxplot_24nt.",acc,".old_new_3newcat.pdf",sep=""),height =3,width = 2.7)
        par(mar=c(6,5,3,2),mgp=c(3,1,0)) 
        a1<-a[a$gene %in% PCgenes_tair10,acc]
        a2<-a[a$gene %in% PCgenes_new_protein_tair10,acc]
        a3<-a[a$gene %in% PCgenes_new_nucl_tair10,acc]
        a4<-a[a$gene %in% PCgenes_new_TElike,acc]
        a5<-a[a$gene %in% TEgenes_tair10,acc]
        a6<-a[a$gene %in% TEgenes_new,acc]
        boxplot(a1,a2,a3,a4,a5,a6,
                outline=F, notch=T, main=paste("CG methylation",acc),names=c("PC_TAIR10","protein similarity","nucleotide similarity","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="CG meth level", col=c("#c4af31","#92de67","#42de67","#a686a0","#542e59","#6e4773"),cex.lab=0.7,cex.names=0.8,cex.main=0.8)
        
        vioplot(a1,a2,a3,a4,a5,a6,names=c("PC_TAIR10","protein similarity","nucleotide similarity","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="CG meth level")
        
        
        #################
        #add p values   #
        #################
        a<-wilcox.test(a1,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=1.5,y=0.2,cex=0.8)
        a<-wilcox.test(a3,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=2.5,y=0.2,cex=0.8)
        a<-wilcox.test(a3,a4)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=3.5,y=0.2,cex=0.8)
        a<-wilcox.test(a4,a5)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=4.5,y=0.2,cex=0.8)
        a<-wilcox.test(a6,a5)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=5.5,y=0.2,cex=0.8)
        ####################################
        dev.off()
        
}
#########################################

#small RNAs ( 24nt TE targeting) - old-new genes
###################################
for (i in 1:14){
  acc=as.character(accessions_miRNA$V1[i])
  pdf(file=paste("figures/miRNA/boxplot_24nt.",acc,".old_new_3newcat.pdf",sep=""),height =3,width = 2.7)
  par(mar=c(6,5,3,2),mgp=c(3,1,0)) 
  a<-get(paste("a",acc,".loci.miRNA",sep=""))
  a1<-a[a$gene %in% PCgenes_tair10,acc]
  a2<-a[a$gene %in% PCgenes_new_proteinsimilarity,acc]
  a3<-a[a$gene %in% PCgenes_new_nucl_tair10,acc]
  a4<-a[a$gene %in% PCgenes_new_TElike,acc]
  a5<-a[a$gene %in% TEgenes_tair10,acc]
  a6<-a[a$gene %in% TEgenes_new,acc]
  boxplot(a1,a2,a3,a4,a5,a6,
          outline=F, notch=T, main=paste("24nt sRNA coverage\nflowers",acc),names=c("PC_TAIR10","protein_similarity","nucleotide_similarity","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="24nt coverage, RPM", col=c("#c4af31","#92de67","#42de67","#a686a0","#542e59","#6e4773"),cex.lab=0.7,cex.names=0.8,cex.main=0.8)
  #################
  #add p values   #
  #################
  a<-wilcox.test(a1,a2)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=0.2,cex=0.8)
  a<-wilcox.test(a3,a2)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=0.2,cex=0.8)
  a<-wilcox.test(a3,a4)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=0.2,cex=0.8)
  a<-wilcox.test(a4,a5)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=4.5,y=0.2,cex=0.8)
  a<-wilcox.test(a6,a5)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=5.5,y=0.2,cex=0.8)
  ####################################
  dev.off()
  
}
#########################################








# 6 new PC subclasses 


#24nt small RNAs
#########################################
for (i in 1:14){
        acc=as.character(accessions_miRNA$V1[i])
        pdf(file=paste("figures/miRNA/boxplot_24nt.",acc,".newgenes_6classes.pdf",sep=""),height =3,width = 2.7)
        par(mar=c(6,5,3,2),mgp=c(3,1,0)) 
        a<-get(paste("a",acc,".loci.miRNA",sep=""))
        a1<-a[a$gene %in% PCgenes_new_protein_tair10,acc]
        a2<-a[a$gene %in% PCgenes_new_protein_thaliana,acc]
        a3<-a[a$gene %in% PCgenes_new_nucl_tair10,acc]
        a4<-a[a$gene %in% PCgenes_new_partial_or_no_homology,acc]
        a5<-a[a$gene %in% PCgenes_new_other_species,acc]
        a6<-a[a$gene %in% PCgenes_new_TElike,acc]
        
        boxplot(a1,a2,a3,a4,a5,a6,
                outline=F, notch=T, main=paste("24nt sRNA coverage\nflowers",acc),names=c("protein_tair10","protein_thaliana","nucl_tair10","partial_or_no_homology","other_species","TE-like"),las=2,ylab="24nt coverage, RPM", col=c("#e5cf36","#e5cf36","#e5cf36","#e5cf36","#e5cf36","#7e30c1"),cex.lab=0.7,cex.names=0.8,cex.main=0.8)
        #################
        #add p values   #
        #################
        a<-wilcox.test(a1,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=1.5,y=0.1,cex=0.8)
        a<-wilcox.test(a3,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=2.5,y=0.1,cex=0.8)
        a<-wilcox.test(a3,a4)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=3.5,y=0.1,cex=0.8)
        a<-wilcox.test(a4,a5)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=4.5,y=0.1,cex=0.8)
        a<-wilcox.test(a6,a5)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=5.5,y=0.1,cex=0.8)
        ####################################
        dev.off()
        
}

#########################################

#K9
#########################################
for (i in 1:6){
        acc=as.character(accessions_chip$V1[i])
        pdf(file=paste("figures/chip/boxplot_K9.",acc,".newgenes_TE_nonTE.pdf",sep=""),height =3,width = 2.7)
        par(mar=c(6,5,3,2),mgp=c(3,1,0)) 
        a<-chipseq.loci.quantstan
        a1<-a[a$gene %in% PCgenes_tair10,paste("K9.",acc,sep="")]
        a2<-a[a$gene %in% PCgenes_new_nonTElike,paste("K9.",acc,sep="")]
        a3<-a[a$gene %in% PCgenes_new_TElike,paste("K9.",acc,sep="")]
        a4<-a[a$gene %in% TEgenes_tair10,paste("K9.",acc,sep="")]
        a5<-a[a$gene %in% TEgenes_new,paste("K9.",acc,sep="")]
        boxplot(a1,a2,a3,a4,a5,
                outline=F, notch=T, main=paste("H3K9me2 level\nleaves",acc),names=c("PC_TAIR10","PC_new","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="normalized ChIP-seq signal", col=c("#f5cf36","#e5cf36","#7e30c1","#7e32a1","#5e32a1",cex.lab=0.7,cex.names=0.8,cex.main=0.8))
        #################
        #add p values   #
        #################
        a<-wilcox.test(a1,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=1.5,y=3,cex=0.8)
        a<-wilcox.test(a3,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=2.5,y=3,cex=0.8)
        a<-wilcox.test(a3,a4)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=3.5,y=3,cex=0.8)
        a<-wilcox.test(a4,a5)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=4.5,y=3,cex=0.8)
        ####################################
        dev.off()
        
}

# 6 new PC subclasses 
#K9
#########################################
for (i in 1:6){
        acc=as.character(accessions_chip$V1[i])
        pdf(file=paste("figures/chip/boxplot_K9.",acc,".newgenes_6classes.pdf",sep=""),height =3,width = 2.7)
        par(mar=c(6,5,3,2),mgp=c(3,1,0)) 
        a<-chipseq.loci.quantstan
        a1<-a[a$gene %in% PCgenes_new_protein_tair10,paste("K9.",acc,sep="")]
        a2<-a[a$gene %in% PCgenes_new_protein_thaliana,paste("K9.",acc,sep="")]
        a3<-a[a$gene %in% PCgenes_new_nucl_tair10,paste("K9.",acc,sep="")]
        a4<-a[a$gene %in% PCgenes_new_partial_or_no_homology,paste("K9.",acc,sep="")]
        a5<-a[a$gene %in% PCgenes_new_other_species,paste("K9.",acc,sep="")]
        a6<-a[a$gene %in% PCgenes_new_TElike,paste("K9.",acc,sep="")]
        
        boxplot(a1,a2,a3,a4,a5,a6,
                outline=F, notch=T, main=paste("24nt sRNA coverage\nflowers",acc),names=c("protein_tair10","protein_thaliana","nucl_tair10","partial_or_no_homology","other_species","TE-like"),las=2,ylab="24nt coverage, RPM", col=c("#e5cf36","#e5cf36","#e5cf36","#e5cf36","#e5cf36","#7e30c1"),cex.lab=0.7,cex.names=0.8,cex.main=0.8)
        #################
        #add p values   #
        #################
        a<-wilcox.test(a1,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=1.5,y=3,cex=0.8)
        a<-wilcox.test(a3,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=2.5,y=3,cex=0.8)
        a<-wilcox.test(a3,a4)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=3.5,y=3,cex=0.8)
        a<-wilcox.test(a4,a5)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=4.5,y=3,cex=0.8)
        a<-wilcox.test(a6,a5)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=5.5,y=3,cex=0.8)
        ####################################
        dev.off()
        
}
#########################################

#K27
#########################################
for (i in 1:6){
        acc=as.character(accessions_chip$V1[i])
        pdf(file=paste("figures/chip/boxplot_K27.",acc,".newgenes_TE_nonTE.pdf",sep=""),height =3,width = 2.7)
        par(mar=c(6,5,3,2),mgp=c(3,1,0)) 
        a<-chipseq.loci.quantstan
        a1<-a[a$gene %in% PCgenes_tair10,paste("K27.",acc,sep="")]
        a2<-a[a$gene %in% PCgenes_new_nonTElike,paste("K27.",acc,sep="")]
        a3<-a[a$gene %in% PCgenes_new_TElike,paste("K27.",acc,sep="")]
        a4<-a[a$gene %in% TEgenes_tair10,paste("K27.",acc,sep="")]
        a5<-a[a$gene %in% TEgenes_new,paste("K27.",acc,sep="")]
        boxplot(a1,a2,a3,a4,a5,
                outline=F, notch=T, main=paste("H3K27me3 level\nleaves",acc),names=c("PC_TAIR10","PC_new","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="normalized ChIP-seq signal", col=c("#f5cf36","#e5cf36","#7e30c1","#7e32a1","#5e32a1",cex.lab=0.7,cex.names=0.8,cex.main=0.8))
        #################
        #add p values   #
        #################
        a<-wilcox.test(a1,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=1.5,y=2,cex=0.8)
        a<-wilcox.test(a3,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=2.5,y=2,cex=0.8)
        a<-wilcox.test(a3,a4)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=3.5,y=2,cex=0.8)
        a<-wilcox.test(a4,a5)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=4.5,y=2,cex=0.8)
        ####################################
        dev.off()
        
}

########################################### 

#K27
#########################################
for (i in 1:6){
        acc=as.character(accessions_chip$V1[i])
        pdf(file=paste("figures/chip/boxplot_K27.",acc,".newgenes_6classes.pdf",sep=""),height =3,width = 2.7)
        par(mar=c(6,5,3,2),mgp=c(3,1,0)) 
        a<-chipseq.loci.quantstan
        a1<-a[a$gene %in% PCgenes_new_protein_tair10,paste("K27.",acc,sep="")]
        a2<-a[a$gene %in% PCgenes_new_protein_thaliana,paste("K27.",acc,sep="")]
        a3<-a[a$gene %in% PCgenes_new_nucl_tair10,paste("K27.",acc,sep="")]
        a4<-a[a$gene %in% PCgenes_new_partial_or_no_homology,paste("K27.",acc,sep="")]
        a5<-a[a$gene %in% PCgenes_new_other_species,paste("K27.",acc,sep="")]
        a6<-a[a$gene %in% PCgenes_new_TElike,paste("K27.",acc,sep="")]
        
        boxplot(a1,a2,a3,a4,a5,a6,
                outline=F, notch=T, main=paste("H3K27me3\nleaves",acc),names=c("protein_tair10","protein_thaliana","nucl_tair10","partial_or_no_homology","other_species","TE-like"),las=2,ylab="normalized ChIP-seq signal", col=c("#e5cf36","#e5cf36","#e5cf36","#e5cf36","#e5cf36","#7e30c1"),cex.lab=0.7,cex.names=0.8,cex.main=0.8)
        #################
        #add p values   #
        #################
        a<-wilcox.test(a1,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=1.5,y=2,cex=0.8)
        a<-wilcox.test(a3,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=2.5,y=2,cex=0.8)
        a<-wilcox.test(a3,a4)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=3.5,y=2,cex=0.8)
        a<-wilcox.test(a4,a5)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=4.5,y=2,cex=0.8)
        a<-wilcox.test(a6,a5)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=5.5,y=2,cex=0.8)
        ####################################
        dev.off()
        
}
#########################################

#K36
#########################################
for (i in 1:6){
        acc=as.character(accessions_chip$V1[i])
        pdf(file=paste("figures/chip/boxplot_K36.",acc,".newgenes_TE_nonTE.pdf",sep=""),height =3,width = 2.7)
        par(mar=c(6,5,3,2),mgp=c(3,1,0)) 
        a<-chipseq.loci.quantstan
        a1<-a[a$gene %in% PCgenes_tair10,paste("K36.",acc,sep="")]
        a2<-a[a$gene %in% PCgenes_new_nonTElike,paste("K36.",acc,sep="")]
        a3<-a[a$gene %in% PCgenes_new_TElike,paste("K36.",acc,sep="")]
        a4<-a[a$gene %in% TEgenes_tair10,paste("K36.",acc,sep="")]
        a5<-a[a$gene %in% TEgenes_new,paste("K36.",acc,sep="")]
        boxplot(a1,a2,a3,a4,a5,
                outline=F, notch=T, main=paste("H3K36me3 level\nleaves",acc),names=c("PC_TAIR10","PC_new","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="normalized ChIP-seq signal", col=c("#f5cf36","#e5cf36","#7e30c1","#7e32a1","#5e32a1",cex.lab=0.7,cex.names=0.8,cex.main=0.8))
        #################
        #add p values   #
        #################
        a<-wilcox.test(a1,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=1.5,y=1.5,cex=0.8)
        a<-wilcox.test(a3,a2)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=2.5,y=1.5,cex=0.8)
        a<-wilcox.test(a3,a4)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=3.5,y=1.5,cex=0.8)
        a<-wilcox.test(a4,a5)
        d<-a$p.value
        if (d<0.0000000001){b="***"} 
        if (d>=0.0000000001 &d<0.00001 ){b="**"} 
        if (d>=0.00001  &d<0.01){b="*"} 
        if (d>=0.01){b="n.s."} 
        text(b,x=4.5,y=1.5,cex=0.8)
        ####################################
        dev.off()
        
}
#########################################


#CG methylation - newPC categories 
#########################################
a<-CG.loci

for (i in 1:12){
  acc=as.character(accessions_methylation$V1[i])
  pdf(file=paste("figures/expression/boxplot_24nt.",acc,".old_new_3newcat.pdf",sep=""),height =3,width = 2.7)
  par(mar=c(8,5,3,2),mgp=c(3,1,0)) 
  a1<-a[a$gene %in% PCgenes_tair10,acc]
  a2<-a[a$gene %in% PCgenes_new_tair10_high_sim,acc]
  a3<-a[a$gene %in% PCgenes_new_tair10_med_sim,acc]
  a31<-a[a$gene %in% PCgenes_new_Brass,acc]
  a32<-a[a$gene %in% PCgenes_new_low_or_nohomology,acc]
  a4<-a[a$gene %in% PCgenes_new_TElike,acc]
  a5<-a[a$gene %in% TEgenes_tair10,acc]
  a6<-a[a$gene %in% TEgenes_new,acc]
  boxplot(a1,a2,a3,a31,a32,a4,a5,a6,
          outline=F, notch=T, main=paste("CG methylation",acc),names=c("PC_TAIR10","new_tair10_high_sim","new_tair10_med_sim","new_Brass","new_low_or_nohomology","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="CG meth level", col=c("#c4af31","#92de67","#42de67","#a686a0","#542e59","#6e4773"),cex.lab=0.7,cex.names=0.8,cex.main=0.8)
  
  vioplot(a1,a2,a3,a31,a32,a4,a5,a6,names=c("PC_TAIR10","new_tair10_high_sim","new_tair10_med_sim","new_Brass","new_low_or_nohomology","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="CG meth")
  #################
  #add p values   #
  #################
  a<-flexible_wilcox(a1,a2)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=1.5,y=0.2,cex=0.8)
  a<-flexible_wilcox(a3,a2)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=2.5,y=0.2,cex=0.8)
  a<-flexible_wilcox(a3,a4)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=3.5,y=0.2,cex=0.8)
  a<-flexible_wilcox(a4,a5)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=4.5,y=0.2,cex=0.8)
  a<-flexible_wilcox(a6,a5)
  d<-a$p.value
  if (d<0.0000000001){b="***"} 
  if (d>=0.0000000001 &d<0.00001 ){b="**"} 
  if (d>=0.00001  &d<0.01){b="*"} 
  if (d>=0.01){b="n.s."} 
  text(b,x=5.5,y=0.2,cex=0.8)
  ####################################
  dev.off()
  
}
#########################################


#new_PCgene categories 
###################################
a<-TPMs.by_full_loci

for (i in 1:27){
  acc=as.character(accessions$V1[i])
  pdf(file=paste("figures/expression/TPM_by_full_loci.",acc,".old_new_3newcat.pdf",sep=""),height =3,width = 2.7)
  acc=as.character(accessions$V1[i])
  a1<-a[a$gene %in% PCgenes_tair10,paste("R.",acc,sep="")]
  a2<-a[a$gene %in% PCgenes_new_proteinsimilarity,paste("R.",acc,sep="")]
  a3<-a[a$gene %in% PCgenes_new_nucl_tair10,paste("R.",acc,sep="")]
  a4<-a[a$gene %in% PCgenes_new_TElike,paste("R.",acc,sep="")]
  a5<-a[a$gene %in% TEgenes_tair10,paste("R.",acc,sep="")]
  a6<-a[a$gene %in% TEgenes_new,paste("R.",acc,sep="")]
  boxplot(a1,a2,a3,a4,a5,a6,
          outline=F, notch=T, main=paste("24nt sRNA coverage\nflowers",acc),names=c("PC_TAIR10","protein_similarity","nucleotide_similarity","TE-like_new","TE_TAIR10","TE_new"),las=2,ylab="24nt coverage, RPM", col=c("#c4af31","#92de67","#42de67","#a686a0","#542e59","#6e4773"),ylim=c(0,5),cex.lab=0.7,cex.names=0.8,cex.main=0.8)
  
  
  







###########################################################

  
  
  
  
  
###########################################################
#########################FIGURE8 - analysis#####
###########################################################

############################ 
# how many genes are wrongly estimated 
# expression correlation - between TPMs calculated on tair10 and own genomes
###########################
correlations_byexons<-data.frame(matrix(nrow = 2,ncol = 2))
names(correlations_byexons)<-c("sample","TPM_correlation")
correlations_byloci<-data.frame(1:2,1:2)
names(correlations_byloci)<-c("sample","TPM_correlation")

a<-counts.by_exons.unprocessed
a<-merge(a,gene_classes[,c("gene","AraportID")],by="gene")
b<-counts.Araport_exons.unprocessed
a<-merge(a,b,by.x="AraportID",by.y="gene")
aa<-cor(a[,grep(".x",names(a))],a[,grep(".y",names(a))],use="complete.obs")
correlations_byexons<-as.data.frame(cbind(row.names(aa),diag(aa)))
names(correlations_byexons)<-c("sample","TPM_correlation")
  
a<-counts.by_full_loci.unprocessed
a<-merge(a,gene_classes[,c("gene","AraportID")],by="gene")
b<-counts.Araport_locus.unprocessed
a<-merge(a,b,by.x="AraportID",by.y="gene")
aa<-cor(a[,grep(".x",names(a))],a[,grep(".y",names(a))],use="complete.obs")
correlations_byloci<-as.data.frame(cbind(row.names(aa),diag(aa)))
names(correlations_byloci)<-c("sample","TPM_correlation")
###############################################


#####################################################
# gene-gene match through annotation 
######################################################

## estimation of wrong gene number own genome vs TAIR10 

# expression genes - through counts
############################################
table_with_wrong_and_good_genes_ownvsTAIR<-eracaps_125samples[,2:3]
table_with_wrong_and_good_genes_ownvsTAIR$cor<-0
table_with_wrong_and_good_genes_ownvsTAIR$N_goodgenes<-0
table_with_wrong_and_good_genes_ownvsTAIR$N_wronggenes<-0

table_with_wrong_and_good_genes_ownvsTAIR$N_wronggenes_PC<-0
table_with_wrong_and_good_genes_ownvsTAIR$N_wronggenes_TE<-0
table_with_wrong_and_good_genes_ownvsTAIR$N_goodgenes_PC<-0
table_with_wrong_and_good_genes_ownvsTAIR$N_goodgenes_TE<-0

table_with_wrong_and_good_genes_ownvsTAIR$N_wronggenes_CNV<-0
table_with_wrong_and_good_genes_ownvsTAIR$N_wronggenes_nonCNV<-0 
table_with_wrong_and_good_genes_ownvsTAIR$N_goodgenes_CNV<-0
table_with_wrong_and_good_genes_ownvsTAIR$N_goodgenes_nonCNV<-0


library(scales)
for (i in 1:125)
{
sample=as.character(eracaps_125samples$sample[i])
acc=as.character(eracaps_125samples$acc[i])
# on own genome 
a<-counts.by_exons.unprocessed[,c("gene",sample)]
a<-merge(a,gene_classes[,c("gene","AraportID")],by="gene")
#length(gene_classes$group[!is.na(gene_classes$AraportID)])
#[1] 27371
#on tair10
b<-counts.Araport_exons.unprocessed[,c("gene",sample)]
a<-merge(a,b,by.x="AraportID",by.y="gene")
a<-merge(a,copyN[,c("group",paste("X",acc,sep=""),"X6909")],by.x="gene",by.y="group")
table_with_wrong_and_good_genes_ownvsTAIR$cor[i]<-cor(a[,3],a[,4],use="complete.obs")
  
names(a)<-c("gene","Araport_ID","on_own","on_tair10","CN_own","CN_6909")
a$max<-apply(a[,3:4],1,max)
a$min<-apply(a[,3:4],1,min)
a$var<-apply(a[,3:4],1,sd)/apply(a[,3:4],1,mean)
a$ratio<-apply(a[,3:4],1,min)/apply(a[,3:4],1,max)
  
#good genes: <30% mistake in TAIR10
#mistakenly estimated genes >30% mistake in TAIR10
#only look at genes with at least 6 reads
wronggenes<-a$gene[a$max>=6 & a$ratio<=0.70 &!is.na(a$max) & !is.na(a$ratio)]
goodgenes<-a$gene[a$max>=6 & a$ratio>0.70 &!is.na(a$max) & !is.na(a$ratio)]


pdf(file=paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/expression_scatter/",sample,".pdf",sep=""),height =4,width = 4)
par(mar=c(6,5,3,2),mgp=c(3,1,0)) 
plot(a$on_own+1,a$on_tair10+1,log="xy",main=paste("exonic read counts\n",sample), xlab="on own genome,log(x+1)",ylab="on TAIR10,log(y+1)",col = alpha("white", 0.2), pch=16)
points(a$on_own[a$max<6 ]+1,a$on_tair10[a$max<6]+1,col = alpha("darkgrey", 0.3), pch=16)
points(a$on_own[a$max>=6 & a$ratio>0.70]+1,a$on_tair10[a$max>=6 & a$ratio>0.70]+1,col = alpha("darkgreen", 0.2), pch=16)
points(a$on_own[a$max>=6 & a$ratio<=0.70]+1,a$on_tair10[a$max>=6 & a$ratio<=0.70]+1,col = alpha("darkred", 0.2), pch=16)
dev.off()


table_with_wrong_and_good_genes_ownvsTAIR$N_goodgenes[i]<-length(goodgenes)
table_with_wrong_and_good_genes_ownvsTAIR$N_wronggenes[i]<-length(wronggenes)


table_with_wrong_and_good_genes_ownvsTAIR$N_wronggenes_PC[i]<-length(wronggenes[wronggenes %in% pcgenes]) 
table_with_wrong_and_good_genes_ownvsTAIR$N_wronggenes_TE[i]<-length(wronggenes[wronggenes %in% tegenes]) 

table_with_wrong_and_good_genes_ownvsTAIR$N_goodgenes_PC[i]<-length(goodgenes[goodgenes %in% pcgenes]) 
table_with_wrong_and_good_genes_ownvsTAIR$N_goodgenes_TE[i]<-length(goodgenes[goodgenes %in% tegenes]) 

table_with_wrong_and_good_genes_ownvsTAIR$N_wronggenes_CNV[i]<-length(wronggenes[wronggenes %in% copyN$group[copyN$sd_CN>0]]) 
table_with_wrong_and_good_genes_ownvsTAIR$N_wronggenes_nonCNV[i]<-length(wronggenes[wronggenes %in% copyN$group[copyN$sd_CN==0]]) 

table_with_wrong_and_good_genes_ownvsTAIR$N_goodgenes_CNV[i]<-length(goodgenes[goodgenes %in% copyN$group[copyN$sd_CN>0]])
table_with_wrong_and_good_genes_ownvsTAIR$N_goodgenes_nonCNV[i]<-length(goodgenes[goodgenes %in% copyN$group[copyN$sd_CN==0]]) 
}

############################################


write.table(table_with_wrong_and_good_genes_ownvsTAIR,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/Expression_wrong_and_good_genes_ownvsTAIR.bed",quote = F,col.names = T,row.names = F,sep="\t")

#CG 
#########################################################
methyl_ownvsTAIR.CG<-accessions_methylation
methyl_ownvsTAIR.CG$acc<-methyl_ownvsTAIR.CG$V1
methyl_ownvsTAIR.CG$cor<-0
methyl_ownvsTAIR.CG$N_goodgenes<-0
methyl_ownvsTAIR.CG$N_wronggenes<-0

methyl_ownvsTAIR.CG$N_wronggenes_PC<-0
methyl_ownvsTAIR.CG$N_wronggenes_TE<-0
methyl_ownvsTAIR.CG$N_goodgenes_PC<-0
methyl_ownvsTAIR.CG$N_goodgenes_TE<-0

methyl_ownvsTAIR.CG$N_wronggenes_CNV<-0
methyl_ownvsTAIR.CG$N_wronggenes_nonCNV<-0 
methyl_ownvsTAIR.CG$N_goodgenes_CNV<-0
methyl_ownvsTAIR.CG$N_goodgenes_nonCNV<-0

for (i in 1:12)
{ 
  
  acc=as.character(methyl_ownvsTAIR.CG$acc[i])
  
  # on own genome 
  a<-CG.loci
  a<-merge(a[,c("gene",acc)],gene_classes[,c("gene","AraportID")],by="gene")
  #on tair10
  b<-get(paste("a",acc,".TAIR10.loci.CG",sep = ""))
  a<-merge(a,b[,c("gene",acc)],by.x="AraportID",by.y="gene")
  a<-merge(a,copyN[,c("group",paste("X",acc,sep=""),"X6909")],by.y="group", by.x="gene")
  methyl_ownvsTAIR.CG$cor[i]<-cor(a[,3],a[,4],use = "complete.obs")
  names(a)<-c("gene","Araport_ID","on_own","on_tair10","CN_own","CN_6909")
  a$max<-apply(a[,3:4],1,max)
  a$min<-apply(a[,3:4],1,min)
  a$var<-apply(a[,3:4],1,sd)/apply(a[,3:4],1,mean)
  a$percentchange<-100*(apply(a[,3:4],1,max)-apply(a[,3:4],1,min))/apply(a[,3:4],1,max)
  #good genes: change in methylation level <50%
  #mistakenly estimated genes change in methylation level >50%
  wronggenes<-a$gene[a$percentchange>50 & !is.na(a$percentchange)& !is.na(a$on_own) & !is.na(a$on_tair10)]
  goodgenes<-a$gene[a$percentchange<=50& !is.na(a$percentchange)&!is.na(a$on_own) & !is.na(a$on_tair10)]
  
  methyl_ownvsTAIR.CG$N_goodgenes[i]<-length(goodgenes)
  methyl_ownvsTAIR.CG$N_wronggenes[i]<-length(wronggenes)
  methyl_ownvsTAIR.CG$N_wronggenes_PC[i]<-length(wronggenes[wronggenes %in% pcgenes]) 
  methyl_ownvsTAIR.CG$N_wronggenes_TE[i]<-length(wronggenes[wronggenes %in% tegenes]) 
  
  methyl_ownvsTAIR.CG$N_goodgenes_PC[i]<-length(goodgenes[goodgenes %in% pcgenes]) 
  methyl_ownvsTAIR.CG$N_goodgenes_TE[i]<-length(goodgenes[goodgenes %in% tegenes]) 
  
  methyl_ownvsTAIR.CG$N_wronggenes_CNV[i]<-length(wronggenes[wronggenes %in% copyN$group[copyN$sd_CN>0]]) 
  
  methyl_ownvsTAIR.CG$N_wronggenes_nonCNV[i]<-length(wronggenes[wronggenes %in% copyN$group[copyN$sd_CN==0]]) 
  
  methyl_ownvsTAIR.CG$N_goodgenes_CNV[i]<-length(goodgenes[goodgenes %in% copyN$group[copyN$sd_CN>0]])
  methyl_ownvsTAIR.CG$N_goodgenes_nonCNV[i]<-length(goodgenes[goodgenes %in% copyN$group[copyN$sd_CN==0]]) 
}
#########################################################

#CHG 
#########################################################
methyl_ownvsTAIR.CHG<-accessions_methylation
methyl_ownvsTAIR.CHG$acc<-methyl_ownvsTAIR.CHG$V1
methyl_ownvsTAIR.CHG$cor<-0
methyl_ownvsTAIR.CHG$N_goodgenes<-0
methyl_ownvsTAIR.CHG$N_wronggenes<-0

methyl_ownvsTAIR.CHG$N_wronggenes_PC<-0
methyl_ownvsTAIR.CHG$N_wronggenes_TE<-0
methyl_ownvsTAIR.CHG$N_goodgenes_PC<-0
methyl_ownvsTAIR.CHG$N_goodgenes_TE<-0

methyl_ownvsTAIR.CHG$N_wronggenes_CNV<-0
methyl_ownvsTAIR.CHG$N_wronggenes_nonCNV<-0 
methyl_ownvsTAIR.CHG$N_goodgenes_CNV<-0
methyl_ownvsTAIR.CHG$N_goodgenes_nonCNV<-0

for (i in 1:12)
{ 
  
  acc=as.character(methyl_ownvsTAIR.CHG$acc[i])
  
  # on own genome 
  a<-CHG.loci
  a<-merge(a[,c("gene",acc)],gene_classes[,c("gene","AraportID")],by="gene")
  #on tair10
  b<-get(paste("a",acc,".TAIR10.loci.CHG",sep = ""))
  a<-merge(a,b[,c("gene",acc)],by.x="AraportID",by.y="gene")
  a<-merge(a,copyN[,c("group",paste("X",acc,sep=""),"X6909")],by.y="group", by.x="gene")
  methyl_ownvsTAIR.CHG$cor[i]<-cor(a[,3],a[,4],use = "complete.obs")
  names(a)<-c("gene","Araport_ID","on_own","on_tair10","CN_own","CN_6909")
  a$max<-apply(a[,3:4],1,max)
  a$min<-apply(a[,3:4],1,min)
  a$var<-apply(a[,3:4],1,sd)/apply(a[,3:4],1,mean)
  a$percentchange<-100*(apply(a[,3:4],1,max)-apply(a[,3:4],1,min))/apply(a[,3:4],1,max)
  #good genes: change in methylation level <50%
  #mistakenly estimated genes change in methylation level >50%
  wronggenes<-a$gene[a$percentchange>50 &!is.na(a$percentchange)& !is.na(a$on_own) & !is.na(a$on_tair10)]
  goodgenes<-a$gene[a$percentchange<=50& !is.na(a$percentchange)&!is.na(a$on_own) & !is.na(a$on_tair10)]
  
  methyl_ownvsTAIR.CHG$N_goodgenes[i]<-length(goodgenes)
  methyl_ownvsTAIR.CHG$N_wronggenes[i]<-length(wronggenes)
  methyl_ownvsTAIR.CHG$N_wronggenes_PC[i]<-length(wronggenes[wronggenes %in% pcgenes]) 
  methyl_ownvsTAIR.CHG$N_wronggenes_TE[i]<-length(wronggenes[wronggenes %in% tegenes]) 
  
  methyl_ownvsTAIR.CHG$N_goodgenes_PC[i]<-length(goodgenes[goodgenes %in% pcgenes]) 
  methyl_ownvsTAIR.CHG$N_goodgenes_TE[i]<-length(goodgenes[goodgenes %in% tegenes]) 
  
  methyl_ownvsTAIR.CHG$N_wronggenes_CNV[i]<-length(wronggenes[wronggenes %in% copyN$group[copyN$sd_CN>0]]) 
  
  methyl_ownvsTAIR.CHG$N_wronggenes_nonCNV[i]<-length(wronggenes[wronggenes %in% copyN$group[copyN$sd_CN==0]]) 
  
  methyl_ownvsTAIR.CHG$N_goodgenes_CNV[i]<-length(goodgenes[goodgenes %in% copyN$group[copyN$sd_CN>0]])
  methyl_ownvsTAIR.CHG$N_goodgenes_nonCNV[i]<-length(goodgenes[goodgenes %in% copyN$group[copyN$sd_CN==0]]) 
}
#########################################################

#CHH 
#########################################################
methyl_ownvsTAIR.CHH<-accessions_methylation
methyl_ownvsTAIR.CHH$acc<-methyl_ownvsTAIR.CHH$V1
methyl_ownvsTAIR.CHH$cor<-0
methyl_ownvsTAIR.CHH$N_goodgenes<-0
methyl_ownvsTAIR.CHH$N_wronggenes<-0

methyl_ownvsTAIR.CHH$N_wronggenes_PC<-0
methyl_ownvsTAIR.CHH$N_wronggenes_TE<-0
methyl_ownvsTAIR.CHH$N_goodgenes_PC<-0
methyl_ownvsTAIR.CHH$N_goodgenes_TE<-0

methyl_ownvsTAIR.CHH$N_wronggenes_CNV<-0
methyl_ownvsTAIR.CHH$N_wronggenes_nonCNV<-0 
methyl_ownvsTAIR.CHH$N_goodgenes_CNV<-0
methyl_ownvsTAIR.CHH$N_goodgenes_nonCNV<-0

for (i in 1:12)
{ 
  
  acc=as.character(methyl_ownvsTAIR.CHH$acc[i])
  
  # on own genome 
  a<-CHH.loci
  a<-merge(a[,c("gene",acc)],gene_classes[,c("gene","AraportID")],by="gene")
  #on tair10
  b<-get(paste("a",acc,".TAIR10.loci.CHH",sep = ""))
  a<-merge(a,b[,c("gene",acc)],by.x="AraportID",by.y="gene")
  a<-merge(a,copyN[,c("group",paste("X",acc,sep=""),"X6909")],by.y="group", by.x="gene")
  methyl_ownvsTAIR.CHH$cor[i]<-cor(a[,3],a[,4],use = "complete.obs")
  names(a)<-c("gene","Araport_ID","on_own","on_tair10","CN_own","CN_6909")
  a$max<-apply(a[,3:4],1,max)
  a$min<-apply(a[,3:4],1,min)
  a$var<-apply(a[,3:4],1,sd)/apply(a[,3:4],1,mean)
  a$percentchange<-100*(apply(a[,3:4],1,max)-apply(a[,3:4],1,min))/apply(a[,3:4],1,max)
  #good genes: change in methylation level <50%
  #mistakenly estimated genes change in methylation level >50%
  wronggenes<-a$gene[a$percentchange>50 & !is.na(a$percentchange)&!is.na(a$on_own) & !is.na(a$on_tair10)]
  goodgenes<-a$gene[a$percentchange<=50& !is.na(a$percentchange)&!is.na(a$on_own) & !is.na(a$on_tair10)]
  
  methyl_ownvsTAIR.CHH$N_goodgenes[i]<-length(goodgenes)
  methyl_ownvsTAIR.CHH$N_wronggenes[i]<-length(wronggenes)
  methyl_ownvsTAIR.CHH$N_wronggenes_PC[i]<-length(wronggenes[wronggenes %in% pcgenes]) 
  methyl_ownvsTAIR.CHH$N_wronggenes_TE[i]<-length(wronggenes[wronggenes %in% tegenes]) 
  
  methyl_ownvsTAIR.CHH$N_goodgenes_PC[i]<-length(goodgenes[goodgenes %in% pcgenes]) 
  methyl_ownvsTAIR.CHH$N_goodgenes_TE[i]<-length(goodgenes[goodgenes %in% tegenes]) 
  
  methyl_ownvsTAIR.CHH$N_wronggenes_CNV[i]<-length(wronggenes[wronggenes %in% copyN$group[copyN$sd_CN>0]]) 
  
  methyl_ownvsTAIR.CHH$N_wronggenes_nonCNV[i]<-length(wronggenes[wronggenes %in% copyN$group[copyN$sd_CN==0]]) 
  
  methyl_ownvsTAIR.CHH$N_goodgenes_CNV[i]<-length(goodgenes[goodgenes %in% copyN$group[copyN$sd_CN>0]])
  methyl_ownvsTAIR.CHH$N_goodgenes_nonCNV[i]<-length(goodgenes[goodgenes %in% copyN$group[copyN$sd_CN==0]]) 
}
#########################################################






#export the table 
write.table(methyl_ownvsTAIR.CG,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/CG.methylation_wrong_and_good_genes_ownvsTAIR.bed",quote = F,col.names = T,row.names = F,sep="\t")
write.table(methyl_ownvsTAIR.CHG,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/CHG.methylation_wrong_and_good_genes_ownvsTAIR.bed",quote = F,col.names = T,row.names = F,sep="\t")
write.table(methyl_ownvsTAIR.CHH,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/CHH.methylation_wrong_and_good_genes_ownvsTAIR.bed",quote = F,col.names = T,row.names = F,sep="\t")


#FIG 8 - plotting 
###########################################
table_with_wrong_and_good_genes_ownvsTAIR$tissue<-unlist(lapply(strsplit(as.character(table_with_wrong_and_good_genes_ownvsTAIR$sample),".", fixed = T), "[",1))


# for plotting - remove 6909 from the tables 

expr<-table_with_wrong_and_good_genes_ownvsTAIR[table_with_wrong_and_good_genes_ownvsTAIR$acc!="6909",]
cg_bias<-methyl_ownvsTAIR.CG[methyl_ownvsTAIR.CG$acc!="6909",]
chg_bias<-methyl_ownvsTAIR.CHG[methyl_ownvsTAIR.CHG$acc!="6909",]
chh_bias<-methyl_ownvsTAIR.CHH[methyl_ownvsTAIR.CHH$acc!="6909",]


#Pearson correlation boxplot
###########################################
pdf(file="Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/boxplot_correlation_onTAIR_vs_on_own.pdf",height =3.5,width = 4)
par(mar=c(6,5,3,2),mgp=c(3,1,0))

set.seed(42)
# Create a boxplot
a1<-expr$cor[expr$tissue=="S"]
a2<-expr$cor[expr$tissue=="R"]
a3<-expr$cor[expr$tissue=="F"]
a4<-expr$cor[expr$tissue=="P"]

a5<-cg_bias$cor
a6<-chg_bias$cor
a7<-chh_bias$cor


boxplot(a1,a2,a3,a4,-1,a5,a6,a7,
        col = "lightgray",
        ylab = "Pearson correlation:\n align to own vs. TAIR10",ylim=c(0.6,1.05),outline = F,notch = F, names=c("seedling","rosette","flowers","pollen","","CG","CHG","CHH"),las=2)

mtext(paste("N=",length(a1),sep = ""),side = 1,line=-1.5,cex=0.7,at = 1,las=1)
mtext(paste("N=",length(a2),sep = ""),side = 1,line=-1.5,cex=0.7,at = 2,las=1)
mtext(paste("N=",length(a3),sep = ""),side = 1,line=-1.5,cex=0.7,at = 3,las=1)
mtext(paste("N=",length(a4),sep = ""),side = 1,line=-1.5,cex=0.7,at = 4,las=1)
mtext(paste("N=",length(a5),sep = ""),side = 1,line=-1.5,cex=0.7,at = 6,las=1)
mtext(paste("N=",length(a6),sep = ""),side = 1,line=-1.5,cex=0.7,at = 7,las=1)
mtext(paste("N=",length(a7),sep = ""),side = 1,line=-1.5,cex=0.7,at = 8,las=1)

mtext("expression",side = 3,line=0,cex=1,at = 2.5)
mtext("methylation",side = 3,line=0,cex=1,at = 7)

# Overlay the data points
stripchart(a1, at = 1,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)
stripchart(a2, at = 2,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)
stripchart(a3, at = 3,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)
stripchart(a4, at = 4,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)
stripchart(a5, at = 6,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)
stripchart(a6, at = 7,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)
stripchart(a7, at = 8,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)

dev.off()

###########################################

# % wrong genes 
###########################################
pdf(file="Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/boxplot_percent_wrong_genes_onTAIR_vs_on_own.pdf",height =3.5,width = 4)
par(mar=c(6,5,3,2),mgp=c(3,1,0))

set.seed(42)
# Create a boxplot
a1<-100*expr$N_wronggenes[expr$tissue=="S"]/(expr$N_wronggenes[expr$tissue=="S"]+expr$N_goodgenes[expr$tissue=="S"])

a2<-100*expr$N_wronggenes[expr$tissue=="R"]/(expr$N_wronggenes[expr$tissue=="R"]+expr$N_goodgenes[expr$tissue=="R"])

a3<-100*expr$N_wronggenes[expr$tissue=="F"]/(expr$N_wronggenes[expr$tissue=="F"]+expr$N_goodgenes[expr$tissue=="F"])
a4<-100*expr$N_wronggenes[expr$tissue=="P"]/(expr$N_wronggenes[expr$tissue=="P"]+expr$N_goodgenes[expr$tissue=="P"])

a5<-100*cg_bias$N_wronggenes/(cg_bias$N_wronggenes+cg_bias$N_goodgenes)
a6<-100*chg_bias$N_wronggenes/(chg_bias$N_wronggenes+chg_bias$N_goodgenes)
a7<-100*chh_bias$N_wronggenes/(chh_bias$N_wronggenes+chh_bias$N_goodgenes)


boxplot(a1,a2,a3,a4,-10,a5,a6,a7,
        col = "lightgray",
        ylab = "% genes wrongly\nestimated on TAIR10",ylim=c(-4,35),outline = F,notch = F, names=c("seedling","rosette","flowers","pollen","","CG","CHG","CHH"),las=2)

mtext(paste("N=",length(a1),sep = ""),side = 1,line=-1,cex=0.7,at = 1,las=1)
mtext(paste("N=",length(a2),sep = ""),side = 1,line=-1,cex=0.7,at = 2,las=1)
mtext(paste("N=",length(a3),sep = ""),side = 1,line=-1,cex=0.7,at = 3,las=1)
mtext(paste("N=",length(a4),sep = ""),side = 1,line=-1,cex=0.7,at = 4,las=1)
mtext(paste("N=",length(a5),sep = ""),side = 1,line=-1,cex=0.7,at = 6,las=1)
mtext(paste("N=",length(a6),sep = ""),side = 1,line=-1,cex=0.7,at = 7,las=1)
mtext(paste("N=",length(a7),sep = ""),side = 1,line=-1,cex=0.7,at = 8,las=1)

mtext("expression",side = 3,line=0,cex=1,at = 2.5)
mtext("methylation",side = 3,line=0,cex=1,at = 7)

# Overlay the data points
stripchart(a1, at = 1,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)
stripchart(a2, at = 2,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)
stripchart(a3, at = 3,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)
stripchart(a4, at = 4,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)
stripchart(a5, at = 6,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)
stripchart(a6, at = 7,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)
stripchart(a7, at = 8,method = "jitter", pch = 16, col = rgb(0, 0, 0, 0.5), vertical = TRUE, add = TRUE, cex=0.7)

dev.off()
###########################################

# % CNV genes 
###########################################
pdf(file="Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/boxplot_percent_CNV_genes_onTAIR_vs_on_own.pdf",height =3.5,width = 6)
par(mar=c(6,5,3,2),mgp=c(3,1,0))

set.seed(42)
# Create a boxplot
a1<-100*expr$N_wronggenes_CNV[expr$tissue=="S"]/expr$N_wronggenes[expr$tissue=="S"]

a11<-100*expr$N_goodgenes_CNV[expr$tissue=="S"]/expr$N_goodgenes[expr$tissue=="S"]


a2<-100*expr$N_wronggenes_CNV[expr$tissue=="R"]/expr$N_wronggenes[expr$tissue=="R"]

a22<-100*expr$N_goodgenes_CNV[expr$tissue=="R"]/expr$N_goodgenes[expr$tissue=="R"]

a3<-100*expr$N_wronggenes_CNV[expr$tissue=="F"]/expr$N_wronggenes[expr$tissue=="F"]

a33<-100*expr$N_goodgenes_CNV[expr$tissue=="F"]/expr$N_goodgenes[expr$tissue=="F"]

a4<-100*expr$N_wronggenes_CNV[expr$tissue=="P"]/expr$N_wronggenes[expr$tissue=="P"]

a44<-100*expr$N_goodgenes_CNV[expr$tissue=="P"]/expr$N_goodgenes[expr$tissue=="P"]

a5<-100*cg_bias$N_wronggenes_CNV/cg_bias$N_wronggenes
a55<-100*cg_bias$N_goodgenes_CNV/cg_bias$N_goodgenes

a6<-100*chg_bias$N_wronggenes_CNV/chg_bias$N_wronggenes
a66<-100*chg_bias$N_goodgenes_CNV/chg_bias$N_goodgenes

a7<-100*chh_bias$N_wronggenes_CNV/chh_bias$N_wronggenes
a77<-100*chh_bias$N_goodgenes_CNV/chh_bias$N_goodgenes


boxplot(a1,a11,a2,a22,a3,a33,a4,a44,-10,a5,a55,a6,a66,a7,a77,
        ylab = "% CNV genes",ylim=c(-2,40),outline = F,notch = F, las=2, col=c("brown1","darkolivegreen3","brown1","darkolivegreen3","brown1","darkolivegreen3","brown1","darkolivegreen3","black","brown1","darkolivegreen3","brown1","darkolivegreen3","brown1","darkolivegreen3"))

#names=c("seedling","rosette","flowers","pollen","","CG","CHG","CHH")
#mtext(paste("N=",length(a1),sep = ""),side = 1,line=-1,cex=0.7,at = 1,las=1)
#mtext(paste("N=",length(a2),sep = ""),side = 1,line=-1,cex=0.7,at = 2,las=1)
#mtext(paste("N=",length(a3),sep = ""),side = 1,line=-1,cex=0.7,at = 3,las=1)
#mtext(paste("N=",length(a4),sep = ""),side = 1,line=-1,cex=0.7,at = 4,las=1)
#mtext(paste("N=",length(a5),sep = ""),side = 1,line=-1,cex=0.7,at = 6,las=1)
#mtext(paste("N=",length(a6),sep = ""),side = 1,line=-1,cex=0.7,at = 7,las=1)
#mtext(paste("N=",length(a7),sep = ""),side = 1,line=-1,cex=0.7,at = 8,las=1)

mtext("expression",side = 3,line=0,cex=1,at = 4.5)
mtext("methylation",side = 3,line=0,cex=1,at = 12.5)

# Overlay the data points
stripchart(a1, at = 1,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a11, at = 2,method = "jitter", pch = 16, col = rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)

stripchart(a2, at = 3,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a22, at = 4,method = "jitter", pch = 16, col =rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)

stripchart(a3, at = 5,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a33, at = 6,method = "jitter", pch = 16, col =rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)

stripchart(a4, at = 7,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a44, at = 8,method = "jitter", pch = 16, col =rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)

stripchart(a5, at = 10,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a55, at = 11,method = "jitter", pch = 16, col =rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)

stripchart(a6, at = 12,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a66, at = 13,method = "jitter", pch = 16, col =rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)

stripchart(a7, at = 14,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a77, at = 15,method = "jitter", pch = 16, col =rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)



#pvalues

wt<-wilcox.test(a1,a11)#p-value = 1.027e-15
wt$p.value
wt<-wilcox.test(a2,a22)# 1.027253e-15
wt$p.value
wt<-wilcox.test(a3,a33)#2.901654e-22
wt$p.value
wt<-wilcox.test(a4,a44) #2.614827e-16
wt$p.value
wt<-wilcox.test(a5,a55)#2.835142e-06
wt$p.value
wt<-wilcox.test(a6,a66)#2.835142e-06
wt$p.value
wt<-wilcox.test(a77,a7) #2.835142e-06
wt$p.value
dev.off()
###########################################

# % TE genes 
###########################################
pdf(file="Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/figures/boxplot_percent_TE_genes_onTAIR_vs_on_own.pdf",height =3.5,width = 6)
par(mar=c(6,5,3,2),mgp=c(3,1,0))

set.seed(42)
# Create a boxplot
a1<-100*expr$N_wronggenes_TE[expr$tissue=="S"]/expr$N_wronggenes[expr$tissue=="S"]

a11<-100*expr$N_goodgenes_TE[expr$tissue=="S"]/expr$N_goodgenes[expr$tissue=="S"]


a2<-100*expr$N_wronggenes_TE[expr$tissue=="R"]/expr$N_wronggenes[expr$tissue=="R"]

a22<-100*expr$N_goodgenes_TE[expr$tissue=="R"]/expr$N_goodgenes[expr$tissue=="R"]

a3<-100*expr$N_wronggenes_TE[expr$tissue=="F"]/expr$N_wronggenes[expr$tissue=="F"]

a33<-100*expr$N_goodgenes_TE[expr$tissue=="F"]/expr$N_goodgenes[expr$tissue=="F"]

a4<-100*expr$N_wronggenes_TE[expr$tissue=="P"]/expr$N_wronggenes[expr$tissue=="P"]

a44<-100*expr$N_goodgenes_TE[expr$tissue=="P"]/expr$N_goodgenes[expr$tissue=="P"]

a5<-100*cg_bias$N_wronggenes_TE/cg_bias$N_wronggenes
a55<-100*cg_bias$N_goodgenes_TE/cg_bias$N_goodgenes

a6<-100*chg_bias$N_wronggenes_TE/chg_bias$N_wronggenes
a66<-100*chg_bias$N_goodgenes_TE/chg_bias$N_goodgenes

a7<-100*chh_bias$N_wronggenes_TE/chh_bias$N_wronggenes
a77<-100*chh_bias$N_goodgenes_TE/chh_bias$N_goodgenes


boxplot(a1,a11,a2,a22,a3,a33,a4,a44,-10,a5,a55,a6,a66,a7,a77,
        ylab = "% TE genes",ylim=c(-0.6,6),outline = F,notch = F, las=2, col=c("brown1","darkolivegreen3","brown1","darkolivegreen3","brown1","darkolivegreen3","brown1","darkolivegreen3","black","brown1","darkolivegreen3","brown1","darkolivegreen3","brown1","darkolivegreen3"))

#names=c("seedling","rosette","flowers","pollen","","CG","CHG","CHH")
#mtext(paste("N=",length(a1),sep = ""),side = 1,line=-1,cex=0.7,at = 1,las=1)
#mtext(paste("N=",length(a2),sep = ""),side = 1,line=-1,cex=0.7,at = 2,las=1)
#mtext(paste("N=",length(a3),sep = ""),side = 1,line=-1,cex=0.7,at = 3,las=1)
#mtext(paste("N=",length(a4),sep = ""),side = 1,line=-1,cex=0.7,at = 4,las=1)
#mtext(paste("N=",length(a5),sep = ""),side = 1,line=-1,cex=0.7,at = 6,las=1)
#mtext(paste("N=",length(a6),sep = ""),side = 1,line=-1,cex=0.7,at = 7,las=1)
#mtext(paste("N=",length(a7),sep = ""),side = 1,line=-1,cex=0.7,at = 8,las=1)

mtext("expression",side = 3,line=0,cex=1,at = 4.5)
mtext("methylation",side = 3,line=0,cex=1,at = 12.5)

# Overlay the data points
stripchart(a1, at = 1,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a11, at = 2,method = "jitter", pch = 16, col = rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)

stripchart(a2, at = 3,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a22, at = 4,method = "jitter", pch = 16, col =rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)

stripchart(a3, at = 5,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a33, at = 6,method = "jitter", pch = 16, col =rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)

stripchart(a4, at = 7,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a44, at = 8,method = "jitter", pch = 16, col =rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)

stripchart(a5, at = 10,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a55, at = 11,method = "jitter", pch = 16, col =rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)

stripchart(a6, at = 12,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a66, at = 13,method = "jitter", pch = 16, col =rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)

stripchart(a7, at = 14,method = "jitter", pch = 16, col = rgb(0.4, 0, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)
stripchart(a77, at = 15,method = "jitter", pch = 16, col =rgb(0, 0.4, 0, 0.5), vertical = TRUE, add = TRUE,cex=0.7)


dev.off()
###########################################

###########################################








# Define a flexible Wilcoxon test function
##########################
flexible_wilcox <- function(a1, a2) {
  # Ensure numeric
  a1 <- as.numeric(a1)
  a2 <- as.numeric(a2)
  
  # Remove NA and infinite values
  a1 <- a1[!is.na(a1) & !is.infinite(a1)]
  a2 <- a2[!is.na(a2) & !is.infinite(a2)]
  
  # Check if exact calculation is possible (no ties)
  use_exact <- (length(unique(a1)) == length(a1) && length(unique(a2)) == length(a2))
  
  # Perform Wilcoxon test
  wilcox.test(a1, a2, exact = use_exact)
}
##########################


