#setwd("/groups/nordborg/user/aleksandra.kornienko/Documents/01_POSTDOC/03_Projects/ERA-CAPS")
#setwd("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/")
#setwd("")


eracaps_125samples<- read.table("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/eracaps_125_RNAseq_samples.txt", quote="\"", comment.char="",header = F)
eracaps_125samples$acc<-unlist(lapply(strsplit(as.character(eracaps_125samples$V1),".", fixed = T), "[",2))
eracaps_125samples$sample<-eracaps_125samples$V1


#import lists of accessions 
##############################################################
accessions_chip <- read.table("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/accessions_chip.txt", quote="\"", comment.char="")
accessions_methylation<- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/accessions_methylation1001.txt", sep="",header = F)
accessions_miRNA<- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/accessions_miRNAseq.txt", sep="",header = F)

accessions <- read.table("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/accessions_annas_order.txt", quote="\"", comment.char="")
accessions_w_220011 <- read.table("Z:/01_POSTDOC/03_Projects/ERA-CAPS/accessions_annas_order_with_220011.txt", quote="\"", comment.char="")
##############################################################

######################################################
#import gene and mRNA frequencies
######################################

freq_loc<-read.delim("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/version_X_eracaps/gff_features/counts_gene.txt")
freq_loc$freq_loc<-apply(freq_loc[2:28],1,function(i) sum(i>0,na.rm=TRUE))
freq_loc$group<-rownames(freq_loc)

freq_mrna<-read.delim("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/version_X_eracaps/gff_features/counts_mrna.txt")
freq_mrna$freq_mrna<-apply(freq_mrna[1:27],1,function(i) sum(i>0,na.rm=TRUE))
freq_mrna$group<-rownames(freq_mrna)

gene_frequency<-merge(freq_loc[,c("group","freq_loc")],freq_mrna[,c("group","freq_mrna")],by="group",all=T)
rm (freq_mrna,freq_loc)
######################################################

#import Annotation-related information 
######################################################
# gene info table from Haim
gene_classes <- read.delim("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/version_X_eracaps/All_1001g_plus.18Dec2024.AtGrp_cat.txt")
gene_classes$group<-gene_classes$AtGrp
gene_classes$AraportID<-gene_classes$TAIR10
gene_classes$te_content<-"NA"
gene_classes$te_content[gene_classes$Pannagram.is.te==0]<-"low"
gene_classes$te_content[gene_classes$Pannagram.is.te==1]<-"high"
gene_classes$te_content[is.na(gene_classes$Pannagram.is.te)]<-"zero"
gene_classes$TEgene<-"NA"
gene_classes$TEgene[gene_classes$is.TE_like==TRUE | gene_classes$te_content=="high"]<- "TRUE"
gene_classes$TEgene[(gene_classes$is.TE_like==FALSE | is.na(gene_classes$is.TE_like)) & gene_classes$te_content!="high"]<- "FALSE"

gene_classes$ancestral_cat<-NA
gene_classes$ancestral_cat[gene_classes$Aly_Ancestral==TRUE]<-"ancestral"
gene_classes$ancestral_cat[gene_classes$Aly_Ancestral==FALSE & gene_classes$Aly_Ancestral_OR_Similar==TRUE]<-"ancestral_onlyseq"
gene_classes$ancestral_cat[gene_classes$Aly_Ancestral_OR_Similar==FALSE]<-"non_ancestral"



gene_classes <- gene_classes[,c("group","AraportID","uniprot_class","is.TE_like"     ,         "TEgene","nestgr_Remove_flag","is.overlap_TAIR10","is.TAIR10_Pseudogene", "is.new" ,"ancestral_cat", "AtGrp_overlap","Pannagram.is.te")   ]

te_content <- read.delim("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/version_X_eracaps/gff_features/te_coverage_mrna_stat.txt", header=T)

gene_classes$remove_final<-NA
gene_classes$remove_final[gene_classes$group %in% filteroutgenes]<-TRUE
gene_classes$remove_final[!(gene_classes$group %in% filteroutgenes)]<-FALSE
gene_classes$gene<-gene_classes$group

write.table(gene_classes[,c("group","AraportID","is.new","remove_final")],"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/Annot_groups_categories.txt",quote = F, sep="\t",col.names = T,row.names = F)


## import categories for NEW PC genes - based on similarity search by haim 


newgenes_classes <- read.delim("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/version_X_eracaps/gene_PC.new.SimCategories.txt")
newgenes_classes$gene<-newgenes_classes$AtGrp
newgenes_classes$group<-newgenes_classes$AtGrp

Brassicaceae       Low or No homology            Other_Species 
817                      421                        9 
TAIR10_High_Similarity TAIR10_Medium_Similarity 
744                      670 

###########################################################


#Import copy number 
###########################################################

copyN_all_to_all<-read.delim("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/version_X_eracaps/copynumber.txt",header = T)
copyN_all_to_all$accessionX<-paste("X",copyN_all_to_all$accession,sep="")
copyN_all_to_all$X220011<-copyN_all_to_all$X22001_mod
copyN_all_to_all$TAIR10<-copyN_all_to_all$X0
copyN_all_to_all$group<-copyN_all_to_all$gene
copyN_all_to_all<-copyN_all_to_all[,c("group","accession","accessionX","length"  ,   "TAIR10"    ,     "X10002" ,    "X10015"    ,  "X10024","X1741","X220011" ,"X22002","X22003","X22004","X22005"    
               , "X22006","X22007","X6024","X6069","X6124","X6244","X6909"
      ,"X6966","X8236","X9075","X9537","X9543","X9638" ,"X9728"     
        ,"X9764"     , "X9888","X9905","X9981"   )]


i=1
acc<-paste("X",as.character(accessions_w_220011$V1[i]),sep="")
a<-copyN_all_to_all[copyN_all_to_all$accessionX==acc,c("group",acc)]
copyN<-a

for (i in 2:27){
  acc<-paste("X",as.character(accessions_w_220011$V1[i]),sep="")
  
  a<-copyN_all_to_all[copyN_all_to_all$accessionX==acc,c("group",acc)]
  copyN<-merge(copyN,a, by="group",all = T)
  }


copyN$max_CN<-apply(copyN[,2:28],1,max,na.rm=T)
copyN$min_CN<-apply(copyN[,2:28],1,min,na.rm=T)
copyN$mean_CN<-apply(copyN[,2:28],1,mean,na.rm=T)
copyN$freq<-apply(copyN[,2:28], 1, function(i) sum(i >=1,na.rm=T))
copyN$sd_CN<-apply(copyN[,2:28],1,sd,na.rm=T)

copyN$CNV<-FALSE
copyN$CNV[copyN$max_CN==copyN$min_CN]<-FALSE
copyN$CNV[copyN$max_CN!=copyN$min_CN]<-TRUE
copyN$enough_0_1_variation<-"non_variable"
copyN$enough_1_multi_variation<-"non_variable"

copyN$enough_0_1_variation[apply(copyN[,2:28], 1, function(i) sum(i ==0))>1 &apply(copyN[,2:28], 1, function(i) sum(i >0))>1 ]<-"variable"
copyN$enough_1_multi_variation[apply(copyN[,2:28], 1, function(i) sum(i ==1))>1 &apply(copyN[,2:28], 1, function(i) sum(i >1))>1 ]<-"variable"

copyN$any_0_1_variation<-"non_variable"
copyN$any_1_multi_variation<-"non_variable"

copyN$any_0_1_variation[apply(copyN[,2:28], 1, function(i) sum(i ==0))>=1 &apply(copyN[,2:28], 1, function(i) sum(i >0))>=1 ]<-"variable"
copyN$any_1_multi_variation[apply(copyN[,2:28], 1, function(i) sum(i ==1))>=1 &apply(copyN[,2:28], 1, function(i) sum(i >1))>=1 ]<-"variable"

#add araport IDs

Annotgr_AraportID<-gene_classes[,c("AraportID","group")]
copyN<-merge(copyN,Annotgr_AraportID,all.x = T,by="group")




write.table(copyN,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/AnnotGroups.CN.27genomes.txt",quote = F, sep="\t",col.names = T,row.names = F)

############################

#import names of genes with tandem dup 
###############################################################

for ( i in 1:27){
  acc=as.character(accessions_w_220011$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/annotations/copies/tandem/genes_with_tandem.DUP.",acc,".txt",sep=""),header = F)
  a<-as.vector(unique(a$V1))
   assign(paste("a",acc,".gen_tandup",sep=""),a)
}

genesnames_tandem_dup<-unique(c(a10002.gen_tandup, a10015.gen_tandup, a10024.gen_tandup ,
a1741.gen_tandup,   a220011.gen_tandup, a22002.gen_tandup ,
 a22003.gen_tandup  ,a22004.gen_tandup , a22005.gen_tandup ,
 a22006.gen_tandup  ,a22007.gen_tandup , a6024.gen_tandup  ,
 a6069.gen_tandup ,  a6124.gen_tandup   ,a6244.gen_tandup  ,
 a6909.gen_tandup  , a6966.gen_tandup,   a8236.gen_tandup  ,
 a9075.gen_tandup   ,a9537.gen_tandup ,  a9543.gen_tandup  ,
 a9638.gen_tandup ,  a9728.gen_tandup  , a9764.gen_tandup  ,
 a9888.gen_tandup,   a9905.gen_tandup  , a9981.gen_tandup  ))

i=1
acc=as.character(accessions$V1[i])
a<-get(paste("a",acc,".gen_tandup",sep=""))
a<-as.data.frame(a)
a$gene<-a$a
a[,acc]<-a$a
tandem.table<-a[,2:3]

for ( i in 2:27){
  acc=as.character(accessions_w_220011$V1[i])
    a<-get(paste("a",acc,".gen_tandup",sep=""))
    a<-as.data.frame(a)
    a$gene<-a$a
    a[,acc]<-a$a
    a<-a[,2:3]
    tandem.table<- merge (tandem.table,a,by="gene", all=T)
}
############################################################################


write.table(tandem.table,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/Tandemdups.AnnotGroups_that_have_a_tandem_copy.27genomes.txt",quote = F, sep="\t",col.names = T,row.names = F)

duplicated_genes<-as.vector(unique(copyN$gene[copyN$max_CN>1]))

############################################################################


######################################################
#import EXPRESSION 
######################################################

#TPMs on own genomes by exons
#on own genomes by exons (mRNAs )
##############################################
i=1
acc=as.character(accessions_w_220011$V1[i])
a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/",acc,".TPMs.by_exons.bed",sep=""),header = T)
assign(paste("a",acc,".TPMs.by_exons",sep=""),a)
TPMs.by_exons<-a

for ( i in 2:27){
  acc=as.character(accessions_w_220011$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/",acc,".TPMs.by_exons.bed",sep=""),header = T)
  assign(paste("a",acc,".TPMs.by_exons",sep=""),a)
  TPMs.by_exons<-merge(TPMs.by_exons,a,by="gene",all = T)
}
TPMs.by_exons.unprocessed<-TPMs.by_exons
# average those samples that have replicates, remove batch names from sample names
TPMs.by_exons$F.10002<-apply(TPMs.by_exons[,grep("F.10002",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$P.10002<-TPMs.by_exons$P.10002.batch2
TPMs.by_exons$F.10015<-apply(TPMs.by_exons[,grep("F.10015",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$P.10015<-apply(TPMs.by_exons[,grep("P.10015",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$S.10015<-apply(TPMs.by_exons[,grep("S.10015",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$R.10015<-apply(TPMs.by_exons[,grep("R.10015",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$F.10024<-apply(TPMs.by_exons[,grep("F.10024",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$P.10024<-TPMs.by_exons$P.10024.batch2
TPMs.by_exons$F.1741<-TPMs.by_exons$F.1741.batchCT
TPMs.by_exons$P.1741<-TPMs.by_exons$P.1741.batchCT

TPMs.by_exons$F.22007<-apply(TPMs.by_exons[,grep("F.22007",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$P.22007<-apply(TPMs.by_exons[,grep("P.22007",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$F.22006<-TPMs.by_exons$F.22006.batchCT
TPMs.by_exons$P.22006<-TPMs.by_exons$P.22006.batchCT
TPMs.by_exons$F.22005<-apply(TPMs.by_exons[,grep("F.22005",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$P.22005<-apply(TPMs.by_exons[,grep("P.22005",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$F.22004<-TPMs.by_exons$F.22004.batchCT
TPMs.by_exons$P.22004<-TPMs.by_exons$P.22004.batchCT
TPMs.by_exons$F.22003<-TPMs.by_exons$F.22003.batchCT
TPMs.by_exons$P.22003<-TPMs.by_exons$P.22003.batchCT
TPMs.by_exons$F.22002<-apply(TPMs.by_exons[,grep("F.22002",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$P.22002<-TPMs.by_exons$P.22002.batch2
TPMs.by_exons$F.220011<-apply(TPMs.by_exons[,grep("F.220011",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$P.220011<-TPMs.by_exons$P.220011.batch2

TPMs.by_exons$F.6024<-TPMs.by_exons$F.6024.batchCT
TPMs.by_exons$P.6024<-TPMs.by_exons$P.6024.batchCT
TPMs.by_exons$F.6244<-TPMs.by_exons$F.6244.batchCT
TPMs.by_exons$P.6244<-TPMs.by_exons$P.6244.batchCT
TPMs.by_exons$F.6909<-apply(TPMs.by_exons[,grep("F.6909",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$P.6909<-TPMs.by_exons$P.6909.batch2
TPMs.by_exons$F.6966<-apply(TPMs.by_exons[,grep("F.6966",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$P.6966<-TPMs.by_exons$P.6966.batch2
TPMs.by_exons$F.8236<-apply(TPMs.by_exons[,grep("F.8236",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$P.8236<-TPMs.by_exons$P.8236.batch2
TPMs.by_exons$F.9075<-TPMs.by_exons$F.9075.batch1
TPMs.by_exons$P.9075<-TPMs.by_exons$P.9075.batchCT
TPMs.by_exons$F.9537<-apply(TPMs.by_exons[,grep("F.9537",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$P.9537<-TPMs.by_exons$P.9537.batch2
TPMs.by_exons$F.9543<-TPMs.by_exons$F.9543.batchCT
TPMs.by_exons$P.9543<-TPMs.by_exons$P.9543.batchCT
TPMs.by_exons$F.9638<-TPMs.by_exons$F.9638.batchCT
TPMs.by_exons$P.9638<-TPMs.by_exons$P.9638.batchCT
TPMs.by_exons$F.9728<-TPMs.by_exons$F.9728.batchCT
TPMs.by_exons$P.9728<-TPMs.by_exons$P.9728.batchCT
TPMs.by_exons$F.9764<-apply(TPMs.by_exons[,grep("F.9764",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$P.9764<-apply(TPMs.by_exons[,grep("P.9764",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$F.9888<-TPMs.by_exons$F.9888.batchCT
TPMs.by_exons$P.9888<-TPMs.by_exons$P.9888.batchCT
TPMs.by_exons$F.9905<-TPMs.by_exons$F.9905.batchCT
TPMs.by_exons$P.9905<-TPMs.by_exons$P.9905.batchCT
TPMs.by_exons$F.9981<-apply(TPMs.by_exons[,grep("F.9981",names(TPMs.by_exons))],1,mean)
TPMs.by_exons$P.9981<-TPMs.by_exons$P.9981.batch2

TPMs.by_exons<-TPMs.by_exons[,grep("batch",names(TPMs.by_exons),invert = T)]
TPMs.by_exons<-TPMs.by_exons[,grep("rep",names(TPMs.by_exons),invert = T)]


TPMs.by_exons$max.flowers<- apply(TPMs.by_exons[,grep("F.",names(TPMs.by_exons))],1,max,na.rm=T)
TPMs.by_exons$max.pollen<- apply(TPMs.by_exons[,grep("P.",names(TPMs.by_exons))],1,max,na.rm=T)
TPMs.by_exons$max.rosette<- apply(TPMs.by_exons[,grep("R.",names(TPMs.by_exons))],1,max,na.rm=T)
TPMs.by_exons$max.seedl<- apply(TPMs.by_exons[,grep("S.",names(TPMs.by_exons))],1,max,na.rm=T)

TPMs.by_exons$mean.flowers<- apply(TPMs.by_exons[,grep("F.",names(TPMs.by_exons))],1,mean,na.rm=T)
TPMs.by_exons$mean.pollen<- apply(TPMs.by_exons[,grep("P.",names(TPMs.by_exons))],1,mean,na.rm=T)
TPMs.by_exons$mean.rosette<- apply(TPMs.by_exons[,grep("R.",names(TPMs.by_exons))],1,mean,na.rm=T)
TPMs.by_exons$mean.seedl<- apply(TPMs.by_exons[,grep("S.",names(TPMs.by_exons))],1,mean,na.rm=T)

TPMs.by_exons$Nacc_expressed.flowers<- apply(TPMs.by_exons[,grep("F.",names(TPMs.by_exons))],1,function(i) sum(i > 0.5,na.rm=T))
TPMs.by_exons$Nacc_expressed.pollen<- apply(TPMs.by_exons[,grep("P.",names(TPMs.by_exons))],1,function(i) sum(i > 0.5,na.rm=T))
TPMs.by_exons$Nacc_expressed.rosette<- apply(TPMs.by_exons[,grep("R.",names(TPMs.by_exons))],1,function(i) sum(i > 0.5,na.rm=T))
TPMs.by_exons$Nacc_expressed.seedl<- apply(TPMs.by_exons[,grep("S.",names(TPMs.by_exons))],1,function(i) sum(i > 0.5,na.rm=T))

TPMs.by_exons$freq<- apply(TPMs.by_exons[,grep("R.",names(TPMs.by_exons))],1,function(i) sum(!is.na(i)))
rownames(TPMs.by_exons)<-TPMs.by_exons$gene

#remove intermediate TPM tables for individual  accessions
rm(list = ls(pattern = "^a[0-9]+\\.TPMs\\.by_exons$"))

############################################################
write.table(TPMs.by_exons,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/TPMs.by_exons.mapped_to_own_genomes.txt",quote = F, sep="\t",col.names = T,row.names = F)

########### TPMs on own genomes by full loci 
#########################################################
#on own genomes  by full loci 
i=1
acc=as.character(accessions_w_220011$V1[i])
a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/",acc,".TPMs.by_full_loci.bed",sep=""),header = T)
assign(paste("a",acc,".TPMs.by_full_loci",sep=""),a)
TPMs.by_full_loci<-a

for ( i in 2:27){
  acc=as.character(accessions_w_220011$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/",acc,".TPMs.by_full_loci.bed",sep=""),header = T)
  assign(paste("a",acc,".TPMs.by_full_loci",sep=""),a)
  TPMs.by_full_loci<-merge(TPMs.by_full_loci,a,by="gene",all = T)
}

TPMs.by_full_loci.unprocessed<-TPMs.by_full_loci
# average those samples that have replicates, remove batch names from sample names
TPMs.by_full_loci$F.10002<-apply(TPMs.by_full_loci[,grep("F.10002",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$P.10002<-TPMs.by_full_loci$P.10002.batch2
TPMs.by_full_loci$F.10015<-apply(TPMs.by_full_loci[,grep("F.10015",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$P.10015<-apply(TPMs.by_full_loci[,grep("P.10015",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$S.10015<-apply(TPMs.by_full_loci[,grep("S.10015",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$R.10015<-apply(TPMs.by_full_loci[,grep("R.10015",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$F.10024<-apply(TPMs.by_full_loci[,grep("F.10024",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$P.10024<-TPMs.by_full_loci$P.10024.batch2
TPMs.by_full_loci$F.1741<-TPMs.by_full_loci$F.1741.batchCT
TPMs.by_full_loci$P.1741<-TPMs.by_full_loci$P.1741.batchCT

TPMs.by_full_loci$F.22007<-apply(TPMs.by_full_loci[,grep("F.22007",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$P.22007<-apply(TPMs.by_full_loci[,grep("P.22007",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$F.22006<-TPMs.by_full_loci$F.22006.batchCT
TPMs.by_full_loci$P.22006<-TPMs.by_full_loci$P.22006.batchCT
TPMs.by_full_loci$F.22005<-apply(TPMs.by_full_loci[,grep("F.22005",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$P.22005<-apply(TPMs.by_full_loci[,grep("P.22005",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$F.22004<-TPMs.by_full_loci$F.22004.batchCT
TPMs.by_full_loci$P.22004<-TPMs.by_full_loci$P.22004.batchCT
TPMs.by_full_loci$F.22003<-TPMs.by_full_loci$F.22003.batchCT
TPMs.by_full_loci$P.22003<-TPMs.by_full_loci$P.22003.batchCT
TPMs.by_full_loci$F.22002<-apply(TPMs.by_full_loci[,grep("F.22002",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$P.22002<-TPMs.by_full_loci$P.22002.batch2
TPMs.by_full_loci$F.220011<-apply(TPMs.by_full_loci[,grep("F.220011",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$P.220011<-TPMs.by_full_loci$P.220011.batch2

TPMs.by_full_loci$F.6024<-TPMs.by_full_loci$F.6024.batchCT
TPMs.by_full_loci$P.6024<-TPMs.by_full_loci$P.6024.batchCT
TPMs.by_full_loci$F.6244<-TPMs.by_full_loci$F.6244.batchCT
TPMs.by_full_loci$P.6244<-TPMs.by_full_loci$P.6244.batchCT
TPMs.by_full_loci$F.6909<-apply(TPMs.by_full_loci[,grep("F.6909",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$P.6909<-TPMs.by_full_loci$P.6909.batch2
TPMs.by_full_loci$F.6966<-apply(TPMs.by_full_loci[,grep("F.6966",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$P.6966<-TPMs.by_full_loci$P.6966.batch2
TPMs.by_full_loci$F.8236<-apply(TPMs.by_full_loci[,grep("F.8236",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$P.8236<-TPMs.by_full_loci$P.8236.batch2
TPMs.by_full_loci$F.9075<-TPMs.by_full_loci$F.9075.batch1
TPMs.by_full_loci$P.9075<-TPMs.by_full_loci$P.9075.batchCT
TPMs.by_full_loci$F.9537<-apply(TPMs.by_full_loci[,grep("F.9537",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$P.9537<-TPMs.by_full_loci$P.9537.batch2
TPMs.by_full_loci$F.9543<-TPMs.by_full_loci$F.9543.batchCT
TPMs.by_full_loci$P.9543<-TPMs.by_full_loci$P.9543.batchCT
TPMs.by_full_loci$F.9638<-TPMs.by_full_loci$F.9638.batchCT
TPMs.by_full_loci$P.9638<-TPMs.by_full_loci$P.9638.batchCT
TPMs.by_full_loci$F.9728<-TPMs.by_full_loci$F.9728.batchCT
TPMs.by_full_loci$P.9728<-TPMs.by_full_loci$P.9728.batchCT
TPMs.by_full_loci$F.9764<-apply(TPMs.by_full_loci[,grep("F.9764",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$P.9764<-apply(TPMs.by_full_loci[,grep("P.9764",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$F.9888<-TPMs.by_full_loci$F.9888.batchCT
TPMs.by_full_loci$P.9888<-TPMs.by_full_loci$P.9888.batchCT
TPMs.by_full_loci$F.9905<-TPMs.by_full_loci$F.9905.batchCT
TPMs.by_full_loci$P.9905<-TPMs.by_full_loci$P.9905.batchCT
TPMs.by_full_loci$F.9981<-apply(TPMs.by_full_loci[,grep("F.9981",names(TPMs.by_full_loci))],1,mean)
TPMs.by_full_loci$P.9981<-TPMs.by_full_loci$P.9981.batch2

TPMs.by_full_loci<-TPMs.by_full_loci[,grep("batch",names(TPMs.by_full_loci),invert = T)]
TPMs.by_full_loci<-TPMs.by_full_loci[,grep("rep",names(TPMs.by_full_loci),invert = T)]


TPMs.by_full_loci$max.flowers<- apply(TPMs.by_full_loci[,grep("F.",names(TPMs.by_full_loci))],1,max,na.rm=T)
TPMs.by_full_loci$max.pollen<- apply(TPMs.by_full_loci[,grep("P.",names(TPMs.by_full_loci))],1,max,na.rm=T)
TPMs.by_full_loci$max.rosette<- apply(TPMs.by_full_loci[,grep("R.",names(TPMs.by_full_loci))],1,max,na.rm=T)
TPMs.by_full_loci$max.seedl<- apply(TPMs.by_full_loci[,grep("S.",names(TPMs.by_full_loci))],1,max,na.rm=T)

TPMs.by_full_loci$mean.flowers<- apply(TPMs.by_full_loci[,grep("F.",names(TPMs.by_full_loci))],1,mean,na.rm=T)
TPMs.by_full_loci$mean.pollen<- apply(TPMs.by_full_loci[,grep("P.",names(TPMs.by_full_loci))],1,mean,na.rm=T)
TPMs.by_full_loci$mean.rosette<- apply(TPMs.by_full_loci[,grep("R.",names(TPMs.by_full_loci))],1,mean,na.rm=T)
TPMs.by_full_loci$mean.seedl<- apply(TPMs.by_full_loci[,grep("S.",names(TPMs.by_full_loci))],1,mean,na.rm=T)

TPMs.by_full_loci$Nacc_expressed.flowers<- apply(TPMs.by_full_loci[,grep("F.",names(TPMs.by_full_loci))],1,function(i) sum(i > 0.25,na.rm=T))
TPMs.by_full_loci$Nacc_expressed.pollen<- apply(TPMs.by_full_loci[,grep("P.",names(TPMs.by_full_loci))],1,function(i) sum(i > 0.25,na.rm=T))
TPMs.by_full_loci$Nacc_expressed.rosette<- apply(TPMs.by_full_loci[,grep("R.",names(TPMs.by_full_loci))],1,function(i) sum(i > 0.25,na.rm=T))
TPMs.by_full_loci$Nacc_expressed.seedl<- apply(TPMs.by_full_loci[,grep("S.",names(TPMs.by_full_loci))],1,function(i) sum(i > 0.25,na.rm=T))

TPMs.by_full_loci$freq<- apply(TPMs.by_full_loci[,grep("R.",names(TPMs.by_full_loci))],1,function(i) sum(!is.na(i)))


rownames(TPMs.by_full_loci)<-TPMs.by_full_loci$gene


#remove intermediate TPM tables for individual  accessions
rm(list = ls(pattern = "^a[0-9]+\\.TPMs\\.by_full_loci$"))
#########################################################

write.table(TPMs.by_full_loci,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/TPMs.by_full_loci.mapped_to_own_genomes.txt",quote = F, sep="\t",col.names = T,row.names = F)


########### counts on own genomes by full loci 
#########################################################
#on own genomes  by full loci 
i=1
acc=as.character(accessions_w_220011$V1[i])
a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/",acc,".counts.by_full_loci.bed",sep=""),header = T)
assign(paste("a",acc,".counts.by_full_loci",sep=""),a)
counts.by_full_loci<-a

for ( i in 2:27){
  acc=as.character(accessions_w_220011$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/",acc,".counts.by_full_loci.bed",sep=""),header = T)
  assign(paste("a",acc,".counts.by_full_loci",sep=""),a)
  counts.by_full_loci<-merge(counts.by_full_loci,a,by="gene",all = T)
}
counts.by_full_loci.unprocessed<-counts.by_full_loci

# average those samples that have replicates, remove batch names from sample names
counts.by_full_loci$F.10002<-apply(counts.by_full_loci[,grep("F.10002",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$P.10002<-counts.by_full_loci$P.10002.batch2
counts.by_full_loci$F.10015<-apply(counts.by_full_loci[,grep("F.10015",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$P.10015<-apply(counts.by_full_loci[,grep("P.10015",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$S.10015<-apply(counts.by_full_loci[,grep("S.10015",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$R.10015<-apply(counts.by_full_loci[,grep("R.10015",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$F.10024<-apply(counts.by_full_loci[,grep("F.10024",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$P.10024<-counts.by_full_loci$P.10024.batch2
counts.by_full_loci$F.1741<-counts.by_full_loci$F.1741.batchCT
counts.by_full_loci$P.1741<-counts.by_full_loci$P.1741.batchCT

counts.by_full_loci$F.22007<-apply(counts.by_full_loci[,grep("F.22007",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$P.22007<-apply(counts.by_full_loci[,grep("P.22007",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$F.22006<-counts.by_full_loci$F.22006.batchCT
counts.by_full_loci$P.22006<-counts.by_full_loci$P.22006.batchCT
counts.by_full_loci$F.22005<-apply(counts.by_full_loci[,grep("F.22005",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$P.22005<-apply(counts.by_full_loci[,grep("P.22005",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$F.22004<-counts.by_full_loci$F.22004.batchCT
counts.by_full_loci$P.22004<-counts.by_full_loci$P.22004.batchCT
counts.by_full_loci$F.22003<-counts.by_full_loci$F.22003.batchCT
counts.by_full_loci$P.22003<-counts.by_full_loci$P.22003.batchCT
counts.by_full_loci$F.22002<-apply(counts.by_full_loci[,grep("F.22002",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$P.22002<-counts.by_full_loci$P.22002.batch2
counts.by_full_loci$F.220011<-apply(counts.by_full_loci[,grep("F.220011",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$P.220011<-counts.by_full_loci$P.220011.batch2

counts.by_full_loci$F.6024<-counts.by_full_loci$F.6024.batchCT
counts.by_full_loci$P.6024<-counts.by_full_loci$P.6024.batchCT
counts.by_full_loci$F.6244<-counts.by_full_loci$F.6244.batchCT
counts.by_full_loci$P.6244<-counts.by_full_loci$P.6244.batchCT
counts.by_full_loci$F.6909<-apply(counts.by_full_loci[,grep("F.6909",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$P.6909<-counts.by_full_loci$P.6909.batch2
counts.by_full_loci$F.6966<-apply(counts.by_full_loci[,grep("F.6966",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$P.6966<-counts.by_full_loci$P.6966.batch2
counts.by_full_loci$F.8236<-apply(counts.by_full_loci[,grep("F.8236",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$P.8236<-counts.by_full_loci$P.8236.batch2
counts.by_full_loci$F.9075<-counts.by_full_loci$F.9075.batch1
counts.by_full_loci$P.9075<-counts.by_full_loci$P.9075.batchCT
counts.by_full_loci$F.9537<-apply(counts.by_full_loci[,grep("F.9537",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$P.9537<-counts.by_full_loci$P.9537.batch2
counts.by_full_loci$F.9543<-counts.by_full_loci$F.9543.batchCT
counts.by_full_loci$P.9543<-counts.by_full_loci$P.9543.batchCT
counts.by_full_loci$F.9638<-counts.by_full_loci$F.9638.batchCT
counts.by_full_loci$P.9638<-counts.by_full_loci$P.9638.batchCT
counts.by_full_loci$F.9728<-counts.by_full_loci$F.9728.batchCT
counts.by_full_loci$P.9728<-counts.by_full_loci$P.9728.batchCT
counts.by_full_loci$F.9764<-apply(counts.by_full_loci[,grep("F.9764",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$P.9764<-apply(counts.by_full_loci[,grep("P.9764",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$F.9888<-counts.by_full_loci$F.9888.batchCT
counts.by_full_loci$P.9888<-counts.by_full_loci$P.9888.batchCT
counts.by_full_loci$F.9905<-counts.by_full_loci$F.9905.batchCT
counts.by_full_loci$P.9905<-counts.by_full_loci$P.9905.batchCT
counts.by_full_loci$F.9981<-apply(counts.by_full_loci[,grep("F.9981",names(counts.by_full_loci))],1,mean)
counts.by_full_loci$P.9981<-counts.by_full_loci$P.9981.batch2

counts.by_full_loci<-counts.by_full_loci[,grep("batch",names(counts.by_full_loci),invert = T)]
counts.by_full_loci<-counts.by_full_loci[,grep("rep",names(counts.by_full_loci),invert = T)]


counts.by_full_loci$max.flowers<- apply(counts.by_full_loci[,grep("F.",names(counts.by_full_loci))],1,max,na.rm=T)
counts.by_full_loci$max.pollen<- apply(counts.by_full_loci[,grep("P.",names(counts.by_full_loci))],1,max,na.rm=T)
counts.by_full_loci$max.rosette<- apply(counts.by_full_loci[,grep("R.",names(counts.by_full_loci))],1,max,na.rm=T)
counts.by_full_loci$max.seedl<- apply(counts.by_full_loci[,grep("S.",names(counts.by_full_loci))],1,max,na.rm=T)

counts.by_full_loci$mean.flowers<- apply(counts.by_full_loci[,grep("F.",names(counts.by_full_loci))],1,mean,na.rm=T)
counts.by_full_loci$mean.pollen<- apply(counts.by_full_loci[,grep("P.",names(counts.by_full_loci))],1,mean,na.rm=T)
counts.by_full_loci$mean.rosette<- apply(counts.by_full_loci[,grep("R.",names(counts.by_full_loci))],1,mean,na.rm=T)
counts.by_full_loci$mean.seedl<- apply(counts.by_full_loci[,grep("S.",names(counts.by_full_loci))],1,mean,na.rm=T)

counts.by_full_loci$Nacc_expressed.flowers<- apply(counts.by_full_loci[,grep("F.",names(counts.by_full_loci))],1,function(i) sum(i > 5,na.rm=T))
counts.by_full_loci$Nacc_expressed.pollen<- apply(counts.by_full_loci[,grep("P.",names(counts.by_full_loci))],1,function(i) sum(i > 5,na.rm=T))
counts.by_full_loci$Nacc_expressed.rosette<- apply(counts.by_full_loci[,grep("R.",names(counts.by_full_loci))],1,function(i) sum(i > 5,na.rm=T))
counts.by_full_loci$Nacc_expressed.seedl<- apply(counts.by_full_loci[,grep("S.",names(counts.by_full_loci))],1,function(i) sum(i > 5,na.rm=T))

counts.by_full_loci$freq<- apply(counts.by_full_loci[,grep("R.",names(counts.by_full_loci))],1,function(i) sum(!is.na(i)))

rownames(counts.by_full_loci)<-counts.by_full_loci$gene


#remove intermediate count tables for individual  accessions
rm(list = ls(pattern = "^a[0-9]+\\.counts\\.by_full_loci$"))


#########################################################

write.table(counts.by_full_loci,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/counts.full_loci.mapped_to_own_genomes.txt",quote = F, sep="\t",col.names = T,row.names = F)


########### counts on own genomes by exons (mRNAs)
#########################################################
#on own genomes  by exons (mRNAs)
i=1
acc=as.character(accessions_w_220011$V1[i])
a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/",acc,".counts.by_exons.bed",sep=""),header = T)
assign(paste("a",acc,".counts.by_exons",sep=""),a)
counts.by_exons<-a

for ( i in 2:27){
  acc=as.character(accessions_w_220011$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/",acc,".counts.by_exons.bed",sep=""),header = T)
  assign(paste("a",acc,".counts.by_exons",sep=""),a)
  counts.by_exons<-merge(counts.by_exons,a,by="gene",all = T)
}
counts.by_exons.unprocessed<-counts.by_exons

# average those samples that have replicates, remove batch names from sample names
counts.by_exons$F.10002<-apply(counts.by_exons[,grep("F.10002",names(counts.by_exons))],1,mean)
counts.by_exons$P.10002<-counts.by_exons$P.10002.batch2
counts.by_exons$F.10015<-apply(counts.by_exons[,grep("F.10015",names(counts.by_exons))],1,mean)
counts.by_exons$P.10015<-apply(counts.by_exons[,grep("P.10015",names(counts.by_exons))],1,mean)
counts.by_exons$S.10015<-apply(counts.by_exons[,grep("S.10015",names(counts.by_exons))],1,mean)
counts.by_exons$R.10015<-apply(counts.by_exons[,grep("R.10015",names(counts.by_exons))],1,mean)
counts.by_exons$F.10024<-apply(counts.by_exons[,grep("F.10024",names(counts.by_exons))],1,mean)
counts.by_exons$P.10024<-counts.by_exons$P.10024.batch2
counts.by_exons$F.1741<-counts.by_exons$F.1741.batchCT
counts.by_exons$P.1741<-counts.by_exons$P.1741.batchCT

counts.by_exons$F.22007<-apply(counts.by_exons[,grep("F.22007",names(counts.by_exons))],1,mean)
counts.by_exons$P.22007<-apply(counts.by_exons[,grep("P.22007",names(counts.by_exons))],1,mean)
counts.by_exons$F.22006<-counts.by_exons$F.22006.batchCT
counts.by_exons$P.22006<-counts.by_exons$P.22006.batchCT
counts.by_exons$F.22005<-apply(counts.by_exons[,grep("F.22005",names(counts.by_exons))],1,mean)
counts.by_exons$P.22005<-apply(counts.by_exons[,grep("P.22005",names(counts.by_exons))],1,mean)
counts.by_exons$F.22004<-counts.by_exons$F.22004.batchCT
counts.by_exons$P.22004<-counts.by_exons$P.22004.batchCT
counts.by_exons$F.22003<-counts.by_exons$F.22003.batchCT
counts.by_exons$P.22003<-counts.by_exons$P.22003.batchCT
counts.by_exons$F.22002<-apply(counts.by_exons[,grep("F.22002",names(counts.by_exons))],1,mean)
counts.by_exons$P.22002<-counts.by_exons$P.22002.batch2
counts.by_exons$F.220011<-apply(counts.by_exons[,grep("F.220011",names(counts.by_exons))],1,mean)
counts.by_exons$P.220011<-counts.by_exons$P.220011.batch2

counts.by_exons$F.6024<-counts.by_exons$F.6024.batchCT
counts.by_exons$P.6024<-counts.by_exons$P.6024.batchCT
counts.by_exons$F.6244<-counts.by_exons$F.6244.batchCT
counts.by_exons$P.6244<-counts.by_exons$P.6244.batchCT
counts.by_exons$F.6909<-apply(counts.by_exons[,grep("F.6909",names(counts.by_exons))],1,mean)
counts.by_exons$P.6909<-counts.by_exons$P.6909.batch2
counts.by_exons$F.6966<-apply(counts.by_exons[,grep("F.6966",names(counts.by_exons))],1,mean)
counts.by_exons$P.6966<-counts.by_exons$P.6966.batch2
counts.by_exons$F.8236<-apply(counts.by_exons[,grep("F.8236",names(counts.by_exons))],1,mean)
counts.by_exons$P.8236<-counts.by_exons$P.8236.batch2
counts.by_exons$F.9075<-counts.by_exons$F.9075.batch1
counts.by_exons$P.9075<-counts.by_exons$P.9075.batchCT
counts.by_exons$F.9537<-apply(counts.by_exons[,grep("F.9537",names(counts.by_exons))],1,mean)
counts.by_exons$P.9537<-counts.by_exons$P.9537.batch2
counts.by_exons$F.9543<-counts.by_exons$F.9543.batchCT
counts.by_exons$P.9543<-counts.by_exons$P.9543.batchCT
counts.by_exons$F.9638<-counts.by_exons$F.9638.batchCT
counts.by_exons$P.9638<-counts.by_exons$P.9638.batchCT
counts.by_exons$F.9728<-counts.by_exons$F.9728.batchCT
counts.by_exons$P.9728<-counts.by_exons$P.9728.batchCT
counts.by_exons$F.9764<-apply(counts.by_exons[,grep("F.9764",names(counts.by_exons))],1,mean)
counts.by_exons$P.9764<-apply(counts.by_exons[,grep("P.9764",names(counts.by_exons))],1,mean)
counts.by_exons$F.9888<-counts.by_exons$F.9888.batchCT
counts.by_exons$P.9888<-counts.by_exons$P.9888.batchCT
counts.by_exons$F.9905<-counts.by_exons$F.9905.batchCT
counts.by_exons$P.9905<-counts.by_exons$P.9905.batchCT
counts.by_exons$F.9981<-apply(counts.by_exons[,grep("F.9981",names(counts.by_exons))],1,mean)
counts.by_exons$P.9981<-counts.by_exons$P.9981.batch2

counts.by_exons<-counts.by_exons[,grep("batch",names(counts.by_exons),invert = T)]
counts.by_exons<-counts.by_exons[,grep("rep",names(counts.by_exons),invert = T)]


counts.by_exons$max.flowers<- apply(counts.by_exons[,grep("F.",names(counts.by_exons))],1,max,na.rm=T)
counts.by_exons$max.pollen<- apply(counts.by_exons[,grep("P.",names(counts.by_exons))],1,max,na.rm=T)
counts.by_exons$max.rosette<- apply(counts.by_exons[,grep("R.",names(counts.by_exons))],1,max,na.rm=T)
counts.by_exons$max.seedl<- apply(counts.by_exons[,grep("S.",names(counts.by_exons))],1,max,na.rm=T)

counts.by_exons$mean.flowers<- apply(counts.by_exons[,grep("F.",names(counts.by_exons))],1,mean,na.rm=T)
counts.by_exons$mean.pollen<- apply(counts.by_exons[,grep("P.",names(counts.by_exons))],1,mean,na.rm=T)
counts.by_exons$mean.rosette<- apply(counts.by_exons[,grep("R.",names(counts.by_exons))],1,mean,na.rm=T)
counts.by_exons$mean.seedl<- apply(counts.by_exons[,grep("S.",names(counts.by_exons))],1,mean,na.rm=T)

counts.by_exons$Nacc_expressed.flowers<- apply(counts.by_exons[,grep("F.",names(counts.by_exons))],1,function(i) sum(i > 5,na.rm=T))
counts.by_exons$Nacc_expressed.pollen<- apply(counts.by_exons[,grep("P.",names(counts.by_exons))],1,function(i) sum(i > 5,na.rm=T))
counts.by_exons$Nacc_expressed.rosette<- apply(counts.by_exons[,grep("R.",names(counts.by_exons))],1,function(i) sum(i > 5,na.rm=T))
counts.by_exons$Nacc_expressed.seedl<- apply(counts.by_exons[,grep("S.",names(counts.by_exons))],1,function(i) sum(i > 5,na.rm=T))

counts.by_exons$freq<- apply(counts.by_exons[,grep("R.",names(counts.by_exons))],1,function(i) sum(!is.na(i)))

rownames(counts.by_exons)<-counts.by_exons$gene


#remove intermediate count tables for individual  accessions
rm(list = ls(pattern = "^a[0-9]+\\.counts\\.by_exons$"))


#########################################################

write.table(counts.by_exons,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/counts.mRNAs.mapped_to_own_genomes.txt",quote = F, sep="\t",col.names = T,row.names = F)

#TPMs on TAIR10 by exons
#on TAIR10 by exons (mRNAs )
##############################################

TPMs.Araport_exons<-read.delim("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/onTAIR10/allsamples.TPMs.by_exons.on_TAIR10.bed",header = T)

TPMs.Araport_exons.unprocessed<-TPMs.Araport_exons

TPMs.Araport_exons$F.10002<-apply(TPMs.Araport_exons[,grep("F.10002",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$P.10002<-TPMs.Araport_exons$P.10002.batch2
TPMs.Araport_exons$F.10015<-apply(TPMs.Araport_exons[,grep("F.10015",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$P.10015<-apply(TPMs.Araport_exons[,grep("P.10015",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$S.10015<-apply(TPMs.Araport_exons[,grep("S.10015",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$R.10015<-apply(TPMs.Araport_exons[,grep("R.10015",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$F.10024<-apply(TPMs.Araport_exons[,grep("F.10024",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$P.10024<-TPMs.Araport_exons$P.10024.batch2
TPMs.Araport_exons$F.1741<-TPMs.Araport_exons$F.1741.batchCT
TPMs.Araport_exons$P.1741<-TPMs.Araport_exons$P.1741.batchCT

TPMs.Araport_exons$F.22007<-apply(TPMs.Araport_exons[,grep("F.22007",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$P.22007<-apply(TPMs.Araport_exons[,grep("P.22007",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$F.22006<-TPMs.Araport_exons$F.22006.batchCT
TPMs.Araport_exons$P.22006<-TPMs.Araport_exons$P.22006.batchCT
TPMs.Araport_exons$F.22005<-apply(TPMs.Araport_exons[,grep("F.22005",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$P.22005<-apply(TPMs.Araport_exons[,grep("P.22005",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$F.22004<-TPMs.Araport_exons$F.22004.batchCT
TPMs.Araport_exons$P.22004<-TPMs.Araport_exons$P.22004.batchCT
TPMs.Araport_exons$F.22003<-TPMs.Araport_exons$F.22003.batchCT
TPMs.Araport_exons$P.22003<-TPMs.Araport_exons$P.22003.batchCT
TPMs.Araport_exons$F.22002<-apply(TPMs.Araport_exons[,grep("F.22002",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$P.22002<-TPMs.Araport_exons$P.22002.batch2
TPMs.Araport_exons$F.220011<-apply(TPMs.Araport_exons[,grep("F.220011",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$P.220011<-TPMs.Araport_exons$P.220011.batch2

TPMs.Araport_exons$F.6024<-TPMs.Araport_exons$F.6024.batchCT
TPMs.Araport_exons$P.6024<-TPMs.Araport_exons$P.6024.batchCT
TPMs.Araport_exons$F.6244<-TPMs.Araport_exons$F.6244.batchCT
TPMs.Araport_exons$P.6244<-TPMs.Araport_exons$P.6244.batchCT
TPMs.Araport_exons$F.6909<-apply(TPMs.Araport_exons[,grep("F.6909",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$P.6909<-TPMs.Araport_exons$P.6909.batch2
TPMs.Araport_exons$F.6966<-apply(TPMs.Araport_exons[,grep("F.6966",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$P.6966<-TPMs.Araport_exons$P.6966.batch2
TPMs.Araport_exons$F.8236<-apply(TPMs.Araport_exons[,grep("F.8236",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$P.8236<-TPMs.Araport_exons$P.8236.batch2
TPMs.Araport_exons$F.9075<-TPMs.Araport_exons$F.9075.batch1
TPMs.Araport_exons$P.9075<-TPMs.Araport_exons$P.9075.batchCT
TPMs.Araport_exons$F.9537<-apply(TPMs.Araport_exons[,grep("F.9537",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$P.9537<-TPMs.Araport_exons$P.9537.batch2
TPMs.Araport_exons$F.9543<-TPMs.Araport_exons$F.9543.batchCT
TPMs.Araport_exons$P.9543<-TPMs.Araport_exons$P.9543.batchCT
TPMs.Araport_exons$F.9638<-TPMs.Araport_exons$F.9638.batchCT
TPMs.Araport_exons$P.9638<-TPMs.Araport_exons$P.9638.batchCT
TPMs.Araport_exons$F.9728<-TPMs.Araport_exons$F.9728.batchCT
TPMs.Araport_exons$P.9728<-TPMs.Araport_exons$P.9728.batchCT
TPMs.Araport_exons$F.9764<-apply(TPMs.Araport_exons[,grep("F.9764",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$P.9764<-apply(TPMs.Araport_exons[,grep("P.9764",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$F.9888<-TPMs.Araport_exons$F.9888.batchCT
TPMs.Araport_exons$P.9888<-TPMs.Araport_exons$P.9888.batchCT
TPMs.Araport_exons$F.9905<-TPMs.Araport_exons$F.9905.batchCT
TPMs.Araport_exons$P.9905<-TPMs.Araport_exons$P.9905.batchCT
TPMs.Araport_exons$F.9981<-apply(TPMs.Araport_exons[,grep("F.9981",names(TPMs.Araport_exons))],1,mean)
TPMs.Araport_exons$P.9981<-TPMs.Araport_exons$P.9981.batch2

TPMs.Araport_exons<-TPMs.Araport_exons[,grep("batch",names(TPMs.Araport_exons),invert = T)]
TPMs.Araport_exons<-TPMs.Araport_exons[,grep("rep",names(TPMs.Araport_exons),invert = T)]


TPMs.Araport_exons$max.flowers<- apply(TPMs.Araport_exons[,grep("F.",names(TPMs.Araport_exons))],1,max,na.rm=T)
TPMs.Araport_exons$max.pollen<- apply(TPMs.Araport_exons[,grep("P.",names(TPMs.Araport_exons))],1,max,na.rm=T)
TPMs.Araport_exons$max.rosette<- apply(TPMs.Araport_exons[,grep("R.",names(TPMs.Araport_exons))],1,max,na.rm=T)
TPMs.Araport_exons$max.seedl<- apply(TPMs.Araport_exons[,grep("S.",names(TPMs.Araport_exons))],1,max,na.rm=T)

TPMs.Araport_exons$mean.flowers<- apply(TPMs.Araport_exons[,grep("F.",names(TPMs.Araport_exons))],1,mean,na.rm=T)
TPMs.Araport_exons$mean.pollen<- apply(TPMs.Araport_exons[,grep("P.",names(TPMs.Araport_exons))],1,mean,na.rm=T)
TPMs.Araport_exons$mean.rosette<- apply(TPMs.Araport_exons[,grep("R.",names(TPMs.Araport_exons))],1,mean,na.rm=T)
TPMs.Araport_exons$mean.seedl<- apply(TPMs.Araport_exons[,grep("S.",names(TPMs.Araport_exons))],1,mean,na.rm=T)

TPMs.Araport_exons$Nacc_expressed.flowers<- apply(TPMs.Araport_exons[,grep("F.",names(TPMs.Araport_exons))],1,function(i) sum(i > 0.5,na.rm=T))
TPMs.Araport_exons$Nacc_expressed.pollen<- apply(TPMs.Araport_exons[,grep("P.",names(TPMs.Araport_exons))],1,function(i) sum(i > 0.5,na.rm=T))
TPMs.Araport_exons$Nacc_expressed.rosette<- apply(TPMs.Araport_exons[,grep("R.",names(TPMs.Araport_exons))],1,function(i) sum(i > 0.5,na.rm=T))
TPMs.Araport_exons$Nacc_expressed.seedl<- apply(TPMs.Araport_exons[,grep("S.",names(TPMs.Araport_exons))],1,function(i) sum(i > 0.5,na.rm=T))

rownames(TPMs.Araport_exons)<-TPMs.Araport_exons$gene

############################################################

#counts on TAIR10 by exons
#on TAIR10 by exons (mRNAs )
##############################################
counts.Araport_exons<-read.delim("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/onTAIR10/allsamples.counts.by_exons.on_TAIR10.bed",header = T)
counts.Araport_exons.unprocessed<-counts.Araport_exons

counts.Araport_exons$F.10002<-apply(counts.Araport_exons[,grep("F.10002",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$P.10002<-counts.Araport_exons$P.10002.batch2
counts.Araport_exons$F.10015<-apply(counts.Araport_exons[,grep("F.10015",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$P.10015<-apply(counts.Araport_exons[,grep("P.10015",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$S.10015<-apply(counts.Araport_exons[,grep("S.10015",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$R.10015<-apply(counts.Araport_exons[,grep("R.10015",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$F.10024<-apply(counts.Araport_exons[,grep("F.10024",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$P.10024<-counts.Araport_exons$P.10024.batch2
counts.Araport_exons$F.1741<-counts.Araport_exons$F.1741.batchCT
counts.Araport_exons$P.1741<-counts.Araport_exons$P.1741.batchCT

counts.Araport_exons$F.22007<-apply(counts.Araport_exons[,grep("F.22007",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$P.22007<-apply(counts.Araport_exons[,grep("P.22007",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$F.22006<-counts.Araport_exons$F.22006.batchCT
counts.Araport_exons$P.22006<-counts.Araport_exons$P.22006.batchCT
counts.Araport_exons$F.22005<-apply(counts.Araport_exons[,grep("F.22005",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$P.22005<-apply(counts.Araport_exons[,grep("P.22005",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$F.22004<-counts.Araport_exons$F.22004.batchCT
counts.Araport_exons$P.22004<-counts.Araport_exons$P.22004.batchCT
counts.Araport_exons$F.22003<-counts.Araport_exons$F.22003.batchCT
counts.Araport_exons$P.22003<-counts.Araport_exons$P.22003.batchCT
counts.Araport_exons$F.22002<-apply(counts.Araport_exons[,grep("F.22002",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$P.22002<-counts.Araport_exons$P.22002.batch2
counts.Araport_exons$F.220011<-apply(counts.Araport_exons[,grep("F.220011",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$P.220011<-counts.Araport_exons$P.220011.batch2

counts.Araport_exons$F.6024<-counts.Araport_exons$F.6024.batchCT
counts.Araport_exons$P.6024<-counts.Araport_exons$P.6024.batchCT
counts.Araport_exons$F.6244<-counts.Araport_exons$F.6244.batchCT
counts.Araport_exons$P.6244<-counts.Araport_exons$P.6244.batchCT
counts.Araport_exons$F.6909<-apply(counts.Araport_exons[,grep("F.6909",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$P.6909<-counts.Araport_exons$P.6909.batch2
counts.Araport_exons$F.6966<-apply(counts.Araport_exons[,grep("F.6966",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$P.6966<-counts.Araport_exons$P.6966.batch2
counts.Araport_exons$F.8236<-apply(counts.Araport_exons[,grep("F.8236",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$P.8236<-counts.Araport_exons$P.8236.batch2
counts.Araport_exons$F.9075<-counts.Araport_exons$F.9075.batch1
counts.Araport_exons$P.9075<-counts.Araport_exons$P.9075.batchCT
counts.Araport_exons$F.9537<-apply(counts.Araport_exons[,grep("F.9537",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$P.9537<-counts.Araport_exons$P.9537.batch2
counts.Araport_exons$F.9543<-counts.Araport_exons$F.9543.batchCT
counts.Araport_exons$P.9543<-counts.Araport_exons$P.9543.batchCT
counts.Araport_exons$F.9638<-counts.Araport_exons$F.9638.batchCT
counts.Araport_exons$P.9638<-counts.Araport_exons$P.9638.batchCT
counts.Araport_exons$F.9728<-counts.Araport_exons$F.9728.batchCT
counts.Araport_exons$P.9728<-counts.Araport_exons$P.9728.batchCT
counts.Araport_exons$F.9764<-apply(counts.Araport_exons[,grep("F.9764",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$P.9764<-apply(counts.Araport_exons[,grep("P.9764",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$F.9888<-counts.Araport_exons$F.9888.batchCT
counts.Araport_exons$P.9888<-counts.Araport_exons$P.9888.batchCT
counts.Araport_exons$F.9905<-counts.Araport_exons$F.9905.batchCT
counts.Araport_exons$P.9905<-counts.Araport_exons$P.9905.batchCT
counts.Araport_exons$F.9981<-apply(counts.Araport_exons[,grep("F.9981",names(counts.Araport_exons))],1,mean)
counts.Araport_exons$P.9981<-counts.Araport_exons$P.9981.batch2

counts.Araport_exons<-counts.Araport_exons[,grep("batch",names(counts.Araport_exons),invert = T)]
counts.Araport_exons<-counts.Araport_exons[,grep("rep",names(counts.Araport_exons),invert = T)]


counts.Araport_exons$max.flowers<- apply(counts.Araport_exons[,grep("F.",names(counts.Araport_exons))],1,max,na.rm=T)
counts.Araport_exons$max.pollen<- apply(counts.Araport_exons[,grep("P.",names(counts.Araport_exons))],1,max,na.rm=T)
counts.Araport_exons$max.rosette<- apply(counts.Araport_exons[,grep("R.",names(counts.Araport_exons))],1,max,na.rm=T)
counts.Araport_exons$max.seedl<- apply(counts.Araport_exons[,grep("S.",names(counts.Araport_exons))],1,max,na.rm=T)

counts.Araport_exons$mean.flowers<- apply(counts.Araport_exons[,grep("F.",names(counts.Araport_exons))],1,mean,na.rm=T)
counts.Araport_exons$mean.pollen<- apply(counts.Araport_exons[,grep("P.",names(counts.Araport_exons))],1,mean,na.rm=T)
counts.Araport_exons$mean.rosette<- apply(counts.Araport_exons[,grep("R.",names(counts.Araport_exons))],1,mean,na.rm=T)
counts.Araport_exons$mean.seedl<- apply(counts.Araport_exons[,grep("S.",names(counts.Araport_exons))],1,mean,na.rm=T)

counts.Araport_exons$Nacc_expressed.flowers<- apply(counts.Araport_exons[,grep("F.",names(counts.Araport_exons))],1,function(i) sum(i > 5,na.rm=T))
counts.Araport_exons$Nacc_expressed.pollen<- apply(counts.Araport_exons[,grep("P.",names(counts.Araport_exons))],1,function(i) sum(i > 5,na.rm=T))
counts.Araport_exons$Nacc_expressed.rosette<- apply(counts.Araport_exons[,grep("R.",names(counts.Araport_exons))],1,function(i) sum(i > 5,na.rm=T))
counts.Araport_exons$Nacc_expressed.seedl<- apply(counts.Araport_exons[,grep("S.",names(counts.Araport_exons))],1,function(i) sum(i > 5,na.rm=T))

rownames(counts.Araport_exons)<-counts.Araport_exons$gene

############################################################

########### TPMs on TAIR10 by full loci 
#########################################################
#on TAIR10  by full loci 
TPMs.Araport_locus<-read.delim("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/onTAIR10/allsamples.TPMs.by_full_loci.on_TAIR10.bed",header = T)
TPMs.Araport_locus.unprocessed<-TPMs.Araport_locus

TPMs.Araport_locus$F.10002<-apply(TPMs.Araport_locus[,grep("F.10002",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$P.10002<-TPMs.Araport_locus$P.10002.batch2
TPMs.Araport_locus$F.10015<-apply(TPMs.Araport_locus[,grep("F.10015",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$P.10015<-apply(TPMs.Araport_locus[,grep("P.10015",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$S.10015<-apply(TPMs.Araport_locus[,grep("S.10015",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$R.10015<-apply(TPMs.Araport_locus[,grep("R.10015",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$F.10024<-apply(TPMs.Araport_locus[,grep("F.10024",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$P.10024<-TPMs.Araport_locus$P.10024.batch2
TPMs.Araport_locus$F.1741<-TPMs.Araport_locus$F.1741.batchCT
TPMs.Araport_locus$P.1741<-TPMs.Araport_locus$P.1741.batchCT

TPMs.Araport_locus$F.22007<-apply(TPMs.Araport_locus[,grep("F.22007",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$P.22007<-apply(TPMs.Araport_locus[,grep("P.22007",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$F.22006<-TPMs.Araport_locus$F.22006.batchCT
TPMs.Araport_locus$P.22006<-TPMs.Araport_locus$P.22006.batchCT
TPMs.Araport_locus$F.22005<-apply(TPMs.Araport_locus[,grep("F.22005",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$P.22005<-apply(TPMs.Araport_locus[,grep("P.22005",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$F.22004<-TPMs.Araport_locus$F.22004.batchCT
TPMs.Araport_locus$P.22004<-TPMs.Araport_locus$P.22004.batchCT
TPMs.Araport_locus$F.22003<-TPMs.Araport_locus$F.22003.batchCT
TPMs.Araport_locus$P.22003<-TPMs.Araport_locus$P.22003.batchCT
TPMs.Araport_locus$F.22002<-apply(TPMs.Araport_locus[,grep("F.22002",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$P.22002<-TPMs.Araport_locus$P.22002.batch2
TPMs.Araport_locus$F.220011<-apply(TPMs.Araport_locus[,grep("F.220011",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$P.220011<-TPMs.Araport_locus$P.220011.batch2

TPMs.Araport_locus$F.6024<-TPMs.Araport_locus$F.6024.batchCT
TPMs.Araport_locus$P.6024<-TPMs.Araport_locus$P.6024.batchCT
TPMs.Araport_locus$F.6244<-TPMs.Araport_locus$F.6244.batchCT
TPMs.Araport_locus$P.6244<-TPMs.Araport_locus$P.6244.batchCT
TPMs.Araport_locus$F.6909<-apply(TPMs.Araport_locus[,grep("F.6909",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$P.6909<-TPMs.Araport_locus$P.6909.batch2
TPMs.Araport_locus$F.6966<-apply(TPMs.Araport_locus[,grep("F.6966",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$P.6966<-TPMs.Araport_locus$P.6966.batch2
TPMs.Araport_locus$F.8236<-apply(TPMs.Araport_locus[,grep("F.8236",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$P.8236<-TPMs.Araport_locus$P.8236.batch2
TPMs.Araport_locus$F.9075<-TPMs.Araport_locus$F.9075.batch1
TPMs.Araport_locus$P.9075<-TPMs.Araport_locus$P.9075.batchCT
TPMs.Araport_locus$F.9537<-apply(TPMs.Araport_locus[,grep("F.9537",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$P.9537<-TPMs.Araport_locus$P.9537.batch2
TPMs.Araport_locus$F.9543<-TPMs.Araport_locus$F.9543.batchCT
TPMs.Araport_locus$P.9543<-TPMs.Araport_locus$P.9543.batchCT
TPMs.Araport_locus$F.9638<-TPMs.Araport_locus$F.9638.batchCT
TPMs.Araport_locus$P.9638<-TPMs.Araport_locus$P.9638.batchCT
TPMs.Araport_locus$F.9728<-TPMs.Araport_locus$F.9728.batchCT
TPMs.Araport_locus$P.9728<-TPMs.Araport_locus$P.9728.batchCT
TPMs.Araport_locus$F.9764<-apply(TPMs.Araport_locus[,grep("F.9764",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$P.9764<-apply(TPMs.Araport_locus[,grep("P.9764",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$F.9888<-TPMs.Araport_locus$F.9888.batchCT
TPMs.Araport_locus$P.9888<-TPMs.Araport_locus$P.9888.batchCT
TPMs.Araport_locus$F.9905<-TPMs.Araport_locus$F.9905.batchCT
TPMs.Araport_locus$P.9905<-TPMs.Araport_locus$P.9905.batchCT
TPMs.Araport_locus$F.9981<-apply(TPMs.Araport_locus[,grep("F.9981",names(TPMs.Araport_locus))],1,mean)
TPMs.Araport_locus$P.9981<-TPMs.Araport_locus$P.9981.batch2

TPMs.Araport_locus<-TPMs.Araport_locus[,grep("batch",names(TPMs.Araport_locus),invert = T)]
TPMs.Araport_locus<-TPMs.Araport_locus[,grep("rep",names(TPMs.Araport_locus),invert = T)]


TPMs.Araport_locus$max.flowers<- apply(TPMs.Araport_locus[,grep("F.",names(TPMs.Araport_locus))],1,max,na.rm=T)
TPMs.Araport_locus$max.pollen<- apply(TPMs.Araport_locus[,grep("P.",names(TPMs.Araport_locus))],1,max,na.rm=T)
TPMs.Araport_locus$max.rosette<- apply(TPMs.Araport_locus[,grep("R.",names(TPMs.Araport_locus))],1,max,na.rm=T)
TPMs.Araport_locus$max.seedl<- apply(TPMs.Araport_locus[,grep("S.",names(TPMs.Araport_locus))],1,max,na.rm=T)

TPMs.Araport_locus$mean.flowers<- apply(TPMs.Araport_locus[,grep("F.",names(TPMs.Araport_locus))],1,mean,na.rm=T)
TPMs.Araport_locus$mean.pollen<- apply(TPMs.Araport_locus[,grep("P.",names(TPMs.Araport_locus))],1,mean,na.rm=T)
TPMs.Araport_locus$mean.rosette<- apply(TPMs.Araport_locus[,grep("R.",names(TPMs.Araport_locus))],1,mean,na.rm=T)
TPMs.Araport_locus$mean.seedl<- apply(TPMs.Araport_locus[,grep("S.",names(TPMs.Araport_locus))],1,mean,na.rm=T)

TPMs.Araport_locus$Nacc_expressed.flowers<- apply(TPMs.Araport_locus[,grep("F.",names(TPMs.Araport_locus))],1,function(i) sum(i > 0.25,na.rm=T))
TPMs.Araport_locus$Nacc_expressed.pollen<- apply(TPMs.Araport_locus[,grep("P.",names(TPMs.Araport_locus))],1,function(i) sum(i > 0.25,na.rm=T))
TPMs.Araport_locus$Nacc_expressed.rosette<- apply(TPMs.Araport_locus[,grep("R.",names(TPMs.Araport_locus))],1,function(i) sum(i > 0.25,na.rm=T))
TPMs.Araport_locus$Nacc_expressed.seedl<- apply(TPMs.Araport_locus[,grep("S.",names(TPMs.Araport_locus))],1,function(i) sum(i > 0.25,na.rm=T))

rownames(TPMs.Araport_locus)<-TPMs.Araport_locus$gene

#########################################################

########### counts on TAIR10 by full loci 
#########################################################
#on TAIR10  by full loci 
counts.Araport_locus<-read.delim("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/onTAIR10/allsamples.counts.by_full_loci.on_TAIR10.bed",header = T)
counts.Araport_locus.unprocessed<-counts.Araport_locus

counts.Araport_locus$F.10002<-apply(counts.Araport_locus[,grep("F.10002",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$P.10002<-counts.Araport_locus$P.10002.batch2
counts.Araport_locus$F.10015<-apply(counts.Araport_locus[,grep("F.10015",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$P.10015<-apply(counts.Araport_locus[,grep("P.10015",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$S.10015<-apply(counts.Araport_locus[,grep("S.10015",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$R.10015<-apply(counts.Araport_locus[,grep("R.10015",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$F.10024<-apply(counts.Araport_locus[,grep("F.10024",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$P.10024<-counts.Araport_locus$P.10024.batch2
counts.Araport_locus$F.1741<-counts.Araport_locus$F.1741.batchCT
counts.Araport_locus$P.1741<-counts.Araport_locus$P.1741.batchCT

counts.Araport_locus$F.22007<-apply(counts.Araport_locus[,grep("F.22007",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$P.22007<-apply(counts.Araport_locus[,grep("P.22007",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$F.22006<-counts.Araport_locus$F.22006.batchCT
counts.Araport_locus$P.22006<-counts.Araport_locus$P.22006.batchCT
counts.Araport_locus$F.22005<-apply(counts.Araport_locus[,grep("F.22005",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$P.22005<-apply(counts.Araport_locus[,grep("P.22005",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$F.22004<-counts.Araport_locus$F.22004.batchCT
counts.Araport_locus$P.22004<-counts.Araport_locus$P.22004.batchCT
counts.Araport_locus$F.22003<-counts.Araport_locus$F.22003.batchCT
counts.Araport_locus$P.22003<-counts.Araport_locus$P.22003.batchCT
counts.Araport_locus$F.22002<-apply(counts.Araport_locus[,grep("F.22002",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$P.22002<-counts.Araport_locus$P.22002.batch2
counts.Araport_locus$F.220011<-apply(counts.Araport_locus[,grep("F.220011",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$P.220011<-counts.Araport_locus$P.220011.batch2

counts.Araport_locus$F.6024<-counts.Araport_locus$F.6024.batchCT
counts.Araport_locus$P.6024<-counts.Araport_locus$P.6024.batchCT
counts.Araport_locus$F.6244<-counts.Araport_locus$F.6244.batchCT
counts.Araport_locus$P.6244<-counts.Araport_locus$P.6244.batchCT
counts.Araport_locus$F.6909<-apply(counts.Araport_locus[,grep("F.6909",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$P.6909<-counts.Araport_locus$P.6909.batch2
counts.Araport_locus$F.6966<-apply(counts.Araport_locus[,grep("F.6966",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$P.6966<-counts.Araport_locus$P.6966.batch2
counts.Araport_locus$F.8236<-apply(counts.Araport_locus[,grep("F.8236",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$P.8236<-counts.Araport_locus$P.8236.batch2
counts.Araport_locus$F.9075<-counts.Araport_locus$F.9075.batch1
counts.Araport_locus$P.9075<-counts.Araport_locus$P.9075.batchCT
counts.Araport_locus$F.9537<-apply(counts.Araport_locus[,grep("F.9537",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$P.9537<-counts.Araport_locus$P.9537.batch2
counts.Araport_locus$F.9543<-counts.Araport_locus$F.9543.batchCT
counts.Araport_locus$P.9543<-counts.Araport_locus$P.9543.batchCT
counts.Araport_locus$F.9638<-counts.Araport_locus$F.9638.batchCT
counts.Araport_locus$P.9638<-counts.Araport_locus$P.9638.batchCT
counts.Araport_locus$F.9728<-counts.Araport_locus$F.9728.batchCT
counts.Araport_locus$P.9728<-counts.Araport_locus$P.9728.batchCT
counts.Araport_locus$F.9764<-apply(counts.Araport_locus[,grep("F.9764",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$P.9764<-apply(counts.Araport_locus[,grep("P.9764",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$F.9888<-counts.Araport_locus$F.9888.batchCT
counts.Araport_locus$P.9888<-counts.Araport_locus$P.9888.batchCT
counts.Araport_locus$F.9905<-counts.Araport_locus$F.9905.batchCT
counts.Araport_locus$P.9905<-counts.Araport_locus$P.9905.batchCT
counts.Araport_locus$F.9981<-apply(counts.Araport_locus[,grep("F.9981",names(counts.Araport_locus))],1,mean)
counts.Araport_locus$P.9981<-counts.Araport_locus$P.9981.batch2

counts.Araport_locus<-counts.Araport_locus[,grep("batch",names(counts.Araport_locus),invert = T)]
counts.Araport_locus<-counts.Araport_locus[,grep("rep",names(counts.Araport_locus),invert = T)]


counts.Araport_locus$max.flowers<- apply(counts.Araport_locus[,grep("F.",names(counts.Araport_locus))],1,max,na.rm=T)
counts.Araport_locus$max.pollen<- apply(counts.Araport_locus[,grep("P.",names(counts.Araport_locus))],1,max,na.rm=T)
counts.Araport_locus$max.rosette<- apply(counts.Araport_locus[,grep("R.",names(counts.Araport_locus))],1,max,na.rm=T)
counts.Araport_locus$max.seedl<- apply(counts.Araport_locus[,grep("S.",names(counts.Araport_locus))],1,max,na.rm=T)

counts.Araport_locus$mean.flowers<- apply(counts.Araport_locus[,grep("F.",names(counts.Araport_locus))],1,mean,na.rm=T)
counts.Araport_locus$mean.pollen<- apply(counts.Araport_locus[,grep("P.",names(counts.Araport_locus))],1,mean,na.rm=T)
counts.Araport_locus$mean.rosette<- apply(counts.Araport_locus[,grep("R.",names(counts.Araport_locus))],1,mean,na.rm=T)
counts.Araport_locus$mean.seedl<- apply(counts.Araport_locus[,grep("S.",names(counts.Araport_locus))],1,mean,na.rm=T)

counts.Araport_locus$Nacc_expressed.flowers<- apply(counts.Araport_locus[,grep("F.",names(counts.Araport_locus))],1,function(i) sum(i > 5,na.rm=T))
counts.Araport_locus$Nacc_expressed.pollen<- apply(counts.Araport_locus[,grep("P.",names(counts.Araport_locus))],1,function(i) sum(i > 5,na.rm=T))
counts.Araport_locus$Nacc_expressed.rosette<- apply(counts.Araport_locus[,grep("R.",names(counts.Araport_locus))],1,function(i) sum(i > 5,na.rm=T))
counts.Araport_locus$Nacc_expressed.seedl<- apply(counts.Araport_locus[,grep("S.",names(counts.Araport_locus))],1,function(i) sum(i > 5,na.rm=T))

rownames(counts.Araport_locus)<-counts.Araport_locus$gene

#########################################################


write.table(TPMs.Araport_locus,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/TPMs.by_full_loci.mapped_to_TAIR10.txt",quote = F, sep="\t",col.names = T,row.names = F)
write.table(TPMs.Araport_exons,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/TPMs.by_exons.mapped_to_TAIR10.txt",quote = F, sep="\t",col.names = T,row.names = F)





# genes and SVs 

########### TPMs on own genomes by full loci - genes and SVs
#########################################################
#on own genomes  by full loci 
i=1
acc=as.character(accessions_w_220011$V1[i])
a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/",acc,".TPMs.by_full_loci_and_SVs.bed",sep=""),header = T)
assign(paste("a",acc,".TPMs.by_full_loci_and_SVs",sep=""),a)
TPMs.by_full_loci_and_SVs<-a

for ( i in 2:27){
  acc=as.character(accessions_w_220011$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/",acc,".TPMs.by_full_loci_and_SVs.bed",sep=""),header = T)
  assign(paste("a",acc,".TPMs.by_full_loci_and_SVs",sep=""),a)
  TPMs.by_full_loci_and_SVs<-merge(TPMs.by_full_loci_and_SVs,a,by="gene",all = T)
}

TPMs.by_full_loci_and_SVs.unprocessed<-TPMs.by_full_loci_and_SVs
# average those samples that have replicates, remove batch names from sample names
TPMs.by_full_loci_and_SVs$F.10002<-apply(TPMs.by_full_loci_and_SVs[,grep("F.10002",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$P.10002<-TPMs.by_full_loci_and_SVs$P.10002.batch2
TPMs.by_full_loci_and_SVs$F.10015<-apply(TPMs.by_full_loci_and_SVs[,grep("F.10015",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$P.10015<-apply(TPMs.by_full_loci_and_SVs[,grep("P.10015",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$S.10015<-apply(TPMs.by_full_loci_and_SVs[,grep("S.10015",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$R.10015<-apply(TPMs.by_full_loci_and_SVs[,grep("R.10015",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$F.10024<-apply(TPMs.by_full_loci_and_SVs[,grep("F.10024",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$P.10024<-TPMs.by_full_loci_and_SVs$P.10024.batch2
TPMs.by_full_loci_and_SVs$F.1741<-TPMs.by_full_loci_and_SVs$F.1741.batchCT
TPMs.by_full_loci_and_SVs$P.1741<-TPMs.by_full_loci_and_SVs$P.1741.batchCT

TPMs.by_full_loci_and_SVs$F.22007<-apply(TPMs.by_full_loci_and_SVs[,grep("F.22007",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$P.22007<-apply(TPMs.by_full_loci_and_SVs[,grep("P.22007",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$F.22006<-TPMs.by_full_loci_and_SVs$F.22006.batchCT
TPMs.by_full_loci_and_SVs$P.22006<-TPMs.by_full_loci_and_SVs$P.22006.batchCT
TPMs.by_full_loci_and_SVs$F.22005<-apply(TPMs.by_full_loci_and_SVs[,grep("F.22005",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$P.22005<-apply(TPMs.by_full_loci_and_SVs[,grep("P.22005",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$F.22004<-TPMs.by_full_loci_and_SVs$F.22004.batchCT
TPMs.by_full_loci_and_SVs$P.22004<-TPMs.by_full_loci_and_SVs$P.22004.batchCT
TPMs.by_full_loci_and_SVs$F.22003<-TPMs.by_full_loci_and_SVs$F.22003.batchCT
TPMs.by_full_loci_and_SVs$P.22003<-TPMs.by_full_loci_and_SVs$P.22003.batchCT
TPMs.by_full_loci_and_SVs$F.22002<-apply(TPMs.by_full_loci_and_SVs[,grep("F.22002",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$P.22002<-TPMs.by_full_loci_and_SVs$P.22002.batch2
TPMs.by_full_loci_and_SVs$F.220011<-apply(TPMs.by_full_loci_and_SVs[,grep("F.220011",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$P.220011<-TPMs.by_full_loci_and_SVs$P.220011.batch2

TPMs.by_full_loci_and_SVs$F.6024<-TPMs.by_full_loci_and_SVs$F.6024.batchCT
TPMs.by_full_loci_and_SVs$P.6024<-TPMs.by_full_loci_and_SVs$P.6024.batchCT
TPMs.by_full_loci_and_SVs$F.6244<-TPMs.by_full_loci_and_SVs$F.6244.batchCT
TPMs.by_full_loci_and_SVs$P.6244<-TPMs.by_full_loci_and_SVs$P.6244.batchCT
TPMs.by_full_loci_and_SVs$F.6909<-apply(TPMs.by_full_loci_and_SVs[,grep("F.6909",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$P.6909<-TPMs.by_full_loci_and_SVs$P.6909.batch2
TPMs.by_full_loci_and_SVs$F.6966<-apply(TPMs.by_full_loci_and_SVs[,grep("F.6966",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$P.6966<-TPMs.by_full_loci_and_SVs$P.6966.batch2
TPMs.by_full_loci_and_SVs$F.8236<-apply(TPMs.by_full_loci_and_SVs[,grep("F.8236",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$P.8236<-TPMs.by_full_loci_and_SVs$P.8236.batch2
TPMs.by_full_loci_and_SVs$F.9075<-TPMs.by_full_loci_and_SVs$F.9075.batch1
TPMs.by_full_loci_and_SVs$P.9075<-TPMs.by_full_loci_and_SVs$P.9075.batchCT
TPMs.by_full_loci_and_SVs$F.9537<-apply(TPMs.by_full_loci_and_SVs[,grep("F.9537",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$P.9537<-TPMs.by_full_loci_and_SVs$P.9537.batch2
TPMs.by_full_loci_and_SVs$F.9543<-TPMs.by_full_loci_and_SVs$F.9543.batchCT
TPMs.by_full_loci_and_SVs$P.9543<-TPMs.by_full_loci_and_SVs$P.9543.batchCT
TPMs.by_full_loci_and_SVs$F.9638<-TPMs.by_full_loci_and_SVs$F.9638.batchCT
TPMs.by_full_loci_and_SVs$P.9638<-TPMs.by_full_loci_and_SVs$P.9638.batchCT
TPMs.by_full_loci_and_SVs$F.9728<-TPMs.by_full_loci_and_SVs$F.9728.batchCT
TPMs.by_full_loci_and_SVs$P.9728<-TPMs.by_full_loci_and_SVs$P.9728.batchCT
TPMs.by_full_loci_and_SVs$F.9764<-apply(TPMs.by_full_loci_and_SVs[,grep("F.9764",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$P.9764<-apply(TPMs.by_full_loci_and_SVs[,grep("P.9764",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$F.9888<-TPMs.by_full_loci_and_SVs$F.9888.batchCT
TPMs.by_full_loci_and_SVs$P.9888<-TPMs.by_full_loci_and_SVs$P.9888.batchCT
TPMs.by_full_loci_and_SVs$F.9905<-TPMs.by_full_loci_and_SVs$F.9905.batchCT
TPMs.by_full_loci_and_SVs$P.9905<-TPMs.by_full_loci_and_SVs$P.9905.batchCT
TPMs.by_full_loci_and_SVs$F.9981<-apply(TPMs.by_full_loci_and_SVs[,grep("F.9981",names(TPMs.by_full_loci_and_SVs))],1,mean)
TPMs.by_full_loci_and_SVs$P.9981<-TPMs.by_full_loci_and_SVs$P.9981.batch2

TPMs.by_full_loci_and_SVs<-TPMs.by_full_loci_and_SVs[,grep("batch",names(TPMs.by_full_loci_and_SVs),invert = T)]
TPMs.by_full_loci_and_SVs<-TPMs.by_full_loci_and_SVs[,grep("rep",names(TPMs.by_full_loci_and_SVs),invert = T)]


TPMs.by_full_loci_and_SVs$max.flowers<- apply(TPMs.by_full_loci_and_SVs[,grep("F.",names(TPMs.by_full_loci_and_SVs))],1,max,na.rm=T)
TPMs.by_full_loci_and_SVs$max.pollen<- apply(TPMs.by_full_loci_and_SVs[,grep("P.",names(TPMs.by_full_loci_and_SVs))],1,max,na.rm=T)
TPMs.by_full_loci_and_SVs$max.rosette<- apply(TPMs.by_full_loci_and_SVs[,grep("R.",names(TPMs.by_full_loci_and_SVs))],1,max,na.rm=T)
TPMs.by_full_loci_and_SVs$max.seedl<- apply(TPMs.by_full_loci_and_SVs[,grep("S.",names(TPMs.by_full_loci_and_SVs))],1,max,na.rm=T)

TPMs.by_full_loci_and_SVs$mean.flowers<- apply(TPMs.by_full_loci_and_SVs[,grep("F.",names(TPMs.by_full_loci_and_SVs))],1,mean,na.rm=T)
TPMs.by_full_loci_and_SVs$mean.pollen<- apply(TPMs.by_full_loci_and_SVs[,grep("P.",names(TPMs.by_full_loci_and_SVs))],1,mean,na.rm=T)
TPMs.by_full_loci_and_SVs$mean.rosette<- apply(TPMs.by_full_loci_and_SVs[,grep("R.",names(TPMs.by_full_loci_and_SVs))],1,mean,na.rm=T)
TPMs.by_full_loci_and_SVs$mean.seedl<- apply(TPMs.by_full_loci_and_SVs[,grep("S.",names(TPMs.by_full_loci_and_SVs))],1,mean,na.rm=T)

TPMs.by_full_loci_and_SVs$Nacc_expressed.flowers<- apply(TPMs.by_full_loci_and_SVs[,grep("F.",names(TPMs.by_full_loci_and_SVs))],1,function(i) sum(i > 0.25,na.rm=T))
TPMs.by_full_loci_and_SVs$Nacc_expressed.pollen<- apply(TPMs.by_full_loci_and_SVs[,grep("P.",names(TPMs.by_full_loci_and_SVs))],1,function(i) sum(i > 0.25,na.rm=T))
TPMs.by_full_loci_and_SVs$Nacc_expressed.rosette<- apply(TPMs.by_full_loci_and_SVs[,grep("R.",names(TPMs.by_full_loci_and_SVs))],1,function(i) sum(i > 0.25,na.rm=T))
TPMs.by_full_loci_and_SVs$Nacc_expressed.seedl<- apply(TPMs.by_full_loci_and_SVs[,grep("S.",names(TPMs.by_full_loci_and_SVs))],1,function(i) sum(i > 0.25,na.rm=T))

TPMs.by_full_loci_and_SVs$freq<- apply(TPMs.by_full_loci_and_SVs[,grep("R.",names(TPMs.by_full_loci_and_SVs))],1,function(i) sum(!is.na(i)))


rownames(TPMs.by_full_loci_and_SVs)<-TPMs.by_full_loci_and_SVs$gene


#remove intermediate TPM tables for individual  accessions
rm(list = ls(pattern = "^a[0-9]+\\.TPMs\\.by_full_loci_and_SVs$"))
#########################################################

write.table(TPMs.by_full_loci_and_SVs,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/TPMs.genes_and_SVs.mapped_to_own_genomes.txt",quote = F, sep="\t",col.names = T,row.names = F)


########### counts on own genomes by full loci - genes and SVs
#########################################################
#on own genomes  by full loci 
i=1
acc=as.character(accessions_w_220011$V1[i])
a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/",acc,".counts.by_full_loci_and_SVs.bed",sep=""),header = T)
assign(paste("a",acc,".counts.by_full_loci_and_SVs",sep=""),a)
counts.by_full_loci_and_SVs<-a

for ( i in 2:27){
  acc=as.character(accessions_w_220011$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/TPMs/",acc,".counts.by_full_loci_and_SVs.bed",sep=""),header = T)
  assign(paste("a",acc,".counts.by_full_loci_and_SVs",sep=""),a)
  counts.by_full_loci_and_SVs<-merge(counts.by_full_loci_and_SVs,a,by="gene",all = T)
}

counts.by_full_loci_and_SVs.unprocessed<-counts.by_full_loci_and_SVs
# average those samples that have replicates, remove batch names from sample names
counts.by_full_loci_and_SVs$F.10002<-apply(counts.by_full_loci_and_SVs[,grep("F.10002",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$P.10002<-counts.by_full_loci_and_SVs$P.10002.batch2
counts.by_full_loci_and_SVs$F.10015<-apply(counts.by_full_loci_and_SVs[,grep("F.10015",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$P.10015<-apply(counts.by_full_loci_and_SVs[,grep("P.10015",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$S.10015<-apply(counts.by_full_loci_and_SVs[,grep("S.10015",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$R.10015<-apply(counts.by_full_loci_and_SVs[,grep("R.10015",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$F.10024<-apply(counts.by_full_loci_and_SVs[,grep("F.10024",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$P.10024<-counts.by_full_loci_and_SVs$P.10024.batch2
counts.by_full_loci_and_SVs$F.1741<-counts.by_full_loci_and_SVs$F.1741.batchCT
counts.by_full_loci_and_SVs$P.1741<-counts.by_full_loci_and_SVs$P.1741.batchCT

counts.by_full_loci_and_SVs$F.22007<-apply(counts.by_full_loci_and_SVs[,grep("F.22007",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$P.22007<-apply(counts.by_full_loci_and_SVs[,grep("P.22007",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$F.22006<-counts.by_full_loci_and_SVs$F.22006.batchCT
counts.by_full_loci_and_SVs$P.22006<-counts.by_full_loci_and_SVs$P.22006.batchCT
counts.by_full_loci_and_SVs$F.22005<-apply(counts.by_full_loci_and_SVs[,grep("F.22005",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$P.22005<-apply(counts.by_full_loci_and_SVs[,grep("P.22005",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$F.22004<-counts.by_full_loci_and_SVs$F.22004.batchCT
counts.by_full_loci_and_SVs$P.22004<-counts.by_full_loci_and_SVs$P.22004.batchCT
counts.by_full_loci_and_SVs$F.22003<-counts.by_full_loci_and_SVs$F.22003.batchCT
counts.by_full_loci_and_SVs$P.22003<-counts.by_full_loci_and_SVs$P.22003.batchCT
counts.by_full_loci_and_SVs$F.22002<-apply(counts.by_full_loci_and_SVs[,grep("F.22002",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$P.22002<-counts.by_full_loci_and_SVs$P.22002.batch2
counts.by_full_loci_and_SVs$F.220011<-apply(counts.by_full_loci_and_SVs[,grep("F.220011",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$P.220011<-counts.by_full_loci_and_SVs$P.220011.batch2

counts.by_full_loci_and_SVs$F.6024<-counts.by_full_loci_and_SVs$F.6024.batchCT
counts.by_full_loci_and_SVs$P.6024<-counts.by_full_loci_and_SVs$P.6024.batchCT
counts.by_full_loci_and_SVs$F.6244<-counts.by_full_loci_and_SVs$F.6244.batchCT
counts.by_full_loci_and_SVs$P.6244<-counts.by_full_loci_and_SVs$P.6244.batchCT
counts.by_full_loci_and_SVs$F.6909<-apply(counts.by_full_loci_and_SVs[,grep("F.6909",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$P.6909<-counts.by_full_loci_and_SVs$P.6909.batch2
counts.by_full_loci_and_SVs$F.6966<-apply(counts.by_full_loci_and_SVs[,grep("F.6966",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$P.6966<-counts.by_full_loci_and_SVs$P.6966.batch2
counts.by_full_loci_and_SVs$F.8236<-apply(counts.by_full_loci_and_SVs[,grep("F.8236",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$P.8236<-counts.by_full_loci_and_SVs$P.8236.batch2
counts.by_full_loci_and_SVs$F.9075<-counts.by_full_loci_and_SVs$F.9075.batch1
counts.by_full_loci_and_SVs$P.9075<-counts.by_full_loci_and_SVs$P.9075.batchCT
counts.by_full_loci_and_SVs$F.9537<-apply(counts.by_full_loci_and_SVs[,grep("F.9537",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$P.9537<-counts.by_full_loci_and_SVs$P.9537.batch2
counts.by_full_loci_and_SVs$F.9543<-counts.by_full_loci_and_SVs$F.9543.batchCT
counts.by_full_loci_and_SVs$P.9543<-counts.by_full_loci_and_SVs$P.9543.batchCT
counts.by_full_loci_and_SVs$F.9638<-counts.by_full_loci_and_SVs$F.9638.batchCT
counts.by_full_loci_and_SVs$P.9638<-counts.by_full_loci_and_SVs$P.9638.batchCT
counts.by_full_loci_and_SVs$F.9728<-counts.by_full_loci_and_SVs$F.9728.batchCT
counts.by_full_loci_and_SVs$P.9728<-counts.by_full_loci_and_SVs$P.9728.batchCT
counts.by_full_loci_and_SVs$F.9764<-apply(counts.by_full_loci_and_SVs[,grep("F.9764",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$P.9764<-apply(counts.by_full_loci_and_SVs[,grep("P.9764",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$F.9888<-counts.by_full_loci_and_SVs$F.9888.batchCT
counts.by_full_loci_and_SVs$P.9888<-counts.by_full_loci_and_SVs$P.9888.batchCT
counts.by_full_loci_and_SVs$F.9905<-counts.by_full_loci_and_SVs$F.9905.batchCT
counts.by_full_loci_and_SVs$P.9905<-counts.by_full_loci_and_SVs$P.9905.batchCT
counts.by_full_loci_and_SVs$F.9981<-apply(counts.by_full_loci_and_SVs[,grep("F.9981",names(counts.by_full_loci_and_SVs))],1,mean)
counts.by_full_loci_and_SVs$P.9981<-counts.by_full_loci_and_SVs$P.9981.batch2

counts.by_full_loci_and_SVs<-counts.by_full_loci_and_SVs[,grep("batch",names(counts.by_full_loci_and_SVs),invert = T)]
counts.by_full_loci_and_SVs<-counts.by_full_loci_and_SVs[,grep("rep",names(counts.by_full_loci_and_SVs),invert = T)]


counts.by_full_loci_and_SVs$max.flowers<- apply(counts.by_full_loci_and_SVs[,grep("F.",names(counts.by_full_loci_and_SVs))],1,max,na.rm=T)
counts.by_full_loci_and_SVs$max.pollen<- apply(counts.by_full_loci_and_SVs[,grep("P.",names(counts.by_full_loci_and_SVs))],1,max,na.rm=T)
counts.by_full_loci_and_SVs$max.rosette<- apply(counts.by_full_loci_and_SVs[,grep("R.",names(counts.by_full_loci_and_SVs))],1,max,na.rm=T)
counts.by_full_loci_and_SVs$max.seedl<- apply(counts.by_full_loci_and_SVs[,grep("S.",names(counts.by_full_loci_and_SVs))],1,max,na.rm=T)

counts.by_full_loci_and_SVs$mean.flowers<- apply(counts.by_full_loci_and_SVs[,grep("F.",names(counts.by_full_loci_and_SVs))],1,mean,na.rm=T)
counts.by_full_loci_and_SVs$mean.pollen<- apply(counts.by_full_loci_and_SVs[,grep("P.",names(counts.by_full_loci_and_SVs))],1,mean,na.rm=T)
counts.by_full_loci_and_SVs$mean.rosette<- apply(counts.by_full_loci_and_SVs[,grep("R.",names(counts.by_full_loci_and_SVs))],1,mean,na.rm=T)
counts.by_full_loci_and_SVs$mean.seedl<- apply(counts.by_full_loci_and_SVs[,grep("S.",names(counts.by_full_loci_and_SVs))],1,mean,na.rm=T)

counts.by_full_loci_and_SVs$Nacc_expressed.flowers<- apply(counts.by_full_loci_and_SVs[,grep("F.",names(counts.by_full_loci_and_SVs))],1,function(i) sum(i > 0.25,na.rm=T))
counts.by_full_loci_and_SVs$Nacc_expressed.pollen<- apply(counts.by_full_loci_and_SVs[,grep("P.",names(counts.by_full_loci_and_SVs))],1,function(i) sum(i > 0.25,na.rm=T))
counts.by_full_loci_and_SVs$Nacc_expressed.rosette<- apply(counts.by_full_loci_and_SVs[,grep("R.",names(counts.by_full_loci_and_SVs))],1,function(i) sum(i > 0.25,na.rm=T))
counts.by_full_loci_and_SVs$Nacc_expressed.seedl<- apply(counts.by_full_loci_and_SVs[,grep("S.",names(counts.by_full_loci_and_SVs))],1,function(i) sum(i > 0.25,na.rm=T))

counts.by_full_loci_and_SVs$freq<- apply(counts.by_full_loci_and_SVs[,grep("R.",names(counts.by_full_loci_and_SVs))],1,function(i) sum(!is.na(i)))


rownames(counts.by_full_loci_and_SVs)<-counts.by_full_loci_and_SVs$gene


#remove intermediate TPM tables for individual  accessions
rm(list = ls(pattern = "^a[0-9]+\\.counts\\.by_full_loci_and_SVs$"))
#########################################################

write.table(counts.by_full_loci_and_SVs,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/counts.genes_and_SVs.mapped_to_own_genomes.txt",quote = F, sep="\t",col.names = T,row.names = F)


############################################################################
# copies TPMs 

############################################################################




######################################################
#import EPIGENETICS 
######################################################

# import ChIPseq coverage 
#######################################################

#loci
a1741.loci.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.1741.loci.log2.bed", sep="",header = T)
a6024.loci.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.6024.loci.log2.bed", sep="",header = T)
a6909.loci.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.6909.loci.log2.bed", sep="",header = T)
a6966.loci.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.6966.loci.log2.bed", sep="",header = T)
a9888.loci.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.9888.loci.log2.bed", sep="",header = T)
a9905.loci.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.9905.loci.log2.bed", sep="",header = T)


chipseq.loci<-merge (a1741.loci.chip[,c(4,6:15)],a6024.loci.chip[,c(4,6:10,13:17)],by="group",all = T)
chipseq.loci<-merge (chipseq.loci,a6909.loci.chip[,c(4,6:10,13:17)],by="group",all = T)
chipseq.loci<-merge (chipseq.loci,a6966.loci.chip[,c(4,6:15)],by="group",all = T)
chipseq.loci<-merge (chipseq.loci,a9888.loci.chip[,c(4,6:15)],by="group",all = T)
chipseq.loci<-merge (chipseq.loci,a9905.loci.chip[,c(4,6:15)],by="group",all = T)


k27<-chipseq.loci[,grep("K27",names(chipseq.loci))]
k4<-chipseq.loci[,grep("K4",names(chipseq.loci))]
h1<-chipseq.loci[,grep("H1",names(chipseq.loci))]
k9<-chipseq.loci[,grep("K9",names(chipseq.loci))]
k36<-chipseq.loci[,grep("K36",names(chipseq.loci))]

cork27<-cor(k27,use = "complete.obs")
corh1<-cor(h1,use = "complete.obs")
cork9<-cor(k9,use = "complete.obs")
cork4<-cor(k4,use = "complete.obs")
cork36<-cor(k36,use = "complete.obs")

#K4 has the worst correlation

mean(apply(cork27,1,mean)) #[1] 0.7700716
mean(apply(corh1,1,mean)) #[1] 0.7742848
mean(apply(cork4,1,mean)) #[1] 0.5796799
mean(apply(cork9,1,mean)) #[1] 0.7941669
mean(apply(cork36,1,mean)) #[1] 0.8213551

cork27[apply(cork27,1,mean)<0.7,]
cork27[apply(cork27,1,mean)<0.5,]
cork9[apply(cork27,1,mean)<0.7,] 
cork9[apply(cork27,1,mean)<0.5,]
cork4[apply(cork27,1,mean)<0.7,] 
cork4[apply(cork27,1,mean)<0.5,]

# exclude: 
# 9888.rep2.H3K4me3
#  9888.rep2.H3K27me3
#6024.rep2.H3K4me3 
#X6966.rep1.H2K27me3 - seems to be contaminated with K9 - 6966.rep1.H2K27me3


a1741.mRNAs.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.1741.mRNAs.log2.bed", sep="",header = T)
a6024.mRNAs.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.6024.mRNAs.log2.bed", sep="",header = T)
a6909.mRNAs.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.6909.mRNAs.log2.bed", sep="",header = T)
a6966.mRNAs.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.6966.mRNAs.log2.bed", sep="",header = T)
a9888.mRNAs.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.9888.mRNAs.log2.bed", sep="",header = T)
a9905.mRNAs.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.9905.mRNAs.log2.bed", sep="",header = T)



chipseq.mRNAs<-merge (a1741.mRNAs.chip[,c(4,6:15)],a6024.mRNAs.chip[,c(4,6:10,13:17)],by="group",all = T)
chipseq.mRNAs<-merge (chipseq.mRNAs,a6909.mRNAs.chip[,c(4,6:10,13:17)],by="group",all = T)
chipseq.mRNAs<-merge (chipseq.mRNAs,a6966.mRNAs.chip[,c(4,6:15)],by="group",all = T)
chipseq.mRNAs<-merge (chipseq.mRNAs,a9888.mRNAs.chip[,c(4,6:15)],by="group",all = T)
chipseq.mRNAs<-merge (chipseq.mRNAs,a9905.mRNAs.chip[,c(4,6:15)],by="group",all = T)


k27<-chipseq.mRNAs[,grep("K27",names(chipseq.mRNAs))]
k4<-chipseq.mRNAs[,grep("K4",names(chipseq.mRNAs))]
h1<-chipseq.mRNAs[,grep("H1",names(chipseq.mRNAs))]
k9<-chipseq.mRNAs[,grep("K9",names(chipseq.mRNAs))]
k36<-chipseq.mRNAs[,grep("K36",names(chipseq.mRNAs))]

cork27<-cor(k27,use = "complete.obs")
corh1<-cor(h1,use = "complete.obs")
cork9<-cor(k9,use = "complete.obs")
cork4<-cor(k4,use = "complete.obs")
cork36<-cor(k36,use = "complete.obs")



mean(apply(cork27,1,mean)) #[1] 0.8194177
mean(apply(corh1,1,mean)) #[1] 0.800482
mean(apply(cork4,1,mean)) #[1] 0.6059178
mean(apply(cork9,1,mean)) #[1] 0.7730113
mean(apply(cork36,1,mean)) #[1] 0.8488054

cork27[apply(cork27,1,mean)<0.7,] 
corh1[apply(corh1,1,mean)<0.7,] 
cork9[apply(cork9,1,mean)<0.7,] 
cork4[apply(cork4,1,mean)<0.7,] 
cork36[apply(cork36,1,mean)<0.7,] 


a1741.SVs.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.1741.svs_v06.log2.bed", sep="",header = T)
a6024.SVs.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.6024.svs_v06.log2.bed", sep="",header = T)
a6909.SVs.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.6909.svs_v06.log2.bed", sep="",header = T)
a6966.SVs.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.6966.svs_v06.log2.bed", sep="",header = T)
a9888.SVs.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.9888.svs_v06.log2.bed", sep="",header = T)
a9905.SVs.chip <- read.csv("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/chipseq/Chipseq_coverage.owngenome.9905.svs_v06.log2.bed", sep="",header = T)

chipseq.SVs<-merge (a1741.SVs.chip[,c(4,6:15)],a6024.SVs.chip[,c(4,6:10,13:17)],by="SV_name",all = T)
chipseq.SVs<-merge (chipseq.SVs,a6909.SVs.chip[,c(4,6:10,13:17)],by="SV_name",all = T)
chipseq.SVs<-merge (chipseq.SVs,a6966.SVs.chip[,c(4,6:15)],by="SV_name",all = T)
chipseq.SVs<-merge (chipseq.SVs,a9888.SVs.chip[,c(4,6:15)],by="SV_name",all = T)
chipseq.SVs<-merge (chipseq.SVs,a9905.SVs.chip[,c(4,6:15)],by="SV_name",all = T)

chipseq.SVs$group<-chipseq.SVs$SV_name
# exclude all K4 
# exclude 9888.rep2.H3K27me3 
# exclude 6966.rep1.H2K27me3 (may be contaminated with K9) 



chipseq.mRNAs<-chipseq.mRNAs[,grep("K4",names(chipseq.mRNAs),invert = T)]
chipseq.mRNAs<-chipseq.mRNAs[,grep("9888.rep2.H3K27me3",names(chipseq.mRNAs),invert = T)]
chipseq.mRNAs<-chipseq.mRNAs[,grep("6966.rep1.H2K27me3",names(chipseq.mRNAs),invert = T)]

chipseq.loci<-chipseq.loci[,grep("K4",names(chipseq.loci),invert = T)]
chipseq.loci<-chipseq.loci[,grep("9888.rep2.H3K27me3",names(chipseq.loci),invert = T)]
chipseq.loci<-chipseq.loci[,grep("6966.rep1.H2K27me3",names(chipseq.loci),invert = T)]

chipseq.SVs<-chipseq.SVs[,grep("K4",names(chipseq.SVs),invert = T)]
chipseq.SVs<-chipseq.SVs[,grep("9888.rep2.H3K27me3",names(chipseq.SVs),invert = T)]
chipseq.SVs<-chipseq.SVs[,grep("6966.rep1.H2K27me3",names(chipseq.SVs),invert = T)]


chipseq.mRNAs$add<-chipseq.mRNAs$group
chipseq.loci$add<-chipseq.loci$group

# exclude: 
# 9888.rep2.H3K4me3
#  9888.rep2.H3K27me3
#6024.rep2.H3K4me3 
#X6966.rep1.H2K27me3
# and H3K4me3 (K4)

#######################################
#define chipseq table processing function
#####################################################################
processchipseq_table <- function(a){
  a<-a[!duplicated(a$gene),]
  rownames(a)<-a$gene
  a$H1.1741<-apply(a[,c( "X1741.rep1.H1" ,"X1741.rep2.H1"  )],1,mean)
  #a$K4.1741<-apply(a[,c( "X1741.rep1.H3K4me3" ,"X1741.rep2.H3K4me3"  )],1,mean)
  a$K9.1741  <-apply(a[,c( "X1741.rep1.H3K9me2" ,"X1741.rep2.H3K9me2"  )],1,mean)
  a$K27.1741<-apply(a[,c( "X1741.rep1.H3K27me3" ,"X1741.rep2.H3K27me3"  )],1,mean)
  a$K36.1741<-apply(a[,c( "X1741.rep1.H3K36me3" ,"X1741.rep2.H3K36me3"  )],1,mean)
    #X6909
  a$H1.6909<-apply(a[,c( "X6909.rep1.H1" ,"X6909.rep2.H1"  )],1,mean)
 # a$K4.6909<-apply(a[,c( "X6909.rep1.H3K4me3" ,"X6909.rep2.H3K4me3"  )],1,mean)
  a$K9.6909  <-apply(a[,c( "X6909.rep1.H3K9me2" ,"X6909.rep2.H3K9me2"  )],1,mean)
  a$K27.6909<-apply(a[,c( "X6909.rep1.H3K27me3" ,"X6909.rep2.H3K27me3"  )],1,mean)
  a$K36.6909<-apply(a[,c( "X6909.rep1.H3K36me3" ,"X6909.rep2.H3K36me3"  )],1,mean)
   #X6966
  a$H1.6966<-apply(a[,c( "X6966.rep1.H1" ,"X6966.rep2.H1"  )],1,mean)
 # a$K4.6966<-apply(a[,c( "X6966.rep1.H3K4me3" ,"X6966.rep2.H3K4me3"  )],1,mean)
  a$K9.6966  <-apply(a[,c( "X6966.rep1.H3K9me2" ,"X6966.rep2.H3K9me2"  )],1,mean)
  #a$K27.6966<-a[,c( "X6966.rep1.H3K27me3" )]
  a$K36.6966<-apply(a[,c( "X6966.rep1.H3K36me3" ,"X6966.rep2.H3K36me3"  )],1,mean)
   #X17414
  a$H1.9888<-apply(a[,c( "X9888.rep1.H1" ,"X9888.rep2.H1"  )],1,mean)
  #a$K4.9888<-a[,c( "X9888.rep1.H3K4me3" )]
  a$K9.9888  <-apply(a[,c( "X9888.rep1.H3K9me2" ,"X9888.rep2.H3K9me2"  )],1,mean)
  a$K27.9888<-a[,c( "X9888.rep1.H3K27me3"  )]
  a$K36.9888<-apply(a[,c( "X9888.rep1.H3K36me3" ,"X9888.rep2.H3K36me3"  )],1,mean)
  #X17415
  a$H1.9905<-apply(a[,c( "X9905.rep1.H1" ,"X9905.rep2.H1"  )],1,mean)
  #a$K4.9905<-apply(a[,c( "X9905.rep1.H3K4me3" ,"X9905.rep2.H3K4me3"  )],1,mean)
  a$K9.9905  <-apply(a[,c( "X9905.rep1.H3K9me2" ,"X9905.rep2.H3K9me2"  )],1,mean)
  a$K27.9905<-apply(a[,c( "X9905.rep1.H3K27me3" ,"X9905.rep2.H3K27me3"  )],1,mean)
  a$K36.9905<-apply(a[,c( "X9905.rep1.H3K36me3" ,"X9905.rep2.H3K36me3"  )],1,mean)
   #X6024
  a$H1.6024<-apply(a[,c( "X6024.rep1.H1" ,"X6024.rep2.H1"  )],1,mean)
  #a$K4.6024<-a[,c( "X6024.rep1.H3K4me3" )]
  a$K9.6024  <-apply(a[,c( "X6024.rep1.H3K9me2" ,"X6024.rep2.H3K9me2"  )],1,mean)
  a$K27.6024<-apply(a[,c( "X6024.rep1.H3K27me3" ,"X6024.rep2.H3K27me3"  )],1,mean)
  a$K36.6024<-apply(a[,c( "X6024.rep1.H3K36me3" ,"X6024.rep2.H3K36me3"  )],1,mean)
  a<-a[,c(50:72)]
  return(a)
}
###################################################

# custom function to implement min max scaling
#################################################################
minMax <- function(x) {  (x - min(x)) / (max(x) - min(x))}
quantile80minmax <- function(x) {  (x - quantile(x,.20,na.rm=T)) / (quantile(x,.80,na.rm=T) - quantile(x,.20,na.rm=T))}
##################################################################

# quantile normalize 
#########################################################
chipseq.mRNAs.quantstan <- cbind(chipseq.mRNAs$group,as.data.frame(lapply(chipseq.mRNAs[,2:48], quantile80minmax)))
chipseq.mRNAs.quantstan$gene<-chipseq.mRNAs.quantstan[,1]
rownames(chipseq.mRNAs.quantstan)<-chipseq.mRNAs.quantstan$gene
chipseq.mRNAs.quantstan<-processchipseq_table(chipseq.mRNAs.quantstan)
chipseq.mRNAs.quantstan$gene<-rownames(chipseq.mRNAs.quantstan)



chipseq.SVs.quantstan <- cbind(chipseq.SVs$group,as.data.frame(lapply(chipseq.SVs[,2:48], quantile80minmax)))
chipseq.SVs.quantstan$gene<-chipseq.SVs.quantstan[,1]
rownames(chipseq.SVs.quantstan)<-chipseq.SVs.quantstan$gene
chipseq.SVs.quantstan<-processchipseq_table(chipseq.SVs.quantstan)
chipseq.SVs.quantstan$gene<-rownames(chipseq.SVs.quantstan)


chipseq.loci.quantstan <- cbind(chipseq.loci$group,as.data.frame(lapply(chipseq.loci[,2:48], quantile80minmax)))
chipseq.loci.quantstan$gene<-chipseq.loci.quantstan[,1]
rownames(chipseq.loci.quantstan)<-chipseq.loci.quantstan$gene
chipseq.loci.quantstan<-processchipseq_table(chipseq.loci.quantstan)
chipseq.loci.quantstan$gene<-rownames(chipseq.loci.quantstan)




#remove intermediate count tables for individual  accessions
rm(list = ls(pattern = "^a.*\\.chip$"))

####################################################


write.table(chipseq.loci.quantstan,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/ChIPseq_coverage.input_and_quantile_normalized.6accessions.genes.txt",quote = F, sep="\t",col.names = T,row.names = F)

write.table(chipseq.SVs.quantstan,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/ChIPseq_coverage.input_and_quantile_normalized.6accessions.svs_v06.txt",quote = F, sep="\t",col.names = T,row.names = F)


######################################################
#import METHYLATION 
######################################################

#import METHYLATION  loci (genes) on own genomes
############################################
i=1
acc=as.character(accessions_methylation$V1[i])
a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CG.own.loci.bed",sep=""),header = T)
names(a)<-c("chr","start","end","gene","score","strand",acc)
assign(paste("a",acc,".loci.CG",sep=""),a)
CG.loci<-a[,c("gene",acc)]

a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CHG.own.loci.bed",sep=""),header = T)
names(a)<-c("chr","start","end","gene","score","strand",acc)
assign(paste("a",acc,".loci.CHG",sep=""),a)
CHG.loci<-a[,c("gene",acc)]

a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CHH.own.loci.bed",sep=""),header = T)
names(a)<-c("chr","start","end","gene","score","strand",acc)
assign(paste("a",acc,".loci.CHH",sep=""),a)
CHH.loci<-a[,c("gene",acc)]

for ( i in 2:12){
  acc=as.character(accessions_methylation$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CG.own.loci.bed",sep=""),header = T)
  names(a)<-c("chr","start","end","gene","score","strand",acc)
  assign(paste("a",acc,".loci.CG",sep=""),a)
  CG.loci<-merge (CG.loci,a[,c("gene",acc)],by="gene",all = T)
  
  acc=as.character(accessions_methylation$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CHG.own.loci.bed",sep=""),header = T)
  names(a)<-c("chr","start","end","gene","score","strand",acc)
  assign(paste("a",acc,".loci.CHG",sep=""),a)
  CHG.loci<-merge (CHG.loci,a[,c("gene",acc)],by="gene",all = T)
  
  acc=as.character(accessions_methylation$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CHH.own.loci.bed",sep=""),header = T)
  names(a)<-c("chr","start","end","gene","score","strand",acc)
  assign(paste("a",acc,".loci.CHH",sep=""),a)
  CHH.loci<-merge (CHH.loci,a[,c("gene",acc)],by="gene",all = T)
  }
############################################

# import METHYLATION SVs on own genomes
############################################
i=1
acc=as.character(accessions_methylation$V1[i])
a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CG.own.SVs.bed",sep=""),header = T)
names(a)<-c("chr","start","end","gene","SV_type","SV_info",acc)
a<-a[a$SV_type!="multi_sv",]
a$length<-lapply(strsplit(as.character(a$SV_info),"len=", fixed = T), "[",2)
a$length<-lapply(strsplit(as.character(a$length),";", fixed = T), "[",1)
a$te_cont<-lapply(strsplit(as.character(a$SV_info),"te.content=", fixed = T), "[",2)
a$te_cont<-lapply(strsplit(as.character(a$te_cont),";", fixed = T), "[",1)
a$te_family<-lapply(strsplit(as.character(a$SV_info),"te.family=", fixed = T), "[",2)
a$te_family<-lapply(strsplit(as.character(a$te_family),";", fixed = T), "[",1)
assign(paste("a",acc,".SVs.CG",sep=""),a)
CG.SVs<-a[,c("gene",acc)]

a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CHG.own.SVs.bed",sep=""),header = T)
names(a)<-c("chr","start","end","gene","SV_type","SV_info",acc)
a<-a[a$SV_type!="multi_sv",]
a$length<-lapply(strsplit(as.character(a$SV_info),"len=", fixed = T), "[",2)
a$length<-lapply(strsplit(as.character(a$length),";", fixed = T), "[",1)
a$te_cont<-lapply(strsplit(as.character(a$SV_info),"te.content=", fixed = T), "[",2)
a$te_cont<-lapply(strsplit(as.character(a$te_cont),";", fixed = T), "[",1)
a$te_family<-lapply(strsplit(as.character(a$SV_info),"te.family=", fixed = T), "[",2)
a$te_family<-lapply(strsplit(as.character(a$te_family),";", fixed = T), "[",1)
assign(paste("a",acc,".SVs.CHG",sep=""),a)
CHG.SVs<-a[,c("gene",acc)]

a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CHH.own.SVs.bed",sep=""),header = T)
names(a)<-c("chr","start","end","gene","SV_type","SV_info",acc)
a<-a[a$SV_type!="multi_sv",]
a$length<-lapply(strsplit(as.character(a$SV_info),"len=", fixed = T), "[",2)
a$length<-lapply(strsplit(as.character(a$length),";", fixed = T), "[",1)
a$te_cont<-lapply(strsplit(as.character(a$SV_info),"te.content=", fixed = T), "[",2)
a$te_cont<-lapply(strsplit(as.character(a$te_cont),";", fixed = T), "[",1)
a$te_family<-lapply(strsplit(as.character(a$SV_info),"te.family=", fixed = T), "[",2)
a$te_family<-lapply(strsplit(as.character(a$te_family),";", fixed = T), "[",1)
assign(paste("a",acc,".SVs.CHH",sep=""),a)
CHH.SVs<-a[,c("gene",acc)]

for ( i in 2:12){
  acc=as.character(accessions_methylation$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CG.own.SVs.bed",sep=""),header = T)
  names(a)<-c("chr","start","end","gene","SV_type","SV_info",acc)
  a<-a[a$SV_type!="multi_sv",]
  a$length<-lapply(strsplit(as.character(a$SV_info),"len=", fixed = T), "[",2)
  a$length<-lapply(strsplit(as.character(a$length),";", fixed = T), "[",1)
  a$te_cont<-lapply(strsplit(as.character(a$SV_info),"te.content=", fixed = T), "[",2)
  a$te_cont<-lapply(strsplit(as.character(a$te_cont),";", fixed = T), "[",1)
  a$te_family<-lapply(strsplit(as.character(a$SV_info),"te.family=", fixed = T), "[",2)
  a$te_family<-lapply(strsplit(as.character(a$te_family),";", fixed = T), "[",1)
  assign(paste("a",acc,".SVs.CG",sep=""),a)
  CG.SVs<-merge (CG.SVs,a[,c("gene",acc)],by="gene",all = T)
  
  acc=as.character(accessions_methylation$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CHG.own.SVs.bed",sep=""),header = T)
  names(a)<-c("chr","start","end","gene","SV_type","SV_info",acc)
  a<-a[a$SV_type!="multi_sv",]
  a$length<-lapply(strsplit(as.character(a$SV_info),"len=", fixed = T), "[",2)
  a$length<-lapply(strsplit(as.character(a$length),";", fixed = T), "[",1)
  a$te_cont<-lapply(strsplit(as.character(a$SV_info),"te.content=", fixed = T), "[",2)
  a$te_cont<-lapply(strsplit(as.character(a$te_cont),";", fixed = T), "[",1)
  a$te_family<-lapply(strsplit(as.character(a$SV_info),"te.family=", fixed = T), "[",2)
  a$te_family<-lapply(strsplit(as.character(a$te_family),";", fixed = T), "[",1)
  assign(paste("a",acc,".SVs.CHG",sep=""),a)
  CHG.SVs<-merge (CHG.SVs,a[,c("gene",acc)],by="gene",all = T)
  
  acc=as.character(accessions_methylation$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CHH.own.SVs.bed",sep=""),header = T)
  names(a)<-c("chr","start","end","gene","SV_type","SV_info",acc)
  a<-a[a$SV_type!="multi_sv",]
  a$length<-lapply(strsplit(as.character(a$SV_info),"len=", fixed = T), "[",2)
  a$length<-lapply(strsplit(as.character(a$length),";", fixed = T), "[",1)
  a$te_cont<-lapply(strsplit(as.character(a$SV_info),"te.content=", fixed = T), "[",2)
  a$te_cont<-lapply(strsplit(as.character(a$te_cont),";", fixed = T), "[",1)
  a$te_family<-lapply(strsplit(as.character(a$SV_info),"te.family=", fixed = T), "[",2)
  a$te_family<-lapply(strsplit(as.character(a$te_family),";", fixed = T), "[",1)
  assign(paste("a",acc,".SVs.CHH",sep=""),a)
  CHH.SVs<-merge (CHH.SVs,a[,c("gene",acc)],by="gene",all = T)
}

############################################

rm(list = ls(pattern = "^a.*\\.CHH$"))
rm(list = ls(pattern = "^a.*\\.CHG$"))
rm(list = ls(pattern = "^a.*\\.CG$"))

write.table(CG.loci,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/CG.methylation_level.12accessions.genes.txt",quote = F, sep="\t",col.names = T,row.names = F)
write.table(CHG.loci,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/CHG.methylation_level.12accessions.genes.txt",quote = F, sep="\t",col.names = T,row.names = F)

write.table(CHH.loci,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/CHH.methylation_level.12accessions.genes.txt",quote = F, sep="\t",col.names = T,row.names = F)

write.table(CG.SVs,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/CG.methylation_level.12accessions.SVs.txt",quote = F, sep="\t",col.names = T,row.names = F)

write.table(CHG.SVs,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/CHG.methylation_level.12accessions.SVs.txt",quote = F, sep="\t",col.names = T,row.names = F)

write.table(CHH.SVs,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/CHH.methylation_level.12accessions.SVs.txt",quote = F, sep="\t",col.names = T,row.names = F)


#import METHYLATION  loci (genes) on TAIR10
############################################
i=1
acc=as.character(accessions_methylation$V1[i])
a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CG.TAIR10.loci.bed",sep=""),header = T)
names(a)<-c("chr","start","end","gene","score","strand",acc)
assign(paste("a",acc,".TAIR10.loci.CG",sep=""),a)
CG.TAIR10.loci<-a[,c("gene",acc)]

a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CHG.TAIR10.loci.bed",sep=""),header = T)
names(a)<-c("chr","start","end","gene","score","strand",acc)
assign(paste("a",acc,".TAIR10.loci.CHG",sep=""),a)
CHG.TAIR10.loci<-a[,c("gene",acc)]

a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CHH.TAIR10.loci.bed",sep=""),header = T)
names(a)<-c("chr","start","end","gene","score","strand",acc)
assign(paste("a",acc,".TAIR10.loci.CHH",sep=""),a)
CHH.TAIR10.loci<-a[,c("gene",acc)]

for ( i in 2:12){
  acc=as.character(accessions_methylation$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CG.TAIR10.loci.bed",sep=""),header = T)
  names(a)<-c("chr","start","end","gene","score","strand",acc)
  assign(paste("a",acc,".TAIR10.loci.CG",sep=""),a)
  CG.TAIR10.loci<-merge (CG.TAIR10.loci,a[,c("gene",acc)],by="gene",all = T)
  
  acc=as.character(accessions_methylation$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CHG.TAIR10.loci.bed",sep=""),header = T)
  names(a)<-c("chr","start","end","gene","score","strand",acc)
  assign(paste("a",acc,".TAIR10.loci.CHG",sep=""),a)
  CHG.TAIR10.loci<-merge (CHG.TAIR10.loci,a[,c("gene",acc)],by="gene",all = T)
  
  acc=as.character(accessions_methylation$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/methylation/",acc,".CHH.TAIR10.loci.bed",sep=""),header = T)
  names(a)<-c("chr","start","end","gene","score","strand",acc)
  assign(paste("a",acc,".TAIR10.loci.CHH",sep=""),a)
  CHH.TAIR10.loci<-merge (CHH.TAIR10.loci,a[,c("gene",acc)],by="gene",all = T)
}
############################################





############## import small RNA 
#####################################################

#  import small RNA  - mRNAs 
############################################
i=1
acc=as.character(accessions_miRNA$V1[i])
a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/mirna/",acc,".24nt.mRNAs.coverage.perbase_calc.bed",sep=""),header = F)
names(a)<-c("chr","start","end","gene","score","strand",acc)
assign(paste("a",acc,".mRNAs.miRNA",sep=""),a)
sRNA.24nt.coverage.mRNAs<-a[,c("gene",acc)]

for ( i in 2:14){
  acc=as.character(accessions_miRNA$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/mirna/",acc,".24nt.mRNAs.coverage.perbase_calc.bed",sep=""),header = F)
  names(a)<-c("chr","start","end","gene","score","strand",acc)
  assign(paste("a",acc,".mRNAs.miRNA",sep=""),a)
  sRNA.24nt.coverage.mRNAs<-merge (sRNA.24nt.coverage.mRNAs,a[,c("gene",acc)],by="gene",all = T)
  
}
############################################

# import small RNA  - loci
############################################
i=1
acc=as.character(accessions_miRNA$V1[i])
a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/mirna/",acc,".24nt.loci.coverage.perbase_calc.bed",sep=""),header = F)
names(a)<-c("chr","start","end","gene","score","strand",acc)
assign(paste("a",acc,".loci.miRNA",sep=""),a)
sRNA.24nt.coverage.loci<-a[,c("gene",acc)]

for ( i in 2:14){
  acc=as.character(accessions_miRNA$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/mirna/",acc,".24nt.loci.coverage.perbase_calc.bed",sep=""),header = F)
  names(a)<-c("chr","start","end","gene","score","strand",acc)
  assign(paste("a",acc,".loci.miRNA",sep=""),a)
  sRNA.24nt.coverage.loci<-merge (sRNA.24nt.coverage.loci,a[,c("gene",acc)],by="gene",all = T)
  
}
############################################

# import small RNA  - SVs
############################################
i=1
acc=as.character(accessions_miRNA$V1[i])
a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/mirna/",acc,".24nt.SVs.coverage.perbase_calc.bed",sep=""),header = F)
names(a)<-c("chr","start","end","gene","SV_type","SV_info ",acc)
assign(paste("a",acc,".SVs.miRNA",sep=""),a)
sRNA.24nt.coverage.SVs<-a[,c("gene",acc)]

for ( i in 2:14){
  acc=as.character(accessions_miRNA$V1[i])
  a<-read.delim(paste("Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/coverage_tables/mirna/",acc,".24nt.SVs.coverage.perbase_calc.bed",sep=""),header = F)
  names(a)<-c("chr","start","end","gene","SV_type","SV_info",acc)
  assign(paste("a",acc,".SVs.miRNA",sep=""),a)
  sRNA.24nt.coverage.SVs<-merge (sRNA.24nt.coverage.SVs,a[,c("gene",acc)],by="gene",all = T)
  
}


############################################



write.table(sRNA.24nt.coverage.loci,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/sRNA.24nt.coverage.14accessions.genes.txt",quote = F, sep="\t",col.names = T,row.names = F)


write.table(sRNA.24nt.coverage.SVs,"Z:/01_POSTDOC/03_Projects/ERA-CAPS/001_paper_analyses/data_for_submission/sRNA.24nt.coverage.14accessions.SVs.txt",quote = F, sep="\t",col.names = T,row.names = F)






