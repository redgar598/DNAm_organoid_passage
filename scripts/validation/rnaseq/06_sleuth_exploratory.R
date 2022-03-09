#'---
#'title: Sleuth explore data
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ### Load Libraries
suppressMessages({
  library(sleuth)
  library(testit)
  library(ggplot2)
  library(RColorBrewer)
  library(reshape2)
  library(here)
  library(rafalib)
  library(grid)
})
options(stringsAsFactors = FALSE)
source(here("scripts","00_pretty_plots.R"))
source(here("scripts","00_heat_scree_plot_generic.R"))


#set input and output dirs
datapath = "data/validation/kallisto"
resultdir = here('data/validation/kallisto/sleuth')

sampleinfo <- read.table("data/validation/sample_info_RNA_seq.txt", header=F, sep=" ", skip=1)
sampleinfo<-sampleinfo[,c(2,7,8:10,12)]
colnames(sampleinfo)<-c("sample","concentration","volume","quantity","ratio","well")

#create a sample to condition metadata description
sample_id = sampleinfo$sample
kal_dirs <- file.path(datapath, sample_id)

sampleinfo$path<-kal_dirs


sampleinfo$individual<-sapply(1:nrow(sampleinfo), function(x) strsplit(sampleinfo$sample[x],"_")[[1]][1])
sampleinfo$sample.site<-sapply(1:nrow(sampleinfo), function(x) strsplit(sampleinfo$sample[x],"_")[[1]][2])
sampleinfo$passage<-as.numeric(sapply(1:nrow(sampleinfo), function(x) gsub("p","",strsplit(sampleinfo$sample[x],"_")[[1]][3])))
sampleinfo$passage_hilo<-sapply(1:nrow(sampleinfo), function(x) if(sampleinfo$passage[x]<5){"low"}else{"high"})
sampleinfo$condition<-sapply(1:nrow(sampleinfo), function(x) strsplit(sampleinfo$sample[x],"_")[[1]][4])
sampleinfo$condition<-as.factor(sampleinfo$condition)
levels(sampleinfo$condition)<-c("D","IFNg","IFNg","TNFa","UD","UT")
sampleinfo$condition<-as.character(sampleinfo$condition)
sampleinfo$comparison<-sapply(1:nrow(sampleinfo), function(x) if(sampleinfo$condition[x]%in%c("UD","D")){"differentiation"}else{"cytokine"})
sampleinfo$treatment<-sapply(1:nrow(sampleinfo), function(x) if(sampleinfo$comparison[x]=="cytokine"){sampleinfo$condition[x]}else{"UT"})
sampleinfo$differentiation<-sapply(1:nrow(sampleinfo), function(x) if(sampleinfo$comparison[x]=="differentiation"){sampleinfo$condition[x]}else{"UD"})




###########
#'# aggregate to gene level
###########
mart <- biomaRt::useMart(biomart = "ensembl",
                         dataset = "hsapiens_gene_ensembl")

ttg <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "transcript_version",
                 "ensembl_gene_id", "external_gene_name", "description",
                 "transcript_biotype"),
  mart = mart)
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id', 'ens_gene', 'ext_gene'))
head(ttg)


#'# run sleuth on the data
so <- sleuth_prep(sampleinfo, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)


#'## PCA my way
mat <- sleuth:::spread_abundance_by(so$obs_norm, "tpm", so$sample_to_covariates$sample)

                              # 
                              # ## ignore outlier
                              # mat_noOutlier<-mat[,which(colnames(mat)!="ERR2271027")]
                              # sampleinfo_noOutlier<-sampleinfo[which(sampleinfo$sample!="ERR2271027"),]
                              # 
                              # save(mat_noOutlier,sampleinfo_noOutlier, file=here("rna_seq/data","ibd_RNAseq_sleuth.RData"))

pca_res <- prcomp(t(mat))
Loadings_meta<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)


#'# PC vs PC plot
Loadings_meta$sample<-rownames(Loadings_meta)
Loadings_meta<-merge(Loadings_meta, sampleinfo, by="sample")

ggplot(Loadings_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=3, color="black")+
  theme_bw()+theme(axis.text = element_text(size=12),
                   axis.title = element_text(size=14),
                   plot.margin = margin(1, 0.1, 1, 0.5, "cm"))+
  fillscale_sampsite+th+xlab("PC1 (42%)")+ylab("PC2 (22%)")
ggsave(here("figs","RNAseq_PC1_PC2_samplesite.pdf"), width = 5, height = 4)
ggsave(here("figs/jpeg","RNAseq_PC1_PC2_samplesite.jpeg"), width = 5, height = 4)

ggplot(Loadings_meta, aes(PC1, PC2, fill=passage_hilo))+geom_point(shape=21,size=3)+theme_bw()+
  th+scale_fill_manual(values=c("#3288BD","#D53E4F"), name="Passage")+xlab("PC1 (42%)")+ylab("PC2 (22%)")
ggsave(here("figs","RNAseq_PC1_PC2_passage.pdf"), width = 7.5, height = 6)
ggsave(here("figs/jpeg","RNAseq_PC1_PC2_passage.jpeg"), width = 7.5, height = 6)

ggplot(Loadings_meta, aes(PC2, PC3, fill=passage_hilo))+
  geom_line(aes(PC2,PC3, group=individual, color=individual))+theme_bw()+
  geom_point(shape=21,size=3)+
  scale_color_manual(values=c("#d9d9d9","#525252","#969696","#737373","#bdbdbd"),name="Individual")+
  th+scale_fill_manual(values=c("#3288BD","#D53E4F"), name="Passage")+
  xlab("PC2 (22%)")+ylab("PC3 (12%)")+theme(axis.text = element_text(size=12),
                                            axis.title = element_text(size=14),
                                            plot.margin = margin(0.6, 1, 0.6, 0.6, "cm"))
ggsave(here("figs","RNAseq_PC2_PC3_passage.pdf"), width = 5.5, height = 4)
ggsave(here("figs/jpeg","RNAseq_PC2_PC3_passage.jpeg"), width = 5.5, height = 4)

cor(Loadings_meta$PC2, Loadings_meta$passage, method="spearman")


ggplot(Loadings_meta, aes(PC1, PC2, fill=treatment, color=differentiation))+geom_point(shape=21,size=3)+theme_bw()+th+
  scale_color_manual(values=c("black","white"))+
  scale_fill_manual(values=c("cornflowerblue","firebrick4","grey80"), name="Treatment")+xlab("PC1 (42%)")+ylab("PC2 (22%)")
ggsave(here("figs","RNAseq_PC1_PC2_condtions.pdf"), width = 7.5, height = 6)
ggsave(here("figs/jpeg","RNAseq_PC1_PC2_conditions.jpeg"), width = 7.5, height = 6)


#Restructure meta
meta_categorical <- sampleinfo[, c(8,9,11,14,15)]  # input column numbers in meta that contain categorical variables
meta_continuous <- sampleinfo[, c(2,3,4,5,10)] # input column numbers in meta that contain continuous variables


ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

Loadings_meta<-as.data.frame(pca_res$x)
suppressWarnings(heat_scree_plot(Loadings_meta, Importance, 3.3, 1.5))

ggsave(here("figs","RNAseq_heat_scree.pdf"), suppressWarnings(heat_scree_plot(Loadings_meta, Importance, 3.3, 1.5)),width = 9, height = 6)
ggsave(here("figs/jpeg","RNAseq_heat_scree.jpeg"), suppressWarnings(heat_scree_plot(Loadings_meta, Importance, 3.3, 1.5)),width = 9, height = 6)



#############
#'# clustering
#############
d <- dist(t(mat))
hc <- hclust(d, method = "complete") #single, complete, average, ward
identical(colnames(mat),sampleinfo$sample)

pdf(here('figs','RNAseq_clustering_condition.pdf'), width=15, h=10)
myplclust(hc, labels=sampleinfo$sample, lab.col=as.fumeric(sampleinfo$differentiation), cex=1)
dev.off()

pdf(here('figs','RNAseq_clustering_passage.pdf'), width=15, h=10)
myplclust(hc, labels=sampleinfo$sample, lab.col=as.fumeric(sampleinfo$passage_hilo), cex=1)
dev.off()





#### #### ####
#### cluster each condition seperate
#### #### ####
sampleinfo_differentiation<-sampleinfo[which(sampleinfo$comparison=="differentiation"),]
mat_differentiation<-mat[,which(colnames(mat)%in%sampleinfo_differentiation$sample)]

identical(colnames(mat_differentiation), sampleinfo_differentiation$sample)

#ibdBetas_cluster<-ibd_combo[complete.cases(ibd_combo),]
d <- dist(t(mat_differentiation))
hc <- hclust(d, method = "complete") #single, complete, average, ward

#pdf(here('rna_seq/figs','cluster_wholeBothArrays_diagnosis.pdf'), width=40, h=6)
myplclust(hc, labels=sampleinfo_differentiation$sample, lab.col=as.fumeric(sampleinfo_differentiation$passage_hilo), cex=1.5)
#dev.off()




#### treatment
sampleinfo_treatment<-sampleinfo[which(sampleinfo$comparison=="cytokine"),]
mat_treatment<-mat[,which(colnames(mat)%in%sampleinfo_treatment$sample)]

identical(colnames(mat_treatment), sampleinfo_treatment$sample)

d <- dist(t(mat_treatment))
hc <- hclust(d, method = "complete") #single, complete, average, ward

#pdf(here('rna_seq/figs','cluster_wholeBothArrays_diagnosis.pdf'), width=40, h=6)
myplclust(hc, labels=sampleinfo_treatment$sample, lab.col=as.fumeric(sampleinfo_treatment$passage_hilo), cex=1.5)
#dev.off()
myplclust(hc, labels=sampleinfo_treatment$sample, lab.col=as.fumeric(sampleinfo_treatment$treatment), cex=1.5)



#'## R Session Info
sessionInfo()
