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
datapath = "data/validation_dataset/kallisto"
resultdir = here('data/validation_dataset/kallisto/sleuth')

sampleinfo <- read.table("data/validation_dataset/sample_info_RNA_seq.txt", header=F, sep=" ", skip=1)
sampleinfo<-sampleinfo[,c(2,7,8:10,12)]
colnames(sampleinfo)<-c("sample","concentration","volume","quantity","ratio","well")

#create a sample to condition metadata description
sample_id = sampleinfo$sample
kal_dirs <- file.path(datapath, sample_id)

sampleinfo$path<-kal_dirs

# sampleinfo$inflammation<-as.factor(sampleinfo$inflammation)
# sampleinfo$sex<-as.factor(sampleinfo$sex)
# sampleinfo$diagnosis<-as.factor(sampleinfo$diagnosis)


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


## ignore outlier
mat_noOutlier<-mat[,which(colnames(mat)!="ERR2271027")]
sampleinfo_noOutlier<-sampleinfo[which(sampleinfo$sample!="ERR2271027"),]

save(mat_noOutlier,sampleinfo_noOutlier, file=here("rna_seq/data","ibd_RNAseq_sleuth.RData"))

pca_res_noOutlier <- prcomp(t(mat_noOutlier))
Loadings_meta_noOutlier<-as.data.frame(pca_res_noOutlier$x)
vars <- pca_res_noOutlier$sdev^2
Importance<-vars/sum(vars)


#Restructure meta
Loadings_meta_noOutlier[order(Loadings_meta_noOutlier$PC1),1:5]
sampleinfo_noOutlier$case.no<-as.factor(sampleinfo_noOutlier$case.no)

meta_categorical <- sampleinfo_noOutlier[, c(1,5,6,11,12)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(sampleinfo_noOutlier[, c(13)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Case No.", "Diagnosis", "Sample Site","Inflammation","Sex")
colnames(meta_continuous) <- c("Age")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings_meta_noOutlier, Importance, 3.3, 1.5))

ggsave(here("rna_seq/figs","heat_scree_RNAseq_noOutlier.pdf"), suppressWarnings(heat_scree_plot(Loadings_meta_noOutlier, Importance, 3.3, 1.5)),width = 9, height = 6)
ggsave(here("rna_seq/figs/jpeg","heat_scree_RNAseq_noOutlier.jpeg"), suppressWarnings(heat_scree_plot(Loadings_meta_noOutlier, Importance, 3.3, 1.5)),width = 9, height = 6)

#'# PC vs PC plot
Loadings_meta_noOutlier$sample<-rownames(Loadings_meta_noOutlier)
Loadings_meta_noOutlier<-merge(Loadings_meta_noOutlier, sampleinfo_noOutlier, by="sample")

ggplot(Loadings_meta_noOutlier, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=3, color="black")+theme_bw()+fillscale_sampsite+th+xlab("PC1 (38%)")+ylab("PC2 (12%)")
ggsave(here("rna_seq/figs","PC1_PC2_RNAseq.pdf"), width = 7.5, height = 6)
ggsave(here("rna_seq/figs/jpeg","PC1_PC2_RNAseq.jpeg"), width = 7.5, height = 6)

ggplot(Loadings_meta_noOutlier, aes(PC1, PC2, fill=diagnosis, color=sample.site))+geom_point(shape=21,size=3)+theme_bw()+
  fillscale_diagnosis+th+scale_color_manual(values=c("grey40","black","white"), name="Sample Site")
ggsave(here("rna_seq/figs","PC1_PC2_RNAseq_diagnosis.pdf"), width = 7.5, height = 6)
ggsave(here("rna_seq/figs/jpeg","PC1_PC2_RNAseq_diagnosis.jpeg"), width = 7.5, height = 6)

  
  #APOA1 (transcripts from emsembl)
  APOA1<-as.data.frame((mat_noOutlier[grep("ENSG00000118137", rownames(mat_noOutlier)),]))
  APOA1$sample<-rownames(APOA1)
  colnames(APOA1)[1]<-"ENSG00000118137"
  
  Loadings_meta_noOutlier_apoa1<-merge(Loadings_meta_noOutlier, APOA1, by="sample")
  Loadings_meta_noOutlier_apoa1$ENSG00000118137_noinf_log2<-log2(Loadings_meta_noOutlier_apoa1$ENSG00000118137)
  Loadings_meta_noOutlier_apoa1$ENSG00000118137_noinf_log2[which(is.infinite(Loadings_meta_noOutlier_apoa1$ENSG00000118137_noinf_log2))]<-min(Loadings_meta_noOutlier_apoa1$ENSG00000118137_noinf_log2[which(is.finite(Loadings_meta_noOutlier_apoa1$ENSG00000118137_noinf_log2))])
  
  ggplot(Loadings_meta_noOutlier_apoa1[,c(81,94)], aes(ENSG00000118137_noinf_log2))+geom_histogram()
  
  ## apoa1 pca
  ggplot(Loadings_meta_noOutlier_apoa1, aes(PC1, PC2, fill=ENSG00000118137_noinf_log2))+
    geom_point(shape=21,size=3)+theme_bw()+scale_fill_gradientn(colours = c("#a6d96a","cornflowerblue", "#195dd6"), values=c(0,0.75,1),
                                                                name="APOA1 transcript \n(ENSG00000118137) \nEstimated Counts (log2)")+th
    
  #EpCAM and PTRC (transcripts from emsembl)
  EPCAM<-as.data.frame((mat_noOutlier[grep("ENSG00000119888", rownames(mat_noOutlier)),]))
  EPCAM$sample<-rownames(EPCAM)
  colnames(EPCAM)[1]<-"ENSG00000119888"
  
  Loadings_meta_noOutlier_EPCAM<-merge(Loadings_meta_noOutlier, EPCAM, by="sample")
  
  Loadings_meta_noOutlier_EPCAM$ENSG00000119888_noinf_log2<-log2(Loadings_meta_noOutlier_EPCAM$ENSG00000119888)
  Loadings_meta_noOutlier_EPCAM$ENSG00000119888_noinf_log2[which(is.infinite(Loadings_meta_noOutlier_EPCAM$ENSG00000119888_noinf_log2))]<-min(Loadings_meta_noOutlier_EPCAM$ENSG00000119888_noinf_log2[which(is.finite(Loadings_meta_noOutlier_EPCAM$ENSG00000119888_noinf_log2))])
  
  ggplot(Loadings_meta_noOutlier_EPCAM, aes(PC1, PC2, fill=ENSG00000119888_noinf_log2))+
    geom_point(shape=21,size=2)+theme_bw()+scale_fill_gradient(low="white", high="cornflowerblue")
  
  PTPRC<-as.data.frame((mat_noOutlier[grep("ENSG00000081237", rownames(mat_noOutlier)),]))
  PTPRC$sample<-rownames(PTPRC)
  colnames(PTPRC)[1]<-"ENSG00000081237"
  
  Loadings_meta_noOutlier_PTPRC<-merge(Loadings_meta_noOutlier, PTPRC, by="sample")

  ggplot(Loadings_meta_noOutlier_PTPRC, aes(PC1, PC2, fill=ENSG00000081237))+
    geom_point(shape=21,size=2)+theme_bw()+scale_fill_gradient(low="white", high="blue")

  
 
  
  
#############  
#'# clustering
#############
d <- dist(t(mat))
hc <- hclust(d, method = "complete") #single, complete, average, ward
identical(colnames(mat),sampleinfo$sample)

sampleinfo$label<-sapply(1:nrow(sampleinfo), function(x) paste(sampleinfo$sample, "_", sampleinfo$sample.site, sep=""))
sampleinfo$diagnosis<-as.character(sampleinfo$diagnosis)

sampleinfo$col_sample_site<-as.factor(sampleinfo$sample.site)
levels(sampleinfo$col_sample_site)<-c( "#1a9850","#a6d96a" ,"cornflowerblue" )
sampleinfo$col_sample_site<-as.character(sampleinfo$col_sample_site)

pdf(here('rna_seq/figs','cluster_RNAseq_site.pdf'), width=40, h=25)
myplclust(hc, labels=sampleinfo$label, lab.col=sampleinfo$col_sample_site, cex=1.5)
dev.off()

pdf(here('rna_seq/figs','cluster_RNAseq_diagnosis.pdf'), width=40, h=25)
myplclust(hc, labels=sampleinfo$label, lab.col=as.fumeric(sampleinfo$diagnosis), cex=1.5)
dev.off()





#### #### #### 
#### cluster each tissues seperate
#### #### #### 
              # mat_noOutlier<-mat[,which(colnames(mat)!="ERR2271027")]
              # sampleinfo_noOutlier<-sampleinfo[which(sampleinfo$sample!="ERR2271027"),]
              # 
              # # mat<-mat_noOutlier
              # # sampleinfo<-sampleinfo_noOutlier

## diagnosis and coloring
sampleinfo_noOutlier$diagnosis<-as.factor(sampleinfo_noOutlier$diagnosis)
sampleinfo_noOutlier$diagnosis<-factor(sampleinfo_noOutlier$diagnosis, levels=c("Control", "CD", "UC"))
sampleinfo_noOutlier$col_diagnosis<-as.factor(sampleinfo_noOutlier$diagnosis)
levels(sampleinfo_noOutlier$col_diagnosis)<-c("lightgrey","dodgerblue3","darkgoldenrod1")
sampleinfo_noOutlier$col_diagnosis<-as.character(sampleinfo_noOutlier$col_diagnosis)



#### TI
identical(colnames(mat_noOutlier), sampleinfo_noOutlier$sample)
sampleinfo_TI<-sampleinfo_noOutlier[which(sampleinfo_noOutlier$sample.site=="TI"),]
TI_est_count<-mat_noOutlier[,which(colnames(mat_noOutlier)%in%sampleinfo_TI$sample)]

identical(colnames(TI_est_count), sampleinfo_TI$sample)

#ibdBetas_cluster<-ibd_combo[complete.cases(ibd_combo),]
d <- dist(t(TI_est_count))
hc <- hclust(d, method = "complete") #single, complete, average, ward

k <- cutree(hc, h = 3e+05)
sampleinfo_TI$k<-k

#pdf(here('rna_seq/figs','cluster_wholeBothArrays_diagnosis.pdf'), width=40, h=6)
myplclust(hc, labels=sampleinfo_TI$case.no, lab.col=sampleinfo_TI$col_diagnosis, cex=1.5)
#dev.off()

#pdf("Rplots.pdf", paper = "USr", height = 8.5, width = 11)
     

dend <- as.dendrogram(hc)
rectanle<-ggplot()+geom_rect(aes(xmin=1:nrow(sampleinfo_TI), xmax=1:nrow(sampleinfo_TI)+1, ymin=0, ymax=1, 
                                 fill=sampleinfo_TI$diagnosis[match(labels(dend),sampleinfo_TI$sample)]), color="black", alpha=0.5) +
  theme_bw()+fillscale_diagnosis+theme(legend.position = "none",
                                       axis.title=element_blank(),
                                       axis.text=element_blank(),
                                       axis.ticks=element_blank(),
                                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line = element_blank(),
                                       panel.border = element_blank())+
  geom_vline(xintercept=2, size=1.5)+geom_vline(xintercept=8, size=1.5)+geom_vline(xintercept=18, size=1.5)+ylim(-0.5,1.5)


pdf(here('rna_seq/figs','cluster_TI_diagnosis_RNAseq.pdf'), width=8, h=6, onefile=F)

par(mfcol = c(1,1), mar=c(5,3,5,3), oma=c(0,0,0,0))
plot.new()
myplclust(hc, labels=sampleinfo_TI$case.no, lab.col=sampleinfo_TI$col_diagnosis, cex=1.5)
vp <- viewport(height = unit(0.1,"npc"), width=unit(0.91, "npc"), 
               just = c("center","top"), y = 0.14, x = 0.5)
print(rectanle, vp = vp)
dev.off()


par(mfcol = c(1,1), mar=c(5,3,5,3), oma=c(0,0,0,0))
plot.new()
myplclust(hc, labels=sampleinfo_TI$case.no, lab.col=sampleinfo_TI$col_diagnosis, cex=1.5)
vp <- viewport(height = unit(0.1,"npc"), width=unit(0.91, "npc"), 
               just = c("center","top"), y = 0.14, x = 0.5)
print(rectanle, vp = vp)
dev.off()





#### SC
identical(colnames(mat_noOutlier), sampleinfo_noOutlier$sample)
sampleinfo_SC<-sampleinfo_noOutlier[which(sampleinfo_noOutlier$sample.site=="SC"),]
SC_est_count<-mat_noOutlier[,which(colnames(mat_noOutlier)%in%sampleinfo_SC$sample)]

identical(colnames(SC_est_count), sampleinfo_SC$sample)

#ibdBetas_cluster<-ibd_combo[complete.cases(ibd_combo),]
d <- dist(t(SC_est_count))
hc <- hclust(d, method = "complete") #single, complete, average, ward

k <- cutree(hc, h = 25e+04)
sampleinfo_SC$k<-k

#pdf(here('rna_seq/figs','cluster_wholeBothArrays_diagnosis.pdf'), width=40, h=6)
myplclust(hc, labels=sampleinfo_SC$case.no, lab.col=sampleinfo_SC$col_diagnosis, cex=1.5)
#dev.off()

#pdf("Rplots.pdf", paper = "USr", height = 8.5, width = 11)
dend <- as.dendrogram(hc)
rectanle<-ggplot()+geom_rect(aes(xmin=1:nrow(sampleinfo_SC), xmax=1:nrow(sampleinfo_SC)+1, ymin=0, ymax=1, 
                                 fill=sampleinfo_SC$diagnosis[match(labels(dend),sampleinfo_SC$sample)]), color="black", alpha=0.5) +
  theme_bw()+fillscale_diagnosis+theme(legend.position = "none",
                                       axis.title=element_blank(),
                                       axis.text=element_blank(),
                                       axis.ticks=element_blank(),
                                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line = element_blank(),
                                       panel.border = element_blank())+
  geom_vline(xintercept=2, size=1.5)+geom_vline(xintercept=25, size=1.5)+ylim(-0.5,1.5)

pdf(here('rna_seq/figs','cluster_SC_diagnosis_RNAseq.pdf'), width=8, h=6, onefile=F)

par(mfcol = c(1,1), mar=c(5,3,5,3), oma=c(0,0,0,0))
plot.new()
myplclust(hc, labels=sampleinfo_SC$case.no, lab.col=sampleinfo_SC$col_diagnosis, cex=1.5)
vp <- viewport(height = unit(0.1,"npc"), width=unit(0.91, "npc"), 
               just = c("center","top"), y = 0.14, x = 0.5)
print(rectanle, vp = vp)
dev.off()

par(mfcol = c(1,1), mar=c(5,3,5,3), oma=c(0,0,0,0))
plot.new()
myplclust(hc, labels=sampleinfo_SC$case.no, lab.col=sampleinfo_SC$col_diagnosis, cex=1.5)
vp <- viewport(height = unit(0.1,"npc"), width=unit(0.91, "npc"), 
               just = c("center","top"), y = 0.14, x = 0.5)
print(rectanle, vp = vp)
dev.off()






#'## R Session Info
sessionInfo()
