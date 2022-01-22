#'---
#'title: Validation in GSE144213
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ### Load Libraries
suppressMessages(library(minfi))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(RColorBrewer))
suppressMessages(library(here))
suppressMessages(library(binom))
suppressMessages(library(limma))
suppressMessages(library(sva))
suppressMessages(library(pamr))
suppressMessages(library(GEOquery))
suppressMessages(library(GEOmetadb))

suppressMessages(library(dplyr))
suppressMessages(library(lmtest))
suppressMessages(library(gridExtra))
suppressMessages(library(gtools))
suppressMessages(library(rafalib))


options(stringsAsFactors = FALSE)


#' ### Load Functions
source(here("scripts","00_pretty_plots.R"))
suppressMessages(source(here("scripts","00_heat_scree_plot_generic.R")))
source(here("scripts","00_EM_array_uniform_background_maximise_betabinom.R"))

#' Analysis of https://www.tandfonline.com/doi/full/10.1080/15592294.2020.1762398

gse <- getGEO("GSE144213", GSEMatrix = TRUE)

GSE144213_meta<-pData(gse[[1]]) 

#cd data/published_organoids/GSE144213
#wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144213/suppl/GSE144213_RAW.tar'
#tar -xvf GSE144213_RAW.tar
path<-"data/published_organoids/GSE144213"

GSE144213_meta$Assay.Name<-sapply(1:nrow(GSE144213_meta), function(x){
  gsub("_Grn","",strsplit(GSE144213_meta$supplementary_file[x],"[/]|[.]")[[1]][13])})

GSE144213_meta$array.id.path <- file.path(here(path), GSE144213_meta$Assay.Name)
rgset_organoid <- read.metharray(GSE144213_meta$array.id.path, verbose = FALSE)


#' ### Normalize DNAm Arrays

# Background and dye bias correction with noob thhrough funnorm implemented in minfi
#http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illumina <- preprocessFunnorm(rgset_organoid, sex=GSE144213_meta$`gender:ch1`)
GSE144213_beta<-getBeta(MSet.illumina)


#' ## Organize meta data
GSE144213_meta$sentrix_ID<-sapply(1:nrow(GSE144213_meta), function(x){
  strsplit(GSE144213_meta$Assay.Name[x],"_")[[1]][2]})
GSE144213_meta<-GSE144213_meta[,c(2,40:47,48,50)] 
colnames(GSE144213_meta)<-gsub(":ch1","",colnames(GSE144213_meta))



#' ### Detection p values across all probes for each sample
avg_detPval <- colMeans(detectionP(rgset_organoid))
GSE144213_meta$det_pval<-avg_detPval

ggplot(GSE144213_meta)+geom_boxplot(aes(as.factor(sentrix_ID), det_pval, fill=as.factor(sentrix_ID)), outlier.shape = NA)+
  geom_point(aes(as.factor(sentrix_ID), det_pval, group=Assay.Name, fill=as.factor(sentrix_ID)), shape=21, color="black",
             position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE)+ylim(0,0.008)

ggsave(here("figs","GSE144213_detection_pvalue.pdf"), width=6, height=5)
ggsave(here("figs/jpeg","GSE144213_detection_pvalue.jpeg"), width=6, height=5)




#' Beta distribution before and after normalization
beta_raw<-getBeta(rgset_organoid)
betas<-getBeta(MSet.illumina)

Beta_melted<- melt(betas)
Beta_melted_raw<- melt(beta_raw)

#remove NAs before plotting (otherwise get many non-inifnite warnings)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
Beta_Plot_raw<-Beta_melted_raw[which(!(is.na(Beta_melted_raw$value))),]

#add meta
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,GSE144213_meta, by.x="ID", by.y="Assay.Name")
colnames(Beta_Plot_raw)<-c("CpG","ID","Beta")
Beta_Plot_raw<-merge(Beta_Plot_raw,GSE144213_meta, by.x="ID", by.y="Assay.Name")

dis1<-ggplot(Beta_Plot, aes(Beta, group=as.character(ID), color=as.character(tissue)))+
  geom_density()+theme_bw()+xlab("DNAm Beta Value")

dis2<-ggplot(Beta_Plot_raw, aes(Beta, group=as.character(ID), color=as.character(tissue)))+
  geom_density()+theme_bw()+xlab("DNAm Beta Value")

ggsave(here("figs","GSE144213_beta_distribution.pdf"), grid.arrange(dis1, dis2),w=8, h=5)
ggsave(here("figs/jpeg","GSE144213_beta_distribution.jpeg"), grid.arrange(dis1, dis2), w=8, h=5)







#' ### Probe Filtering 
GSE144213_beta <- GSE144213_beta[!grepl("rs",rownames(GSE144213_beta)), ]
print(paste("Samples available: ",ncol(GSE144213_beta),"\nProbes available: ",nrow(GSE144213_beta),sep=""))

#' #### 450K annotation from illumina
# https://emea.support.illumina.com/downloads/humanmethylation450_15017482_v1-2_product_files.html
anno_450k<-read.csv(here("data","HumanMethylation450_15017482_v1-2.csv"), skip=7)
anno_450k<-anno_450k[match(rownames(GSE144213_beta),anno_450k$IlmnID),]
identical(rownames(GSE144213_beta), anno_450k$IlmnID)

#' #### Sex Chromosomes 
GSE144213_beta<-GSE144213_beta[which(!(anno_450k$CHR%in%c('X','Y'))),] 
filt_sex<-nrow(GSE144213_beta)
print(paste("Samples available: ",ncol(GSE144213_beta),"Probes available: ",nrow(GSE144213_beta),sep=""))


#' #### Cross-hybridizing probes and polymorphic probes. 
#' Some probes have been found to cross-hybridize with other chromosomes (Price et al. 2013 *Epigenetics*).
#' https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16304
price<-read.table(here("data","GPL16304-47833.txt"), sep='\t', header=T, skip=22)
price<-price[match(rownames(GSE144213_beta),price$ID),]

cross_hyb<-price[which(price$XY_Hits=="XY_YES" | price$Autosomal_Hits=="A_YES"),]
GSE144213_beta<-GSE144213_beta[which(!(rownames(GSE144213_beta)%in%cross_hyb$ID)),]
filt_cross<-nrow(GSE144213_beta)
print(paste("Samples available: ",ncol(GSE144213_beta),"\nProbes available: ",nrow(GSE144213_beta),sep=""))


#' Polymorphic probes
SnpatCpG<-price[which(price$Target.CpG.SNP!=""),]
GSE144213_beta<-GSE144213_beta[which(!(rownames(GSE144213_beta)%in%SnpatCpG$ID)),]
filt_poly<-nrow(GSE144213_beta)
print(paste("Samples available: ",ncol(GSE144213_beta),"\nProbes available: ",nrow(GSE144213_beta),sep=""))

#' #### Probe filtering based on detection pvalue and detection over background (NA)
#' Remove probes with high NA count
na_count_probe <-sapply(1:nrow(GSE144213_beta), function(y) length(which(is.na(GSE144213_beta[y,]))))
na_count_probe_good<-which(na_count_probe<(ncol(GSE144213_beta)*0.05))
GSE144213_beta<-GSE144213_beta[na_count_probe_good,]
filt_bead<-nrow(GSE144213_beta)
print(paste("Samples available: ",ncol(GSE144213_beta),"\nProbes available: ",nrow(GSE144213_beta),sep=""))

#' Remove probes with high detection p value across samples, and any samples with many high detection p value probes
detP <- detectionP(rgset_organoid)
failed <- detP>0.05
bad_det_p<-names(which(rowMeans(failed)>0.01))
bad_det_psamp<-names(which(colMeans(failed)>0.01))

GSE144213_beta<-GSE144213_beta[which(!(rownames(GSE144213_beta)%in%bad_det_p)),]
GSE144213_beta<-GSE144213_beta[,which(!(colnames(GSE144213_beta)%in%bad_det_psamp))]

filt_detp<-nrow(GSE144213_beta)
print(paste("Samples available: ",ncol(GSE144213_beta),"\nProbes available: ",nrow(GSE144213_beta),sep=""))
colnames(GSE144213_meta)<-gsub(" ", "_",colnames(GSE144213_meta))


#' ### PCA
pca_res <- prcomp(t(GSE144213_beta))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

GSE144213_meta$sentrix_ID<-as.factor(GSE144213_meta$sentrix_ID)
GSE144213_meta$age_at_acquisition<-as.factor(GSE144213_meta$age_at_acquisition)

meta_categorical <- GSE144213_meta[, c(4,5,6,8,9,11)]  # input column numbers in meta that contain categorical variables
meta_continuous <- GSE144213_meta[, c(2,12)]  # input column numbers in meta that contain continuous variables

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
PCs_to_view<-10
suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7))

## PC vs PC plot
Loadings$Assay.Name<-rownames(Loadings)
Loadings_meta<-merge(Loadings, GSE144213_meta, by="Assay.Name")

ggplot(Loadings_meta, aes(PC1, PC2, fill=primary_site))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab(paste("PC1 (",round(Importance[1]*100,0),"%)", sep=""))+ylab(paste("PC2 (",round(Importance[2]*100,0),"%)", sep=""))+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=c("#cb181d","#1d91c0","#df65b0","#fc8d59","#c7e9b4","#41ab5d"))



save(GSE144213_beta, GSE144213_meta, file=paste(here("data"),"/GSE144213_beta_organoids.RData",sep=""))
#load(here("data", "GSE144213_beta_organoids.RData"))


#' #### Variable Beta Distribution (not fetal)
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
Mval<-function(beta) log2(beta/(1-beta))

GSE144213_mval= apply(GSE144213_beta, 1, Mval) # need mvalues for combat
GSE144213_mval = as.data.frame(GSE144213_mval)
GSE144213_mval = t(GSE144213_mval)

ref_range_dnam<-sapply(1:nrow(GSE144213_mval), function(x) Variation(GSE144213_mval[x,]))
GSE144213_beta_VeryVariable<-GSE144213_beta[which(ref_range_dnam>=2.75),]

print(paste("There are ",nrow(GSE144213_beta_VeryVariable), " variable CpGs (10th-90th quantile range in M value >2.75)",sep=""))

dim(GSE144213_beta_VeryVariable<-GSE144213_beta[rev(order(ref_range_dnam)),])
#' Include the same number of vairable CpGs as organoids varible. So take the 71384 most variable
dim(GSE144213_beta_VeryVariable<-GSE144213_beta_VeryVariable[1:71384 ,])

## Beta distribution plot
Beta_melted<- melt(GSE144213_beta_VeryVariable)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,GSE144213_meta, by.x="ID", by.y="Assay.Name")

ggplot(Beta_Plot, aes(Beta,color=primary_site))+
  geom_density(size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=c("#cb181d","#1d91c0","#df65b0","#fc8d59","#c7e9b4","#41ab5d"))+
  th+theme(legend.text = element_text(size=7),
           legend.title = element_text(size=10),
           legend.key.size = unit(0.7,"line"))


ggplot()+
  geom_density(data=Beta_Plot, aes(Beta,color=primary_site, group=ID),size=0.75)+
  theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  th+theme(legend.text = element_text(size=7),
           legend.title = element_text(size=10),
           legend.key.size = unit(0.7,"line"))+
  scale_color_manual(values=c("#cb181d","#1d91c0","#df65b0","#fc8d59","#c7e9b4","#41ab5d"), guide=F)+
  facet_wrap(~primary_site)
ggsave(here("figs","GSE144213_variable_CpGs_not_split_tri.pdf"),width = 7, height = 5)
ggsave(here("figs/jpeg","GSE144213_variable_CpGs_not_split_tri.jpeg"), width = 7, height = 5)






#' ## Thresholding trimodal samples

GSE144213_meta$thresholded_prior_ratio<-sapply(1:nrow(GSE144213_meta), function(x){
  converted<-as.numeric(round(GSE144213_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  passage_threshold_params(converted, counts, res)
})

GSE144213_meta$thresholded_ratio_max<-"Bimodal"
GSE144213_meta$thresholded_ratio_max[which(GSE144213_meta$thresholded_prior_ratio>1)]<-"Trimodal"



## complicated color palette
GSE144213_meta$color<-as.factor(GSE144213_meta$primary_site)
levels(GSE144213_meta$color)<-c("#cb181d","#1d91c0","#5aae61","#fccde5","#41b6c4","#454514")

color_df<-do.call(rbind, lapply(unique(GSE144213_meta$color), function(x) {
  colgroup<-GSE144213_meta$color[which(GSE144213_meta$color==x)]
  data.frame(ID=GSE144213_meta$Assay.Name[which(GSE144213_meta$color==x)], color=sapply(1:length(colgroup), function(y){ muted(colgroup[y], l=y*10)}))}))

GSE144213_meta_color<-color_df$color
names(GSE144213_meta_color)<-color_df$ID

GSE144213_meta_color[5:6]<-c("#fd8d3c","#d94801")
GSE144213_meta_color[15]<-c("#762a83")
GSE144213_meta_color[23:25]<-c("#3690c0","#67a9cf","#a6bddb")

labels<-as.data.frame(table(GSE144213_meta$primary_site, GSE144213_meta$thresholded_ratio_max))
colnames(labels)<-c("primary_site","thresholded_ratio_max","count")


Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,GSE144213_meta, by.x="ID", by.y="Assay.Name")



ggplot()+
  geom_density(data=Beta_Plot, aes(Beta,color=ID, group=ID),size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  th+theme(legend.text = element_text(size=7),
           legend.title = element_text(size=10),
           legend.key.size = unit(0.7,"line"))+
  scale_color_manual(values=GSE144213_meta_color, guide=F)+
  facet_grid(primary_site~thresholded_ratio_max)+
  geom_text(aes(x= 0.5, y=4,label=paste("n=",count, sep="")),data=labels)



ggsave(here("figs","GSE144213_variable_CpGs.pdf"),width = 7.5, height = 15)
ggsave(here("figs/jpeg","GSE144213_variable_CpGs.jpeg"), width = 7.5, height = 15)



## plot all samples
plts_paired<-lapply(1:nrow(GSE144213_meta), function(x){
  passage<-paste("Primary site:",GSE144213_meta$primary_site[x],"\nIndividual: ", GSE144213_meta$patient_id[x],"\nRatio I/H: " ,round(GSE144213_meta$thresholded_prior_ratio[x],2), sep="")
  converted<-as.numeric(round(GSE144213_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  draw_fit_params_gg(converted, counts, res,passage)
})


pdf(here("figs","GSE144213_organoids_thresholding_all_samples.pdf"))
plts_paired
dev.off()





#'## R Session Info
sessionInfo()
