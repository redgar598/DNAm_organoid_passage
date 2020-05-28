#'---
#'title: Validation in E-MTAB-4957
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

options(stringsAsFactors = FALSE)


#' ### Load Functions
source(here("scripts","00_pretty_plots.R"))
suppressMessages(source(here("scripts","00_heat_scree_plot_generic.R")))
source(here("scripts","00_EM_array_uniform_background_maximise_betabinom.R"))




#' ### Download Data
# wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4957/E-MTAB-4957.raw.1.zip
# wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4957/E-MTAB-4957.raw.2.zip
# wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4957/E-MTAB-4957.raw.3.zip
# wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4957/E-MTAB-4957.raw.4.zip
# wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4957/E-MTAB-4957.raw.5.zip
# unzip E-MTAB-4957.raw.1.zip 
# unzip E-MTAB-4957.raw.2.zip 
# unzip E-MTAB-4957.raw.3.zip 
# unzip E-MTAB-4957.raw.4.zip 
# unzip E-MTAB-4957.raw.5.zip 

#wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4957/E-MTAB-4957.sdrf.txt
path<-"data/published_organoids/E_MTAB_4957"

sampleinfo_organoid <- read.table(here(path, "E-MTAB-4957.sdrf.txt"), header=T, sep="\t")
sampleinfo_organoid<-sampleinfo_organoid[,c(1:15,30:34,38:41)]

#sample info cleanup
sampleinfo_organoid$sentrix<-sapply(1:nrow(sampleinfo_organoid), function(x) strsplit(as.character(sampleinfo_organoid$Assay.Name)[x], "_")[[1]][1])

sampleinfo_organoid$array.id.path <- file.path(here(path), sampleinfo_organoid$Assay.Name)
rgset_organoid <- read.metharray(sampleinfo_organoid$array.id.path, verbose = FALSE)


#' ### Normalize DNAm Arrays

# Background and dye bias correction with noob thhrough funnorm implemented in minfi
#http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illumina <- preprocessFunnorm(rgset_organoid)
MTAB_organoid_beta<-getBeta(MSet.illumina)



#' ### Detection p values across all probes for each sample
avg_detPval <- colMeans(detectionP(rgset_organoid))
sampleinfo_organoid$det_pval<-avg_detPval

print("detection pval")

ggplot(sampleinfo_organoid)+geom_boxplot(aes(as.factor(sentrix), det_pval, fill=as.factor(sentrix)), outlier.shape = NA)+
  geom_point(aes(as.factor(sentrix), det_pval, group=Assay.Name, fill=as.factor(sentrix)), shape=21, color="black",
             position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE)+ylim(0,0.008)

ggsave(here("figs","MTAB4957_detection_pvalue.pdf"), width=6, height=5)
ggsave(here("figs/jpeg","MTAB4957_detection_pvalue.jpeg"), width=6, height=5)




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
Beta_Plot<-merge(Beta_Plot,sampleinfo_organoid, by.x="ID", by.y="Assay.Name")
colnames(Beta_Plot_raw)<-c("CpG","ID","Beta")
Beta_Plot_raw<-merge(Beta_Plot_raw,sampleinfo_organoid, by.x="ID", by.y="Assay.Name")

dis1<-ggplot(Beta_Plot, aes(Beta, group=as.character(ID), color=as.character(Characteristics.sampling.site.)))+
  geom_density()+theme_bw()+xlab("DNAm Beta Value")

dis2<-ggplot(Beta_Plot_raw, aes(Beta, group=as.character(ID), color=as.character(Characteristics.sampling.site.)))+
  geom_density()+theme_bw()+xlab("DNAm Beta Value")

ggsave(here("figs","MTAB4957_beta_distribution.pdf"), grid.arrange(dis1, dis2),w=20, h=5)
ggsave(here("figs/jpeg","MTAB4957_beta_distribution.jpeg"), grid.arrange(dis1, dis2), w=20, h=5)







#' ### Probe Filtering 
MTAB_organoid_beta <- MTAB_organoid_beta[!grepl("rs",rownames(MTAB_organoid_beta)), ]
print(paste("Samples available: ",ncol(MTAB_organoid_beta),"\nProbes available: ",nrow(MTAB_organoid_beta),sep=""))

#' #### 450K annotation from illumina
# https://emea.support.illumina.com/downloads/humanmethylation450_15017482_v1-2_product_files.html
anno_450k<-read.csv(here("data","HumanMethylation450_15017482_v1-2.csv"), skip=7)
anno_450k<-anno_450k[match(rownames(MTAB_organoid_beta),anno_450k$IlmnID),]
identical(rownames(MTAB_organoid_beta), anno_450k$IlmnID)

#' #### Sex Chromosomes 
MTAB_organoid_beta<-MTAB_organoid_beta[which(!(anno_450k$CHR%in%c('X','Y'))),] 
filt_sex<-nrow(MTAB_organoid_beta)
print(paste("Samples available: ",ncol(MTAB_organoid_beta),"Probes available: ",nrow(MTAB_organoid_beta),sep=""))


#' #### Cross-hybridizing probes and polymorphic probes. 
#' Some probes have been found to cross-hybridize with other chromosomes (Price et al. 2013 *Epigenetics*).
#' https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16304
price<-read.table(here("data","GPL16304-47833.txt"), sep='\t', header=T, skip=22)
price<-price[match(rownames(MTAB_organoid_beta),price$ID),]

cross_hyb<-price[which(price$XY_Hits=="XY_YES" | price$Autosomal_Hits=="A_YES"),]
MTAB_organoid_beta<-MTAB_organoid_beta[which(!(rownames(MTAB_organoid_beta)%in%cross_hyb$ID)),]
filt_cross<-nrow(MTAB_organoid_beta)
print(paste("Samples available: ",ncol(MTAB_organoid_beta),"\nProbes available: ",nrow(MTAB_organoid_beta),sep=""))


#' Polymorphic probes
SnpatCpG<-price[which(price$Target.CpG.SNP!=""),]
MTAB_organoid_beta<-MTAB_organoid_beta[which(!(rownames(MTAB_organoid_beta)%in%SnpatCpG$ID)),]
filt_poly<-nrow(MTAB_organoid_beta)
print(paste("Samples available: ",ncol(MTAB_organoid_beta),"\nProbes available: ",nrow(MTAB_organoid_beta),sep=""))

#' #### Probe filtering based on detection pvalue and detection over background (NA)
#' Remove probes with high NA count
na_count_probe <-sapply(1:nrow(MTAB_organoid_beta), function(y) length(which(is.na(MTAB_organoid_beta[y,]))))
na_count_probe_good<-which(na_count_probe<(ncol(MTAB_organoid_beta)*0.05))
MTAB_organoid_beta<-MTAB_organoid_beta[na_count_probe_good,]
filt_bead<-nrow(MTAB_organoid_beta)
print(paste("Samples available: ",ncol(MTAB_organoid_beta),"\nProbes available: ",nrow(MTAB_organoid_beta),sep=""))

#' Remove probes with high detection p value across samples, and any samples with many high detection p value probes
detP <- detectionP(rgset_organoid)
failed <- detP>0.05
bad_det_p<-names(which(rowMeans(failed)>0.01))
bad_det_psamp<-names(which(colMeans(failed)>0.01))

MTAB_organoid_beta<-MTAB_organoid_beta[which(!(rownames(MTAB_organoid_beta)%in%bad_det_p)),]
MTAB_organoid_beta<-MTAB_organoid_beta[,which(!(colnames(MTAB_organoid_beta)%in%bad_det_psamp))]

filt_detp<-nrow(MTAB_organoid_beta)
print(paste("Samples available: ",ncol(MTAB_organoid_beta),"\nProbes available: ",nrow(MTAB_organoid_beta),sep=""))







#' ### Tidy meta data available
sampleinfo_organoid<-sampleinfo_organoid[,c(1:16,24:27)]
identical(colnames(MTAB_organoid_beta),sampleinfo_organoid$Assay.Name)
head(sampleinfo_organoid)


#' #### Seperate out primary pediatric samples
sampleinfo_ped_primary<-sampleinfo_organoid[which(sampleinfo_organoid$Characteristics.developmental.stage.=="juvenile stage" & sampleinfo_organoid$Characteristics.biosource.type.=="purified cell"),]
print(paste("Primary pediatric samples available: ",nrow(sampleinfo_ped_primary), sep=""))

# Three samples are the non-epithelial portion
sampleinfo_ped_primary<-sampleinfo_ped_primary[-grep("non-epithelial",sampleinfo_ped_primary$Characteristics.cell.type.),]
organoid_ped_primary<-MTAB_organoid_beta[,which(colnames(MTAB_organoid_beta)%in%sampleinfo_ped_primary$Assay.Name)]
identical(colnames(organoid_ped_primary),sampleinfo_ped_primary$Assay.Name)
print(paste("Primary samples available: ",nrow(sampleinfo_ped_primary), sep=""))


#' #### Seperate out primary fetal samples
sampleinfo_fetal_primary<-sampleinfo_organoid[which(sampleinfo_organoid$Characteristics.developmental.stage.=="fetal stage" & sampleinfo_organoid$Characteristics.biosource.type.=="purified cell"),]
print(paste("Primary fetal samples available: ",nrow(sampleinfo_fetal_primary), sep=""))

organoid_fetal_primary<-MTAB_organoid_beta[,which(colnames(MTAB_organoid_beta)%in%sampleinfo_fetal_primary$Assay.Name)]
identical(colnames(organoid_fetal_primary),sampleinfo_fetal_primary$Assay.Name)
print(paste("Primary samples available: ",nrow(sampleinfo_fetal_primary), sep=""))




#' #### Only organoids
sampleinfo_organoid<-sampleinfo_organoid[which(sampleinfo_organoid$Characteristics.biosource.type.=="organoid"),]
MTAB_organoid_beta<-MTAB_organoid_beta[,which(colnames(MTAB_organoid_beta)%in%sampleinfo_organoid$Assay.Name)]
identical(colnames(MTAB_organoid_beta),sampleinfo_organoid$Assay.Name)

print(paste("Organoid samples available: ",nrow(sampleinfo_organoid), sep=""))


#' #### Numeric passage
sampleinfo_organoid$passage.or.rescope.no_numeric<-as.factor(as.character(sampleinfo_organoid$Characteristics.passage.))
levels(sampleinfo_organoid$passage.or.rescope.no_numeric)<-c(1,11,12,14,18,2,21,23,3,4,5,6,7,8,9)
sampleinfo_organoid$passage.or.rescope.no_numeric<-as.numeric(as.character(sampleinfo_organoid$passage.or.rescope.no_numeric))

#' sample_ID to include patient and tissue
sampleinfo_organoid$sample_ID<-paste(sampleinfo_organoid$Characteristics.individual., sampleinfo_organoid$Characteristics.sampling.site.)


#' #### Extracting information from source name and associated paper

#' remove AP as it seems to indicate a mixture of passages (i.e P1AP6)? These aen't mentioned in manscript?
sampleinfo_organoid<-sampleinfo_organoid[-grep("AP",sampleinfo_organoid$Source.Name),]
#' Unclear what GM and IM are but it seems like only 281 was under these conditions, so will exclude to be safe
sampleinfo_organoid<-sampleinfo_organoid[-grep("281",sampleinfo_organoid$Source.Name),]

#' tidy condition
sampleinfo_organoid$condition<-NA
#' need to grep for WT and Clone B
sampleinfo_organoid$condition[grep("WT",sampleinfo_organoid$Source.Name)]<-"WT"
sampleinfo_organoid$condition[grep("Clone B|CloneB",sampleinfo_organoid$Source.Name)]<-"KO"
#' DM seems to be differenitated and SCM is maybe undifferentiated
sampleinfo_organoid$condition[grep("SCM",sampleinfo_organoid$Source.Name)]<-"SCM"
sampleinfo_organoid$condition[grep("DM",sampleinfo_organoid$Source.Name)]<-"DM4d"

print(paste("Organoid samples available: ",nrow(sampleinfo_organoid), sep=""))

# match the DNAm data
MTAB_organoid_beta<-MTAB_organoid_beta[,which(colnames(MTAB_organoid_beta)%in%sampleinfo_organoid$Assay.Name)]
identical(colnames(MTAB_organoid_beta),sampleinfo_organoid$Assay.Name)






#' ##### Filter samples later run in both dataset 
#' Individuals "212" "223" "224" "225" "369" are in both studies. 
#' The organoids from these individuals will be removed from the differential DNAm analysis, but visulized later to show the patterns occur across array type

#' Saving them for later integration
sampleinfo_organoid_paired<-sampleinfo_organoid[grep("212|223|224",sampleinfo_organoid$Characteristics.individual.),]
MTAB4957_beta_for_paird<-MTAB_organoid_beta[,which(colnames(MTAB_organoid_beta)%in%sampleinfo_organoid_paired$Assay.Name)]
identical(colnames(MTAB4957_beta_for_paird), sampleinfo_organoid_paired$Assay.Name)
print(paste("Organoid samples available: ",nrow(sampleinfo_organoid_paired), sep=""))


#' load the EPIC organoids
load(here("data","beta_organoids.RData"))

intersect(sampleinfo_organoid$Characteristics.individual., epic.organoid$case.no)

#' 369  and 225 are run at identical passages in both cohorts so all samples will be removed. 
sampleinfo_organoid[grep("369",sampleinfo_organoid$Characteristics.individual.),c(1,4,8,9,13,21,22)]
epic.organoid[grep("369",epic.organoid$case.no),c(1,5,17, 9, 10)]
sampleinfo_organoid<-sampleinfo_organoid[-grep("369|225",sampleinfo_organoid$Characteristics.individual.),]

#' 212,223 and 224 have organoids at passages only used in MTAB-4957. So organoids that are included at the same passage will be removed. 
remove_sampleID<-c("212 SC P6","212 TI P6","223 SC P2","223 TI P2","224 SC P3","224 TI P2")
sampleinfo_organoid<-sampleinfo_organoid[which(!(sampleinfo_organoid$Source.Name%in%remove_sampleID)),]

print(paste("Organoid samples available: ",nrow(sampleinfo_organoid), sep=""))

# match the DNAm data
MTAB_organoid_beta<-MTAB_organoid_beta[,which(colnames(MTAB_organoid_beta)%in%sampleinfo_organoid$Assay.Name)]
identical(colnames(MTAB_organoid_beta),sampleinfo_organoid$Assay.Name)


#' ### Fetal, pedatric and adult organoids together
pca_res <- prcomp(t(MTAB_organoid_beta))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

meta_categorical <- sampleinfo_organoid[, c(4,7,8,13,17,18)]  # input column numbers in meta that contain categorical variables
meta_continuous <- sampleinfo_organoid[, c(9,21)]  # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Individual", "Developmental Stage","Sample Site","Sex","Block","Sentrix ID")
colnames(meta_continuous) <- c("Age", "Passage")
meta_continuous$Age<-as.numeric(meta_continuous$Age)
meta_categorical$Block<-as.factor(meta_categorical$Block)


ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
PCs_to_view<-10
suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7))

## PC vs PC plot
Loadings$Assay.Name<-rownames(Loadings)
Loadings_meta<-merge(Loadings, sampleinfo_organoid, by="Assay.Name")

ggplot(Loadings_meta, aes(PC1, PC2, fill=Characteristics.developmental.stage.))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (29%)")+ylab("PC2 (10%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))

ggplot(Loadings_meta, aes(PC1, PC2, fill=Characteristics.sampling.site.))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (29%)")+ylab("PC2 (10%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))

ggplot(Loadings_meta, aes(PC2, PC3, fill=Characteristics.sampling.site.))+geom_point(shape=21,size=3, color="black")+theme_bw()+
    xlab("PC2 (10%)")+ylab("PC3 (8%)")+th+theme(axis.text = element_text(size=12),
                                                 axis.title = element_text(size=14),
                                                 plot.margin = margin(1, 0.1, 1, 1, "cm"))



pc_plt<-ggplot(Loadings_meta, aes(PC2, PC3, fill=as.factor(passage.or.rescope.no_numeric),color=Characteristics.developmental.stage.))+geom_line(aes(PC2,PC3, group=sample_ID), color="lightgrey")+#, color=sampling.time.point
  geom_point(shape=21,size=3)+#
  theme_bw()+xlab("PC2 (10%)")+ylab("PC3 (8%)")+th+theme(axis.text = element_text(size=12),axis.title = element_text(size=14))+
  scale_fill_manual(values=pass_col,name="Passage\nNumber")+scale_color_manual(values=c("black","white","black"))

legend<-ggplot(sampleinfo_organoid, aes(as.factor(-passage.or.rescope.no_numeric), fill=as.factor(passage.or.rescope.no_numeric)))+geom_bar(color="black")+
  theme_bw()+theme(legend.position = "none", axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   legend.title=element_text(size=10),
                   legend.text=element_text(size=8))+
  coord_flip()+
  scale_fill_manual(values=pass_col,name="Passage\nNumber")+th

r <- ggplot() + theme_void()

grid.arrange(pc_plt,arrangeGrob(r,legend,r, heights=c(0.6,1.25,0.4)), ncol=2, widths=c(7,1))


#' High passage fetal toward older organoids?
ggplot(Loadings_meta, aes(PC1, PC2, fill=as.factor(passage.or.rescope.no_numeric), color=Characteristics.developmental.stage.))+geom_line(aes(PC1,PC2, group=sample_ID), color="lightgrey")+#, color=sampling.time.point
  geom_point(shape=21,size=3)+#
  theme_bw()+xlab("PC2 (10%)")+ylab("PC3 (8%)")+th+theme(axis.text = element_text(size=12),axis.title = element_text(size=14))+
  scale_fill_manual(values=pass_col,name="Passage\nNumber")+scale_color_manual(values=c("black","white","black"))









#' # Non Fetal samples
sampleinfo_organoid_notfetal<-sampleinfo_organoid[which(sampleinfo_organoid$Characteristics.biosource.type.=="organoid" & sampleinfo_organoid$Characteristics.developmental.stage.!="fetal stage"),]
MTAB_organoid_beta_notfetal<-MTAB_organoid_beta[,which(colnames(MTAB_organoid_beta)%in%sampleinfo_organoid_notfetal$Assay.Name)]
identical(colnames(MTAB_organoid_beta_notfetal),sampleinfo_organoid_notfetal$Assay.Name)


#' ### PCA 
pca_res <- prcomp(t(MTAB_organoid_beta_notfetal))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

meta_categorical <- sampleinfo_organoid_notfetal[, c(4,7,8,13,17,18)]  # input column numbers in meta that contain categorical variables
meta_continuous <- sampleinfo_organoid_notfetal[, c(9,21)]  # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Individual", "Developmental Stage","Sample Site","Sex","Block","Sentrix ID")
colnames(meta_continuous) <- c("Age", "Passage")
meta_continuous$Age<-as.numeric(meta_continuous$Age)
meta_categorical$Block<-as.factor(meta_categorical$Block)

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
PCs_to_view<-10
suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7))


## PC vs PC plot
Loadings$Assay.Name<-rownames(Loadings)
Loadings_meta<-merge(Loadings, sampleinfo_organoid_notfetal, by="Assay.Name")

#' Sample Site
ggplot(Loadings_meta, aes(PC1, PC2, fill=Characteristics.sampling.site.))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (19%)")+ylab("PC2 (13%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))
#' Condition
ggplot(Loadings_meta, aes(PC1, PC2, fill=condition))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (19%)")+ylab("PC2 (13%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))
#' Sample Site PC2/3
ggplot(Loadings_meta, aes(PC2, PC3, fill=Characteristics.sampling.site.))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC2 (13%)")+ylab("PC3 (10%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))
#' Sample Site PC3/4
ggplot(Loadings_meta, aes(PC3, PC4, fill=Characteristics.sampling.site.))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC3 (10%)")+ylab("PC4 (6%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))



pc_plt<-ggplot(Loadings_meta, aes(PC2, PC3, fill=as.factor(passage.or.rescope.no_numeric)))+geom_line(aes(PC2,PC3, group=sample_ID), color="lightgrey")+#, color=sampling.time.point
  geom_point(shape=21,size=3)+#
  theme_bw()+xlab("PC2 (10%)")+ylab("PC3 (8%)")+th+theme(axis.text = element_text(size=12),axis.title = element_text(size=14))+
  scale_fill_manual(values=pass_col,name="Passage\nNumber")+scale_color_manual(values=c("black","white","black"))

legend<-ggplot(sampleinfo_organoid_notfetal, aes(as.factor(-passage.or.rescope.no_numeric), fill=as.factor(passage.or.rescope.no_numeric)))+geom_bar(color="black")+
  theme_bw()+theme(legend.position = "none", axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   legend.title=element_text(size=10),
                   legend.text=element_text(size=8))+
  coord_flip()+
  scale_fill_manual(values=pass_col,name="Passage\nNumber")+th

r <- ggplot() + theme_void()

grid.arrange(pc_plt,arrangeGrob(r,legend,r, heights=c(0.6,1.25,0.4)), ncol=2, widths=c(7,1))



#' #### Variable Beta Distribution (not fetal)
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
Mval<-function(beta) log2(beta/(1-beta))

MTAB4957_mval= apply(MTAB_organoid_beta_notfetal, 1, Mval) # need mvalues for combat
MTAB4957_mval = as.data.frame(MTAB4957_mval)
MTAB4957_mval = t(MTAB4957_mval)

ref_range_dnam<-sapply(1:nrow(MTAB4957_mval), function(x) Variation(MTAB4957_mval[x,]))
dim(MTAB4957_beta_VeryVariable<-MTAB_organoid_beta_notfetal[which(ref_range_dnam>=2.75),])#  47924

## Beta distribution plot
Beta_melted<- melt(MTAB4957_beta_VeryVariable)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,sampleinfo_organoid_notfetal, by.x="ID", by.y="Assay.Name")

Beta_Plot$passage.or.rescope.no_numeric.factor <- factor(Beta_Plot$passage.or.rescope.no_numeric, levels = c(11,9,8,6,4,3,2,1))

ggplot(Beta_Plot, aes(Beta,color=passage.or.rescope.no_numeric.factor))+
  geom_density(size=1)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=pass_col, name="Passage\nNumber")+th+theme(legend.text = element_text(size=7),
                                                      legend.title = element_text(size=10),
                                                      legend.key.size = unit(0.7,"line"))



#' To view the beta distributions we will also include a line for the 32 primary samples to compare each passage to primary 
MTAB4957_mval= apply(organoid_ped_primary, 1, Mval) # need mvalues for combat
MTAB4957_mval = as.data.frame(MTAB4957_mval)
MTAB4957_mval = t(MTAB4957_mval)

ref_range_dnam_primary<-sapply(1:nrow(MTAB4957_mval), function(x) Variation(MTAB4957_mval[x,]))
dim(organoid_ped_primary_VeryVariable<-organoid_ped_primary[rev(order(ref_range_dnam_primary)),])
#' Include the same number of vairable CpGs as organoids varible. So take the 47924 most variable
dim(organoid_ped_primary_VeryVariable<-organoid_ped_primary_VeryVariable[1:47924 ,])

Beta_melted_MTAB_primary<- melt(organoid_ped_primary_VeryVariable)
Beta_Plot_MTAB_primary<-Beta_melted_MTAB_primary[which(!(is.na(Beta_melted_MTAB_primary$value))),]
colnames(Beta_Plot_MTAB_primary)<-c("CpG","ID","Beta")
Beta_Plot_MTAB_primary<-merge(Beta_Plot_MTAB_primary,sampleinfo_ped_primary, by.x="ID", by.y="Assay.Name")
Beta_plot_primary<-Beta_Plot_MTAB_primary[,c(1:3)]
Beta_plot_primary$passage.or.rescope.no_numeric<-0
Beta_Plot<-Beta_Plot[,c(1:3,23)]

Beta_Plot_combined<-rbind(Beta_plot_primary,Beta_Plot)

Beta_Plot_combined$passage.or.rescope.no_numeric.factor <- factor(Beta_Plot_combined$passage.or.rescope.no_numeric, levels = c(11,9,8,6,4,3,2,1,0))

ggplot(Beta_Plot_combined, aes(Beta,color=as.factor(passage.or.rescope.no_numeric.factor)))+
  geom_density(size=1)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=pass_col, name="Passage\nNumber")+th+theme(legend.text = element_text(size=7),
                                                                       legend.title = element_text(size=10),
                                                                       legend.key.size = unit(0.7,"line"))


ggsave(here("figs","MTAB4957_beta_notfetal_with_primary.pdf"),width = 3.75, height = 2.5)
ggsave(here("figs/jpeg","MTAB4957_beta_notfetal_with_primary.jpeg"), width = 3.75, height = 2.5)





#' #### Paired Samples in beta plots
#' These are just the from the individuals not shared across studies
sampleinfo_organoid_notfetal$sample.site<-as.factor(sampleinfo_organoid_notfetal$Characteristics.sampling.site.)
levels(sampleinfo_organoid_notfetal$sample.site)<-c("SC","TI")
sampleinfo_organoid_notfetal$sample_ID<-paste(sampleinfo_organoid_notfetal$Characteristics.individual., sampleinfo_organoid_notfetal$sample.site)
sampleinfo_organoid_paired_unique<-sampleinfo_organoid_notfetal[which(sampleinfo_organoid_notfetal$sample_ID%in%sampleinfo_organoid_notfetal$sample_ID[duplicated(sampleinfo_organoid_notfetal$sample_ID)]),]


MTAB4957.organoid_paired<-do.call(rbind,lapply(1:length(unique(sampleinfo_organoid_paired_unique$sample_ID)), function(x){
  sample<-unique(sampleinfo_organoid_paired_unique$sample_ID)[x]
  samp<-sampleinfo_organoid_paired_unique[sampleinfo_organoid_paired_unique$sample_ID==sample,]
  samp<-samp[order(samp$passage.or.rescope.no_numeric),]
  samp$hilo<-as.factor(samp$passage.or.rescope.no_numeric)

  if(length(levels(samp$hilo))==2){levels(samp$hilo)<-c("lower","higher")}else{
    if(length(levels(samp$hilo))==3){levels(samp$hilo)<-c("lower","higher","highest")}else{
      if(length(levels(samp$hilo))==4){levels(samp$hilo)<-c("lowest","lower","higher","highest")}else{samp$hilo<-NA}
    }
  }
  samp
}))

MTAB4957.organoid_paired<-MTAB4957.organoid_paired[which(!is.na(MTAB4957.organoid_paired$hilo)),]
## will include 224 in the combined plot with its passage 2,3 samples that were on the epic
MTAB4957.organoid_paired<-MTAB4957.organoid_paired[-(grep("224", MTAB4957.organoid_paired$Source.Name)),]

MTAB4957_beta_VeryVariable_paird<-MTAB4957_beta_VeryVariable[,which(sampleinfo_organoid_notfetal$sample_ID%in%MTAB4957.organoid_paired$sample_ID)]

Beta_melted<- melt(MTAB4957_beta_VeryVariable_paird)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,MTAB4957.organoid_paired, by.x="ID", by.y="Assay.Name")


labels<-as.data.frame(tapply(MTAB4957.organoid_paired$passage.or.rescope.no_numeric, MTAB4957.organoid_paired$sample_ID, function(x) paste(x, collapse=", ")))
colnames(labels)<-"passge"
labels$sample_ID<-rownames(labels)

ggplot()+
  geom_density(aes(Beta,color=hilo, group=ID),Beta_Plot, size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values = c ("#9ecae1", "#225ea8", "#081d58"), name="Relative\nPassage\nLevel within\nPatient")+facet_wrap(~sample_ID, nrow=1)+
  geom_text(aes(0.75, 2.75, label=passge), data=labels, color="grey20")+th+theme(strip.text = element_text(size = 10),
                                                                                 axis.text=element_text(size=4),
                                                                                 panel.spacing = unit(0.7, "lines"))+th+
  scale_x_continuous(breaks = c(0,0.5,1))

ggsave(here("figs","MTAB4957_beta_paired_notfetal.pdf"),width = 5, height = 2.2)
ggsave(here("figs/jpeg","MTAB4957_beta_paired_notfetal.jpeg"), width = 5, height = 2.2)




#' #### Paired not fetal EPIC and 450K
#' these are individuals present in both studies
intersect(sampleinfo_organoid_paired$Characteristics.individual., epic.organoid$case.no)

organoid_Mval = apply(organoid_beta, 1, Mval) # need mvalues for combat
organoid_Mval = as.data.frame(organoid_Mval)
organoid_Mval = t(organoid_Mval)

ref_range_dnam<-sapply(1:nrow(organoid_Mval), function(x) Variation(organoid_Mval[x,]))
organoid_beta_VeryVariable<-organoid_beta[which(ref_range_dnam>=2.75),]#  51545


#' "212" "223" "224"  in both studes
epic.organoid_paired<-epic.organoid[grep("212|223|224",epic.organoid$case.no),]

MTAB4957_beta_VeryVariable_paird<-MTAB4957_beta_for_paird[which(rownames(MTAB4957_beta_for_paird)%in%rownames(MTAB4957_beta_VeryVariable)),]
identical(colnames(MTAB4957_beta_VeryVariable_paird), as.character(sampleinfo_organoid_paired$Assay.Name))


epic.organoid_beta_VeryVariable_paird<-organoid_beta_VeryVariable[,which(colnames(organoid_beta_VeryVariable)%in%epic.organoid_paired$array.id)]
epic.organoid_beta_VeryVariable_paird<-epic.organoid_beta_VeryVariable_paird[,match(epic.organoid_paired$array.id,colnames(epic.organoid_beta_VeryVariable_paird))]
identical(colnames(epic.organoid_beta_VeryVariable_paird), as.character(epic.organoid_paired$array.id))

epic.organoid_paired$Study<-"Original Organoids"

sampleinfo_organoid_paired$Study<-"E-MTAB-4957"
sampleinfo_organoid_paired$sample_ID<-gsub(" sigmoid colon", " SC",sampleinfo_organoid_paired$sample_ID)
sampleinfo_organoid_paired$sample_ID<-gsub(" terminal ileum", " TI",sampleinfo_organoid_paired$sample_ID)
colnames(sampleinfo_organoid_paired)[c(4,16)]<-c("case.no","array.id")

samplinfo_paired_combined<-rbind(epic.organoid_paired[,c(1,2,14,17,18)],sampleinfo_organoid_paired[,c(4,16,22,21,24)])

variable_beta_combined<-merge(epic.organoid_beta_VeryVariable_paird,MTAB4957_beta_VeryVariable_paird, by="row.names")
variable_beta_combined$Row.names<-NULL
identical(colnames(variable_beta_combined), as.character(samplinfo_paired_combined$array.id))


samplinfo_paired_combined<-do.call(rbind,lapply(1:length(unique(samplinfo_paired_combined$sample_ID)), function(x){
  sample<-unique(samplinfo_paired_combined$sample_ID)[x]
  samp<-samplinfo_paired_combined[samplinfo_paired_combined$sample_ID==sample,]
  samp<-samp[order(samp$passage.or.rescope.no_numeric),]
  samp$hilo<-as.factor(samp$passage.or.rescope.no_numeric)
  
  if(length(levels(samp$hilo))==2){levels(samp$hilo)<-c("lower","higher")}else{
    if(length(levels(samp$hilo))==3){levels(samp$hilo)<-c("lower","higher","highest")}else{
      if(length(levels(samp$hilo))==4){levels(samp$hilo)<-c("lowest","lower","higher","highest")}else{samp$hilo<-NA}
    }
  }
  samp
}))


Beta_melted<- melt(variable_beta_combined)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("ID","Beta")
Beta_Plot<-merge(Beta_Plot,samplinfo_paired_combined, by.x="ID", by.y="array.id")

Beta_Plot$hilo_epic<-paste(Beta_Plot$hilo,"\n", Beta_Plot$Study, sep="")
Beta_Plot$hilo_epic<-factor(Beta_Plot$hilo_epic, levels=c("lower\nOriginal Organoids", "lower\nE-MTAB-4957", "higher\nE-MTAB-4957", "highest\nE-MTAB-4957"))

labels<-as.data.frame(tapply(samplinfo_paired_combined$passage.or.rescope.no_numeric, list(samplinfo_paired_combined$sample_ID,samplinfo_paired_combined$Study), function(x) paste(x, collapse=", ")))
labels$sample_ID<-rownames(labels)
colnames(labels)<-c("MTAB","OG","sample_ID")

labels$MTAB<-paste("E-MTAB-4957: ", labels$MTAB, sep="")
labels$OG<-paste("Original Organoid: ", labels$OG, sep="")



ggplot()+
  geom_density(aes(Beta,color=hilo_epic, group=ID),Beta_Plot, size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values = c ("#ef3b2c", "#9ecae1","#225ea8","#081d58"), name="Relative\nPassage\nLevel within\nPatient")+facet_wrap(~sample_ID, nrow=3)+
  geom_text(aes(0.5, 2.65, label=MTAB), data=labels, color="grey20", size=2.75)+
  geom_text(aes(0.5, 2.35, label=OG), data=labels, color="grey20", size=2.75)+
  th+theme(strip.text = element_text(size = 12),
           axis.text=element_text(size=10),
           panel.spacing = unit(0.7, "lines"),
        legend.text=element_text(size=9.5),
        legend.title=element_text(size=12))+ scale_x_continuous(breaks = c(0,0.5,1))

ggsave(here("figs","MTAB4957_EPIC_beta_paired.pdf"),width = 6, height = 6)
ggsave(here("figs/jpeg","MTAB4957_EPIC_beta_paired.jpeg"), width = 6, height = 6)





#' ## Thresholding trimodal samples

sampleinfo_organoid_notfetal$thresholded_prior_ratio<-sapply(1:nrow(sampleinfo_organoid_notfetal), function(x){
  print(x)
  converted<-as.numeric(round(MTAB4957_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  passage_threshold_params(converted, counts, res)
})

ggplot(sampleinfo_organoid_notfetal, aes(as.numeric(as.character(passage.or.rescope.no_numeric)), thresholded_prior_ratio))+
  geom_point(size=2,shape=21,color="black",aes(fill=as.factor(passage.or.rescope.no_numeric)))+xlab("Passage")+
  ylab("Intermediate Peak Prior")+theme_bw()+theme(axis.title = element_text(size=12))+
  #geom_text(aes(label=count, vjust=vjust, hjust=hjust), color="grey40", size=3)+
  scale_x_continuous(breaks=c(1,2,3,4,6,7,8,2,4,10,11,14,16))+ scale_fill_manual(values=pass_col,name="Passage\nNumber", guide=F)

ggsave(here("figs","MTAB4957_mixture_model_ratio_maximize.pdf"), width=3, height=2)


sampleinfo_organoid_notfetal$thresholded_ratio_max<-F
sampleinfo_organoid_notfetal$thresholded_ratio_max[which(sampleinfo_organoid_notfetal$thresholded_prior_ratio>1)]<-T

percent_passing<-round((tapply(sampleinfo_organoid_notfetal$thresholded_ratio_max, sampleinfo_organoid_notfetal$passage.or.rescope.no_numeric, sum)/tapply(sampleinfo_organoid_notfetal$array.id, sampleinfo_organoid_notfetal$passage.or.rescope.no_numeric, length))*100,2)
passed_num<-tapply(sampleinfo_organoid_notfetal$thresholded_ratio_max, sampleinfo_organoid_notfetal$passage.or.rescope.no_numeric, sum)
org_numer<-tapply(sampleinfo_organoid_notfetal$array.id, sampleinfo_organoid_notfetal$passage.or.rescope.no_numeric, length)

df<-data.frame(passage=names(percent_passing), passing=percent_passing, pro_passing=percent_passing/100, count=org_numer, passed_num=passed_num)
df$passage.factor <- factor(df$passage, levels = c(11,9,8,6,4,3,2,1))


df<-cbind(df,(binom.confint(df$passed_num, df$count, method="exact", conf.level=0.95)))
df$upper<-df$upper*100
df$lower<-df$lower*100

ggplot(df, aes(as.numeric(as.character(passage)), passing))+
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="grey70", width=.25)+
  geom_line(color="grey20")+geom_point(size=1.25,shape=21,color="black",aes(fill=passage.factor))+xlab("Passage")+
  ylab("Samples with Trimodal\nDistribution (%)")+theme_bw()+theme(axis.title = element_text(size=10))+
  scale_x_continuous(breaks=c(1,2,3,4,6,7,8,2,4,10,11,14,16))+ scale_fill_manual(values=pass_col,name="Passage\nNumber", guide=F)

ggsave(here("figs","MTAB4957_mixture_model_ratio_threshold_maximize.pdf"), width=3, height=2)



## plot all samples
plts_paired<-lapply(1:nrow(sampleinfo_organoid_notfetal), function(x){
  print(x)
  passage<-paste("passage: ",sampleinfo_organoid_notfetal$passage.or.rescope.no_numeric[x],"\nIndividual: ", sampleinfo_organoid_notfetal$Characteristics.individual.[x],"\nRatio I/H: " ,round(epic.organoid$thresholded_prior_ratio[x],2), sep="")
  converted<-as.numeric(round(MTAB4957_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  draw_fit_params_gg(converted, counts, res,passage)
})

plts_paired_order<-plts_paired[order(sampleinfo_organoid_notfetal$passage.or.rescope.no_numeric)]

pdf(here("figs","MTAB4957_organoids_thresholding_all_samples.pdf"))
plts_paired_order
dev.off()








#' # Fetal Organoids

sampleinfo_organoid_fetal<-sampleinfo_organoid[which(sampleinfo_organoid$Characteristics.biosource.type.=="organoid" & sampleinfo_organoid$Characteristics.developmental.stage.=="fetal stage"),]
sampleinfo_organoid_fetal<-sampleinfo_organoid_fetal[-grep("KO", sampleinfo_organoid_fetal$condition),]
MTAB_organoid_beta_fetal<-MTAB_organoid_beta[,which(colnames(MTAB_organoid_beta)%in%sampleinfo_organoid_fetal$Assay.Name)]
identical(colnames(MTAB_organoid_beta_fetal),sampleinfo_organoid_fetal$Assay.Name)


# ' ### PCA organoids
pca_res <- prcomp(t(MTAB_organoid_beta_fetal))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

meta_categorical <- sampleinfo_organoid_fetal[, c(4,8,13,17,18)]  # input column numbers in meta that contain categorical variables
meta_continuous <- sampleinfo_organoid_fetal[, c(9,21)]  # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Individual", "Sample Site","Sex","Block","Sentrix ID")
colnames(meta_continuous) <- c("Age", "Passage")
meta_continuous$Age<-as.numeric(meta_continuous$Age)
meta_categorical$Block<-as.factor(meta_categorical$Block)

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
PCs_to_view<-10
suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7))


## PC vs PC plot
Loadings$Assay.Name<-rownames(Loadings)
Loadings_meta<-merge(Loadings, sampleinfo_organoid_fetal, by="Assay.Name")

#' Sample Site
ggplot(Loadings_meta, aes(PC1, PC2, fill=Characteristics.sampling.site.))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (20%)")+ylab("PC2 (13%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))




pc_plt<-ggplot(Loadings_meta, aes(PC1, PC2, fill=as.factor(passage.or.rescope.no_numeric)))+geom_line(aes(PC1,PC2, group=sample_ID), color="lightgrey")+#, color=sampling.time.point
  geom_point(shape=21,size=3)+#
  theme_bw()+xlab("PC1 (36%)")+ylab("PC1 (14%)")+th+theme(axis.text = element_text(size=12),axis.title = element_text(size=14))+
  scale_fill_manual(values=c(colorRampPalette(brewer.pal(11, "Spectral"))(11), "#544791","#4a3e80", "#40366f","#221d3c"),name="Passage\nNumber")+scale_color_manual(values=c("black","white","black"))

legend<-ggplot(sampleinfo_organoid_fetal, aes(as.factor(-passage.or.rescope.no_numeric), fill=as.factor(passage.or.rescope.no_numeric)))+geom_bar(color="black")+
  theme_bw()+theme(legend.position = "none", axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   legend.title=element_text(size=10),
                   legend.text=element_text(size=8))+
  coord_flip()+
  scale_fill_manual(values=c(colorRampPalette(brewer.pal(11, "Spectral"))(11), "#544791","#4a3e80", "#40366f","#221d3c"),name="Passage\nNumber")+th

r <- ggplot() + theme_void()

grid.arrange(pc_plt,arrangeGrob(r,legend,r, heights=c(0.6,1.25,0.4)), ncol=2, widths=c(7,1))




#' #### Variable Beta Distribution (fetal)
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
Mval<-function(beta) log2(beta/(1-beta))

MTAB4957_mval= apply(MTAB_organoid_beta_fetal, 1, Mval) # need mvalues for combat
MTAB4957_mval = as.data.frame(MTAB4957_mval)
MTAB4957_mval = t(MTAB4957_mval)


ref_range_dnam<-sapply(1:nrow(MTAB4957_mval), function(x) Variation(MTAB4957_mval[x,]))
dim(MTAB4957_beta_VeryVariable<-MTAB_organoid_beta_fetal[which(ref_range_dnam>=2.75),])#  28418

## Beta distribution plot
Beta_melted<- melt(MTAB4957_beta_VeryVariable)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,sampleinfo_organoid_fetal, by.x="ID", by.y="Assay.Name")

Beta_Plot$passage.or.rescope.no_numeric.factor <- factor(Beta_Plot$passage.or.rescope.no_numeric, levels = c(23,21,14,12,9,7,6,5,3,2,1))

ggplot(Beta_Plot, aes(Beta,color=passage.or.rescope.no_numeric.factor))+
  geom_density(size=1)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=pass_col, name="Passage\nNumber")+th+theme(legend.text = element_text(size=7),
                                                                       legend.title = element_text(size=10),
                                                                       legend.key.size = unit(0.7,"line"))

#' To view the beta distributions we will also include a line for the 11 primary fetal samples to compare each passage to primary 
identical(colnames(organoid_fetal_primary),sampleinfo_fetal_primary$Assay.Name)

MTAB4957_mval= apply(organoid_fetal_primary, 1, Mval) # need mvalues for combat
MTAB4957_mval = as.data.frame(MTAB4957_mval)
MTAB4957_mval = t(MTAB4957_mval)

ref_range_dnam_fetal<-sapply(1:nrow(MTAB4957_mval), function(x) Variation(MTAB4957_mval[x,]))
dim(organoid_fetal_primary_VeryVariable<-organoid_fetal_primary[rev(order(ref_range_dnam_fetal)),])
dim(organoid_fetal_primary_VeryVariable<-organoid_fetal_primary_VeryVariable[1:28418 ,])# same number as MTAB organoid varible


Beta_melted_MTAB_primary<- melt(organoid_fetal_primary_VeryVariable)
Beta_Plot_MTAB_primary<-Beta_melted_MTAB_primary[which(!(is.na(Beta_melted_MTAB_primary$value))),]
colnames(Beta_Plot_MTAB_primary)<-c("CpG","ID","Beta")
Beta_Plot_MTAB_primary<-merge(Beta_Plot_MTAB_primary,sampleinfo_fetal_primary, by.x="ID", by.y="Assay.Name")
Beta_plot_primary<-Beta_Plot_MTAB_primary[,c(1:3)]
Beta_plot_primary$passage.or.rescope.no_numeric<-0
Beta_Plot<-Beta_Plot[,c(1:3,23)]

Beta_Plot_combined<-rbind(Beta_plot_primary,Beta_Plot)

Beta_Plot_combined$passage.or.rescope.no_numeric.factor <- factor(Beta_Plot_combined$passage.or.rescope.no_numeric, levels = c(23,21,14,12,9,7,6,5,3,2,1,0))

ggplot(Beta_Plot_combined, aes(Beta,color=as.factor(passage.or.rescope.no_numeric.factor)))+
  geom_density(size=1)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=pass_col, name="Passage\nNumber")+th+theme(legend.text = element_text(size=7),
                                                                       legend.title = element_text(size=10),
                                                                       legend.key.size = unit(0.7,"line"))


ggsave(here("figs","MTAB4957_beta_fetal_with_primary.pdf"),width = 3.75, height = 2.5)
ggsave(here("figs/jpeg","MTAB4957_beta_fetal_with_primary.jpeg"), w=5, h=3)



#' #### Paired Samples in beta plots
sampleinfo_organoid_fetal$sample_ID<-paste(sampleinfo_organoid_fetal$Characteristics.individual., sampleinfo_organoid_fetal$Characteristics.sampling.site.)
sampleinfo_organoid_paired<-sampleinfo_organoid_fetal[which(sampleinfo_organoid_fetal$sample_ID%in%sampleinfo_organoid_fetal$sample_ID[duplicated(sampleinfo_organoid_fetal$sample_ID)]),]

MTAB4957.organoid_paired<-do.call(rbind,lapply(1:length(unique(sampleinfo_organoid_paired$sample_ID)), function(x){
  sample<-unique(sampleinfo_organoid_paired$sample_ID)[x]
  samp<-sampleinfo_organoid_paired[sampleinfo_organoid_paired$sample_ID==sample,]
  samp<-samp[order(samp$passage.or.rescope.no_numeric),]
  samp$hilo<-as.factor(samp$passage.or.rescope.no_numeric)

  if(length(levels(samp$hilo))==2){levels(samp$hilo)<-c("lower","higher")}else{
    if(length(levels(samp$hilo))==3){levels(samp$hilo)<-c("lower","higher","highest")}else{
      if(length(levels(samp$hilo))==4){levels(samp$hilo)<-c("lowest","lower","higher","highest")}else{samp$hilo<-NA}
    }
  }
  samp
}))

MTAB4957.organoid_paired<-MTAB4957.organoid_paired[which(!is.na(MTAB4957.organoid_paired$hilo)),]
MTAB4957.organoid_paired$hilo<-factor(MTAB4957.organoid_paired$hilo, c("lowest","lower","higher","highest"))

identical(colnames(MTAB4957_beta_VeryVariable), sampleinfo_organoid_fetal$Assay.Name)

MTAB4957_beta_VeryVariable_paird<-MTAB4957_beta_VeryVariable[,which(sampleinfo_organoid_fetal$sample_ID%in%MTAB4957.organoid_paired$sample_ID)]



Beta_melted<- melt(MTAB4957_beta_VeryVariable_paird)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,MTAB4957.organoid_paired, by.x="ID", by.y="Assay.Name")


labels<-as.data.frame(tapply(MTAB4957.organoid_paired$passage.or.rescope.no_numeric, MTAB4957.organoid_paired$sample_ID, function(x) paste(x, collapse=", ")))
colnames(labels)<-"passge"
labels$sample_ID<-rownames(labels)

ggplot()+
  geom_density(aes(Beta,color=hilo, group=ID),Beta_Plot, size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values = c ("#9ecae1","#4292c6", "#225ea8", "#081d58"), name="Relative\nPassage\nLevel within\nPatient")+facet_wrap(~sample_ID, nrow=2)+
  geom_text(aes(0.75, 2.75, label=passge), data=labels, color="grey20")+th+theme(strip.text = element_text(size = 10),
                                                                                 axis.text=element_text(size=4),
                                                                                 panel.spacing = unit(0.7, "lines"))+
  scale_x_continuous(breaks = c(0,0.5,1))

ggsave(here("figs","MTAB4957_beta_paired_fetal.pdf"),w=6, height = 3.75)
ggsave(here("figs/jpeg","MTAB4957_beta_paired_fetal.jpeg"), w=6, height = 3.75)



#' ## Thresholding trimodal samples
sampleinfo_organoid_fetal$thresholded_prior_ratio<-sapply(1:nrow(sampleinfo_organoid_fetal), function(x){
  print(x)
  converted<-as.numeric(round(MTAB4957_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  passage_threshold_params(converted, counts, res)
})

ggplot(sampleinfo_organoid_fetal, aes(as.numeric(as.character(passage.or.rescope.no_numeric)), thresholded_prior_ratio))+
  geom_point(size=2,shape=21,color="black",aes(fill=as.factor(passage.or.rescope.no_numeric)))+xlab("Passage")+
  ylab("Intermediate Peak Prior")+theme_bw()+theme(axis.title = element_text(size=12))+
  #geom_text(aes(label=count, vjust=vjust, hjust=hjust), color="grey40", size=3)+
  scale_x_continuous(breaks=c(1,2,3,4,6,7,8,2,4,10,11,14,16))+ scale_fill_manual(values=pass_col,name="Passage\nNumber", guide=F)

ggsave(here("figs","MTAB4957_fetal_mixture_model_ratio_maximize.pdf"), width=3, height=2)


sampleinfo_organoid_fetal$thresholded_ratio_max<-F
sampleinfo_organoid_fetal$thresholded_ratio_max[which(sampleinfo_organoid_fetal$thresholded_prior_ratio>1)]<-T

percent_passing<-round((tapply(sampleinfo_organoid_fetal$thresholded_ratio_max, sampleinfo_organoid_fetal$passage.or.rescope.no_numeric, sum)/tapply(sampleinfo_organoid_fetal$array.id, sampleinfo_organoid_fetal$passage.or.rescope.no_numeric, length))*100,2)
passed_num<-tapply(sampleinfo_organoid_fetal$thresholded_ratio_max, sampleinfo_organoid_fetal$passage.or.rescope.no_numeric, sum)
org_numer<-tapply(sampleinfo_organoid_fetal$array.id, sampleinfo_organoid_fetal$passage.or.rescope.no_numeric, length)

df<-data.frame(passage=names(percent_passing), passing=percent_passing, pro_passing=percent_passing/100, count=org_numer, passed_num=passed_num)
df$passage.factor <- factor(df$passage, levels = c(23,21,14,12,9,7,6,5,3,2,1))


df<-cbind(df,(binom.confint(df$passed_num, df$count, method="exact", conf.level=0.95)))
df$upper<-df$upper*100
df$lower<-df$lower*100

ggplot(df, aes(as.numeric(as.character(passage)), passing))+
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="grey70", width=.25)+
  geom_line(color="grey20")+geom_point(size=1.25,shape=21,color="black",aes(fill=passage.factor))+xlab("Passage")+
  ylab("Samples with Trimodal\nDistribution (%)")+theme_bw()+theme(axis.title = element_text(size=10))+
  scale_x_continuous(breaks=c(1,2,3,4,6,7,8,2,4,10,11,14,16))+ scale_fill_manual(values=pass_col,name="Passage\nNumber", guide=F)

ggsave(here("figs","MTAB4957_fetal_mixture_model_ratio_threshold_maximize.pdf"), width=3, height=2)


## plot all samples
plts_paired<-lapply(1:nrow(sampleinfo_organoid_fetal), function(x){
  print(x)
  passage<-paste("passage: ",sampleinfo_organoid_fetal$passage.or.rescope.no_numeric[x],"\nIndividual: ", sampleinfo_organoid_fetal$Characteristics.individual.[x],"\nRatio I/H: " ,round(epic.organoid$thresholded_prior_ratio[x],2), sep="")
  converted<-as.numeric(round(MTAB4957_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  draw_fit_params_gg(converted, counts, res,passage)
})

plts_paired_order<-plts_paired[order(sampleinfo_organoid_fetal$passage.or.rescope.no_numeric)]

pdf(here("figs","MTAB4957_fetal_organoids_thresholding_all_samples.pdf"))
plts_paired_order
dev.off()










#' ## Differential methylation with passage not fetal
mod<-model.matrix(~ 0 + passage.or.rescope.no_numeric, data=sampleinfo_organoid_notfetal)
fit <- lmFit(MTAB_organoid_beta_notfetal, mod)
ebfit <- eBayes(fit)

# covariate adjusted beta values
beta<-MTAB_organoid_beta_notfetal

passage_db<-sapply(1:nrow(beta), function(x){
  sampleinfo_cpg<-sampleinfo_organoid_notfetal
  sampleinfo_cpg$beta<-as.numeric(beta[x,])
  
  fit<-lm(beta ~ passage.or.rescope.no_numeric, data=sampleinfo_cpg)
  pval<-summary(fit)$coef["passage.or.rescope.no_numeric","Pr(>|t|)"]
  slope<-fit$coefficients[2]
  
  (min(sampleinfo_organoid_notfetal$passage.or.rescope.no_numeric)*slope) - (max(sampleinfo_organoid_notfetal$passage.or.rescope.no_numeric)*slope)})

passage_MTAB<-data.frame(p.value=ebfit$p.value[,"passage.or.rescope.no_numeric"], CpG=rownames(beta), db=passage_db)

# Adjust P values
passage_MTAB$p_adjusted<-p.adjust(passage_MTAB$p.value, method="BH")

diff_CpG_dbMTAB<-passage_MTAB[which(passage_MTAB$p_adjusted<0.05 & abs(passage_MTAB$db)>0.15),] #25086
diff_CpG_db_hypoMTAB<-diff_CpG_dbMTAB$CpG[which((diff_CpG_dbMTAB$db)>=0.15)] #  17419
diff_CpG_db_hyperMTAB<-diff_CpG_dbMTAB$CpG[which((diff_CpG_dbMTAB$db)<=(-0.15))] #  7667


#' ### load original organoid passage CpGs
load(here("data","beta_organoids.RData"))
pvals_long<-read.csv(here("data","Heteroskedactiy_pvalues_FDR_1000iter.csv"), header=T)
pvals_long[,1]<-NULL
colnames(pvals_long)<-c("BP_count","diff_count","mean_db","fdr_BP","fdr_diff")
pvals_long$CpG<-rownames(organoid_beta)
pvals_long$BP_pval<-((1000-pvals_long$BP_count)+1)/1001
pvals_long$diff_pval<-((1000-pvals_long$diff_count)+1)/1001
pvals_long$BP_fdr<-((1000-pvals_long$fdr_BP)+1)/1001
pvals_long$diff_fdr<-((1000-pvals_long$fdr_diff)+1)/1001

hetero_CpG<-rownames(organoid_beta)[which(pvals_long$BP_fdr<0.05)]
print(paste("CpGs with significant (adjusted p<0.05) heteroskedactiy: ", length(hetero_CpG), sep=""))
diff_CpG<-rownames(organoid_beta)[which(pvals_long$diff_fdr<0.05)]
diff_CpG_db<-pvals_long[which(pvals_long$diff_fdr<0.05 & abs(pvals_long$mean_db)>0.15),]
print(paste("CpGs with significant (adjusted p<0.05; delta beta >0.05) differential methylation: ", nrow(diff_CpG_db), sep=""))

diff_CpG_db_hypo<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)>=0.15)] #  11772
diff_CpG_db_hyper<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)<=(-0.15))] #  5214


#' How many differential CpGs could overlap? 450K vs EPIC
print(paste("Of the ", nrow(diff_CpG_db), " CpGs differential with passage in the original organoids, ",
            length(diff_CpG_db$CpG[which(diff_CpG_db$CpG%in%rownames(MTAB_organoid_beta_notfetal))]), " are in the 450K data", sep=""))

diff_CpG_db_hypo_overlap<-diff_CpG_db_hypo[which(diff_CpG_db_hypo%in%rownames(MTAB_organoid_beta_notfetal))]
diff_CpG_db_hyper_overlap<-diff_CpG_db_hyper[which(diff_CpG_db_hyper%in%rownames(MTAB_organoid_beta_notfetal))]

diff_CpG_db_hypoMTAB_overlap<-diff_CpG_db_hypoMTAB[which(diff_CpG_db_hypoMTAB%in%pvals_long$CpG)]
diff_CpG_db_hyperMTAB_overlap<-diff_CpG_db_hyperMTAB[which(diff_CpG_db_hyperMTAB%in%pvals_long$CpG)]

print(paste("Of the ",length(diff_CpG_db_hypo_overlap)," hypo CpGs also on the 450K ",
            length(intersect(diff_CpG_db_hypoMTAB_overlap, diff_CpG_db_hypo_overlap))," are also hypo in the MTAB-4957 cohort (",
            round((length(intersect(diff_CpG_db_hypoMTAB_overlap, diff_CpG_db_hypo_overlap))/length(diff_CpG_db_hypo_overlap))*100,2),"%)",sep=""))

print(paste("Of the ",length(diff_CpG_db_hyper_overlap)," hyper CpGs also on the 450K ",
            length(intersect(diff_CpG_db_hyperMTAB_overlap, diff_CpG_db_hyper_overlap))," are also hypo in the MTAB-4957 cohort (",
            round((length(intersect(diff_CpG_db_hyperMTAB_overlap, diff_CpG_db_hyper_overlap))/length(diff_CpG_db_hyper_overlap))*100,2),"%)",sep=""))


#' ### delta beta directionality plot
plt_db_direction<-merge(pvals_long[,c(3,6)], passage_MTAB, by="CpG")# 380776

plt_db_direction$sig<-"Not Significant"
plt_db_direction$sig[which(plt_db_direction$CpG%in%c(intersect(diff_CpG_db_hypoMTAB_overlap, diff_CpG_db_hypo_overlap),intersect(diff_CpG_db_hyperMTAB_overlap, diff_CpG_db_hyper_overlap)))]<-"Significant\nIn Both\nCohorts"

ggplot(plt_db_direction, aes(mean_db, db))+geom_point(aes(color=sig, alpha=sig),shape=19)+th+theme_bw()+
  scale_color_manual(values=c("lightgrey", "cornflowerblue"), name="Significant\nWith Passage")+
  scale_alpha_manual(values=c(0.25,1), guide=F)+
  geom_hline(yintercept=c(-0.15,0.15), color="grey60")+geom_vline(xintercept=c(-0.15,0.15), color="grey60")+
  ylim(-0.8,0.8)+xlim(-0.8,0.8)+xlab("Original Organoid\nPassage Delta Beta")+ylab("MTAB-4957 Organoid\nPassage Delta Beta")+
  stat_smooth(method="lm", se=F, color="black")


#ggsave(here("figs","MTAB_db_directionality.pdf"), width=5, height=3.75)
ggsave(here("figs/jpeg","MTAB_db_directionality.jpeg"), width=5, height=3.75)

print(paste("Correlation of delta betas between cohorts: ", round(cor(plt_db_direction$db, plt_db_direction$mean_db),2), sep=""))


#' ### representative CpGs
epic.organoid_minimal<-epic.organoid[,c(2, 14, 17)]
colnames(epic.organoid_minimal)[1]<-"Assay.Name"
epic.organoid_minimal$cohort<-"Original Organoids"

sampleinfo_organoid_notfetal_minimal<-sampleinfo_organoid_notfetal[,c(16,22,21)]
sampleinfo_organoid_notfetal_minimal$cohort<-"MTAB-4957 Organoids"

sample_info_both<-rbind(sampleinfo_organoid_notfetal_minimal,epic.organoid_minimal)


plt_hetero_MTAB<-function(CpGs, legend, axislab, title){
  betas<-melt(cbind(MTAB_organoid_beta_notfetal[CpGs,],organoid_beta[CpGs,]))
  organoid_plt<-merge(sample_info_both, betas, by.x="Assay.Name",by.y="Var2")
  
  p<-ggplot(organoid_plt, aes(passage.or.rescope.no_numeric,value))+
    geom_line(aes(group=sample_ID),color="lightgrey")+
    stat_smooth(method="lm", color="grey30", size=0.7, se=F)+th+theme_bw()+
    geom_point(aes(fill=as.factor(passage.or.rescope.no_numeric)),shape=21, size=1.25)+
    scale_fill_manual(values=pass_col,name="Passage\nNumber", drop=T)+
    facet_grid(cohort~Var1)+
    ylab("DNAm Beta")+xlab("Passage Number")+ylim(0,1)+
    theme(plot.margin = margin(0.5, 0.15, 0.5, 0.15, "cm"),plot.title = element_text(size=12))+
    xlim(1,16)
  
  if(missing(legend) & missing(axislab) & missing(title)){p}else{
    if(legend=="N" & axislab=="N"){p + theme(legend.position = "none",axis.title.y=element_blank(),
                                             axis.text.y=element_blank(),
                                             axis.ticks.y=element_blank())+ ggtitle(title)}else{
                                               if(legend=="N" & axislab=="Y"){p + theme(legend.position = "none") + ggtitle(title)}}}}

plt_hetero_MTAB(c("cg25402228","cg09146328"))

ggsave(here("figs","Passage_differential_CpGs_MTAB4957.pdf"),width = 4.75, height = 4)
ggsave(here("figs/jpeg","Passage_differential_CpGs_MTAB4957.jpeg"), width = 4.75, height = 4)






#' ## Differential methylation with passage fetal
mod<-model.matrix(~ 0 + passage.or.rescope.no_numeric, data=sampleinfo_organoid_fetal)
fit <- lmFit(MTAB_organoid_beta_fetal, mod)
ebfit <- eBayes(fit)

# covariate adjusted beta values
beta<-MTAB_organoid_beta_fetal

passage_db<-sapply(1:nrow(beta), function(x){
  sampleinfo_cpg<-sampleinfo_organoid_fetal
  sampleinfo_cpg$beta<-as.numeric(beta[x,])
  
  fit<-lm(beta ~ passage.or.rescope.no_numeric, data=sampleinfo_cpg)
  pval<-summary(fit)$coef["passage.or.rescope.no_numeric","Pr(>|t|)"]
  slope<-fit$coefficients[2]
  
  (min(sampleinfo_organoid_fetal$passage.or.rescope.no_numeric)*slope) - (max(sampleinfo_organoid_fetal$passage.or.rescope.no_numeric)*slope)})

passage_MTAB<-data.frame(p.value=ebfit$p.value[,"passage.or.rescope.no_numeric"], CpG=rownames(beta), db=passage_db)

# p adjust
passage_MTAB$p_adjusted<-p.adjust(passage_MTAB$p.value, method="BH")

diff_CpG_dbMTAB<-passage_MTAB[which(passage_MTAB$p_adjusted<0.05 & abs(passage_MTAB$db)>0.15),] #58958
diff_CpG_db_hypoMTAB<-diff_CpG_dbMTAB$CpG[which((diff_CpG_dbMTAB$db)>=0.15)] #  45098
diff_CpG_db_hyperMTAB<-diff_CpG_dbMTAB$CpG[which((diff_CpG_dbMTAB$db)<=(-0.15))] #  13860


#' How many differential CpGs could overlap? 450K vs EPIC
print(paste("Of the ", nrow(diff_CpG_db), " CpGs differential with passage in the original organoids, ",
            length(diff_CpG_db$CpG[which(diff_CpG_db$CpG%in%rownames(MTAB_organoid_beta_fetal))]), " are in the 450K data", sep=""))

diff_CpG_db_hypo_overlap<-diff_CpG_db_hypo[which(diff_CpG_db_hypo%in%rownames(MTAB_organoid_beta_fetal))]
diff_CpG_db_hyper_overlap<-diff_CpG_db_hyper[which(diff_CpG_db_hyper%in%rownames(MTAB_organoid_beta_fetal))]

diff_CpG_db_hypoMTAB_overlap<-diff_CpG_db_hypoMTAB[which(diff_CpG_db_hypoMTAB%in%pvals_long$CpG)]
diff_CpG_db_hyperMTAB_overlap<-diff_CpG_db_hyperMTAB[which(diff_CpG_db_hyperMTAB%in%pvals_long$CpG)]

print(paste("Of the ",length(diff_CpG_db_hypo_overlap)," hypo CpGs also on the 450K ",
            length(intersect(diff_CpG_db_hypoMTAB_overlap, diff_CpG_db_hypo_overlap))," are also hypo in the MTAB-4957 cohort (",
            round((length(intersect(diff_CpG_db_hypoMTAB_overlap, diff_CpG_db_hypo_overlap))/length(diff_CpG_db_hypo_overlap))*100,2),"%)",sep=""))

print(paste("Of the ",length(diff_CpG_db_hyper_overlap)," hyper CpGs also on the 450K ",
            length(intersect(diff_CpG_db_hyperMTAB_overlap, diff_CpG_db_hyper_overlap))," are also hypo in the MTAB-4957 cohort (",
            round((length(intersect(diff_CpG_db_hyperMTAB_overlap, diff_CpG_db_hyper_overlap))/length(diff_CpG_db_hyper_overlap))*100,2),"%)",sep=""))



### delta beta directionality plot
plt_db_direction<-merge(pvals_long[,c(3,6)], passage_MTAB, by="CpG")

plt_db_direction$sig<-"Not Significant"
plt_db_direction$sig[which(plt_db_direction$CpG%in%c(intersect(diff_CpG_db_hypoMTAB_overlap, diff_CpG_db_hypo_overlap),intersect(diff_CpG_db_hyperMTAB_overlap, diff_CpG_db_hyper_overlap)))]<-"Significant\nIn Both\nCohorts"

ggplot(plt_db_direction, aes(mean_db, db))+geom_point(aes(color=sig, alpha=sig),shape=19)+th+theme_bw()+
  scale_color_manual(values=c("lightgrey", "cornflowerblue"), name="Significant\nWith Passage")+
  scale_alpha_manual(values=c(0.25,1), guide=F)+
  geom_hline(yintercept=c(-0.15,0.15), color="grey60")+geom_vline(xintercept=c(-0.15,0.15), color="grey60")+
  ylim(-0.8,0.8)+xlim(-0.8,0.8)+xlab("Original Organoid\nPassage Delta Beta")+ylab("MTAB-4957 Fetal Organoid\nPassage Delta Beta")+
  stat_smooth(method="lm", se=F, color="black")


#ggsave(here("figs","MTAB_db_directionality_fetal.pdf"), width=5, height=3.75)
ggsave(here("figs/jpeg","MTAB_db_directionality_fetal.jpeg"), width=5, height=3.75)

print(paste("Correlation of delta betas between cohorts: ", round(cor(plt_db_direction$db, plt_db_direction$mean_db),2), sep=""))


#' ### representative CpGs
sampleinfo_organoid_fetal_minimal<-sampleinfo_organoid_fetal[,c(16,22,21)]
sampleinfo_organoid_fetal_minimal$cohort<-"MTAB-4957 Fetal Organoids"

sample_info_both<-rbind(sampleinfo_organoid_fetal_minimal,epic.organoid_minimal)

plt_hetero_MTAB<-function(CpGs, legend, axislab, title){
  betas<-melt(cbind(MTAB_organoid_beta_fetal[CpGs,],organoid_beta[CpGs,]))
  organoid_plt<-merge(sample_info_both, betas, by.x="Assay.Name",by.y="Var2")
  
  p<-ggplot(organoid_plt, aes(passage.or.rescope.no_numeric,value))+
    geom_line(aes(group=sample_ID),color="lightgrey")+
    stat_smooth(method="lm", color="grey30", size=0.7, se=F)+th+theme_bw()+
    geom_point(aes(fill=as.factor(passage.or.rescope.no_numeric)),shape=21, size=1.25)+
    scale_fill_manual(values=pass_col,name="Passage\nNumber", drop=T)+
    facet_grid(cohort~Var1)+
    ylab("DNAm Beta")+xlab("Passage Number")+ylim(0,1)+
    theme(plot.margin = margin(0.5, 0.15, 0.5, 0.15, "cm"),plot.title = element_text(size=12))+
    xlim(1,23)
  
  if(missing(legend) & missing(axislab) & missing(title)){p}else{
    if(legend=="N" & axislab=="N"){p + theme(legend.position = "none",axis.title.y=element_blank(),
                                             axis.text.y=element_blank(),
                                             axis.ticks.y=element_blank())+ ggtitle(title)}else{
                                               if(legend=="N" & axislab=="Y"){p + theme(legend.position = "none") + ggtitle(title)}}}}

plt_hetero_MTAB(c("cg25402228","cg09146328"))

ggsave(here("figs","Passage_differential_CpGs_MTAB4957_fetal.pdf"),width = 4.75, height = 4)
ggsave(here("figs/jpeg","Passage_differential_CpGs_MTAB4957_fetal.jpeg"), width = 4.75, height = 4)





#'## R Session Info
sessionInfo()
