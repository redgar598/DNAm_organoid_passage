#'---
#'title: Validation organoid data
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


#' ### Load Data
path<-"data/validation/DNAm"
epic.organoid1<-read.csv(here(path, "METH_K_Nayak_20035_SS_Plate1.csv"), skip=7)
epic.organoid1$plate_path<-paste(epic.organoid1$Sample_Plate,"_96_Samples",sep="")
epic.organoid2<-read.csv(here(path, "METH_K_Nayak_20035_SS_Plate2.csv"), skip=7)
epic.organoid2$plate_path<-paste(epic.organoid2$Sample_Plate,"_96_samples",sep="")
epic.organoid3<-read.csv(here(path, "METH_K_Nayak_20035_SS_Plate3.csv"), skip=7)
epic.organoid3$plate_path<-paste(epic.organoid3$Sample_Plate,"_16_samples",sep="")

epic.organoid<-rbind(epic.organoid1,epic.organoid2,epic.organoid3)

#'### add gender
sampleInfo<-read.csv(here("data/validation", "Final samples for submission 2021.csv"))
epic.organoid<-merge(epic.organoid,sampleInfo[c("Sample.Name", "Age", "Diagnosis", "Gender", "GEPadGI", "BOX", "Segment","Source", "Wnt.type","Biobank.Rachel.Replicates")], by.x="Sample_Name", by.y="Sample.Name")

epic.organoid<-epic.organoid[which(epic.organoid$Biobank.Rachel.Replicates=="Rachel"),]

#' ### Normalize DNAm Arrays
here(path)
epic.organoid$array.id.path <- file.path(here(path,"Data"),epic.organoid$plate_path, epic.organoid$Sentrix_ID, paste(epic.organoid$Sentrix_ID, epic.organoid$Sentrix_Position, sep="_"))
epic.organoid$array.id<-paste(epic.organoid$Sentrix_ID, epic.organoid$Sentrix_Position, sep="_")
epic.organoid$individual<-sapply(1:nrow(epic.organoid), function(x) strsplit(epic.organoid$Sample_Name[x]," ")[[1]][1])

epic.organoid$passage<-as.numeric(sapply(1:nrow(epic.organoid), function(x) gsub("p","",strsplit(epic.organoid$Sample_Name[x]," ")[[1]][3])))
epic.organoid$passage_hilo<-sapply(1:nrow(epic.organoid), function(x) if(epic.organoid$passage[x]<5){"low"}else{"high"})
epic.organoid$condition<-sapply(1:nrow(epic.organoid), function(x) strsplit(epic.organoid$Sample_Name[x]," ")[[1]][4])
epic.organoid$condition<-as.factor(epic.organoid$condition)
levels(epic.organoid$condition)<-c("D","IFNg","IFNg","TNFa","UD","UT")
epic.organoid$condition<-as.character(epic.organoid$condition)
epic.organoid$comparison<-sapply(1:nrow(epic.organoid), function(x) if(epic.organoid$condition[x]%in%c("UD","D")){"differentiation"}else{"cytokine"})
epic.organoid$treatment<-sapply(1:nrow(epic.organoid), function(x) if(epic.organoid$comparison[x]=="cytokine"){epic.organoid$condition[x]}else{"UT"})
epic.organoid$differentiation<-sapply(1:nrow(epic.organoid), function(x) if(epic.organoid$comparison[x]=="differentiation"){epic.organoid$condition[x]}else{"UD"})



# multiple DMAP files common with epic so need to force https://support.bioconductor.org/p/97773/
rgset_organoid <- read.metharray(epic.organoid$array.id.path, verbose = FALSE,force=TRUE)

# Background and dye bias correction with noob thhrough funnorm implemented in minfi
#http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illumina <- preprocessFunnorm(rgset_organoid, sex=epic.organoid$Gender)
organoid_beta<-getBeta(MSet.illumina)

print(paste("Samples available: ",ncol(organoid_beta),"; Probes available: ",nrow(organoid_beta),sep=""))

#' ### Detection p values across all probes for each sample
avg_detPval <- colMeans(detectionP(rgset_organoid))
epic.organoid$det_pval<-avg_detPval

ggplot(epic.organoid)+geom_boxplot(aes(as.factor(Sentrix_ID), det_pval, fill=as.factor(Sentrix_ID)), outlier.shape = NA)+
  geom_point(aes(as.factor(Sentrix_ID), det_pval, group=Sample_Name, fill=as.factor(Sentrix_ID)), shape=21, color="black",
             position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE)

ggsave(here("figs","validation_detection_pvalue_organoids.pdf"), width=6, height=5)
ggsave(here("figs/jpeg","validation_detection_pvalue_organoids.jpeg"), width=6, height=5)




#' Beta distribution before and after normalization

# extract raw beta values for plotting
beta_raw<-getBeta(rgset_organoid)
identical(colnames(beta_raw),epic.organoid$array.id)

Beta_melted<- melt(organoid_beta)
Beta_melted_raw<- melt(beta_raw)

# Remove NAs before plotting (otherwise get many non-inifnite warnings)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
Beta_Plot_raw<-Beta_melted_raw[which(!(is.na(Beta_melted_raw$value))),]

# Add meta data
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid, by.x="ID", by.y="array.id")
colnames(Beta_Plot_raw)<-c("CpG","ID","Beta")
Beta_Plot_raw<-merge(Beta_Plot_raw,epic.organoid, by.x="ID", by.y="array.id")

beta_dis_EPIC_raw<-ggplot(Beta_Plot_raw, aes(Beta, group=as.character(ID), color=as.character(Segment)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm Beta Value")

beta_dis_EPIC<-ggplot(Beta_Plot, aes(Beta, group=as.character(ID), color=as.character(Segment)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm Beta Value")


grid.arrange(beta_dis_EPIC_raw,beta_dis_EPIC)

ggsave(here("figs","validation_beta_distribution_organoid.pdf"),grid.arrange(beta_dis_EPIC_raw,beta_dis_EPIC),  w=10, h=5)
ggsave(here("figs/jpeg","validation_beta_distribution_organoid.jpeg"),grid.arrange(beta_dis_EPIC_raw,beta_dis_EPIC), w=10, h=5)



#'#### Confirm individuals ID with SNPs probes and clustering by DNAm

# remove rows with NAs
Betas_cluster<-organoid_beta[complete.cases(organoid_beta),]

d <- dist(t(Betas_cluster))
hc <- hclust(d, method = "complete") #single, complete, average, ward
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(epic.organoid$Segment), cex=1.5)

pdf(here("figs","validation_cluster_wholeEPIC_organoid.pdf"), width=30)
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(epic.organoid$Segment), cex=1.5)
dev.off()


#' #### Genotyping Probes
SNPs <- getSnpBeta(rgset_organoid)
SNPs<-SNPs[complete.cases(SNPs),]# 65 cause one was all NA

SNPs<-SNPs[,which(colnames(SNPs)%in%epic.organoid$array.id)]
identical(colnames(SNPs),epic.organoid$array.id)

d <- dist(t(SNPs))
hc <- hclust(d, method = "complete") #single, complete, average, ward
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(as.character(epic.organoid$individual)), cex=1.5)

pdf(here("figs","validation_cluster_snps_EPIC_organoid.pdf"), width=30)
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(as.character(epic.organoid$individual)), cex=1.5)
dev.off()



#' #### Sex clustering
#' Using the cg ID to chromosome annotation from illumina 
#' https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html
anno_EPIC<-read.csv(here("data", "MethylationEPIC_v-1-0_B4.csv"), skip=7)

organoid_beta<-organoid_beta[which(rownames(organoid_beta)%in%anno_EPIC$IlmnID),]
anno_EPIC<-anno_EPIC[match(rownames(organoid_beta),anno_EPIC$IlmnID),]
identical(rownames(organoid_beta),anno_EPIC$IlmnID)

organoid_beta_sex<-organoid_beta[which(anno_EPIC$CHR%in%c('X','Y')),]

d <- dist(t(organoid_beta_sex))
hc <- hclust(d, method = "complete") #single, complete, average, ward
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(epic.organoid$Gender), cex=1.5)

pdf(here("figs","validation_cluster_sex_EPIC_organoid.pdf"), width=30)
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(epic.organoid$Gender), cex=1.5)
dev.off()



#' #### Remove samples which do not cluster correctly
#' none to remove here though

#' ### Probe Filtering 
# SNP probes should already be removed
organoid_beta <- organoid_beta[!grepl("rs",rownames(organoid_beta)), ]
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))

#' #### Sex Chromosomes 
anno_EPIC<-anno_EPIC[anno_EPIC$IlmnID%in%rownames(organoid_beta),]
identical(rownames(organoid_beta),anno_EPIC$IlmnID)
organoid_beta <- organoid_beta[!anno_EPIC$CHR%in%c("X", "Y"), ]

filt_sex<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))


#' #### Cross-hybridizing probes and polymorphic probes. 
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1
#' "43,254 cross-reactive probes with ≥ 47 bp homology with an off-target site, of which 15,782 (36.5 %) are new to the EPIC platform"
#' They include this annotated list in their supplement.
#' wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM2_ESM.csv
cross_reactive<-read.csv(here("data", "13059_2016_1066_MOESM2_ESM.csv"), stringsAsFactors = F)
organoid_beta<-organoid_beta[which(!(rownames(organoid_beta)%in%cross_reactive$PROBE)),]

filt_cross<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))


#'For polymorphic probes I will The Pidsley annotation aswell for "Probes overlapping genetic variants at targeted CpG sites." and "Probes overlapping genetic variants at single base extension sites for Infinium Type I probes" but NOT "Probes with genetic variants overlapping the body of the probe: 48 base pairs for Infinium Type I probes and 49 base pairs for Infinium Type II probes."

#wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM4_ESM.csv
#wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM5_ESM.csv

polymorphic<-read.csv(here("data", "13059_2016_1066_MOESM4_ESM.csv"), stringsAsFactors = F)
print(paste("Filtering ",length(unique(polymorphic$PROBE))," polymorphic probes (genetic variants at targeted CpG sites).", sep=""))

baseext<-read.csv(here("data", "13059_2016_1066_MOESM5_ESM.csv"), stringsAsFactors = F)
print(paste("Filtering ",length(unique(baseext$PROBE))," polymorphic probes (single base extension sites for Infinium Type I probes).", sep=""))

organoid_beta<-organoid_beta[which(!(rownames(organoid_beta)%in%c(polymorphic$PROBE, baseext$PROBE))),]

filt_poly<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))



#' #### Probe filtering based on detection pvalue and detection over background (NA)

#' Remove probes with high NA count
na_count_probe <-sapply(1:nrow(organoid_beta), function(y) length(which(is.na(organoid_beta[y,]))))
na_count_probe_good<-which(na_count_probe<(ncol(organoid_beta)*0.05))
organoid_beta<-organoid_beta[na_count_probe_good,]

filt_bead<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))


#' Remove probes with high detection p value across samples, and any samples with many high detection p value probes
detP <- detectionP(rgset_organoid)
detP<-detP[,which(colnames(detP)%in%epic.organoid$array.id)]
identical(colnames(detP),epic.organoid$array.id)

failed <- detP>0.05
bad_det_p<-names(which(rowMeans(failed)>0.01))
bad_det_psamp<-names(which(colMeans(failed)>0.01))

organoid_beta<-organoid_beta[which(!(rownames(organoid_beta)%in%bad_det_p)),]
organoid_beta<-organoid_beta[,which(!(colnames(organoid_beta)%in%bad_det_psamp))]
identical(colnames(organoid_beta), as.character(epic.organoid$array.id))

filt_detp<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))



#' #### Probe attrition plot
df<-data.frame(sample_num_remaining=c(866238,865918,865859,filt_sex,filt_cross,filt_poly,filt_bead,filt_detp),
               filter=c("EPIC Probe Number","Missing Annotation Data","Removal of SNP Probes",
                        "Removal of X and Y chromosome probes","Removal of Cross Reactive Probes",
                        "Removal of Polymorphic Probes", "Removal of Probes with Beadcount <3\nin 5 % of Samples",
                        "Removal of Probes with 1 % of samples\nwith a detection p-value greater than 0.05"))
df$sample_num_lost<-c(0,sapply(2:nrow(df), function(x) df$sample_num_remaining[x-1]-df$sample_num_remaining[x]))

df$filter<-factor(df$filter, rev(df$filter))

ggplot(df)+
  geom_bar(aes(filter,-sample_num_remaining), stat="identity", fill="grey70", color="black")+
  geom_bar(aes(filter,sample_num_lost), stat="identity",fill="darkred", color="black")+
  geom_text(aes(x=filter, y=-min(sample_num_remaining)/2,  label=comma(sample_num_remaining)))+
  geom_text(aes(x=filter, y=max(sample_num_lost)/1.5,  label=comma(sample_num_lost)))+
  geom_hline(yintercept=0)+
  coord_flip()+theme_bw()+ylab("")+xlab("")+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "grey20", size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(position = "top")

ggsave(here("figs","validation_probe_attrition.pdf"), width = 8, height = 3)
ggsave(here("figs/jpeg","validation_probe_attrition.jpeg"), width = 8, height = 3)

epic.organoid<-epic.organoid[,which(!(colnames(epic.organoid)%in%c("array.id.path","plate_path","Biobank.Rachel.Replicates","GEPadGI","BOX","Source","Wnt.type","Sample_Group", "Pool_ID")))]

validation_organoid_beta<-organoid_beta
validation_epic.organoid<-epic.organoid

save(validation_organoid_beta, validation_epic.organoid, file=here("data/validation/DNAm","validation_betas_normalized.RData"))








load(file=here("data/validation/DNAm","validation_betas_normalized.RData"))


#' ### Principal Component Analysis (PCA)
pca_res <- prcomp(t(validation_organoid_beta))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

validation_epic.organoid$Sentrix_ID<-as.factor(as.character(validation_epic.organoid$Sentrix_ID))

meta_categorical <- validation_epic.organoid[, c(3,4, 8,9,11,13,16,17)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(validation_epic.organoid[, c(6, 12,18)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Plate","Sentrix ID","Gender","Segment","Individual", "passage_hilo", "treatment", "differentiation")
colnames(meta_continuous) <- c( "age","passage_numeric", "det_pval")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8))

ggsave(here("figs/validation_heat_scree.pdf"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)
ggsave(here("figs/jpeg","validation_heat_scree.jpeg"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)




## PC vs PC plot
Loadings$array.id<-rownames(Loadings)
Loadings_meta<-merge(Loadings, validation_epic.organoid, by="array.id")

ggplot(Loadings_meta, aes(PC1, PC2, fill=Segment))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (42%)")+ylab("PC2 (12%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=c("#1a9850", "cornflowerblue","#e1e1ff","#1347a4","#000192"))
ggsave(here("figs","validation_PC12_site.pdf"), width = 5, height = 4)
ggsave(here("figs/jpeg","validation_PC12_site.jpeg"), width = 5, height = 4)




ggplot(Loadings_meta, aes(PC2, PC3, fill=passage_hilo))+
  geom_line(aes(PC2,PC3, group=individual, color=individual))+theme_bw()+
  geom_point(shape=21,size=3, color="black")+
  scale_color_manual(values=c("#d9d9d9","#525252","#969696","#737373","#bdbdbd"),name="Individual")+
  scale_fill_manual(values=c("#3288BD","#D53E4F"), name="Passage")+
  xlab("PC2 (12%)")+ylab("PC3 (11%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 1, 1, 1, "cm"))
ggsave(here("figs","validation_PC23_passage.pdf"), width = 5.5, height = 4)
ggsave(here("figs/jpeg","validation_PC23_passage.jpeg"), width = 5.5, height = 4)



ggplot(Loadings_meta, aes(PC2, PC3, fill=treatment, color=differentiation))+geom_point(shape=21,size=3)+
  scale_color_manual(values=c("black","white"))+scale_fill_manual(values=c("cornflowerblue","firebrick4","grey80"), name="Treatment")+
  theme_bw()+
  xlab("PC2 (12%)")+ylab("PC3 (11%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 1, 1, 1, "cm"))
ggsave(here("figs","validation_PC23_condition.pdf"), width = 6, height = 4)
ggsave(here("figs/jpeg","validation_PC23_condtion.jpeg"), width = 6, height = 4)





#' ## Overall Variance Across most Variable CpGs with Passage
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
Mval<-function(beta) log2(beta/(1-beta))

organoid_mval = apply(validation_organoid_beta, 1, Mval)
organoid_mval = as.data.frame(organoid_mval)
organoid_mval = t(organoid_mval)

ref_range_dnam<-sapply(1:nrow(organoid_mval), function(x) Variation(organoid_mval[x,]))
validation_organoid_beta_VeryVariable<-validation_organoid_beta[which(ref_range_dnam>=2.75),]

print(paste("There are ",nrow(validation_organoid_beta_VeryVariable), " variable CpGs (10th-90th quantile range in M value >2.75)",sep=""))


# beta plot variable CpGs
Beta_melted<- melt(validation_organoid_beta_VeryVariable)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,validation_epic.organoid, by.x="ID", by.y="array.id")
Beta_Plot$passage.or.numeric.factor <- factor(Beta_Plot$passage, levels = c(12,9,8,4,3,2))

ggplot(Beta_Plot, aes(Beta,  color=passage.or.numeric.factor))+
  geom_density(size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=rev(c("#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4","#4CA5B1","#5E4FA2")), name="Passage\nNumber")

ggsave(here("figs/validation_variable_CpGs.pdf"),width = 3.75, height = 2.5)
ggsave(here("figs/jpeg","validation_variable_CpGs.jpeg"), width = 7.5, height = 5)

ggplot(Beta_Plot, aes(Beta,  color=passage.or.numeric.factor, group=ID))+
  geom_density(size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+facet_grid(condition~passage.or.numeric.factor)+
  scale_color_manual(values=rev(c("#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4","#4CA5B1","#5E4FA2")), name="Passage\nNumber")

ggplot(Beta_Plot[which(Beta_Plot$differentiation=="UD" & Beta_Plot$treatment=="UT"),], aes(Beta,  color=passage.or.numeric.factor))+
  geom_density(size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=rev(c("#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4","#4CA5B1","#5E4FA2")), name="Passage\nNumber")



ggplot()+
  geom_density(aes(Beta,color=passage_hilo, group=ID),Beta_Plot, size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values = c ("#9ecae1", "#225ea8", "#081d58"), name="Relative\nPassage\nLevel within\nPatient")+facet_grid(condition~individual)+
  th+theme(strip.text = element_text(size = 10), axis.text=element_text(size=4),panel.spacing = unit(0.7, "lines"))+th+
  scale_x_continuous(breaks = c(0,0.5,1))

Beta_Plot$ID_nopass<-paste(Beta_Plot$individual, Beta_Plot$Segment,Beta_Plot$condition)
ggplot()+
  geom_density(aes(Beta,color=passage_hilo, group=ID),Beta_Plot, size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values = c ("#081d58","#9ecae1"), name="Relative\nPassage\nLevel within\nPatient")+facet_wrap(~ID_nopass)+
  th+theme(strip.text = element_text(size = 10), axis.text=element_text(size=4),panel.spacing = unit(0.7, "lines"))+th+
  scale_x_continuous(breaks = c(0,0.5,1))
ggsave(here("figs","validation_paired_beta.pdf"),width = 5, height = 2.2)
ggsave(here("figs/jpeg","validation_paired_beta.jpeg"),width = 10, height = 10)

ggplot()+
  geom_density(aes(Beta,color=passage.or.numeric.factor, group=ID),Beta_Plot, size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=rev(c("#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4","#4CA5B1","#5E4FA2")))+facet_wrap(~ID_nopass)+
  th+theme(strip.text = element_text(size = 10), axis.text=element_text(size=4),panel.spacing = unit(0.7, "lines"))+th+
  scale_x_continuous(breaks = c(0,0.5,1))

  
  

#' # Thresholding Trimodality
validation_epic.organoid$thresholded_prior_ratio<-sapply(1:nrow(validation_epic.organoid), function(x){
  print(x)
  converted<-as.numeric(round(validation_organoid_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  passage_threshold_params(converted, counts, res)
})

validation_epic.organoid$thresholded_ratio_max<-F
validation_epic.organoid$thresholded_ratio_max[which(validation_epic.organoid$thresholded_prior_ratio>1)]<-T

percent_passing<-round((tapply(validation_epic.organoid$thresholded_ratio_max, validation_epic.organoid$passage_hilo, sum)/tapply(validation_epic.organoid$array.id, validation_epic.organoid$passage_hilo, length))*100,2)
passed_num<-tapply(validation_epic.organoid$thresholded_ratio_max, validation_epic.organoid$passage_hilo, sum)
org_numer<-tapply(validation_epic.organoid$array.id, validation_epic.organoid$passage_hilo, length)

df<-data.frame(passage=names(percent_passing), passing=percent_passing, pro_passing=percent_passing/100, count=org_numer, passed_num=passed_num)
df<-cbind(df,(binom.confint(df$passed_num, df$count, method="exact", conf.level=0.95)))
df$upper<-df$upper*100
df$lower<-df$lower*100

print(df)



#' ### Plot all samples
plts_paired<-lapply(1:nrow(validation_epic.organoid), function(x){#1:nrow(validation_epic.organoid)
  print(x)
  passage<-paste("passage: ",validation_epic.organoid$passage[x],"\nIndividual: ", validation_epic.organoid$Sample_Name[x],"\nRatio I/H: ",round(validation_epic.organoid$thresholded_prior_ratio[x],2), sep="")
  converted<-as.numeric(round(validation_organoid_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  draw_fit_params_gg(converted, counts, res,passage)
})

plts_paired_order<-plts_paired[order(validation_epic.organoid$passage)]

pdf(here("figs","validation_organoids_thresholding_all_samples.pdf"))
plts_paired_order
dev.off()





#' ## Differential methylation with passage
mod<-model.matrix(~ 0 + passage, data=validation_epic.organoid)
fit <- lmFit(validation_organoid_beta, mod)
ebfit <- eBayes(fit)

# covariate adjusted beta values
beta<-validation_organoid_beta

passage_db<-sapply(1:nrow(beta), function(x){
  sampleinfo_cpg<-validation_epic.organoid
  sampleinfo_cpg$beta<-as.numeric(beta[x,])

  fit<-lm(beta ~ passage, data=sampleinfo_cpg)
  pval<-summary(fit)$coef["passage","Pr(>|t|)"]
  slope<-fit$coefficients[2]

  (min(validation_epic.organoid$passage)*slope) - (max(validation_epic.organoid$passage)*slope)})

passage_validation<-data.frame(p.value=ebfit$p.value[,"passage"], CpG=rownames(beta), db=passage_db)

# Adjust P values
passage_validation$p_adjusted<-p.adjust(passage_validation$p.value, method="BH")

diff_CpG_dbvalidation<-passage_validation[which(passage_validation$p_adjusted<0.05 & abs(passage_validation$db)>0.15),] #25086
diff_CpG_db_hypovalidation<-diff_CpG_dbvalidation$CpG[which((diff_CpG_dbvalidation$db)>=0.15)] #  30061
diff_CpG_db_hypervalidation<-diff_CpG_dbvalidation$CpG[which((diff_CpG_dbvalidation$db)<=(-0.15))] #  2859

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
            length(diff_CpG_db$CpG[which(diff_CpG_db$CpG%in%rownames(validation_organoid_beta))]), " are in the validation data", sep=""))

diff_CpG_db_hypo_overlap<-diff_CpG_db_hypo[which(diff_CpG_db_hypo%in%rownames(validation_organoid_beta))]
diff_CpG_db_hyper_overlap<-diff_CpG_db_hyper[which(diff_CpG_db_hyper%in%rownames(validation_organoid_beta))]

diff_CpG_db_hypovalidation_overlap<-diff_CpG_db_hypovalidation[which(diff_CpG_db_hypovalidation%in%pvals_long$CpG)]
diff_CpG_db_hypervalidation_overlap<-diff_CpG_db_hypervalidation[which(diff_CpG_db_hypervalidation%in%pvals_long$CpG)]

print(paste("Of the ",length(diff_CpG_db_hypo_overlap)," hypo CpGs also on the 450K ",
            length(intersect(diff_CpG_db_hypovalidation_overlap, diff_CpG_db_hypo_overlap))," are also hypo in the validation cohort (",
            round((length(intersect(diff_CpG_db_hypovalidation_overlap, diff_CpG_db_hypo_overlap))/length(diff_CpG_db_hypo_overlap))*100,2),"%)",sep=""))

print(paste("Of the ",length(diff_CpG_db_hyper_overlap)," hypo CpGs also on the 450K ",
            length(intersect(diff_CpG_db_hypervalidation_overlap, diff_CpG_db_hyper_overlap))," are also hypo in the validation cohort (",
            round((length(intersect(diff_CpG_db_hypervalidation_overlap, diff_CpG_db_hyper_overlap))/length(diff_CpG_db_hyper_overlap))*100,2),"%)",sep=""))


#' ### delta beta directionality plot
plt_db_direction<-merge(pvals_long[,c(3,6)], passage_validation, by="CpG")# 380776

plt_db_direction$sig<-"Not Significant"
plt_db_direction$sig[which(plt_db_direction$CpG%in%c(intersect(diff_CpG_db_hypovalidation_overlap, diff_CpG_db_hypo_overlap),intersect(diff_CpG_db_hypervalidation_overlap, diff_CpG_db_hyper_overlap)))]<-"Significant\nIn Both\nCohorts"

ggplot(plt_db_direction, aes(mean_db, db))+geom_point(aes(color=sig, alpha=sig),shape=19)+th+theme_bw()+
  scale_color_manual(values=c("lightgrey", "cornflowerblue"), name="Significant\nWith Passage")+
  scale_alpha_manual(values=c(0.25,1), guide=F)+
  geom_hline(yintercept=c(-0.15,0.15), color="grey60")+geom_vline(xintercept=c(-0.15,0.15), color="grey60")+
  ylim(-0.8,0.8)+xlim(-0.8,0.8)+xlab("Cohort 1 Organoid\nPassage Delta Beta")+ylab("Validation Organoid\nPassage Delta Beta")+
  stat_smooth(method="lm", se=F, color="black")


#ggsave(here("figs","validation_db_directionality.pdf"), width=5, height=3.75)
ggsave(here("figs/jpeg","validation_db_directionality.jpeg"), width=5, height=3.75)

print(paste("Correlation of delta betas between cohorts: ", round(cor(plt_db_direction$db, plt_db_direction$mean_db),2), sep=""))


#' ### representative CpGs
epic.organoid_minimal<-epic.organoid[,c(2, 14, 17)]
colnames(epic.organoid_minimal)[1]<-"Assay.Name"
epic.organoid_minimal$cohort<-"Cohort 1 Organoids"

validation_epic.organoid_minimal<-validation_epic.organoid[,c(10, 1, 12)]
colnames(validation_epic.organoid_minimal)[1]<-"Assay.Name"
validation_epic.organoid_minimal$cohort<-"Validation Organoids"
colnames(validation_epic.organoid_minimal)[2:3]<-c("sample_ID","passage.or.rescope.no_numeric")

sample_info_both<-rbind(validation_epic.organoid_minimal,epic.organoid_minimal)
sample_info_both$passage.or.rescope.no_numeric<-as.numeric(as.character(sample_info_both$passage.or.rescope.no_numeric))

plt_hetero_validation<-function(CpGs, legend, axislab, title){
  betas<-melt(cbind(validation_organoid_beta[CpGs,],organoid_beta[CpGs,]))
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

plt_hetero_validation(c("cg25402228","cg22009751"))

ggsave(here("figs","Passage_differential_CpGs_validation.pdf"),width = 4.75, height = 4)
ggsave(here("figs/jpeg","Passage_differential_CpGs_validation.jpeg"), width = 4.75, height = 4)



#' ### CpG to gene associations
EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194

diff_genes_db_hypovalidation<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypovalidation)] ) #11442
diff_genes_db_hypervalidation<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypervalidation)] ) # 2084

write.table(diff_genes_db_hypovalidation, file=here("data/validation/DNAm/","validation_genes_hypomethylation.txt"), quote=F, row.names = F, col.names = F)
write.table(diff_genes_db_hypervalidation, file=here("data/validation/DNAm/","validation_genes_hypermethylation.txt"), quote=F, row.names = F, col.names = F)

#'### Genes differential in original and validation
diff_genes_db_hypovalidation_original<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypo_overlap)] ) # 7875
diff_genes_db_hypervalidation_original<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hyper_overlap)] ) # 4281
write.table(diff_genes_db_hypovalidation_original, file=here("data/validation/DNAm/","validation_original_genes_hypomethylation.txt"), quote=F, row.names = F, col.names = F)
write.table(diff_genes_db_hypervalidation_original, file=here("data/validation/DNAm/","validation_original_genes_hypermethylation.txt"), quote=F, row.names = F, col.names = F)


#'## R Session Info
sessionInfo()
