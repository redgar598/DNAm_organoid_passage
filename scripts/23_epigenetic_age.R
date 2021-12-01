#'---
#'title: Epigenetic age
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

suppressMessages(library(minfi))
suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
library(here)


source(here("scripts","00_pretty_plots.R"))
options(stringsAsFactors = FALSE)

# # filter to CpGs needed for clock and add probe IDs with NAs for missing 19 probes of 393 (5079 of dat mini)
#wget https://horvath.genetics.ucla.edu/html/dnamage/datMiniAnnotation.csv
datmini<-read.csv("../ibd/data/datMiniAnnotation.csv")


#' ### Load Data
# Load the cohort meta data Felicity Payne shared
path<-"data/original_organoids"
epic.organoid<-read.csv(here(path, "epic_organoid_MethylationArraySamples_31Jul19.csv"))

#' ### Normalize DNAm Arrays
here(path)
epic.organoid$array.id.path <- file.path(here(path), epic.organoid$array.id)
# multiple DMAP files common with epic so need to force https://support.bioconductor.org/p/97773/
rgset_organoid <- read.metharray(epic.organoid$array.id.path, verbose = FALSE,force=TRUE)

# Background and dye bias correction with noob thhrough funnorm implemented in minfi
#http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illumina <- preprocessFunnorm(rgset_organoid, sex=epic.organoid$sex)
organoid_beta<-as.data.frame(getBeta(MSet.illumina))
#' 
#' avg_detPval <- colMeans(detectionP(rgset_organoid))
#' epic.organoid$det_pval<-avg_detPval
#' 
#' print(paste("Number of samples: ", nrow(epic.organoid), sep=""))
#' epic.organoid<-epic.organoid[which(epic.organoid$det_pval<0.005),]
#' epic.organoid[(which(epic.organoid$det_pval>0.005)),]
#' 
#' #'Number of samples after removal of high detection p value: `r nrow(epic.organoid)`
#' #'
#' #' Normalize raw again but this time without bad samples
#' # multiple DMAP files common with epic so need to force https://support.bioconductor.org/p/97773/
#' rgset_organoid <- read.metharray(epic.organoid$array.id.path, verbose = FALSE,force=TRUE)
#' MSet.illumina <- preprocessFunnorm(rgset_organoid, sex=epic.organoid$sex)
#' organoid_beta<-as.data.frame(getBeta(MSet.illumina))

dim(organoid_beta)

organoid_beta<-organoid_beta[which(rownames(organoid_beta)%in%datmini$Name),]#24443
datNA<-datmini[which(!(datmini$Name%in%rownames(organoid_beta))),]
organoid_beta[datNA$Name,]<-NA
organoid_beta$ProbeID<-rownames(organoid_beta)
dim(organoid_beta)
organoid_beta<-organoid_beta[,c(83, 1:82)]

#format for clock
write.csv(organoid_beta, file="data/organoid_clock.csv", row.names = F, quote = F)

colnames(epic.organoid)[2]<-"Sample_ID"
colnames(epic.organoid)[10]<-"Age"
epic.organoid$Tissue<-"Gastric"
colnames(epic.organoid)[9]<-"Female"
epic.organoid$Female<-as.factor(epic.organoid$Female)
levels(epic.organoid$Female)<-c("1","0")
sampleinfo_min<-epic.organoid[,c("Sample_ID","Age","Female","Tissue")]
sampleinfo_min$Age<-as.numeric(sampleinfo_min$Age)
sampleinfo_min$Female<-as.character(sampleinfo_min$Female)

write.csv(sampleinfo_min, file="data/sampleinfo_organoid_clock_format.csv", row.names = F, quote = F)
# 
# # this was then put into https://dnamage.genetics.ucla.edu/new with clock and fun normalization and horvath normalization
# dnamage<-read.csv("../../output/clock/ibd_beta_organoid_funnorm.output.csv")
# dnamage$SampleID<-gsub("X","",dnamage$SampleID)
# 
# 
# # passage as numeric
# epic.organoid$passage.or.rescope.no_numeric<-as.factor(epic.organoid$passage.or.rescope.no)
# levels(epic.organoid$passage.or.rescope.no_numeric)<-c(1,10,11,14,16,2,3,4,6,7,8,2,4)
# epic.organoid$passage.or.rescope.no_numeric<-as.numeric(as.character(epic.organoid$passage.or.rescope.no_numeric))
# 
# 
# sampleinfo<-merge(dnamage, epic.organoid, by.x="Sample_ID", by.y="Sample_ID")
# 
# 
# 
# 
# #overall correlation plot
# ggplot(sampleinfo, aes(Age.x, DNAmAge))+geom_abline(color="lightgrey")+geom_point(size=2, color="black", fill="grey", shape=21)+
#   stat_smooth(method="lm", se=F, color="black")+th+theme_bw()+xlab("Chronological Age")+
#   annotate("text", x=20, y=5, label = paste("r = ", round(cor(sampleinfo$Age.x, sampleinfo$DNAmAge), 2)), size = 5)+
#   xlim(0,35)+ylim(0,35)
# 
# ggsave("../../figs/AgeCor_organoid.pdf", width = 4, height = 4)
# ggsave("../../figs/jpeg/AgeCor_organoid.jpeg", width = 4, height = 4)
# 
# mae(sampleinfo$Age.x, sampleinfo$DNAmAge)
# median(abs(sampleinfo$Age.x-sampleinfo$DNAmAge))
# 
# 
# ggplot(sampleinfo, aes(Age.x, DNAmAge))+geom_abline(color="lightgrey")+geom_point(aes(fill=as.factor(passage.or.rescope.no_numeric)),size=2, color="black",shape=21)+
#   stat_smooth(method="lm", se=F, color="black")+th+theme_bw()+xlab("Chronological Age")+
#   annotate("text", x=20, y=5, label = paste("r = ", round(cor(sampleinfo$Age.x, sampleinfo$DNAmAge), 2)), size = 5)+
#   xlim(0,35)+ylim(0,35)+scale_fill_brewer(palette = "Spectral",name="Passage\nNumber")
# 
# ggsave("../../figs/AgeCor_organoid_passage.pdf", width = 5, height = 4)
# ggsave("../../figs/jpeg/AgeCor_organoid_passage.jpeg", width = 5, height = 4)
# 
# 
# 
# # More simple grouping of samplesite
# sampleinfo$sample.site_grouped<-as.factor(sampleinfo$sample.site)
# levels(sampleinfo$sample.site_grouped)<-c("colon","ileum")
# site<-ggplot(sampleinfo, aes(sample.site_grouped, AgeAccelerationResidual, fill=sample.site_grouped))+geom_boxplot(color="black", outlier.shape=NA)+
#   geom_point(shape=21,position=position_jitter(w=0.15))+theme_bw()+scale_fill_manual(values=c("darkgreen","cornflowerblue"), name="Sample\nSite")+th+xlab("")+ylab("Age Acceleration Residual ")
# tapply(sampleinfo$AgeAccelerationResidual, sampleinfo$sample.site_grouped, mean)
# 
# summary(aov(sampleinfo$AgeAccelerationResidual~sampleinfo$sample.site_grouped))
# 
# sampleinfo$diagnosis_grouped<-as.factor(sampleinfo$diagnosis)
# levels(sampleinfo$diagnosis_grouped)<-c("CD","Control","IBD-U","Other.GI","UC")
# # 
# diagnosis<-ggplot(sampleinfo, aes(diagnosis_grouped, AgeAccelerationResidual, fill=diagnosis_grouped))+geom_boxplot(color="black", outlier.shape=NA)+
#   geom_point(shape=21,position=position_jitter(w=0.15))+theme_bw()+fillscale_diagnosis+th+xlab("")+ylab("Age Acceleration Residual")+
#   theme(axis.text.x=element_text(angle=45, hjust=1))+facet_wrap(~sample.site_grouped)
# 
# summary(aov(sampleinfo$AgeAccelerationResidual~sampleinfo$diagnosis_grouped))
# summary(aov(sampleinfo$AgeAccelerationResidual~sampleinfo$diagnosis_grouped+sampleinfo$sample.site_grouped))
# 
# 
# 
# passage<-ggplot(sampleinfo, aes(passage.or.rescope.no_numeric,AgeAccelerationResidual, 
#                                 fill=as.factor(passage.or.rescope.no_numeric)))+
#   geom_line(aes(passage.or.rescope.no_numeric,AgeAccelerationResidual, group=sample_ID), color="lightgrey")+
#   geom_point(shape=21, size=3)+scale_fill_brewer(palette = "Spectral",name="Passage\nNumber")+th+theme_bw()
# 
# cor(sampleinfo$passage.or.rescope.no_numeric,sampleinfo$AgeAccelerationResidual)
# 
# 
# ggsave("../../figs/Ageresid_organoid.pdf", grid.arrange(site, diagnosis, passage), width = 6, height = 12)
# ggsave("../../figs/jpeg/Ageresid_organoid.jpeg", grid.arrange(site, diagnosis, passage), width = 6, height = 12)
# 
# 
# 
# ####################################
# ### Organoids - PedClock
# ####################################
# 
# load("../../data/ibd_beta_organoids.RData")
# 
# dat0<-ibd_EPIC_organoid_beta
# 
# datM=t(dat0)
# colnames(datM)=as.character(rownames(dat0))
# anti.trafo= function(x,adult.age=20) {ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
# 
# datClock=read.csv("../../data/public_other/organoid/datcoefInteresting94.csv") #based on ALL samples
# selectCpGsClock=is.element(dimnames(datM)[[2]], as.character(datClock[,1][-1]))
# datMethClock0=data.frame(datM[,selectCpGsClock])
# datMethClock= data.frame(datMethClock0[as.character(datClock[,1][-1])])
# PedBE_age=as.numeric(anti.trafo(datClock[1,2]+as.numeric(as.matrix(datMethClock)%*%as.numeric(datClock[,2][-1]))))
# 
# identical(colnames(ibd_EPIC_organoid_beta), epic.organoid$array.id)
# 
# # passage as numeric
# epic.organoid$passage.or.rescope.no_numeric<-as.factor(epic.organoid$passage.or.rescope.no)
# levels(epic.organoid$passage.or.rescope.no_numeric)<-c(1,10,11,14,16,2,3,4,6,7,8,2,4)
# epic.organoid$passage.or.rescope.no_numeric<-as.numeric(as.character(epic.organoid$passage.or.rescope.no_numeric))
# 
# epic.organoid$PedBE_age<-PedBE_age
# epic.organoid$ped_age_acceleration<- residuals(lm(epic.organoid$PedBE_age~epic.organoid$age))
# 
# 
# 
# #overall correlation plot
# ggplot(epic.organoid, aes(age, PedBE_age))+geom_abline(color="lightgrey")+geom_point(size=2, color="black", fill="grey", shape=21)+
#   stat_smooth(method="lm", se=F, color="black")+th+theme_bw()+xlab("Chronological Age")+
#   annotate("text", x=20, y=5, label = paste("r = ", round(cor(epic.organoid$age, epic.organoid$PedBE_age), 2)), size = 5)+
#   xlim(0,35)+ylim(0,35)
# 
# ggsave("../../figs/AgeCor_organoid_PedBE.pdf", width = 4, height = 4)
# ggsave("../../figs/jpeg/AgeCor_organoid_PedBE.jpeg", width = 4, height = 4)
# 
# mae(epic.organoid$age, epic.organoid$PedBE_age)
# median(abs(epic.organoid$age-epic.organoid$PedBE_age))
# 
# 
# # More simple grouping of samplesite
# epic.organoid$sample.site_grouped<-as.factor(epic.organoid$sample.site)
# levels(epic.organoid$sample.site_grouped)<-c("colon","ileum")
# site<-ggplot(epic.organoid, aes(sample.site_grouped, ped_age_acceleration, fill=sample.site_grouped))+geom_boxplot(color="black", outlier.shape=NA)+
#   geom_point(shape=21,position=position_jitter(w=0.15))+theme_bw()+scale_fill_manual(values=c("darkgreen","cornflowerblue"), name="Sample\nSite")+th+xlab("")+ylab("Age Acceleration Residual ")
# tapply(epic.organoid$ped_age_acceleration, epic.organoid$sample.site_grouped, mean)
# 
# summary(aov(epic.organoid$ped_age_acceleration~epic.organoid$sample.site_grouped))
# 
# epic.organoid$diagnosis_grouped<-as.factor(epic.organoid$diagnosis)
# levels(epic.organoid$diagnosis_grouped)<-c("CD","Control","IBD-U","Other.GI","UC")
# # 
# diagnosis<-ggplot(epic.organoid, aes(diagnosis_grouped, ped_age_acceleration, fill=diagnosis_grouped))+geom_boxplot(color="black", outlier.shape=NA)+
#   geom_point(shape=21,position=position_jitter(w=0.15))+theme_bw()+fillscale_diagnosis+th+xlab("")+ylab("Age Acceleration Residual")+
#   theme(axis.text.x=element_text(angle=45, hjust=1))+facet_wrap(~sample.site_grouped)
# 
# summary(aov(epic.organoid$ped_age_acceleration~epic.organoid$diagnosis_grouped))
# summary(aov(epic.organoid$ped_age_acceleration~epic.organoid$diagnosis_grouped+epic.organoid$sample.site_grouped))
# 
# 
# 
# passage<-ggplot(epic.organoid, aes(passage.or.rescope.no_numeric,ped_age_acceleration, 
#                                    fill=as.factor(passage.or.rescope.no_numeric)))+
#   geom_line(aes(passage.or.rescope.no_numeric,ped_age_acceleration, group=sample_ID), color="lightgrey")+
#   geom_point(shape=21, size=3)+scale_fill_brewer(palette = "Spectral",name="Passage\nNumber")+th+theme_bw()
# 
# cor(epic.organoid$passage.or.rescope.no_numeric,epic.organoid$ped_age_acceleration)
# 
# ggsave("../../figs/Ageresid_organoid_PedBE.pdf", grid.arrange(site, diagnosis, passage), width = 6, height = 12)
# ggsave("../../figs/jpeg/Ageresid_organoid_PedBE.jpeg", grid.arrange(site, diagnosis, passage), width = 6, height = 12)
# 
# 
# 
