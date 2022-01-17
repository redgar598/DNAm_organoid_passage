suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(lmtest)
library(gridExtra)
library(gtools)
library(scales)
library(here)

source(here("scripts","00_pretty_plots.R"))

options(stringsAsFactors = FALSE)
options(scipen = 999)

load(here("data","beta_organoids.RData"))



###################
# Subset Sample site
###################
epic.organoid_TI<-epic.organoid[which(epic.organoid$sample.site=="TI"),]
organoid_TI<-organoid_beta[,which(colnames(organoid_beta)%in%epic.organoid_TI$array.id)]
identical(colnames(organoid_TI),epic.organoid_TI$array.id)

epic.organoid_SC<-epic.organoid[which(epic.organoid$sample.site=="SC"),]
organoid_SC<-organoid_beta[,which(colnames(organoid_beta)%in%epic.organoid_SC$array.id)]
identical(colnames(organoid_SC),epic.organoid_SC$array.id)




#################
# Heteroskedasticity with passage
#################

#' Plotting function for seeing the trend as individual CpGs
plt_hetero_segment<-function(CpGs, legend, axislab, title){
  betas<-melt(organoid_beta[CpGs,])
  epic.organoid_plt<-merge(epic.organoid, betas, by.x="array.id",by.y="X2")
  
  p<-ggplot(epic.organoid_plt, aes(passage.or.rescope.no_numeric,value))+
    geom_line(aes(group=sample_ID),color="lightgrey")+
    stat_smooth(method="lm", color="grey30", size=0.7, se=F)+th+theme_bw()+
    geom_point(aes(fill=as.factor(passage.or.rescope.no_numeric)),shape=21, size=1.25)+
    scale_fill_brewer(palette = "Spectral",name="Passage\nNumber")+facet_grid(sample.site~X1)+
    ylab("DNAm Beta")+xlab("Passage Number")+ylim(0,1)+
    theme(plot.margin = margin(0.5, 0.15, 0.5, 0.15, "cm"),plot.title = element_text(size=12))
  
  if(missing(legend) & missing(axislab) & missing(title)){p}else{
    if(legend=="N" & axislab=="N"){p + theme(legend.position = "none",axis.title.y=element_blank(),
                                             axis.text.y=element_blank(),
                                             axis.ticks.y=element_blank())+ ggtitle(title)}else{
                                               if(legend=="N" & axislab=="Y"){p + theme(legend.position = "none") + ggtitle(title)}}}}




###################
# Heteroskedacity from python
###################
#all organoids
load(here("data","Heteroskedactiy_pvalues_FDR_1000iter_w_CpG.Rdata"))
     
pvals_long$BP_pval<-((1000-pvals_long$BP_count)+1)/1001
pvals_long$diff_pval<-((1000-pvals_long$diff_count)+1)/1001
pvals_long$BP_fdr<-((1000-pvals_long$fdr_BP)+1)/1001
pvals_long$diff_fdr<-((1000-pvals_long$fdr_diff)+1)/1001


# split by segment
imports_pvals<-function(tissue){
  pvals_long<-read.csv(paste("data/Heteroskedactiy_pvalues_FDR_1000iter",tissue,".csv", sep=""), header=T)
  pvals_long[,1]<-NULL
  colnames(pvals_long)<-c("BP_count","diff_count","mean_db","fdr_BP","fdr_diff")
  pvals_long$CpG<-rownames(organoid_beta)
  
  pvals_long$BP_pval<-((1000-pvals_long$BP_count)+1)/1001
  pvals_long$diff_pval<-((1000-pvals_long$diff_count)+1)/1001
  pvals_long$BP_fdr<-((1000-pvals_long$fdr_BP)+1)/1001
  pvals_long$diff_fdr<-((1000-pvals_long$fdr_diff)+1)/1001
  pvals_long}
  
pvals_long_TI<-imports_pvals("TI")
pvals_long_SC<-imports_pvals("SC")

# combined
length(hetero_CpG<-rownames(organoid_beta)[which(pvals_long$BP_fdr<0.05)])
nrow(diff_CpG<-pvals_long[which(pvals_long$diff_fdr<0.05),])
nrow(diff_CpG_db<-pvals_long[which(pvals_long$diff_fdr<0.05 & abs(pvals_long$mean_db)>0.15),])
length(diff_CpG_db_hypo<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)>=0.15)]) #  17352
length(diff_CpG_db_hyper<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)<=(-0.15))]) #  6414


# TI
length(hetero_CpG_TI<-rownames(organoid_beta)[which(pvals_long_TI$BP_fdr<0.05)])#51474
nrow(diff_CpG_db_TI<-pvals_long_TI[which(pvals_long_TI$diff_fdr<0.05 & abs(pvals_long_TI$mean_db)>0.15),]) #33921
length(diff_CpG_db_hypo_TI<-diff_CpG_db_TI$CpG[which((diff_CpG_db_TI$mean_db)>=0.15)]) #  24103
length(diff_CpG_db_hyper_TI<-diff_CpG_db_TI$CpG[which((diff_CpG_db_TI$mean_db)<=(-0.15))]) #  9818

# SC
length(hetero_CpG_SC<-rownames(organoid_beta)[which(pvals_long_SC$BP_fdr<0.05)])#34416
nrow(diff_CpG_db_SC<-pvals_long_SC[which(pvals_long_SC$diff_fdr<0.05 & abs(pvals_long_SC$mean_db)>0.15),]) #11846
length(diff_CpG_db_hypo_SC<-diff_CpG_db_SC$CpG[which((diff_CpG_db_SC$mean_db)>=0.15)]) #  6923
length(diff_CpG_db_hyper_SC<-diff_CpG_db_SC$CpG[which((diff_CpG_db_SC$mean_db)<=(-0.15))]) #  4923
                      


length(intersect(hetero_CpG, hetero_CpG_TI))/41852

length(intersect(diff_CpG_db_hypo, diff_CpG_db_hypo_TI))/17352
length(intersect(diff_CpG_db_hypo, diff_CpG_db_hypo_SC))/17352
length(intersect(diff_CpG_db_hypo_TI, diff_CpG_db_hypo_SC))
length(intersect(intersect(diff_CpG_db_hypo_TI, diff_CpG_db_hypo_SC), diff_CpG_db_hypo))/17352


length(intersect(diff_CpG_db_hyper, diff_CpG_db_hyper_TI))/6414
length(intersect(diff_CpG_db_hyper, diff_CpG_db_hyper_SC))/6414
length(intersect(intersect(diff_CpG_db_hyper_TI, diff_CpG_db_hyper_SC), diff_CpG_db_hyper))/6414


## all differential
all_diff<-c(diff_CpG_db_hypo, diff_CpG_db_hyper)
length(unique(all_diff))

length(intersect(unique(c(diff_CpG_db_hyper_TI, diff_CpG_db_hypo_TI)), unique(all_diff)))
length(intersect(unique(c(diff_CpG_db_hyper_SC, diff_CpG_db_hypo_SC)), unique(all_diff)))

length(intersect(unique(c(diff_CpG_db_hyper_SC, diff_CpG_db_hypo_SC)), unique(c(diff_CpG_db_hyper_TI, diff_CpG_db_hypo_TI))))

length(unique(c(diff_CpG_db_hyper_SC, diff_CpG_db_hypo_SC)))
length(unique(c(diff_CpG_db_hyper_TI, diff_CpG_db_hypo_TI)))

## for venndiagram


plt_hetero_segment(c(intersect(diff_CpG_db_hyper_TI, diff_CpG_db_hyper_SC)[1:2],intersect(diff_CpG_db_hypo_TI, diff_CpG_db_hypo_SC)[1:2]))
# ggsave("../../figs/Passage_Heteroskedacity_CpGs_FDR.pdf",width = 4.25, height = 3.9)
# ggsave("../../figs/jpeg/Passage_Heteroskedacity_CpGs_FDR.jpeg", width = 4.25, height = 3.9)

