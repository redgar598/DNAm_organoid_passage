#'---
#'title: Exploration of CpGs effected by passage
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' #### Load Libraries
suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(lmtest)
library(gridExtra)
library(gtools)
library(scales)
library(here)
library(binom)


options(stringsAsFactors = FALSE)
options(scipen = 999)


#' #### Load Functions
source(here("scripts","00_pretty_plots.R"))
source(here("scripts","00_EM_array_uniform_background_maximise_betabinom.R"))


#' #### Load Normalized Data
load(here("data","beta_organoids.RData"))



#' Plotting function for seeing the trend as individual CpGs
plt_hetero<-function(CpGs, legend, axislab, title){
  betas<-melt(organoid_beta[CpGs,])
  epic.organoid_plt<-merge(epic.organoid, betas, by.x="array.id",by.y="X2")
  
  p<-ggplot(epic.organoid_plt, aes(passage.or.rescope.no_numeric,value))+
    geom_line(aes(group=sample_ID),color="lightgrey")+
    stat_smooth(method="lm", color="grey30", size=0.7, se=F)+th+theme_bw()+
    geom_point(aes(fill=as.factor(passage.or.rescope.no_numeric)),shape=21, size=1.25)+
    scale_fill_brewer(palette = "Spectral",name="Passage\nNumber")+facet_wrap(~X1)+
    ylab("DNAm Beta")+xlab("Passage Number")+ylim(0,1)+
    theme(plot.margin = margin(0.5, 0.15, 0.5, 0.15, "cm"),plot.title = element_text(size=12))
  
  if(missing(legend) & missing(axislab) & missing(title)){p}else{
    if(legend=="N" & axislab=="N"){p + theme(legend.position = "none",axis.title.y=element_blank(),
                                             axis.text.y=element_blank(),
                                             axis.ticks.y=element_blank())+ ggtitle(title)}else{
                                               if(legend=="N" & axislab=="Y"){p + theme(legend.position = "none") + ggtitle(title)}}}}




#' ### Import differential DNAm and heteroskedacity statstics
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
print(paste("CpGs with significant (adjusted p<0.05; absolute delta beta >0.15) differential methylation: ", nrow(diff_CpG_db), sep=""))

#' ##### Plot representative CpGs
#' First CpGs which are significantly heteroskedastic but not differentially methylated
plt_hetero(c("cg16867657"))

EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194
plt_hetero(unique(EPIC_genes[which(EPIC_genes$Gene.name=="ELOVL2"),]$IlmnID))
EPIC_genes[which(EPIC_genes$IlmnID=="cg16867657"),]

plt_hetero(unique(EPIC_genes[which(EPIC_genes$Gene.name=="ELOVL2"),]$IlmnID))
pvals_long[which(pvals_long$CpG%in%unique(EPIC_genes[which(EPIC_genes$Gene.name=="ELOVL2"),]$IlmnID)),]

EPIC_genes[which(EPIC_genes$Gene.name=="EVL" & EPIC_genes$CpG_in=="promoter"),]
plt_hetero(unique(EPIC_genes[which(EPIC_genes$Gene.name=="EVL" & EPIC_genes$CpG_in=="promoter"),]$IlmnID))
pvals_long[which(pvals_long$CpG%in%unique(EPIC_genes[which(EPIC_genes$Gene.name=="EVL" & EPIC_genes$CpG_in=="promoter"),]$IlmnID)),]
plt_hetero(c("cg02314758","cg02121312","cg27060175"))

EPIC_genes[which(EPIC_genes$Gene.name=="ESR1"),]
plt_hetero(unique(EPIC_genes[which(EPIC_genes$Gene.name=="ESR1" & EPIC_genes$CpG_in=="promoter"),]$IlmnID))
pvals_long[which(pvals_long$CpG%in%unique(EPIC_genes[which(EPIC_genes$Gene.name=="ESR1"),]$IlmnID)),]



#' ## Enrichment in aging lists
#' ### Lists of CpGs to test for enrichment
#' Using the cg ID to chromosome and coordinte annotation from illumina
#' https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html
anno_EPIC<-read.csv(here("data", "MethylationEPIC_v-1-0_B4.csv"), skip=7)

background<-anno_EPIC[which(anno_EPIC$IlmnID%in%rownames(organoid_beta)), c('IlmnID', 'CHR', 'MAPINFO')]
hetero_CpG<-anno_EPIC[which(anno_EPIC$IlmnID%in%hetero_CpG), c('IlmnID', 'CHR', 'MAPINFO')] #14867

# split diff into hypo and hyper
diff_CpG_db_hypo<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)>=0.15)] #  11772
diff_CpG_db_hyper<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)<=(-0.15))] #  5214
diff_CpG_hypo<-anno_EPIC[which(anno_EPIC$IlmnID%in%diff_CpG_db_hypo), c('IlmnID', 'CHR', 'MAPINFO')] #
diff_CpG_hyper<-anno_EPIC[which(anno_EPIC$IlmnID%in%diff_CpG_db_hyper), c('IlmnID', 'CHR', 'MAPINFO')] #



#' ### Backgorund CpG List to compare to
#' CpGs which have an passage=0 intercept <0.15 need to be excluded from the hypomethylayed background. Since they start <0.15 they can decrease by >0.15 as there is a lower limit of 0.
#' CpGs which have an passage=0 intercept >0.85 need to be excluded from the hypermethylayed background. Since they start >0.15 they can increase by >0.15 as there is a upper limit of 1.

intercepts<-sapply(1:nrow(organoid_beta), function(x) as.numeric(lm(organoid_beta[x,]~epic.organoid$passage.or.rescope.no_numeric)$coefficients[1]))

rm_hypo<-rownames(organoid_beta)[which(intercepts<=0.15)]
rm_hyper<-rownames(organoid_beta)[which(intercepts>=0.85)]

length(rm_hypo)
length(rm_hyper)

plt_hetero(rm_hypo[1:4])
plt_hetero(rm_hyper[1:4])

background_hypo<-background[which(!(background$IlmnID%in%rm_hypo)),]
background_hyper<-background[which(!(background$IlmnID%in%rm_hyper)),]

#' Load age CpGs
ageCG<-read.csv(here("data", "Day_etal_gen_biol_2013_Additional_file_2.csv"))

ageCG$direction<-"positive"
ageCG$direction[which(ageCG$slope.age<0)]<-"negative"
ageCG$tissue_direction<-paste(ageCG$direction, ageCG$TISSUE)

#' ### Permutation function
age_permutation_enrichment<-function(CpG_list,background.probes, permutation_number, plotmin,plotmax, axislab){
  ageCG$tissue_direction<-as.factor(ageCG$tissue_direction)
  
  ageCG_hits<-ageCG[which(ageCG$probe%in%CpG_list),] 
  ageCG_hits_featureMeans<-tapply(ageCG_hits$probe, ageCG_hits$tissue_direction, length)
  ageCG_hits_featureMeans[is.na(ageCG_hits_featureMeans)]<-0
  ageCG_hits_featureMeans<-data.frame(Probe_Count=as.numeric(ageCG_hits_featureMeans),
                                     Feature=names(ageCG_hits_featureMeans))
  ageCG_hits_featureMeans$Feature<-as.factor(ageCG_hits_featureMeans$Feature)
  
  # random sampling (to see if hits more in feature than expected)
  bootstrap_age<-lapply(1:permutation_number, function(x){
    set.seed(x)
    Hit_number<-length(CpG_list)
    rnd_CpGs<-background.probes[sample(1:length(background.probes),Hit_number)]
    ageCG_rnd<-ageCG[which(ageCG$probe%in%rnd_CpGs),]
    ageCG_rnd_featureMeans<-tapply(ageCG_rnd$probe, ageCG_rnd$tissue_direction, length)
    ageCG_rnd_featureMeans[is.na(ageCG_rnd_featureMeans)]<-0
    ageCG_rnd_featureMeans
  })
  bootstrap_age<-do.call(rbind,bootstrap_age)
  

  age_results<-lapply(1:nrow(ageCG_hits_featureMeans), function(x){
    real_CpG_in_region<-ageCG_hits_featureMeans$Probe_Count[x]
    Adjusted_enrich_p<-round(p.adjust((length(which(bootstrap_age[,x]>=real_CpG_in_region))+1)/(permutation_number+1), method="fdr", n=nrow(ageCG_hits_featureMeans)),3)
    Adjusted_depletion_p<-round(p.adjust((length(which(bootstrap_age[,x]<=real_CpG_in_region))+1)/(permutation_number+1), method="fdr", n=nrow(ageCG_hits_featureMeans)),3)
    print(paste("Feature: ", ageCG_hits_featureMeans$Feature[x], "   Enrichment: ", Adjusted_enrich_p, "; Depletion: ", Adjusted_depletion_p,  sep=""))
    
    data.frame(Feature = ageCG_hits_featureMeans$Feature[x], 
               Enrichment = Adjusted_enrich_p, 
               Depletion = Adjusted_depletion_p, 
               Sig = if(Adjusted_enrich_p<0.05 | Adjusted_depletion_p<0.05){"*"}else{""})
  })
  
  age_results<-do.call(rbind, age_results)
  
  # Plot the fold change
  ageCG_hits_featureMeans$Fold_change<-sapply(1:nrow(ageCG_hits_featureMeans), function(x) mean(foldchange(ageCG_hits_featureMeans$Probe_Count[x]+1, bootstrap_age[,x]+1)))# +1 cause there are 0
  ageCG_hits_featureMeans$Fold_change[is.infinite(ageCG_hits_featureMeans$Fold_change)]<-NA
  
  se <- function(x) sd(x)/sqrt(length(x))
  ageCG_hits_featureMeans$Fold_change_se<-sapply(1:nrow(ageCG_hits_featureMeans), function(x) se(foldchange(ageCG_hits_featureMeans$Probe_Count[x]+1, bootstrap_age[,x]+1)))
  ageCG_hits_featureMeans$Fold_change_se[is.infinite(ageCG_hits_featureMeans$Fold_change_se)]<-NA
  
  ageCG_hits_featureMeans<-merge(ageCG_hits_featureMeans, age_results, by="Feature")
  
  ageCG_hits_featureMeans$Feature<-as.factor(ageCG_hits_featureMeans$Feature)

  p<-ggplot(ageCG_hits_featureMeans, aes(Feature, Fold_change, fill=Feature))+ # don't know what just DMR is
    geom_bar(position=position_dodge(width=0.9),stat="identity", color="grey25")+theme_bw()+
    geom_errorbar(aes(ymax = Fold_change+Fold_change_se, ymin=Fold_change-Fold_change_se),width=0.25,position=position_dodge(width=0.9))+
    ylab("Fold Change")+xlab("")+ylim(plotmin,plotmax)+
    theme(legend.position="none",plot.margin = margin(0.5, 0.1, 0.1, 0.5, "cm"),
          axis.text = element_text(size =10, angle=45,color="black", hjust=1),
          axis.title = element_text(size =12),
          legend.text = element_text(size =12),
          legend.title = element_text(size =12),
          strip.text.x = element_text(size = 12))+
    geom_text(aes(label = Sig, y=7), size=6, color="grey40")+
    scale_fill_manual(values=c("#a6d96a","#92c5de","blue","pink","#a6d96a","#92c5de","blue","pink"))
  
  if(missing(axislab)){p}else{
    if(axislab=="N"){p + theme(axis.title.y=element_blank(),
                               axis.text.y=element_blank(),
                               axis.ticks.y=element_blank())}else{p}}
}


age_enrichment_hetero<-age_permutation_enrichment(hetero_CpG$IlmnID,background$IlmnID, 1000, -20,20)
age_enrichment_hetero
#ggsave(here("figs","Passage_hetero_cDMR_CpGs_FDR.pdf"),width = 3, height = 6)
#ggsave(here("figs/jpeg","Passage_hetero_cDMR_CpGs_FDR.jpeg"), width = 3, height = 6)

## with direction of effect backgrounds
age_enrichment_diff_hypo<-age_permutation_enrichment(diff_CpG_hypo$IlmnID,background_hypo$IlmnID, 1000,-20,20)
age_enrichment_diff_hypo
#ggsave(here("figs","Passage_diffhypo_cDMR_CpGs_FDR_background015.pdf"),width = 3, height = 6)
#ggsave(here("figs/jpeg","Passage_diffhypo_cDMR_CpGs_FDR_background015.jpeg"), width = 3, height = 6)
age_enrichment_diff_hyper<-age_permutation_enrichment(diff_CpG_hyper$IlmnID,background_hyper$IlmnID, 1000,-20,20)
age_enrichment_diff_hyper
#ggsave(here("figs","Passage_diffhyper_cDMR_CpGs_FDR_background015.pdf"),width = 3, height = 6)
#ggsave(here("figs/jpeg","Passage_diffhyper_cDMR_CpGs_FDR_background015.jpeg"), width = 3, height = 6)



## Age data Heyn
age_dmr<-read.csv(here("data","Heyn_etal_PNAS_2010_data_S2.csv"), skip=3)

#' ### Permutation function
agedmr_permutation_enrichment<-function(CpG_list,background.probes, permutation_number, plotmin,plotmax, axislab){

  ageCG_hits<-age_dmr[which(age_dmr$Target.ID%in%CpG_list),] 
  ageCG_hits_featureMeans<-tapply(ageCG_hits$Target.ID, ageCG_hits$Direction.in.Y103, length)
  ageCG_hits_featureMeans[is.na(ageCG_hits_featureMeans)]<-0
  ageCG_hits_featureMeans<-data.frame(Probe_Count=as.numeric(ageCG_hits_featureMeans),
                                      Feature=names(ageCG_hits_featureMeans))
  ageCG_hits_featureMeans$Feature<-as.factor(ageCG_hits_featureMeans$Feature)
  
  # random sampling (to see if hits more in feature than expected)
  bootstrap_age<-lapply(1:permutation_number, function(x){
    set.seed(x)
    Hit_number<-length(CpG_list)
    rnd_CpGs<-background.probes[sample(1:length(background.probes),Hit_number)]
    ageCG_rnd<-age_dmr[which(age_dmr$Target.ID%in%rnd_CpGs),]
    ageCG_rnd_featureMeans<-tapply(ageCG_rnd$Target.ID, ageCG_rnd$Direction.in.Y103, length)
    ageCG_rnd_featureMeans[is.na(ageCG_rnd_featureMeans)]<-0
    ageCG_rnd_featureMeans
  })
  bootstrap_age<-do.call(rbind,bootstrap_age)
  
  
  age_results<-lapply(1:nrow(ageCG_hits_featureMeans), function(x){
    real_CpG_in_region<-ageCG_hits_featureMeans$Probe_Count[x]
    Adjusted_enrich_p<-round(p.adjust((length(which(bootstrap_age[,x]>=real_CpG_in_region))+1)/(permutation_number+1), method="fdr", n=nrow(ageCG_hits_featureMeans)),3)
    Adjusted_depletion_p<-round(p.adjust((length(which(bootstrap_age[,x]<=real_CpG_in_region))+1)/(permutation_number+1), method="fdr", n=nrow(ageCG_hits_featureMeans)),3)
    print(paste("Feature: ", ageCG_hits_featureMeans$Feature[x], "   Enrichment: ", Adjusted_enrich_p, "; Depletion: ", Adjusted_depletion_p,  sep=""))
    
    data.frame(Feature = ageCG_hits_featureMeans$Feature[x], 
               Enrichment = Adjusted_enrich_p, 
               Depletion = Adjusted_depletion_p, 
               Sig = if(Adjusted_enrich_p<0.05 | Adjusted_depletion_p<0.05){"*"}else{""})
  })
  
  age_results<-do.call(rbind, age_results)
  
  # Plot the fold change
  ageCG_hits_featureMeans$Fold_change<-sapply(1:nrow(ageCG_hits_featureMeans), function(x) mean(foldchange(ageCG_hits_featureMeans$Probe_Count[x]+1, bootstrap_age[,x]+1)))# +1 cause there are 0
  ageCG_hits_featureMeans$Fold_change[is.infinite(ageCG_hits_featureMeans$Fold_change)]<-NA
  
  se <- function(x) sd(x)/sqrt(length(x))
  ageCG_hits_featureMeans$Fold_change_se<-sapply(1:nrow(ageCG_hits_featureMeans), function(x) se(foldchange(ageCG_hits_featureMeans$Probe_Count[x]+1, bootstrap_age[,x]+1)))
  ageCG_hits_featureMeans$Fold_change_se[is.infinite(ageCG_hits_featureMeans$Fold_change_se)]<-NA
  
  ageCG_hits_featureMeans<-merge(ageCG_hits_featureMeans, age_results, by="Feature")
  
  ageCG_hits_featureMeans$Feature<-as.factor(ageCG_hits_featureMeans$Feature)
  
  p<-ggplot(ageCG_hits_featureMeans, aes(Feature, Fold_change, fill=Feature))+ # don't know what just DMR is
    geom_bar(position=position_dodge(width=0.9),stat="identity", color="grey25")+theme_bw()+
    geom_errorbar(aes(ymax = Fold_change+Fold_change_se, ymin=Fold_change-Fold_change_se),width=0.25,position=position_dodge(width=0.9))+
    ylab("Fold Change")+xlab("")+ylim(plotmin,plotmax)+
    theme(legend.position="none",plot.margin = margin(0.5, 0.1, 0.1, 0.5, "cm"),
          axis.text = element_text(size =10, angle=45,color="black", hjust=1),
          axis.title = element_text(size =12),
          legend.text = element_text(size =12),
          legend.title = element_text(size =12),
          strip.text.x = element_text(size = 12))+
    geom_text(aes(label = Sig, y=7), size=6, color="grey40")+
    scale_fill_manual(values=c("#a6d96a","#92c5de","blue","pink"))
  
  if(missing(axislab)){p}else{
    if(axislab=="N"){p + theme(axis.title.y=element_blank(),
                               axis.text.y=element_blank(),
                               axis.ticks.y=element_blank())}else{p}}
}


age_enrichment_hetero<-agedmr_permutation_enrichment(hetero_CpG$IlmnID,background$IlmnID, 1000, -20,20)
age_enrichment_hetero
#ggsave(here("figs","Passage_hetero_cDMR_CpGs_FDR.pdf"),width = 3, height = 6)
#ggsave(here("figs/jpeg","Passage_hetero_cDMR_CpGs_FDR.jpeg"), width = 3, height = 6)

## with direction of effect backgrounds
age_enrichment_diff_hypo<-agedmr_permutation_enrichment(diff_CpG_hypo$IlmnID,background_hypo$IlmnID, 1000,-20,20)
age_enrichment_diff_hypo
#ggsave(here("figs","Passage_diffhypo_cDMR_CpGs_FDR_background015.pdf"),width = 3, height = 6)
#ggsave(here("figs/jpeg","Passage_diffhypo_cDMR_CpGs_FDR_background015.jpeg"), width = 3, height = 6)
age_enrichment_diff_hyper<-agedmr_permutation_enrichment(diff_CpG_hyper$IlmnID,background_hyper$IlmnID, 1000,-20,20)
age_enrichment_diff_hyper
#ggsave(here("figs","Passage_diffhyper_cDMR_CpGs_FDR_background015.pdf"),width = 3, height = 6)
#ggsave(here("figs/jpeg","Passage_diffhyper_cDMR_CpGs_FDR_background015.jpeg"), width = 3, height = 6)

