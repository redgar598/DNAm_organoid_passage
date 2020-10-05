#'---
#'title: Enrichment of CpGs effected by passage in genome features
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' #### Load Libraries
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(lmtest))
suppressMessages(library(gridExtra))
suppressMessages(library(gtools))
suppressMessages(library(scales))
suppressMessages(library(grid))
suppressMessages(library(here))

options(stringsAsFactors = FALSE)
options(scipen = 999)


#' #### Load Functions
source(here("scripts","00_pretty_plots.R"))

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


#' ## Import differential DNAm and heteroskedacity statstics
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




#' ## IBD CpGs
IBD_CpG<-read.csv(here("data","Agliata_IBD_CpGs.csv"))
IBD_CpG<-IBD_CpG$probe


IBD_overlap<-function(CpG_list,background.probes,permutation_number){
  real_CpG_IBD<-length(intersect(IBD_CpG, CpG_list))
  bootstrap_IBD<-sapply(1:permutation_number, function(x){
    set.seed(x)
    Hit_number<-length(CpG_list)
    rnd_CpGs<-background.probes[sample(1:length(background.probes),Hit_number)]
    length(intersect(IBD_CpG, rnd_CpGs))
  })
  (length(which(bootstrap_IBD>=real_CpG_IBD))+1)/(permutation_number+1)}

IBD_overlap(diff_CpG_db_hypo, background_hypo$IlmnID, 1000)
IBD_overlap(diff_CpG_db_hyper, background_hyper$IlmnID, 1000)
IBD_overlap(hetero_CpG$IlmnID,background$IlmnID, 1000)


#' top CpGs from IBD study
plt_hetero(c("cg16465027","cg19269426","cg07839457","cg16240683"))






library(limma)
#' ## Differential methylation with samplesite 
mod<-model.matrix(~ 0 + sample.site, data=epic.organoid)
fit <- lmFit(organoid_beta, mod)
ebfit <- eBayes(fit)

cont.matrix <- makeContrasts(SCvsTI=sample.siteSC-sample.siteTI, levels=mod)
fit2 <- contrasts.fit(fit, cont.matrix)
ebfit2 <- eBayes(fit2)
topTable(ebfit2, adjust="BH")


# delta beta
samplesite_db<-sapply(1:nrow(organoid_beta), function(x){
  sampleinfo_cpg<-epic.organoid
  sampleinfo_cpg$beta<-as.numeric(organoid_beta[x,])
  
  mns<-tapply(sampleinfo_cpg$beta, sampleinfo_cpg$sample.site, mean)
  mns[1]-mns[2]
})

sample_site_CpGs<-data.frame(CpG=rownames(organoid_beta), p.value=ebfit2$p.value[,1], db=samplesite_db)


# Adjust P values
sample_site_CpGs$p_adjusted<-p.adjust(sample_site_CpGs$p.value, method="BH")
sig_sample_site_CpG<-sample_site_CpGs[which(sample_site_CpGs$p_adjusted<0.05 & abs(sample_site_CpGs$db)>0.15),] 


# hypo
length(intersect(diff_CpG_db_hypo, sig_sample_site_CpG$CpG))

expected_overlap<-sapply(1:1000, function(x){
  set.seed(x)
  rnd_cpg<-rownames(organoid_beta)[sample(1:nrow(organoid_beta), length(diff_CpG_db_hypo))]
  length(intersect(rnd_cpg,sig_sample_site_CpG$CpG))
})

(length(which(expected_overlap>length(intersect(diff_CpG_db_hypo, sig_sample_site_CpG$CpG))))+1)/1001


# hyper
length(intersect(diff_CpG_db_hyper, sig_sample_site_CpG$CpG))

expected_overlap<-sapply(1:1000, function(x){
  set.seed(x)
  rnd_cpg<-rownames(organoid_beta)[sample(1:nrow(organoid_beta), length(diff_CpG_db_hyper))]
  length(intersect(rnd_cpg,sig_sample_site_CpG$CpG))
})

(length(which(expected_overlap>length(intersect(diff_CpG_db_hyper, sig_sample_site_CpG$CpG))))+1)/1001



# hetero
length(intersect(hetero_CpG, sig_sample_site_CpG$CpG))

expected_overlap<-sapply(1:1000, function(x){
  set.seed(x)
  rnd_cpg<-rownames(organoid_beta)[sample(1:nrow(organoid_beta), length(hetero_CpG))]
  length(intersect(rnd_cpg,sig_sample_site_CpG$CpG))
})

(length(which(expected_overlap>length(intersect(hetero_CpG, sig_sample_site_CpG$CpG))))+1)/1001


plt_sample_site<-function(CpGs){
  betas<-melt(organoid_beta[CpGs,])
  epic.organoid_plt<-merge(epic.organoid, betas, by.x="array.id",by.y="X2")
  
  ggplot(epic.organoid_plt, aes(sample.site,value, fill=sample.site))+th+theme_bw()+
    geom_boxplot()+
    geom_point(aes(fill=as.factor(sample.site)),shape=21, size=1.25)+
    fillscale_sampsite+facet_wrap(~X1)+
    ylab("DNAm Beta")+ylim(0,1)
}

plt_sample_site(sig_sample_site_CpG$CpG[1:4])





#' ## c-DMR enrichment of passage CpGs
#' c-DMRs are from # https://www.nature.com/articles/ng.298#Sec25
#' supplementary data 2
#' The paper does not list the genome build but they DNAm was measured on the Nimblegen CHARM array. The GEO page for the Nimblegen CHARM array used does say is it hg18 (GRch36)

cDMR<-read.csv(here("data","EPIC_Features_cDMR36.csv"))
dim(cDMR)
#' Most CpGs don't have an associated cDMR, but here is an example with one.
cDMR[110:120,]

#' remove []
cDMR$X1<-gsub("\\[|\\]","", cDMR$X1)
print(paste("There are ",length(irizarry_sup_cDMR<-cDMR$X0[which(cDMR$X1!="")])," CpGs overlapping a cDMR",sep="") )
 
anno_EPIC<-anno_EPIC[which(anno_EPIC$IlmnID%in%cDMR$X0),]

#' Illumina annotation includes a CDMR column, so will compare to see my annotation is close
#' But they don't have direction, just yes no
illumina_cDMR<-anno_EPIC$IlmnID[which(anno_EPIC$DMR=="CDMR")] # 6581
print(paste("There are ",length(illumina_cDMR)," CpGs overlapping a cDMR",sep="") )
print(paste("Of the ",length(illumina_cDMR)," CpGs overlapping a cDMR, my annotation captures ",length(intersect(irizarry_sup_cDMR, illumina_cDMR))," of them.",sep="") )

#' Why do I have more? which of my cDMR CpG are not in cDMR according to illumina
x<-irizarry_sup_cDMR[which(!(irizarry_sup_cDMR%in%illumina_cDMR))]
nrow(x)
x<-anno_EPIC[which(anno_EPIC$IlmnID%in%x),c("IlmnID","DMR","Chromosome_36", "Coordinate_36")]
table(x$DMR)

#' Of the 721 cDMR CpG not called in cDMR by illumna
#' 713 of these are because it looks like illumina did not allow for multiple entries in the DMR column. So there is no CpG annotated as both RDMR CDMR and therefore there are less CDMR than should be based on irizarry supplement
#' And there are 8 cases where the CpG is on DMR boundary and I called as in and Illumina as out

#' tidy up the cDMR file
cDMR<-cDMR[,3:4]
colnames(cDMR)<-c("CpG", "cancer_db")
# 15 CpGs are in 2 cDMR. Fortunately in all cases the delta beta for these cDMRs are in the same direction so a mean will be taken
cDMR$cancer_db_mean<-sapply(1:nrow(cDMR), function(x){
  mean(as.numeric(strsplit(cDMR$cancer_db[x], ",")[[1]]))
})

#' the change in methylation are from the cDMR paper where "DeltaM is cancer minus normal"
cDMR$cDMR<-sapply(1:nrow(cDMR), function(x){
  if(is.finite(cDMR$cancer_db_mean[x])){
    if(cDMR$cancer_db_mean[x]>=0){"cDMR with more DNAm in cancer"}else{
      if(cDMR$cancer_db_mean[x]<=0){"cDMR with less DNAm in cancer"}}}else("Not cDMR")
})

# view the same CpG again
cDMR[110:120,]


#' Visulize the overall methylation level of cDMR across all samples
CpG_mnBeta<-data.frame(CpG = names(rowMeans(organoid_beta)), mnBeta = rowMeans(organoid_beta))

CpG_mnBeta_dmr<-merge(cDMR, CpG_mnBeta, by="CpG")

ggplot(CpG_mnBeta_dmr, aes(cDMR, mnBeta)) + geom_violin(fill="lightgrey", color="lightgrey", width=1)+
  geom_boxplot(aes(fill=cDMR), width=0.1, outlier.size=0.5)+th+theme_bw() + scale_fill_manual(values=c("#a6d96a","#92c5de","lightgrey"))+
  theme(legend.position="none",
        axis.text = element_text(size =12, color="black"),
        axis.title = element_text(size =15),
        legend.text = element_text(size =14),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust=1))

ggsave(here("figs","cDMR_CpGs_beta.pdf"),width = 4, height = 6)
ggsave(here("figs/jpeg","cDMR_CpGs_beta.jpeg"), width = 4, height = 6)


#' #### Monte carlo function on background of CpGs
DMR_permutation_enrichment<-function(CpG_list,background.probes, permutation_number, plotmin,plotmax, axislab){
  cDMR_DMR<-cDMR[which(cDMR$cDMR!="Not cDMR"),]
  cDMR_DMR$cDMR<-as.factor(cDMR_DMR$cDMR)
  
  cDMR_hits<-cDMR_DMR[which(cDMR_DMR$CpG%in%CpG_list),] 
  cDMR_hits_featureMeans<-tapply(cDMR_hits$CpG, cDMR_hits$cDMR, length)
  cDMR_hits_featureMeans[is.na(cDMR_hits_featureMeans)]<-0
  cDMR_hits_featureMeans<-data.frame(Probe_Count=as.numeric(cDMR_hits_featureMeans),
                                     Feature=names(cDMR_hits_featureMeans))
  cDMR_hits_featureMeans$Feature<-as.factor(cDMR_hits_featureMeans$Feature)

  # random sampling (to see if hits more in feature than expected)
  bootstrap_cDMR<-lapply(1:permutation_number, function(x){
    set.seed(x)
    Hit_number<-length(CpG_list)
    rnd_CpGs<-background.probes[sample(1:length(background.probes),Hit_number)]
    cDMR_rnd<-cDMR_DMR[which(cDMR_DMR$CpG%in%rnd_CpGs),]
    cDMR_rnd_featureMeans<-tapply(cDMR_rnd$CpG, cDMR_rnd$cDMR, length)
    cDMR_rnd_featureMeans[is.na(cDMR_rnd_featureMeans)]<-0
    cDMR_rnd_featureMeans
  })
  bootstrap_cDMR<-do.call(rbind,bootstrap_cDMR)
  
  cDMR_results<-lapply(1:nrow(cDMR_hits_featureMeans), function(x){
    real_CpG_in_region<-cDMR_hits_featureMeans$Probe_Count[x]
    Adjusted_enrich_p<-round(p.adjust((length(which(bootstrap_cDMR[,x]>=real_CpG_in_region))+1)/(permutation_number+1), method="fdr", n=nrow(cDMR_hits_featureMeans)),3)
    Adjusted_depletion_p<-round(p.adjust((length(which(bootstrap_cDMR[,x]<=real_CpG_in_region))+1)/(permutation_number+1), method="fdr", n=nrow(cDMR_hits_featureMeans)),3)
    print(paste("Feature: ", cDMR_hits_featureMeans$Feature[x], "   Enrichment: ", Adjusted_enrich_p, "; Depletion: ", Adjusted_depletion_p,  sep=""))
    
    data.frame(Feature = cDMR_hits_featureMeans$Feature[x], 
               Enrichment = Adjusted_enrich_p, 
               Depletion = Adjusted_depletion_p, 
               Sig = if(Adjusted_enrich_p<0.05 | Adjusted_depletion_p<0.05){"*"}else{""})
  })
  
  cDMR_results<-do.call(rbind, cDMR_results)
  
  # Plot the fold change
  cDMR_hits_featureMeans$Fold_change<-sapply(1:nrow(cDMR_hits_featureMeans), function(x) mean(foldchange(cDMR_hits_featureMeans$Probe_Count[x], bootstrap_cDMR[,x])))
  cDMR_hits_featureMeans$Fold_change[is.infinite(cDMR_hits_featureMeans$Fold_change)]<-NA
  
  se <- function(x) sd(x)/sqrt(length(x))
  cDMR_hits_featureMeans$Fold_change_se<-sapply(1:nrow(cDMR_hits_featureMeans), function(x) se(foldchange(cDMR_hits_featureMeans$Probe_Count[x], bootstrap_cDMR[,x])))
  cDMR_hits_featureMeans$Fold_change_se[is.infinite(cDMR_hits_featureMeans$Fold_change_se)]<-NA
  
  cDMR_hits_featureMeans<-merge(cDMR_hits_featureMeans, cDMR_results, by="Feature")
  
  cDMR_hits_featureMeans$Feature<-as.factor(cDMR_hits_featureMeans$Feature)
  levels(cDMR_hits_featureMeans$Feature)<-c("cDMR with less\nDNAm in cancer", "cDMR with more\nDNAm in cancer")
  
  p<-ggplot(cDMR_hits_featureMeans, aes(Feature, Fold_change, fill=Feature))+ # don't know what just DMR is
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
          scale_fill_manual(values=c("#a6d96a","#92c5de"))
  
  if(missing(axislab)){p}else{
    if(axislab=="N"){p + theme(axis.title.y=element_blank(),
                               axis.text.y=element_blank(),
                               axis.ticks.y=element_blank())}else{p}}
}


cDMR_enrichment_hetero<-DMR_permutation_enrichment(hetero_CpG$IlmnID,background$IlmnID, 1000, -20,20)
cDMR_enrichment_hetero
ggsave(here("figs","Passage_hetero_cDMR_CpGs_FDR.pdf"),width = 3, height = 6)
ggsave(here("figs/jpeg","Passage_hetero_cDMR_CpGs_FDR.jpeg"), width = 3, height = 6)

## with direction of effect backgrounds
cDMR_enrichment_diff_hypo<-DMR_permutation_enrichment(diff_CpG_hypo$IlmnID,background_hypo$IlmnID, 1000,-20,20)
cDMR_enrichment_diff_hypo
ggsave(here("figs","Passage_diffhypo_cDMR_CpGs_FDR_background015.pdf"),width = 3, height = 6)
ggsave(here("figs/jpeg","Passage_diffhypo_cDMR_CpGs_FDR_background015.jpeg"), width = 3, height = 6)
cDMR_enrichment_diff_hyper<-DMR_permutation_enrichment(diff_CpG_hyper$IlmnID,background_hyper$IlmnID, 1000,-20,20)
cDMR_enrichment_diff_hyper
ggsave(here("figs","Passage_diffhyper_cDMR_CpGs_FDR_background015.pdf"),width = 3, height = 6)
ggsave(here("figs/jpeg","Passage_diffhyper_cDMR_CpGs_FDR_background015.jpeg"), width = 3, height = 6)

t <- textGrob("* FDR\nadjusted\np < 0.05")

grid.arrange(DMR_permutation_enrichment(hetero_CpG$IlmnID,background$IlmnID, 1000, -20,20, "Y"),
             DMR_permutation_enrichment(diff_CpG_hypo$IlmnID,background_hypo$IlmnID, 1000,-20,20, "N"),
             DMR_permutation_enrichment(diff_CpG_hyper$IlmnID,background_hyper$IlmnID, 1000,-20,20, "N"), t, ncol=4, widths=c(1.3,1,1,0.5))

ggsave(here("figs","Passage_hetero_hypo_hyper_DMR_CpGs_FDR_background015.pdf"),
       grid.arrange(DMR_permutation_enrichment(hetero_CpG$IlmnID,background$IlmnID, 1000, -20,20, "Y"),
                    DMR_permutation_enrichment(diff_CpG_hypo$IlmnID,background_hypo$IlmnID, 1000,-20,20, "N"),
                    DMR_permutation_enrichment(diff_CpG_hyper$IlmnID,background_hyper$IlmnID, 1000,-20,20, "N"), t, ncol=4, widths=c(1.3,1,1,0.5)),width = 8, height = 3)


ggsave(here("figs/jpeg","Passage_hetero_hypo_hyper_DMR_CpGs_FDR_background015.jpeg"),
       grid.arrange(DMR_permutation_enrichment(hetero_CpG$IlmnID,background$IlmnID, 1000, -20,20, "Y"),
                    DMR_permutation_enrichment(diff_CpG_hypo$IlmnID,background_hypo$IlmnID, 1000,-20,20, "N"),
                    DMR_permutation_enrichment(diff_CpG_hyper$IlmnID,background_hyper$IlmnID, 1000,-20,20, "N"), t, ncol=4, widths=c(1.3,1,1,0.5)),width = 8, height = 3)

                     






#' ## Regulatory region enrichment
EPIC_ensembl_annotation<-read.csv(here("data","EPIC_ensembl_reg_annotation.csv"), header=T)
head(EPIC_ensembl_annotation)

EPIC_ensembl_annotation<-EPIC_ensembl_annotation[,3:4]
colnames(EPIC_ensembl_annotation)<-c("CpG", "Feature")

#' Remove special characters and make one feature-CpG per row
EPIC_ensembl_annotation_unique<-lapply(1:nrow(EPIC_ensembl_annotation), function(x) {
  Features=gsub("'","",strsplit(gsub("[^0-9A-Za-z///' ]","",as.character(EPIC_ensembl_annotation$Feature[x])), "' '")[[1]])
  if(length(Features)==0){Features<-"None"}else{}
  data.frame(CpG=EPIC_ensembl_annotation$CpG[x], Feature=Features)
})

EPIC_ensembl_annotation_unique<-do.call(rbind, EPIC_ensembl_annotation_unique)

head(EPIC_ensembl_annotation_unique)
print(paste("There are ",length(unique(EPIC_ensembl_annotation_unique[which(EPIC_ensembl_annotation_unique$Feature!="None"),]$CpG)),
            "/800383 CpGs in regulatory features, and ",sum(duplicated(EPIC_ensembl_annotation_unique$CpG))," in multiple features (mostly CTCF and Promoter)", sep=""))

#' ### Base EPIC array distribution
EPIC_featureCount<-tapply(EPIC_ensembl_annotation_unique$CpG, EPIC_ensembl_annotation_unique$Feature, length)
EPIC_featureCount<-data.frame(CpG_Count=as.numeric(EPIC_featureCount),
                              Feature=names(EPIC_featureCount))
EPIC_featureCount$Feature<-factor(EPIC_featureCount$Feature,
                                  levels=c("Open chromatin","TF binding site","CTCF Binding Site",
                                           "Enhancer","Promoter Flanking Region",
                                           "Promoter", "None"))
levels(EPIC_featureCount$Feature)<-c("Open\nchromatin","TF\nBinding\nsite","CTCF\nBinding\nSite",
                                     "Enhancer","Promoter\nFlanking\nRegion",
                                     "Promoter", "None")
EPIC_featureCount$percent<-(EPIC_featureCount$CpG_Count/sum(EPIC_featureCount$CpG_Count))*100

ggplot(EPIC_featureCount, aes(Feature, CpG_Count, fill=Feature))+
  geom_bar(position=position_dodge(width=0.9),stat="identity", color="grey25")+theme_bw()+
  ylab("CpG Count")+xlab("")+th+
  theme(legend.position="none",
        axis.text = element_text(size =10, color="black"),
        axis.title = element_text(size =15),
        legend.text = element_text(size =14),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))+
  geom_text(aes(label = paste (round(percent, 0), "%"), y=CpG_Count, vjust=-0.5), size=4, color="grey30")+
  scale_fill_manual(values=c("#d8d8d8","#cb95cc","#44dfcf","#fac900","#ff6969","#ff0000","grey97"))+
  scale_y_continuous(label=comma)

ggsave(here("figs","EPIC_distribution_reg_CpGs.pdf"),width = 6, height = 6)
ggsave(here("figs/jpeg","EPIC_distribution_reg_CpGs.jpeg"), width = 6, height = 6)



# Overall methylation level of these regions across all samples
# Are some regions depleted for hypo cause they are already unmethylated

CpG_mnBeta<-data.frame(CpG = names(rowMeans(organoid_beta)), mnBeta = rowMeans(organoid_beta))
CpG_mnBeta_feature<-merge(EPIC_ensembl_annotation_unique, CpG_mnBeta, by = "CpG")
CpG_mnBeta_feature$Feature<-factor(CpG_mnBeta_feature$Feature,
                                   levels=c("Open chromatin","TF binding site","CTCF Binding Site",
                                            "Enhancer","Promoter Flanking Region",
                                            "Promoter", "None"))
levels(CpG_mnBeta_feature$Feature)<-c("Open\nchromatin","TF\nBinding\nsite","CTCF\nBinding\nSite",
                                      "Enhancer","Promoter\nFlanking\nRegion",
                                      "Promoter", "None")

ggplot(CpG_mnBeta_feature, aes(Feature, mnBeta)) + geom_violin(fill="lightgrey", color="lightgrey", width=1.5)+
  geom_boxplot(aes(fill=Feature), width=0.1, outlier.size=0.5)+th+theme_bw() + scale_fill_manual(values=c("#44dfcf","#fac900","#d8d8d8","#ff0000","#ff6969","#cb95cc","grey97"))+
  theme(legend.position="none",
        axis.text = element_text(size =12, color="black"),
        axis.title = element_text(size =15),
        legend.text = element_text(size =14),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust=1))

ggsave(here("figs","EPIC_distribution_reg_CpGs_beta.pdf"),width = 8, height = 6)
ggsave(here("figs/jpeg","EPIC_distribution_reg_CpGs_beta.jpeg"), width = 8, height = 6)


#' #### Monte carlo function on background of CpGs
regulation_permutation_enrichment<-function(CpG_list,background.probes, permutation_number, plotmin,plotmax, axislab){
  
  reg_hits<-EPIC_ensembl_annotation_unique[which(EPIC_ensembl_annotation_unique$CpG%in%CpG_list),] 
  reg_hits_featureMeans<-tapply(reg_hits$CpG, reg_hits$Feature, length)
  reg_hits_featureMeans[is.na(reg_hits_featureMeans)]<-0
  reg_hits_featureMeans<-data.frame(Probe_Count=as.numeric(reg_hits_featureMeans),
                                    Feature=names(reg_hits_featureMeans))
  
  # random sampling (to see if hits more in feature than expected)
  bootstrap_reg<-lapply(1:permutation_number, function(x){
    set.seed(x)
    Hit_number<-length(CpG_list)
    rnd_CpGs<-background.probes[sample(1:length(background.probes),Hit_number)]
    reg_rnd<-EPIC_ensembl_annotation_unique[which(EPIC_ensembl_annotation_unique$CpG%in%rnd_CpGs),]
    reg_rnd_featureMeans<-tapply(reg_rnd$CpG, reg_rnd$Feature, length)
    reg_rnd_featureMeans[is.na(reg_rnd_featureMeans)]<-0
    reg_rnd_featureMeans
  })
  bootstrap_reg<-do.call(rbind,bootstrap_reg)
  
  reg_results<-lapply(1:nrow(reg_hits_featureMeans), function(x){
    real_CpG_in_region<-reg_hits_featureMeans$Probe_Count[x]
    
    Adjusted_enrich_p<-round(p.adjust((length(which(bootstrap_reg[,x]>=real_CpG_in_region))+1)/(permutation_number+1), method="fdr", n=nrow(reg_hits_featureMeans)),3)
    Adjusted_depletion_p<-round(p.adjust((length(which(bootstrap_reg[,x]<=real_CpG_in_region))+1)/(permutation_number+1), method="fdr", n=nrow(reg_hits_featureMeans)),3)
    print(paste("Feature: ", reg_hits_featureMeans$Feature[x], "   Enrichment: ", Adjusted_enrich_p, "; Depletion: ", Adjusted_depletion_p,  sep=""))
    
    data.frame(Feature = reg_hits_featureMeans$Feature[x], 
               Enrichment = Adjusted_enrich_p, 
               Depletion = Adjusted_depletion_p, 
               Sig = if(Adjusted_enrich_p<0.05 | Adjusted_depletion_p<0.05){"*"}else{""})
  })
  reg_results<-do.call(rbind, reg_results)
  
  ## fold change
  reg_hits_featureMeans$Fold_change<-sapply(1:nrow(reg_hits_featureMeans), function(x) mean(foldchange(reg_hits_featureMeans$Probe_Count[x], bootstrap_reg[,x])))
  reg_hits_featureMeans$Fold_change[is.infinite(reg_hits_featureMeans$Fold_change)]<-NA
  se <- function(x) sd(x)/sqrt(length(x))
  reg_hits_featureMeans$Fold_change_se<-sapply(1:nrow(reg_hits_featureMeans), function(x) se(foldchange(reg_hits_featureMeans$Probe_Count[x], bootstrap_reg[,x])))
  reg_hits_featureMeans$Fold_change_se[is.infinite(reg_hits_featureMeans$Fold_change_se)]<-NA
  reg_hits_featureMeans<-merge(reg_hits_featureMeans, reg_results, by="Feature")
  
  reg_hits_featureMeans$Feature<-factor(reg_hits_featureMeans$Feature,
                                        levels=c("Open chromatin","TF binding site","CTCF Binding Site",
                                                 "Enhancer","Promoter Flanking Region",
                                                 "Promoter", "None"))
  
  levels(reg_hits_featureMeans$Feature)<-c("Open chromatin","TF Binding site","CTCF Binding Site",
                                           "Enhancer","Promoter Flanking Region",
                                           "Promoter", "None")
  
  p<-ggplot(reg_hits_featureMeans, aes(Feature, Fold_change, fill=Feature))+ 
    geom_hline(yintercept = reg_hits_featureMeans$Fold_change[which(reg_hits_featureMeans$Feature=="None")], color="grey25", size=0.5)+
          geom_bar(position=position_dodge(width=0.9),stat="identity", color="grey25")+theme_bw()+
          geom_errorbar(aes(ymax = Fold_change+Fold_change_se, ymin=Fold_change-Fold_change_se),width=0.25,position=position_dodge(width=0.9))+
          ylab("Fold Change")+xlab("")+ylim(plotmin,plotmax)+
          theme(legend.position="none",plot.margin = margin(0.5, 0.1, 0.1, 0.5, "cm"),
                axis.title = element_text(size =12),
                legend.text = element_text(size =14),
                legend.title = element_text(size =12),
                strip.text.x = element_text(size = 12),
                axis.text.x = element_text(size = 10, angle = 45, hjust=1, color="black"))+
          geom_text(aes(label = Sig, y=2.5), size=6, color="grey40")+
          scale_fill_manual(values=c("#d8d8d8","#cb95cc","#44dfcf","#fac900","#ff6969","#ff0000","grey97"))#+ geom_text(aes(label = Probe_Count, y=Fold_change, vjust=-0.8, hjust=0.1),angle = 45, size=3, color="grey30")#
    
  
  if(missing(axislab)){p}else{
    if(axislab=="N"){p + theme(axis.title.y=element_blank(),
                               axis.text.y=element_blank(),
                               axis.ticks.y=element_blank())}else{p}}
}



reg_enrichment_hetero<-regulation_permutation_enrichment(hetero_CpG$IlmnID,background$IlmnID, 1000, -5,5)
reg_enrichment_hetero
ggsave(here("figs","Passage_hetero_reg_CpGs_FDR.pdf"),width = 3, height = 6)
ggsave(here("figs/jpeg","Passage_hetero_reg_CpGs_FDR.jpeg"), width = 3, height = 6)

reg_enrichment_diff_hypo<-regulation_permutation_enrichment(diff_CpG_hypo$IlmnID,background_hypo$IlmnID, 1000,-5,5)
reg_enrichment_diff_hypo
ggsave(here("figs","Passage_diff_hypo_reg_CpGs_FDR_background015.pdf"),width = 3, height = 6)
ggsave(here("figs/jpeg","Passage_diff_hypo_reg_CpGs_FDR_background015.jpeg"), width = 3, height = 6)

reg_enrichment_diff_hyper<-regulation_permutation_enrichment(diff_CpG_hyper$IlmnID,background_hyper$IlmnID, 1000,-5,5)
reg_enrichment_diff_hyper
ggsave(here("figs","Passage_diff_hyper_reg_CpGs_FDR_background015.pdf"),width = 3, height = 6)
ggsave(here("figs/jpeg","Passage_diff_hyper_reg_CpGs_FDR_background015.jpeg"), width = 3, height = 6)



## multipanel
t <- textGrob("* FDR\nadjusted\np < 0.05")



grid.arrange(regulation_permutation_enrichment(hetero_CpG$IlmnID,background$IlmnID, 1000, -3,3, "Y"), 
             regulation_permutation_enrichment(diff_CpG_hypo$IlmnID,background_hypo$IlmnID, 1000,-3,3, "N"), 
             regulation_permutation_enrichment(diff_CpG_hyper$IlmnID,background_hyper$IlmnID, 1000,-3,3, "N"), t, ncol=4, widths=c(1.3,1,1,0.5))


ggsave(here("figs","Passage_hetero_hypo_hyper_reg_CpGs_FDR_background015.pdf"),
       grid.arrange(regulation_permutation_enrichment(hetero_CpG$IlmnID,background$IlmnID, 1000, -3,3, "Y"), 
                    regulation_permutation_enrichment(diff_CpG_hypo$IlmnID,background_hypo$IlmnID, 1000,-3,3, "N"), 
                    regulation_permutation_enrichment(diff_CpG_hyper$IlmnID,background_hyper$IlmnID, 1000,-3,3, "N"), t, ncol=4, widths=c(1.3,1,1,0.5)),width = 8, height = 4)

ggsave(here("figs/jpeg","Passage_hetero_hypo_hyper_reg_CpGs_FDR_background015.jpeg"),
       grid.arrange(regulation_permutation_enrichment(hetero_CpG$IlmnID,background$IlmnID, 1000, -3,3, "Y"), 
                    regulation_permutation_enrichment(diff_CpG_hypo$IlmnID,background_hypo$IlmnID, 1000,-3,3, "N"), 
                    regulation_permutation_enrichment(diff_CpG_hyper$IlmnID,background_hyper$IlmnID, 1000,-3,3, "N"), t, ncol=4, widths=c(1.3,1,1,0.5)),width = 8, height = 4)




#'## R Session Info
sessionInfo()
