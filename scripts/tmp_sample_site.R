#'---
#'title: Sample Site CpGs
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ### Load Libraries
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(RColorBrewer))
suppressMessages(library(here))
suppressMessages(library(limma))


suppressMessages(library(dplyr))
suppressMessages(library(lmtest))
suppressMessages(library(gridExtra))
suppressMessages(library(gtools))


options(stringsAsFactors = FALSE)


#' ### Load Functions
source(here("scripts","00_pretty_plots.R"))

#' #### Load Normalized Data
load(here("data","beta_organoids.RData"))


#' ## Differential methylation with passage 
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

#Passage sites
pvals_long<-read.csv(here("data","Heteroskedactiy_pvalues_FDR_1000iter.csv"), header=T)
pvals_long[,1]<-NULL
colnames(pvals_long)<-c("BP_count","diff_count","mean_db","fdr_BP","fdr_diff")
pvals_long$CpG<-rownames(organoid_beta)
pvals_long$BP_pval<-((1000-pvals_long$BP_count)+1)/1001
pvals_long$diff_pval<-((1000-pvals_long$diff_count)+1)/1001
pvals_long$BP_fdr<-((1000-pvals_long$fdr_BP)+1)/1001
pvals_long$diff_fdr<-((1000-pvals_long$fdr_diff)+1)/1001

hetero_CpG<-rownames(organoid_beta)[which(pvals_long$BP_fdr<0.05)]
diff_CpG<-rownames(organoid_beta)[which(pvals_long$diff_fdr<0.05)]
diff_CpG_db<-pvals_long[which(pvals_long$diff_fdr<0.05 & abs(pvals_long$mean_db)>0.15),]

# differential
length(intersect(diff_CpG_db$CpG, sig_sample_site_CpG$CpG))

expected_overlap<-sapply(1:1000, function(x){
  set.seed(x)
  rnd_cpg<-rownames(organoid_beta)[sample(1:nrow(organoid_beta), length(diff_CpG_db$CpG))]
  length(intersect(rnd_cpg,sig_sample_site_CpG$CpG))
})

(length(which(expected_overlap>length(intersect(diff_CpG_db$CpG, sig_sample_site_CpG$CpG))))+1)/1001


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
