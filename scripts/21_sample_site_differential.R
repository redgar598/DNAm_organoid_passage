#'---
#'title: Exploration of CpGs different with sample site and effects of passage
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' #### Load Libraries
suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
library(limma)
suppressMessages(library(here))
library(cowplot)


options(stringsAsFactors = FALSE)
options(scipen = 999)

#' #### Load Functions
source(here("scripts","00_pretty_plots.R"))
source(here("scripts","00_EM_array_uniform_background_maximise_betabinom.R"))


#' #### Load Normalized Data
load(here("data","beta_organoids.RData"))

# convert to M values
Mval<-function(beta) log2(beta/(1-beta))
organoid_Mval = apply(organoid_beta, 1, Mval) # need mvalues for combat
organoid_Mval = as.data.frame(organoid_Mval)
organoid_Mval = t(organoid_Mval)


## model
print(paste("Testing in ", ncol(organoid_beta), " individuals", sep=""))
# p value
mod<-model.matrix(~as.factor(sample.site), data=epic.organoid)
fit <- lmFit(organoid_Mval, mod)
ebfit <- eBayes(fit)

delta_beta<-sapply(1:nrow(organoid_beta), function(x) {
  group_means<-tapply(organoid_beta[x,], epic.organoid$sample.site, mean)
  group_means[1]-group_means[2] # SC-TI
})

limmares_organoid <- data.frame(p.value=ebfit$p.value[,"as.factor(sample.site)TI"], delta_beta=delta_beta)
limmares_organoid$adjusted_p<-p.adjust(limmares_organoid$p.value, method="fdr", n=nrow(limmares_organoid))



organoid_hits<-limmares_organoid[which(limmares_organoid$adjusted_p<=0.0001),]
organoid_hits<-organoid_hits[which(abs(organoid_hits$delta_beta)>=0.2),]
nrow(organoid_hits)

organoid_hits_top<-organoid_hits[(order(abs(organoid_hits$p.value))),]
organoid_hits_top<-organoid_hits_top[1:500,]

organoid_beta_topdiff<-organoid_beta[which(rownames(organoid_beta)%in%rownames(organoid_hits_top)),]

heat_plot_df<-melt(organoid_beta_topdiff)
heat_plot_df<-merge(heat_plot_df, epic.organoid, by.x="X2", by.y="array.id")

heat_plot_df$X2<-factor(heat_plot_df$X2, levels=epic.organoid$array.id[order(epic.organoid$sample.site, epic.organoid$passage.or.rescope.no_numeric)])

db_top<-organoid_hits[which(rownames(organoid_hits)%in%rownames(organoid_beta_topdiff)),]
CpG_order<-rownames(db_top)[order(db_top$delta_beta)]

heat_plot_df$X1<-factor(heat_plot_df$X1, levels=CpG_order)

beta_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))

plot_grid(
  ggplot()+geom_tile(aes(X2, 1, fill=sample.site),heat_plot_df)+theme_void()+fillscale_sampsite,
  ggplot()+geom_tile(aes(X2, 1, fill=as.factor(passage.or.rescope.no_numeric)),heat_plot_df)+ 
    scale_fill_manual(values=pass_col,name="Passage\nNumber")+th+theme_void(),
  ggplot()+geom_tile(aes(X2, X1, fill=value),heat_plot_df)+
    ylab("CpG")+xlab("Individual")+
    scale_fill_gradientn(colours = beta_palette(100), limits=c(0, 1))+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()),
  ncol=1, align="v", rel_heights = c(1,1,15))



#### overlap with passage hits
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

diff_CpG_db_hypo<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)>=0.15)]
diff_CpG_db_hyper<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)<=(-0.15))]
print(paste("CpGs with significant (adjusted p<0.05; delta beta >0.15) hypomethylation: ", length(diff_CpG_db_hypo), sep=""))
print(paste("CpGs with significant (adjusted p<0.05; delta beta < -0.15) hypermethylation: ", length(diff_CpG_db_hyper), sep=""))


length(intersect(diff_CpG_db_hypo, rownames(organoid_hits)))
length(intersect(diff_CpG_db_hyper, rownames(organoid_hits)))
length(intersect(hetero_CpG, rownames(organoid_hits)))


## hypo heatplot
organoid_beta_topdiff<-organoid_beta[which(rownames(organoid_beta)%in%intersect(diff_CpG_db_hypo, rownames(organoid_hits))),]

heat_plot_df<-melt(organoid_beta_topdiff)
heat_plot_df<-merge(heat_plot_df, epic.organoid, by.x="X2", by.y="array.id")

heat_plot_df$X2<-factor(heat_plot_df$X2, levels=epic.organoid$array.id[order(epic.organoid$sample.site, epic.organoid$passage.or.rescope.no_numeric)])

db_top<-organoid_hits[which(rownames(organoid_hits)%in%rownames(organoid_beta_topdiff)),]
CpG_order<-rownames(db_top)[order(db_top$delta_beta)]

heat_plot_df$X1<-factor(heat_plot_df$X1, levels=CpG_order)

beta_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))

plot_grid(
  ggplot()+geom_tile(aes(X2, 1, fill=sample.site),heat_plot_df)+theme_void()+fillscale_sampsite,
  ggplot()+geom_tile(aes(X2, 1, fill=as.factor(passage.or.rescope.no_numeric)),heat_plot_df)+ 
    scale_fill_manual(values=pass_col,name="Passage\nNumber")+th+theme_void(),
  ggplot()+geom_tile(aes(X2, X1, fill=value),heat_plot_df)+
    ylab("CpG")+xlab("Individual")+
    scale_fill_gradientn(colours = beta_palette(100), limits=c(0, 1))+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()),
  ncol=1, align="v", rel_heights = c(1,1,15))



## hyper heatplot
organoid_beta_topdiff<-organoid_beta[which(rownames(organoid_beta)%in%intersect(diff_CpG_db_hyper, rownames(organoid_hits))),]

heat_plot_df<-melt(organoid_beta_topdiff)
heat_plot_df<-merge(heat_plot_df, epic.organoid, by.x="X2", by.y="array.id")

heat_plot_df$X2<-factor(heat_plot_df$X2, levels=epic.organoid$array.id[order(epic.organoid$sample.site, epic.organoid$passage.or.rescope.no_numeric)])

db_top<-organoid_hits[which(rownames(organoid_hits)%in%rownames(organoid_beta_topdiff)),]
CpG_order<-rownames(db_top)[order(db_top$delta_beta)]

heat_plot_df$X1<-factor(heat_plot_df$X1, levels=CpG_order)

beta_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))

plot_grid(
  ggplot()+geom_tile(aes(X2, 1, fill=sample.site),heat_plot_df)+theme_void()+fillscale_sampsite,
  ggplot()+geom_tile(aes(X2, 1, fill=as.factor(passage.or.rescope.no_numeric)),heat_plot_df)+ 
    scale_fill_manual(values=pass_col,name="Passage\nNumber")+th+theme_void(),
  ggplot()+geom_tile(aes(X2, X1, fill=value),heat_plot_df)+
    ylab("CpG")+xlab("Individual")+
    scale_fill_gradientn(colours = beta_palette(100), limits=c(0, 1))+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()),
  ncol=1, align="v", rel_heights = c(1,1,15))


## hetero heatplot
organoid_beta_topdiff<-organoid_beta[which(rownames(organoid_beta)%in%intersect(hetero_CpG, rownames(organoid_hits))),]

heat_plot_df<-melt(organoid_beta_topdiff)
heat_plot_df<-merge(heat_plot_df, epic.organoid, by.x="X2", by.y="array.id")

heat_plot_df$X2<-factor(heat_plot_df$X2, levels=epic.organoid$array.id[order(epic.organoid$sample.site, epic.organoid$passage.or.rescope.no_numeric)])

db_top<-organoid_hits[which(rownames(organoid_hits)%in%rownames(organoid_beta_topdiff)),]
CpG_order<-rownames(db_top)[order(db_top$delta_beta)]

heat_plot_df$X1<-factor(heat_plot_df$X1, levels=CpG_order)

beta_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))

plot_grid(
  ggplot()+geom_tile(aes(X2, 1, fill=sample.site),heat_plot_df)+theme_void()+fillscale_sampsite,
  ggplot()+geom_tile(aes(X2, 1, fill=as.factor(passage.or.rescope.no_numeric)),heat_plot_df)+ 
    scale_fill_manual(values=pass_col,name="Passage\nNumber")+th+theme_void(),
  ggplot()+geom_tile(aes(X2, X1, fill=value),heat_plot_df)+
    ylab("CpG")+xlab("Individual")+
    scale_fill_gradientn(colours = beta_palette(100), limits=c(0, 1))+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()),
  ncol=1, align="v", rel_heights = c(1,1,15))
