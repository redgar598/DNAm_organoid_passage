suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(lmtest)
library(gridExtra)
library(gtools)
library(scales)
library(here)


#' #### Load Functions
source(here("scripts","00_pretty_plots.R"))

#' #### Load Normalized Data
load(here("../ibd/data/","ibd_beta_organoids.RData"))

epic.organoid$passage.or.rescope.no_numeric<-as.factor(epic.organoid$passage.or.rescope.no)
levels(epic.organoid$passage.or.rescope.no_numeric)<-c(1,10,11,14,16,2,3,4,6,7,8,2,4)
epic.organoid$passage.or.rescope.no_numeric<-as.numeric(as.character(epic.organoid$passage.or.rescope.no_numeric))

organoid_beta<-ibd_EPIC_organoid_beta
#' Controls only
epic.organoid_ctrl<-epic.organoid[which(epic.organoid$diagnosis%in%c("Control","Other.GI")),]
organoid_beta_ctrl<-organoid_beta[,which(colnames(organoid_beta)%in%epic.organoid_ctrl$array.id)]

pass_col_ctrls<-c("#9E0142",  "#F46D43","#FEC776", "#E6F598", "#ABDDA4","#4CA5B1","#762a83")
names(pass_col_ctrls)<-c(1,2,3,4,6,7,8)


table(epic.organoid_ctrl$passage.or.rescope.no_numeric)

#' ### Principal Component Analysis (PCA)
pca_res <- prcomp(t(organoid_beta_ctrl))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)


#' ### PC vs PC plot
Loadings$array.id<-rownames(Loadings)
Loadings_meta<-merge(Loadings, epic.organoid_ctrl, by="array.id")

ggplot(Loadings_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=3, color="black")+theme_bw()+fillscale_sampsite+xlab("PC1 (29%)")+ylab("PC2 (9%)")+th+theme(axis.text = element_text(size=12),axis.title = element_text(size=14))

pc_plt<-ggplot(Loadings_meta, aes(PC2, PC3, fill=as.factor(passage.or.rescope.no_numeric)))+geom_line(aes(PC2,PC3, group=sample_ID), color="grey")+#, color=sampling.time.point
  geom_point(shape=21,size=3, color="black")+#
  theme_bw()+xlab("PC2 (9%)")+ylab("PC3 (3%)")+th+theme(axis.text = element_text(size=12),axis.title = element_text(size=14))+
  scale_fill_manual(values=pass_col_ctrls,name="Passage\nNumber")#+ #scale_color_manual(values=c("white", "black"))

legend<-ggplot(epic.organoid_ctrl, aes(as.factor(-passage.or.rescope.no_numeric), fill=as.factor(passage.or.rescope.no_numeric)))+geom_bar(color="black")+
  theme_bw()+theme(legend.position = "none", axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   legend.title=element_text(size=10),
                   legend.text=element_text(size=8))+
  coord_flip()+
  scale_fill_manual(values=pass_col_ctrls)+ylab("Sample\nNumber")+th

r <- ggplot() + theme_void()

grid.arrange(pc_plt,arrangeGrob(r,legend,r, heights=c(0.6,1.25,0.4)), ncol=2, widths=c(7,1))


#' ###### Overall Variance Across most Variable CpGs with Passage
Mval<-function(beta) log2(beta/(1-beta))
organoid_Mval = apply(organoid_beta_ctrl, 1, Mval) # need mvalues for combat
organoid_Mval = as.data.frame(organoid_Mval)
organoid_Mval = t(organoid_Mval)

# variable mvalues only
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
ref_range_dnam<-sapply(1:nrow(organoid_Mval), function(x) Variation(organoid_Mval[x,]))
organoid_beta_VeryVariable<-organoid_beta_ctrl[which(ref_range_dnam>=2.75),]#  51545

print(paste("There are ",nrow(organoid_beta_VeryVariable), " variable CpGs (10th-90th quantile range in M value >2.75)",sep=""))

# beta plot variable CpGs
Beta_melted<- melt(organoid_beta_VeryVariable)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid_ctrl, by.x="ID", by.y="array.id")
Beta_Plot$passage.or.rescope.no_numeric.factor <- factor(Beta_Plot$passage.or.rescope.no_numeric, levels = c(16,14,11,10,8,7,6,4,3,2,1))





ggplot(Beta_Plot, aes(Beta,color=passage.or.rescope.no_numeric.factor))+
  geom_density(size=0.5)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=pass_col_ctrls,name="Passage\nNumber")+th+theme(legend.text = element_text(size=7),
                                                                            legend.title = element_text(size=10),
                                                                            legend.key.size = unit(0.7,"line"))




#' ##### Beta distributions paired by individual
epic.organoid_paired<-epic.organoid_ctrl[which(epic.organoid_ctrl$sample_ID%in%epic.organoid_ctrl$sample_ID[duplicated(epic.organoid_ctrl$sample_ID)]),]

#'Paired samples for beta distributions: `r unique(epic.organoid_paired$case.no)`
#'From passage: `r unique(epic.organoid_paired$passage.or.rescope.no_numeric)`


