suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(lmtest)
library(gridExtra)
library(gtools)
library(scales)

source("pretty_plots.R")

options(stringsAsFactors = FALSE)
options(scipen = 999)

load("../../data/ibd_beta_organoids.RData")
#Restructure meta
epic.organoid$sentrix_ID<-as.factor(epic.organoid$sentrix_ID)
epic.organoid$case.no<-as.factor(epic.organoid$case.no)
epic.organoid$inflammation<-as.factor(epic.organoid$inflammation)
epic.organoid$sex<-as.factor(epic.organoid$sex)
epic.organoid$sample.site<-as.factor(epic.organoid$sample.site)
epic.organoid$diagnosis<-as.factor(epic.organoid$diagnosis)

epic.organoid$passage.or.rescope.no_numeric<-as.factor(epic.organoid$passage.or.rescope.no)
levels(epic.organoid$passage.or.rescope.no_numeric)<-c(1,10,11,14,16,2,3,4,6,7,8,2,4)
epic.organoid$passage.or.rescope.no_numeric<-as.numeric(as.character(epic.organoid$passage.or.rescope.no_numeric))

identical(epic.organoid$array.id, colnames(ibd_EPIC_organoid_beta))


###################
# Subset Sample site
###################
epic.organoid_TI<-epic.organoid[which(epic.organoid$sample.site=="TI"),]
ibd_organoid_TI<-ibd_EPIC_organoid_beta[,which(colnames(ibd_EPIC_organoid_beta)%in%epic.organoid_TI$array.id)]
identical(colnames(ibd_organoid_TI),epic.organoid_TI$array.id)

epic.organoid_SC<-epic.organoid[which(epic.organoid$sample.site=="SC"),]
ibd_organoid_SC<-ibd_EPIC_organoid_beta[,which(colnames(ibd_EPIC_organoid_beta)%in%epic.organoid_SC$array.id)]
identical(colnames(ibd_organoid_SC),epic.organoid_SC$array.id)




#################
# Heteroskedasticity with passage
#################
## Down sampled > 4 passages to 3 samples

##plt and sig hetero

plt_hetero<-function(CpGs, legend, axislab, title){
  betas<-melt(ibd_EPIC_organoid_beta[CpGs,])
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


imports_pvals<-function(tissue){
  pvals_long<-read.csv(paste("../../output/Heteroskedactiy_pvalues_FDR_1000iter",tissue,".csv", sep=""), header=T)
  pvals_long[,1]<-NULL
  colnames(pvals_long)<-c("BP_count","diff_count","mean_db","fdr_BP","fdr_diff")
  pvals_long$CpG<-rownames(ibd_EPIC_organoid_beta)
  
  pvals_long$BP_pval<-((1000-pvals_long$BP_count)+1)/1001
  pvals_long$diff_pval<-((1000-pvals_long$diff_count)+1)/1001
  pvals_long$BP_fdr<-((1000-pvals_long$fdr_BP)+1)/1001
  pvals_long$diff_fdr<-((1000-pvals_long$fdr_diff)+1)/1001
  pvals_long}
  
pvals_long_TI<-imports_pvals("TI")
pvals_long_SC<-imports_pvals("SC")
pvals_long<-imports_pvals("")

# combined
length(hetero_CpG<-rownames(ibd_EPIC_organoid_beta)[which(pvals_long$BP_fdr<0.05)])#19490
nrow(diff_CpG<-pvals_long[which(pvals_long$diff_fdr<0.05),]) #19064
nrow(diff_CpG_db<-pvals_long[which(pvals_long$diff_fdr<0.05 & abs(pvals_long$mean_db)>0.15),]) #16986
length(diff_CpG_db_hypo<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)>=0.15)]) #  11772
length(diff_CpG_db_hyper<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)<=(-0.15))]) #  5214


# TI
length(hetero_CpG_TI<-rownames(ibd_EPIC_organoid_beta)[which(pvals_long_TI$BP_fdr<0.05)])#6305
nrow(diff_CpG_db_TI<-pvals_long_TI[which(pvals_long_TI$diff_fdr<0.05 & abs(pvals_long_TI$mean_db)>0.15),]) #18286
length(diff_CpG_db_hypo_TI<-diff_CpG_db_TI$CpG[which((diff_CpG_db_TI$mean_db)>=0.15)]) #  11802
length(diff_CpG_db_hyper_TI<-diff_CpG_db_TI$CpG[which((diff_CpG_db_TI$mean_db)<=(-0.15))]) #  6484

# SC
length(hetero_CpG_SC<-rownames(ibd_EPIC_organoid_beta)[which(pvals_long_SC$BP_fdr<0.05)])#0
nrow(diff_CpG_db_SC<-pvals_long_SC[which(pvals_long_SC$diff_fdr<0.05 & abs(pvals_long_SC$mean_db)>0.15),]) #3069
length(diff_CpG_db_hypo_SC<-diff_CpG_db_SC$CpG[which((diff_CpG_db_SC$mean_db)>=0.15)]) #  1310
length(diff_CpG_db_hyper_SC<-diff_CpG_db_SC$CpG[which((diff_CpG_db_SC$mean_db)<=(-0.15))]) #  1759
                      


length(intersect(hetero_CpG, hetero_CpG_TI))/6305

length(intersect(diff_CpG_db_hypo, diff_CpG_db_hypo_TI))/11802
length(intersect(diff_CpG_db_hypo, diff_CpG_db_hypo_SC))/1310
length(intersect(diff_CpG_db_hypo_TI, diff_CpG_db_hypo_SC))


length(intersect(diff_CpG_db_hyper, diff_CpG_db_hyper_TI))/6484
length(intersect(diff_CpG_db_hyper, diff_CpG_db_hyper_SC))/1759
length(intersect(diff_CpG_db_hyper_TI, diff_CpG_db_hyper_SC))



                      
plt_hetero(hetero_CpG[which(!(hetero_CpG%in%diff_CpG))][c(2,73,20,9)]) #73 2 20 9
# ggsave("../../figs/Passage_Heteroskedacity_CpGs_FDR.pdf",width = 4.25, height = 3.9)
# ggsave("../../figs/jpeg/Passage_Heteroskedacity_CpGs_FDR.jpeg", width = 4.25, height = 3.9)

# most differential but not hetero
diff_not_hetero<-diff_CpG_db[which(!(diff_CpG_db$CpG%in%hetero_CpG)),]#13152
diff_not_hetero_order<-diff_not_hetero[rev(order(abs(diff_not_hetero$mean_db))),]

plt_hetero(c("cg16152117","cg16182572","cg20869635","cg24157274"))
# ggsave("../../figs/Passage_differential_CpGs_FDR_hypo.pdf",width = 4.25, height = 3.9)
# ggsave("../../figs/jpeg/Passage_differential_CpGs_FDR_hypo.jpeg", width = 4.25, height = 3.9)


plt_hetero(c("cg03106584","cg09898122","cg22009751","cg24882644"))
# ggsave("../../figs/Passage_differential_CpGs_FDR_hyper.pdf",width = 4.25, height = 3.9)
# ggsave("../../figs/jpeg/Passage_differential_CpGs_FDR_hyper.jpeg", width = 4.25, height = 3.9)
                      # 
                      # 
                      # # mostly things become unmethylated with passage (delta beta are positive)
                      # ggplot(pvals_long, aes(mean_db))+geom_histogram(binwidth=0.025)+th+theme_bw()+geom_vline(xintercept=c(-0.2, 0.2))
                      # ggplot(pvals_long[which(abs(pvals_long$mean_db)>0.2),], aes(mean_db))+geom_histogram(binwidth=0.025)+th+theme_bw()+geom_vline(xintercept=c(-0.2, 0.2))
                      # 
                      # # grab the legend for all plots
                      # tmp <- ggplot_gtable(ggplot_build(plt_hetero(diff_not_hetero_order[1:4,"CpG"])))
                      # leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
                      # legend <- tmp$grobs[[leg]]
                      # 
                      # ggsave("../../figs/Passage_representative_CpGs_FDR.pdf",
                      #        grid.arrange(plt_hetero(hetero_CpG[which(!(hetero_CpG%in%diff_CpG))][c(2,66,14,8)],"N","Y","Heteroskedastic CpGs"),
                      #                     plt_hetero(diff_not_hetero_order[which(diff_not_hetero_order$mean_db>0),"CpG"][c(8,16,23,31)],"N","N","Hypomethylated CpGs"),
                      #                     plt_hetero(diff_not_hetero_order[which(diff_not_hetero_order$mean_db<0),"CpG"][c(1,2,22,32)],"N","N","Hypermethylated CpGs"),legend, nrow=1, widths=c(4.3,3.85,3.85,1)),
                      #        width = 9.75, height = 3.9)
                      # 
                      # ggsave("../../figs/jpeg/Passage_representative_CpGs_FDR.jpeg", 
                      #        grid.arrange(plt_hetero(hetero_CpG[which(!(hetero_CpG%in%diff_CpG))][c(2,66,14,8)],"N","Y"),
                      #                     plt_hetero(diff_not_hetero_order[which(diff_not_hetero_order$mean_db>0),"CpG"][c(8,16,23,31)],"N","N"),
                      #                     plt_hetero(diff_not_hetero_order[which(diff_not_hetero_order$mean_db<0),"CpG"][c(1,2,22,32)],"N","N"),legend, nrow=1, widths=c(4.3,3.85,3.85,1)),
                      #        width = 6.5, height = 2.5)
                      # 
                      # 

##########
## CpG Background to annotate in python
############

anno_EPIC<-read.csv("../../data/MethylationEPIC_v-1-0_B4.csv", skip=7)

background<-anno_EPIC[which(anno_EPIC$IlmnID%in%rownames(ibd_EPIC_organoid_beta)), c('IlmnID', 'CHR', 'MAPINFO')]
write.csv(background, file="../../output/passage_background.csv", row.names = F)

background_36<-anno_EPIC[which(anno_EPIC$IlmnID%in%rownames(ibd_EPIC_organoid_beta)), c('IlmnID', 'Chromosome_36', 'Coordinate_36')]
write.csv(background_36, file="../../output/passage_background_build36.csv", row.names = F)









