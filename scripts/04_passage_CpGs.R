#'---
#'title: Exploration of CpGs effected by passage
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
suppressMessages(library(here))
suppressMessages(library(binom))


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

save(pvals_long, file=here("data","Heteroskedactiy_pvalues_FDR_1000iter_w_CpG.Rdata"))

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
plt_hetero(hetero_CpG[which(!(hetero_CpG%in%diff_CpG))][c(31,37,132,18)])


diff_CpG_db_hypo<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)>=0.15)]
diff_CpG_db_hyper<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)<=(-0.15))]
print(paste("CpGs with significant (adjusted p<0.05; delta beta >0.15) hypomethylation: ", length(diff_CpG_db_hypo), sep=""))
print(paste("CpGs with significant (adjusted p<0.05; delta beta < -0.15) hypermethylation: ", length(diff_CpG_db_hyper), sep=""))


#' Then CpGs which are significantly differentially methylated but not heteroskedastic
diff_not_hetero<-diff_CpG_db[which(!(diff_CpG_db$CpG%in%hetero_CpG)),]
diff_not_hetero_order<-diff_not_hetero[rev(order(abs(diff_not_hetero$mean_db))),]

#' Examples which lose DNAm with passage
plt_hetero(diff_not_hetero_order[which(diff_not_hetero_order$mean_db>0),"CpG"][c(41,17,39,80)])
#' Examples which gain DNAm with passage
plt_hetero(diff_not_hetero_order[which(diff_not_hetero_order$mean_db<0),"CpG"][c(1,46,3,20)])


#' Combine these into a multi panel plot
tmp <- suppressWarnings(ggplot_gtable(ggplot_build(plt_hetero(diff_not_hetero_order[1:4,"CpG"]))))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]

ggsave(here("figs","Passage_representative_CpGs_FDR.pdf"),
       suppressWarnings(grid.arrange(plt_hetero(hetero_CpG[which(!(hetero_CpG%in%diff_CpG))][c(31,37,132,18)],"N","Y","Heteroskedastic CpGs"),
                    plt_hetero(diff_not_hetero_order[which(diff_not_hetero_order$mean_db>0),"CpG"][c(41,17,39,80)],"N","N","Hypomethylated CpGs"),
                    plt_hetero(diff_not_hetero_order[which(diff_not_hetero_order$mean_db<0),"CpG"][c(1,46,3,20)],"N","N","Hypermethylated CpGs"),legend, nrow=1, widths=c(4.3,3.85,3.85,1))),
       width = 9.75, height = 3.9)

ggsave(here("figs/jpeg","Passage_representative_CpGs_FDR.jpeg"),
       suppressWarnings(grid.arrange(plt_hetero(hetero_CpG[which(!(hetero_CpG%in%diff_CpG))][c(31,37,132,18)],"N","Y","Heteroskedastic CpGs"),
                                     plt_hetero(diff_not_hetero_order[which(diff_not_hetero_order$mean_db>0),"CpG"][c(41,17,39,80)],"N","N","Hypomethylated CpGs"),
                                     plt_hetero(diff_not_hetero_order[which(diff_not_hetero_order$mean_db<0),"CpG"][c(1,46,3,20)],"N","N","Hypermethylated CpGs"),legend, nrow=1, widths=c(4.3,3.85,3.85,1))),
       width = 9.75, height = 3.9)



#' ##### Export CpG background for python annotations

#' Using the cg ID to chromosome and coordinte annotation from illumina
#' https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html
anno_EPIC<-read.csv(here("data", "MethylationEPIC_v-1-0_B4.csv"), skip=7)

#' Will use the columns 'CHR' and 'MAPINFO' when possible as ther are GRCh37
background<-anno_EPIC[which(anno_EPIC$IlmnID%in%rownames(organoid_beta)), c('IlmnID', 'CHR', 'MAPINFO')]
write.csv(background, here("data", "passage_background.csv"), row.names = F)

#' The cDMR published coordinates are GRCh36 so will export those as well
background_36<-anno_EPIC[which(anno_EPIC$IlmnID%in%rownames(organoid_beta)), c('IlmnID', 'Chromosome_36', 'Coordinate_36')]
write.csv(background_36, here("data", "passage_background_build36.csv"), row.names = F)


### bar plot percent CpGs
CpG_proportions<-data.frame(CpG=c("Heteroskedastic","Hypomethylated","Hypermethylated","Not Passage\nAssociated"),
           Number=c(length(hetero_CpG),length(diff_CpG_db_hypo), length(diff_CpG_db_hyper), (nrow(pvals_long)-sum(length(hetero_CpG),length(diff_CpG_db_hypo), length(diff_CpG_db_hyper))) ))

CpG_proportions$Proportion_CpG<-(CpG_proportions$Number/nrow(pvals_long))*100

ggplot(CpG_proportions, aes(fill=CpG, y=Proportion_CpG, x=1)) + 
  geom_bar(position="stack", stat="identity")

CpG_proportions<-data.frame(CpG=c("Heteroskedastic","Hypomethylated","Hypermethylated"),
                            Number=c(length(hetero_CpG),length(diff_CpG_db_hypo), length(diff_CpG_db_hyper) ))

CpG_proportions$Proportion_CpG<-(CpG_proportions$Number/(sum(length(hetero_CpG),length(diff_CpG_db_hypo), length(diff_CpG_db_hyper))))*100

CpG_proportions$CpG<-factor(CpG_proportions$CpG, levels=rev(c("Heteroskedastic","Hypomethylated","Hypermethylated")))
CpG_proportions$cumsum<-cumsum(CpG_proportions$Proportion_CpG)

ggplot(CpG_proportions, aes(fill=CpG, y=Proportion_CpG, x=1)) + 
  geom_bar(position="stack", stat="identity", color="black")+coord_flip()+
  geom_text(aes(x = 1, y = cumsum-(Proportion_CpG/2), label = paste(round(Proportion_CpG),"%", sep="")),colour = "black")+
  th_present+theme_void()+scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

ggsave(here("figs","Proportions_Passage_all_CpGs.pdf"),width = 6, height = 0.5)
ggsave(here("figs/jpeg","Proportions_Passage_all_CpGs.jpeg"), width = 6, height = 0.5)



#' ##### Beta value distribution change with passage

# beta plot all CpGs
Beta_melted<- melt(organoid_beta)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid, by.x="ID", by.y="array.id")
Beta_Plot$passage.or.rescope.no_numeric.factor <- factor(Beta_Plot$passage.or.rescope.no_numeric, levels = c(16,14,11,10,8,7,6,4,3,2,1))

ggplot(Beta_Plot, aes(Beta,color=passage.or.rescope.no_numeric.factor))+
  geom_density(size=0.5)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=pass_col,name="Passage\nNumber")+th+theme(legend.text = element_text(size=7),
                                                                      legend.title = element_text(size=10),
                                                                      legend.key.size = unit(0.7,"line"))


ggsave(here("figs","Passage_all_CpGs.pdf"),width = 3.75, height = 2.5)
ggsave(here("figs/jpeg","Passage_all_CpGs.jpeg"), width = 3.75, height = 2.5)





#' ###### Overall Variance Across most Variable CpGs with Passage
Mval<-function(beta) log2(beta/(1-beta))
organoid_Mval = apply(organoid_beta, 1, Mval) # need mvalues for combat
organoid_Mval = as.data.frame(organoid_Mval)
organoid_Mval = t(organoid_Mval)

# variable mvalues only
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
ref_range_dnam<-sapply(1:nrow(organoid_Mval), function(x) Variation(organoid_Mval[x,]))
organoid_beta_VeryVariable<-organoid_beta[which(ref_range_dnam>=2.75),]#  51545

print(paste("There are ",nrow(organoid_beta_VeryVariable), " variable CpGs (10th-90th quantile range in M value >2.75)",sep=""))

# beta plot variable CpGs
Beta_melted<- melt(organoid_beta_VeryVariable)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid, by.x="ID", by.y="array.id")
Beta_Plot$passage.or.rescope.no_numeric.factor <- factor(Beta_Plot$passage.or.rescope.no_numeric, levels = c(16,14,11,10,8,7,6,4,3,2,1))

ggplot(Beta_Plot, aes(Beta,color=passage.or.rescope.no_numeric.factor))+
  geom_density(size=0.5)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=pass_col,name="Passage\nNumber")+th+theme(legend.text = element_text(size=7),
                                                                  legend.title = element_text(size=10),
                                                                  legend.key.size = unit(0.7,"line"))

ggsave(here("figs","Passage_variable_CpGs.pdf"),width = 3.75, height = 2.5)
ggsave(here("figs/jpeg","Passage_variable_CpGs.jpeg"), width = 3.75, height = 2.5)









#' ##### Beta distributions paired by individual
epic.organoid_paired<-epic.organoid[which(epic.organoid$sample_ID%in%epic.organoid$sample_ID[duplicated(epic.organoid$sample_ID)]),]
organoid_beta_VeryVariable_paird<-organoid_beta_VeryVariable[,which(epic.organoid$sample_ID%in%epic.organoid$sample_ID[duplicated(epic.organoid$sample_ID)])]

epic.organoid_paired<-do.call(rbind,lapply(1:length(unique(epic.organoid_paired$sample_ID)), function(x){
  sample<-unique(epic.organoid_paired$sample_ID)[x]
  samp<-epic.organoid_paired[epic.organoid_paired$sample_ID==sample,]
  samp<-samp[order(samp$passage.or.rescope.no_numeric),]
  samp$hilo<-1:nrow(samp)
  samp}))

epic.organoid_paired$hilo<-as.factor(epic.organoid_paired$hilo)
levels(epic.organoid_paired$hilo)<-c("low","high","highest")

Beta_melted<- melt(organoid_beta_VeryVariable_paird)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid_paired, by.x="ID", by.y="array.id")

labels<-as.data.frame(tapply(epic.organoid_paired$passage.or.rescope.no_numeric, epic.organoid_paired$sample_ID, function(x) paste(x, collapse=", ")))
colnames(labels)<-"passge"
labels$sample_ID<-rownames(labels)

Beta_Plot$sample_ID_long<-gsub("_", " ", Beta_Plot$sample_ID)
labels$sample_ID_long<-gsub("_", " ", labels$sample_ID)

Beta_Plot$sample_ID_long<-as.factor(Beta_Plot$sample_ID_long)
Beta_Plot$sample_ID_long<-factor(Beta_Plot$sample_ID_long, levels = c("259 SC","T116 SC", "298 SC","T091 SC","T095 TI", "270 TI","374 TI","268 TI","T046 SC","T046 TI","T049 SC","T049 TI" ))
labels$sample_ID_long<-factor(labels$sample_ID_long, levels = c("259 SC","T116 SC", "298 SC","T091 SC","T095 TI","270 TI","374 TI","268 TI","T046 SC","T046 TI","T049 SC","T049 TI" ))


ggplot()+
  geom_density(aes(Beta,color=hilo),Beta_Plot, size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values = c ("#41b6c4", "#225ea8", "#081d58"), name="Relative\nPassage\nLevel within\nPatient")+facet_wrap(~sample_ID_long, nrow=3)+
  geom_text(aes(0.75, 2.75, label=passge), data=labels, color="grey20")+th+
  theme(strip.text = element_text(size = 10),axis.text=element_text(size=4),panel.spacing = unit(0.7, "lines"))+th+
  scale_x_continuous(breaks = c(0,0.5,1))

ggsave(here("figs","Passage_paried_beta.pdf"),width = 8, height = 5.625)
ggsave(here("figs/jpeg","Passage_paried_beta.jpeg"), width = 8, height = 5.625)
    
 




#' ###### Are heteroskedastic CpGs contributing highly to the trimodial? -yup
organoid_beta_hetero<-organoid_beta[which(rownames(organoid_beta)%in%hetero_CpG),]

# density hi low
epic.organoid$passage_hilo<-"High Passage\n(>3 Passages)"
epic.organoid$passage_hilo[which(epic.organoid$passage.or.rescope.no_numeric<=3)]<-"Low Passage\n(<4 Passages)"


Beta_melted<- melt(organoid_beta_hetero)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid, by.x="ID", by.y="array.id")
Beta_Plot$passage.or.rescope.no_numeric.factor <- factor(Beta_Plot$passage.or.rescope.no_numeric, levels = c(16,14,11,10,8,7,6,4,3,2,1))

# all passage
ggplot(Beta_Plot, aes(Beta,color=passage.or.rescope.no_numeric.factor))+
  geom_density(size=1)+theme_bw()+xlab("DNAm Beta Value")+
  scale_color_manual(values=pass_col,name="Passage\nNumber")+th

# hilo
ggplot(Beta_Plot, aes(Beta,color=passage_hilo))+
  geom_density(size=1)+theme_bw()+xlab("DNAm Beta Value")+
  scale_color_manual(values=c ("#41b6c4", "#081d58"))+th



#' ###### Are differential CpGs contributing highly to the trimodial? -no
diff_CpG_db_hypo<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)>=0.15)]
diff_CpG_db_hyper<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)<=(-0.15))]

#' Hypomethylated
organoid_beta_diffhypo<-organoid_beta[which(rownames(organoid_beta)%in%diff_CpG_db_hypo),]
Beta_melted<- melt(organoid_beta_diffhypo)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid, by.x="ID", by.y="array.id")
Beta_Plot$passage.or.rescope.no_numeric.factor <- factor(Beta_Plot$passage.or.rescope.no_numeric, levels = c(16,14,11,10,8,7,6,4,3,2,1))

## all passage
ggplot(Beta_Plot, aes(Beta,color=passage.or.rescope.no_numeric.factor))+
  geom_density(size=1)+theme_bw()+xlab("DNAm Beta Value")+
  scale_color_manual(values=pass_col,name="Passage\nNumber")+th

## hilo
ggplot(Beta_Plot, aes(Beta,color=passage_hilo))+
  geom_density(size=1)+theme_bw()+xlab("DNAm Beta Value")+
  scale_color_manual(values=c ("#41b6c4", "#081d58"))+th

#' Hypermethylated
organoid_beta_diffhyper<-organoid_beta[which(rownames(organoid_beta)%in%diff_CpG_db_hyper),]
Beta_melted<- melt(organoid_beta_diffhyper)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid, by.x="ID", by.y="array.id")
Beta_Plot$passage.or.rescope.no_numeric.factor <- factor(Beta_Plot$passage.or.rescope.no_numeric, levels = c(16,14,11,10,8,7,6,4,3,2,1))

## all passage
ggplot(Beta_Plot, aes(Beta,color=passage.or.rescope.no_numeric.factor))+
  geom_density(size=1)+theme_bw()+xlab("DNAm Beta Value")+
  scale_color_manual(values=pass_col,name="Passage\nNumber")+th

## hilo
ggplot(Beta_Plot, aes(Beta,color=passage_hilo))+
  geom_density(size=1)+theme_bw()+xlab("DNAm Beta Value")+
  scale_color_manual(values=c ("#41b6c4", "#081d58"))+th







#' # Thresholding Trimodality
epic.organoid$thresholded_prior_ratio<-sapply(1:nrow(epic.organoid), function(x){
  print(x)
  converted<-as.numeric(round(organoid_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  passage_threshold_params(converted, counts, res)
})


ggplot(epic.organoid, aes(as.numeric(as.character(passage.or.rescope.no_numeric)), thresholded_prior_ratio))+
  geom_point(size=2,shape=21,color="black",aes(fill=as.factor(passage.or.rescope.no_numeric)))+xlab("Passage")+
  ylab("I/H Ratio")+theme_bw()+theme(axis.title = element_text(size=12))+
  #geom_text(aes(label=count, vjust=vjust, hjust=hjust), color="grey40", size=3)+
  scale_x_continuous(breaks=c(1,2,3,4,6,7,8,2,4,10,11,14,16))+ scale_fill_manual(values=pass_col,name="Passage\nNumber", guide=F)


ggsave(here("figs","Mixture_model_ratio_maximize.pdf"), width=3, height=2)


epic.organoid$thresholded_ratio_max<-F
epic.organoid$thresholded_ratio_max[which(epic.organoid$thresholded_prior_ratio>1)]<-T

percent_passing<-round((tapply(epic.organoid$thresholded_ratio_max, epic.organoid$passage.or.rescope.no_numeric, sum)/tapply(epic.organoid$array.id, epic.organoid$passage.or.rescope.no_numeric, length))*100,2)
passed_num<-tapply(epic.organoid$thresholded_ratio_max, epic.organoid$passage.or.rescope.no_numeric, sum)
org_numer<-tapply(epic.organoid$array.id, epic.organoid$passage.or.rescope.no_numeric, length)

df<-data.frame(passage=names(percent_passing), passing=percent_passing, pro_passing=percent_passing/100, count=org_numer, passed_num=passed_num)
df$vjust<-c(1.75,1.75,1.75,-0.75,1.75,-0.75,-0.75,0.5,1.75,1.75)
df$hjust<-c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,2,2,0.5)# for the text label you removed anyway

df$passage.factor <- factor(df$passage, levels = c(16,14,11,10,8,7,6,4,3,2,1))

df<-cbind(df,(binom.confint(df$passed_num, df$count, method="exact", conf.level=0.95)))
df$upper<-df$upper*100
df$lower<-df$lower*100

print(df)

ggplot(df, aes(as.numeric(as.character(passage)), passing))+
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="grey70", width=.25)+
  geom_line(color="grey20")+geom_point(size=1.25,shape=21,color="black",aes(fill=passage.factor))+xlab("Passage")+
  ylab("Samples with Trimodal\nDistribution (%)")+theme_bw()+theme(axis.title = element_text(size=10))+
  #geom_text(aes(label=count, vjust=vjust, hjust=hjust), color="grey40", size=3)+
  scale_x_continuous(breaks=c(1,2,3,4,6,7,8,11,14,16))+ scale_fill_manual(values=pass_col,name="Passage\nNumber", guide=F)

ggsave(here("figs","Mixture_model_ratio_threshold_maximize.pdf"), width=3, height=2)




plts_paired<-lapply(1:nrow(epic.organoid), function(x){#1:nrow(epic.organoid)
  print(x)
  passage<-paste("passage: ",epic.organoid$passage.or.rescope.no_numeric[x],"\nIndividual: ", epic.organoid$case.no[x],"\nRatio I/H: ",round(epic.organoid$thresholded_prior_ratio[x],2), sep="")
  converted<-as.numeric(round(organoid_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  draw_fit_params_gg(converted, counts, res,passage)
})

plts_paired_order<-plts_paired[order(epic.organoid$passage.or.rescope.no_numeric)]

pdf(here("figs","Original_organoids_thresholding_all_samples.pdf"))
plts_paired_order
dev.off()



#' ### paired samples for EM explination
plts_paired_example<-lapply(c(rev(grep("T046 SC", epic.organoid$sample_ID))), function(x){
  print(x)
  passage<-paste("passage: ",epic.organoid$passage.or.rescope.no_numeric[x],"\nIndividual: ", epic.organoid$case.no[x],"\nRatio I/H: ",round(epic.organoid$thresholded_prior_ratio[x],2), sep="")
  converted<-as.numeric(round(organoid_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  draw_fit_params_gg(converted, counts, res,passage)
})


ggsave(here("figs","Mix_model_Paired.pdf"),do.call(grid.arrange,plts_paired_example), width = 5, height=8)
ggsave(here("figs/jpeg","Mix_model_Paired.jpeg"),do.call(grid.arrange,plts_paired_example), width = 5, height=8)






#'## R Session Info
sessionInfo()
