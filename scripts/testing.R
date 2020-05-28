#'---
#'title: Exploration of CpGs effected by passage
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

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
#source(here("scripts","00_EM_array_uniform_background.R"))


#' #### Load Normalized Data
load(here("data","beta_organoids.RData"))




Mval<-function(beta) log2(beta/(1-beta))
organoid_Mval = apply(organoid_beta, 1, Mval) # need mvalues for combat
organoid_Mval = as.data.frame(organoid_Mval)
organoid_Mval = t(organoid_Mval)

# variable mvalues only
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
ref_range_dnam<-sapply(1:nrow(organoid_Mval), function(x) Variation(organoid_Mval[x,]))
organoid_beta_VeryVariable<-organoid_beta[which(ref_range_dnam>=2.75),]#  51545

print(paste("There are ",nrow(organoid_beta_VeryVariable), " variable CpGs (10th-90th quantile range in M value >2.75)",sep=""))


# epic.organoid$thresholded_prior_I<-sapply(1:nrow(epic.organoid), function(x){
#   print(x)
#   converted<-as.numeric(round(organoid_beta_VeryVariable[,x]*1000,0))
#   counts<-rep(1000, length(converted))
#   res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
#   res$prior_I
# })
# 
# 
# ggplot(epic.organoid, aes(as.numeric(as.character(passage.or.rescope.no_numeric)), thresholded_prior_I))+
#   geom_point(size=2,shape=21,color="black",aes(fill=as.factor(passage.or.rescope.no_numeric)))+xlab("Passage")+
#   ylab("Total Likelihood")+theme_bw()+theme(axis.title = element_text(size=12))+
#   #geom_text(aes(label=count, vjust=vjust, hjust=hjust), color="grey40", size=3)+
#   scale_x_continuous(breaks=c(1,2,3,4,6,7,8,2,4,10,11,14,16))+ scale_fill_manual(values=pass_col,name="Passage\nNumber", guide=F)
# 
# ggsave(here("figs","Mixture_model_prior_I.pdf"), width=3, height=2)
# 
# 
# epic.organoid$thresholded_prior10<-F
# epic.organoid$thresholded_prior10[which(epic.organoid$thresholded_prior_I>0.10)]<-T
# 
# percent_passing<-round((tapply(epic.organoid$thresholded_prior10, epic.organoid$passage.or.rescope.no_numeric, sum)/tapply(epic.organoid$array.id, epic.organoid$passage.or.rescope.no_numeric, length))*100,2)
# passed_num<-tapply(epic.organoid$thresholded_prior10, epic.organoid$passage.or.rescope.no_numeric, sum)
# org_numer<-tapply(epic.organoid$array.id, epic.organoid$passage.or.rescope.no_numeric, length)
# 
# df<-data.frame(passage=names(percent_passing), passing=percent_passing, pro_passing=percent_passing/100, count=org_numer, passed_num=passed_num)
# df$vjust<-c(1.75,1.75,1.75,-0.75,1.75,-0.75,-0.75,0.5,1.75,1.75)
# df$hjust<-c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,2,2,0.5)# for the text label you removed anyway
# 
# df$passage.factor <- factor(df$passage, levels = c(16,14,11,10,8,7,6,4,3,2,1))
# 
# df<-cbind(df,(binom.confint(df$passed_num, df$count, method="exact", conf.level=0.95)))
# df$upper<-df$upper*100
# df$lower<-df$lower*100
# 
# ggplot(df, aes(as.numeric(as.character(passage)), passing))+
#   geom_errorbar(aes(ymin=lower, ymax=upper), colour="grey70", width=.25)+
#   geom_line(color="grey20")+geom_point(size=1.25,shape=21,color="black",aes(fill=passage.factor))+xlab("Passage")+
#   ylab("Samples with Trimodal\nDistribution (%)")+theme_bw()+theme(axis.title = element_text(size=10))+
#   #geom_text(aes(label=count, vjust=vjust, hjust=hjust), color="grey40", size=3)+
#   scale_x_continuous(breaks=c(1,2,3,4,6,7,8,11,14,16))+ scale_fill_manual(values=pass_col,name="Passage\nNumber", guide=F)
# 
# ggsave(here("figs","Mixture_model_prior_I_threshold.pdf"), width=3, height=2)
# 
# 
# 
# plts_paired<-lapply(1:nrow(epic.organoid), function(x){#1:nrow(epic.organoid)
#   print(x)
#   passage<-paste("passage: ",epic.organoid$passage.or.rescope.no_numeric[x],"\nIndividual: ", epic.organoid$case.no[x],"\nPrior I: ",round(epic.organoid$thresholded_prior_I[x],2), sep="")
#   converted<-as.numeric(round(organoid_beta_VeryVariable[,x]*1000,0))
#   counts<-rep(1000, length(converted))
#   res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
#   draw_fit_params_gg(converted, counts, res,passage)
# })
# 
# plts_paired_order<-plts_paired[order(epic.organoid$passage.or.rescope.no_numeric)]
# 
# pdf(here("figs","Original_organoids_thresholding_all_samples.pdf"))
# plts_paired_order
# dev.off()



# epic.organoid$L_fits<-sapply(1:nrow(epic.organoid), function(x){
#   converted<-as.numeric(round(organoid_beta_VeryVariable[,x]*1000,0))
#   counts<-rep(1000, length(converted))
#   res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
#   fit_params_L(converted, counts, res)
# })
# 




source(here("scripts","00_EM_array_uniform_background_maximise_betabinom.R"))

epic.organoid$thresholded_prior_ratio<-sapply(1:nrow(epic.organoid), function(x){
  print(x)
  converted<-as.numeric(round(organoid_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .419, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  passage_threshold_params(converted, counts, res)
})


# epic.organoid$thresholded_prior_I_max<-sapply(1:nrow(epic.organoid), function(x){
#   print(x)
#   converted<-as.numeric(round(organoid_beta_VeryVariable[,x]*1000,0))
#   counts<-rep(1000, length(converted))
#   res = em(converted, counts, .419, .31, .27, 0.001, .1, .1, .90, .03, .5, .05)
#   res$prior_I
# })
# 
# 
ggplot(epic.organoid, aes(as.numeric(as.character(passage.or.rescope.no_numeric)), thresholded_prior_ratio))+
  geom_point(size=2,shape=21,color="black",aes(fill=as.factor(passage.or.rescope.no_numeric)))+xlab("Passage")+
  ylab("I/H Ratio")+theme_bw()+theme(axis.title = element_text(size=12))+
  #geom_text(aes(label=count, vjust=vjust, hjust=hjust), color="grey40", size=3)+
  scale_x_continuous(breaks=c(1,2,3,4,6,7,8,2,4,10,11,14,16))+ scale_fill_manual(values=pass_col,name="Passage\nNumber", guide=F)

ggsave(here("figs","Mixture_model_ratio_maximize.pdf"), width=3, height=2)

save(epic.organoid, file=paste(here("data"),"/tmp_epic.organoid_ratio.RData",sep=""))

load(here("data", "tmp_epic.organoid_ratio.RData"))

# 
# 
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

ggplot(df, aes(as.numeric(as.character(passage)), passing))+
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="grey70", width=.25)+
  geom_line(color="grey20")+geom_point(size=1.25,shape=21,color="black",aes(fill=passage.factor))+xlab("Passage")+
  ylab("Samples with Trimodal\nDistribution (%)")+theme_bw()+theme(axis.title = element_text(size=10))+
  #geom_text(aes(label=count, vjust=vjust, hjust=hjust), color="grey40", size=3)+
  scale_x_continuous(breaks=c(1,2,3,4,6,7,8,11,14,16))+ scale_fill_manual(values=pass_col,name="Passage\nNumber", guide=F)

ggsave(here("figs","Mixture_model_ratio_threshold_maximize.pdf"), width=3, height=2)


# 
plts_paired<-lapply(1:nrow(epic.organoid), function(x){#1:nrow(epic.organoid)
  print(x)
  passage<-paste("passage: ",epic.organoid$passage.or.rescope.no_numeric[x],"\nIndividual: ", epic.organoid$case.no[x],"\nRatio I/H: ",round(epic.organoid$thresholded_prior_ratio[x],2), sep="")
  converted<-as.numeric(round(organoid_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .419, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  draw_fit_params_gg(converted, counts, res,passage)
})
# 
plts_paired_order<-plts_paired[order(epic.organoid$passage.or.rescope.no_numeric)]

pdf(here("figs","Original_organoids_thresholding_all_samples_maximize_ratio.pdf"))
plts_paired_order
dev.off()
# 
# 
# epic.organoid$L_fits_max<-sapply(1:nrow(epic.organoid), function(x){
#   converted<-as.numeric(round(organoid_beta_VeryVariable[,x]*1000,0))
#   counts<-rep(1000, length(converted))
#   res = em(converted, counts, .419, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
#   fit_params_L(converted, counts, res)
# })
# 
# 
# epic.organoid$improvement<-epic.organoid$L_fits_max-epic.organoid$L_fits
# epic.organoid$MIP<-""
# epic.organoid$MIP[which(epic.organoid$improvement>1000)]<-as.character(epic.organoid$case.no[which(epic.organoid$improvement>1000)])
# 
# ggplot(epic.organoid, aes(x = reorder(sample_ID, passage.or.rescope.no_numeric))) +geom_segment(aes(xend = sample_ID, yend = L_fits, y = L_fits_max), color="grey")+
#   geom_point(aes(y=L_fits), fill="grey", shape=21, size=1)+geom_point(aes(y=L_fits_max, fill=as.factor(passage.or.rescope.no_numeric)), shape=21, size=2)+
#   theme_bw()+ geom_text(aes(y=L_fits_max+250, label=MIP))+
#   scale_fill_manual(values=pass_col,name="Passage\nNumber")+
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank())+xlab("Sample")+ylab("Total Likelihood")
# 
# ggsave(here("figs","Mixture_model_likelihood_compare_maximize.pdf"), width=12, height=5)

