#'---
#'title: Public colon cancer DNAm distributions
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---


#' #### Load Libraries
suppressMessages(library(GEOquery))
suppressMessages(library(GEOmetadb))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(lmtest))
suppressMessages(library(gridExtra))
suppressMessages(library(gtools))
suppressMessages(library(grid))
suppressMessages(library(here))

options(stringsAsFactors = FALSE)
options(scipen = 999)


#' #### Load Functions
source(here("scripts","00_pretty_plots.R"))
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
Mval<-function(beta) log2(beta/(1-beta))


#' ## Collect cancer colon
getSQLiteFile(destdir=here("data/public_cancer"))
con <- dbConnect(SQLite(), paste(here("data/public_cancer"),"/GEOmetadb.sqlite", sep=""))
meta<-dbGetQuery(con, "select title,description,series_id,gsm,source_name_ch1,characteristics_ch1 from gsm where gpl='GPL13534' OR gpl='GPL21145'")
colon_data<-meta[grep("colon", meta$characteristics_ch1),]


#' ## GSE42752
GSE42752<-meta[grep("GSE42752", meta$series_id),]
GSE42752_beta<- as.data.frame(exprs(getGEO("GSE42752")[[1]]))
identical(colnames(GSE42752_beta), GSE42752$gsm)
cancer_Mval = apply(GSE42752_beta, 1, Mval) # need mvalues for combat
cancer_Mval = as.data.frame(cancer_Mval)
cancer_Mval = t(cancer_Mval)
cancer_Mval<-as.data.frame(cancer_Mval)

ref_range_dnam_mval<-sapply(1:nrow(cancer_Mval), function(x) Variation(cancer_Mval[x,]))
dim(GSE42752_beta_VeryVariable<-GSE42752_beta[rev(order(ref_range_dnam_mval)),])
dim(GSE42752_beta_VeryVariable<-GSE42752_beta_VeryVariable[1:71384 ,])# same number as organoid varible

Beta_melted<- melt(GSE42752_beta_VeryVariable)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]

colnames(Beta_Plot)<-c("ID","Beta")
Beta_Plot<-merge(Beta_Plot,GSE42752, by.x="ID", by.y="gsm")

GSE42752_plt_df<-Beta_Plot



#' ## GSE48684
GSE48684<-meta[grep("GSE48684", meta$series_id),]
GSE48684_beta<- as.data.frame(exprs(getGEO("GSE48684")[[1]]))
identical(colnames(GSE48684_beta), GSE48684$gsm)

GSE48684$description_simple<-sapply(1:nrow(GSE48684), function(x) strsplit(GSE48684$description[x],";\t|-")[[1]][1])

GSE48684_beta<-GSE48684_beta
GSE48684_cancer<-GSE48684
identical(colnames(GSE48684_beta), GSE48684_cancer$gsm)

cancer_Mval = apply(GSE48684_beta, 1, Mval) # need mvalues for combat
cancer_Mval = as.data.frame(cancer_Mval)
cancer_Mval = t(cancer_Mval)
cancer_Mval<-as.data.frame(cancer_Mval)

ref_range_dnam_mval<-sapply(1:nrow(cancer_Mval), function(x) Variation(cancer_Mval[x,]))

dim(GSE48684_beta_VeryVariable<-GSE48684_beta[rev(order(ref_range_dnam_mval)),])
dim(GSE48684_beta_VeryVariable<-GSE48684_beta_VeryVariable[1:71384 ,])# same number as organoid varible

Beta_melted<- melt(GSE48684_beta_VeryVariable)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]

colnames(Beta_Plot)<-c("ID","Beta")
Beta_Plot<-merge(Beta_Plot,GSE48684_cancer, by.x="ID", by.y="gsm")

GSE48684_plt_df<-Beta_Plot





#' ## Plot beta distrbutions of variable CpGs

GSE48684_plt_df$description_simple<-as.factor(GSE48684_plt_df$description)
levels(GSE48684_plt_df$description_simple)<-c(rep("adenoma",25), "cancer","normal colon\nfrom CRC patient","normal colon\nfrom healthy control" )

GSE42752_plt_df$diagnosis<-as.factor(GSE42752_plt_df$source_name_ch1)
levels(GSE42752_plt_df$diagnosis)<-c("adenocarcinoma","normal colon\nfrom CRC patient","normal colon\nfrom healthy control")

head(GSE48684_plt_df)

GSE42752_plt<-ggplot()+
  geom_density(aes(Beta, color=as.character(diagnosis)), GSE42752_plt_df, size=1)+theme_bw()+xlab("DNAm Beta Value")+
  scale_color_manual(values=c("#b10026","#1d91c0","#225ea8"), name="Diagnosis")+guides(color=guide_legend(ncol=2))+
  theme(legend.position="bottom",
        plot.margin = margin(0.5, 0.5, 0.25, 0.5, "cm"),
        legend.text=element_text(size=4))+
  ggtitle("GSE42752")+th

GSE42752_plt

GSE42752_grob<-arrangeGrob(GSE42752_plt, top = textGrob("A", x = unit(0, "npc")
                                                        , y  = unit(0, "npc"), just=c("left","top"),
                                                        gp=gpar(col="black", fontsize=14)))

GSE48684_plt<-ggplot()+
  geom_density(aes(Beta, color=as.character(description_simple)), GSE48684_plt_df, size=1)+theme_bw()+
  scale_color_manual(values=c("#6a51a3","#b10026","#1d91c0","#225ea8"), name="Diagnosis")+
  theme_bw()+xlab("DNAm Beta Value")+ guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom",
        plot.margin = margin(0.5, 0.5, 0.25, 0.5, "cm"),
        legend.text=element_text(size=4))+
  ggtitle("GSE48684")+th

GSE48684_plt

GSE48684_grob<-arrangeGrob(GSE48684_plt, top = textGrob("B", x = unit(0, "npc")
                                                        , y  = unit(0, "npc"), just=c("left","top"),
                                                        gp=gpar(col="black", fontsize=14)))

ggsave(here("figs","cancer_public_same_variable_cpg_num.pdf"),
       grid.arrange(arrangeGrob(GSE42752_grob,GSE48684_grob, ncol=2)),
       width=10, height=5)


    

#'## R Session Info
sessionInfo()        