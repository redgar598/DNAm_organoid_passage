#'---
#'title: CpGs hyper DMR, in wnt signalling
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


#' #### Load Normalized Data
load(here("data","beta_organoids.RData"))


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
diff_CpG<-rownames(organoid_beta)[which(pvals_long$diff_fdr<0.05)]
diff_CpG_db<-pvals_long[which(pvals_long$diff_fdr<0.05 & abs(pvals_long$mean_db)>0.15),]
diff_CpG_db_hypo<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)>=0.15)]
diff_CpG_db_hyper<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)<=(-0.15))]


#' Using the cg ID to chromosome and coordinte annotation from illumina
#' https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html
anno_EPIC<-read.csv(here("data", "MethylationEPIC_v-1-0_B4.csv"), skip=7)

#' ## c-DMR enrichment of passage CpGs
#' c-DMRs are from # https://www.nature.com/articles/ng.298#Sec25
#' supplementary data 2
#' The paper does not list the genome build but they DNAm was measured on the Nimblegen CHARM array. The GEO page for the Nimblegen CHARM array used does say is it hg18 (GRch36)

cDMR<-read.csv(here("data","EPIC_Features_cDMR36.csv"))

#' remove []
cDMR$X1<-gsub("\\[|\\]","", cDMR$X1)
print(paste("There are ",length(irizarry_sup_cDMR<-cDMR$X0[which(cDMR$X1!="")])," CpGs overlapping a cDMR",sep="") )

anno_EPIC<-anno_EPIC[which(anno_EPIC$IlmnID%in%cDMR$X0),]

#' Illumina annotation includes a CDMR column, so will compare to see my annotation is close
#' But they don't have direction, just yes no
illumina_cDMR<-anno_EPIC$IlmnID[which(anno_EPIC$DMR=="CDMR")] # 6581
print(paste("There are ",length(illumina_cDMR)," CpGs overlapping a cDMR",sep="") )
print(paste("Of the ",length(illumina_cDMR)," CpGs overlapping a cDMR, my annotation captures ",length(intersect(irizarry_sup_cDMR, illumina_cDMR))," of them.",sep="") )

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



#' #### Monte carlo function on background of CpGs
cDMR_DMR<-cDMR[which(cDMR$cDMR!="Not cDMR"),]
cDMR_DMR$cDMR<-as.factor(cDMR_DMR$cDMR)
cDMR_DMR<-cDMR_DMR[which(cDMR_DMR$cDMR!="cDMR with more DNAm in cancer"),]


cDMR_hyper_hits<-cDMR_DMR[which(cDMR_DMR$CpG%in%diff_CpG_db_hyper),] 


#' ### CpG to gene associations
EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194
EPIC_genes[which(EPIC_genes$IlmnID%in%cDMR_hyper_hits$CpG),]

#AXIN2: cg09231741


cDMR_coor<-read.csv(here("data","41588_2009_BFng298_MOESM18_ESM.csv"), skip=3)

load(here("data","ensembl_transcripts.RData"))


#' #' ## DMR sytle plot of a gene region
#'   CpG_OI<-EPIC_genes[which(EPIC_genes$Gene.name==gene_name),]$IlmnID
#' 
#'   
#'   hypo_CpG_gene<-diff_CpG_db_hypo[which(diff_CpG_db_hypo%in%CpG_OI)]
#'   hyper_CpG_gene<-diff_CpG_db_hyper[which(diff_CpG_db_hyper%in%CpG_OI)]
#' 
#'   hypo_CpG_gene<-anno_EPIC[which(anno_EPIC$IlmnID%in%hypo_CpG_gene),c("IlmnID","MAPINFO")]
#'   hyper_CpG_gene<-anno_EPIC[which(anno_EPIC$IlmnID%in%hyper_CpG_gene),c("IlmnID","MAPINFO")]
#' 
#'   if(nrow(hypo_CpG_gene)>0) {hypo_CpG_gene$direction<-"Hypo"} else{}
#'   hyper_CpG_gene$direction<-"Hyper"
#' 
#'   sig<-do.call(rbind,list(hypo_CpG_gene,hyper_CpG_gene))
#'   
#' gene<-ensembl_transcripts[which(ensembl_transcripts$Gene.name=="AXIN2"),][1,]
#' 
#' hyper_cdmr<-cDMR_coor[which(cDMR_coor$delta.M==cDMR_hyper_hits[which(cDMR_hyper_hits$CpG=="cg09231741"),"cancer_db"]),]
#' 
#' 
#'   betas_org<-as.data.frame(organoid_beta[which(rownames(organoid_beta)%in%CpG_OI),])
#'   betas_org$CpG<-rownames(betas_org)
#'   betas_org<-melt(betas_org)
#'   betas_org$sample.type<-"organoid"
#'   
#'   # betas_pri<-as.data.frame(ibd_combo_og[which(rownames(ibd_combo_og)%in%CpG_OI),])
#'   # betas_pri$CpG<-rownames(betas_pri)
#'   # betas_pri<-melt(betas_pri)
#'   # betas_pri$sample.type<-"primary"
#'   betas<-betas_org
#'   #betas<-rbind(betas_org,betas_pri)
#'   
#'   betas<-merge(betas, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")
#'   plt<-merge(betas, epic.organoid, by.x="variable", by.y="array.id")
#'   
#'   sum_stat<-melt(tapply(plt$value, list(plt$MAPINFO,plt$passage.or.rescope.no_numeric), mean))
#'   colnames(sum_stat)<-c("MAPINFO","passage.or.rescope.no_numeric","value")
#'   sum_stat<-sum_stat[which(!(is.na(sum_stat$value))),]
#'   
#'   
#'   
#'   ggplot()+ geom_rect(aes(xmin = Gene.start..bp., xmax = Gene.end..bp., ymin = -0.05, ymax = 0), fill = '#f2f500', alpha = 1, data = gene)+
#'     geom_rect(aes(xmin = promoter_start, xmax = promoter_end, ymin = -0.05, ymax = 0), fill = 'blue', alpha = 1, data = gene)+
#'     geom_rect(aes(xmin = start, xmax = end, ymin = -0.10, ymax = -0.05), fill = 'red', alpha = 1, data = hyper_cdmr)+
#'     geom_vline(data = sig, aes(xintercept = MAPINFO), size=1)+
#'     geom_point(aes(MAPINFO, value, fill=as.factor(passage.or.rescope.no_numeric)), plt, size=1.5,shape=21, color="white")+
#'     geom_line(aes(MAPINFO,value, color=as.factor(passage.or.rescope.no_numeric)),sum_stat,size=0.5)+
#'     theme_bw()+th+ylim(-0.1,1)+ylab("DNAm Beta Value")+scale_fill_manual(values=pass_col, name="passage")+scale_color_manual(values=pass_col, name="passage")
#'   
  


### alternate
gene_name<-"AXIN2"


CpG_OI<-EPIC_genes[which(EPIC_genes$Gene.name==gene_name),]
hyper_CpG_gene<-diff_CpG_db_hyper[which(diff_CpG_db_hyper%in%CpG_OI$IlmnID)]
CpG_OI_sig<-hyper_CpG_gene
CpG_in_sig_range<-CpG_OI[which(CpG_OI$MAPINFO>=min(EPIC_genes[which(EPIC_genes$IlmnID%in%CpG_OI_sig),]$MAPINFO) & CpG_OI$MAPINFO<=max(EPIC_genes[which(EPIC_genes$IlmnID%in%CpG_OI_sig),]$MAPINFO)),]
CpG_OI<-CpG_in_sig_range$IlmnID


betas_org<-as.data.frame(organoid_beta[which(rownames(organoid_beta)%in%CpG_OI),])
betas_org$CpG<-rownames(betas_org)
betas_org<-melt(betas_org)
betas_org$sample.type<-"organoid"

# betas_pri<-as.data.frame(ibd_combo_og[which(rownames(ibd_combo_og)%in%CpG_OI),])
# betas_pri$CpG<-rownames(betas_pri)
# betas_pri<-melt(betas_pri)
# betas_pri$sample.type<-"primary"
betas<-betas_org
#betas<-rbind(betas_org,betas_pri)

betas<-merge(betas, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")
plt<-merge(betas, epic.organoid, by.x="variable", by.y="array.id")

sum_stat<-melt(tapply(plt$value, list(plt$MAPINFO,plt$passage.or.rescope.no_numeric), mean))
colnames(sum_stat)<-c("MAPINFO","passage.or.rescope.no_numeric","value")
sum_stat<-sum_stat[which(!(is.na(sum_stat$value))),]

ggplot()+ 
  geom_point(aes(MAPINFO, value, fill=as.factor(passage.or.rescope.no_numeric)), plt, size=2,shape=21, color="white")+
  geom_line(aes(MAPINFO,value, color=as.factor(passage.or.rescope.no_numeric)),sum_stat,size=0.5)+
  theme_bw()+th+ylim(-0.1,1)+ylab("DNAm Beta Value")+scale_fill_manual(values=pass_col, name="passage")+scale_color_manual(values=pass_col, name="passage")


## load cancer primary
suppressMessages(library(GEOquery))
suppressMessages(library(GEOmetadb))

#' ## Collect cancer colon
getSQLiteFile(destdir=here("data/public_cancer"))
con <- dbConnect(SQLite(), paste(here("data/public_cancer"),"/GEOmetadb.sqlite", sep=""))
meta<-dbGetQuery(con, "select title,description,series_id,gsm,source_name_ch1,characteristics_ch1 from gsm where gpl='GPL13534' OR gpl='GPL21145'")
colon_data<-meta[grep("colon", meta$characteristics_ch1),]


#' ## GSE42752
GSE42752<-meta[grep("GSE42752", meta$series_id),]
GSE42752_beta<- as.data.frame(exprs(getGEO("GSE42752")[[1]]))
identical(colnames(GSE42752_beta), GSE42752$gsm)

#' ## GSE48684
GSE48684<-meta[grep("GSE48684", meta$series_id),]
GSE48684_beta<- as.data.frame(exprs(getGEO("GSE48684")[[1]]))
identical(colnames(GSE48684_beta), GSE48684$gsm)

GSE48684$description_simple<-sapply(1:nrow(GSE48684), function(x) strsplit(GSE48684$description[x],";\t|-")[[1]][1])

GSE48684_beta<-GSE48684_beta
GSE48684_cancer<-GSE48684
identical(colnames(GSE48684_beta), GSE48684_cancer$gsm)

GSE48684_cancer$description_simple<-as.factor(GSE48684_cancer$description)
levels(GSE48684_cancer$description_simple)<-c(rep("adenoma",25), "cancer","normal colon\nfrom CRC patient","normal colon\nfrom healthy control" )

GSE42752$diagnosis<-as.factor(GSE42752$source_name_ch1)
levels(GSE42752$diagnosis)<-c("adenocarcinoma","normal colon\nfrom CRC patient","normal colon\nfrom healthy control")
