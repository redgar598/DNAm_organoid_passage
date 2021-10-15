#'---
#'title: origin of replication
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ### Load Libraries
suppressMessages(library(minfi))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(RColorBrewer))
suppressMessages(library(here))
suppressMessages(library(binom))
suppressMessages(library(limma))
suppressMessages(library(sva))
suppressMessages(library(pamr))
suppressMessages(library(GEOquery))
suppressMessages(library(GEOmetadb))

suppressMessages(library(dplyr))
suppressMessages(library(lmtest))
suppressMessages(library(gridExtra))
suppressMessages(library(gtools))
suppressMessages(library(rafalib))


options(stringsAsFactors = FALSE)


#' ### Load Functions
source(here("scripts","00_pretty_plots.R"))
suppressMessages(source(here("scripts","00_heat_scree_plot_generic.R")))
source(here("scripts","00_EM_array_uniform_background_maximise_betabinom.R"))



#' ### load original organoid passage CpGs
load(here("data","beta_organoids.RData"))
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

diff_CpG_db_hypo<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)>=0.15)] #  11772
diff_CpG_db_hyper<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)<=(-0.15))] #  5214


#'## hg19 ORC
ORC_coor<-read.csv(here("data/Miotto_ORC_hg19_suppT1.csv"),skip=2, header=F)
colnames(ORC_coor)<-c("CHR","start","end")
ORC_coor$CHR<-gsub("chr","",ORC_coor$CHR)

#'# distance to closest ORC
#'###https://www.pnas.org/content/pnas/113/33/E4810.full.pdf
#'## CpG Coordiantes from illumina annotation
anno_EPIC<-read.csv(here("data", "MethylationEPIC_v-1-0_B4.csv"), skip=7)
anno_EPIC_minimal<-anno_EPIC[,c("IlmnID","CHR", "MAPINFO")]
anno_EPIC_minimal<-anno_EPIC_minimal[which(!(is.na(anno_EPIC_minimal$MAPINFO))),]


x=1

anno_EPIC_minimal$distance_to_ORC<-sapply(1:nrow(anno_EPIC_minimal), function(x){
  chr_ORC<- ORC_coor[which(ORC_coor$CHR==anno_EPIC_minimal$CHR[x]),]
  start_dis<-abs(chr_ORC$start-anno_EPIC_minimal$MAPINFO[x])
  end_dis<-abs(chr_ORC$end-anno_EPIC_minimal$MAPINFO[x])
  min(start_dis, end_dis)})


hypo_distance<-anno_EPIC_minimal$distance_to_ORC[which(anno_EPIC_minimal$IlmnID%in%diff_CpG_db_hypo)]
hyper_distance<-anno_EPIC_minimal$distance_to_ORC[which(anno_EPIC_minimal$IlmnID%in%diff_CpG_db_hyper)]
hetero_distance<-anno_EPIC_minimal$distance_to_ORC[which(anno_EPIC_minimal$IlmnID%in%hetero_CpG)]

mean(hypo_distance)
mean(hyper_distance)
mean(hetero_distance)
mean(anno_EPIC_minimal$distance_to_ORC)


#' ### Backgorund CpG List to compare to
#' CpGs which have an passage=0 intercept <0.15 need to be excluded from the hypomethylayed background. Since they start <0.15 they can decrease by >0.15 as there is a lower limit of 0.
#' CpGs which have an passage=0 intercept >0.85 need to be excluded from the hypermethylayed background. Since they start >0.15 they can increase by >0.15 as there is a upper limit of 1.
background<-anno_EPIC[which(anno_EPIC$IlmnID%in%rownames(organoid_beta)), c('IlmnID', 'CHR', 'MAPINFO')]
background_distance<-anno_EPIC_minimal$distance_to_ORC[which(anno_EPIC_minimal$IlmnID%in%background$IlmnID)]
mean(background_distance)

intercepts<-sapply(1:nrow(organoid_beta), function(x) as.numeric(lm(organoid_beta[x,]~epic.organoid$passage.or.rescope.no_numeric)$coefficients[1]))

rm_hypo<-rownames(organoid_beta)[which(intercepts<=0.15)]
rm_hyper<-rownames(organoid_beta)[which(intercepts>=0.85)]

length(rm_hypo)
length(rm_hyper)

background_hypo<-background[which(!(background$IlmnID%in%rm_hypo)),]
background_hyper<-background[which(!(background$IlmnID%in%rm_hyper)),]

### Random distances
dis_means<-function(background_CpGs, original_list, perm){
  sapply(1:perm, function(x) {
    set.seed(x)
    rnd_cpg<-sample(background_CpGs$IlmnID, length(original_list))
    rnd_distance<-anno_EPIC_minimal$distance_to_ORC[which(anno_EPIC_minimal$IlmnID%in%rnd_cpg)]
    mean(rnd_distance)})}

perm=10000
hypo_rnd<-dis_means(background_hypo, hypo_distance, perm)
(length(which(hypo_rnd>mean(hypo_distance)))+1)/(perm+1)
(length(which(hypo_rnd<mean(hypo_distance)))+1)/(perm+1)

hyper_rnd<-dis_means(background_hyper, hyper_distance, perm)
(length(which(hyper_rnd>mean(hyper_distance)))+1)/(perm+1)
(length(which(hyper_rnd<mean(hyper_distance)))+1)/(perm+1)

hetero_rnd<-dis_means(background, hetero_distance, perm)
(length(which(hetero_rnd>mean(hetero_distance)))+1)/(perm+1)
(length(which(hetero_rnd<mean(hetero_distance)))+1)/(perm+1)

plt_df<-data.frame(distance=c(hypo_rnd,hyper_rnd,hetero_rnd), list_cpg=rep(c("hypo", "hyper","hetero"), each = perm))
plt_lines<-data.frame(distance=c(mean(hypo_distance), mean(hyper_distance), mean(hetero_distance)), list_cpg=c("hypo", "hyper","hetero"))

ggplot()+geom_density(aes(distance),plt_df, fill="lightgrey", color="grey40")+
  geom_vline(xintercept=mean(background_distance), color="black")+
  geom_vline(aes(xintercept=distance),plt_lines, color="red")+
  theme_bw()+th+facet_wrap(~list_cpg, ncol=1, scales="free_y")
ggsave(here("figs","ORC_CpGs_passage.pdf"),width = 4, height = 6)
ggsave(here("figs/jpeg","ORC_CpGs_passage.jpeg"), width = 4, height = 6)
