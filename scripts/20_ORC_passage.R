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
suppressMessages(library(cowplot))



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


### compared to eachother
plt_df_passage<-data.frame(distance=c(hypo_distance,hyper_distance,hetero_distance), list_cpg=rep(c("Hypomethylated", "Hypermethylated","Heteroskedastic"), times = c(length(hypo_distance),length(hyper_distance),length(hetero_distance))))
ggplot(plt_df_passage, aes(log(distance), color=list_cpg))+geom_density(size=0.75)+
  scale_color_manual(values=c("#8da0cb","#66c2a5","#fc8d62"), name="CpG Type")+
  theme_bw()+th+xlab("Distance to ORC Site (log)")

ggsave(here("figs","Compare_to_each_other_ORC_CpGs_passage.pdf"),width = 5, height = 2)
ggsave(here("figs/jpeg","Compare_to_each_other_ORC_CpGs_passage.jpeg"), width = 5, height = 2)


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

perm=1000
hypo_rnd<-dis_means(background_hypo, hypo_distance, perm)
(length(which(hypo_rnd>mean(hypo_distance)))+1)/(perm+1)
(length(which(hypo_rnd<mean(hypo_distance)))+1)/(perm+1)

hyper_rnd<-dis_means(background_hyper, hyper_distance, perm)
(length(which(hyper_rnd>mean(hyper_distance)))+1)/(perm+1)
(length(which(hyper_rnd<mean(hyper_distance)))+1)/(perm+1)

hetero_rnd<-dis_means(background, hetero_distance, perm)
(length(which(hetero_rnd>mean(hetero_distance)))+1)/(perm+1)
(length(which(hetero_rnd<mean(hetero_distance)))+1)/(perm+1)

plt_df_background<-data.frame(distance=c(hypo_rnd,hyper_rnd,hetero_rnd), list_cpg=rep(c("hypo", "hyper","hetero"), each = perm))
plt_lines<-data.frame(distance=c(mean(hypo_distance), mean(hyper_distance), mean(hetero_distance)), list_cpg=c("hypo", "hyper","hetero"))

plt_df_passage<-data.frame(distance=c(hypo_distance,hyper_distance,hetero_distance), list_cpg=rep(c("hypo", "hyper","hetero"), times = c(length(hypo_distance),length(hyper_distance),length(hetero_distance))))


mean(hypo_distance)
mean(hyper_distance)
mean(hetero_distance)


ggplot()+geom_density(aes(distance),plt_df_background, fill="lightgrey", color="grey40")+
  geom_vline(xintercept=mean(background_distance), color="black")+
  geom_vline(aes(xintercept=distance),plt_lines, color="red")+
  theme_bw()+th+facet_wrap(~list_cpg, ncol=1, scales="free_y")
ggsave(here("figs","ORC_CpGs_passage.pdf"),width = 4, height = 6)
ggsave(here("figs/jpeg","ORC_CpGs_passage.jpeg"), width = 4, height = 6)


ggplot()+#geom_density(aes(distance),plt_df_background, fill="lightgrey", color="grey40")+
  geom_density(aes(distance, color=list_cpg),plt_df_passage)+
  geom_vline(xintercept=mean(background_distance), color="black")+
  geom_vline(aes(xintercept=distance),plt_lines, color="red")+
  theme_bw()+th

ggplot()+#geom_density(aes(distance),plt_df_background, fill="lightgrey", color="grey40")+
  geom_density(aes(log(distance), color=list_cpg),plt_df_passage)+
  #geom_vline(xintercept=mean(background_distance), color="black")+
  geom_vline(aes(xintercept=log(distance), color=list_cpg),plt_lines)+
  theme_bw()+th

ggplot()+#geom_density(aes(distance),plt_df_background, fill="lightgrey", color="grey40")+
  geom_histogram(aes(log(distance), color=list_cpg),plt_df_passage)+
  geom_vline(xintercept=mean(background_distance), color="black")+
  geom_vline(aes(xintercept=distance),plt_lines, color="red")+
  theme_bw()+th




######
# plot single random versus actual
######
# save the background of permutations
plt_df_background_permutated<-plt_df_background
plt_lines_background<-as.data.frame(tapply(plt_df_background_permutated$distance, plt_df_background_permutated$list_cpg, mean))
plt_lines_background$list_cpg<-rownames(plt_lines_background)
colnames(plt_lines_background)<-c("mean_distance", "list_cpg")
plt_lines_background$list_cpg<-c("Heteroskedastic","Hypermethylated","Hypomethylated")

set.seed(1)
hypo_back<-sample(background_hypo$IlmnID, length(hypo_distance))
hypo_back_distance<-anno_EPIC_minimal$distance_to_ORC[which(anno_EPIC_minimal$IlmnID%in%hypo_back)]
plt_df_background_hypo<-data.frame(distance=c(hypo_back_distance,hypo_distance), 
                              rnd_og=rep(c("Background", "Passage CpGs"), each = length(hypo_distance)),
                              list_cpg=rep(c("Hypomethylated"), each = length(hypo_distance)*2))
set.seed(1)
hyper_back<-sample(background_hyper$IlmnID, length(hyper_distance))
hyper_back_distance<-anno_EPIC_minimal$distance_to_ORC[which(anno_EPIC_minimal$IlmnID%in%hyper_back)]
plt_df_background_hyper<-data.frame(distance=c(hyper_back_distance,hyper_distance), 
                                   rnd_og=rep(c("Background", "Passage CpGs"), each = length(hyper_distance)),
                                   list_cpg=rep(c("Hypermethylated"), each = length(hyper_distance)*2))
set.seed(1)
hetero_back<-sample(background$IlmnID, length(hetero_distance))
hetero_back_distance<-anno_EPIC_minimal$distance_to_ORC[which(anno_EPIC_minimal$IlmnID%in%hetero_back)]
plt_df_background_hetero<-data.frame(distance=c(hetero_back_distance,hetero_distance), 
                                    rnd_og=rep(c("Background", "Passage CpGs"), each = length(hetero_distance)),
                                    list_cpg=rep(c("Heteroskedastic"), each = length(hetero_distance)*2))

plt_df_background<-rbind(plt_df_background_hypo,plt_df_background_hyper, plt_df_background_hetero)
plt_df_background$rnd_og<-factor(plt_df_background$rnd_og, levels=c("Background","Passage CpGs"))
plt_lines<-data.frame(distance=c(mean(hypo_distance), mean(hyper_distance), mean(hetero_distance)), 
                      list_cpg=c("Hypomethylated", "Hypermethylated","Heteroskedastic"))

ggplot()+
  geom_density(aes((distance), color=rnd_og, fill=rnd_og),plt_df_background)+
  geom_vline(aes(xintercept=mean_distance),plt_lines_background, color="grey30")+
  geom_vline(aes(xintercept=distance),plt_lines, color="#3676e8")+
  theme_bw()+th+facet_wrap(~list_cpg, ncol=1)+
  scale_fill_manual(values=c("lightgrey",NA))+scale_color_manual(values=c("lightgrey","cornflowerblue"))+
  xlab("Distance to ORC Peak")+guides(fill=guide_legend(title="CpGs"),color=guide_legend(title="CpGs"))
ggsave(here("figs","ORC_CpGs_passage_distribution.pdf"),width = 8, height = 6)
ggsave(here("figs/jpeg","ORC_CpGs_passage_distribution.jpeg"), width = 8, height = 6)


quantile(plt_df_background$distance, seq(0,1,0.1))
ggplot()+
  geom_density(aes((distance), color=rnd_og, fill=rnd_og),plt_df_background)+
  geom_vline(aes(xintercept=mean_distance),plt_lines_background, color="grey30")+
  geom_vline(aes(xintercept=distance),plt_lines, color="#3676e8")+
  theme_bw()+th+facet_wrap(~list_cpg, ncol=1)+
  scale_fill_manual(values=c("lightgrey",NA))+scale_color_manual(values=c("lightgrey","cornflowerblue"))+
  xlim(0,300000)+
  xlab("Distance to ORC Peak")+guides(fill=guide_legend(title="CpGs"),color=guide_legend(title="CpGs"))
ggsave(here("figs","ORC_CpGs_passage_distributionzoom.pdf"),width = 8, height = 6)
ggsave(here("figs/jpeg","ORC_CpGs_passage_distributionzoom.jpeg"), width = 8, height = 6)


ggplot()+
  geom_density(aes(log(distance), color=rnd_og, fill=rnd_og),plt_df_background)+
  geom_vline(aes(xintercept=log(mean_distance)),plt_lines_background, color="grey30")+
  geom_vline(aes(xintercept=log(distance)),plt_lines, color="#2125A5", size=1)+
  theme_bw()+th+facet_wrap(~list_cpg, ncol=1)+
  scale_fill_manual(values=c("lightgrey",NA))+scale_color_manual(values=c("lightgrey","#2125A5"))+
  xlab("Distance to ORC Site (log)")+guides(fill=guide_legend(title="CpGs"),color=guide_legend(title="CpGs"))
ggsave(here("figs","ORC_CpGs_passage_distributionlog.pdf"),width = 8, height = 6)
ggsave(here("figs/jpeg","ORC_CpGs_passage_distributionlog.jpeg"), width = 8, height = 6)


plt_df_background$list_cpg<-factor(plt_df_background$list_cpg, levels=c())
ggplot()+
  geom_density(aes(log(distance), color=rnd_og, fill=rnd_og),plt_df_background)+
  geom_vline(aes(xintercept=log(mean_distance)),plt_lines_background, color="grey30", size=0.75)+
  geom_vline(aes(xintercept=log(distance)),plt_lines, color="#2125A5", size=0.75)+
  theme_bw()+th+facet_wrap(~list_cpg, ncol=3)+
  scale_fill_manual(values=c("lightgrey",NA))+scale_color_manual(values=c("lightgrey","#2125A5"))+
  xlab("Distance to ORC Site (log)")+guides(fill=guide_legend(title="CpGs"),color=guide_legend(title="CpGs"))
ggsave(here("figs","ORC_CpGs_passage_distributionlog_wide.pdf"),width = 12, height = 3)
ggsave(here("figs/jpeg","ORC_CpGs_passage_distributionlog_wide.jpeg"), width = 12, height = 3)



### Random distances median
dis_med<-function(background_CpGs, original_list, perm){
  sapply(1:perm, function(x) {
    set.seed(x)
    rnd_cpg<-sample(background_CpGs$IlmnID, length(original_list))
    rnd_distance<-anno_EPIC_minimal$distance_to_ORC[which(anno_EPIC_minimal$IlmnID%in%rnd_cpg)]
    median(rnd_distance)})}

perm=1000
hypo_rnd_median<-dis_means(background_hypo, hypo_distance, perm)
(length(which(hypo_rnd_median>median(hypo_distance)))+1)/(perm+1)
(length(which(hypo_rnd_median<median(hypo_distance)))+1)/(perm+1)

hyper_rnd_median<-dis_means(background_hyper, hyper_distance, perm)
(length(which(hyper_rnd_median>median(hyper_distance)))+1)/(perm+1)
(length(which(hyper_rnd_median<median(hyper_distance)))+1)/(perm+1)

hetero_rnd_median<-dis_means(background, hetero_distance, perm)
(length(which(hetero_rnd_median>median(hetero_distance)))+1)/(perm+1)
(length(which(hetero_rnd_median<median(hetero_distance)))+1)/(perm+1)

plt_df_background_median<-data.frame(distance=c(hypo_rnd_median,hyper_rnd_median,hetero_rnd_median), list_cpg=rep(c("hypo", "hyper","hetero"), each = perm))
plt_lines_median<-data.frame(distance=c(median(hypo_distance), median(hyper_distance), median(hetero_distance)), list_cpg=c("hypo", "hyper","hetero"))

plt_df_passage_median<-data.frame(distance=c(hypo_distance,hyper_distance,hetero_distance), list_cpg=rep(c("hypo1", "hyper1","hetero1"), times = c(length(hypo_distance),length(hyper_distance),length(hetero_distance))))


median(hypo_distance)
median(hyper_distance)
median(hetero_distance)



############
#'### continuous
############

pvals_long_distance<-merge(pvals_long, anno_EPIC_minimal, by.x="CpG", by.y="IlmnID")
pvals_long_distance$direction<-"gain"
pvals_long_distance$direction[which((pvals_long_distance$mean_db)>=0)]<-"lose"

ggplot(pvals_long_distance, aes(log(distance_to_ORC), mean_db))+geom_bin2d()
ggplot(pvals_long_distance, aes(log(distance_to_ORC), mean_db))+geom_hex()
ggplot(pvals_long_distance, aes(distance_to_ORC, mean_db))+geom_bin2d()


ggplot(pvals_long_distance[sample(1:nrow(pvals_long_distance), 1000),], aes(distance_to_ORC, mean_db)) +
  geom_density_2d(show.legend = FALSE) +
  coord_cartesian(expand = FALSE)

ggplot(pvals_long_distance, aes(log(distance_to_ORC), mean_db))+geom_smooth()

ggplot(pvals_long_distance, aes(distance_to_ORC, mean_db))+geom_smooth()+facet_wrap(~direction, scales="free_y", ncol=1)
ggplot(pvals_long_distance, aes(-log10(distance_to_ORC), mean_db))+geom_smooth()+facet_wrap(~direction, scales="free_y", ncol=1)


ggplot(pvals_long_distance, aes(distance_to_ORC, mean_db))+geom_bin2d()+geom_smooth()+facet_wrap(~direction, scales="free_y", ncol=1)
ggplot(pvals_long_distance, aes(distance_to_ORC, mean_db))+geom_bin2d()+geom_smooth()+facet_wrap(~direction, scales="free_y")+ylim()
ggplot(pvals_long_distance, aes(distance_to_ORC, mean_db))+geom_bin2d()+geom_smooth()+facet_wrap(~direction, scales="free_y")+ylim()

ggplot(pvals_long_distance, aes(-log10(distance_to_ORC), mean_db))+geom_bin2d(bins=100)




cor(pvals_long_distance[which(pvals_long_distance$direction>0),]$mean_db, pvals_long_distance[which(pvals_long_distance$direction>0),]$distance_to_ORC)
cor(pvals_long_distance[which(pvals_long_distance$direction<0),]$mean_db, pvals_long_distance[which(pvals_long_distance$direction<0),]$distance_to_ORC)

plot_grid(ggplot(pvals_long_distance, aes(distance_to_ORC, mean_db))+geom_smooth()+facet_wrap(~direction, scales="free_y", ncol=1),
          ggplot(pvals_long_distance, aes(mean_db))+geom_density(fill="lightgrey")+facet_wrap(~direction, scales="free_y", ncol=1)+coord_flip(),
          ggplot(pvals_long_distance, aes(distance_to_ORC))+geom_density(fill="lightgrey")+theme_bw()+th)

placeholder<-ggplot() + theme_void()
plot_grid(ggplot(pvals_long_distance, aes(distance_to_ORC))+geom_density(fill="lightgrey")+theme_bw()+th,placeholder,
          ggplot(pvals_long_distance[which(pvals_long_distance$direction=="gain"),], aes(distance_to_ORC, mean_db))+geom_smooth()+ylim(0,-0.5),
          ggplot(pvals_long_distance[which(pvals_long_distance$direction=="gain"),], aes(mean_db))+geom_density(fill="lightgrey")+coord_flip()+xlim(0,-0.5),
          ggplot(pvals_long_distance[which(pvals_long_distance$direction=="lose"),], aes(distance_to_ORC, mean_db))+geom_smooth()+ylim(0,0.5),
          ggplot(pvals_long_distance[which(pvals_long_distance$direction=="lose"),], aes(mean_db))+geom_density(fill="lightgrey")+coord_flip()+xlim(0,0.5), 
          ncol=2,align = "hv")


plot_grid(ggplot(pvals_long_distance, aes(distance_to_ORC))+geom_density(fill="lightgrey")+theme_bw()+th,
          ggplot(pvals_long_distance[which(pvals_long_distance$direction=="gain"),], aes(distance_to_ORC, mean_db))+geom_smooth(),
          ggplot(pvals_long_distance[which(pvals_long_distance$direction=="lose"),], aes(distance_to_ORC, mean_db))+geom_smooth(),
          ncol=1,align = "hv")


## volcano style plot
pvals_long_distance$Interesting_CpG3<-sapply(1:nrow(pvals_long_distance), function(x) if(pvals_long_distance$diff_fdr[x]<=0.05){
  if(abs(pvals_long_distance$mean_db[x])>0.15){
    if(pvals_long_distance$mean_db[x]>0.15){"Lose Methylation\n(Significant)"}else{"Gain Methylation\n (Significant)"}
  }else{if(pvals_long_distance$mean_db[x]>0){"Lose Methylation"}else{"Gain Methylation"}}}else{"Not Significantly Different"})


myColors <- c(muted("red", l=80, c=30),"red",muted("blue", l=70, c=40),"blue", "grey")

color_possibilities<-c("Lose Methylation",
                       "Lose Methylation\n(Significant)",
                       "Gain Methylation",
                       "Gain Methylation\n (Significant)",
                       "Not Significantly Different")

names(myColors) <- color_possibilities
colscale <- scale_color_manual(name = "Passage Association",
                               values = myColors, drop = FALSE)


#omg
ggplot(pvals_long_distance[1:20000,], aes(mean_db, -log10(distance_to_ORC), color=Interesting_CpG3))+
  geom_point(shape=19, size=1, alpha=0.1)+theme_bw()+
  colscale+
  geom_vline(xintercept=c(-0.15,0.15), color="grey60")+
  ylab("-log 10 distance_to_ORC")+xlab("Delta Beta")+
  theme(plot.margin=unit(c(1,1,1,2),"cm"))+ th+
  guides(color = guide_legend(override.aes = list(size = 4)))

ggplot(pvals_long_distance[1:200000,], aes(mean_db, (distance_to_ORC), color=Interesting_CpG3))+
  geom_point(shape=19, size=1, alpha=0.1)+theme_bw()+
  colscale+
  geom_vline(xintercept=c(-0.15,0.15), color="grey60")+
  ylab("distance_to_ORC")+xlab("Delta Beta")+
  theme(plot.margin=unit(c(1,1,1,2),"cm"))+ th+
  guides(color = guide_legend(override.aes = list(size = 4)))

hist(-log10(pvals_long_distance$distance_to_ORC))
hist((pvals_long_distance$mean_db))
