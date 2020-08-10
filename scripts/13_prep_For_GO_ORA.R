#'---
#'title: Prep lists For GO ORA
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---
#'

#' #### Load Libraries
suppressMessages(library(here))
suppressMessages(library(dplyr))
options(stringsAsFactors = FALSE)


#' #### Load Normalized Data
load(here("data","beta_organoids.RData"))

#' ### CpG to gene associations
EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194


#' Gene IDs to GO gene sets
# erminej gene to GO annotation file
# file from https://gemma.msl.ubc.ca/annots/Generic_human_ensemblIds_noParents.an.txt.gz

generic_human_go<-read.delim(here("data","Generic_human_ensemblIds_noParents.an.txt"), skip=7, sep="\t", header = T)
head(generic_human_go)
print(paste("Of the ",length(unique(EPIC_genes$Gene.stable.ID))," genes associated to EPIC CpG ",length(which(unique(EPIC_genes$Gene.stable.ID)%in%generic_human_go$ProbeName)), " have an annotated GO group",sep="" ))

#'  So 21,151 emsembl IDs dont have a GO label?  confirmed a few in the browser
head(EPIC_genes$Gene.stable.ID[which(!(EPIC_genes$Gene.stable.ID%in%generic_human_go$ProbeName))])

generic_human_go[which(generic_human_go$ProbeName=="ENSG00000149658"),]

#' seems correct for sampled genes. So will save the file of genes with Go groups for the ORA backgroud. 
generic_human_go_EPIC<-generic_human_go[which(generic_human_go$ProbeName%in%unique(EPIC_genes$Gene.stable.ID)),]
write.table(generic_human_go_EPIC, file=paste(here("data"),"/Generic_human_ensemblIds_noParents_EPIC_array.an.txt", sep=""), sep="\t", quote=F, row.names = F)



#' ### Import differential DNAm and heteroskedasticity statstics
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

diff_CpG_db_hypo<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)>=0.15)]
diff_CpG_db_hyper<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)<=(-0.15))]

hetero_genes<-EPIC_genes[which(EPIC_genes$IlmnID%in%hetero_CpG),]
hetero_genes<-unique(hetero_genes$Gene.stable.ID)
hyper_genes<-EPIC_genes[which(EPIC_genes$IlmnID%in%diff_CpG_db_hyper),]
hyper_genes<-unique(hyper_genes$Gene.stable.ID)
hypo_genes<-EPIC_genes[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypo),]
hypo_genes<-unique(hypo_genes$Gene.stable.ID)


#' Passage CpGs to gene lists for enrichment
hetero_genes<-hetero_genes[which(hetero_genes%in%generic_human_go_EPIC$ProbeName)]
head(hetero_genes)
write.csv(hetero_genes, file=paste(here("data"),"/hetero_genes.csv", sep=""), row.names = F, quote = F)
print(paste("There are ",length(hetero_genes), " genes associated with heteroskedastic CpGs", sep=""))

hyper_db<-pvals_long[which(pvals_long$CpG%in%diff_CpG_db_hyper),]
hyper_genes<-merge(hyper_db[,c(6,3,10)], EPIC_genes[,c(3,6,7,9)], by.x="CpG", by.y="IlmnID")
head(hyper_genes)
#' for genes with multiple differential CpGs a mean will be taken of their delta betas
hyper_genes_mean_db<-hyper_genes[,c(4,2,5)] %>%
  group_by(Gene.stable.ID) %>%
  dplyr::summarize(Mean = mean(mean_db, na.rm=TRUE))

hyper_genes_mean_db<-as.data.frame(hyper_genes_mean_db)
hyper_genes_mean_db<-hyper_genes_mean_db[which(hyper_genes_mean_db$Gene.stable.ID%in%generic_human_go_EPIC$ProbeName),]
write.csv(hyper_genes_mean_db, file=paste(here("data"),"/hyper_genes.csv", sep=""), row.names = F, quote = F)
print(paste("There are ",nrow(hyper_genes_mean_db), " genes associated with hypermethylated CpGs", sep=""))


hypo_db<-pvals_long[which(pvals_long$CpG%in%diff_CpG_db_hypo),]
hypo_genes<-merge(hypo_db[,c(6,3,10)], EPIC_genes[,c(3,6,7,9)], by.x="CpG", by.y="IlmnID")
hypo_genes_mean_db<-hypo_genes[,c(4,2,5)] %>%
  group_by(Gene.stable.ID) %>%
  dplyr::summarize(Mean = mean(mean_db, na.rm=TRUE))

hypo_genes_mean_db<-as.data.frame(hypo_genes_mean_db)
hypo_genes_mean_db<-hypo_genes_mean_db[which(hypo_genes_mean_db$Gene.stable.ID%in%generic_human_go_EPIC$ProbeName),]
write.csv(hypo_genes_mean_db, file=paste(here("data"),"/hypo_genes.csv", sep=""), row.names = F, quote = F)
print(paste("There are ",nrow(hypo_genes_mean_db), " genes associated with hypomethylated CpGs", sep=""))

#' ## GO ORA in erminej run as below:
#' Gene annotation file: Generic_human_ensemblIds_noParents_EPIC_array.an.txt
#' Analysis -> Run Analysis -> ORA
#' ORA Paste genes lists as "Quick list" (did not use a score)
#' Gene lists to paste: hetero_genes.csv hyper_genes.csv hypo_genes.csv
#' Biological process only
#' max set size=400 min=5
#' no log of score, larger scores are not better gene score threshold 0.0
#' analysis-> save analysis -> chose each ORA run seperately -> include all genes in output (yes)
#' save hyper.erminej.txt hypo.erminej.txt hetero.erminej.txt

#'## R Session Info
sessionInfo()