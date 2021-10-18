###########
#' Table of genes
###########
library(dplyr)

genes<-c(diff_genes_db_hypovalidation_original$V1[1:150],"HHAT")
stats_table<-sleuth_significant_low


load(here("data/validation/CpG_validated.RData"))
EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194


stats<-stats_table[which(stats_table$ext_gene%in%genes),c("ext_gene","target_id","pval","qval")]

EPIC_genes_CpGs<-EPIC_genes[which(EPIC_genes$Gene.name%in%stats$ext_gene),]
EPIC_genes_CpGs_hypo<-EPIC_genes_CpGs[which(EPIC_genes_CpGs$IlmnID%in%CpG_hypo_validated),]
EPIC_genes_CpGs_hypo<-as.data.frame(EPIC_genes_CpGs_hypo %>% 
  group_by(Gene.name) %>% 
  mutate(hypo_DNAm_CpG = paste0(unique(IlmnID), collapse = ", ")) )
EPIC_genes_CpGs_hypo<-EPIC_genes_CpGs_hypo[,c("Gene.name","hypo_DNAm_CpG")]
EPIC_genes_CpGs_hypo<-EPIC_genes_CpGs_hypo[!duplicated(EPIC_genes_CpGs_hypo),]
EPIC_genes_CpGs_hypo

EPIC_genes_CpGs_hyper<-EPIC_genes_CpGs[which(EPIC_genes_CpGs$IlmnID%in%CpG_hyper_validated),]
EPIC_genes_CpGs_hyper<-as.data.frame(EPIC_genes_CpGs_hyper %>% 
                                      group_by(Gene.name) %>% 
                                      mutate(hyper_DNAm_CpG = paste0(unique(IlmnID), collapse = ", ")) )
EPIC_genes_CpGs_hyper<-EPIC_genes_CpGs_hyper[,c("Gene.name","hyper_DNAm_CpG")]
EPIC_genes_CpGs_hyper<-EPIC_genes_CpGs_hyper[!duplicated(EPIC_genes_CpGs_hyper),]
EPIC_genes_CpGs_hyper

EPIC_genes_CpGs<-merge(EPIC_genes_CpGs_hyper, EPIC_genes_CpGs_hypo, by="Gene.name", all=T)

stats_cpgs<-merge(stats, EPIC_genes_CpGs, by.x="ext_gene", by.y="Gene.name")
