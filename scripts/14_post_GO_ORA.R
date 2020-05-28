#'---
#'title: Prep lists fro Go ORA
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---
#'

#' #### Load Libraries
library(here)
options(stringsAsFactors = FALSE)

#' ### CpG to gene associations
EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194

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
print(paste("CpGs with significant (adjusted p<0.05) heteroskedactiy: ", length(hetero_CpG), sep=""))

diff_CpG<-rownames(organoid_beta)[which(pvals_long$diff_fdr<0.05)]
diff_CpG_db<-pvals_long[which(pvals_long$diff_fdr<0.05 & abs(pvals_long$mean_db)>0.15),]
print(paste("CpGs with significant (adjusted p<0.05; delta beta >0.05) differential methylation: ", nrow(diff_CpG_db), sep=""))

diff_CpG_db_hypo<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)>=0.15)]
diff_CpG_db_hyper<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)<=(-0.15))]

hetero_genes<-EPIC_genes[which(EPIC_genes$IlmnID%in%hetero_CpG),]
hetero_gene_name<-unique(hetero_genes$Gene.name)
hyper_genes<-EPIC_genes[which(EPIC_genes$IlmnID%in%diff_CpG_db_hyper),]
hyper_gene_name<-unique(hyper_genes$Gene.name)
hypo_genes<-EPIC_genes[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypo),]
hypo_gene_name<-unique(hypo_genes$Gene.name)



#'##  Tidy the output of erimineJ
#' ### hyper
hyperGO<-read.delim(here("data","hyper.erminej.txt"), skip=33, header=F)
colnames(hyperGO)<-c("","Name",	"ID",	"NumProbes"	,"NumGenes",	"RawScore",	"Pval",	"CorrectedPvalue"	,"MFPvalue",	"CorrectedMFPvalue"	,"Multifunctionality","SameAs","GeneMembers")
hyperGO[,1]<-NULL
hyperGO[,13]<-NULL

#' Limit the GeneMemberDiff to just genes in the GO gene set also associated to passage
hyperGO$GeneMemberDiff<-sapply(1:nrow(hyperGO), function(x){
  genemember<-strsplit(hyperGO$GeneMembers[x],"\\|")[[1]]
  paste(genemember[which(genemember%in%hyper_gene_name)],collapse = "|")
})

sighyper<-hyperGO[which(hyperGO$CorrectedPvalue<0.05),c(2,1,4,7,10,13)]
sighyper$CorrectedPvalue<-signif(sighyper$CorrectedPvalue, 1)
sighyper$GeneMemberDiff<-sapply(1:nrow(sighyper), function(x) gsub("\\|", ", ",sighyper$GeneMemberDiff[x]))

head(sighyper)

write.csv(sighyper, file=paste(here("output"), "/Supplementary_Table3_Hyper.csv", sep=""), row.names=F)





#' ### hypo
hypoGO<-read.delim(here("data","hypo.erminej.txt"), skip=33, header=F)
colnames(hypoGO)<-c("","Name",	"ID",	"NumProbes"	,"NumGenes",	"RawScore",	"Pval",	"CorrectedPvalue"	,"MFPvalue",	"CorrectedMFPvalue"	,"Multifunctionality","SameAs","GeneMembers")
hypoGO[,1]<-NULL
hypoGO[,13]<-NULL

#' Limit the GeneMemberDiff to just genes in the GO gene set also associated to passage
hypoGO$GeneMemberDiff<-sapply(1:nrow(hypoGO), function(x){
  genemember<-strsplit(hypoGO$GeneMembers[x],"\\|")[[1]]
  paste(genemember[which(genemember%in%hypo_gene_name)],collapse = "|")
})

sighypo<-hypoGO[which(hypoGO$CorrectedPvalue<0.05),c(2,1,4,7,10,13)]
sighypo$CorrectedPvalue<-signif(sighypo$CorrectedPvalue, 1)
sighypo$GeneMemberDiff<-sapply(1:nrow(sighypo), function(x) gsub("\\|", ", ",sighypo$GeneMemberDiff[x]))

head(sighypo)

write.csv(sighypo, file=paste(here("output"), "/Supplementary_Table2_Hypo.csv", sep=""), row.names=F)




#' ### hetero
heteroGO<-read.delim(here("data","hetero.erminej.txt"), skip=33, header=F)
colnames(heteroGO)<-c("","Name",	"ID",	"NumProbes"	,"NumGenes",	"RawScore",	"Pval",	"CorrectedPvalue"	,"MFPvalue",	"CorrectedMFPvalue"	,"Multifunctionality","SameAs","GeneMembers")
heteroGO[,1]<-NULL
heteroGO[,13]<-NULL

#' Limit the GeneMemberDiff to just genes in the GO gene set also associated to passage
heteroGO$GeneMemberDiff<-sapply(1:nrow(heteroGO), function(x){
  genemember<-strsplit(heteroGO$GeneMembers[x],"\\|")[[1]]
  paste(genemember[which(genemember%in%hetero_gene_name)],collapse = "|")
})

sighetero<-heteroGO[which(heteroGO$CorrectedPvalue<0.05),c(2,1,4,7,10,13)]
sighetero$CorrectedPvalue<-signif(sighetero$CorrectedPvalue, 1)
sighetero$GeneMemberDiff<-sapply(1:nrow(sighetero), function(x) gsub("\\|", ", ",sighetero$GeneMemberDiff[x]))

head(sighetero)

write.csv(sighetero, file=paste(here("output"), "/Supplementary_Table1_Hetero.csv", sep=""), row.names=F)


#' ## Overlap of enriched Go Groups
print(paste("Of the ",length(sighetero$ID)," and ",length(sighypo$ID)," significantly enriched GO gene sets (FDR < 0.05) in heteroskedastic and hypomethylated CpGs, respectively, ",length(intersect(sighetero$ID, sighypo$ID))," were overlapping",sep=""))
print(paste("Of the ",length(sighetero$ID)," and ",length(sighyper$ID)," significantly enriched GO gene sets (FDR < 0.05) in heteroskedastic and hypomethylated CpGs, respectively, ",length(intersect(sighetero$ID, sighyper$ID))," were overlapping",sep=""))
print(paste("Of the ",length(sighyper$ID)," and ",length(sighypo$ID)," significantly enriched GO gene sets (FDR < 0.05) in heteroskedastic and hypomethylated CpGs, respectively, ",length(intersect(sighyper$ID, sighypo$ID))," were overlapping",sep=""))


#'## R Session Info
sessionInfo()