#' ## Differential methylation with passage
#' Only in UT UD
validation_epic.organoid_UT_UD<-validation_epic.organoid[which(validation_epic.organoid$differentiation=="UD" & validation_epic.organoid$treatment=="UT"),]
table(validation_epic.organoid_UT_UD$passage_hilo)

validation_organoid_beta_UT_UD<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%validation_epic.organoid_UT_UD$array.id)]
identical(validation_epic.organoid_UT_UD$array.id, colnames(validation_organoid_beta_UT_UD))


mod<-model.matrix(~ 0 + passage + individual, data=validation_epic.organoid_UT_UD)
fit <- lmFit(validation_organoid_beta_UT_UD, mod)
ebfit <- eBayes(fit)

# covariate adjusted beta values
beta<-validation_organoid_beta_UT_UD

passage_db<-sapply(1:nrow(beta), function(x){
  sampleinfo_cpg<-validation_epic.organoid_UT_UD
  sampleinfo_cpg$beta<-as.numeric(beta[x,])
  
  fit<-lm(beta ~ passage + individual, data=sampleinfo_cpg)
  pval<-summary(fit)$coef["passage","Pr(>|t|)"]
  slope<-fit$coefficients[2]
  
  (min(validation_epic.organoid_UT_UD$passage)*slope) - (max(validation_epic.organoid_UT_UD$passage)*slope)})

passage_validation<-data.frame(p.value=ebfit$p.value[,"passage"], CpG=rownames(beta), db=passage_db)

# Adjust P values
passage_validation$p_adjusted<-p.adjust(passage_validation$p.value, method="BH")

diff_CpG_dbvalidation<-passage_validation[which(passage_validation$p_adjusted<0.05 & abs(passage_validation$db)>0.15),] #21098
diff_CpG_db_hypovalidation<-diff_CpG_dbvalidation$CpG[which((diff_CpG_dbvalidation$db)>=0.15)] #  20559
diff_CpG_db_hypervalidation<-diff_CpG_dbvalidation$CpG[which((diff_CpG_dbvalidation$db)<=(-0.15))] #  539


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


#' How many differential CpGs could overlap? 450K vs EPIC
print(paste("Of the ", nrow(diff_CpG_db), " CpGs differential with passage in the original organoids, ",
            length(diff_CpG_db$CpG[which(diff_CpG_db$CpG%in%rownames(validation_organoid_beta))]), " are in the validation data", sep=""))

diff_CpG_db_hypo_overlap<-diff_CpG_db_hypo[which(diff_CpG_db_hypo%in%rownames(validation_organoid_beta))]
diff_CpG_db_hyper_overlap<-diff_CpG_db_hyper[which(diff_CpG_db_hyper%in%rownames(validation_organoid_beta))]

diff_CpG_db_hypovalidation_overlap<-diff_CpG_db_hypovalidation[which(diff_CpG_db_hypovalidation%in%pvals_long$CpG)]
diff_CpG_db_hypervalidation_overlap<-diff_CpG_db_hypervalidation[which(diff_CpG_db_hypervalidation%in%pvals_long$CpG)]

print(paste("Of the ",length(diff_CpG_db_hypo_overlap)," hypo CpGs also on the 450K ",
            length(intersect(diff_CpG_db_hypovalidation_overlap, diff_CpG_db_hypo_overlap))," are also hypo in the validation cohort (",
            round((length(intersect(diff_CpG_db_hypovalidation_overlap, diff_CpG_db_hypo_overlap))/length(diff_CpG_db_hypo_overlap))*100,2),"%)",sep=""))

print(paste("Of the ",length(diff_CpG_db_hyper_overlap)," hypo CpGs also on the 450K ",
            length(intersect(diff_CpG_db_hypervalidation_overlap, diff_CpG_db_hyper_overlap))," are also hypo in the validation cohort (",
            round((length(intersect(diff_CpG_db_hypervalidation_overlap, diff_CpG_db_hyper_overlap))/length(diff_CpG_db_hyper_overlap))*100,2),"%)",sep=""))



#' ### CpG to gene associations
EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194

diff_genes_db_hypovalidation<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypovalidation)] ) #8689
diff_genes_db_hypervalidation<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypervalidation)] ) # 478

write.table(diff_genes_db_hypovalidation, file=here("data/validation/DNAm/","validation_genes_hypomethylation_UTUD.txt"), quote=F, row.names = F, col.names = F)
write.table(diff_genes_db_hypervalidation, file=here("data/validation/DNAm/","validation_genes_hypermethylation_UTUD.txt"), quote=F, row.names = F, col.names = F)

#'### Genes differential in original and validation
diff_genes_db_hypovalidation_original<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%intersect(diff_CpG_db_hypovalidation_overlap, diff_CpG_db_hypo_overlap))] ) # 4789
diff_genes_db_hypervalidation_original<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%intersect(diff_CpG_db_hypervalidation_overlap, diff_CpG_db_hyper_overlap))] ) # 410
write.table(diff_genes_db_hypovalidation_original, file=here("data/validation/DNAm/","validation_original_genes_hypomethylation_UTUD.txt"), quote=F, row.names = F, col.names = F)
write.table(diff_genes_db_hypervalidation_original, file=here("data/validation/DNAm/","validation_original_genes_hypermethylation_UTUD.txt"), quote=F, row.names = F, col.names = F)





#' ## Overall Variance Across most Variable CpGs with Passage
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
Mval<-function(beta) log2(beta/(1-beta))

organoid_mval = apply(validation_organoid_beta_UT_UD, 1, Mval)
organoid_mval = as.data.frame(organoid_mval)
organoid_mval = t(organoid_mval)

ref_range_dnam<-sapply(1:nrow(organoid_mval), function(x) Variation(organoid_mval[x,]))
validation_organoid_beta_VeryVariable<-validation_organoid_beta_UT_UD[which(ref_range_dnam>=2.75),]

print(paste("There are ",nrow(validation_organoid_beta_VeryVariable), " variable CpGs (10th-90th quantile range in M value >2.75)",sep=""))


# beta plot variable CpGs
Beta_melted<- melt(validation_organoid_beta_VeryVariable)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,validation_epic.organoid_UT_UD, by.x="ID", by.y="array.id")
Beta_Plot$passage.or.numeric.factor <- factor(Beta_Plot$passage, levels = c(12,9,8,4,3,2))

ggplot(Beta_Plot, aes(Beta,  color=passage.or.numeric.factor))+
  geom_density(size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=rev(c("#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4","#4CA5B1","#5E4FA2")), name="Passage\nNumber")
