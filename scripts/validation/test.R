validation_epic.organoid_UT_UD<-validation_epic.organoid[which(validation_epic.organoid$differentiation=="UD" & validation_epic.organoid$treatment=="UT"),]
table(validation_epic.organoid_UT_UD$passage_hilo)

validation_organoid_beta_UT_UD<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%validation_epic.organoid_UT_UD$array.id)]
identical(validation_epic.organoid_UT_UD$array.id, colnames(validation_organoid_beta_UT_UD))


#' ## Differential methylation with passage (high low)
mod<-model.matrix(~ passage_hilo, data=validation_epic.organoid_UT_UD)
fit <- lmFit(validation_organoid_beta_UT_UD, mod)
ebfit <- eBayes(fit)

# covariate adjusted beta values
beta<-validation_organoid_beta_UT_UD

passage_hilo_db<-sapply(1:nrow(beta), function(x){
  sampleinfo_cpg<-validation_epic.organoid_UT_UD
  sampleinfo_cpg$beta<-as.numeric(beta[x,])
  mns<-tapply(sampleinfo_cpg$beta,  sampleinfo_cpg$passage_hilo, mean)
  mns[2]-mns[1]# low minus high
  })

passage_validation_hilo<-data.frame(p.value=ebfit$p.value[,"passage_hilolow"], CpG=rownames(beta), db=passage_hilo_db)

# Adjust P values
passage_validation_hilo$p_adjusted<-p.adjust(passage_validation_hilo$p.value, method="BH")

diff_CpG_dbvalidation<-passage_validation_hilo[which(passage_validation_hilo$p_adjusted<0.05 & abs(passage_validation_hilo$db)>0.15),] #5532
diff_CpG_db_hypovalidation<-diff_CpG_dbvalidation$CpG[which((diff_CpG_dbvalidation$db)>=0.15)] #  5458
diff_CpG_db_hypervalidation<-diff_CpG_dbvalidation$CpG[which((diff_CpG_dbvalidation$db)<=(-0.15))] #  74

EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194

diff_genes_db_hypovalidation<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypovalidation)] ) #11442
diff_genes_db_hypervalidation<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypervalidation)] ) # 2084

write.table(diff_genes_db_hypovalidation, file=here("data/validation/DNAm/","validation_genes_hypomethylation_hilo.txt"), quote=F, row.names = F, col.names = F)
write.table(diff_genes_db_hypervalidation, file=here("data/validation/DNAm/","validation_genes_hypermethylation_hilo.txt"), quote=F, row.names = F, col.names = F)





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
