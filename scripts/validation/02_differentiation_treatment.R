#'---
#'title: Differentiation and treatment organoids
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




load(file=here("data/validation/DNAm","validation_betas_normalized.RData"))



################################################################################
#' ## Differential methylation with differentiation 
################################################################################
#########
#' ## low passage
#########
validation_epic.organoid_low<-validation_epic.organoid[which(validation_epic.organoid$passage_hilo=="low" & validation_epic.organoid$comparison=="differentiation"),]
table(validation_epic.organoid_low$differentiation, validation_epic.organoid_low$individual)

validation_organoid_beta_low<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%validation_epic.organoid_low$array.id)]
identical(validation_epic.organoid_low$array.id, colnames(validation_organoid_beta_low))

validation_epic.organoid_low$individual<-as.factor(validation_epic.organoid_low$individual)

mod<-model.matrix(~ differentiation + Segment + individual, data=validation_epic.organoid_low)
fit <- lmFit(validation_organoid_beta_low, mod)
ebfit <- eBayes(fit)

# Delta beta
low_pass_diff_db<-sapply(1:nrow(validation_organoid_beta_low), function(x){
  sampleinfo_cpg<-validation_epic.organoid_low
  sampleinfo_cpg$beta<-as.numeric(validation_organoid_beta_low[x,])
  mns<-tapply(sampleinfo_cpg$beta,  sampleinfo_cpg$differentiation, mean)
  mns[2]-mns[1]# UD minus D
})


differnetiation_low_stats<-data.frame(p.value=ebfit$p.value[,"differentiationUD"], CpG=rownames(validation_organoid_beta_low), db=low_pass_diff_db)

# Adjust P values
differnetiation_low_stats$p_adjusted<-p.adjust(differnetiation_low_stats$p.value, method="BH")

diff_CpG_lowdiff<-differnetiation_low_stats[which(differnetiation_low_stats$p_adjusted<0.05 & abs(differnetiation_low_stats$db)>0.15),] #
diff_CpG_db_lowdiff<-diff_CpG_lowdiff$CpG[which((diff_CpG_lowdiff$db)>=0.15)] #  
diff_CpG_db_lowdiff<-diff_CpG_lowdiff$CpG[which((diff_CpG_lowdiff$db)<=(-0.15))] #  



#########
#' ## high passage
#########
validation_epic.organoid_high<-validation_epic.organoid[which(validation_epic.organoid$passage_hilo=="high" & validation_epic.organoid$comparison=="differentiation"),]
table(validation_epic.organoid_high$differentiation, validation_epic.organoid_high$individual)

validation_organoid_beta_high<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%validation_epic.organoid_high$array.id)]
identical(validation_epic.organoid_high$array.id, colnames(validation_organoid_beta_high))

validation_epic.organoid_high$individual<-as.factor(validation_epic.organoid_high$individual)

mod<-model.matrix(~ differentiation + Segment + individual, data=validation_epic.organoid_high)
fit <- lmFit(validation_organoid_beta_high, mod)
ebfit <- eBayes(fit)

# Delta beta
high_pass_diff_db<-sapply(1:nrow(validation_organoid_beta_high), function(x){
  sampleinfo_cpg<-validation_epic.organoid_high
  sampleinfo_cpg$beta<-as.numeric(validation_organoid_beta_high[x,])
  mns<-tapply(sampleinfo_cpg$beta,  sampleinfo_cpg$differentiation, mean)
  mns[2]-mns[1]# UD minus D
})


differnetiation_high_stats<-data.frame(p.value=ebfit$p.value[,"differentiationUD"], CpG=rownames(validation_organoid_beta_high), db=high_pass_diff_db)

# Adjust P values
differnetiation_high_stats$p_adjusted<-p.adjust(differnetiation_high_stats$p.value, method="BH")

diff_CpG_highdiff<-differnetiation_high_stats[which(differnetiation_high_stats$p_adjusted<0.05 & abs(differnetiation_high_stats$db)>0.15),] #
diff_CpG_db_highdiff<-diff_CpG_highdiff$CpG[which((diff_CpG_highdiff$db)>=0.15)] #  
diff_CpG_db_highdiff<-diff_CpG_highdiff$CpG[which((diff_CpG_highdiff$db)<=(-0.15))] #  




#' ### CpG to gene associations
EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194

diff_genes_db_hypovalidation<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypovalidation)] ) #11442
diff_genes_db_hypervalidation<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypervalidation)] ) # 2084

write.table(diff_genes_db_hypovalidation, file=here("data/validation/DNAm/","validation_genes_hypomethylation.txt"), quote=F, row.names = F, col.names = F)
write.table(diff_genes_db_hypervalidation, file=here("data/validation/DNAm/","validation_genes_hypermethylation.txt"), quote=F, row.names = F, col.names = F)



#######################################################################
#' ## Differential methylation with IFNg 
#######################################################################
#' ## low passage
validation_epic.organoid_low_IFNg<-validation_epic.organoid[which(validation_epic.organoid$passage_hilo=="low" & validation_epic.organoid$comparison=="cytokine" & validation_epic.organoid$treatment!="TNFa"),]
table(validation_epic.organoid_low_IFNg$treatment, validation_epic.organoid_low_IFNg$individual)

validation_organoid_beta_low_IFNg<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%validation_epic.organoid_low_IFNg$array.id)]
identical(validation_epic.organoid_low_IFNg$array.id, colnames(validation_organoid_beta_low_IFNg))

validation_epic.organoid_low_IFNg$individual<-as.factor(validation_epic.organoid_low_IFNg$individual)

mod<-model.matrix(~ treatment  + individual, data=validation_epic.organoid_low_IFNg)
fit <- lmFit(validation_organoid_beta_low_IFNg, mod)
ebfit <- eBayes(fit)

# Delta beta
low_pass_IFNg_db<-sapply(1:nrow(validation_organoid_beta_low_IFNg), function(x){
  sampleinfo_cpg<-validation_epic.organoid_low_IFNg
  sampleinfo_cpg$beta<-as.numeric(validation_organoid_beta_low_IFNg[x,])
  mns<-tapply(sampleinfo_cpg$beta,  sampleinfo_cpg$treatment, mean)
  mns[2]-mns[1]# UT minus IFNg
})


IFNg_low_stats<-data.frame(p.value=ebfit$p.value[,"treatmentUT"], CpG=rownames(validation_organoid_beta_low_IFNg), db=low_pass_IFNg_db)

# Adjust P values
IFNg_low_stats$p_adjusted<-p.adjust(IFNg_low_stats$p.value, method="BH")

diff_CpG_lowIFNg<-IFNg_low_stats[which(IFNg_low_stats$p_adjusted<0.05 & abs(differnetiation_low_stats$db)>0.15),] #
diff_CpG_lowIFNg_hypo<-diff_CpG_lowIFNg$CpG[which((diff_CpG_lowIFNg$db)>=0.15)] #  
diff_CpG_lowIFNg_hyper<-diff_CpG_lowIFNg$CpG[which((diff_CpG_lowIFNg$db)<=(-0.15))] #  


#########
#' ## high passage
#########
validation_epic.organoid_high_IFNg<-validation_epic.organoid[which(validation_epic.organoid$passage_hilo=="high" & validation_epic.organoid$comparison=="cytokine" & validation_epic.organoid$treatment!="TNFa"),]
table(validation_epic.organoid_high_IFNg$treatment, validation_epic.organoid_high_IFNg$individual)

validation_organoid_beta_high_IFNg<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%validation_epic.organoid_high_IFNg$array.id)]
identical(validation_epic.organoid_high_IFNg$array.id, colnames(validation_organoid_beta_high_IFNg))

validation_epic.organoid_high_IFNg$individual<-as.factor(validation_epic.organoid_high_IFNg$individual)

mod<-model.matrix(~ treatment  + individual, data=validation_epic.organoid_high_IFNg)
fit <- lmFit(validation_organoid_beta_high_IFNg, mod)
ebfit <- eBayes(fit)

# Delta beta
high_pass_IFNg_db<-sapply(1:nrow(validation_organoid_beta_high_IFNg), function(x){
  sampleinfo_cpg<-validation_epic.organoid_high_IFNg
  sampleinfo_cpg$beta<-as.numeric(validation_organoid_beta_high_IFNg[x,])
  mns<-tapply(sampleinfo_cpg$beta,  sampleinfo_cpg$treatment, mean)
  mns[2]-mns[1]# UT minus IFNg
})


IFNg_high_stats<-data.frame(p.value=ebfit$p.value[,"treatmentUT"], CpG=rownames(validation_organoid_beta_high_IFNg), db=high_pass_IFNg_db)

# Adjust P values
IFNg_high_stats$p_adjusted<-p.adjust(IFNg_high_stats$p.value, method="BH")

diff_CpG_highIFNg<-IFNg_high_stats[which(IFNg_high_stats$p_adjusted<0.05 & abs(differnetiation_high_stats$db)>0.15),] #
diff_CpG_highIFNg_hypo<-diff_CpG_highIFNg$CpG[which((diff_CpG_highIFNg$db)>=0.15)] #  
diff_CpG_highIFNg_hyper<-diff_CpG_highIFNg$CpG[which((diff_CpG_highIFNg$db)<=(-0.15))] #  
