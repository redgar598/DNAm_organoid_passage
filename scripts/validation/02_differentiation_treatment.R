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

mod<-model.matrix(~ differentiation + individual, data=validation_epic.organoid_low)
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

mod<-model.matrix(~ differentiation  + individual, data=validation_epic.organoid_high)
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



#' 
#' #' ### CpG to gene associations
#' EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194
#' 
#' diff_genes_db_hypovalidation<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypovalidation)] ) #11442
#' diff_genes_db_hypervalidation<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypervalidation)] ) # 2084
#' 
#' write.table(diff_genes_db_hypovalidation, file=here("data/validation/DNAm/","validation_genes_hypomethylation.txt"), quote=F, row.names = F, col.names = F)
#' write.table(diff_genes_db_hypervalidation, file=here("data/validation/DNAm/","validation_genes_hypermethylation.txt"), quote=F, row.names = F, col.names = F)



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

diff_CpG_lowIFNg<-IFNg_low_stats[which(IFNg_low_stats$p_adjusted<0.05 & abs(IFNg_low_stats$db)>0.15),] #
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

diff_CpG_highIFNg<-IFNg_high_stats[which(IFNg_high_stats$p_adjusted<0.05 & abs(IFNg_high_stats$db)>0.15),] #
diff_CpG_highIFNg_hypo<-diff_CpG_highIFNg$CpG[which((diff_CpG_highIFNg$db)>=0.15)] #  
diff_CpG_highIFNg_hyper<-diff_CpG_highIFNg$CpG[which((diff_CpG_highIFNg$db)<=(-0.15))] #  





#######################################################################
#' ## Differential methylation with TNFa 
#######################################################################
#' ## low passage
validation_epic.organoid_low_TNFa<-validation_epic.organoid[which(validation_epic.organoid$passage_hilo=="low" & validation_epic.organoid$comparison=="cytokine" & validation_epic.organoid$treatment!="IFNg"),]
table(validation_epic.organoid_low_TNFa$treatment, validation_epic.organoid_low_TNFa$individual)

validation_organoid_beta_low_TNFa<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%validation_epic.organoid_low_TNFa$array.id)]
identical(validation_epic.organoid_low_TNFa$array.id, colnames(validation_organoid_beta_low_TNFa))

validation_epic.organoid_low_TNFa$individual<-as.factor(validation_epic.organoid_low_TNFa$individual)

mod<-model.matrix(~ treatment  + individual, data=validation_epic.organoid_low_TNFa)
fit <- lmFit(validation_organoid_beta_low_TNFa, mod)
ebfit <- eBayes(fit)

# Delta beta
low_pass_TNFa_db<-sapply(1:nrow(validation_organoid_beta_low_TNFa), function(x){
  sampleinfo_cpg<-validation_epic.organoid_low_TNFa
  sampleinfo_cpg$beta<-as.numeric(validation_organoid_beta_low_TNFa[x,])
  mns<-tapply(sampleinfo_cpg$beta,  sampleinfo_cpg$treatment, mean)
  mns[2]-mns[1]# UT minus TNFa
})


TNFa_low_stats<-data.frame(p.value=ebfit$p.value[,"treatmentUT"], CpG=rownames(validation_organoid_beta_low_TNFa), db=low_pass_TNFa_db)

# Adjust P values
TNFa_low_stats$p_adjusted<-p.adjust(TNFa_low_stats$p.value, method="BH")

diff_CpG_lowTNFa<-TNFa_low_stats[which(TNFa_low_stats$p_adjusted<0.05 & abs(TNFa_low_stats$db)>0.15),] #
diff_CpG_lowTNFa_hypo<-diff_CpG_lowTNFa$CpG[which((diff_CpG_lowTNFa$db)>=0.15)] #  
diff_CpG_lowTNFa_hyper<-diff_CpG_lowTNFa$CpG[which((diff_CpG_lowTNFa$db)<=(-0.15))] #  


#########
#' ## high passage
#########
validation_epic.organoid_high_TNFa<-validation_epic.organoid[which(validation_epic.organoid$passage_hilo=="high" & validation_epic.organoid$comparison=="cytokine" & validation_epic.organoid$treatment!="IFNg"),]
table(validation_epic.organoid_high_TNFa$treatment, validation_epic.organoid_high_TNFa$individual)

validation_organoid_beta_high_TNFa<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%validation_epic.organoid_high_TNFa$array.id)]
identical(validation_epic.organoid_high_TNFa$array.id, colnames(validation_organoid_beta_high_TNFa))

validation_epic.organoid_high_TNFa$individual<-as.factor(validation_epic.organoid_high_TNFa$individual)

mod<-model.matrix(~ treatment  + individual, data=validation_epic.organoid_high_TNFa)
fit <- lmFit(validation_organoid_beta_high_TNFa, mod)
ebfit <- eBayes(fit)

# Delta beta
high_pass_TNFa_db<-sapply(1:nrow(validation_organoid_beta_high_TNFa), function(x){
  sampleinfo_cpg<-validation_epic.organoid_high_TNFa
  sampleinfo_cpg$beta<-as.numeric(validation_organoid_beta_high_TNFa[x,])
  mns<-tapply(sampleinfo_cpg$beta,  sampleinfo_cpg$treatment, mean)
  mns[2]-mns[1]# UT minus IFNg
})


TNFa_high_stats<-data.frame(p.value=ebfit$p.value[,"treatmentUT"], CpG=rownames(validation_organoid_beta_high_TNFa), db=high_pass_TNFa_db)

# Adjust P values
TNFa_high_stats$p_adjusted<-p.adjust(TNFa_high_stats$p.value, method="BH")

diff_CpG_highIFNg<-TNFa_high_stats[which(TNFa_high_stats$p_adjusted<0.05 & abs(TNFa_high_stats$db)>0.15),] #
diff_CpG_highTNFa_hypo<-diff_CpG_highIFNg$CpG[which((diff_CpG_highIFNg$db)>=0.15)] #  
diff_CpG_highTNFa_hyper<-diff_CpG_highIFNg$CpG[which((diff_CpG_highIFNg$db)<=(-0.15))] #  




### plot whatever gene DMR
EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194

gene_DNAm<-function(gene, condition){
    CpG_goi<-EPIC_genes[which(EPIC_genes$Gene.name%in%gene & EPIC_genes$CpG_in=="promoter"),]
    CpGs<-unique(CpG_goi$IlmnID)
    
    if(condition=="differentiation"){
      sample_info<-validation_epic.organoid[which(validation_epic.organoid$comparison=="differentiation"),]
      betas<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%sample_info$array.id)]
      columnX<-"condition"
      fills<-scale_fill_manual(values=c("#abd9e9","#a1d99b"), name="Differentiation")
      }else{
        if(condition=="TNFa"){
          sample_info<-validation_epic.organoid[which(validation_epic.organoid$comparison=="cytokine" & validation_epic.organoid$treatment!="IFNg"),]
          betas<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%sample_info$array.id)]
          columnX<-"treatment"
          fills<-scale_fill_manual(values=c("grey80","firebrick4"), name="Treatment")}else{
            if(condition=="IFNg"){
              sample_info<-validation_epic.organoid[which(validation_epic.organoid$comparison=="cytokine" & validation_epic.organoid$treatment!="TNFa"),]
              betas<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%sample_info$array.id)]
              columnX<-"treatment"
              fills<-scale_fill_manual(values=c("grey80","cornflowerblue"), name="Treatment")
              }}}
    
    betas<-betas[which(rownames(betas)%in%CpGs),]
    betas<-melt(betas)
    organoid_plt<-merge(sample_info, betas, by.x="array.id",by.y="Var2")
    organoid_plt<-merge(organoid_plt, CpG_goi, by.x="Var1",by.y="IlmnID")
    
    ggplot(organoid_plt, aes(Var1, value, fill=organoid_plt[,which(colnames(organoid_plt)==columnX)]))+
      geom_boxplot()+
      geom_point(size=1.5,shape=21, color="black",position = position_dodge(width=0.75))+
      th+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      facet_wrap(~passage_hilo)+cols+fills+
      ylab("DNAm Beta")+xlab("")+ylim(0,1)}

gene<-"LGR5"
condition<-"IFNg"

gene_DNAm("LGR5","differentiation")
gene_DNAm("FABP1","differentiation")
gene_DNAm("KRT19","differentiation")

gene_DNAm("NLRC5","IFNg")
gene_DNAm("TAP1","IFNg")
gene_DNAm("HLA-A","IFNg")
gene_DNAm("PSMB9","IFNg")



# DMR
# sum_stat<-melt(tapply(organoid_plt$value, list(organoid_plt$MAPINFO,organoid_plt$condition,organoid_plt$passage_hilo), mean))
# colnames(sum_stat)<-c("MAPINFO","condition","passage_hilo","value")
# sum_stat<-sum_stat[which(!(is.na(sum_stat$value))),]
# 
# ggplot()+
#   geom_point(aes(MAPINFO, value, fill=condition), organoid_plt, size=1.5,shape=21, color="white")+
#   geom_line(aes(MAPINFO,value,color=condition),sum_stat, size=1)+th+theme_bw()+
#   facet_wrap(~passage_hilo)+cols+fills+
#   ylab("DNAm Beta")+xlab(xaxis)+ylim(0,1)
# 
