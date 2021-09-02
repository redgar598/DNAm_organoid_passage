#'---
#'title: Differential expresion gene level
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---



#load sleuth library
suppressMessages({
  library(sleuth)
  library(testit)
  library(ggplot2)
  library(RColorBrewer)
  library(reshape2)
  library(gridExtra)
  library(here)
})
options(stringsAsFactors = FALSE)
source(here("general_functions/00_pretty_plots.R"))

num_cores=1

#set input and output dirs
datapath = here("../data/public_rna_seq/trimmed/kallisto")
resultdir = here('../data/public_rna_seq/trimmed/kallisto/sleuth')

sampleinfo <- read.table(here("../data/public_rna_seq/Zerbino_RNASeq_SampleInfo.txt"), header=T, sep="\t")
# ont need ftp location
sampleinfo<-sampleinfo[,1:13]

#'###create a sample to condition metadata description
sample_id = sampleinfo$ena.id
kal_dirs <- file.path(datapath, sample_id)

sampleinfo$path<-kal_dirs
colnames(sampleinfo)[colnames(sampleinfo) =="ena.id"]<-"sample"

sampleinfo$inflammation<-as.factor(sampleinfo$inflammation)
sampleinfo$sex<-as.factor(sampleinfo$sex)
sampleinfo$diagnosis<-as.factor(sampleinfo$diagnosis)

#' ignore outlier
sampleinfo_noOutlier<-sampleinfo[which(sampleinfo$sample!="ERR2271027"),]



###########################
#' ## Cell type counts from DNAm
###########################
sampleinfo_rnaseq<-sampleinfo_noOutlier
sampleinfo_rnaseq$sample_ID<-paste(sampleinfo_rnaseq$case.no, sampleinfo_rnaseq$sample.site)

load(file=here("DNAm/data/DNAm_organoid_counts.RData"))
sampleinfo_DNAm_countsConstrained<-sampleinfo_DNAm_countsConstrained[which(is.na(sampleinfo_DNAm_countsConstrained$passage.or.rescope.no)),]
sampleinfo_DNAm_countsConstrained$sample_ID<-paste(sampleinfo_DNAm_countsConstrained$case.no, sampleinfo_DNAm_countsConstrained$sample.site)


sampleinfo_rnaseq<-merge(sampleinfo_rnaseq, sampleinfo_DNAm_countsConstrained[,c("sample_ID","organoid", "tcell")], by="sample_ID")

sampleinfo_rnaseq_TI<-sampleinfo_rnaseq[which(sampleinfo_rnaseq$sample.site=="TI"),]
sampleinfo_rnaseq_TI$diagnosis_grouped<-sampleinfo_rnaseq_TI$diagnosis
levels(sampleinfo_rnaseq_TI$diagnosis_grouped)<-c("CD","Control_UC","Control_UC")

sampleinfo_rnaseq_SC<-sampleinfo_rnaseq[which(sampleinfo_rnaseq$sample.site=="SC"),]
sampleinfo_rnaseq_SC$diagnosis_grouped<-sampleinfo_rnaseq_SC$diagnosis
levels(sampleinfo_rnaseq_SC$diagnosis_grouped)<-c("IBD","Control","IBD")

###########
#'# aggergate to gene level
###########
mart <- biomaRt::useMart(biomart = "ensembl",
                         dataset = "hsapiens_gene_ensembl")

ttg <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "transcript_version",
                 "ensembl_gene_id", "external_gene_name", "description",
                 "transcript_biotype"),
  mart = mart)
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id', 'ens_gene', 'ext_gene'))
head(ttg)

#############
#'# run sleuth on the tissue seperated data
#############

#'## TI
so_TI <- sleuth_prep(sampleinfo_rnaseq_TI, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)
#so_TI <- sleuth_prep(sampleinfo_rnaseq_TI, extra_bootstrap_summary =TRUE) # transcript level

so_TI <- sleuth_fit(so_TI, ~0+diagnosis_grouped+organoid+age+sex+inflammation, 'full')
so_TI <- sleuth_fit(so_TI, ~0+organoid+age+sex+inflammation, 'reduced')
so_TI <- sleuth_lrt(so_TI, 'reduced', 'full')
models(so_TI)


#' ### summarize the sleuth results and view 20 most significant DE transcripts
sleuth_table_TI <- sleuth_results(so_TI, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_TI <- dplyr::filter(sleuth_table_TI, qval <= 0.05)
head(sleuth_significant_TI, 20)


#'### plot an example DE gene result
plot_bootstrap(so_TI, target_id ="ENSG00000140853", units = "scaled_reads_per_base", color_by = "diagnosis")+fillscale_diagnosis
plot_bootstrap(so_TI, target_id ="ENSG00000159212", units = "scaled_reads_per_base", color_by = "diagnosis")+fillscale_diagnosis

sleuth_table_TI[which(sleuth_table_TI$target_id=="ENSG00000140853"),]
sleuth_table_TI[which(sleuth_table_TI$ext_gene=="TAP1"),]
MHCI_NLRC5 = c("TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")

sleuth_table_TI[which(sleuth_table_TI$ext_gene%in%MHCI_NLRC5),]
sleuth_significant_TI[which(sleuth_significant_TI$ext_gene%in%MHCI_NLRC5),]

mat <- sleuth:::spread_abundance_by(so_TI$obs_norm, "scaled_reads_per_base",  so_TI$sample_to_covariates$sample)


#'## gene plots
gene_exp_plot<-function(gene){
  goi<-as.data.frame(mat[which(rownames(mat)%in%(unique(ttg$ens_gene[which(ttg$ext_gene==gene)]))),])
  if(ncol(goi)==1){
    goi$sample_ID<-rownames(goi)
    colnames(goi)[1]<-"scaled_reads_per_base"
    plt<-merge(sampleinfo_rnaseq_TI,goi, by.x="sample",by.y="sample_ID")
    ggplot(plt, aes(diagnosis, scaled_reads_per_base,fill=diagnosis))+geom_boxplot()+geom_point(aes(fill=diagnosis), shape=21, color="black")+theme_bw()+th+fillscale_diagnosis
  }else{
    goi$gene_ID<-rownames(goi)
    goi<-melt(goi)
    colnames(goi)[3]<-"scaled_reads_per_base"
    plt<-merge(sampleinfo_rnaseq_TI,goi, by.x="sample",by.y="variable")
    ggplot(plt, aes(diagnosis, scaled_reads_per_base,fill=diagnosis))+geom_boxplot()+geom_point(aes(fill=diagnosis), shape=21, color="black")+facet_wrap(~gene_ID)+theme_bw()+th+fillscale_diagnosis
    }}

# howell saw diff exp  "IRF1"  "NLRC5" "TAP1"  "PSMB8" "PSMB9" i.e respectively "ENSG00000125347" "ENSG00000140853" "ENSG00000168394" "ENSG00000204264" "ENSG00000240065"
gene_exp_plot("IRF1")
gene_exp_plot("NLRC5")
gene_exp_plot("TAP1")
gene_exp_plot("PSMB8")
gene_exp_plot("PSMB9")
## not sig howell but in dnam
gene_exp_plot("B2M")


geneID_exp_plot<-function(gene_ID){
  goi<-as.data.frame(mat[which(rownames(mat)%in%(unique(ttg$ens_gene[which(ttg$ens_gene==gene_ID)]))),])
  goi$sample_ID<-rownames(goi)
  colnames(goi)[1]<-"scaled_reads_per_base"
  plt<-merge(sampleinfo_rnaseq_TI,goi, by.x="sample",by.y="sample_ID")
  label<-unique(ttg$ext_gene[which(ttg$ens_gene==gene_ID)])
  plt$diagnosis<-factor(plt$diagnosis, levels = c("Control","UC","CD"))
  ggplot(plt, aes(diagnosis, scaled_reads_per_base,fill=diagnosis))+geom_boxplot()+
    geom_point(shape=21, size=2, color="black")+theme_bw()+th+fillscale_diagnosis+ggtitle(label)}
    
grid.arrange(geneID_exp_plot("ENSG00000125347"),#IRF1
  geneID_exp_plot("ENSG00000140853"),#NLRC5
  geneID_exp_plot("ENSG00000168394"),#TAP1
  geneID_exp_plot("ENSG00000204264"),#PSMB8
  geneID_exp_plot("ENSG00000240065"))#PSMB9


# one big plot
goi<-as.data.frame(mat[which(rownames(mat)%in%(  unique(ttg$ens_gene[which(ttg$ens_gene%in%c("ENSG00000125347", "ENSG00000140853", "ENSG00000168394", "ENSG00000204264", "ENSG00000240065","ENSG00000166710"))]))),])
  goi$gene_ID<-rownames(goi)
  goi<-melt(goi)
  colnames(goi)[3]<-"scaled_reads_per_base"
  plt<-merge(sampleinfo_rnaseq_TI,goi, by.x="sample",by.y="variable")
  ttg_min<-ttg[,2:3][!duplicated(ttg[,2:3]),]
  ttg_min$label<-paste(ttg_min$ext_gene, " - ", ttg_min$ens_gene, sep="")
  plt<-merge(plt, ttg_min, by.x="gene_ID","ens_gene")
  
  ggplot(plt, aes(diagnosis, scaled_reads_per_base,fill=diagnosis))+
    geom_boxplot()+geom_point(aes(fill=diagnosis), shape=21, color="black")+
    facet_wrap(~label, scale="free_y")+theme_bw()+th+fillscale_diagnosis+
    ylab("Scaled reads per base")+xlab("Diagnosis")


ggsave(file=here("rna_seq/figs/Bulk_RNA_seq_MHC1_TI.pdf"), w=12, h=6)
ggsave(file=here("rna_seq/figs/jpeg/Bulk_RNA_seq_MHC1_TI.jpeg"), w=12, h=6)


#' 
#' #'## SC
#' so_SC <- sleuth_prep(sampleinfo_rnaseq_SC, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)
#' 
#' so_SC <- sleuth_fit(so_SC, ~0+diagnosis_grouped+organoid+age+sex+inflammation, 'full')
#' so_SC <- sleuth_fit(so_SC, ~0+organoid+age+sex+inflammation, 'reduced')
#' so_SC <- sleuth_lrt(so_SC, 'reduced', 'full')
#' models(so_SC)
#' 
#' #' ### summarize the sleuth results and view 20 most significant DE transcripts
#' sleuth_table_SC <- sleuth_results(so_SC, 'reduced:full', 'lrt', show_all = FALSE)
#' sleuth_significant_SC <- dplyr::filter(sleuth_table_SC, qval <= 0.05)
#' head(sleuth_significant_SC, 20)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #################
#' ## Without cell covariate
#' ################
#' #so_TI <- sleuth_prep(sampleinfo_rnaseq_TI, extra_bootstrap_summary =TRUE) 
#' so_TI_nocell <- sleuth_prep(sampleinfo_rnaseq_TI, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)
#' so_TI_nocell <- sleuth_fit(so_TI_nocell, ~0+diagnosis_grouped+age+sex+inflammation, 'full')
#' so_TI_nocell <- sleuth_fit(so_TI_nocell, ~0+age+sex+inflammation, 'reduced')
#' so_TI_nocell <- sleuth_lrt(so_TI_nocell, 'reduced', 'full')
#' models(so_TI_nocell)
#' 
#' 
#' # #summarize the sleuth results and view 20 most significant DE transcripts
#' sleuth_table_TI_nocell <- sleuth_results(so_TI_nocell, 'reduced:full', 'lrt', show_all = FALSE)
#' sleuth_significant_TI_nocell <- dplyr::filter(sleuth_table_TI_nocell, qval <= 0.05)
#' head(sleuth_significant_TI_nocell, 20)
#' 
#' sleuth_table_TI_nocell[which(sleuth_table_TI_nocell$ext_gene%in%MHCI_NLRC5),]



#'## R Session Info
sessionInfo()


