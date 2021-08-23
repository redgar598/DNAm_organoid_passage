#load sleuth library
suppressMessages({
  library("sleuth")
  library(testit)
  library(ggplot2)
  library(RColorBrewer)
  library(reshape2)
  library(here)
  library(cowplot)
})
options(stringsAsFactors = FALSE)
source(here("general_functions/00_pretty_plots.R"))

#'## All controls are presented in MTAB5015
#'##For the IBD part will use MTAB5464

#set input and output dirs
datapath = here("rna_seq/data/MTAB5015/trimmed/kallisto")
resultdir = here('rna_seq/data/MTAB5015/trimmed/kallisto/sleuth')

sampleinfo <- read.table("../data/public_other/organoid/MTAB5015/E-MTAB-5015.sdrf.txt", header=T, sep="\t")
sampleinfo<-sampleinfo[,c(1:3,6:11,15,16,28,29,31,35,38,40,41,45)]
sample_id = sampleinfo$Comment.ENA_RUN.
kal_dirs <- file.path(datapath, sample_id)


sampleinfo$path<-kal_dirs
colnames(sampleinfo)[colnames(sampleinfo) =="Comment.ENA_RUN."]<-"sample"


table(sampleinfo$Characteristics.sampling.site.,sampleinfo$Characteristics.biosource.type., sampleinfo$Characteristics.developmental.stage.)

## maybe need to remove duplicated rows?
sampleinfo<-sampleinfo[!duplicated(sampleinfo$sample),]



###########
#'# aggregate to gene level
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

###############
#'# run sleuth on the data
###############
so <- sleuth_prep(sampleinfo, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)


#'## PCA my way
mat <- sleuth:::spread_abundance_by(so$obs_norm, "tpm", so$sample_to_covariates$sample)



pca_res <- prcomp(t(mat))
Loadings_meta<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

#'# PC vs PC plot
Loadings_meta$sample<-rownames(Loadings_meta)
Loadings_meta<-merge(Loadings_meta, sampleinfo, by="sample")

ggplot(Loadings_meta, aes(PC1, PC2, fill=Characteristics.sampling.site.))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  th+xlab("PC1 (60%)")+ylab("PC2 (13%)")
ggplot(Loadings_meta, aes(PC1, PC2, fill=Characteristics.developmental.stage.))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  th+xlab("PC1 (60%)")+ylab("PC2 (13%)")
ggplot(Loadings_meta, aes(PC1, PC2, fill=Characteristics.biosource.type.))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  th+xlab("PC1 (60%)")+ylab("PC2 (13%)")

# ggsave(here("rna_seq/figs","PC1_PC2_RNAseq.pdf"), width = 7.5, height = 6)
# ggsave(here("rna_seq/figs/jpeg","PC1_PC2_RNAseq.jpeg"), width = 7.5, height = 6)



#'## might be best to just take the highest expressed transcript for comparisons
mat_top<-as.data.frame(mat)
mat_top$transcript<-rownames(mat_top)
mat_top<-merge(mat_top, ttg,by.x="transcript",by.y="ens_gene")
transcript_means<-data.frame(transcript=mat_top$transcript, gene=mat_top$ext_gene, mn=rowMeans(select(mat_top, -c(1,52:53))))
transcript_max_expressed<-transcript_means %>% group_by(gene) %>% slice(which.max(mn))

mat_top<-mat[which(rownames(mat)%in%transcript_max_expressed$transcript),]


identical(colnames(mat_top),sampleinfo$sample)

ttg_unique<-ttg[!duplicated(ttg[,c("ens_gene","ext_gene")]),c("ens_gene","ext_gene")]

gene_exp_plot<-function(gene){
  goi<-as.data.frame(mat_top[which(rownames(mat_top)%in%(unique(ttg$ens_gene[which(ttg$ext_gene%in%gene)]))),])
  goi$gene_ID<-rownames(goi)
  goi<-melt(goi)
  colnames(goi)[3]<-"scaled_reads_per_base"
  plt<-merge(sampleinfo,goi, by.x="sample",by.y="variable")
  plt<-merge(plt, ttg_unique, by.x="gene_ID", by.y="ens_gene")
  plt$label<-paste(plt$ext_gene, "\n(", plt$gene_ID,")", sep="")
  plt$sampletype<-paste(plt$Characteristics.developmental.stage., "\n", plt$Characteristics.biosource.type.,"", sep="")
  
  plt_fetal<-plt[which(plt$Characteristics.developmental.stage.=="fetal stage"),]
  plt_fetal$Characteristics.sampling.site.<-as.factor(plt_fetal$Characteristics.sampling.site.)
  levels(plt_fetal$Characteristics.sampling.site.)<-c("distal","proximal")
  plt_juvanile<-plt[which(plt$Characteristics.developmental.stage.=="juvenile stage"),]
  
  fetal<-ggplot(plt_fetal, aes(Characteristics.sampling.site., scaled_reads_per_base,fill=Characteristics.sampling.site.))+
    geom_boxplot()+geom_point(aes(fill=Characteristics.sampling.site.), shape=21, color="black")+
    facet_grid(sampletype~label)+theme_bw()+th+ylab("Scaled Reads Per Base (TPM)")+xlab("")+theme(axis.text = element_text(size=12))+
    fillscale_sampsite
  
  juvanile<-ggplot(plt_juvanile, aes(Characteristics.sampling.site., scaled_reads_per_base,fill=Characteristics.sampling.site.))+
    geom_boxplot()+geom_point(aes(fill=Characteristics.sampling.site.), shape=21, color="black")+
    facet_grid(sampletype~label)+theme_bw()+th+ylab("Scaled Reads Per Base (TPM)")+xlab("")+theme(axis.text = element_text(size=12))+
    fillscale_sampsite
  
  plot_grid(fetal, juvanile, align="hv", ncol=1)
  }


# howell saw diff exp  "IRF1"  "NLRC5" "TAP1"  "PSMB8" "PSMB9" i.e respectively "ENSG00000125347" "ENSG00000140853" "ENSG00000168394" "ENSG00000204264" "ENSG00000240065"
gene_exp_plot(c("IRF1",  "NLRC5", "TAP1" , "PSMB8", "PSMB9"))

gene_exp_plot(c("NLRC5", "TAP1" ,"TAP2" , "PSMB8", "PSMB9","B2M","EPCAM"))

gene_exp_plot(c("NLRC5", "TAP1" ,"TAP2" , "PSMB8", "PSMB9","TNF","KRT18"))


