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
  library(dplyr)
})
options(stringsAsFactors = FALSE)
source(here("general_functions/00_pretty_plots.R"))

#################
#'# load adjusted RNAseq pediatric data
#################
load(file=here("rna_seq/data/ibd_adjusted_rnaseq.RData"))
sampleinfo_rnaseq$case.no<-as.character(sampleinfo_rnaseq$case.no)
sampleinfo_rnaseq$path<-NULL

#' gene ID convert
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ttg <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version",
                                     "ensembl_gene_id", "external_gene_name", "description",
                                     "transcript_biotype"),mart = mart)
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id', 'ens_gene', 'ext_gene'))
ttg$transcript_IDversion<-paste(ttg$target_id, ttg$transcript_version, sep=".")
head(ttg)

#'## might be best to just take the highest expressed transcript for comparisons
ibd_rnaseq_adjusted_top<-as.data.frame(ibd_rnaseq_adjusted)
ibd_rnaseq_adjusted_top$transcript<-rownames(ibd_rnaseq_adjusted)
ibd_rnaseq_adjusted_top<-merge(ibd_rnaseq_adjusted_top, ttg,by.x="transcript",by.y="ens_gene")
transcript_means<-data.frame(transcript=ibd_rnaseq_adjusted_top$transcript, gene=ibd_rnaseq_adjusted_top$ext_gene, mn=rowMeans(select(ibd_rnaseq_adjusted_top, -c(1,80:82))))
transcript_max_expressed<-transcript_means %>% group_by(gene) %>% slice(which.max(mn))

ibd_rnaseq_adjusted_top<-ibd_rnaseq_adjusted[which(rownames(ibd_rnaseq_adjusted)%in%transcript_max_expressed$transcript),]

sampleinfo_rnaseq_controls<-sampleinfo_rnaseq[which(sampleinfo_rnaseq$diagnosis=="Control"),]
ibd_rnaseq_adjusted_top_controls<-ibd_rnaseq_adjusted_top[,which(colnames(ibd_rnaseq_adjusted_top)%in%sampleinfo_rnaseq_controls$sample)]
identical(colnames(ibd_rnaseq_adjusted_top_controls),sampleinfo_rnaseq_controls$sample)

ttg_unique<-ttg[!duplicated(ttg[,c("ens_gene","ext_gene")]),c("ens_gene","ext_gene")]


gene_exp_plot<-function(gene){
  goi<-as.data.frame(ibd_rnaseq_adjusted_top_controls[which(rownames(ibd_rnaseq_adjusted_top_controls)%in%(unique(ttg$ens_gene[which(ttg$ext_gene%in%gene)]))),])
    goi$gene_ID<-rownames(goi)
    goi<-melt(goi)
    colnames(goi)[3]<-"scaled_reads_per_base"
    plt<-merge(sampleinfo_rnaseq_controls,goi, by.x="sample",by.y="variable")
    plt<-merge(plt, ttg_unique, by.x="gene_ID", by.y="ens_gene")
    plt$label<-paste(plt$ext_gene, "\n(", plt$gene_ID,")", sep="")
    ggplot(plt, aes(sample.site, scaled_reads_per_base,fill=sample.site))+
      geom_boxplot()+geom_point(aes(fill=sample.site), shape=21, color="black")+
      facet_wrap(~label, nrow=1)+theme_bw()+th+fillscale_sampsite+
      ylab("Scaled Reads Per Base (TPM)")+xlab("")+theme(axis.text = element_text(size=12))}


# howell saw diff exp  "IRF1"  "NLRC5" "TAP1"  "PSMB8" "PSMB9" i.e respectively "ENSG00000125347" "ENSG00000140853" "ENSG00000168394" "ENSG00000204264" "ENSG00000240065"
gene_exp_plot(c("IRF1",  "NLRC5", "TAP1" , "PSMB8", "PSMB9"))

gene_exp_plot(c("NLRC5", "TAP1" ,"TAP2" , "PSMB8", "PSMB9"))

gene_exp_plot(c("NLRC5", "TAP1" ,"TAP2" , "PSMB8", "PSMB9","HLA-DOA"))


# ggsave(file=here("rna_seq/figs/Bulk_RNA_seq_MHC1_TI.pdf"), w=12, h=6)
# ggsave(file=here("rna_seq/figs/jpeg/Bulk_RNA_seq_MHC1_TI.jpeg"), w=12, h=6)





#################
#' Correlation NLRC5 and MHCI
#################
ibd_rnaseq_adjusted_top_controls



gene_cor_plot<-function(gene){
  goi<-as.data.frame(ibd_rnaseq_adjusted_top_controls[which(rownames(ibd_rnaseq_adjusted_top_controls)%in%(unique(ttg$ens_gene[which(ttg$ext_gene%in%gene)]))),])
  goi<-as.data.frame(t(goi))
  colnames(goi)[1]<-paste(unique(ttg$ext_gene[grep(colnames(goi)[1], ttg$ens_gene)]),"\n(", colnames(goi)[1], ")",sep="")
  colnames(goi)[2]<-paste(unique(ttg$ext_gene[grep(colnames(goi)[2], ttg$ens_gene)]),"\n(", colnames(goi)[2], ")",sep="")
  goi$sample<-rownames(goi)
  
  plt<-merge(sampleinfo_rnaseq_controls,goi, by="sample")

  ggplot(plt, aes(plt[,17], plt[,18]))+
    geom_point(aes(fill=sample.site), shape=21, color="black")+
    theme_bw()+th+fillscale_sampsite+
    theme(axis.text = element_text(size=12))+
    xlab(colnames(plt)[17])+ylab(colnames(plt)[18])+
    stat_smooth(se=F, method="lm",color="#1f66e5")
  ggsave(file=paste(here("rna_seq/figs/jpeg/"),gene[1],"_",gene[2],"_controlsOnly_correlation_bulk.jpeg", sep=""), w=5, h=3.5)
  ggsave(file=paste(here("rna_seq/figs/"),gene[1],"_",gene[2],"_controlsOnly_correlation_bulk.pdf", sep=""), w=5, h=3.5)  }

gene_cor_plot(c("NLRC5",  "TAP1"))



## Correlation heatmap
MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","MR1","CD1D","IRF1","NLRC5")

gene_combos<-combn(MHCI, 2)

gene_correlation_values<-do.call(rbind, lapply(1:ncol(gene_combos), function(x){
  gene_name1=gene_combos[1,x]
  gene_name2=gene_combos[2,x]
  plt_gene<-as.data.frame(t(as.data.frame(ibd_rnaseq_adjusted_top_controls[which(rownames(ibd_rnaseq_adjusted_top_controls)%in%(unique(ttg$ens_gene[which(ttg$ext_gene%in%c(gene_name1, gene_name2))]))),])))
  
  allcor<-cor(plt_gene[,1],plt_gene[,2], method="spearman")

  data.frame(gene_name1=gene_name1, gene_name2=gene_name2, allcor=allcor)
}))

gene_correlation_values$gene_name1<-factor(gene_correlation_values$gene_name1, levels=rev(MHCI))
gene_correlation_values$gene_name2<-factor(gene_correlation_values$gene_name2, levels=rev(MHCI))

selected_labels<-gene_correlation_values[which(gene_correlation_values$gene_name2=="NLRC5"),]

ggplot()+geom_tile(aes(gene_name1, gene_name2, fill=allcor),gene_correlation_values, color="black")+
  geom_text(aes(gene_name1, gene_name2, label=round(allcor,2)), selected_labels, color="grey10", size=4)+
  scale_fill_gradient2(high="#2171b5",mid="#f3f7fc",low="red",
                       na.value="white", midpoint=0, name="Correlation")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size =12, color="black"),
        legend.text = element_text(size =10),
        legend.title = element_text(size =11),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_y_discrete(position = "right") 


ggsave(file=here("rna_seq/figs/jpeg","MHCI_score_controlsOnly_correlation_bulk_pediatric.jpeg"), w=10, h=9)
ggsave(file=here("rna_seq/figs","MHCI_score_controlsOnly_correlation_bulk_pediatric.pdf"), w=10, h=9)



#'## R Session Info
sessionInfo()


