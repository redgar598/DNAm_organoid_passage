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
source(here("scripts/00_pretty_plots.R"))

num_cores=1

#set input and output dirs
datapath = "data/validation/kallisto"
resultdir = here('data/validation/kallisto/sleuth')

sampleinfo <- read.table("data/validation/sample_info_RNA_seq.txt", header=F, sep=" ", skip=1)
sampleinfo<-sampleinfo[,c(2,7,8:10,12)]
colnames(sampleinfo)<-c("sample","concentration","volume","quantity","ratio","well")

#create a sample to condition metadata description
sample_id = sampleinfo$sample
kal_dirs <- file.path(datapath, sample_id)

sampleinfo$path<-kal_dirs


sampleinfo$individual<-sapply(1:nrow(sampleinfo), function(x) strsplit(sampleinfo$sample[x],"_")[[1]][1])
sampleinfo$sample.site<-sapply(1:nrow(sampleinfo), function(x) strsplit(sampleinfo$sample[x],"_")[[1]][2])
sampleinfo$passage<-as.numeric(sapply(1:nrow(sampleinfo), function(x) gsub("p","",strsplit(sampleinfo$sample[x],"_")[[1]][3])))
sampleinfo$passage_hilo<-sapply(1:nrow(sampleinfo), function(x) if(sampleinfo$passage[x]<5){"low"}else{"high"})
sampleinfo$condition<-sapply(1:nrow(sampleinfo), function(x) strsplit(sampleinfo$sample[x],"_")[[1]][4])
sampleinfo$condition<-as.factor(sampleinfo$condition)
levels(sampleinfo$condition)<-c("D","IFNg","IFNg","TNFa","UD","UT")
sampleinfo$condition<-as.character(sampleinfo$condition)
sampleinfo$comparison<-sapply(1:nrow(sampleinfo), function(x) if(sampleinfo$condition[x]%in%c("UD","D")){"differentiation"}else{"cytokine"})
sampleinfo$treatment<-sapply(1:nrow(sampleinfo), function(x) if(sampleinfo$comparison[x]=="cytokine"){sampleinfo$condition[x]}else{"UT"})
sampleinfo$differentiation<-sapply(1:nrow(sampleinfo), function(x) if(sampleinfo$comparison[x]=="differentiation"){sampleinfo$condition[x]}else{"UD"})


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




#############
#'# run sleuth on passage high low all undifferentiated
#############
sampleinfo_UD_site<-sampleinfo[which(sampleinfo$differentiation=="UD" & sampleinfo$treatment=="UT"),]
table(sampleinfo_UD_site$sample.site)
table(sampleinfo_UD_site$sample.site, sampleinfo_UD_site$individual)


#'## TI
so_UD_site <- sleuth_prep(sampleinfo_UD_site, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)
so_UD_site <- sleuth_fit(so_UD_site, ~sample.site + individual , 'full')
so_UD_site <- sleuth_fit(so_UD_site, ~1 + individual, 'reduced')
so_UD_site <- sleuth_lrt(so_UD_site, 'reduced', 'full')
models(so_UD_site)
tests(so_UD_site)


#' ### summarize the sleuth results and view 20 most significant DE transcripts
sleuth_table_UD_site <- sleuth_results(so_UD_site, 'reduced:full', 'lrt', show_all = TRUE)  #show_all = F
sleuth_significant_UD_site <- dplyr::filter(sleuth_table_UD_site, qval <= 0.05)
head(sleuth_significant_UD_site, 20)

sleuth_table_UD_site[grep("WNT5A|LRP5|LRP6",sleuth_table_UD_site$ext_gene),]
sleuth_significant_UD_site[grep("WNT5A|LRP5|LRP6",sleuth_significant_UD_site$ext_gene),]


#'## gene plots
mat <- sleuth:::spread_abundance_by(so_UD_site$obs_norm, "scaled_reads_per_base",  so_UD_site$sample_to_covariates$sample)

gene_exp_plot<-function(gene){
  #if(length(gene)==1){
    goi<-as.data.frame(mat[which(rownames(mat)%in%(unique(ttg$ens_gene[which(ttg$ext_gene==gene)]))),])
    goi$sample_ID<-rownames(goi)
    colnames(goi)[1]<-"scaled_reads_per_base"
    plt<-merge(sampleinfo,goi, by.x="sample",by.y="sample_ID")
    ggplot(plt, aes(sample.site, scaled_reads_per_base))+geom_boxplot()+geom_point(aes(fill=sample.site), shape=21, color="black")+theme_bw()+th+
      scale_fill_manual(values=c("#a6d96a","cornflowerblue"), name="Sample site")#SC TI
  # }else{
  #   goi<-as.data.frame(mat[which(rownames(mat)%in%(unique(ttg$ens_gene[which(ttg$ext_gene%in%gene)]))),])
  #   gene_ID<-ttg[which(ttg$ext_gene%in%gene),2:3]
  #   gene_ID<-gene_ID[!duplicated(gene_ID),]
  #   goi$gene_ID<-rownames(goi)
  #   goi<-melt(goi)
  #   colnames(goi)[3]<-"scaled_reads_per_base"
  #   goi<-merge(goi, gene_ID, by.x="gene_ID", by.y="ens_gene")
  #   goi$label<-paste(goi$ext_gene, "\n(",goi$gene_ID, ")",sep="")
  #   plt<-merge(sampleinfo,goi, by.x="sample",by.y="variable")
  #   plt$Sample_ID<-paste(plt$individual, plt$sample.site, plt$condition)
  #   
  #   ggplot(plt, aes(passage, scaled_reads_per_base))+
  #     geom_line(aes(group=Sample_ID),color="lightgrey")+
  #     stat_smooth(method="lm", color="grey30", size=0.7, se=F)+ylab("Scaled reads per base")+xlab("Passage Number")+
  #     geom_point(aes(fill=as.factor(passage)), shape=21, color="black")+
  #     facet_wrap(~ext_gene, scale="free_y", nrow=1)+theme_bw()+th+
  #     scale_fill_manual(values=pass_col,name="Passage\nNumber", drop=T)+scale_x_continuous(breaks = c(2, 4, 6,8,10,12))+
  #     theme(plot.margin = margin(0.5, 0.15, 0.5, 0.15, "cm"),plot.title = element_text(size=12))  }
    }

# gene level sumary
gene_exp_plot("LRP5")
ggsave(here("../ibd/figs/single_gene_IBD/jpeg","LRP5_differenital_expression_SCTI_passageUTUD.jpeg"), width = 3, height = 2.5)

gene_exp_plot("LRP6")

gene_exp_plot("WNT5A")
ggsave(here("../ibd/figs/single_gene_IBD/jpeg","WNT5A_differenital_expression_SCTI_passageUTUD.jpeg"), width = 3, height = 2.5)


goi<-as.data.frame(mat[which(rownames(mat)%in%(unique(ttg$ens_gene[which(ttg$ext_gene=="LRP6")]))),])
goi<-as.data.frame(t(goi[1,]))
goi$sample_ID<-rownames(goi)
colnames(goi)[1]<-"scaled_reads_per_base"
plt<-merge(sampleinfo,goi, by.x="sample",by.y="sample_ID")
ggplot(plt, aes(sample.site, scaled_reads_per_base))+geom_boxplot()+geom_point(aes(fill=sample.site), shape=21, color="black")+theme_bw()+th+
  scale_fill_manual(values=c("#a6d96a","cornflowerblue"), name="Sample site")#SC TI
ggsave(here("../ibd/figs/single_gene_IBD/jpeg","LRP6_differenital_expression_SCTI_passageUTUD.jpeg"), width = 3, height = 2.5)


