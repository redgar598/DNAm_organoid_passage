


#'---
#'title: Differential expression gene levl sample site split
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
  library(cowplot)
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



#############################
#'# run sleuth on differentiation
#############################
diff_genes<-c("ASCL2","MKI67","LGR5","CA2","OLFM4","LYZ","MUC2","MUC1","CYP3A4","PLA2G2A","FABP1","KRT19","HELLS","SPINK4","FCGBP","NEAT1","TOP2A","DEFA5","DEFA6")#"TTF3",

## Fold change

diff_function<-function(condtion, segment){
  
  sampleinfo_segment<-sampleinfo[which(sampleinfo$sample.site==segment),]
  sampleinfo_diff<-sampleinfo_segment[which(sampleinfo_segment$passage_hilo==condtion & sampleinfo_segment$comparison=="differentiation"),]
  
  so_diff <- sleuth_prep(sampleinfo_diff, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)
  
  so_diff <- sleuth_fit(so_diff, ~differentiation + individual, 'full')
  so_diff <- sleuth_fit(so_diff, ~1 + individual, 'reduced')
  so_diff <- sleuth_lrt(so_diff, 'reduced', 'full')
  
  sleuth_table <- sleuth_results(so_diff, 'reduced:full', 'lrt', show_all = FALSE)
  
  so_diff <- sleuth_wt(so_diff, paste0('differentiationUD'))
  sleuth_table_wt <- sleuth_results(so_diff, 'differentiationUD', 'wt', show_all = FALSE)
  #b (Wald only): 'beta' value (effect size). Technically a biased estimator of the fold change. Only seen with Wald test results.
  sleuth_table<-merge(sleuth_table, sleuth_table_wt[,c("target_id","b")], by="target_id")
  sleuth_table
}

sleuth_table_low_TI<-diff_function("low","TI")
sleuth_table_low_SC<-diff_function("low","SC")
sleuth_table_high_TI<-diff_function("high","TI")
sleuth_table_high_SC<-diff_function("high","SC")

sleuth_significant_low_TI <- dplyr::filter(sleuth_table_low_TI, qval <= 0.05)
sleuth_significant_low_SC <- dplyr::filter(sleuth_table_low_SC, qval <= 0.05)
sleuth_significant_high_TI <- dplyr::filter(sleuth_table_high_TI, qval <= 0.05)
sleuth_significant_high_SC <- dplyr::filter(sleuth_table_high_SC, qval <= 0.05)

#key genes

# TI low
sleuth_significant_low_TI[which(sleuth_significant_low_TI$ext_gene%in%c("FABP6","MUC5B")),]
#SC low
sleuth_significant_low_SC[which(sleuth_significant_low_SC$ext_gene%in%c("FABP6","MUC5B")),]
#TI high
sleuth_significant_high_TI[which(sleuth_significant_high_TI$ext_gene%in%c("FABP6","MUC5B")),]
#SC high
sleuth_significant_high_SC[which(sleuth_significant_high_SC$ext_gene%in%c("FABP6","MUC5B")),]

# TI low
sleuth_table_low_TI[which(sleuth_table_low_TI$ext_gene%in%c("FABP6","MUC5B")),]
#SC low
sleuth_table_low_SC[which(sleuth_table_low_SC$ext_gene%in%c("FABP6","MUC5B")),]
#TI high
sleuth_table_high_TI[which(sleuth_table_high_TI$ext_gene%in%c("FABP6","MUC5B")),]
#SC high
sleuth_table_high_SC[which(sleuth_table_high_SC$ext_gene%in%c("FABP6","MUC5B")),]


# gene level sumary
sampleinfo_diff<-sampleinfo[which(sampleinfo$comparison=="differentiation"),]
so <- sleuth_prep(sampleinfo_diff, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)
mat <- sleuth:::spread_abundance_by(so$obs_norm, "scaled_reads_per_base",  so$sample_to_covariates$sample)

gene_exp_plot_differentiation<-function(gene, pltrows, dupens_rm, ordered){
  gene_ID<-ttg[which(ttg$ext_gene%in%gene),2:3]
  gene_ID<-gene_ID[!duplicated(gene_ID),]
  
  if(missing(dupens_rm)){}else{if(dupens_rm==T){
    dup_genes<-names(table(gene_ID$ext_gene))[which(table(gene_ID$ext_gene)>1)]
    pval_dup<-sleuth_table_low[which(sleuth_table_low$ext_gene%in%dup_genes),]
    pval_dup<-pval_dup %>% group_by(ext_gene) %>% slice(which.min(pval))
    gene_ID<-gene_ID[which(!gene_ID$ext_gene%in%dup_genes | gene_ID$ens_gene%in%pval_dup$target_id),]
  }}
  
  label<-paste(gene_ID$ext_gene,"\n(",gene_ID$ens_gene,")",sep="")
  
  goi<-as.data.frame(mat[which(rownames(mat)%in%(unique(ttg$ens_gene[which(ttg$ext_gene==gene)]))),])
  goi$sample_ID<-rownames(goi)
  colnames(goi)[1]<-"scaled_reads_per_base"
  plt<-merge(sampleinfo,goi, by.x="sample",by.y="sample_ID")
  plt$passage_hilo<-factor(plt$passage_hilo, levels=c("low","high"))
  plt$differentiation<-factor(plt$differentiation, levels=c("UD","D"))
  
  ggplot(plt, aes(sample.site, scaled_reads_per_base,fill=differentiation))+geom_boxplot()+
    geom_point(aes(fill=differentiation), shape=21, color="black", position=position_dodge(width=0.75))+theme_bw()+th_present+facet_wrap(~passage_hilo)+
    scale_fill_manual(values=c("#abd9e9","#a1d99b"), name="Differentiation")+ggtitle(label)+
    xlab("Scaled reads per base")+ylab("Segment")}

gene_exp_plot_differentiation("MUC5B")
gene_exp_plot_differentiation("FABP6")

gene_exp_plot_differentiation("NR1H4")
gene_exp_plot_differentiation("TFF3")

gene_exp_plot_differentiation("LGR5")
gene_exp_plot_differentiation("FABP1")

legend <- get_legend(
  # create some space to the left of the legend
  gene_exp_plot_differentiation("FABP1") + theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot_grid(gene_exp_plot_differentiation("LGR5")+ theme(legend.position="none"),
          gene_exp_plot_differentiation("FABP1")+ theme(legend.position="none"), legend,
          gene_exp_plot_differentiation("MUC5B")+ theme(legend.position="none"),
          gene_exp_plot_differentiation("FABP6")+ theme(legend.position="none"),
          ncol=3, rel_widths=c(1,1,0.5))

ggsave(here("figs","differenital_expression_differentiation_genes_segment_specific.pdf"),width = 7, height = 6)
ggsave(here("figs/jpeg","differenital_expression_differentiation_genes_segment_specific.jpeg"), width = 7, height = 6)


########################################## interaction
#'---
#'title: Differential expression gene levl sample site split
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



#############################
#'# run sleuth on differentiation
#############################
diff_genes<-c("ASCL2","MKI67","LGR5","CA2","OLFM4","LYZ","MUC2","MUC1","CYP3A4","PLA2G2A","FABP1","KRT19","HELLS","SPINK4","FCGBP","NEAT1","TOP2A","DEFA5","DEFA6")#"TTF3",

#'## low
sampleinfo_diff_low<-sampleinfo[which(sampleinfo$passage_hilo=="low" & sampleinfo$comparison=="differentiation"),]
table(sampleinfo_diff_low$differentiation, sampleinfo_diff_low$sample.site)

so_diff_low <- sleuth_prep(sampleinfo_diff_low, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)

# so_diff_low <- sleuth_fit(so_diff_low, ~differentiation + sample.site + individual, 'full')
# so_diff_low <- sleuth_fit(so_diff_low, ~1 + individual, 'reduced')
# so_diff_low <- sleuth_lrt(so_diff_low, 'reduced', 'full')
# models(so_diff_low)
# tests(so_diff_low)

#Null is that the differentiation effect is the same in both segments 
so_diff_low <- sleuth_fit(so_diff_low, ~differentiation + sample.site + differentiation:sample.site + individual, 'full')
so_diff_low <- sleuth_fit(so_diff_low, ~1 + differentiation + sample.site + individual, 'reduced')
so_diff_low <- sleuth_lrt(so_diff_low, 'reduced', 'full')
models(so_diff_low)
tests(so_diff_low)


sleuth_table_low <- sleuth_results(so_diff_low, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_low <- dplyr::filter(sleuth_table_low, qval <= 0.05)
head(sleuth_significant_low, 20)

#key genes
diff_markers_low<-sleuth_significant_low[which(sleuth_significant_low$ext_gene%in%diff_genes),]
sleuth_significant_low[which(sleuth_significant_low$ext_gene%in%c("FABP6","MUC5B")),]
sleuth_table_low[which(sleuth_table_low$ext_gene%in%c("FABP6","MUC5B")),]



#'## high
sampleinfo_diff_high<-sampleinfo[which(sampleinfo$passage_hilo=="high" & sampleinfo$comparison=="differentiation"),]
table(sampleinfo_diff_high$differentiation, sampleinfo_diff_high$sample.site)

so_diff_high <- sleuth_prep(sampleinfo_diff_high, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)

#Null is that the differentiation effect is the same in both segments 
so_diff_high <- sleuth_fit(so_diff_high, ~differentiation + sample.site + differentiation:sample.site + individual, 'full')
so_diff_high <- sleuth_fit(so_diff_high, ~1 + differentiation + sample.site + individual, 'reduced')
so_diff_high <- sleuth_lrt(so_diff_high, 'reduced', 'full')
models(so_diff_high)
tests(so_diff_high)

sleuth_table_high <- sleuth_results(so_diff_high, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_high <- dplyr::filter(sleuth_table_high, qval <= 0.05)
head(sleuth_significant_high, 20)

#key genes
diff_markers_high<-sleuth_significant_high[which(sleuth_significant_high$ext_gene%in%diff_genes),]
sleuth_significant_high[which(sleuth_significant_high$ext_gene%in%c("FABP6","MUC5B")),]
sleuth_table_high[which(sleuth_table_high$ext_gene%in%c("FABP6","MUC5B")),]

## low not high
diff_markers_low[which(!(diff_markers_low$ext_gene%in%diff_markers_high$ext_gene)),]
diff_markers_low[which(!(diff_markers_low$target_id%in%diff_markers_high$target_id)),]


# gene level sumary
mat_high <- sleuth:::spread_abundance_by(so_diff_high$obs_norm, "scaled_reads_per_base",  so_diff_high$sample_to_covariates$sample)
mat_low <- sleuth:::spread_abundance_by(so_diff_low$obs_norm, "scaled_reads_per_base",  so_diff_low$sample_to_covariates$sample)

gene_exp_plot_differentiation<-function(gene, pltrows, dupens_rm, ordered){
  if(length(gene)==1){
    goi_low<-as.data.frame(mat_low[which(rownames(mat_low)%in%(unique(ttg$ens_gene[which(ttg$ext_gene==gene)]))),])
    goi_high<-as.data.frame(mat_high[which(rownames(mat_high)%in%(unique(ttg$ens_gene[which(ttg$ext_gene==gene)]))),])
    goi_low$passage<-"low"
    goi_low$sample_ID<-rownames(goi_low)
    colnames(goi_low)[1]<-"scaled_reads_per_base"
    goi_high$passage<-"high"
    goi_high$sample_ID<-rownames(goi_high)
    colnames(goi_high)[1]<-"scaled_reads_per_base"
    goi<-rbind(goi_low, goi_high)
    plt<-merge(sampleinfo,goi, by.x="sample",by.y="sample_ID")
    plt$passage_hilo<-factor(plt$passage_hilo, levels=c("low","high"))
    plt$differentiation<-factor(plt$differentiation, levels=c("UD","D"))
    ggplot(plt, aes(differentiation, scaled_reads_per_base,fill=differentiation))+geom_boxplot()+
      geom_point(aes(fill=differentiation), shape=21, color="black")+theme_bw()+th+facet_grid(sample.site~passage_hilo)+
      scale_fill_manual(values=c("#abd9e9","#a1d99b"), name="Differentiation")
  }else{
    goi_low<-as.data.frame(mat_low[which(rownames(mat_low)%in%(unique(ttg$ens_gene[which(ttg$ext_gene%in%gene)]))),])
    goi_high<-as.data.frame(mat_high[which(rownames(mat_high)%in%(unique(ttg$ens_gene[which(ttg$ext_gene%in%gene)]))),])
    gene_ID<-ttg[which(ttg$ext_gene%in%gene),2:3]
    gene_ID<-gene_ID[!duplicated(gene_ID),]
    
    if(missing(dupens_rm)){}else{if(dupens_rm==T){
      dup_genes<-names(table(gene_ID$ext_gene))[which(table(gene_ID$ext_gene)>1)]
      pval_dup<-sleuth_table_low[which(sleuth_table_low$ext_gene%in%dup_genes),]
      pval_dup<-pval_dup %>% group_by(ext_gene) %>% slice(which.min(pval))
      gene_ID<-gene_ID[which(!gene_ID$ext_gene%in%dup_genes | gene_ID$ens_gene%in%pval_dup$target_id),]
    }}
    
    gene_ID$label<-paste(gene_ID$ext_gene,"\n(",gene_ID$ens_gene,")",sep="")
    goi_high$gene_ID<-rownames(goi_high)
    goi_low$gene_ID<-rownames(goi_low)
    goi_high<-melt(goi_high)
    goi_low<-melt(goi_low)
    goi_low$passage<-"low"
    goi_high$passage<-"high"
    
    goi<-rbind(goi_low, goi_high)
    goi<-merge(goi,gene_ID, by.x="gene_ID", by.y="ens_gene")
    plt<-merge(sampleinfo,goi, by.x="sample",by.y="variable")
    plt$passage_hilo<-factor(plt$passage_hilo, levels=c("low","high"))
    plt$differentiation<-factor(plt$differentiation, levels=c("UD","D"))
    
    if(missing(ordered)){}else{if(ordered==T){
      uni_gene<-plt[!duplicated(plt[,c("ext_gene","label")]),]
      plt$label<-factor(plt$label, levels=uni_gene$label[match(gene,uni_gene$ext_gene)])
    }}
    
    ggplot(plt, aes(passage_hilo, value,fill=differentiation))+geom_boxplot(outlier.shape=NA)+
      geom_point(aes(fill=differentiation), shape=21, color="black",position = position_dodge(width=0.75))+theme_bw()+th+
      facet_wrap(~label, scales="free_y", nrow=pltrows)+ylab("Scaled reads per base")+xlab("Passage")+
      scale_fill_manual(values=c("#abd9e9","#a1d99b"), name="Differentiation")
  }}
gene_exp_plot_differentiation(c("LGR5","ASCL2","HELLS","FABP1","KRT19","MUC2"),2,T,T)

