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
      geom_point(aes(fill=differentiation), shape=21, color="black", position=position_dodge(width=0.75))+theme_bw()+th+facet_wrap(~passage_hilo)+
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

ggsave(here("figs","differenital_expression_differentiation_genes_segment_specific.pdf"),width = 10, height = 6)
ggsave(here("figs/jpeg","differenital_expression_differentiation_genes_segment_specific.jpeg"), width = 10, height = 6)



#   }else{
#     goi<-as.data.frame(mat[which(rownames(mat)%in%(unique(ttg$ens_gene[which(ttg$ext_gene%in%gene)]))),])
#     gene_ID<-ttg[which(ttg$ext_gene%in%gene),2:3]
#     gene_ID<-gene_ID[!duplicated(gene_ID),]
# 
#     if(missing(dupens_rm)){}else{if(dupens_rm==T){
#       dup_genes<-names(table(gene_ID$ext_gene))[which(table(gene_ID$ext_gene)>1)]
#       pval_dup<-sleuth_table_low[which(sleuth_table_low$ext_gene%in%dup_genes),]
#       pval_dup<-pval_dup %>% group_by(ext_gene) %>% slice(which.min(pval))
#       gene_ID<-gene_ID[which(!gene_ID$ext_gene%in%dup_genes | gene_ID$ens_gene%in%pval_dup$target_id),]
#     }}
# 
#     gene_ID$label<-paste(gene_ID$ext_gene,"\n(",gene_ID$ens_gene,")",sep="")
#     goi$gene_ID<-rownames(goi)
#     goi<-melt(goi)
# 
#     goi<-merge(goi,gene_ID, by.x="gene_ID", by.y="ens_gene")
#     plt<-merge(sampleinfo,goi, by.x="sample",by.y="variable")
#     plt$passage_hilo<-factor(plt$passage_hilo, levels=c("low","high"))
#     plt$differentiation<-factor(plt$differentiation, levels=c("UD","D"))
# 
#     if(missing(ordered)){}else{if(ordered==T){
#       uni_gene<-plt[!duplicated(plt[,c("ext_gene","label")]),]
#       plt$label<-factor(plt$label, levels=uni_gene$label[match(gene,uni_gene$ext_gene)])
#     }}
# 
#     ggplot(plt, aes(passage_hilo, value,fill=differentiation))+geom_boxplot(outlier.shape=NA)+
#       geom_point(aes(fill=differentiation), shape=21, color="black",position = position_dodge(width=0.75))+theme_bw()+th+
#       facet_wrap(~label, scales="free_y", nrow=pltrows)+ylab("Scaled reads per base")+xlab("Passage")+
#       scale_fill_manual(values=c("#abd9e9","#a1d99b"), name="Differentiation")
#   }}
# 
# ggplot(plt, aes(sample.site, scaled_reads_per_base,fill=differentiation))+geom_boxplot()+
#   geom_point(aes(fill=differentiation), shape=21, color="black", position=position_dodge(width=0.75))+theme_bw()+th+facet_wrap(~passage_hilo)+
#   scale_fill_manual(values=c("#abd9e9","#a1d99b"), name="Differentiation")}
# 
# 
# gene_exp_plot_differentiation(c("LGR5","ASCL2","HELLS","FABP1","KRT19","MUC2"),2,T,T)

gene_exp_plot_differentiation("LYZ",1)
gene_exp_plot_differentiation(diff_genes,4)
gene_exp_plot_differentiation(c("LGR5","ASCL2","HELLS","FABP1","KRT19","MUC2"),2,T,T)
ggsave(here("figs","differenital_expression_differentiation_genes.pdf"),width = 10, height = 5)
ggsave(here("figs/jpeg","differenital_expression_differentiation_genes.jpeg"), width = 10, height = 5)




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

gene_exp_plot_differentiation("LYZ",1)
gene_exp_plot_differentiation(diff_genes,4)
gene_exp_plot_differentiation(c("LGR5","ASCL2","HELLS","FABP1","KRT19","MUC2"),2,T,T)
ggsave(here("figs","differenital_expression_differentiation_genes.pdf"),width = 10, height = 5)
ggsave(here("figs/jpeg","differenital_expression_differentiation_genes.jpeg"), width = 10, height = 5)


gene_exp_plot_differentiation("ASCL2",1)
gene_exp_plot_differentiation("FABP1",1)
gene_exp_plot_differentiation("KRT19",1)
gene_exp_plot_differentiation("TOP2A",1)

## differential between passage and DNAm diff with passage
# but actually tested in high
length(unique(sleuth_significant_high$ext_gene))
length(unique(sleuth_significant_low$ext_gene))

sleuth_significant_low_tested<-sleuth_significant_low[which(sleuth_significant_low$target_id%in%sleuth_table_high$target_id),]
sleuth_significant_high_tested<-sleuth_significant_high[which(sleuth_significant_high$target_id%in%sleuth_table_low$target_id),]

length(unique(sleuth_significant_high_tested$ext_gene))
length(unique(sleuth_significant_low_tested$ext_gene))

length(unique(sleuth_table_high$ext_gene))
length(unique(sleuth_table_low$ext_gene))

low_not_high_diff<-sleuth_significant_low_tested[which(!(sleuth_significant_low_tested$ext_gene%in%sleuth_significant_high_tested$ext_gene)),]
length(unique(low_not_high_diff$ext_gene))


#' ## DNAm passage genes original 80 and validation
diff_genes_db_hypovalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypomethylation_UTUD.txt"))
diff_genes_db_hypervalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypermethylation_UTUD.txt"))
sleuth_sig_DNAm<-low_not_high_diff[which(low_not_high_diff$ext_gene%in%c(diff_genes_db_hypovalidation_original$V1, diff_genes_db_hypervalidation_original$V1)),]
length(unique(sleuth_sig_DNAm$ext_gene))


gene_exp_plot_differentiation(c("DISC1","BCL3","NR1D1","LYZ"))


write.table(unique(low_not_high_diff$ext_gene), file=here("data/validation","low_not_high_differentiationgenes.txt"), quote=F, row.names = F, col.names = F)

low_not_high_diff[which(low_not_high_diff$ext_gene%in%diff_genes),]
low_not_high_diff[grep("DNMT|TET|TLR",low_not_high_diff$ext_gene),]


gene_exp_plot_differentiation("HSPA1A",1)
gene_exp_plot_differentiation("LYZ",1)
### different in low and high LYZ 

## interesting for GO enrichment
gene_exp_plot_differentiation("FUS",1)


#'### Gene Table
load(here("data/validation/CpG_validated.RData"))
EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194


table_low_not_high<-function(genes, stats_table_low, stats_table_high){
  stats_low<-stats_table_low[which(stats_table_low$ext_gene%in%genes),c("ext_gene","target_id","pval","qval")]
  stats_high<-stats_table_high[which(stats_table_high$ext_gene%in%genes),c("ext_gene","target_id","pval","qval")]
  
  EPIC_genes_CpGs<-EPIC_genes[which(EPIC_genes$Gene.name%in%genes),]
  
  EPIC_genes_CpGs_hypo<-EPIC_genes_CpGs[which(EPIC_genes_CpGs$IlmnID%in%CpG_hypo_validated),]
  EPIC_genes_CpGs_hypo<-as.data.frame(EPIC_genes_CpGs_hypo %>% 
                                        group_by(Gene.name) %>% 
                                        mutate(hypo_DNAm_CpG = paste0(unique(IlmnID), collapse = ", ")) )
  EPIC_genes_CpGs_hypo<-EPIC_genes_CpGs_hypo[,c("Gene.name","hypo_DNAm_CpG")]
  EPIC_genes_CpGs_hypo<-EPIC_genes_CpGs_hypo[!duplicated(EPIC_genes_CpGs_hypo),]
  EPIC_genes_CpGs_hypo
  
  EPIC_genes_CpGs_hyper<-EPIC_genes_CpGs[which(EPIC_genes_CpGs$IlmnID%in%CpG_hyper_validated),]
  EPIC_genes_CpGs_hyper<-as.data.frame(EPIC_genes_CpGs_hyper %>% 
                                         group_by(Gene.name) %>% 
                                         mutate(hyper_DNAm_CpG = paste0(unique(IlmnID), collapse = ", ")) )
  EPIC_genes_CpGs_hyper<-EPIC_genes_CpGs_hyper[,c("Gene.name","hyper_DNAm_CpG")]
  EPIC_genes_CpGs_hyper<-EPIC_genes_CpGs_hyper[!duplicated(EPIC_genes_CpGs_hyper),]
  EPIC_genes_CpGs_hyper
  
  EPIC_genes_CpGs<-merge(EPIC_genes_CpGs_hyper, EPIC_genes_CpGs_hypo, by="Gene.name", all=T)
  
  colnames(stats_low)<-c("Gene","Ensembl gene name","Low Passage p value","Low Passage q value")
  colnames(stats_high)<-c("Gene","Ensembl gene name","High Passage p value","High Passage q value")
  
  stats<-merge(stats_low, stats_high, by=c("Gene","Ensembl gene name"))
  stats_cpgs<-merge(stats, EPIC_genes_CpGs, by.x="Gene", by.y="Gene.name")
  colnames(stats_cpgs)<-c("Gene","Ensembl gene name", "Low Passage p value","Low Passage q value",
                          "High Passage p value", "High Passage q value","Hypermethylated CpGs", "Hypomethylated CpGs")
  stats_cpgs$`High Passage q value`<-signif(stats_cpgs$`High Passage q value`, 2)
  stats_cpgs$`High Passage p value`<-signif(stats_cpgs$`High Passage p value`, 2)
  stats_cpgs$`Low Passage q value`<-signif(stats_cpgs$`Low Passage q value`, 2)
  stats_cpgs$`Low Passage p value`<-signif(stats_cpgs$`Low Passage p value`, 2)
  stats_cpgs[order(stats_cpgs$`Low Passage p value`),]}

stats_cpgs_diff<-table_low_not_high(unique(sleuth_sig_DNAm$ext_gene), sleuth_sig_DNAm,sleuth_table_high)
head(stats_cpgs_diff)
write.csv(stats_cpgs_diff, file="reports_figs/validation_update/Supplement_Table5_differentiation_genes.csv", row.names = F)


########
#'# wnt genes
########
wnt_genes<-read.table(file="/home/redgar/Documents/ibd/data/GO_0060070.txt", header=T, sep="\t")

wnt_high<-sleuth_significant_high[which(sleuth_significant_high$ext_gene%in%wnt_genes$Element),]
wnt_low<-sleuth_significant_low[which(sleuth_significant_low$ext_gene%in%wnt_genes$Element),]

wnt_low[which(!(wnt_low$ext_gene%in%wnt_high$ext_gene)),]


gene_exp_plot_differentiation(wnt_lowsleuth_sig_DNAm[which(!(wnt_low$ext_gene%in%wnt_high$ext_gene)),]$ext_gene,2)

sleuth_sig_DNAm[which(sleuth_sig_DNAm$ext_gene%in%wnt_genes$Element),]
gene_exp_plot_differentiation(c("DISC1","WNT11"),1)

#############
#'## volcano
#############
source(here("scripts/validation/rnaseq/00_volcano.R"))

so_diff_high <- sleuth_wt(so_diff_high, paste0('differentiationUD'))
sleuth_table_high_wt <- sleuth_results(so_diff_high, 'differentiationUD', 'wt', show_all = FALSE)

vol_high<-makeVolcano(sleuth_table_high$pval, sleuth_table_high_wt$b, 2,0.015, "Difference Between\nDifferentiated and \nUndifferentiated", 6,13,11000)


so_diff_low <- sleuth_wt(so_diff_low, paste0('differentiationUD'))
sleuth_table_low_wt <- sleuth_results(so_diff_low, 'differentiationUD', 'wt', show_all = FALSE)

vol_low<-makeVolcano(sleuth_table_low$pval, sleuth_table_low_wt$b, 2,0.015, "Difference Between\nDifferentiated and \nUndifferentiated", 6,13,11000)

grid.arrange(vol_low, vol_high, ncol=2)
ggsave("figs/jpeg/RNAseq_volcano_differentiation.jpeg",grid.arrange(vol_low, vol_high, ncol=2),  width = 15, height = 6)


###############
#### Genes with no induction in high (or low)
###############
#ie not significant in high but also low FC

low_not_high_highFC<-sleuth_table_high_wt[which(sleuth_table_high_wt$target_id%in%low_not_high_diff$target_id),]
head(low_not_high_highFC[order(abs(low_not_high_highFC$b)),])

# top low p
not_induced<-low_not_high_highFC[which(abs(low_not_high_highFC$b)<0.01),]
nrow(not_induced)
low_not_high_top<-low_not_high_diff[order(low_not_high_diff$pval),][1:1000,]
intersect(low_not_high_top$ext_gene, not_induced$ext_gene)
gene_exp_plot_differentiation(intersect(low_not_high_top$ext_gene, not_induced$ext_gene))

# top low b
not_induced<-low_not_high_highFC[order(abs(low_not_high_highFC$b)),][1:200,]
sleuth_table_low_wt_sig<-sleuth_table_low_wt[which(sleuth_table_low_wt$target_id%in%low_not_high_diff$target_id),]
low_not_high_topb<-sleuth_table_low_wt_sig[rev(order(abs(sleuth_table_low_wt_sig$b))),][1:1000,]
intersect(low_not_high_topb$ext_gene, not_induced$ext_gene)
gene_exp_plot_differentiation(intersect(low_not_high_topb$ext_gene, not_induced$ext_gene))

gene_exp_plot_differentiation(c("RC3H2","TAOK1"))

# top low b (negative low b)
not_induced<-low_not_high_highFC[order(abs(low_not_high_highFC$b)),][1:100,]
sleuth_table_low_wt_sig<-sleuth_table_low_wt[which(sleuth_table_low_wt$target_id%in%low_not_high_diff$target_id),]
low_not_high_topb<-sleuth_table_low_wt_sig[which(sleuth_table_low_wt_sig$b<0),]
intersect(low_not_high_topb$ext_gene, not_induced$ext_gene)
gene_exp_plot_differentiation(intersect(low_not_high_topb$ext_gene, not_induced$ext_gene))

gene_exp_plot_differentiation(c("RC3H2","TAOK1","OR10AH1P","PEAR1","SGMS1","ZNF321P","VEGFB","STK39"),4)

gene_exp_plot_differentiation(c("RC3H2","PEAR1","ZNF321P","VEGFB"),2)

# but also diff DNAm with passage
intersect(intersect(low_not_high_topb$ext_gene, not_induced$ext_gene),sleuth_sig_DNAm$ext_gene)
gene_exp_plot_differentiation(intersect(intersect(low_not_high_topb$ext_gene, not_induced$ext_gene),sleuth_sig_DNAm$ext_gene))

## representaive genes
gene_exp_plot_differentiation(unique(sleuth_sig_DNAm$ext_gene)[151:175],5)
gene_exp_plot_differentiation(c("SP3","RAB3B","DISC1","STK39","FAM13A","ZNF702P","BEGAIN","COL1A1","THRB","FKBP5", "WNT11"),4)



gene_exp_plot_differentiation(c("FKBP5","STK39","COL1A1","WNT11"),1)
ggsave(here("figs","differenital_expression_differentiation_genes_difflowhigh.pdf"),width = 10, height = 2.5)
ggsave(here("figs/jpeg","differenital_expression_differentiation_genes_difflowhigh.jpeg"), width = 10, height = 2.5)



#' 
#' 
#' #############################
#' #'# run sleuth on differentiation interaction terms
#' #############################
#' diff_genes<-c("ASCL2","MKI67","LGR5","CA2","OLFM4","LYZ","MUC2","MUC1","CYP3A4","PLA2G2A","FABP1","KRT19","HELLS","SPINK4","FCGBP","NEAT1","TOP2A")#"TTF3",
#' 
#' #'## low
#' sampleinfo_diff<-sampleinfo[which(sampleinfo$comparison=="differentiation"),]
#' table(sampleinfo_diff$differentiation)
#' 
#' so_diff <- sleuth_prep(sampleinfo_diff, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)
#' 
#' so_diff <- sleuth_fit(so_diff, ~differentiation + passage_hilo + differentiation:passage_hilo +individual, 'full')
#' so_diff <- sleuth_fit(so_diff, ~differentiation + passage_hilo +individual, 'reduced')
#' so_diff <- sleuth_lrt(so_diff, 'reduced', 'full')
#' models(so_diff)
#' tests(so_diff)
#' 
#' sleuth_table <- sleuth_results(so_diff, 'reduced:full', 'lrt', show_all = FALSE)
#' sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
#' head(sleuth_significant, 20)
#' 
#' #sleuth_table[which(sleuth_table$ext_gene%in%diff_genes),]
#' sleuth_significant[which(sleuth_significant$ext_gene%in%diff_genes),]
#' 
#' sleuth_significant <- dplyr::filter(sleuth_table, pval <= 0.01)
#' head(sleuth_significant, 20)
#' 
#' 
#' gene_exp_plot_differentiation("SEMA3B")
#' 





###############################################################################################################################################
#############################
#'# run sleuth on IFNg
#############################
MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","MR1","CD1D","IRF1","NLRC5")

#'## low
sampleinfo_IFNg_low<-sampleinfo[which(sampleinfo$passage_hilo=="low" &  sampleinfo$comparison=="cytokine" & sampleinfo$treatment!="TNFa"),]
table(sampleinfo_IFNg_low$treatment)

so_IFNg_low <- sleuth_prep(sampleinfo_IFNg_low, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)

so_IFNg_low <- sleuth_fit(so_IFNg_low, ~treatment + individual, 'full')
so_IFNg_low <- sleuth_fit(so_IFNg_low, ~1 + individual, 'reduced')
so_IFNg_low <- sleuth_lrt(so_IFNg_low, 'reduced', 'full')
models(so_IFNg_low)
tests(so_IFNg_low)

sleuth_table_low <- sleuth_results(so_IFNg_low, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_low <- dplyr::filter(sleuth_table_low, qval <= 0.05)
head(sleuth_significant_low, 20)

#sleuth_table_low[which(sleuth_table_low$ext_gene%in%diff_genes),]
sleuth_significant_low[which(sleuth_significant_low$ext_gene%in%diff_genes),]
MHCI_low<-sleuth_significant_low[which(sleuth_significant_low$ext_gene%in%MHCI),]



#'## high
sampleinfo_IFNg_high<-sampleinfo[which(sampleinfo$passage_hilo=="high" &  sampleinfo$comparison=="cytokine" & sampleinfo$treatment!="TNFa"),]
table(sampleinfo_IFNg_high$treatment)

so_IFNg_high <- sleuth_prep(sampleinfo_IFNg_high, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)

so_IFNg_high <- sleuth_fit(so_IFNg_high, ~treatment + individual, 'full')
so_IFNg_high <- sleuth_fit(so_IFNg_high, ~1 + individual, 'reduced')
so_IFNg_high <- sleuth_lrt(so_IFNg_high, 'reduced', 'full')
models(so_IFNg_high)
tests(so_IFNg_high)

sleuth_table_high <- sleuth_results(so_IFNg_high, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_high <- dplyr::filter(sleuth_table_high, qval <= 0.05)
head(sleuth_significant_high, 20)

#sleuth_table_high[which(sleuth_table_high$ext_gene%in%diff_genes),]
sleuth_significant_high[which(sleuth_significant_high$ext_gene%in%diff_genes),]
MHCI_high<-sleuth_significant_high[which(sleuth_significant_high$ext_gene%in%MHCI),]


#low not high
MHCI_low[which(!(MHCI_low$ext_gene%in%MHCI_high$ext_gene)),]


# gene level sumary
mat_high <- sleuth:::spread_abundance_by(so_IFNg_high$obs_norm, "scaled_reads_per_base",  so_IFNg_high$sample_to_covariates$sample)
mat_low <- sleuth:::spread_abundance_by(so_IFNg_low$obs_norm, "scaled_reads_per_base",  so_IFNg_low$sample_to_covariates$sample)

gene_exp_plot_treatment<-function(gene, sig){
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
    plt$treatment<-factor(plt$treatment, levels=c("UT","IFNg"))
    ggplot(plt, aes(treatment, scaled_reads_per_base,fill=treatment))+geom_boxplot()+
      geom_point(aes(fill=treatment), shape=21, color="black")+theme_bw()+th+facet_wrap(~passage_hilo)+
      scale_fill_manual(values=c("grey80","cornflowerblue"), name="Treatment")
  }else{
    goi_low<-as.data.frame(mat_low[which(rownames(mat_low)%in%(unique(ttg$ens_gene[which(ttg$ext_gene%in%gene)]))),])
    goi_high<-as.data.frame(mat_high[which(rownames(mat_high)%in%(unique(ttg$ens_gene[which(ttg$ext_gene%in%gene)]))),])
    gene_ID<-ttg[which(ttg$ext_gene%in%gene),2:3]
    if(missing(sig)){}else{
      gene_ID<-gene_ID[which(gene_ID$ens_gene%in%sig),]
    }
    gene_ID<-gene_ID[!duplicated(gene_ID),]
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
    plt$treatment<-factor(plt$treatment, levels=c("UT","IFNg"))
    ggplot(plt, aes(passage_hilo, value,fill=treatment))+geom_boxplot(outlier.shape=NA)+
      geom_point(aes(fill=treatment), shape=21, color="black",position = position_dodge(width=0.75))+theme_bw()+th+
      facet_wrap(~label, scales="free_y", nrow=1)+ylab("Scaled reads per base")+xlab("Passage")+
      scale_fill_manual(values=c("grey80","cornflowerblue"), name="Treatment")
  }}

gene_exp_plot_treatment(MHCI)
gene_exp_plot_treatment(c("NLRC5","TAP1"))

top_sig<-sleuth_significant_low %>% group_by(ext_gene) %>% slice(which.min(pval))
gene_exp_plot_treatment(MHCI,top_sig$target_id)
gene_exp_plot_treatment(c("NLRC5","TAP1","PSMB9","PSMB8","B2M","HLA-A"),top_sig$target_id)
ggsave(here("figs","differenital_expression_MHCI_genes_treatmentlowhigh.pdf"),width = 10, height = 5)
ggsave(here("figs/jpeg","differenital_expression_MHCI_genes_treatmentlowhigh.jpeg"), width = 10, height = 5)


## differential between passage and DNAm diff with passage
# but actually tested in high
length(unique(sleuth_significant_high$ext_gene))
length(unique(sleuth_significant_low$ext_gene))

sleuth_significant_low_tested<-sleuth_significant_low[which(sleuth_significant_low$target_id%in%sleuth_table_high$target_id),]
sleuth_significant_high_tested<-sleuth_significant_high[which(sleuth_significant_high$target_id%in%sleuth_table_low$target_id),]

length(unique(sleuth_significant_high_tested$ext_gene))
length(unique(sleuth_significant_low_tested$ext_gene))

length(unique(sleuth_table_high$ext_gene))
length(unique(sleuth_table_low$ext_gene))


low_not_high_treatment<-sleuth_significant_low_tested[which(!(sleuth_significant_low_tested$ext_gene%in%sleuth_significant_high_tested$ext_gene)),]
length(unique(low_not_high_treatment$ext_gene))

low_not_high_treatment[which(low_not_high_treatment$ext_gene%in%MHCI),]
low_not_high_treatment[grep("DNMT|TET|TLR",low_not_high_treatment$ext_gene),]

gene_exp_plot_treatment("NDUFA8")
gene_exp_plot_treatment("NDUFA7")
gene_exp_plot_treatment("SEC24B")

write.table(unique(low_not_high_treatment$ext_gene), file=here("data/validation","low_not_high_treatmentgenes.txt"), quote=F, row.names = F, col.names = F)




#' ## DNAm passage genes original 80 and validation
diff_genes_db_hypovalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypomethylation_UTUD.txt"))
diff_genes_db_hypervalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypermethylation_UTUD.txt"))
sleuth_sig_DNAm<-low_not_high_treatment[which(low_not_high_treatment$ext_gene%in%c(diff_genes_db_hypovalidation_original$V1, diff_genes_db_hypervalidation_original$V1)),]
length(unique(sleuth_sig_DNAm$ext_gene))

gene_exp_plot_treatment(c("IDH1","TGFB1","LYZ","MSH2"))

gene_exp_plot_treatment(unique(sleuth_sig_DNAm$ext_gene)[169:194])
gene_exp_plot_treatment(c("ARHGAP19","CADM1","CDKN1A","LNX1","MSH2","ACSL6","DPYD","NR1D1",
                          "RDX","ABCG1","PSTPIP2","ACVR1","FOXO1","CLK3","CLCN6"))

gene_exp_plot_treatment(c("LNX1","MSH2","CLK3","CADM1"))
ggsave(here("figs","differenital_expression_treatment_genes_treatmentlowhigh.pdf"),width = 10, height = 2.5)
ggsave(here("figs/jpeg","differenital_expression_treatment_genes_treatmentlowhigh.jpeg"), width = 10, height = 2.5)

#'### Gene Table
stats_cpgs_IFNg<-table_low_not_high(unique(sleuth_sig_DNAm$ext_gene), sleuth_sig_DNAm,sleuth_table_high)
head(stats_cpgs_IFNg)
write.csv(stats_cpgs_IFNg, file="reports_figs/validation_update/Supplement_Table6_IFNg_genes.csv", row.names = F)



######
#'## volcano
######
source(here("scripts/validation/rnaseq/00_volcano.R"))

so_IFNg_high <- sleuth_wt(so_IFNg_high, paste0('treatmentUT'))
sleuth_table_high_wt <- sleuth_results(so_IFNg_high, 'treatmentUT', 'wt', show_all = FALSE)

sleuth_table_high_wt<-sleuth_table_high_wt[match(sleuth_table_high$target_id,sleuth_table_high_wt$target_id),]
identical(sleuth_table_high_wt$target_id, sleuth_table_high$target_id)
vol_high<-makeVolcano(sleuth_table_high$pval, -sleuth_table_high_wt$b, 2,0.015, "Difference Between\nUntreated and \nTreated", 9,11, 6000)


so_IFNg_low <- sleuth_wt(so_IFNg_low, paste0('treatmentUT'))
sleuth_table_low_wt <- sleuth_results(so_IFNg_low, 'treatmentUT', 'wt', show_all = FALSE)


sleuth_table_low_wt<-sleuth_table_low_wt[match(sleuth_table_low$target_id,sleuth_table_low_wt$target_id),]
identical(sleuth_table_low_wt$target_id, sleuth_table_low$target_id)
vol_low<-makeVolcano(sleuth_table_low$pval, -sleuth_table_low_wt$b, 2,0.015, "Difference Between\nUntreated and \nTreated", 9,11,6000)

grid.arrange(vol_low, vol_high, ncol=2)
ggsave("figs/jpeg/RNAseq_volcano_treatment.jpeg",grid.arrange(vol_low, vol_high, ncol=2),  width = 15, height = 6)

gene_exp_plot_treatment("WARS1")


## low no high with fold change
low_fc<-sleuth_table_low_wt[which(abs(sleuth_table_low_wt$b)>2),]
sleuth_significant_low_fc<-sleuth_significant_low[which(sleuth_significant_low$target_id%in%low_fc$target_id),]
high_fc<-sleuth_table_high_wt[which(abs(sleuth_table_high_wt$b)>2),]
sleuth_significant_high_fc<-sleuth_significant_high[which(sleuth_significant_high$target_id%in%high_fc$target_id),]

low_not_high_treatment_fc<-sleuth_significant_low_fc[which(!(sleuth_significant_low_fc$ext_gene%in%sleuth_significant_high_fc$ext_gene)),]
gene_exp_plot_treatment("BTN3A3")
gene_exp_plot_treatment("JAK2")
gene_exp_plot_treatment("TRIM69")


#' ### heatmap
sampleinfo_IFNg<-sampleinfo[which(sampleinfo$comparison=="cytokine" & sampleinfo$treatment!="TNFa"),]
table(sampleinfo_IFNg$treatment)
so_IFNg <- sleuth_prep(sampleinfo_IFNg, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)

plot_transcript_heatmap(so_IFNg,transcripts =sleuth_significant_low$target_id[which(sleuth_significant_low$ext_gene%in%MHCI)] , 
                        annotation_cols = c("individual","passage_hilo", "treatment"))
dev.off()

het_low<-plot_transcript_heatmap(so_IFNg_low,transcripts =sleuth_significant_low$target_id[which(sleuth_significant_low$ext_gene%in%MHCI)] , 
                                 annotation_cols = c("individual","passage_hilo", "treatment"))
ggsave("figs/jpeg/RNAseq_heatmap_low_MHCI.jpeg",het_low,  width = 10, height = 10)
dev.off()

het_high<-plot_transcript_heatmap(so_IFNg_high,transcripts =sleuth_significant_low$target_id[which(sleuth_significant_low$ext_gene%in%MHCI)] , 
                                  annotation_cols = c("individual","passage_hilo", "treatment"))
ggsave("figs/jpeg/RNAseq_heatmap_high_MHCI.jpeg",het_high,  width = 10, height = 10)
dev.off()



###############################################################################################################################################
#############################
#'# run sleuth on TNFa
#############################

#'## low
sampleinfo_TNFa_low<-sampleinfo[which(sampleinfo$passage_hilo=="low" &  sampleinfo$comparison=="cytokine" & sampleinfo$treatment!="IFNg"),]
table(sampleinfo_TNFa_low$treatment)

so_TNFa_low <- sleuth_prep(sampleinfo_TNFa_low, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)

so_TNFa_low <- sleuth_fit(so_TNFa_low, ~treatment + individual, 'full')
so_TNFa_low <- sleuth_fit(so_TNFa_low, ~1 + individual, 'reduced')
so_TNFa_low <- sleuth_lrt(so_TNFa_low, 'reduced', 'full')
models(so_TNFa_low)
tests(so_TNFa_low)

sleuth_table_low <- sleuth_results(so_TNFa_low, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_low <- dplyr::filter(sleuth_table_low, qval <= 0.05)
head(sleuth_significant_low, 20)

#sleuth_table_low[which(sleuth_table_low$ext_gene%in%diff_genes),]
sleuth_significant_low[which(sleuth_significant_low$ext_gene%in%MHCI),]



#'## high
sampleinfo_TNFa_high<-sampleinfo[which(sampleinfo$passage_hilo=="high" &  sampleinfo$comparison=="cytokine" & sampleinfo$treatment!="IFNg"),]
table(sampleinfo_TNFa_high$treatment)

so_TNFa_high <- sleuth_prep(sampleinfo_TNFa_high, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)

so_TNFa_high <- sleuth_fit(so_TNFa_high, ~treatment + individual, 'full')
so_TNFa_high <- sleuth_fit(so_TNFa_high, ~1 + individual, 'reduced')
so_TNFa_high <- sleuth_lrt(so_TNFa_high, 'reduced', 'full')
models(so_TNFa_high)
tests(so_TNFa_high)

sleuth_table_high <- sleuth_results(so_TNFa_high, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_high <- dplyr::filter(sleuth_table_high, qval <= 0.05)
head(sleuth_significant_high, 20)

#sleuth_table_high[which(sleuth_table_high$ext_gene%in%diff_genes),]
sleuth_significant_high[which(sleuth_significant_high$ext_gene%in%MHCI),]



# gene level sumary
mat_high <- sleuth:::spread_abundance_by(so_TNFa_high$obs_norm, "scaled_reads_per_base",  so_TNFa_high$sample_to_covariates$sample)
mat_low <- sleuth:::spread_abundance_by(so_TNFa_low$obs_norm, "scaled_reads_per_base",  so_TNFa_low$sample_to_covariates$sample)

gene_exp_plot_treatment<-function(gene){
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
    plt$treatment<-factor(plt$treatment, levels=c("UT","TNFa"))
    ggplot(plt, aes(treatment, scaled_reads_per_base,fill=treatment))+geom_boxplot()+
      geom_point(aes(fill=treatment), shape=21, color="black")+theme_bw()+th+facet_wrap(~passage_hilo)+
      scale_fill_manual(values=c("grey80","firebrick4"), name="Treatment")
  }else{
    goi_low<-as.data.frame(mat_low[which(rownames(mat_low)%in%(unique(ttg$ens_gene[which(ttg$ext_gene%in%gene)]))),])
    goi_high<-as.data.frame(mat_high[which(rownames(mat_high)%in%(unique(ttg$ens_gene[which(ttg$ext_gene%in%gene)]))),])
    gene_ID<-ttg[which(ttg$ext_gene%in%gene),2:3]
    gene_ID<-gene_ID[!duplicated(gene_ID),]
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
    plt$treatment<-factor(plt$treatment, levels=c("UT","TNFa"))
    ggplot(plt, aes(passage_hilo, value,fill=treatment))+geom_boxplot(outlier.shape=NA)+
      geom_point(aes(fill=treatment), shape=21, color="black",position = position_dodge(width=0.75))+theme_bw()+th+
      facet_wrap(~label, scales="free_y", nrow=1)+ylab("Scaled reads per base")+xlab("Passage")+
      scale_fill_manual(values=c("grey80","firebrick4"), name="Treatment")
  }}


## differential between passage and DNAm diff with passage
# but actually tested in high
length(unique(sleuth_significant_high$ext_gene))
length(unique(sleuth_significant_low$ext_gene))

sleuth_significant_low_tested<-sleuth_significant_low[which(sleuth_significant_low$target_id%in%sleuth_table_high$target_id),]
sleuth_significant_high_tested<-sleuth_significant_high[which(sleuth_significant_high$target_id%in%sleuth_table_low$target_id),]

length(unique(sleuth_significant_high_tested$ext_gene))
length(unique(sleuth_significant_low_tested$ext_gene))

length(unique(sleuth_table_high$ext_gene))
length(unique(sleuth_table_low$ext_gene))

low_not_high_treatment<-sleuth_significant_low_tested[which(!(sleuth_significant_low_tested$ext_gene%in%sleuth_significant_high_tested$ext_gene)),]
length(unique(low_not_high_treatment$ext_gene))


low_not_high_treatment[grep("DNMT|TET|TLR",low_not_high_treatment$ext_gene),]

gene_exp_plot_treatment("PGLYRP4")
gene_exp_plot_treatment("ARRDC3")
gene_exp_plot_treatment("DUOXA2")

write.table(unique(low_not_high_treatment$ext_gene), file=here("data/validation","low_not_high_treatmentgenes_TNFa.txt"), quote=F, row.names = F, col.names = F)

## GO genes
gene_exp_plot_treatment("CXCL5")
gene_exp_plot_treatment("LTF")
gene_exp_plot_treatment("BCL3")
gene_exp_plot_treatment("GSN")


#' ## DNAm passage genes original 80 and validation
diff_genes_db_hypovalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypomethylation_UTUD.txt"))
diff_genes_db_hypervalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypermethylation_UTUD.txt"))
sleuth_sig_DNAm<-low_not_high_treatment[which(low_not_high_treatment$ext_gene%in%c(diff_genes_db_hypovalidation_original$V1, diff_genes_db_hypervalidation_original$V1)),]
length(unique(sleuth_sig_DNAm$ext_gene))

gene_exp_plot_treatment(c("FOXO1","ENPP6","PDK2","CNTNAP1"))
gene_exp_plot_treatment(unique(sleuth_sig_DNAm$ext_gene)[50:59])
gene_exp_plot_treatment(c("FOXO1","LDLRAD4","JAZF1","UNC5D","ENPP6","PDK2","CNTNAP1","CIITA"))

gene_exp_plot_treatment(c("FOXO1","CIITA","JAZF1","UNC5D"))
ggsave(here("figs","differenital_expression_treatment_genes_treatmentlowhigh_TNFa.pdf"),width = 10, height = 2.5)
ggsave(here("figs/jpeg","differenital_expression_treatment_genes_treatmentlowhigh_TNFa.jpeg"), width = 10, height = 2.5)


######
#'## volcano
######
source(here("scripts/validation/rnaseq/00_volcano.R"))

so_TNFa_high <- sleuth_wt(so_TNFa_high, paste0('treatmentUT'))
sleuth_table_high_wt <- sleuth_results(so_TNFa_high, 'treatmentUT', 'wt', show_all = FALSE)

sleuth_table_high_wt<-sleuth_table_high_wt[match(sleuth_table_high$target_id,sleuth_table_high_wt$target_id),]
identical(sleuth_table_high_wt$target_id, sleuth_table_high$target_id)
vol_high<-makeVolcano(sleuth_table_high$pval, -sleuth_table_high_wt$b, 2,0.015, "Difference Between\nUntreated and \nTreated", 9,10, 3000)


so_TNFa_low <- sleuth_wt(so_TNFa_low, paste0('treatmentUT'))
sleuth_table_low_wt <- sleuth_results(so_TNFa_low, 'treatmentUT', 'wt', show_all = FALSE)


sleuth_table_low_wt<-sleuth_table_low_wt[match(sleuth_table_low$target_id,sleuth_table_low_wt$target_id),]
identical(sleuth_table_low_wt$target_id, sleuth_table_low$target_id)
vol_low<-makeVolcano(sleuth_table_low$pval, -sleuth_table_low_wt$b, 2,0.015, "Difference Between\nUntreated and \nTreated", 9,10,3000)

grid.arrange(vol_low, vol_high, ncol=2)
ggsave("figs/jpeg/RNAseq_volcano_treatment_TNFa.jpeg",grid.arrange(vol_low, vol_high, ncol=2),  width = 15, height = 6)

gene_exp_plot_treatment("SOX7")


## low no high with fold change
low_fc<-sleuth_table_low_wt[which(abs(sleuth_table_low_wt$b)>2),]
sleuth_significant_low_fc<-sleuth_significant_low[which(sleuth_significant_low$target_id%in%low_fc$target_id),]
high_fc<-sleuth_table_high_wt[which(abs(sleuth_table_high_wt$b)>2),]
sleuth_significant_high_fc<-sleuth_significant_high[which(sleuth_significant_high$target_id%in%high_fc$target_id),]

low_not_high_treatment_fc<-sleuth_significant_low_fc[which(!(sleuth_significant_low_fc$ext_gene%in%sleuth_significant_high_fc$ext_gene)),]
gene_exp_plot_treatment("DUOX2")
gene_exp_plot_treatment("SLC2A6")
gene_exp_plot_treatment("IL1A")


#'### Gene Table
stats_cpgs_TNFa<-table_low_not_high(unique(sleuth_sig_DNAm$ext_gene), sleuth_sig_DNAm,sleuth_table_high)
head(stats_cpgs_TNFa)
write.csv(stats_cpgs_TNFa, file="reports_figs/validation_update/Supplement_Table7_TNFa_genes.csv", row.names = F)




#' ### heatmap

plot_transcript_heatmap(so_TNFa_high, low_not_high_treatment_fc$target_id)

plot_transcript_heatmap(so_TNFa_high,transcripts = sleuth_significant_high$target_id[1:20])




#'## R Session Info
sessionInfo()


