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
source(here("scripts/00_pretty_plots.R"))

num_cores=1

#set input and output dirs
datapath = "data/validation_dataset/kallisto"
resultdir = here('data/validation_dataset/kallisto/sleuth')

sampleinfo <- read.table("data/validation_dataset/sample_info_RNA_seq.txt", header=F, sep=" ", skip=1)
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
sampleinfo_UD_hilo<-sampleinfo[which(sampleinfo$differentiation=="UD" & sampleinfo$treatment=="UT"),]
table(sampleinfo_UD_hilo$passage_hilo)

#'## TI
so_UD_hilo <- sleuth_prep(sampleinfo_UD_hilo, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)


so_UD_hilo <- sleuth_fit(so_UD_hilo, ~passage_hilo, 'full')
#so_UD_hilo <- sleuth_fit(so_UD_hilo, ~passage, 'full')

so_UD_hilo <- sleuth_fit(so_UD_hilo, ~1, 'reduced')
so_UD_hilo <- sleuth_lrt(so_UD_hilo, 'reduced', 'full')
models(so_UD_hilo)
tests(so_UD_hilo)


#' ### summarize the sleuth results and view 20 most significant DE transcripts
sleuth_table <- sleuth_results(so_UD_hilo, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

sleuth_significant <- dplyr::filter(sleuth_table, pval <= 0.05)
dim(sleuth_significant)
head(sleuth_significant, 20)
head(sleuth_significant[order(sleuth_significant$pval),])

sleuth_table[grep("WNT",sleuth_table$ext_gene),]
sleuth_table[grep("LGR",sleuth_table$ext_gene),]
sleuth_table[grep("NOTCH",sleuth_table$ext_gene),]

sleuth_table[grep("SOX9",sleuth_table$ext_gene),]


plot_bootstrap(so_UD_hilo, 
               target_id = "ENSG00000159884", 
               units = "scaled_reads_per_base", 
               color_by = "passage_hilo")

plot_bootstrap(so_UD_hilo, #lgr5
               target_id = "ENSG00000139292", 
               units = "scaled_reads_per_base", 
               color_by = "passage_hilo")

plot_bootstrap(so_UD_hilo, #sox9
               target_id = "ENSG00000125398", 
               units = "scaled_reads_per_base", 
               color_by = "passage_hilo")

#'### plot an example DE gene result
plot_bootstrap(so_TI, target_id ="ENSG00000140853", units = "scaled_reads_per_base", color_by = "diagnosis")+fillscale_diagnosis
plot_bootstrap(so_TI, target_id ="ENSG00000159212", units = "scaled_reads_per_base", color_by = "diagnosis")+fillscale_diagnosis

sleuth_table_TI[which(sleuth_table_TI$target_id=="ENSG00000140853"),]
sleuth_table_TI[which(sleuth_table_TI$ext_gene=="TAP1"),]
MHCI_NLRC5 = c("TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")

sleuth_table_TI[which(sleuth_table_TI$ext_gene%in%MHCI_NLRC5),]
sleuth_significant_TI[which(sleuth_significant_TI$ext_gene%in%MHCI_NLRC5),]

mat <- sleuth:::spread_abundance_by(so_UD_hilo$obs_norm, "scaled_reads_per_base",  so_UD_hilo$sample_to_covariates$sample)


#'## gene plots
gene_exp_plot<-function(gene){
  goi<-as.data.frame(mat[which(rownames(mat)%in%(unique(ttg$ens_gene[which(ttg$ext_gene==gene)]))),])
  if(ncol(goi)==1){
    goi$sample_ID<-rownames(goi)
    colnames(goi)[1]<-"scaled_reads_per_base"
    plt<-merge(sampleinfo_UD_hilo,goi, by.x="sample",by.y="sample_ID")
    ggplot(plt, aes(passage_hilo, scaled_reads_per_base,fill=passage_hilo))+geom_boxplot()+geom_point(aes(fill=passage_hilo), shape=21, color="black")+theme_bw()+th
  }else{
    goi$gene_ID<-rownames(goi)
    goi<-melt(goi)
    colnames(goi)[3]<-"scaled_reads_per_base"
    plt<-merge(sampleinfo_UD_hilo,goi, by.x="sample",by.y="variable")
    ggplot(plt, aes(passage_hilo, scaled_reads_per_base,fill=passage_hilo))+geom_boxplot()+geom_point(aes(fill=passage_hilo), shape=21, color="black")+facet_wrap(~gene_ID)+theme_bw()+th
    }}

# gene level sumary
gene_exp_plot("LGR5")
gene_exp_plot("SOX9")


gene_exp_plot("PSMB8")
gene_exp_plot("PSMB9")
## not sig howell but in dnam
gene_exp_plot("B2M")


geneID_exp_plot<-function(gene_ID){
  goi<-as.data.frame(mat[which(rownames(mat)%in%(unique(ttg$ens_gene[which(ttg$ens_gene==gene_ID)]))),])
  goi$sample_ID<-rownames(goi)
  colnames(goi)[1]<-"scaled_reads_per_base"
  plt<-merge(sampleinfo_UD_hilo,goi, by.x="sample",by.y="sample_ID")
  label<-unique(ttg$ext_gene[which(ttg$ens_gene==gene_ID)])
  ggplot(plt, aes(passage_hilo, scaled_reads_per_base,fill=passage_hilo))+geom_boxplot()+
    geom_point(shape=21, size=2, color="black")+theme_bw()+th+ggtitle(label)}
    
grid.arrange(geneID_exp_plot("ENSG00000001626"),
  geneID_exp_plot("ENSG00000001630"))


##################
#'# Gene expression varibility in each group
##################
mat <- sleuth:::spread_abundance_by(so_UD_hilo$obs_norm, "scaled_reads_per_base",  so_UD_hilo$sample_to_covariates$sample)

sampleinfo_low_UD<-sampleinfo_UD_hilo[which(sampleinfo_UD_hilo$passage_hilo=="low"),]
sampleinfo_high_UD<-sampleinfo_UD_hilo[which(sampleinfo_UD_hilo$passage_hilo=="high"),]

mat_low_UD<-mat[,which(colnames(mat)%in%sampleinfo_low_UD$sample)]
mat_high_UD<-mat[,which(colnames(mat)%in%sampleinfo_high_UD$sample)]


#' ## Clustering on variable probes
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}

mean_count<-rowMeans(mat)

ref_range_low<-sapply(1:nrow(mat_low_UD), function(x) Variation(mat_low_UD[x,]))
ref_range_high<-sapply(1:nrow(mat_high_UD), function(x) Variation(mat_high_UD[x,]))

plt_refrange<-data.frame(gene=rownames(mat_low_UD),lowrr=ref_range_low, highrr=ref_range_high)
ggplot(plt_refrange, aes(lowrr, highrr))+geom_point()


mad_low<-sapply(1:nrow(mat_low_UD), function(x) mad(mat_low_UD[x,]))
mad_high<-sapply(1:nrow(mat_high_UD), function(x) mad(mat_high_UD[x,]))

plt_mad<-data.frame(gene=rownames(mat_low_UD),lowmad=mad_low, highmad=mad_high)
plt_mad<-plt_mad[which(mean_count!=0),]

ggplot(plt_mad, aes(lowmad, highmad))+geom_point()

plt_mad$difference<-plt_mad$lowmad-plt_mad$highmad
ggplot(plt_mad, aes(difference))+geom_histogram()+xlim(-1000,1000)

length(which(plt_mad$difference>0))
length(which(plt_mad$difference<0))


###############################################################################################################################################
#############################
#'# run sleuth on differentiation
#############################
diff_genes<-c("ASCL2","MKI67","LGR5","CA2","OLFM4","LYZ","MUC2","MUC1","CYP3A4","PLA2G2A","FABP1","KRT19","HELLS","SPINK4","FCGBP","NEAT1","TOP2A")#"TTF3",

#'## low
sampleinfo_diff_low<-sampleinfo[which(sampleinfo$passage_hilo=="low" & sampleinfo$comparison=="differentiation"),]
table(sampleinfo_diff_low$differentiation)

so_diff_low <- sleuth_prep(sampleinfo_diff_low, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)

so_diff_low <- sleuth_fit(so_diff_low, ~differentiation, 'full')
so_diff_low <- sleuth_fit(so_diff_low, ~1, 'reduced')
so_diff_low <- sleuth_lrt(so_diff_low, 'reduced', 'full')
models(so_diff_low)
tests(so_diff_low)

sleuth_table_low <- sleuth_results(so_diff_low, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_low <- dplyr::filter(sleuth_table_low, qval <= 0.05)
head(sleuth_significant_low, 20)

#sleuth_table_low[which(sleuth_table_low$ext_gene%in%diff_genes),]
sleuth_significant_low[which(sleuth_significant_low$ext_gene%in%diff_genes),]



#'## high
sampleinfo_diff_high<-sampleinfo[which(sampleinfo$passage_hilo=="high" & sampleinfo$comparison=="differentiation"),]
table(sampleinfo_diff_high$differentiation)

so_diff_high <- sleuth_prep(sampleinfo_diff_high, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)

so_diff_high <- sleuth_fit(so_diff_high, ~differentiation, 'full')
so_diff_high <- sleuth_fit(so_diff_high, ~1, 'reduced')
so_diff_high <- sleuth_lrt(so_diff_high, 'reduced', 'full')
models(so_diff_high)
tests(so_diff_high)

sleuth_table_high <- sleuth_results(so_diff_high, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_high <- dplyr::filter(sleuth_table_high, qval <= 0.05)
head(sleuth_significant_high, 20)

#sleuth_table_high[which(sleuth_table_high$ext_gene%in%diff_genes),]
sleuth_significant_high[which(sleuth_significant_high$ext_gene%in%diff_genes),]



# gene level sumary
mat_high <- sleuth:::spread_abundance_by(so_diff_high$obs_norm, "scaled_reads_per_base",  so_diff_high$sample_to_covariates$sample)
mat_low <- sleuth:::spread_abundance_by(so_diff_low$obs_norm, "scaled_reads_per_base",  so_diff_low$sample_to_covariates$sample)

gene_exp_plot_differentiation<-function(gene){
  goi_low<-as.data.frame(mat_low[which(rownames(mat_low)%in%(unique(ttg$ens_gene[which(ttg$ext_gene==gene)]))),])
  goi_high<-as.data.frame(mat_high[which(rownames(mat_high)%in%(unique(ttg$ens_gene[which(ttg$ext_gene==gene)]))),])
  
  if(ncol(goi_high)==1){
    goi_low$passage<-"low"
    goi_low$sample_ID<-rownames(goi_low)
    colnames(goi_low)[1]<-"scaled_reads_per_base"
    goi_high$passage<-"high"
    goi_high$sample_ID<-rownames(goi_high)
    colnames(goi_high)[1]<-"scaled_reads_per_base"
    goi<-rbind(goi_low, goi_high)
    plt<-merge(sampleinfo,goi, by.x="sample",by.y="sample_ID")
    ggplot(plt, aes(differentiation, scaled_reads_per_base,fill=differentiation))+geom_boxplot()+geom_point(aes(fill=differentiation), shape=21, color="black")+theme_bw()+th+facet_wrap(~passage_hilo)
    }else{
      goi_high$gene_ID<-rownames(goi_high)
      goi_low$gene_ID<-rownames(goi_low)
      goi_high<-melt(goi_high)
      goi_low<-melt(goi_low)
      goi_low$passage<-"low"
      goi_high$passage<-"high"
      
      goi<-rbind(goi_low, goi_high)
      plt<-merge(sampleinfo,goi, by.x="sample",by.y="variable")
      ggplot(plt, aes(differentiation, value,fill=differentiation))+geom_boxplot()+geom_point(aes(fill=differentiation), shape=21, color="black")+theme_bw()+th+facet_grid(gene_ID~passage_hilo)
      }}


gene_exp_plot_differentiation("LYZ")
gene_exp_plot_differentiation("LGR5")
gene_exp_plot_differentiation("ASCL2")
gene_exp_plot_differentiation("FABP1")
gene_exp_plot_differentiation("KRT19")
gene_exp_plot_differentiation("TOP2A")


low_not_high_diff<-sleuth_significant_low[which(!(sleuth_significant_low$ext_gene%in%sleuth_significant_high$ext_gene)),]

gene_exp_plot_differentiation("HELLS")
gene_exp_plot_differentiation("LYZ")


### different in low and high LYZ and HELLS


########
#'# wnt genes
########
wnt_genes<-read.table(file="/home/redgar/Documents/ibd/data/GO_0060070.txt", header=T, sep="\t")

wnt_high<-sleuth_significant_high[which(sleuth_significant_high$ext_gene%in%wnt_genes$Element),]
wnt_low<-sleuth_significant_low[which(sleuth_significant_low$ext_gene%in%wnt_genes$Element),]

wnt_low[which(!(wnt_low$ext_gene%in%wnt_high$ext_gene)),]

gene_exp_plot_differentiation("LRP6")
gene_exp_plot_differentiation("LRP5")
gene_exp_plot_differentiation("LGR4")
gene_exp_plot_differentiation("PTEN")
gene_exp_plot_differentiation("WNT8B")




#############
#'## volcano
#############
source(here("scripts/validation/rnaseq/00_volcano.R"))

so_diff_high <- sleuth_wt(so_diff_high, paste0('differentiationUD'))
sleuth_table_high_wt <- sleuth_results(so_diff_high, 'differentiationUD', 'wt', show_all = FALSE)

vol_high<-makeVolcano(sleuth_table_high$pval, sleuth_table_high_wt$b, 2,0.015, "Difference Between\nDifferentiated and \nUndifferentiated", 5,11)


so_diff_low <- sleuth_wt(so_diff_low, paste0('differentiationUD'))
sleuth_table_low_wt <- sleuth_results(so_diff_low, 'differentiationUD', 'wt', show_all = FALSE)

vol_low<-makeVolcano(sleuth_table_low$pval, sleuth_table_low_wt$b, 2,0.015, "Difference Between\nDifferentiated and \nUndifferentiated", 5,11)

grid.arrange(vol_low, vol_high, ncol=2)
ggsave("figs/jpeg/RNAseq_volcano_differentiation.jpeg",grid.arrange(vol_low, vol_high, ncol=2),  width = 15, height = 6)








#############################
#'# run sleuth on differentiation interaction terms
#############################
diff_genes<-c("ASCL2","MKI67","LGR5","CA2","OLFM4","LYZ","MUC2","MUC1","CYP3A4","PLA2G2A","FABP1","KRT19","HELLS","SPINK4","FCGBP","NEAT1","TOP2A")#"TTF3",

#'## low
sampleinfo_diff<-sampleinfo[which(sampleinfo$comparison=="differentiation"),]
table(sampleinfo_diff$differentiation)

so_diff <- sleuth_prep(sampleinfo_diff, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)

so_diff <- sleuth_fit(so_diff, ~differentiation + passage_hilo + differentiation:passage_hilo, 'full')
so_diff <- sleuth_fit(so_diff, ~differentiation + passage_hilo, 'reduced')
so_diff <- sleuth_lrt(so_diff, 'reduced', 'full')
models(so_diff)
tests(so_diff)

sleuth_table <- sleuth_results(so_diff, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

#sleuth_table[which(sleuth_table$ext_gene%in%diff_genes),]
sleuth_significant[which(sleuth_significant$ext_gene%in%diff_genes),]

sleuth_significant <- dplyr::filter(sleuth_table, pval <= 0.01)
head(sleuth_significant, 20)

gene_exp_plot_differentiation("CIDEC")










#####

#'  #'## gene plots
#' gene_exp_plot_differentiation<-function(gene, mat, sampleinf){
#'   goi<-as.data.frame(mat[which(rownames(mat)%in%(unique(ttg$ens_gene[which(ttg$ext_gene==gene)]))),])
#'   if(ncol(goi)==1){
#'     goi$sample_ID<-rownames(goi)
#'     colnames(goi)[1]<-"scaled_reads_per_base"
#'     plt<-merge(sampleinf,goi, by.x="sample",by.y="sample_ID")
#'     ggplot(plt, aes(differentiation, scaled_reads_per_base,fill=differentiation))+geom_boxplot()+geom_point(aes(fill=differentiation), shape=21, color="black")+theme_bw()+th
#'   }else{
#'     goi$gene_ID<-rownames(goi)
#'     goi<-melt(goi)
#'     colnames(goi)[3]<-"scaled_reads_per_base"
#'     plt<-merge(sampleinf,goi, by.x="sample",by.y="variable")
#'     ggplot(plt, aes(differentiation, scaled_reads_per_base,fill=differentiation))+geom_boxplot()+geom_point(aes(fill=differentiation), shape=21, color="black")+facet_wrap(~gene_ID)+theme_bw()+th
#'   }}
#' 
#' # gene level sumary
#' mat_high <- sleuth:::spread_abundance_by(so_diff_high$obs_norm, "scaled_reads_per_base",  so_diff_high$sample_to_covariates$sample)
#' mat_low <- sleuth:::spread_abundance_by(so_diff_low$obs_norm, "scaled_reads_per_base",  so_diff_low$sample_to_covariates$sample)

sleuth_significant <- dplyr::filter(sleuth_table, pval <= 0.05)
dim(sleuth_significant)
head(sleuth_significant, 20)
head(sleuth_significant[order(sleuth_significant$pval),])

sleuth_table[grep("WNT",sleuth_table$ext_gene),]
sleuth_table[grep("LGR",sleuth_table$ext_gene),]
sleuth_table[grep("NOTCH",sleuth_table$ext_gene),]

sleuth_table[grep("SOX9",sleuth_table$ext_gene),]




#'## R Session Info
sessionInfo()


