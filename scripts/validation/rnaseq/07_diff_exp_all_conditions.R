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

        #' #############
        #' #'# run sleuth on passage high low all undifferentiated
        #' #############
        #' sampleinfo_UD_hilo<-sampleinfo[which(sampleinfo$differentiation=="UD" & sampleinfo$treatment=="UT"),]
        #' table(sampleinfo_UD_hilo$passage_hilo)
        #' 
        #' #'## TI
        #' so_UD_hilo <- sleuth_prep(sampleinfo_UD_hilo, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)
        #' 
        #' 
        #' so_UD_hilo <- sleuth_fit(so_UD_hilo, ~passage_hilo, 'full')
        #' #so_UD_hilo <- sleuth_fit(so_UD_hilo, ~passage, 'full')
        #' 
        #' so_UD_hilo <- sleuth_fit(so_UD_hilo, ~1, 'reduced')
        #' so_UD_hilo <- sleuth_lrt(so_UD_hilo, 'reduced', 'full')
        #' models(so_UD_hilo)
        #' tests(so_UD_hilo)
        #' 
        #' 
        #' #' ### summarize the sleuth results and view 20 most significant DE transcripts
        #' sleuth_table <- sleuth_results(so_UD_hilo, 'reduced:full', 'lrt', show_all = FALSE)
        #' sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
        #' head(sleuth_significant, 20)
        #' 
        #' sleuth_significant <- dplyr::filter(sleuth_table, pval <= 0.05)
        #' dim(sleuth_significant)
        #' head(sleuth_significant, 20)
        #' head(sleuth_significant[order(sleuth_significant$pval),])
        #' 
        #' sleuth_table[grep("WNT",sleuth_table$ext_gene),]
        #' sleuth_table[grep("LGR",sleuth_table$ext_gene),]
        #' sleuth_table[grep("NOTCH",sleuth_table$ext_gene),]
        #' 
        #' sleuth_table[grep("SOX9",sleuth_table$ext_gene),]
        #' 
        #' 
        #' plot_bootstrap(so_UD_hilo, 
        #'                target_id = "ENSG00000159884", 
        #'                units = "scaled_reads_per_base", 
        #'                color_by = "passage_hilo")
        #' 
        #' plot_bootstrap(so_UD_hilo, #lgr5
        #'                target_id = "ENSG00000139292", 
        #'                units = "scaled_reads_per_base", 
        #'                color_by = "passage_hilo")
        #' 
        #' plot_bootstrap(so_UD_hilo, #sox9
        #'                target_id = "ENSG00000125398", 
        #'                units = "scaled_reads_per_base", 
        #'                color_by = "passage_hilo")
        #' 
        #' #'### plot an example DE gene result
        #' plot_bootstrap(so_TI, target_id ="ENSG00000140853", units = "scaled_reads_per_base", color_by = "diagnosis")+fillscale_diagnosis
        #' plot_bootstrap(so_TI, target_id ="ENSG00000159212", units = "scaled_reads_per_base", color_by = "diagnosis")+fillscale_diagnosis
        #' 
        #' sleuth_table_TI[which(sleuth_table_TI$target_id=="ENSG00000140853"),]
        #' sleuth_table_TI[which(sleuth_table_TI$ext_gene=="TAP1"),]
        #' MHCI_NLRC5 = c("TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")
        #' 
        #' sleuth_table_TI[which(sleuth_table_TI$ext_gene%in%MHCI_NLRC5),]
        #' sleuth_significant_TI[which(sleuth_significant_TI$ext_gene%in%MHCI_NLRC5),]
        #' 
        #' mat <- sleuth:::spread_abundance_by(so_UD_hilo$obs_norm, "scaled_reads_per_base",  so_UD_hilo$sample_to_covariates$sample)
        #' 
        #' 
        #' #'## gene plots
        #' gene_exp_plot<-function(gene){
        #'   goi<-as.data.frame(mat[which(rownames(mat)%in%(unique(ttg$ens_gene[which(ttg$ext_gene==gene)]))),])
        #'   if(ncol(goi)==1){
        #'     goi$sample_ID<-rownames(goi)
        #'     colnames(goi)[1]<-"scaled_reads_per_base"
        #'     plt<-merge(sampleinfo_UD_hilo,goi, by.x="sample",by.y="sample_ID")
        #'     ggplot(plt, aes(passage_hilo, scaled_reads_per_base,fill=passage_hilo))+geom_boxplot()+geom_point(aes(fill=passage_hilo), shape=21, color="black")+theme_bw()+th
        #'   }else{
        #'     goi$gene_ID<-rownames(goi)
        #'     goi<-melt(goi)
        #'     colnames(goi)[3]<-"scaled_reads_per_base"
        #'     plt<-merge(sampleinfo_UD_hilo,goi, by.x="sample",by.y="variable")
        #'     ggplot(plt, aes(passage_hilo, scaled_reads_per_base,fill=passage_hilo))+geom_boxplot()+geom_point(aes(fill=passage_hilo), shape=21, color="black")+facet_wrap(~gene_ID)+theme_bw()+th
        #'     }}
        #' 
        #' # gene level sumary
        #' gene_exp_plot("LGR5")
        #' gene_exp_plot("SOX9")
        #' 
        #' 
        #' 
        #'           # geneID_exp_plot<-function(gene_ID){
        #'           #   goi<-as.data.frame(mat[which(rownames(mat)%in%(unique(ttg$ens_gene[which(ttg$ens_gene==gene_ID)]))),])
        #'           #   goi$sample_ID<-rownames(goi)
        #'           #   colnames(goi)[1]<-"scaled_reads_per_base"
        #'           #   plt<-merge(sampleinfo_UD_hilo,goi, by.x="sample",by.y="sample_ID")
        #'           #   label<-unique(ttg$ext_gene[which(ttg$ens_gene==gene_ID)])
        #'           #   ggplot(plt, aes(passage_hilo, scaled_reads_per_base,fill=passage_hilo))+geom_boxplot()+
        #'           #     geom_point(shape=21, size=2, color="black")+theme_bw()+th+ggtitle(label)}
        #'           #     
        #'           # grid.arrange(geneID_exp_plot("ENSG00000001626"),
        #'           #   geneID_exp_plot("ENSG00000001630"))
        #' 
        #' 
        #' #' ## DNAm passage genes
        #' diff_genes_db_hypovalidation<-read.table(file=here("data/validation/DNAm/","validation_genes_hypomethylation_hilo.txt"))
        #' diff_genes_db_hypervalidation<-read.table(file=here("data/validation/DNAm/","validation_genes_hypermethylation_hilo.txt"))
        #' 
        #' hypo_diffexp<-sleuth_significant[which(sleuth_significant$ext_gene%in%diff_genes_db_hypovalidation$V1),]
        #' hyper_diffexp<-sleuth_significant[which(sleuth_significant$ext_gene%in%diff_genes_db_hypervalidation$V1),]
        #' 
        #' head(hypo_diffexp[order(hypo_diffexp$pval),])
        #' head(hyper_diffexp[order(hyper_diffexp$pval),])
        #' 
        #' gene_exp_plot("KANK1")
        #' gene_exp_plot("AQP1")
        #' gene_exp_plot("NOTCH1")

        #' ##################
        #' #'# Gene expression varibility in each group
        #' ##################
        #' mat <- sleuth:::spread_abundance_by(so_UD_hilo$obs_norm, "scaled_reads_per_base",  so_UD_hilo$sample_to_covariates$sample)
        #' 
        #' sampleinfo_low_UD<-sampleinfo_UD_hilo[which(sampleinfo_UD_hilo$passage_hilo=="low"),]
        #' sampleinfo_high_UD<-sampleinfo_UD_hilo[which(sampleinfo_UD_hilo$passage_hilo=="high"),]
        #' 
        #' mat_low_UD<-mat[,which(colnames(mat)%in%sampleinfo_low_UD$sample)]
        #' mat_high_UD<-mat[,which(colnames(mat)%in%sampleinfo_high_UD$sample)]
        #' 
        #' 
        #' #' ## Clustering on variable probes
        #' Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
        #' 
        #' mean_count<-rowMeans(mat)
        #' 
        #' ref_range_low<-sapply(1:nrow(mat_low_UD), function(x) Variation(mat_low_UD[x,]))
        #' ref_range_high<-sapply(1:nrow(mat_high_UD), function(x) Variation(mat_high_UD[x,]))
        #' 
        #' plt_refrange<-data.frame(gene=rownames(mat_low_UD),lowrr=ref_range_low, highrr=ref_range_high)
        #' ggplot(plt_refrange, aes(lowrr, highrr))+geom_point()
        #' 
        #' 
        #' mad_low<-sapply(1:nrow(mat_low_UD), function(x) mad(mat_low_UD[x,]))
        #' mad_high<-sapply(1:nrow(mat_high_UD), function(x) mad(mat_high_UD[x,]))
        #' 
        #' plt_mad<-data.frame(gene=rownames(mat_low_UD),lowmad=mad_low, highmad=mad_high)
        #' plt_mad<-plt_mad[which(mean_count!=0),]
        #' 
        #' ggplot(plt_mad, aes(lowmad, highmad))+geom_point()
        #' 
        #' plt_mad$difference<-plt_mad$lowmad-plt_mad$highmad
        #' ggplot(plt_mad, aes(difference))+geom_histogram()+xlim(-1000,1000)
        #' 
        #' length(which(plt_mad$difference>0))
        #' length(which(plt_mad$difference<0))
        #' 

#############
#'# run sleuth on passage all organoids (and continious so most comparable to DNAm)
#############
table(sampleinfo$passage_hilo, sampleinfo$individual)

so_all_hilo <- sleuth_prep(sampleinfo, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)

so_all_hilo <- sleuth_fit(so_all_hilo, ~passage + individual , 'full')
so_all_hilo <- sleuth_fit(so_all_hilo, ~1 + individual, 'reduced')
so_all_hilo <- sleuth_lrt(so_all_hilo, 'reduced', 'full')
models(so_all_hilo)
tests(so_all_hilo)


#' ### summarize the sleuth results and view 20 most significant DE transcripts
sleuth_table <- sleuth_results(so_all_hilo, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

# sleuth_significant <- dplyr::filter(sleuth_table, pval <= 0.05)
# dim(sleuth_significant)
# head(sleuth_significant, 20)
# head(sleuth_significant[order(sleuth_significant$pval),])

sleuth_table[grep("WNT",sleuth_table$ext_gene),]
sleuth_table[grep("LGR",sleuth_table$ext_gene),]
sleuth_table[grep("NOTCH",sleuth_table$ext_gene),]



#'## gene plots
mat <- sleuth:::spread_abundance_by(so_all_hilo$obs_norm, "scaled_reads_per_base",  so_all_hilo$sample_to_covariates$sample)

gene_exp_plot<-function(gene){
  goi<-as.data.frame(mat[which(rownames(mat)%in%(unique(ttg$ens_gene[which(ttg$ext_gene==gene)]))),])
  if(ncol(goi)==1){
    goi$sample_ID<-rownames(goi)
    colnames(goi)[1]<-"scaled_reads_per_base"
    plt<-merge(sampleinfo,goi, by.x="sample",by.y="sample_ID")
    ggplot(plt, aes(passage, scaled_reads_per_base,fill=passage))+geom_point(aes(fill=passage_hilo), shape=21, color="black")+theme_bw()+th+
      stat_smooth(method="lm", se=F, color="grey70")
  }else{
    goi$gene_ID<-rownames(goi)
    goi<-melt(goi)
    colnames(goi)[3]<-"scaled_reads_per_base"
    plt<-merge(sampleinfo,goi, by.x="sample",by.y="variable")
    ggplot(plt, aes(passage, scaled_reads_per_base,fill=passage))+geom_point(aes(fill=passage_hilo), shape=21, color="black")+facet_wrap(~gene_ID)+theme_bw()+th+
      stat_smooth(method="lm", se=F, color="grey70")
  }}

# gene level sumary
gene_exp_plot("SPART")
gene_exp_plot("FASTKD5")




#' ## DNAm passage genes original 80 and validation
diff_genes_db_hypovalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypomethylation.txt"))
diff_genes_db_hypervalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypermethylation.txt"))

hypo_diffexp_original<-sleuth_significant[which(sleuth_significant$ext_gene%in%diff_genes_db_hypovalidation_original$V1),]
hyper_diffexp_original<-sleuth_significant[which(sleuth_significant$ext_gene%in%diff_genes_db_hypervalidation_original$V1),]

head(hypo_diffexp_original[order(hypo_diffexp_original$pval),])
head(hyper_diffexp_original[order(hyper_diffexp_original$pval),])

gene_exp_plot("PIK3R3")



## fdr on just diff DNAM
sleuth_sig_DNAm<-sleuth_table[which(sleuth_table$ext_gene%in%c(diff_genes_db_hypovalidation_original$V1, diff_genes_db_hypervalidation_original$V1)),]
pdjust_dnam<-p.adjust(sleuth_sig_DNAm$pval, method="fdr", n=nrow(sleuth_sig_DNAm))


## wnt genes
wnt_genes<-read.table(file="/home/redgar/Documents/ibd/data/GO_0060070.txt", header=T, sep="\t")
hypo_diffexp_original[which(hypo_diffexp_original$ext_gene%in%wnt_genes$Symbol),]
hyper_diffexp_original[which(hyper_diffexp_original$ext_gene%in%wnt_genes$Symbol),]

sleuth_significant[which(sleuth_significant$ext_gene%in%wnt_genes$Symbol),]
gene_exp_plot("PTEN")
gene_exp_plot("RECK")


sleuth_table[which(sleuth_table$ext_gene%in%wnt_genes$Symbol),]
# RECK is differentially expressed and hypermethylated, PTEN is just differenitall expressed



## fold change
so_all_hilo <- sleuth_wt(so_all_hilo, paste0('passage'))
sleuth_table_wt <- sleuth_results(so_all_hilo, 'passage', 'wt', show_all = FALSE)

              # hist(sleuth_table_wt[which(sleuth_table_wt$target_id%in%hypo_diffexp_original$target_id),]$b)
              # head(sleuth_table_wt[which(sleuth_table_wt$target_id%in%hypo_diffexp_original$target_id),])
              # 
              # hist(sleuth_table_wt[which(sleuth_table_wt$target_id%in%hyper_diffexp_original$target_id),]$b)
              # head(sleuth_table_wt[which(sleuth_table_wt$target_id%in%hyper_diffexp_original$target_id),])
              
              # gene_exp_plot("CLK1")
              # gene_exp_plot("TESK2")
              # gene_exp_plot("TESK2")


#'### more overlapping than chance?
EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194

hypo_rnd<-nrow(diff_genes_db_hypovalidation_original)
EPIC_genes_unique<-unique(EPIC_genes$Gene.name)

rnd_overlap_hypo<-sapply(1:100, function(x){
  set.seed(x)
  rnd_genes<-EPIC_genes_unique[sample(1:length(EPIC_genes_unique), hypo_rnd)]
  length(intersect(rnd_genes, sleuth_significant$ext_gene))})

nrow(hypo_diffexp_original)

(length(which(rnd_overlap_hypo>nrow(hypo_diffexp_original)))+1)/(100+1)
(length(which(rnd_overlap_hypo<nrow(hypo_diffexp_original)))+1)/(100+1)


hyper_rnd<-nrow(diff_genes_db_hypervalidation_original)

rnd_overlap_hyper<-sapply(1:100, function(x){
  set.seed(x)
  rnd_genes<-EPIC_genes_unique[sample(1:length(EPIC_genes_unique), hyper_rnd)]
  length(intersect(rnd_genes, sleuth_significant$ext_gene))})

nrow(hyper_diffexp_original)

(length(which(rnd_overlap_hyper>nrow(hyper_diffexp_original)))+1)/(100+1)
(length(which(rnd_overlap_hyper<nrow(hyper_diffexp_original)))+1)/(100+1)



###############################################################################################################################################
#############################
#'# run sleuth on differentiation
#############################
diff_genes<-c("ASCL2","MKI67","LGR5","CA2","OLFM4","LYZ","MUC2","MUC1","CYP3A4","PLA2G2A","FABP1","KRT19","HELLS","SPINK4","FCGBP","NEAT1","TOP2A")#"TTF3",

#'## low
sampleinfo_diff_low<-sampleinfo[which(sampleinfo$passage_hilo=="low" & sampleinfo$comparison=="differentiation"),]
table(sampleinfo_diff_low$differentiation, sampleinfo_diff_low$individual)

so_diff_low <- sleuth_prep(sampleinfo_diff_low, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)

so_diff_low <- sleuth_fit(so_diff_low, ~differentiation +individual, 'full')
so_diff_low <- sleuth_fit(so_diff_low, ~1 + individual, 'reduced')
so_diff_low <- sleuth_lrt(so_diff_low, 'reduced', 'full')
models(so_diff_low)
tests(so_diff_low)

sleuth_table_low <- sleuth_results(so_diff_low, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_low <- dplyr::filter(sleuth_table_low, qval <= 0.05)
head(sleuth_significant_low, 20)

#sleuth_table_low[which(sleuth_table_low$ext_gene%in%diff_genes),]
diff_markers_low<-sleuth_significant_low[which(sleuth_significant_low$ext_gene%in%diff_genes),]



#'## high
sampleinfo_diff_high<-sampleinfo[which(sampleinfo$passage_hilo=="high" & sampleinfo$comparison=="differentiation"),]
table(sampleinfo_diff_high$differentiation)

so_diff_high <- sleuth_prep(sampleinfo_diff_high, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, gene_mode = TRUE)

so_diff_high <- sleuth_fit(so_diff_high, ~differentiation +individual, 'full')
so_diff_high <- sleuth_fit(so_diff_high, ~1 +individual, 'reduced')
so_diff_high <- sleuth_lrt(so_diff_high, 'reduced', 'full')
models(so_diff_high)
tests(so_diff_high)

sleuth_table_high <- sleuth_results(so_diff_high, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_high <- dplyr::filter(sleuth_table_high, qval <= 0.05)
head(sleuth_significant_high, 20)

#sleuth_table_high[which(sleuth_table_high$ext_gene%in%diff_genes),]
diff_markers_high<-sleuth_significant_high[which(sleuth_significant_high$ext_gene%in%diff_genes),]


## low not high
diff_markers_low[which(!(diff_markers_low$ext_gene%in%diff_markers_high$ext_gene)),]
diff_markers_low[which(!(diff_markers_low$target_id%in%diff_markers_high$target_id)),]


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
gene_exp_plot_differentiation("GPR135")
gene_exp_plot_differentiation("KATNAL2")

write.table(unique(low_not_high_diff$ext_gene), file=here("data/validation","low_not_high_differentiationgenes.txt"), quote=F, row.names = F, col.names = F)

low_not_high_diff[which(low_not_high_diff$ext_gene%in%diff_genes),]

gene_exp_plot_differentiation("HSPA1A")
gene_exp_plot_differentiation("LYZ")
### different in low and high LYZ 

## interesting for GO enrichment
gene_exp_plot_differentiation("FUS")



########
#'# wnt genes
########
wnt_genes<-read.table(file="/home/redgar/Documents/ibd/data/GO_0060070.txt", header=T, sep="\t")

wnt_high<-sleuth_significant_high[which(sleuth_significant_high$ext_gene%in%wnt_genes$Element),]
wnt_low<-sleuth_significant_low[which(sleuth_significant_low$ext_gene%in%wnt_genes$Element),]

wnt_low[which(!(wnt_low$ext_gene%in%wnt_high$ext_gene)),]

gene_exp_plot_differentiation("WNT8B")
gene_exp_plot_differentiation("RECK")
gene_exp_plot_differentiation("AMER1")
gene_exp_plot_differentiation("FZD8")
gene_exp_plot_differentiation("DISC1")




#############
#'## volcano
#############
source(here("scripts/validation/rnaseq/00_volcano.R"))

so_diff_high <- sleuth_wt(so_diff_high, paste0('differentiationUD'))
sleuth_table_high_wt <- sleuth_results(so_diff_high, 'differentiationUD', 'wt', show_all = FALSE)

vol_high<-makeVolcano(sleuth_table_high$pval, sleuth_table_high_wt$b, 2,0.015, "Difference Between\nDifferentiated and \nUndifferentiated", 5,11,9000)


so_diff_low <- sleuth_wt(so_diff_low, paste0('differentiationUD'))
sleuth_table_low_wt <- sleuth_results(so_diff_low, 'differentiationUD', 'wt', show_all = FALSE)

vol_low<-makeVolcano(sleuth_table_low$pval, sleuth_table_low_wt$b, 2,0.015, "Difference Between\nDifferentiated and \nUndifferentiated", 5,11,9000)

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

so_diff <- sleuth_fit(so_diff, ~differentiation + passage_hilo + differentiation:passage_hilo +individual, 'full')
so_diff <- sleuth_fit(so_diff, ~differentiation + passage_hilo +individual, 'reduced')
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


gene_exp_plot_differentiation("SEMA3B")






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

gene_exp_plot_treatment<-function(gene){
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
    ggplot(plt, aes(treatment, scaled_reads_per_base,fill=treatment))+geom_boxplot()+geom_point(aes(fill=treatment), shape=21, color="black")+theme_bw()+th+facet_wrap(~passage_hilo)
  }else{
    goi_high$gene_ID<-rownames(goi_high)
    goi_low$gene_ID<-rownames(goi_low)
    goi_high<-melt(goi_high)
    goi_low<-melt(goi_low)
    goi_low$passage<-"low"
    goi_high$passage<-"high"
    
    goi<-rbind(goi_low, goi_high)
    plt<-merge(sampleinfo,goi, by.x="sample",by.y="variable")
    ggplot(plt, aes(treatment, value,fill=treatment))+geom_boxplot()+geom_point(aes(fill=treatment), shape=21, color="black")+theme_bw()+th+facet_grid(gene_ID~passage_hilo)
  }}


gene_exp_plot_treatment("HLA-B")
gene_exp_plot_treatment("NLRC5")
gene_exp_plot_treatment("TAP1")
gene_exp_plot_treatment("MR1")
gene_exp_plot_treatment("PSMB8")
gene_exp_plot_treatment("CD1D")




low_not_high_treatment<-sleuth_significant_low[which(!(sleuth_significant_low$ext_gene%in%sleuth_significant_high$ext_gene)),]
gene_exp_plot_treatment("NDUFA8")
gene_exp_plot_treatment("NDUFA7")
gene_exp_plot_treatment("SEC24B")

write.table(unique(low_not_high_treatment$ext_gene), file=here("data/validation","low_not_high_treatmentgenes.txt"), quote=F, row.names = F, col.names = F)


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
    ggplot(plt, aes(treatment, scaled_reads_per_base,fill=treatment))+geom_boxplot()+geom_point(aes(fill=treatment), shape=21, color="black")+theme_bw()+th+facet_wrap(~passage_hilo)
  }else{
    goi_high$gene_ID<-rownames(goi_high)
    goi_low$gene_ID<-rownames(goi_low)
    goi_high<-melt(goi_high)
    goi_low<-melt(goi_low)
    goi_low$passage<-"low"
    goi_high$passage<-"high"
    
    goi<-rbind(goi_low, goi_high)
    plt<-merge(sampleinfo,goi, by.x="sample",by.y="variable")
    ggplot(plt, aes(treatment, value,fill=treatment))+geom_boxplot()+geom_point(aes(fill=treatment), shape=21, color="black")+theme_bw()+th+facet_grid(gene_ID~passage_hilo)
  }}



low_not_high_treatment<-sleuth_significant_low[which(!(sleuth_significant_low$ext_gene%in%sleuth_significant_high$ext_gene)),]
gene_exp_plot_treatment("PGLYRP4")
gene_exp_plot_treatment("ARRDC3")
gene_exp_plot_treatment("DUOXA2")

write.table(unique(low_not_high_treatment$ext_gene), file=here("data/validation","low_not_high_treatmentgenes_TNFa.txt"), quote=F, row.names = F, col.names = F)

## GO genes
gene_exp_plot_treatment("CXCL5")
gene_exp_plot_treatment("LTF")
gene_exp_plot_treatment("BCL3")
gene_exp_plot_treatment("GSN")


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




#' ### heatmap

plot_transcript_heatmap(so_TNFa_high, low_not_high_treatment_fc$target_id)

plot_transcript_heatmap(so_TNFa_high,transcripts = sleuth_significant_high$target_id[1:20])




#'## R Session Info
sessionInfo()


