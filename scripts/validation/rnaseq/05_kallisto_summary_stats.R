library(ggplot2)
library(gridExtra)
library(here)

options(stringsAsFactors = FALSE)
source("scripts/00_pretty_plots.R")

datapath = "data/validation/kallisto"
sampleinfo <- read.table("data/validation/sample_info_RNA_seq.txt", header=F, sep=" ", skip=1)
sampleinfo<-sampleinfo[,c(2,7,8:10,12)]
colnames(sampleinfo)<-c("sample","concentration","volume","quantity","ratio","well")
sample_id = sampleinfo$sample
kal_dirs <- file.path(datapath, sample_id)

summary_stats<-do.call(rbind,lapply(1:length(sample_id), function(x){
  json<-read.table(paste(kal_dirs[x],"/run_info.json", sep=""), skip=1, nrow = 10, sep="\t")
  fragments<-strsplit(json$V2[grep("n_processed", json$V2)]," |,")[[1]][2]
  p_pseudoaligned<-strsplit(json$V2[grep("p_pseudoaligned", json$V2)]," |,")[[1]][2]
  n_pseudoaligned<-strsplit(json$V2[grep("n_pseudoaligned", json$V2)]," |,")[[1]][2]
  data.frame(sample_id=sample_id[x],fragments=as.numeric(fragments),p_pseudoaligned=as.numeric(p_pseudoaligned), n_pseudoaligned=as.numeric(n_pseudoaligned))
  }))


processed<-ggplot(summary_stats, aes(fragments/1000000))+geom_histogram(fill="lightgrey", color="black")+th+theme_bw()+xlab("Reads Processed (Millions)")
mapped<-ggplot(summary_stats, aes(p_pseudoaligned))+geom_histogram(fill="lightgrey", color="black")+th+theme_bw()+xlab("Percent Reads Pseudoaligned")

panel<-grid.arrange(processed, mapped, nrow=1)

ggsave("figs/pseudoalignment_stats.pdf", panel,w=10, h=5)
ggsave("figs/jpeg/pseudoalignment_stats.jpeg",panel, w=10, h=5)




counts<-read.table(here("data/validation/merged/fastqc_sequence_counts_plot.tsv"), sep="\t", header=T)
counts$total<-rowSums(counts[,2:3])
counts$perdup<-counts$Duplicate.Reads/counts$total

counts$sample_id<-gsub("_merged","",counts$Category)

counts_align<-merge(summary_stats,counts, by="sample_id")

ggplot(counts_align, aes(total, p_pseudoaligned))+geom_point()+theme_bw()
ggplot(counts_align, aes(perdup, p_pseudoaligned))+geom_point()+theme_bw()


