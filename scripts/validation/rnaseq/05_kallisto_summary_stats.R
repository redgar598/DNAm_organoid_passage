library(ggplot2)
library(gridExtra)
library(here)

options(stringsAsFactors = FALSE)
source("scripts/00_pretty_plots.R")

datapath = "data/validation_dataset/kallisto"
sampleinfo <- read.table("data/validation_dataset/sample_info_RNA_seq.txt", header=F, sep=" ", skip=1)
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




### correlate to multipqc statistics
overrepresented<-read.table("../../data/public_rna_seq/multiqc_report_raw_data/mqc_fastqc_overrepresented_sequencesi_plot_1.txt", sep="\t", header=T)
overrepresented$all_over<-rowSums(overrepresented[,2:3])
overrepresented$Sample_run<-gsub(".*_","",overrepresented$Sample)
overrepresented$Sample<-gsub("_1|_2","",overrepresented$Sample)
over_align<-merge(summary_stats,overrepresented, by.x="sample_id", by.y="Sample")
ggplot(over_align, aes(all_over, p_pseudoaligned))+geom_point()+theme_bw()+xlim(0,50)+ylim(25,75)+xlab("Percent Overrepresented Sequences")

ggsave("../../figs/pseudoalignment_overrepresented_stats.pdf",w=5.2, h=5)
ggsave("../../figs/jpeg/pseudoalignment_overrepresented_stats.jpeg", w=5.2, h=5)


##just the overrepresented bar plot
ggplot(overrepresented, aes(reorder(Sample, all_over),all_over, fill=Sample_run))+geom_bar(position="dodge",stat="identity")+
  theme_bw()+theme(axis.text.x=element_text(angle=-90))+
  scale_fill_manual(values=c("#9ecae1","#6baed6"))+xlab("Sample")+ylim(0,50)
  
ggsave("../../figs/overrepresented_stats.pdf",w=10, h=5)
ggsave("../../figs/jpeg/overrepresented_stats.jpeg", w=10, h=5)



## number aligned also associated to over represented
ggplot(over_align, aes(all_over, n_pseudoaligned))+geom_point()+theme_bw()




counts<-read.table("../../data/public_rna_seq/multiqc_report_raw_data/mqc_fastqc_sequence_counts_plot_1.txt", sep="\t", header=T)
counts$total<-rowSums(counts[,2:3])
counts$perdup<-counts$Duplicate.Reads/counts$total

counts$Sample<-gsub("_1|_2","",overrepresented$Sample)

counts_align<-merge(summary_stats,counts, by.x="sample_id", by.y="Sample")

ggplot(counts_align, aes(total, p_pseudoaligned))+geom_point()+theme_bw()
ggplot(counts_align, aes(perdup, p_pseudoaligned))+geom_point()+theme_bw()


