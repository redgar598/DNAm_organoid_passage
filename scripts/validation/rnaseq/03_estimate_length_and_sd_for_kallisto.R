length_reads<-read.table("data/validation_dataset/multiQC/run1/fastqc_sequence_length_distribution_plot.tsv", sep="\t", header=T)

# from https://math.stackexchange.com/questions/857566/how-to-get-the-standard-deviation-of-a-given-histogram-image
run1<-do.call(rbind, lapply(2:ncol(length_reads), function(samp) {
  mean_read_length<-(1/sum(length_reads[,samp]))*sum(sapply(1:nrow(length_reads), function(x) length_reads[x,1]*length_reads[x,samp]))
  sd_squrd<-(1/sum(length_reads$T317_SC_p12_D_S15_R1_001))*sum(sapply(1:nrow(length_reads), function(x) length_reads[x,samp]*((length_reads[x,1]-mean_read_length)^2)))
  data.frame(sample=colnames(length_reads)[samp], mean_read_length=mean_read_length, sd=sqrt(sd_squrd))}))


length_reads<-read.table("data/validation_dataset/multiQC/run2/fastqc_sequence_length_distribution_plot.tsv", sep="\t", header=T)

# from https://math.stackexchange.com/questions/857566/how-to-get-the-standard-deviation-of-a-given-histogram-image
run2<-do.call(rbind, lapply(2:ncol(length_reads), function(samp) {
  mean_read_length<-(1/sum(length_reads[,samp]))*sum(sapply(1:nrow(length_reads), function(x) length_reads[x,1]*length_reads[x,samp]))
  sd_squrd<-(1/sum(length_reads$T317_SC_p12_D_S15_R1_001))*sum(sapply(1:nrow(length_reads), function(x) length_reads[x,samp]*((length_reads[x,1]-mean_read_length)^2)))
  data.frame(sample=colnames(length_reads)[samp], mean_read_length=mean_read_length, sd=sqrt(sd_squrd))}))

mean(c(run1$mean_read_length, run2$mean_read_length))
mean(c(run1$sd, run2$sd))
