export RNA_DATA_DIR1=/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/run1_fastqs
export RNA_DATA_DIR2=/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/run2_fastqs

export RNA_REFS_DIR=../../data/public_rna_seq/refs

conda activate ibd_dnam_rnaseq


# ERR2270966 issues

##  QC
cd $RNA_DATA_DIR1
fastqc *.fastq.gz
multiqc --filename $RNA_DATA_DIR1/multiqc_report_raw.html $RNA_DATA_DIR1

##  QC
cd $RNA_DATA_DIR2
fastqc *.fastq.gz
multiqc --filename $RNA_DATA_DIR2/multiqc_report_raw.html $RNA_DATA_DIR2



# ##adapter trimming
# #http://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
# #mkdir trimmed

# #cutadapt adapters from truseq illumina site and -quality filter 3' ends
# ##http://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html

# awk 'NR>1{ print $3 }' $RNA_DATA_DIR/Zerbino_RNASeq_SampleInfo.txt | while read samplename 
# do
# 	cutadapt --cores=0 -q 10 -m 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $RNA_DATA_DIR/trimmed/${samplename}_adapter_trimmed_1.fastq.gz -p $RNA_DATA_DIR/trimmed/${samplename}_adapter_trimmed_2.fastq.gz $RNA_DATA_DIR/${samplename}_1.fastq.gz $RNA_DATA_DIR/${samplename}_2.fastq.gz
# done



# ##  QC
# fastqc $RNA_DATA_DIR/trimmed/*.fastq.gz
# multiqc -f --filename $RNA_DATA_DIR/trimmed/multiqc_report_trimmed.html $RNA_DATA_DIR/trimmed









