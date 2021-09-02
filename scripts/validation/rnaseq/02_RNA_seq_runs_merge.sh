export RNA_DATA_DIR1=/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/run1_fastqs
export RNA_DATA_DIR2=/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/run2_fastqs

export RNA_REFS_DIR=/nfs/research1/zerbino/redgar/ibd/data/public_rna_seq/refs
export RNA_DATA_MERGED=/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/merged


awk 'NR>1{ print $2 }' ../sample_info_RNA_seq.txt | while read samplename 
do
	cat $RNA_DATA_DIR1/${samplename}_*_R1_001.fastq.gz $RNA_DATA_DIR2/${samplename}_*_R1_001.fastq.gz > $RNA_DATA_MERGED/${samplename}_merged.fastq.gz
done


##  QC
cd $RNA_DATA_MERGED
fastqc *.fastq.gz
multiqc --filename $RNA_DATA_MERGED/multiqc_report_raw.html $RNA_DATA_MERGED

