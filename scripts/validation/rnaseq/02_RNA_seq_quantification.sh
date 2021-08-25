export RNA_DATA_DIR1=/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/run1_fastqs
export RNA_DATA_DIR2=/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/run2_fastqs

export RNA_REFS_DIR=/nfs/research1/zerbino/redgar/ibd/data/public_rna_seq/refs
export RNA_DATA_QUANT=/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/kallisto

#conda activate ibd_dnam_rnaseq

#####################
## Kallisto pipeline
#####################
#wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -P $RNA_REFS_DIR
######wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz -P $RNA_REFS_DIR


## Buiding the kallisto index from the human cDNA file
kallisto index -i $RNA_REFS_DIR/human_transcriptome_index.idx $RNA_REFS_DIR/Homo_sapiens.GRCh38.cdna.all.fa.gz 


#cd $RNA_DATA_TRIM
#mkdir $RNA_DATA_TRIM/kallisto
cd $RNA_DATA_QUANT

awk 'NR>1{ print $2 }' ../sample_info_RNA_seq.txt | while read samplename 
do
	kallisto quant -b 100 -t 8 -i $RNA_REFS_DIR/human_transcriptome_index.idx -o ${samplename} $RNA_DATA_DIR1/${samplename}_R1_001.fastq.gz $RNA_DATA_DIR2/${samplename}_R1_001.fastq.gz
done



					# #cd $RNA_DATA_TRIM
					# #mkdir $RNA_DATA_TRIM/kallisto
					# cd $RNA_DATA_TRIM/kallisto

					# awk 'NR>1{ print $3 }' ../../Zerbino_RNASeq_SampleInfo.txt | while read samplename 
					# do
					# 	kallisto quant -b 100 -t 8 -i ../../refs/human_transcriptome_index.idx -o ${samplename} ../${samplename}_adapter_trimmed_1.fastq.gz ../${samplename}_adapter_trimmed_2.fastq.gz
					# done

