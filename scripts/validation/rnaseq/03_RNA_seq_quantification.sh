export RNA_DATA_DIR=../../data/public_rna_seq
export RNA_REFS_DIR=../../data/public_rna_seq/refs
export RNA_DATA_TRIM=../../data/public_rna_seq/trimmed

conda activate ibd_dnam_rnaseq

#####################
## Kallisto pipeline
#####################
#wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -P $RNA_REFS_DIR
######wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz -P $RNA_REFS_DIR


## Buiding the kallisto index from the human cDNA file
kallisto index -i $RNA_REFS_DIR/human_transcriptome_index.idx $RNA_REFS_DIR/Homo_sapiens.GRCh38.cdna.all.fa.gz 
#kallisto index -i human_transcriptome_index_GRCh37.idx Homo_sapiens.GRCh37.75.cdna.all.fa.gz

# ## Get the human GTF file
# wget ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz
# gunzip -v -f ./Homo_sapiens.GRCh38.95.gtf.gz


# paired end
#cd $RNA_DATA_DIR
#mkdir kallisto
#cd kallisto
#mkdir sleuth

# awk 'NR>1{ print $3 }' ../../Zerbino_RNASeq_SampleInfo.txt | while read samplename 
# do
# 	kallisto quant --fr-stranded -i ../../refs/human_transcriptome_index.idx -o ${samplename} ../${samplename}_adapter_trimmed_1.fastq.gz ../${samplename}_adapter_trimmed_2.fastq.gz
# done


#awk 'NR>1{ print $3 }' ../Zerbino_RNASeq_SampleInfo.txt | while read samplename 
#do
#	kallisto quant -b 10 -t 8 -i ../refs/human_transcriptome_index.idx -o ${samplename} ../${samplename}_1.fastq.gz ../${samplename}_2.fastq.gz
#done


													# #single test
													# kallisto quant -i refs/human_transcriptome_index.idx -o output --single -l 126 -s 20 ERR2271024_1.fastq.gz ERR2271030_1.fastq.gz

											# awk 'NR>1{ print $3 }' ../../Zerbino_RNASeq_SampleInfo.txt | while read samplename 
											# do
											# 	kallisto quant -i ../../refs/human_transcriptome_index.idx -o ${samplename} -b 100 ../${samplename}_adapter_trimmed_1.fastq.gz ../${samplename}_adapter_trimmed_2.fastq.gz
											# done


											# kallisto quant -i ../../refs/human_transcriptome_index.idx -o test -b 100 ../ERR2270963_adapter_trimmed_1.fastq.gz ../ERR2270963_adapter_trimmed_2.fastq.gz
											# kallisto quant -i ../../refs/human_transcriptome_index.idx -o test -b 100 ../../ERR2270963_1.fastq.gz ../../ERR2270963_2.fastq.gz


### on cutadapt trimmed instead

#cd $RNA_DATA_TRIM
#mkdir $RNA_DATA_TRIM/kallisto
cd $RNA_DATA_TRIM/kallisto

awk 'NR>1{ print $3 }' ../../Zerbino_RNASeq_SampleInfo.txt | while read samplename 
do
	kallisto quant -b 100 -t 8 -i ../../refs/human_transcriptome_index.idx -o ${samplename} ../${samplename}_adapter_trimmed_1.fastq.gz ../${samplename}_adapter_trimmed_2.fastq.gz
done