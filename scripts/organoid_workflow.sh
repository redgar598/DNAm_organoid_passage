#### Organoid pipeline
conda activate org_pass
cd /nfs/research1/zerbino/redgar/DNAm_organoid_passage
git pull


bsub -M 80000 -R "rusage[mem=80000]" Rscript scripts/organoid_workflow_Rmarkdown.R
bsub -M 50000 -R "rusage[mem=50000]" Rscript scripts/testing.R






##############
## Normalize and sample quality control
##############
bsub -M 80000 -R "rusage[mem=80000]" Rscript scripts/01_organoids_preprocessing.R


##############
## Permutate the heteroskedacity and differential pvalue calculations
##############
for iter in $(seq 1 10 1000); do
	bsub -M 15000 -R "rusage[mem=15000]" python scripts/02_heteroskedasicity_CpGs.py $iter
done

#combine outputs into one file
bsub -M 25000 -R "rusage[mem=25000]" python scripts/03_heteroskedascity_CpGs_combine.py



##############
### Check out the CpGs
##############
bsub -M 60000 -R "rusage[mem=60000]" Rscript scripts/04_passage_CpGs.R





#############
# Annotate the EPIC with cDMR directionality
#############

# split into 80 jobs of 10,000 CpGs
for iter in $(seq 0 10000 800000); do
	bsub -M 15000 -R "rusage[mem=15000]" python scripts/05_EPIC_annotation_cDMR.py $iter
done

#combine 80 outputs into one file
bsub -M 25000 -R "rusage[mem=25000]" python scripts/06_EPIC_annotation_cDMR_combine.py





#############
# Annotate the EPIC with ensembl regulatory features
#############

# split into 80 jobs of 10,000 CpGs
for iter in $(seq 0 10000 800000); do
	bsub -M 15000 -R "rusage[mem=15000]" python scripts/07_EPIC_ensembl_reg_annotation.py $iter
done

#combine 80 outputs into one file
bsub -M 25000 -R "rusage[mem=25000]" python scripts/08_EPIC_ensembl_reg_annotation_combine.py




##############
### Run enrichment of CpGs over background
##############
bsub -M 55000 -R "rusage[mem=55000]" Rscript scripts/09_enrichment_heteroskedasticity_CpG.R



##############
### EPIC annotation CpG to gene
##############
bsub -M 15000 -R "rusage[mem=15000]" Rscript scripts/10_ensembl_genes_EPIC_annotation.R 


for iter in $(seq 1 10000 870000); do
	bsub -M 40000 -R "rusage[mem=40000]" Rscript scripts/11_EPIC_gene_annotation.R $iter
done

bsub -M 25000 -R "rusage[mem=25000]" python scripts/12_EPIC_gene_annotation_combine.py



##############
### GO Over Representation Analysis
##############
bsub -M 20000 -R "rusage[mem=20000]" Rscript scripts/13_prep_For_GO_ORA.R

############ go enrichment
sudo apt -y install icedtea-netx icedtea-plugin
javaws data/tools/ermineJ.jnlp
# Gene annotation file: Generic_human_ensemblIds_noParents_EPIC_array.an.txt
# Analysis -> Run Analysis -> ORA
# ORA Paste genes lists as "Quick list" (did not use a score)
# Gene lists to paste: hetero_genes.csv hyper_genes.csv hypo_genes.csv
# Biological process only
# max set size=400 min=5
# no log of score, larger scores are not better gene score threshold 0.0

# analysis-> save analysis -> chose each ORA run seperately -> include all genes in output (yes)
# save hyper.erminej.txt hypo.erminej.txt hetero.erminej.txt

bsub -M 20000 -R "rusage[mem=20000]" Rscript scripts/14_post_GO_ORA.R



##############
### Cancer beta distributions
##############
bsub -M 20000 -R "rusage[mem=20000]" Rscript scripts/15_cancer_beta_simplified.R



###############
### Validation in other datasets
##############
bsub -M 80000 -R "rusage[mem=80000]" Rscript scripts/16_MTAB4957_passage.R
bsub -M 80000 -R "rusage[mem=80000]" Rscript scripts/17_GSE141256_organoids.R



bsub -M 40000 -R "rusage[mem=40000]" Rscript scripts/18_GSE142213_organoids.R



######################
# rnaseq pipeline
######################
cd /nfs/research1/zerbino/redgar/DNAm_organoid_passage
conda activate ibd_dnam_rnaseq
git pull


bsub -M 10000 -R "rusage[mem=10000]" -n 8 bash scripts/validation/rnaseq/01_RNA_seq_processing.sh
bsub -M 10000 -R "rusage[mem=10000]" -n 8 bash scripts/validation/rnaseq/02_RNA_seq_runs_merge.sh
bsub -M 10000 -R "rusage[mem=10000]" -n 8 bash scripts/validation/rnaseq/04_RNA_seq_quantification.sh

# bsub -M 10000 -R "rusage[mem=10000]" -n 8 bash RNAseq/rna_seq_quantification.sh
# bsub -M 80000 -R "rusage[mem=40000]" Rscript RNAseq/sleuth_exploratory.R
# bsub -M 80000 -R "rusage[mem=20000]" Rscript RNAseq/rna_seq_countadjustment.R


# bsub -M 20000 -R "rusage[mem=20000]" Rscript RNAseq/testing.R

scp redgar@ebi:/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/run1_fastqs/multiqc_report_raw.html /home/redgar/Documents/DNAm_organoid_passage/data/validation_dataset
scp redgar@ebi:/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/run2_fastqs/multiqc_report_raw.html /home/redgar/Documents/DNAm_organoid_passage/data/validation_dataset

scp redgar@ebi:/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/run1_fastqs/T317_SC_p12_D_S15_R1_001_fastqc.html /home/redgar/Documents/DNAm_organoid_passage/data/validation_dataset
scp /home/redgar/Documents/DNAm_organoid_passage/data/validation_dataset/sample_info_RNA_seq.csv redgar@ebi:/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage

scp -r redgar@ebi:/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/kallisto /home/redgar/Documents/DNAm_organoid_passage/data/validation_dataset



				scp redgar@ebi:/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/run1_fastqs/T317_SC_p12_D_S15_R1_001* /home/redgar/Documents/DNAm_organoid_passage/data/validation_dataset/test

scp redgar@ebi:/nfs/research1/zerbino/redgar/ibd/data/raw/RNAseq_passage/merged/multiqc_report_raw.html /home/redgar/Documents/DNAm_organoid_passage/data/validation_dataset/merged

scp redgar@ebi:/nfs/research1/zerbino/redgar/DNAm_organoid_passage/data/validation/DNAm/*.csv /home/redgar/Documents/DNAm_organoid_passage/data/validation/DNAm
scp /home/redgar/Documents/DNAm_organoid_passage/data/validation/Final\ samples\ for\ submission\ 2021.csv redgar@ebi:/nfs/research1/zerbino/redgar/DNAm_organoid_passage/data/validation

scp redgar@ebi:/nfs/research1/zerbino/redgar/DNAm_organoid_passage/figs/validation_detection_pvalue_organoids.pdf /home/redgar/Documents/DNAm_organoid_passage/figs
scp redgar@ebi:/nfs/research1/zerbino/redgar/DNAm_organoid_passage/figs/jpeg/validation_detection_pvalue_organoids.jpeg /home/redgar/Documents/DNAm_organoid_passage/figs/jpeg


scp redgar@ebi:/nfs/research1/zerbino/redgar/DNAm_organoid_passage/output/01_validation_organoids.html /home/redgar/Documents/DNAm_organoid_passage/output


scp redgar@ebi:/nfs/research1/zerbino/redgar/DNAm_organoid_passage/output/data/validation/DNAm/validation_betas_normalized.RData /home/redgar/Documents/DNAm_organoid_passage/data/validation/DNAm


