# Organoid Passage DNAm
Analysis of the effect of passage on organoid DNAm 

<a href="https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/01_organoids_preprocessing.html" title="Precrocessing">About Me</a>

# Workflow
The scripts were run in the order below (also in [organoid_workflow.sh](https://github.com/redgar598/DNAm_organoid_passage/tree/master/scripts/organoid_workflow.sh))


| Script                                                                                                                                                          | Action                                                                                                                                                     | Figures in Paper |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------|
| [01_organoids_preprocessing.R](https://htmlpreview.github.io/?https://github.com/redgar598/DNAm_organoid_passage/blob/master/output/01_organoids_preprocessing.html)                           | Loads IDAT files normalizes, probe QC and sample QC                                                                                                        |       1, S1          |
| [02_heteroskedasicity_CpGs.py](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/02_heteroskedasicity_CpGs.ipynb)                           | Subsampling of the low passage samples in 1,000 iterations. Calculates differential DNAm and heteroskedastic p values                                      |        -          |
| [03_heteroskedascity_CpGs_combine.py](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/03_heteroskedascity_CpGs_combine.ipynb)             | Combines differential DNAm and heteroskedastic p values from each iteration into one table                                                                 |         -         |
| [04_passage_CpGs.R](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/04_passage_CpGs.html)                                                 | Plots representative passage CpGs, makes the variable CpG beta distributions                                                                               |          2, 3a, S2,S3        |
| [05_EPIC_annotation_cDMR.py](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/05_EPIC_annotation_cDMR.ipynb)                               | Annotates all CpGs on the EPIC array with [cDMR](https://www.nature.com/articles/ng.298#Sec25) and direction of change of cDMR between cancer and healthy. |        -          |
| [06_EPIC_annotation_cDMR_combine.py](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/06_EPIC_annotation_cDMR_combine.ipynb)               | The output of annotating the EPIC with cDMRs is split into 80 jobs each with 10,000 CpGs so these are then combined into one output.                       |        -          |
| [07_EPIC_ensembl_reg_annotation.py](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/07_EPIC_ensembl_reg_annotation.ipynb)                 | Annotates all CpGs on the EPIC array with the [ensembl regulatory build](http://grch37.ensembl.org/info/genome/funcgen/regulatory_build.html).             |         -         |
| [08_EPIC_ensembl_reg_annotation_combine.py](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/08_EPIC_ensembl_reg_annotation_combine.ipynb) | The output of annotating the EPIC with the regulatory build is split into 80 jobs each with 10,000 CpGs so these are then combined into one output.        |        -          |
| [09_enrichment_heteroskedasticity_CpG.R](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/09_enrichment_heteroskedasticity_CpG.html)       | Calculate enrichment statistics for passage CpGs in cDMRs and regulatory features.                                                                         |        3bc, S8         |
| [10_ensembl_genes_EPIC_annotation.R](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/09_enrichment_heteroskedasticity_CpG.html)           | Annotate human genes with promoter, gene body and 3’UTR regions for later annotation of CpGs on EPIC into these regions.                                   |       -           |
| [11_EPIC_gene_annotation.R](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/11_EPIC_gene_annotation.html)                                 | Annotates all CpGs on the EPIC array with which gene and gene feature they are in.                                                                         |        -          |
| [12_EPIC_gene_annotation_combine.py](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/12_EPIC_gene_annotation_combine.ipynb)               | The output of annotating the EPIC with the gene features into 80 jobs each with 10,000 CpGs so these are then combined into one output.                    |         -         |
| [13_prep_For_GO_ORA.R](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/13_prep_For_GO_ORA.html)                                           | Annotate the genes on the EPIC array with [GO sene sets](https://gemma.msl.ubc.ca/annots/Generic_human_ensemblIds_noParents.an.txt.gz) for a background gene list to compared to the genes associated with passage.                                                                                                                                                         |         -         |
| [14_post_GO_ORA.R](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/13_prep_For_GO_ORA.html)                                               | Tidy the output of erminej to include in the supplement of enriched GO groups.                                                                                                                                                         |         Tables S1-3         |
| [15_cancer_beta_simplified.R](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/15_cancer_beta_simplified.html)                             | View the distributions of the top variable CpGs in primary colorectal cancer samples [GSE42752](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) and [GSE48684](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48684). |        S9          |
| [16_MTAB4957_passage.R](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/16_MTAB4957_passage.html)                                         | Validate passage effect in pediatric and fetal organoids from [E-MTAB-4957](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-4957/) |         S4, S5, S7         |
| [17_GSE141256_organoids.R](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/17_GSE141256_organoids.html)                                   | Validate passage effect in organoid samples from [GSE141256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141256)  |        S6          |
| [18_GSE142213_organoids.R](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/18_GSE142213_organoids.html)                                   | Explore the trimodality of cancer organoids [GSE142213](https://www.ncbi.nlm.nih.gov/gds/?term=GSE142213[Accession]) |        S10          |



## Functions imported throughout workflow above

| Script                                                                                                                                                                            | Functions                                                                                     |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| [00_EM_array_uniform_background_maximise_betabinom.R](https://github.com/redgar598/organoid_passage_DNAm/tree/master/scripts/00_EM_array_uniform_background_maximise_betabinom.R) | Expectation Maximization fitting of beta-binomial distributions to methylation distribution   |
| [00_heat_scree_plot_generic.R](https://github.com/redgar598/organoid_passage_DNAm/tree/master/scripts/00_heat_scree_plot_generic.R)                                               | Plot the associations between PCA loadings and meta data                                      |
| [00_pretty_plots.R](https://github.com/redgar598/organoid_passage_DNAm/tree/master/scripts/00_pretty_plots.R)                                                                     | Colour palettes and themes for all plots                                                      |


| Script                                                                                                                                                                                                                                                           | Figures in Paper |
|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------|
| [01_organoids_preprocessing.R](https://htmlpreview.github.io/?https://github.com/redgar598/DNAm_organoid_passage/blob/master/output/01_organoids_preprocessing.html) <br /> Loads IDAT files normalizes, probe QC and sample QC                                    | 1, S1            |
| [02_heteroskedasicity_CpGs.py](https://github.com/redgar598/DNAm_organoid_passage/tree/master/output/02_heteroskedasicity_CpGs.ipynb) <br /> Subsampling of the low passage samples in 1,000 iterations. Calculates differential DNAm and heteroskedastic p values | -                |
----



# Environment used for analysis
The envrionment used is org_pass and can be created through the environment from yml file.
```
conda env create -f org_pass.yml
```

Unfortunately some key packages did not get exported to yml. So they can be installed through conda and R.
```
conda activate org_pass

conda install -c bioconda bioconductor-bumphunter
BiocManager::install("rtracklayer") # and update none
BiocManager::install("GEOmetadb") # and update none
install.packages("here")
```

For these packages  used versions

bumphunter_1.24.5

rtracklayer_1.42.1

GEOmetadb_1.44.0


Also I installed "here" but just to manage the pathways for rendering of these markdown documents.
```
install.packages("here")
```
Version here_0.1

The data folder provided on github does not include data. All files used include links to where the data is deposited pubically.
