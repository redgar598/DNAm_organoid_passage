# Organoid Passage DNAm
Analysis of the effect of passage on organoid DNAm

# Workflow
The scripts were run in the order below (also in [organoid_workflow.sh](https://github.com/redgar598/organoid_passage_DNAm/tree/master/scripts/organoid_workflow.sh))

| Script                                                                                                                                                          | Action                                                                                                                                                                                   | Figures in Paper |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------|
| [01_organoids_preprocessing.R](https://github.com/redgar598/organoid_passage_DNAm/tree/master/output/01_organoids_preprocessing.html)                           | Loads IDAT files normalizes, probe QC and sample QC                                                                                                                                      |                  |
| [02_heteroskedasicity_CpGs.py](https://github.com/redgar598/organoid_passage_DNAm/tree/master/output/02_heteroskedasicity_CpGs.ipynb)                           | Subsampling of the low passage samples in 1,000 iterations. Calculates differential DNAm and heteroskedastic p values                                                                    |                  |
| [03_heteroskedascity_CpGs_combine.py](https://github.com/redgar598/organoid_passage_DNAm/tree/master/output/03_heteroskedascity_CpGs_combine.ipynb)             | Combines differential DNAm and heteroskedastic p values from each iteration into one table                                                                                               |                  |
| [04_passage_CpGs.R](https://github.com/redgar598/organoid_passage_DNAm/tree/master/output/04_passage_CpGs.html)                                                 | Plots representative passage CpGs, makes the variable CpG beta distributions                                                                                                             |                  |
| [05_EPIC_annotation_cDMR.py](https://github.com/redgar598/organoid_passage_DNAm/tree/master/output/05_EPIC_annotation_cDMR.ipynb)                               | Annotates all CpGs on the EPIC array with [cDMR](https://www.nature.com/articles/ng.298#Sec25) and direction of change of cDMR between cancer and healthy.                               |                  |
| [06_EPIC_annotation_cDMR_combine.py](https://github.com/redgar598/organoid_passage_DNAm/tree/master/output/06_EPIC_annotation_cDMR_combine.ipynb)               | The output of annotating the EPIC with cDMRs is split into 80 jobs each with 10,000 CpGs so these are then combined into on output.                                                      |                  |
| [07_EPIC_ensembl_reg_annotation.py](https://github.com/redgar598/organoid_passage_DNAm/tree/master/output/07_EPIC_ensembl_reg_annotation.ipynb)                 | Annotates all CpGs on the EPIC array with the [ensembl regulatory build](https://europepmc.org/articles/PMC4407537 http://grch37.ensembl.org/info/genome/funcgen/regulatory_build.html). |                  |
| [08_EPIC_ensembl_reg_annotation_combine.py](https://github.com/redgar598/organoid_passage_DNAm/tree/master/output/08_EPIC_ensembl_reg_annotation_combine.ipynb) | The output of annotating the EPIC with the regulatory build is split into 80 jobs each with 10,000 CpGs so these are then combined into on output.                                       |                  |
| [09_enrichment_heteroskedasticity_CpG.py](https://github.com/redgar598/organoid_passage_DNAm/tree/master/output/09_enrichment_heteroskedasticity_CpG.html)      | Calculate enrichment statistics for passage CpGs in cDMRs and regulatory features.                                                                                                       |                  |
|                                                                                                                                                                 |                                                                                                                                                                                          |                  |
|                                                                                                                                                                 |                                                                                                                                                                                          |                  |
|                                                                                                                                                                 |                                                                                                                                                                                          |                  |
|                                                                                                                                                                 |                                                                                                                                                                                          |                  |
|                                                                                                                                                                 |                                                                                                                                                                                          |                  |
|                                                                                                                                                                 |                                                                                                                                                                                          |                  |
|                                                                                                                                                                 |                                                                                                                                                                                          |                  |
|                                                                                                                                                                 |                                                                                                                                                                                          |                  |
|                                                                                                                                                                 |                                                                                                                                                                                          |                  |
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
