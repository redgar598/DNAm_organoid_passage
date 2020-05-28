# Organoid Passage DNAm
Analysis of the effect of passage on organoid DNAm

# Workflow
The scripts were run in the order below (also in [organoid_workflow.sh](https://github.com/redgar598/organoid_passage_DNAm/tree/master/scripts/organoid_workflow.sh))


### [01_organoids_preprocessing.R](https://github.com/redgar598/organoid_passage_DNAm/tree/master/output/01_organoids_preprocessing.html)
*  Loads IDAT files normalizes, probe QC and sample QC
### [02_heteroskedasicity_CpGs.py](https://github.com/redgar598/organoid_passage_DNAm/tree/master/output/02_heteroskedasicity_CpGs.ipynb)
*  Run in iterations of 1000 random samples. Does subsampling of the low passage samples and calculated differential DNAm and heteroskedastic p values
### [03_heteroskedascity_CpGs_combine.py](https://github.com/redgar598/organoid_passage_DNAm/tree/master/output/03_heteroskedascity_CpGs_combine.ipynb)
*  Combines differential DNAm and heteroskedastic p values from each iteration into one table
### [04_passage_CpGs.R](https://github.com/redgar598/organoid_passage_DNAm/tree/master/output/04_passage_CpGs.html)
*  Plots representative passage CpGs, makes the variable CpG beta distributions


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
