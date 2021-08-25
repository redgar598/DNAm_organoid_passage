########### ibd_dnam (RNA-seq)
conda create --name ibd_dnam_rnaseq

conda activate ibd_dnam_rnaseq

conda install -c bioconda -c conda-forge kallisto fastqc multiqc

# R
conda install -c r r-devtools
conda install -c conda-forge r-gplots 
conda install -c r r-dplyr r-ggplot2 
conda install -c bioconda r-sleuth 
conda install -c r r-testit
conda install -c r r-reshape 
conda install -c conda-forge r-gridextra 
conda install -c conda-forge r-lme4
conda install -c conda-forge r-rafalib 


# # to run rstudio on right version 
# export RSTUDIO_WHICH_R=/home/redgar/anaconda3/envs/ibd_dnam_rnaseq/bin/R

#conda install -c biobuilds fastx-toolkit 
conda install -c bioconda cutadapt

