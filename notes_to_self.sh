conda create --name org_pass --clone ibd_dnam_epic

## on server
conda env create -f org_pass.yml

conda activate org_pass

conda install -c bioconda bioconductor-bumphunter
BiocManager::install("rtracklayer") # and update none
BiocManager::install("GEOmetadb") # and update none
install.packages("here")


## maybe I dont need any of this for organoids?


					conda install -c bioconda bioconductor-methylumi bioconductor-watermelon 
					conda install -c conda-forge r-rpmm r-ggsci r-rcolorbrewer r-lme4 r-gridextra
					conda install -c bioconda bioconductor-minfi 
					conda install -c bioconda bioconductor-illuminahumanmethylation450kmanifest 
					conda install -c conda-forge r-rafalib 
					conda install -c r r-tidyverse
					conda install -c bioconda bioconductor-quantro
					conda install -c r r-gridbase 
					conda install -c conda-forge r-pvclust
					conda install -c bioconda bioconductor-geoquery 
					conda install -c bioconda bioconductor-geometadb 
					conda install -c conda-forge r-statmod 
					conda install -c bioconda bioconductor-arrayexpress 
					conda install -c omnia quadprog 
					conda install -c r r-car 
					conda install -c conda-forge r-lmtest 
					conda install -c conda-forge r-mixtools 
					conda install -c conda-forge numpy 
					conda install -c conda-forge statsmodels 
					conda install -c conda-forge r-metrics 
					conda install -c conda-forge r-ggridges 
					conda install -c bioconda bioconductor-watermelon 
					conda install -c conda-forge r-vgam
					conda install -c r r-xtable
					conda install -c conda-forge r-binom
					conda install -c conda-forge r-cowplot




# conda install -c bioconda bioconductor-illuminahumanmethylation450kmanifest bioconductor-illuminahumanmethylationepicmanifest 
BiocManager::install("IlluminaHumanMethylation450kmanifest", version = "3.8")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19", version = "3.8")
install.packages("inferr")

awk -F "\"*,\"*" 'NR>1{print $2}' ../../../../../organoid_passage_DNAm/data/original_organoids/epic_organoid_MethylationArraySamples_31Jul19.csv | while read samplename 
do
	cp ${samplename}* /nfs/research1/zerbino/redgar/organoid_passage_DNAm/data/original_organoids/
done








## sample selection


library(here)

load(here("../ibd/data/","raw/all_methylation_share/AllMethylationArraySamples_31Jul19_Phenotypes.RData"))
head(epic.organoid)

epic.organoid$sample_ID<-paste(epic.organoid$case.no, epic.organoid$sample.site)
epic.organoid$passage.or.rescope.no_numeric<-as.factor(epic.organoid$passage.or.rescope.no)
levels(epic.organoid$passage.or.rescope.no_numeric)<-c(1,10,11,14,16,2,3,4,6,7,8,2,4)
epic.organoid$passage.or.rescope.no_numeric<-as.numeric(as.character(epic.organoid$passage.or.rescope.no_numeric))

IBD_samples<-epic.organoid[which(!(epic.organoid$diagnosis%in%c("Control","Other.GI"))),]
IBD_samples_paired<-IBD_samples[which(IBD_samples$sample_ID%in%IBD_samples$sample_ID[duplicated(IBD_samples$sample_ID)]),]

epic.organoid_ctrl<-epic.organoid[which(epic.organoid$diagnosis%in%c("Control","Other.GI")),]


#Sample IDs I can pull 20 IBD samples

ibd_select<-c("T091 SC", "T049 SC", "T049 TI", "374 TI", "374 TI", "268 TI", "T046 SC", "T046 TI", "T088 TI", "279 SC", "T095 TI", "T013 SC")#T091 SC - 3, T049 - 4, T046 - 4, 374 TI - 2, 268 TI - 2, "T088 TI - 1  279 SC - 1 ,T095 TI -2 T013 SC - 1 20

epic.organoid_samples<-epic.organoid[which(epic.organoid$diagnosis%in%c("Control","Other.GI") | epic.organoid$sample_ID%in%ibd_select),]

table(epic.organoid_samples$diagnosis, epic.organoid_samples$passage.or.rescope.no_numeric)
table(epic.organoid_samples$passage.or.rescope.no_numeric)
epic.organoid_samples_paired<-epic.organoid_samples[which(epic.organoid_samples$sample_ID%in%epic.organoid_samples$sample_ID[duplicated(epic.organoid_samples$sample_ID)]),]

epic.organoid_samples$sentrix_ID<-gsub("_", "", epic.organoid_samples$sentrix.id)

epic.organoid_samples<-epic.organoid_samples[,c(1:8, 14,15, 32, 18, 29, 30)]
epic.organoid_samples$diagnosis<-as.factor(epic.organoid_samples$diagnosis)
levels(epic.organoid_samples$diagnosis)<-c("IBD","Control","Other.GI","IBD")

write.csv(epic.organoid_samples, file="data/original_organoids/epic_organoid_MethylationArraySamples_31Jul19.csv", row.names=F, quote=F)




load(here("../ibd/data/","ibd_beta_organoids.RData"))
epic.organoid_samples<-epic.organoid_samples[which(epic.organoid_samples$array.id%in%epic.organoid$array.id),]


table(epic.organoid_samples$diagnosis)
table(epic.organoid_samples$passage.or.rescope.no_numeric)
table(epic.organoid_samples$diagnosis, epic.organoid_samples$passage.or.rescope.no_numeric)

