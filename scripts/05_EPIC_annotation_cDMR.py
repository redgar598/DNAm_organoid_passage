#### cDMR annotation of EPIC
import pandas as pd
from numpy import genfromtxt
from itertools import chain
import sys
from collections import Counter
import functools


# https://www.nature.com/articles/ng.298#Sec25
# supplementary data 2
Feature_bed = 'data/41588_2009_BFng298_MOESM18_ESM.csv'
background = 'data/passage_background_build36.csv'

Feature_bed = pd.read_csv(Feature_bed, header=1,skiprows=2)
CpG_background_pd = pd.read_csv(background)

Feature_bed['chr'] = Feature_bed['chr'].str.replace('chr','')

CpG_start = int(sys.argv[1])
CpG_end = int(sys.argv[1])+10000

# subset to system arguments
CpG_background_pd = CpG_background_pd[CpG_start:CpG_end]

### Actual background SNPs in TFBS
EPIC_features=[]

for snp, row in CpG_background_pd.iterrows():
    # match chromosome
    chr_features = Feature_bed.loc[Feature_bed['chr'] == str(row['Chromosome_36'])]

    features_snp_in = []
    for feat in range(0,len(chr_features.index)):#chr_TFBS.size
        if chr_features['start'].iloc[feat] <= int(row['Coordinate_36']) <= chr_features['end'].iloc[feat]:
            features_snp_in.append(chr_features['delta M'].iloc[feat])

    EPIC_features.append([row['IlmnID'], features_snp_in])

df = pd.DataFrame(EPIC_features)
df.to_csv("data/annotation_split/EPIC_Features_cDMR36_" + sys.argv[1] + ".csv") 