### EPIC annotation with Reg feature
import pandas as pd
from numpy import genfromtxt
from itertools import chain
import sys
from collections import Counter
import functools

#The regulatory build (https://europepmc.org/articles/PMC4407537 http://grch37.ensembl.org/info/genome/funcgen/regulatory_build.html) was downloaded using biomart
Feature_bed = 'data/human_regulatory_features_GRCh37p13.txt'
background = 'data/passage_background.csv'

Feature_bed = pd.read_csv(Feature_bed, header=None, names=['chr','start','end','Feature'],skiprows=1)
CpG_background_pd = pd.read_csv(background)


CpG_start = int(sys.argv[1])
CpG_end = int(sys.argv[1])+10000

# subset to system arguments
CpG_background_pd = CpG_background_pd[CpG_start:CpG_end]


# make merge object to fill missing TFs to 0 
features = Feature_bed
features['count'] = 0
features = features[['Feature','count']]
features = pd.DataFrame.drop_duplicates(features)

### Actual background SNPs in TFBS
EPIC_features=[]

for snp, row in CpG_background_pd.iterrows():
    # match chromosome
    chr_features = Feature_bed.loc[Feature_bed['chr'] == str(row['CHR'])]

    features_snp_in = []
    for feat in range(1,len(chr_features.index)):#chr_TFBS.size
        if chr_features['start'].iloc[feat] <= row['MAPINFO'] <= chr_features['end'].iloc[feat]:
            features_snp_in.append(chr_features['Feature'].iloc[feat])

    EPIC_features.append([row['IlmnID'], features_snp_in])

df = pd.DataFrame(EPIC_features)
df.to_csv("data/annotation_split/EPIC_Features" + sys.argv[1] + ".csv") 