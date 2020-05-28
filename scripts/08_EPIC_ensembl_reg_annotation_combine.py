import numpy as np
import statsmodels
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.stats.api as sms
import sys

EPIC_annot = []

for i in range(0,810000,10000): #for i in range(0,800000,10000):
    EPICnext = pd.read_csv("data/annotation_split/EPIC_Features" + str(i) + ".csv")
    EPIC_annot.append(EPICnext)


EPIC_annot = pd.concat(EPIC_annot)
EPIC_annot.to_csv("data/EPIC_ensembl_reg_annotation.csv") 
