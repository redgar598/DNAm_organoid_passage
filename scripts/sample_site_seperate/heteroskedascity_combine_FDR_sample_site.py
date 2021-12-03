import numpy as np
import statsmodels
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.stats.api as sms
import sys

pvalues10 = pd.read_csv("../../../../../output/Heteroskedactiy_pvalues_FDR" +  sys.argv[1] + "1.csv") 

#pvalues10 = pd.read_csv("../../../../output/Heteroskedactiy_pvalues_FDR1.csv") 

pval_BP = pd.DataFrame(pvalues10["0"])
pval_diff = pd.DataFrame(pvalues10["1"])
mean_db = pd.DataFrame(pvalues10["2"])
fdr_BP = pd.DataFrame(pvalues10["3"])
fdr_diff = pd.DataFrame(pvalues10["4"])

for i in range(11,1000,10): #
    pvalues = pd.read_csv("../../../../../output/Heteroskedactiy_pvalues_FDR" +  sys.argv[1] + str(i) + ".csv") 
    pval_BP["" + str(i)] = pd.to_numeric(pvalues["0"])
    pval_diff["" + str(i)] = pd.to_numeric(pvalues["1"])
    mean_db["" + str(i)] = pd.DataFrame(pvalues["2"])
    fdr_BP["" + str(i)] = pd.DataFrame(pvalues["3"])
    fdr_diff["" + str(i)] = pd.DataFrame(pvalues["4"])



pval_BP = pval_BP.sum(axis = 1, skipna = True) 
pval_diff = pval_diff.sum(axis = 1, skipna = True) 
mean_db = mean_db.mean(axis = 1, skipna = True) 
fdr_BP = fdr_BP.sum(axis = 1, skipna = True) 
fdr_diff = fdr_diff.sum(axis = 1, skipna = True) 


pval_BP_df = pd.DataFrame([pval_BP, pval_diff,mean_db,fdr_BP,fdr_diff])
pval_BP_df = pval_BP_df.transpose() 
#pval_BP_df.to_csv("../../../../output/Heteroskedactiy_pvalues_FDR_1000iter.csv") 
pval_BP_df.to_csv("../../../../../output/Heteroskedactiy_pvalues_FDR_1000iter" + sys.argv[1] + ".csv") 