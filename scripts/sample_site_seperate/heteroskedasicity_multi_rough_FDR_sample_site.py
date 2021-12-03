
# coding: utf-8

# In[196]:


import numpy as np
import statsmodels
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.stats.api as sms
import sys
import statistics 


# In[197]:
beta = pd.read_csv('data/beta_organoids.csv')


# In[198]:
meta = pd.read_csv('data/meta_organoids.csv')


# In[7]:


# prepare passage column from linear modelling
meta.rename(columns={"passage.or.rescope.no": "passage", "sample.type": "sampletype"}, inplace=True)

#df['score_num'] = df['score'].apply(score_to_numeric)
meta['passage'] = meta['passage'].str.replace('P','')
meta['passage'] = meta['passage'].str.replace('RE1.','')
meta['passage'] = pd.to_numeric(meta['passage'])

meta = meta[meta['sample.site'] == sys.argv[2]]

# In[192]:


import random

permstart = int(sys.argv[1])
permend = int(sys.argv[1])+10
CpGnum = beta.shape[0]


pval_all_BP = []
pval_all_diff = []
db_all_diff = []
fdr_all_BP = []
fdr_all_diff = []



for n in range(permstart,permend):
    
    random.seed(n)


    ## Sample the cohort in the lower passage number samples. 
    #Pull 3 random samples from passage 2 as it is over represented

    meta_sampled_high_passage = meta[meta['passage'] != 2]

    meta_sampled = meta[meta['passage'] == 2]

    meta_sampled_grouped = meta_sampled.groupby('passage')

    meta_sampled_subset = []
    for name, group in meta_sampled_grouped:
        meta_sampled_subset.append(group.sample(3))

    meta_sampled_subset = pd.concat([pd.concat(meta_sampled_subset),meta_sampled_high_passage])


    ## collect a p value for each CpG

    beta_sampled = beta[meta_sampled_subset['array.id'].values.tolist()]

    CpG_pval_passage_subset = []
    CpG_pval_BP_subset = []
    CpG_db_passage_subset = []


    for cpg in range(0, CpGnum): #beta_sampled.shape[0]

        meta_sampled_subset['beta'] = beta_sampled.iloc[cpg,0:45].values.tolist()
        meta_sampled_subset['constant'] = 1

        reg = smf.ols('beta ~ passage', data=meta_sampled_subset).fit()
        # Differential p value is interesting as well
        pval_passage = reg.pvalues[1]
        db = (reg.params[1]*1)-(reg.params[1]*16)

        pred_val = reg.fittedvalues.copy()
        true_val = meta_sampled_subset['beta'].values.copy()
        residual = true_val - pred_val

        #BP heteroskedacity test
        _, pval_BP, __, f_pval = statsmodels.stats.diagnostic.het_breuschpagan(residual, meta_sampled_subset[['passage','constant']])
        # studentized or not (p vs f) values do match the ones from bptest in R

        CpG_pval_BP_subset.append(pval_BP)
        CpG_pval_passage_subset.append(pval_passage)
        CpG_db_passage_subset.append(db)

        
    pval_all_BP.append(CpG_pval_BP_subset)
    pval_all_diff.append(CpG_pval_passage_subset)
    db_all_diff.append(CpG_db_passage_subset)
    fdr_all_BP.append(statsmodels.stats.multitest.multipletests(CpG_pval_BP_subset, method='fdr_bh', is_sorted=False, returnsorted=False)[1])
    fdr_all_diff.append(statsmodels.stats.multitest.multipletests(CpG_pval_passage_subset, method='fdr_bh', is_sorted=False, returnsorted=False)[1])
   



# In[193]:


pval_BP_df = pd.DataFrame(pval_all_BP)
pval_diff_df = pd.DataFrame(pval_all_diff)
db_all_diff = pd.DataFrame(db_all_diff)
fdr_all_BP = pd.DataFrame(fdr_all_BP)
fdr_all_diff = pd.DataFrame(fdr_all_diff)


# In[194]:


sig_BP = []
for cpg in range(0, CpGnum): #beta_sampled.shape[0]
    sig = sum(pval_BP_df.iloc[:,cpg] < 0.05)
    sig_BP.append(sig)
    
sig_diff = []
for cpg in range(0, CpGnum): #beta_sampled.shape[0]
    sig = sum(pval_diff_df.iloc[:,cpg] < 0.05)
    sig_diff.append(sig)

mn_db = []
for cpg in range(0, CpGnum): #beta_sampled.shape[0]
    mn = statistics.mean(db_all_diff.iloc[:,cpg])
    mn_db.append(mn)
    
sig_BP_fdr = []
for cpg in range(0, CpGnum): #beta_sampled.shape[0]
    sig = sum(fdr_all_BP.iloc[:,cpg] < 0.05)
    sig_BP_fdr.append(sig)

sig_diff_fdr = []
for cpg in range(0, CpGnum): #beta_sampled.shape[0]
    sig = sum(fdr_all_diff.iloc[:,cpg] < 0.05)
    sig_diff_fdr.append(sig)
    
    
pval_BP_df = pd.DataFrame([sig_BP, sig_diff,mn_db, sig_BP_fdr,sig_diff_fdr])

pval_BP_df = pval_BP_df.transpose() 
pval_BP_df.to_csv("data/passage_CpG_iterations/Heteroskedactiy_pvalues_FDR"+  sys.argv[2] + sys.argv[1] + ".csv") 
#pval_BP_df.to_csv("../../../../output/Heteroskedactiy_pvalues_iter.csv") 
