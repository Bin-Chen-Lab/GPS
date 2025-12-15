# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 21:04:58 2024

@author: xingjin1
"""
import pandas as pd
import numpy as np
import os
from scipy.stats import ranksums, spearmanr
import matplotlib
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
import seaborn as sns


wrkdir = './data_for_code/GSAR/'
rnaseq_dir = wrkdir + 'rna-seq/'
cpd_sigs = pd.read_csv(wrkdir + 'HCC_cmpds_gps_2types.csv', index_col=0)
cpd_info = pd.read_csv(wrkdir + 'DataSx_HCC_SAR_2types_cluster.csv', index_col=0)
cpd_info = cpd_info.loc[cpd_info['Good_pred'] == 'Yes']
cpd_info = cpd_info.sort_values(by=['Chemical Type', 'pIC50 Summary'], ascending=[True, False])
ovlp = [cpd for cpd in cpd_info.index if cpd in cpd_sigs.columns]
cpd_info = cpd_info.loc[ovlp]
cpd_sigs = cpd_sigs.loc[:, ovlp]
'''
keep = cpd_sigs.where(cpd_sigs != 0).count(axis=1)
keep = keep.loc[keep >= 2].index
cpd_sigs = cpd_sigs.loc[keep]
'''
print(cpd_sigs.shape)

tst = ['PB56874852', '5260420', '45302']
#tst += ['44443110', '45300']
trn = [cpd for cpd in cpd_info.index if not cpd in tst]
cpd_sigs_trn = cpd_sigs.loc[:, trn]
eff_trn = cpd_info.loc[trn, 'pIC50 Summary']

result = []
for gene, row in cpd_sigs_trn.iterrows():
    ct = row.value_counts()
    es, p, d, good = np.nan, np.nan, 0, False
    if len(ct) == 3:    
        es, p = spearmanr(row, eff_trn, nan_policy='omit')
        if p < 0.005:
        #if round(p, 2) <= 0.01:
            d = 1 if es > 0 else -1
            good = True
    elif len(ct) == 2:
        t = row.loc[row == ct.index.min()].index
        a = eff_trn.loc[t]
        t = row.loc[row == ct.index.max()].index
        b = eff_trn.loc[t]
        es, p = ranksums(a, b, nan_policy='omit')
        #if round(p, 2) <= 0.01:
        if p < 0.005:
            d = 1 if es < 0 else -1
            good = True
    else: continue
    result.append([gene, len(ct), es, p, good, d])
        
result = pd.DataFrame(data=result, columns=['Gene', 'Cats', 'ES', 'P', 'Selected', 'Direction'])
result.index = result['Gene']
result = result.iloc[:, 1:]
result = result.sort_values(by=['P'])
sele = result.loc[result['Selected'] == True].index.to_list()
print(len(sele))

viz = []
file_list = ['DE_MSU45302.csv', 'DE_4852_10uM.csv', 'DE_4852_4uM.csv']
for fn in file_list:
    treat = fn.replace('DE_', '').replace('.csv', '').replace('_', ' ')
    if 'Niclosamide' in treat: continue
    if '5260420' in treat: continue
    de = pd.read_csv(rnaseq_dir + fn, index_col='Symbol')
    de = de.loc[~de.index.isna()]
    de['abs'] = de['log2FoldChange'].abs()
    de = de.sort_values(by='abs', ascending=False)
    de = de.loc[~de.index.duplicated()]
    #de['log2FoldChange'] = [round(lfc, 2) if p < 0.2 else 0 for lfc, p in zip(de['log2FoldChange'], de['padj'])]
    
    cal = result.loc[(result.index.isin(de.index)) & (result['Selected'] == True), ['Direction']]
    cal['log2FoldChange'] = de.loc[cal.index, 'log2FoldChange']
    cal['decomposition'] = cal['log2FoldChange'] * cal['Direction']
    score = cal['decomposition'].sum() / len(sele)
    
    #tmp = cal['log2FoldChange'].reindex(sele).fillna(0).round(2)
    tmp = de.reindex(sele)[['pvalue', 'log2FoldChange']]
    t = []
    for p, lfc in zip(tmp['pvalue'], tmp['log2FoldChange']):
        if np.isnan(p):
            t.append(0)
        elif lfc > 0:
            t.append(-np.log10(p))
        else:
            t.append(np.log10(p))
    tmp[treat] = t
    viz.append(tmp[treat].round(2))
    
    print(treat, score)
viz = pd.concat(viz, axis=1)
viz = viz.sort_values(by='MSU45302')
#viz['Required'] = np.log10(result.loc[viz.index, 'P'])

FIG = plt.figure(figsize=(3, 6))
sns.heatmap(viz, vmin=-2, center=0, vmax=2, cmap='coolwarm', cbar=False, annot=True, linewidths=1, linecolor='grey')
plt.ylabel('Efficacy responsive genes', fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.title('Direction & Significance\n(RNA-seq)')
plt.tight_layout()
plt.savefig(wrkdir + 'Heatmap_GSAR-sig_302-4852.png', transparent=True)


gm = cpd_sigs.loc[sele, cpd_info.index].T
mp = {-1: 'Down', 0: 'None', 1: 'Up'}
for col in gm.columns:
    if col == 'pIC50 Summary': continue
    gm[col] = [mp[r] for r in gm[col]]
gm = pd.concat([cpd_info['pIC50 Summary'], gm], axis=1)
gm = gm.loc[gm['pIC50 Summary'] > 0]


genes2show = ['WDR75', 'KIF23', 'UHRF1', 'MCM6']
clr_code = {'Down': 'cornflowerblue', 'None': 'silver', 'Up': 'Salmon'}
for g in genes2show:
    gm = gm.sort_values(by=g)
    fig = plt.figure(figsize=(2.4, 2.4))
    sns.boxplot(data=gm, x=g, y='pIC50 Summary', palette=clr_code)
    sns.swarmplot(data=gm, x=g, y='pIC50 Summary', color='black')
    #plt.ylabel(None)
    plt.xlabel(None)
    plt.title(g + ': P = %.4f'%result.loc[g, 'P'])
    plt.tight_layout()
    plt.savefig(wrkdir + 'pIC50_separation_%s.png'%g, transparent=True)


de = pd.read_csv('../data/GSAR/rna-seq/DE_MSU45302.csv', index_col='Symbol')
de = de.loc[~de.index.isna()]
de['abs'] = de['log2FoldChange'].abs()
de = de.sort_values(by='abs', ascending=False)
c302 = de.loc[~de.index.duplicated()]
de = pd.read_csv('../data/GSAR/rna-seq/DE_4852_10uM.csv', index_col='Symbol')
de = de.loc[~de.index.isna()]
de['abs'] = de['log2FoldChange'].abs()
de = de.sort_values(by='abs', ascending=False)
c4852 = de.loc[~de.index.duplicated()]
com = pd.concat([c302[['log2FoldChange', 'padj']], c4852[['log2FoldChange', 'padj']]], axis=1)
com.columns = ['log2FoldChange_MSU45302', 'padj_MSU45302', 'log2FoldChange_4852', 'padj_4852']
com = com.reindex(cpd_sigs.index).dropna()
tmp, thr = [], 0.1
for a, b in zip(com['padj_MSU45302'], com['padj_4852']):
    tmp.append(sum([i < thr for i in [a, b]]))
com['Significant'] = tmp
com = com.sort_values(by='log2FoldChange_MSU45302', ascending=False)

clr_code = {0: 'grey', 1: 'salmon', 2: 'chocolate'}
fig = plt.figure(figsize=(4, 4))
sns.scatterplot(data=com, x='log2FoldChange_MSU45302', y='log2FoldChange_4852', \
                hue='Significant', palette=clr_code, alpha=0.75, size='Significant')
for key, row in com.loc[com['Significant'] == 2].iterrows():
    a, b = row['log2FoldChange_MSU45302'], row['log2FoldChange_4852']
    plt.text(a, b, key)
plt.xlim(-1, 1)
plt.ylim(-3, 3)
plt.show()

    




