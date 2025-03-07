# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 15:08:21 2024

@author: xingjin1
"""
import pandas as pd
from sklearn.manifold import TSNE
import matplotlib
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
import seaborn as sns


wrkdir = '../data/GSAR/'
cpd_sigs = pd.read_csv(wrkdir + 'HCC_cmpds_gps_2types.csv', index_col=0)
cpd_info = pd.read_csv(wrkdir + 'HCC_SAR_2types_cluster.csv', index_col=0)
cpd_info = cpd_info.loc[cpd_info['pIC50 Summary'] > 0]
cpd_info = cpd_info.loc[cpd_info['Good_pred'] == 'Yes']
cpd_info = cpd_info.sort_values(by=['Chemical Type', 'pIC50 Summary'], ascending=[True, False])
ovlp = [cpd for cpd in cpd_info.index if cpd in cpd_sigs.columns]
cpd_info = cpd_info.loc[ovlp]
cpd_sigs = cpd_sigs.loc[:, ovlp]
keep = cpd_sigs.where(cpd_sigs != 0).count(axis=1)
keep = keep.loc[keep >= 1].index
cpd_sigs = cpd_sigs.loc[keep]
print(cpd_sigs.shape)

model = TSNE(random_state=1, perplexity=6, metric='euclidean')
embd = model.fit_transform(cpd_sigs.T)
cpd_info['x'] = embd[:, 0]
cpd_info['y'] = embd[:, 1]

FIG = plt.figure(figsize=(3.5, 3))
sns.scatterplot(x='x', y='y', hue='pIC50 Summary', style='Chemical Type', data=cpd_info, palette='coolwarm')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=8)
'''
for cpd in ['45302', 'PB56874852', 'Niclosamide', '5260420']:
    a, b  = cpd_info.loc[cpd, ['x', 'y']]
    plt.text(a, b, cpd, ha='left')
'''
#plt.axis('off')
plt.tight_layout()
#plt.savefig(wrkdir + 'TSNE_HCC_cmpds_GPS_chemtypes.png', transparent=True)






