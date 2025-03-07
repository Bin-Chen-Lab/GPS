# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 15:16:39 2023

@author: Jing
"""

import pandas as pd
import os
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["font.family"] = "Arial"


cpd_sele = pd.read_csv('../data/cmpd_profile_robustness/selected_multiCell_cmpds-1.csv', index_col='StructID')
profiles_dir = '../data/cmpd_profile_robustness/LINCS2020_Level5_gold_cmpd_24h_10uM/'
cell_sele = ['PC3', 'MCF7', 'VCAP', 'A549', 'A375', 'HA1E', 'HT29', 'NPC', \
             'HCC515', 'ASC', 'HELA', 'SKB', 'YAPC', 'HEPG2', 'MCF10A', 'NEU', \
             'PHH', 'MDAMB231', 'XC.L10', 'THP1']

sele_profiles = []
for fn in os.listdir(profiles_dir):
    if not fn.endswith('.csv'): continue
    cell = fn.replace('.csv', '').split('_')[-1]
    if not cell in cell_sele: continue
    tmp = pd.read_csv(profiles_dir + fn, index_col=0)
    tmp = tmp.loc[:, tmp.columns.isin(cpd_sele.index)].T
    print(cell, tmp.shape)
    if len(tmp) == 0: continue
    tmp['Compound'] = tmp.index
    tmp['Cell line'] = [cell for i in range(tmp.shape[0])]
    sele_profiles.append(tmp)
sele_profiles = pd.concat(sele_profiles, ignore_index=True)
tmp = dict(zip(cpd_sele.index, [' '.join(nt) for nt in cpd_sele[['cmap_name', 'target']].to_numpy()]))
sele_profiles['Compound'] = [tmp[name] for name in sele_profiles['Compound']]

model = umap.UMAP(random_state=1, n_neighbors=200, metric='cosine')
Xe = model.fit_transform(sele_profiles.iloc[:, :-2])
Xe = pd.DataFrame(data=Xe, columns=['Feature 0', 'Feature 1'])
Xe['Compound'] = sele_profiles['Compound']
Xe['Cell line'] = sele_profiles['Cell line']

FIG = plt.figure(figsize=(7, 7))
xx = FIG.add_subplot()
sns.scatterplot(x='Feature 0', y='Feature 1', hue='Compound', style='Cell line', \
                s=120, data=Xe, ax=xx)
sns.move_legend(xx, "upper center", bbox_to_anchor=(0.5, -0.1), ncol=5, \
                frameon=False, fontsize='small')
plt.tight_layout()
#plt.savefig('../doc_fig_table/umap_cmpd_profile_robustness.pdf', transparent=True)










