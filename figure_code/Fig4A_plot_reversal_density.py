#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 16:22:25 2022

@author: jingxing
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["font.size"] = 12
import seaborn as sns

wrkdir = '../data/HCC_efficacy_rges/'
dzsig_file = wrkdir + 'DZSIG__HCC_opt.csv'
dzsig_title = 'HCC Signature'
profiles_outfile = wrkdir + 'HCC_opt_candidate_GE_profiles.csv'
reversal_outfile = '../doc_fig_table/HCC_5candidates_reversal_density.pdf'
M = {'ZINC000067803790': '44443110\n(-9.83)', \
     'ZINC000004631199': '4896-3969\n(-9.65)', \
     'ZINC000000086363': 'PB56874852\n(-9.51)', \
     'ZINC000100424631': '1C-001_4622\n(-9.36)', \
     'ZINC000169294501': 'STOCK7S-74424\n(-9.10)'}

druglib = pd.read_csv(profiles_outfile, index_col=0)
if M != None:
    druglib = druglib.loc[:, M.keys()]
    druglib.columns = M.values()

dzsig = pd.read_csv(dzsig_file, index_col=0)
dzsig = dzsig.sort_values(by='Value')
dzsig['RANK'] = dzsig['Value'].rank(pct=True)
druglib = druglib.reindex(dzsig.index)
plot_df = []
for cpd in druglib.columns:
    profile = druglib[cpd]
    up_genes = profile.loc[profile == 1].index
    down_genes = profile.loc[profile == -1].index
    tmp = [[dzsig.loc[g, 'RANK'], cpd, 'Up'] for g in up_genes]
    plot_df += tmp
    tmp = [[dzsig.loc[g, 'RANK'], cpd, 'Down'] for g in down_genes]
    plot_df += tmp
plot_df = pd.DataFrame(data=plot_df, columns=['Value', 'Compound', 'Change'])

my_pal = {'Up': 'tomato', 'Down': 'cornflowerblue'}
FIG, axes = plt.subplots(2, 1, figsize=(5, 4), gridspec_kw={'height_ratios': [0.75, druglib.shape[1]]})
sns.heatmap(dzsig[['Value']].T, center=0, vmax=3, vmin=-3, cmap='coolwarm', cbar=False, yticklabels=False, ax=axes[0])
axes[0].set_yticks([0.5])
axes[0].set_yticklabels(['log2 Fold Change'])
axes[0].set_xticks([])
axes[0].set_xlabel(None)
axes[0].set_title(dzsig_title)
sns.violinplot(y='Compound', x='Value', hue='Change', split=True, orient='h', kde_kws={'kernel': 'triw'}, \
               data=plot_df, inner='stick', palette=my_pal, ax=axes[1])
axes[1].set_xlabel(None)
axes[1].get_xaxis().set_ticks([])
axes[1].set_ylabel(None)
axes[1].set_xlim(0, 1)
axes[1].legend(bbox_to_anchor=(1, 0), ncol=2)
plt.tight_layout()
plt.show()
FIG.savefig(reversal_outfile, transparent=True)


