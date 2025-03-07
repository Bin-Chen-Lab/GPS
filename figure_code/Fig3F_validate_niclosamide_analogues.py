#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  8 18:41:27 2022

@author: jingxing
"""

import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams['figure.dpi'] = 300
import seaborn as sns


wrkdir = '../data/HCC_efficacy_rges/'
df = pd.read_csv(wrkdir + 'niclosamideAnalogs_HCC_ZRGES_efficacy.csv', index_col=0)

FIG, axes = plt.subplots(2, 1, figsize=(4, 6))
X, Y, names = df['HepG2_pIC50'], df['Z_RGES'], df.index
axes[0].plot(X, Y, 'o')
for x, y, name in zip(X, Y, names):
    axes[0].text(x, y, name, color='grey')
cor0, p0 = spearmanr(X, Y)
axes[0].text(4.5, -6, 'HepG2\nCor = %.2f\nP = %.3f'%(cor0, p0), weight='bold', size=13)
#axes[0].set_xlabel('HepG2 pIC50')
axes[0].set_ylabel('Z-RGES')
axes[0].set_xlim(3.8, 7.5)

X, Y, names = df['Huh7_pIC50'], df['Z_RGES'], df.index
axes[1].plot(X, Y, 'o')
for x, y, name in zip(X, Y, names):
    axes[1].text(x, y, name, color='grey')
cor1, p1 = spearmanr(X, Y)
axes[1].text(4.5, -6, 'Huh7\nCor = %.2f\nP = %.3f'%(cor1, p1), weight='bold', size=13)
axes[1].set_xlabel('pIC50')
axes[1].set_ylabel('Z-RGES')
axes[1].set_xlim(3.8, 7.5)

plt.tight_layout()
#FIG.savefig('../doc_fig_table/HCC_rges_niclosamide_analogues.pdf', transparent=True)

s1 = df[['Z_RGES', 'HepG2_pIC50']]
s1.columns = ['Z_RGES', 'pIC50']
s1['Cell'] = ['HepG2' for i in range(s1.shape[0])]
s1['Cpd'] = s1.index
s2 = df[['Z_RGES', 'Huh7_pIC50']]
s2.columns = ['Z_RGES', 'pIC50']
s2['Cell'] = ['Huh7' for i in range(s2.shape[0])]
s2['Cpd'] = s2.index
df1 = pd.concat([s1, s2], ignore_index=True)

FIG = plt.figure(figsize=(3.5, 3))
sns.scatterplot(data=df1, x='pIC50', y='Z_RGES', style='Cell')
plt.text(4, -6, 'HepG2 Cor = %.2f\nHuh7 Cor = %.2f'%(cor0, cor1), weight='bold', size=12)
#plt.legend(loc='center', bbox_to_anchor=(0.5, 1.1), ncols=2, frameon=False)
plt.tight_layout()
FIG.savefig('../doc_fig_table/HCC_rges_niclosamide_analogues-2.png', transparent=True)




