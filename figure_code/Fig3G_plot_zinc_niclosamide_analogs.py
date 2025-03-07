#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 16:15:10 2022

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
df = pd.read_csv(wrkdir + 'ZINC_niclosamide_analogs_sim0.3_Z2.csv', index_col=0)
df['S_uM'] = 10 ** df['LogS'] * 10 ** 6

mol_rm = ['ZINC000003874496']
mol_sele = ['ZINC000003874496', 'ZINC000004595567', 'ZINC000004017999', 'ZINC000000246015', 'ZINC000004978017']

#df = df.loc[~df.index.isin(mol_rm)]
FIG = plt.figure(figsize=(5, 4))
#cmp = matplotlib.cm.get_cmap('Blues').reversed()
plt.scatter(df['Z_RGES'], df['Similarity'], s=10, c=df['S_uM'], cmap='Blues', vmin=1, vmax=50)
plt.colorbar()
plt.scatter(df.loc[mol_sele, 'Z_RGES'], df.loc[mol_sele, 'Similarity'], s=100, c='red', marker='$O$')
for k, row in df.loc[mol_sele].iterrows():
    X, Y, S = row['Z_RGES'], row['Similarity'], row['S_uM']
    plt.text(X+0.05, Y+0.01, '%.1f'%S)
plt.xlabel('Z-RGES')
plt.ylabel('Similarity to niclosamide')
plt.tight_layout()
FIG.savefig('../doc_fig_table/HCC_ZINC_niclosamide_analogs.pdf', transparent=True)

