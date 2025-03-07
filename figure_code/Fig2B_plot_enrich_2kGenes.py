#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:06:43 2021

@author: jingxing
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["font.family"] = "Arial"
import seaborn as sns


wrkdir = '/Users/jingxing/Documents/myc_reverse/DL/manuscript_cpdProfilePred/doc_fig_table/'
df = pd.read_excel(wrkdir + 'summary_Gorilla_GO_enrichment_2kGenes.xlsx', index_col=0)
df = pd.concat([df.loc[df.Category == t].iloc[:5] for t in ['BP', 'CC', 'MF']])
df['x'] = range(df.shape[0], 0, -1)
df['Description'] = [s[0].upper() + s[1:] for s in df.Description]

FIG = plt.figure(figsize=(6, 6))
sc = sns.scatterplot(x='logQ', y='x', size='Enrichment', data=df, sizes=(80, 250))
plt.setp(sc.get_legend().get_title(), fontsize='12')
plt.yticks(ticks=df.x, labels=df.Description, size=14)
plt.ylabel(None)
plt.xlabel('-log10 q value', size=14)
plt.xlim(0, 30)
plt.tight_layout()
FIG.savefig(wrkdir + 'summary_GOenrich_2kGenes.pdf', transparent=True)



