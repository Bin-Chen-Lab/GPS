#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 15:34:06 2020

@author: jingxing
"""

import pandas as pd
import random
random.seed(1)
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
import seaborn as sns


inpdir = './data_for_code/predictable_genes/'
outdir = './data_for_code/'
CCLE = pd.read_hdf('/Users/jingxing/Documents/Data_Sources/CCLE/CCLE_expression_20Q1.hdf5')
lmg = set([l.strip() for l in open(inpdir + 'lincs_landmark_genes.txt')])

cell_genesele = dict([('MCF7', 'selected_genes_MCF7_t1.csv'), \
                      ('PC3', 'selected_genes_PC3_t1.csv'), \
                      ('VCAP', 'selected_genes_VCAP_t1.csv'), \
                      ('HEPG2', 'selected_genes_HEPG2_t0.csv')])
for k in cell_genesele:
    tmp = [l.strip() for l in open(inpdir + cell_genesele[k])][1:]
    cell_genesele[k] = set(tmp)


FIG, axes = plt.subplots(len(cell_genesele), 2, figsize=(12, len(cell_genesele) * 4))
for num, (cell, gset) in enumerate(cell_genesele.items()):
    other_union = set.union(*([v for k, v in cell_genesele.items() if k != cell]))
    spgset = set.difference(gset, other_union)
    
    S = CCLE.loc[cell, CCLE.columns.isin(spgset)]
    A = CCLE.loc[cell, CCLE.columns.isin(gset)]
    B = CCLE.loc[cell, CCLE.columns.isin(set.difference(lmg, gset))]
    s, p = stats.ranksums(A, B)
    p_with_rand = [p]
    for i in range(49):
        randset = set(random.sample(list(lmg), len(gset)))
        A1 = CCLE.loc[cell, CCLE.columns.isin(randset)]
        B1 = CCLE.loc[cell, CCLE.columns.isin(set.difference(lmg, randset))]
        _, p1 = stats.ranksums(A1, B1)
        p_with_rand.append(p1)
    fdr = multipletests(p_with_rand, method='fdr_bh')[1][0]
    sns.kdeplot(A, label='Predictable genes\nP = %.2E\nFDR = %.2E'%(p, fdr), ax=axes[num][0])
    sns.kdeplot(B, label='Other landmark genes', ax=axes[num][0])
    #sns.kdeplot(S, label='%s Only'%cell, ax=axes[num][0])
    axes[num][0].set_title(cell)
    axes[num][0].set_xlabel('TPM')
    axes[num][0].legend()
    
    rand_cell = random.sample([c for c in CCLE.index.to_list() if c != cell], 49)
    for i, rc in enumerate(rand_cell):
        X = CCLE.loc[rc, CCLE.columns.isin(spgset)]
        l = 'Random cell lines' if i == 0 else None
        sns.kdeplot(X, ax=axes[num][1], color='grey', linestyle='dotted', label=l)
    sns.kdeplot(S, label=cell, linewidth=3, ax=axes[num][1])
    axes[num][1].legend()
    axes[num][1].set_title('%s Only Predictable Genes'%cell)
    axes[num][1].set_xlabel('TPM')
plt.tight_layout()
FIG.savefig(outdir + 'predictable_gene_expr_distribution.pdf')






