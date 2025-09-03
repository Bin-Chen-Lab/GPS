#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 21:54:55 2021

@author: jingxing
"""

import argparse
import os, sys
import pandas as pd
from scipy.stats import zscore
import numpy as np
np.random.seed(1)
from multiprocessing import Pool


pathname = os.path.dirname(sys.argv[0])        
path = os.path.abspath(pathname)
GATE = path.split('code')[0]

parser = argparse.ArgumentParser()
parser.add_argument('--dzSigFile', type = str, default='../data/dzsig/DZSIG__TestJob1.csv')
parser.add_argument('--cmpdLibFile', type = str, default='../data/profile_pred/MEDIAN/preds_TestJob0/TestJob0_MEDIAN_GeneExpressionChange.csv')
parser.add_argument('--cpu_num', type = int, help = 'Number of cpu cores to use', default=10)
parser.add_argument('--RGES_bgrd_file', type=str, default='../data/dzsig/BGRD__TestJob1.pkl')
args = parser.parse_args()
job_prefix = os.path.basename(args.dzSigFile).replace('DZSIG__', '').replace('.csv', '')


def score_norm_rges(gene_up, gene_down, dz_de, bgrd):
    up_tags_position = dz_de.loc[gene_up].sort_values(by='RANK').RANK
    down_tags_position = dz_de.loc[gene_down].sort_values(by='RANK').RANK
    if len(gene_up) >= 2:
        a_up = max(map(lambda x: float(x)/len(gene_up) - up_tags_position[x-1]/float(dz_de.shape[0]), \
                       range(1, len(gene_up)+1)))
        b_up = max(map(lambda x: up_tags_position[x-1]/float(dz_de.shape[0]) - float(x-1)/len(gene_up), \
                       range(1, len(gene_up)+1)))
        ks_up = a_up if a_up > b_up else -b_up
    else:
        ks_up = 0
    if len(gene_down) >= 2:
        a_down = max(map(lambda x: float(x)/len(gene_down) - down_tags_position[x-1]/float(dz_de.shape[0]), \
                       range(1, len(gene_down)+1)))
        b_down = max(map(lambda x: down_tags_position[x-1]/float(dz_de.shape[0]) - float(x-1)/len(gene_down), \
                       range(1, len(gene_down)+1)))
        ks_down = a_down if a_down > b_down else -b_down
    else:
        ks_down = 0
    

    # z-score normalization
    bgrd_up, bgrd_down = bgrd
    z_up = zscore(np.array([ks_up] + bgrd_up))[0]
    z_down = zscore(np.array([ks_down] + bgrd_down))[0]
    score = z_up - z_down
    
    return score


def profile_score(params):
    cpd, profile, dzde, permt_bgrd, score_fun = params
    up_genes = set(profile.loc[profile == 1].index)
    down_genes = set(profile.loc[profile == -1].index)
    
    if not (permt_bgrd is None):
        bins_up = [int(b.split('_')[1]) for b in permt_bgrd.index if b.startswith('Up')]
        bins_down = [int(b.split('_')[1]) for b in permt_bgrd.index if b.startswith('Down')]
        mark_up = min(bins_up, key=lambda x: abs(x - len(up_genes)))
        mark_down = min(bins_down, key=lambda x: abs(x - len(down_genes)))
        bgrd_up = permt_bgrd.loc['Up_%s'%mark_up].to_list()
        bgrd_down = permt_bgrd.loc['Down_%s'%mark_down].to_list()
        rges_bgrd = (bgrd_up, bgrd_down)
        rev_score = score_fun(up_genes, down_genes, dzde, rges_bgrd)
    else:
        rev_score = score_fun(up_genes, down_genes, dzde)
    
    return cpd, rev_score


### Load compound profiles, score backgrounds and disease signatures
if args.cmpdLibFile == 'ZINC':
    tmp = np.load(GATE + 'data/profile_pred/MEDIAN/ZINC_strong.npz')
    colnames = tmp['columns'].tolist()
    rownames = tmp['index'].tolist()
    library = pd.DataFrame(data=tmp['X'], index=rownames, columns=colnames).T
    del tmp, rownames, colnames
elif args.cmpdLibFile == 'HTS':
    tmp = np.load(GATE + 'data/profile_pred/MEDIAN/ENAMINE_HTS_strong.npz')
    colnames = tmp['columns'].tolist()
    rownames = tmp['index'].tolist()
    library = pd.DataFrame(data=tmp['X'], index=rownames, columns=colnames).T
    del tmp, rownames, colnames
elif args.cmpdLibFile == 'CNS':
    tmp = np.load(GATE + 'data/profile_pred/MEDIAN/ENAMINE_CNS_strong.npz')
    colnames = tmp['columns'].tolist()
    rownames = tmp['index'].tolist()
    library = pd.DataFrame(data=tmp['X'], index=rownames, columns=colnames).T
    del tmp, rownames, colnames
elif args.cmpdLibFile == 'DIPG':
    tmp = np.load(GATE + 'data/profile_pred/MEDIAN/DIPG_ALL.npz')
    colnames = tmp['columns'].tolist()
    rownames = tmp['index'].tolist()
    library = pd.DataFrame(data=tmp['X'], index=rownames, columns=colnames).T
    del tmp, rownames, colnames
elif args.cmpdLibFile == 'LINCS':
    tmp = np.load(GATE + 'data/profile_pred/MEDIAN/LINCS.npz')
    colnames = tmp['columns'].tolist()
    rownames = tmp['index'].tolist()
    library = pd.DataFrame(data=tmp['X'], index=rownames, columns=colnames).T
    del tmp, rownames, colnames
else:
    library = pd.read_csv(args.cmpdLibFile, index_col=0).T
    library.index = library.index.astype(str)

dzde = pd.read_csv(args.dzSigFile, index_col='GeneSymbol')
dzde['RANK'] = dzde['Value'].rank(ascending=False)

permt_bgrd = pd.read_pickle(args.RGES_bgrd_file)


### Align compound profile genes to dzsig
rdc_transcrptm = set.intersection(set(library.columns), set(dzde.index))
library = library.loc[:, library.columns.isin(rdc_transcrptm)]


### Prepare iterator and run
def param_iter():
    for cpd, profile in library.iterrows():
        yield cpd, profile, dzde, permt_bgrd, score_norm_rges

pl = Pool(args.cpu_num)

result = []
for cpd, rev_score in pl.imap(profile_score, param_iter()):
    result.append([cpd, rev_score])
    if len(result) % 500 == 0: print(len(result) / float(library.shape[0]), 'finished')
pl.close()

result = pd.DataFrame(data=result, columns=['ID', 'Z_RGES'])
result = result.sort_values(by='Z_RGES')
result.to_csv(GATE + 'data/reversal_score/%s_RGES_norm.csv'%job_prefix, index=False)




