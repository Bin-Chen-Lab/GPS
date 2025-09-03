#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 17:08:11 2021

@author: jingxing

Given a disease signature, permute RGES of any number of up/down regulated genes a compound induces.
"""

import argparse
import pandas as pd
import os, sys
import numpy as np
np.random.seed(10)
from multiprocessing import Pool


pathname = os.path.dirname(sys.argv[0])        
path = os.path.abspath(pathname)
GATE = path.split('code')[0]

parser = argparse.ArgumentParser()
parser.add_argument('--dzSigFile', type = str, default='../data/dzsig/DZSIG__TestJob1.csv')
parser.add_argument('--gene_input', type=str, help='Input list of gene symbols', default='preselect')
parser.add_argument('--cpu_num', type = int, help = 'Number of cpu cores to use', default=10)
parser.add_argument('--out_name_prefix', type=str, default='BGRD__')
args = parser.parse_args()



def score_rges_up(gene_up, dz_de):
    thrs = (2, 2)
    gene_up = set.intersection(gene_up, set(dz_de.index))
    up_tags_position = dz_de.loc[gene_up].sort_values(by='RANK').RANK
    if len(gene_up) >= thrs[0]:
        a_up = max(map(lambda x: float(x)/len(gene_up) - up_tags_position[x-1]/float(dz_de.shape[0]), \
                       range(1, len(gene_up)+1)))
        b_up = max(map(lambda x: up_tags_position[x-1]/float(dz_de.shape[0]) - float(x-1)/len(gene_up), \
                       range(1, len(gene_up)+1)))
        ks_up = a_up if a_up > b_up else -b_up
    else:
        ks_up = 0
    return ks_up

def score_rges_down(gene_down, dz_de):
    thrs = (2, 2)
    gene_down = set.intersection(gene_down, set(dz_de.index))
    down_tags_position = dz_de.loc[gene_down].sort_values(by='RANK').RANK
    if len(gene_down) >= thrs[1]:
        a_down = max(map(lambda x: float(x)/len(gene_down) - down_tags_position[x-1]/float(dz_de.shape[0]), \
                       range(1, len(gene_down)+1)))
        b_down = max(map(lambda x: down_tags_position[x-1]/float(dz_de.shape[0]) - float(x-1)/len(gene_down), \
                       range(1, len(gene_down)+1)))
        ks_down = a_down if a_down > b_down else -b_down
    else:
        ks_down = 0
    return ks_down


def single_run(params):
    num_genes, dzde, NP, rdc_transcrptm = params
    tmp_up, tmp_down = [], []
    for i in range(NP):
        gene_sele = set(np.random.choice(rdc_transcrptm, num_genes, replace=False))
        score_up = score_rges_up(gene_sele, dzde)
        score_down = score_rges_down(gene_sele, dzde)
        tmp_up.append(score_up)
        tmp_down.append(score_down)
    return tmp_up, tmp_down, num_genes



job_prefix = os.path.basename(args.dzSigFile).replace('DZSIG__', '').replace('.csv', '')
wrkdir = GATE + 'data/dzsig/'
B, NP = 1, 1500
fn = GATE + 'data/input_gene_features/selected_genes_2198.csv' if args.gene_input == 'preselect' else args.gene_input
gene_sele = [l.strip() for l in open(fn)]
dzde = pd.read_csv(args.dzSigFile, index_col='GeneSymbol')

dzde['RANK'] = dzde['Value'].rank(ascending=False)
rdc_transcrptm = [g for g in gene_sele if g in dzde.index]
print(len(rdc_transcrptm))

bins = range(2, len(rdc_transcrptm), B)

def param_gen(bins):
    for num in bins:
        yield num, dzde, NP, rdc_transcrptm


pl = Pool(args.cpu_num)
mtx, title = [], []
for tmp_up, tmp_down, num in pl.imap(single_run, param_gen(bins)):
    mtx.append(tmp_up)
    mtx.append(tmp_down)
    title += ['Up_%s'%num, 'Down_%s'%num]
    if num / B % 5 == 0: print(num)
pl.close()    


mtx = pd.DataFrame(data=mtx, index=title, columns=['RGES_Rnd%s'%i for i in range(NP)])
mtx.to_pickle(wrkdir + args.out_name_prefix + job_prefix + '.pkl')










