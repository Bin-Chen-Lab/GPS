#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 19:32:26 2020

@author: jingxing
"""

import argparse
import pandas as pd
import numpy as np
import os, sys
from multiprocessing import Pool

pathname = os.path.dirname(sys.argv[0])        
path = os.path.abspath(pathname)
GATE = path.split('code')[0]

parser = argparse.ArgumentParser()
parser.add_argument('--result_dir', type = str, help = 'dir to save result txt files', default = GATE + 'data/profile_pred/')
parser.add_argument('--cmpd_input', type=str, help='Input csv for compound ID and SMILES', default='../data/input_cmpd_gene/cmpd__TestJob0.csv')
parser.add_argument('--cpu_num', type = int, help = 'Number of cpu cores to use', default=10)
args = parser.parse_args()
prefix = os.path.basename(args.cmpd_input).replace('.csv', '').replace('cmpd__', '')
cpunum = args.cpu_num


def catg_assign(abc):
    if np.max(abc) < 0.95: return 0
    return np.argmax(abc) - 1

def single_run0(files2merge):
    out_name = files2merge[0].split('/')[-1]
    tmp = [np.loadtxt(fn) for fn in files2merge]
    tmp = np.median(np.stack(tmp), axis=0) if len(tmp) > 1 else tmp[0]
    return tmp, out_name

def single_run1(filename):
    #print(filename)
    probs = np.loadtxt(filename)
    preds = np.apply_along_axis(catg_assign, 1, probs)
    return preds


if __name__ == '__main__':
    ### Merge probabilities
    wrkdir = args.result_dir if args.result_dir.endswith('/') else args.result_dir + '/'
    ds = 'prob_' + prefix
    outdir = wrkdir + 'MEDIAN/%s/'%ds
    if not os.path.exists(outdir): os.mkdir(outdir)

    prob_dirs = ['HEPG2_t0/', 'MCF7_t1/', 'PC3_t1/', 'VCAP_t1/']
    geneset = set()
    for d in prob_dirs:
        for fn in os.listdir(wrkdir + d + ds):
            if fn.startswith('.'): continue
            geneset.add(fn.split('_')[0])
    files2merge_list = [['%s%s%s/%s_probs.txt'%(wrkdir, c, ds, gene) for c in prob_dirs \
                         if os.path.exists('%s%s%s/%s_probs.txt'%(wrkdir, c, ds, gene))] \
                        for gene in geneset]

    p = Pool(cpunum)
    for mrgd_mtx, out_name in p.imap(single_run0, files2merge_list):
        np.savetxt(outdir + out_name, mrgd_mtx)
    p.close()
    
    ### Probs to preds
    wrkdir = args.result_dir + 'MEDIAN/'
    # outdir = wrkdir + 'preds_%s/'%prefix
    outdir = wrkdir + 'preds_all/'
    inpdir = wrkdir + 'prob_%s/'%prefix
    if not os.path.exists(outdir): os.mkdir(outdir)
    
    drug_info = pd.read_csv(args.cmpd_input, index_col='ID')
    filename_list = [inpdir + fn for fn in os.listdir(inpdir) if fn.endswith('_probs.txt')]
    gene_info = [fn.replace(inpdir, '').replace('_probs.txt', '') for fn in filename_list]
    
    result = np.empty(shape=(0, drug_info.shape[0]))
    for fn in filename_list:
        preds = single_run1(fn)
        result = np.row_stack((result, preds))

    '''
    p = Pool(cpunum)
    result = np.empty(shape=(0, drug_info.shape[0]))
    for preds in p.imap(single_run1, filename_list):
        result = np.row_stack((result, preds))
    p.close()    
    '''
    result = pd.DataFrame(data=result, index=gene_info, columns=drug_info.index)
    result.to_csv(outdir + '%s_MEDIAN_GeneExpressionChange.csv'%prefix)


    print('WORK DONE: MERGE')
