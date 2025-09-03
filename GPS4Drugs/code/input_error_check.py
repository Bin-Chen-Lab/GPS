#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 30 19:45:08 2021

@author: jingxing
"""

import os, sys
import pandas as pd
from rdkit import Chem

pathname = os.path.dirname(sys.argv[0])        
path = os.path.abspath(pathname)
GATE = path.split('code')[0]

def CheckError(option, value):
    warnmsg = ''
    
    if option == 'cmpd_input':
        if not os.path.exists(value): raise Exception('ERROR: Compound input file not found!')
        if os.stat(value).st_size > 10 ** 7: raise Exception('ERROR: Compound input file too large (over 10M)!')
        
        table = pd.read_csv(value)
        if not 'ID' in table.columns: raise Exception('ERROR: Please specify "ID"!')
        if not 'SMILES' in table.columns: raise Exception('ERROR: Please specify "SMILES"!')
        
        shape_org, changed = table.shape[0], False
        if shape_org > 1000: raise Exception('ERROR: reached compound number limitation: 1000!')
        
        table = table[(~table.ID.isna()) & (~table.SMILES.isna())]
        if table.shape[0] != shape_org:
            shape_org, changed = table.shape[0], True
            tmp = 'WARNING: Found NA in ID or SMILES and removed.\n'
            warnmsg += tmp
            print(tmp)
        
        table = table[~table.ID.duplicated(keep='first')]
        if table.shape[0] != shape_org:
            shape_org, changed = table.shape[0], True
            tmp = 'WARNING: Found duplicated IDs and kept the first.\n'
            warnmsg += tmp
            print(tmp)
        
        vid = []
        for i, s in zip(table.index, table.SMILES):
            mol = Chem.MolFromSmiles(s)
            if mol: vid.append(i)
        table = table.loc[vid]
        if table.shape[0] != shape_org:
            changed = True
            tmp = 'WARNING: Found invalid SMILES and removed.\n'
            warnmsg += tmp
            print(tmp)
        
        if changed:
            os.rename(value, value + '_ERR.csv')
            table.to_csv(value, index=False)
    
    
    elif option == 'gene_input':
        if value != 'preselect':
            if not os.path.exists(value): raise Exception('ERROR: Gene input file not found!')
            if os.stat(value).st_size > 10 ** 7: raise Exception('ERROR: Gene input file too large (over 10M)!')
            
            gl = set([l.strip() for l in open(value)])
            df = pd.read_csv(GATE + 'data/input_gene_features/go_fingerprints_allGenesExt.csv', index_col=0)
            fail_genes = '|'.join(set.difference(gl, set(df.index)))
            df = df.loc[df.index.isin(gl)]
            task_num = df.shape[0]
            if task_num == 0: raise Exception('ERROR: Not able to generate features for any of the input genes!')
            if task_num != len(gl):
                tmp = 'WARNING: Removed unpredictable gene(s): %s\n'%fail_genes
                warnmsg += tmp
                print(tmp)
                os.rename(value, value + '_ERR.csv')
                with open(value, 'w') as f:
                    f.write('\n'.join(df.index) + '\n')
            del df
            
            
    elif option == 'cpu_num':
        if any([type(value) is not int, value < 1, value > 14]):
            raise Exception('ERROR: Please specify a valid number of cores (1 - 14), or leave it as default!')
    
    
    elif option == 'dzSigFile':
        if not os.path.exists(value): raise Exception('ERROR: Disease signature input file not found!')
        if os.stat(value).st_size > 10 ** 8: raise Exception('ERROR: Disease signature input file too large (over 100M)!')
        
        table = pd.read_csv(value)
        if not 'GeneSymbol' in table.columns: raise Exception('ERROR: Please specify "GeneSymbol"!')
        if not 'Value' in table.columns: raise Exception('ERROR: Please specify "Value"! Positive values indicate up-regulation, and negative indicate down-regulation.')
        
        shape_org, changed = table.shape[0], False
        table = table[(~table.GeneSymbol.isna()) & (~table.Value.isna())]
        if table.shape[0] != shape_org:
            shape_org, changed = table.shape[0], True
            tmp = 'WARNING: Found NA in GeneSymbol or Value and removed.\n'
            warnmsg += tmp
            print(tmp)
        
        table1 = table.groupby(by='GeneSymbol').mean()
        if table1.shape[0] != shape_org:
            shape_org, changed = table1.shape[0], True
            tmp = 'WARNING: Found duplicated gene symbols and merged as average.'
            warnmsg += tmp
            print(tmp)
            table = table1.copy()
            table['GeneSymbol'] = table.index
            del table1
        
        # Keep at most 2000 genes each side
        table = table.sort_values(by='Value')
        U, D = table[table.Value > 0], table[table.Value < 0]
        if U.shape[0] > 2000 or D.shape[0] > 2000:
            changed = True
            table = pd.concat([D.iloc[:2000], U.iloc[-2000:]])
            tmp = 'WARNING: Disease signature too long. Trimmed to 2000 each side.\n'
            warnmsg += tmp
            print(tmp)
        
        if changed:
            os.rename(value, value + '_ERR.csv')
            table.to_csv(value, index=False)
    
    
    elif option == 'cmpdLibID':
        if value != 'ZINC':
            file_dir = GATE + 'data/profile_pred/MEDIAN/preds_%s/%s_MEDIAN_GeneExpressionChange.csv'%(value, value)
            if not os.path.exists(file_dir):
                raise Exception('ERROR: Gene expression profile with ID "%s" input file not found!'%value)
    
    
    elif option == 'RGES_bgrd_ID':
        if value != 'NONE':
            file_dir = GATE + 'data/dzsig/BGRD__%s.pkl'%value
            if not os.path.exists(file_dir):
                raise Exception('ERROR: RGES background file with ID "%s" was not found!'%value)
    
    
    else:
        raise Exception('Unknown option!')


    return warnmsg



        
            
