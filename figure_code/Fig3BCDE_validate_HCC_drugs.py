#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 21:07:49 2020

@author: jingxing
"""

import pandas as pd
from scipy.stats import ranksums, spearmanr, zscore, rankdata, beta
import numpy as np
np.random.seed(1)
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["font.family"] = "Arial"
import seaborn as sns
from multiprocessing import Pool

    

def score_norm_rges(gene_up, gene_down, dz_de, bgrd):
    bgrd_up, bgrd_down = bgrd
    up_tags_position = dz_de.loc[gene_up].sort_values(by='RANK').RANK
    down_tags_position = dz_de.loc[gene_down].sort_values(by='RANK').RANK
    
    if len(gene_up) >= 2:
        a_up = max(map(lambda x: x/len(gene_up) - up_tags_position[x-1]/dz_de.shape[0], \
                       range(1, len(gene_up)+1)))
        b_up = max(map(lambda x: up_tags_position[x-1]/dz_de.shape[0] - (x-1)/len(gene_up), \
                       range(1, len(gene_up)+1)))
        ks_up = a_up if a_up > b_up else -b_up
    else:
        ks_up = 0
    if len(gene_down) >= 2:
        a_down = max(map(lambda x: x/len(gene_down) - down_tags_position[x-1]/dz_de.shape[0], \
                       range(1, len(gene_down)+1)))
        b_down = max(map(lambda x: down_tags_position[x-1]/dz_de.shape[0] - (x-1)/len(gene_down), \
                       range(1, len(gene_down)+1)))
        ks_down = a_down if a_down > b_down else -b_down
    else:
        ks_down = 0
    z_up = zscore([ks_up] + bgrd_up)[0]
    z_down = zscore([ks_down] + bgrd_down)[0]
    score = z_up - z_down
    return score, z_up, z_down

def score_rges(gene_up, gene_down, dz_de):
    up_tags_position = dz_de.loc[gene_up].sort_values(by='RANK').RANK
    down_tags_position = dz_de.loc[gene_down].sort_values(by='RANK').RANK
    
    if len(gene_up) >= 2:
        a_up = max(map(lambda x: x/len(gene_up) - up_tags_position[x-1]/dz_de.shape[0], \
                       range(1, len(gene_up)+1)))
        b_up = max(map(lambda x: up_tags_position[x-1]/dz_de.shape[0] - (x-1)/len(gene_up), \
                       range(1, len(gene_up)+1)))
        ks_up = a_up if a_up > b_up else -b_up
    else:
        ks_up = 0
    if len(gene_down) >= 2:
        a_down = max(map(lambda x: x/len(gene_down) - down_tags_position[x-1]/dz_de.shape[0], \
                       range(1, len(gene_down)+1)))
        b_down = max(map(lambda x: down_tags_position[x-1]/dz_de.shape[0] - (x-1)/len(gene_down), \
                       range(1, len(gene_down)+1)))
        ks_down = a_down if a_down > b_down else -b_down
    else:
        ks_down = 0
    return ks_up - ks_down, ks_up, ks_down



def penalty_pct(profile, dz_de_ori):
    dz_de = dz_de_ori[dz_de_ori.index.isin(profile.index)]
    dz_gene_up = set(dz_de[dz_de['Value'] > 0].index)
    dz_gene_down = set(dz_de[dz_de['Value'] < 0].index)
    gene_up = set(profile.loc[profile == 1].index)
    gene_down = set(profile.loc[profile == -1].index)
    p_up = len(set.intersection(gene_up, dz_gene_up)) / len(dz_gene_up)
    p_down = len(set.intersection(gene_down, dz_gene_down)) / len(dz_gene_down)
    return (p_up + p_down) / 2, len(gene_up), len(gene_down)


def profile_score(params):
    cpd, profile, dzde, transcrptm, permt_rges = params
    up_genes = sorted(set(profile.loc[profile == 1].index))
    down_genes = sorted(set(profile.loc[profile == -1].index))
    
    bins_up = [int(b.split('_')[-1]) for b in permt_rges.index if b.startswith('Up_')]
    bins_down = [int(b.split('_')[-1]) for b in permt_rges.index if b.startswith('Down_')]
    mark_up = min(bins_up, key=lambda x: abs(x - len(up_genes)))
    mark_down = min(bins_down, key=lambda x: abs(x - len(down_genes)))
    bgrd_up = permt_rges.loc['Up_%s'%mark_up].to_list()
    bgrd_down = permt_rges.loc['Down_%s'%mark_down].to_list()
    rges_bgrd = (bgrd_up, bgrd_down)
    zrges, z_up, z_down = score_norm_rges(up_genes, down_genes, dzde, rges_bgrd)
    penalty, num_up, num_down = penalty_pct(profile, dzde)
    rges, ks_up, ks_down = score_rges(up_genes, down_genes, dzde)
    return cpd, zrges, z_up, z_down, rges, penalty, num_up, num_down


        
### Working directory and input files
### Use high-throughput dataset to evaluate enrichment, manually curated dataset to evalutate correlation.
wrkdir = './data_for_code/HCC_efficacy_rges/'

com_roc = pd.read_hdf(wrkdir + 'HCC_CTRP_MEDIAN_mlp_2kGenes_preds_0.95.hdf5').T
com_cor = pd.read_hdf(wrkdir + 'HCC_drugs_MEDIAN_mlp_2kGenes_preds_0.95.hdf5').T

dzde = pd.read_csv(wrkdir + 'DZSIG__HCC_opt.csv', index_col='GeneSymbol')
dzde['RANK'] = dzde['Value'].rank(ascending=False)
permt_bgrd = pd.read_csv(wrkdir + 'BGRD__HCC_opt.csv', index_col=0)

meta_cor = pd.read_csv(wrkdir + 'HCC_drug_efficacy.csv', index_col=0)
meta_roc = pd.read_csv(wrkdir + 'CTRP_HEPG2_efficacy.csv', index_col='cpd_name')


### Define the transcriptome
rdc_transcrptm = set.intersection(set(com_cor.columns), set(com_roc.columns), set(dzde.index))

### Align compound profiles, activity information and disease signatures
com_cor = com_cor.loc[:, com_cor.columns.isin(rdc_transcrptm)]
com_roc = com_roc.loc[:, com_roc.columns.isin(rdc_transcrptm)]
com_roc = com_roc.loc[com_roc.index.isin(meta_roc.index)]
meta_cor = meta_cor.loc[com_cor.index]
meta_roc = meta_roc.loc[com_roc.index]
meta_cor['pIC50'] = [-np.log10(i * 10 ** (-9)) for i in meta_cor['standard_value (nm)']]
meta_cor['Active'] = [True if i > 5 else False for i in meta_cor['pIC50']]

meta_roc['Active'] = [True if auc <= 9 else False for auc in meta_roc.area_under_curve]
meta_roc = meta_roc.loc[(meta_roc.area_under_curve <= 9) | (meta_roc.area_under_curve >= 12)]
meta_roc.loc['niclosamide', 'Active'] = True
com_roc = com_roc.loc[com_roc.index.isin(meta_roc.index)]


### Calculate reversal scores of the correlation dataset
eva_cor = []
for cpd, profile in com_cor.iterrows():
    params = (cpd, profile, dzde, rdc_transcrptm, permt_bgrd)
    res = profile_score(params)
    eva_cor.append(list(res))
eva_cor = pd.DataFrame(data=eva_cor, columns=['DrugName', 'Z_RGES', 'Z_Up', 'Z_Down', 'Raw_RGES', 'Penalty', 'Num_Up', 'Num_Down'])
eva_cor.index = eva_cor.iloc[:, 0]
eva_cor = eva_cor.iloc[:, 1:]
eva_cor = pd.concat([meta_cor, eva_cor], axis=1)


### Calculate reversal scores of the enrichment dataset
eva_roc = []
for cpd, profile in com_roc.iterrows():
    params = (cpd, profile, dzde, rdc_transcrptm, permt_bgrd)
    res = profile_score(params)
    eva_roc.append(list(res))
eva_roc = pd.DataFrame(data=eva_roc, columns=['DrugName', 'Z_RGES', 'Z_Up', 'Z_Down', 'Raw_RGES', 'Penalty', 'Num_Up', 'Num_Down'])
eva_roc.index = eva_roc.iloc[:, 0]
eva_roc = eva_roc.iloc[:, 1:]
eva_roc = pd.concat([meta_roc, eva_roc], axis=1)
eva_roc['Active'] = eva_roc['Active'].fillna(False)
#eva_roc = eva_roc[(eva_roc['Num_Up'] >= 5) & (eva_roc['Num_Down'] >= 5)]

FS = (3.5, 3.5)
### Evaluate correlation: pIC50 ~ reversal scores
eva_cor_tstd = eva_cor.loc[~eva_cor.pIC50.isna()]
Xs, Ys = eva_cor_tstd['pIC50'], eva_cor_tstd['Z_RGES']
FIG = plt.figure(figsize=FS)
plt.plot(Xs, Ys, 'o', color='black')
plt.xlabel('pIC50')
plt.ylabel('Z-RGES')
r, p = spearmanr(Xs, Ys)
plt.text(7.5, -2.5, 'Cor = %.3f\nP = %.4f'%(r, p), size=12, weight='bold')
for drug, row in eva_cor_tstd.iterrows():
    x, y = row['pIC50'], row['Z_RGES']
    plt.text(x, y, drug, color='grey')
plt.tight_layout()
#FIG.savefig('./data_for_code/HCCinSilico_zrges_pIC50.pdf', transparent=True)

Ys = eva_cor_tstd['Raw_RGES']
FIG = plt.figure(figsize=FS)
plt.plot(Xs, Ys, 'o', color='cornflowerblue')
plt.xlabel('pIC50')
plt.ylabel('Raw RGES')
r, p = spearmanr(Xs, Ys)
plt.text(7.5, -1.1, 'Cor = %.3f\nP = %.4f'%(r, p), size=12, weight='bold')
for drug, row in eva_cor_tstd.iterrows():
    x, y = row['pIC50'], row['Raw_RGES']
    plt.text(x, y, drug, color='grey')
plt.tight_layout()
#FIG.savefig('./data_for_code/HCCinSilico_rges_pIC50.pdf', transparent=True)


### Evaluate AUC-ROC: area_under_curve ~ reversal scores
y_true = eva_roc.Active.to_list()
ratio = sum(y_true) / len(y_true)
score_zrges = -eva_roc['Z_RGES']
auc_zrges = roc_auc_score(y_true, score_zrges)
fprz, tprz, thrz = roc_curve(y_true, score_zrges)
score_rges = -eva_roc['Raw_RGES']
auc_rges = roc_auc_score(y_true, score_rges)
fpr, tpr, thr = roc_curve(y_true, score_rges)

FIG = plt.figure(figsize=FS)    
plt.plot(fprz, tprz, color='black', label='Z-RGES AU-ROC = %.3f'%auc_zrges)
plt.plot(fpr, tpr, color='cornflowerblue', label='Raw RGES AU-ROC = %.3f'%auc_rges)
plt.text(0, 0.9, 'CTRP HepG2\n%s compounds'%eva_roc.shape[0], size=12, weight='bold')
plt.legend()
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.tight_layout()
#FIG.savefig('./data_for_code/HCCinSilico_ROC_zrges_rges.pdf', transparent=True)

topn = 100        
ranking_zrges = sorted([(s, yt) for s, yt in zip(score_zrges, y_true)], reverse=True)
prcs_zrges = [sum([r[1] for r in ranking_zrges[:i]]) / i for i in range(1, topn+1)]
scr_zrges = range(1, topn+1)
ranking_rges = sorted([(s, yt) for s, yt in zip(score_rges, y_true)], reverse=True)
prcs_rges = [sum([r[1] for r in ranking_rges[:i]]) / i for i in range(1, topn+1)]
scr_rges = range(1, topn+1)

FIG = plt.figure(figsize=FS)
plt.plot(scr_zrges, prcs_zrges, color='black', label='Z-RGES')
plt.plot(scr_rges, prcs_rges, color='cornflowerblue', label='Raw RGES')
plt.plot([0.0, topn], [ratio, ratio], '--', color='grey', label='Overal positive ratio')
plt.legend()
plt.xlabel('Rank')
plt.ylabel('Hit Rate')
plt.ylim(0.0, 1.05)
plt.text(20, 0.5, 'CTRP HepG2\n%s compounds'%eva_roc.shape[0], size=12, weight='bold')
plt.tight_layout()
#FIG.savefig('./data_for_code/HCCinSilico_HitRate_zrges_rges.pdf', transparent=True)
    
print('Number of actives: %s\nActive ratio: %s'%(sum(y_true), ratio))


### Compare hit rate of octad compounds and unseen compounds in PRISM
eff = pd.read_csv(wrkdir + 'PRISM_HUH7_efficacy.csv', index_col=0)
eff = eff.loc[~eff['Active'].isna()]
srges = pd.read_csv(wrkdir + 'sRGES_HCC_octad.csv', index_col='pert_iname')
zrges = pd.read_csv(wrkdir + 'PRISM_HCC_opt_RGES_norm.csv', index_col=0)
eff['sRGES'] = srges.reindex(eff.index)['sRGES']
eff['Z_RGES'] = zrges.reindex(eff.index)['Z_RGES']
eff_sub_s = eff.loc[~eff['sRGES'].isna()]
eff_sub_z = eff.loc[(~eff.index.isin(eff_sub_s.index)) & (~eff['Z_RGES'].isna())]

topn = 100
ranking_srges = sorted([(s, yt) for s, yt in zip(eff_sub_s['sRGES'], eff_sub_s['Active'])])
prcs_srges = [sum([r[1] for r in ranking_srges[:i]]) / i for i in range(1, topn+1)]
scr_srges = range(1, topn+1)
ratio_srges = eff_sub_s['Active'].mean()
ranking_zrges = sorted([(s, yt) for s, yt in zip(eff_sub_z['Z_RGES'], eff_sub_z['Active'])])
prcs_zrges = [sum([r[1] for r in ranking_zrges[:i]]) / i for i in range(1, topn+1)]
scr_zrges = range(1, topn+1)
ratio_zrges = eff_sub_z['Active'].mean()

FIG = plt.figure(figsize=FS)
plt.plot(scr_srges, prcs_srges, color='cornflowerblue', label='sRGES %s OCTAD cpds'%eff_sub_s.shape[0])
plt.plot(scr_zrges, prcs_zrges, color='black', label='Z-RGES %s other cpds'%eff_sub_z.shape[0])
plt.plot([0.0, topn], [ratio_srges, ratio_srges], '--', color='cornflowerblue', \
         label='OCTAD cpds pos ratio')
plt.plot([0.0, topn], [ratio_zrges, ratio_zrges], '--', color='black', \
         label='Other cpds pos ratio')
plt.text(topn*0.6, 0.5, 'PRISM Huh7', size=12, weight='bold')
plt.legend()
plt.xlabel('Rank')
plt.ylabel('Hit Rate')
plt.ylim(0.0, 1.05)
plt.tight_layout()
#FIG.savefig('./data_for_code/HCCinSilico_HitRate_zrges_srges.pdf', transparent=True)










