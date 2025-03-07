# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 16:34:39 2023

@author: Jing
"""

import pandas as pd
from scipy.stats import zscore, ranksums
import numpy as np
#import qnorm
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["font.family"] = "Arial"


def score_rges(gene_up, gene_down, dz_de):
    dz_de = pd.DataFrame(dz_de)
    dz_de['RANK'] = dz_de.rank(ascending=False)
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
    score = ks_up - ks_down
    return score

def score_rs(gene_up, gene_down, dz_de):
    a = dz_de.loc[dz_de.index.isin(gene_up)]
    b = dz_de.loc[dz_de.index.isin(gene_down)]
    n = dz_de.loc[(~dz_de.index.isin(a.index)) & (~dz_de.index.isin(b.index))]
    if len(a) > 2:
        _, pa = ranksums(a, n, alternative='greater')
    else:
        pa = 1
    la = -np.log10(pa)
    if len(b) > 2:
        _, pb = ranksums(b, n, alternative='less')
    else:
        pb = 1
    lb = -np.log10(pb)
    return la + lb
    

meta = pd.read_csv('../data/example_target_pred/v0/summary20_target_examples.csv', index_col='ID')
cpd_profiles = pd.read_csv('../data/example_target_pred/v0/exampleTars_MEDIAN_GeneExpressionChange.csv', index_col=0)
tar_profiles = pd.read_csv('../data/example_target_pred/example_target_profiles_AveStd_12k.csv', index_col=0)
tar_profiles = tar_profiles.fillna(0)

ovlp = sorted(set.intersection(set(cpd_profiles.index), set(tar_profiles.index)))
print(len(ovlp))
cpd_profiles = cpd_profiles.loc[ovlp]


result = []
for cpd in cpd_profiles.columns:
    vec = cpd_profiles[cpd]
    gene_up = vec.loc[vec == 1].index
    gene_down = vec.loc[vec == -1].index
    tmp = []
    for tar in tar_profiles.columns:
        #if tar != 'KDM4A': continue
        vec = tar_profiles[tar]
        '''
        mu, sigma = vec.mean(), vec.std()
        vec = vec.loc[(vec <= mu - 0.25 * sigma) | (vec >= mu + 0.25 * sigma)]
        gene_up_ = list(set.intersection(set(gene_up), set(vec.index)))
        gene_down_ = list(set.intersection(set(gene_down), set(vec.index)))
        score = score_rs(gene_up_, gene_down_, vec)
        '''
        score = score_rs(gene_up, gene_down, vec)
        tmp.append(score)
    result.append([cpd] + tmp)
result = pd.DataFrame(data=result, columns=['Cmpd'] + tar_profiles.columns.to_list())
result.index = result.iloc[:, 0]
result = result.iloc[:, 1:]
print(result.shape)
#result_norm = qnorm.quantile_normalize(result, axis=0)
#result_norm = zscore(result, axis=0)
#result_norm = zscore(result_norm, axis=1)
#result_norm = result.rank(axis=1)
result_norm = result.copy()

viz = []
for cpd in result_norm.index:
    qtar = meta.loc[cpd, 'Target']
    for tar in result_norm.columns:
        score = result_norm.loc[cpd, tar]
        on = 'Yes' if tar == qtar else 'No'
        viz.append([cpd, qtar, tar, score, on])
viz = pd.DataFrame(data=viz, columns=['Cmpd', 'Target_True', 'Target', 'shRNA & Cmpd. Profile Sim.', 'Is inhibitor'])
viz = pd.concat([viz.loc[viz['Target'] == tar] for tar in sorted(result.columns)])
viz = viz.sort_values(by=['Target', 'Is inhibitor'])

tar_signf = dict()
for tar, subviz in viz.groupby(by='Target'):
    A = subviz.loc[subviz['Is inhibitor'] == 'Yes', 'shRNA & Cmpd. Profile Sim.']
    B = subviz.loc[subviz['Is inhibitor'] == 'No', 'shRNA & Cmpd. Profile Sim.']
    s, p = ranksums(A, B, alternative='greater')
    l = ''
    if p < 0.05:
        l = '*'
        if p < 0.01:
            l = '**'
            if p < 0.001:
                l = '***'
    tar_signf[tar] = l
tar_signf = [(tar, tar_signf[tar]) for tar in sorted(result.columns)]
mp = dict(tar_signf)
tarlist = sorted(set(viz['Target']))
tarlist_small = ['ABCC1', 'ABCG2', 'ADRB2', 'ALDH1A1', 'APP', 'CACNA1C', 'EHMT2', 'GABRA5', 'MAPK1', 'VDR']
tarlist_large = [t for t in tarlist if not t in tarlist_small]

clrs = {'Yes': 'cornflowerblue', 'No': 'lightgray'}
FIG = plt.figure(figsize=(5, 5))
xx = FIG.add_subplot()
sns.boxplot(x='Target', y='shRNA & Cmpd. Profile Sim.', hue='Is inhibitor', \
            fliersize=2, linewidth=1, data=viz.loc[viz['Target'].isin(tarlist_large)], \
            palette=clrs, ax=xx)
tklbs = ['%s %s'%(v, k) for k, v in tar_signf]
plt.xticks(ticks=range(len(tarlist_small)), labels=[' '.join([mp[tar], tar]) for tar in tarlist_large], rotation=45, ha='right')
sns.move_legend(xx, "upper left", frameon=False)
xx.set_ylabel('Similarity score')
xx.set_ylim(0, 30)
#xx.set_title('Similarity between cmpd. profile (GPS) & target shRNA profile (LINCS)')
plt.tight_layout()
plt.savefig('../doc_fig_table/boxplot_example_target_pred-B.pdf', transparent=True)
'''

FIG, axes = plt.subplots(nrows=10, ncols=1, figsize=(5, 6), sharex=True)
for ai, tar in enumerate(tarlist_small):
    sub = viz.loc[viz['Target'] == tar]
    sns.boxplot(data=sub, y='Is inhibitor', x='shRNA & Cmpd. Profile Sim.', fliersize=2, ax=axes[ai])
    axes[ai].set_ylabel(None)
    axes[ai].set_yticks([0.5])
    axes[ai].set_yticklabels([mp[tar] + ' ' + tar])
    axes[ai].set_xlabel(None)
    axes[ai].set_xlim(0, 12)
    #axes[ai].set_xticklabels(axes[ai].get_xticklabels(), size=8)
plt.tight_layout()
plt.show()
'''






