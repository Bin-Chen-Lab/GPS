# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 20:46:23 2023

@author: Jing
"""

import pandas as pd
import os
from scipy.stats import spearmanr
from itertools import combinations, product
import random
random.seed(0)
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["font.family"] = "Arial"


wrkdir = './data_for_code/cmpd_profile_robustness/'
profiles_dir = wrkdir + 'LINCS2020_Level5_gold_cmpd_24h_10uM/'
cpd_cell_smry = pd.read_csv(wrkdir + 'cmpd_cell_line.csv', index_col=0)
cpd_moas = pd.read_csv(wrkdir + 'cmpd_MoAs.csv', index_col=0)
gps_cells = ['MCF7', 'PC3', 'VCAP', 'HEPG2']
#gps_cells = ['MCF7', 'PC3', 'HEPG2']


df_cells = dict()
for fn in os.listdir(profiles_dir):
    if not fn.endswith('.csv'): continue
    cell = fn.replace('.csv', '').split('_')[-1]
    tmp = pd.read_csv(profiles_dir + fn, index_col=0)
    cnt = tmp.where(tmp.abs() >= 2).count()
    tmp = tmp.sort_index()
    df_cells[cell] = tmp
    print(cell, tmp.shape)
    
moa_dict = dict([(cpd, set(s.split('|'))) for cpd, s in \
                 zip(cpd_moas.index, cpd_moas['MoAs'])])

profile_corr_gps = []
for cell1, cell2 in combinations(gps_cells, 2):
    print(cell1, cell2)
    dfa, dfb = df_cells[cell1], df_cells[cell2]
    ovlp = set.intersection(set(dfa.columns), set(dfb.columns))
    for cpd in ovlp:
        cor, p = spearmanr(dfa[cpd], dfb[cpd])
        profile_corr_gps.append([cor, 'Same cmpd. %s-%s'%(cell1, cell2)])
print(len(profile_corr_gps))

rn = 1000
for cell, df in df_cells.items():
    if not cell in gps_cells: continue
    candt = set.intersection(set(moa_dict.keys()), set(df.columns))
    combos = [c for c in combinations(candt, 2)]
    combos = [(a, b) for (a, b) in combos if len(set.intersection(moa_dict[a], moa_dict[b])) == 0]
    if len(combos) > rn:
        combos = random.sample(combos, rn)
    corrs = [spearmanr(df[a], df[b])[0] for a, b in combos]
    profile_corr_gps += [[cor, 'Diff. cmpd. within %s'%cell] for cor in corrs]
    print(cell, len(corrs))

profile_corr_gps = pd.DataFrame(data=profile_corr_gps, columns=['LINCS profile correlation', 'Group'])
render_order = ['Same cmpd. MCF7-HEPG2', 'Same cmpd. MCF7-PC3', 'Same cmpd. PC3-HEPG2', \
                'Same cmpd. PC3-VCAP', 'Same cmpd. MCF7-VCAP', 'Same cmpd. VCAP-HEPG2', \
                'Diff. cmpd. within HEPG2', 'Diff. cmpd. within MCF7', 'Diff. cmpd. within PC3', \
                'Diff. cmpd. within VCAP']

clr = dict([(k, 'royalblue' if k.startswith('Same') else 'orange') for k in \
            profile_corr_gps['Group'].unique()])
nums = [profile_corr_gps.loc[profile_corr_gps['Group'] == grp].shape[0] for grp in render_order]
FIG = plt.figure(figsize=(6, 6))
xx = plt.subplot()
sns.boxenplot(x='LINCS profile correlation', y='Group', order=render_order, \
              palette=clr, data=profile_corr_gps, ax=xx)
xx.set_xlim(-0.9, 0.9)
xx.set_ylabel(None)
xx.set_yticklabels(['%s\n(%s pairs)'%(a, b) for a, b in zip(render_order, nums)])
plt.tight_layout()
#plt.savefig('./data_for_code/boxplot_cmpd_profiles_correlation_A.pdf', transparent=True)


cpd_gps_cell = cpd_cell_smry[gps_cells]
cpd_other_cell = cpd_cell_smry.loc[:, ~cpd_cell_smry.columns.isin(gps_cells)]
cpd_gps_cell = cpd_gps_cell.loc[cpd_gps_cell.sum(axis=1) >= len(gps_cells) - 1]
cpd_other_cell = cpd_other_cell.loc[cpd_other_cell.sum(axis=1) >= 5]
ovlp = set.intersection(set(cpd_gps_cell.index), set(cpd_other_cell.index))

profile_gps_sum = []
for cpd, row in cpd_gps_cell.iterrows():
    cells = row[row == 1].index
    cpd_df = pd.concat([df_cells[c][cpd] for c in cells], axis=1)
    cpd_sum = cpd_df.mean(axis=1) / cpd_df.std(axis=1)
    cpd_sum.name = cpd
    profile_gps_sum.append(cpd_sum)
profile_gps_sum = pd.concat(profile_gps_sum, axis=1)

profile_other_sum = []
for cpd, row in cpd_other_cell.iterrows():
    cells = row[row == 1].index
    cpd_df = pd.concat([df_cells[c][cpd] for c in cells], axis=1)
    cpd_sum = cpd_df.mean(axis=1) / cpd_df.std(axis=1)
    cpd_sum.name = cpd
    profile_other_sum.append(cpd_sum)
profile_other_sum = pd.concat(profile_other_sum, axis=1)

profile_corr_sum = []
for cpd, row in cpd_gps_cell.iterrows():
    cells = row[row == 1].index
    cpd_profile_sum = profile_gps_sum[cpd]
    profile_corr_sum += [[spearmanr(cpd_profile_sum, df_cells[c][cpd])[0], \
                          'Same cmpd. aggr. vs. %s'%c] for c in cells]

for cpd in ovlp:
    cor, p = spearmanr(profile_gps_sum[cpd], profile_other_sum[cpd])
    profile_corr_sum.append([cor, 'Same cmpd. aggr. GPS\nvs. other cells'])

combos = list(product(profile_gps_sum.columns, profile_other_sum.columns))
combos = [c for c in combos if all([x in moa_dict for x in c])]
combos = [(a, b) for a, b in combos if len(set.intersection(moa_dict[a], moa_dict[b])) == 0]
if len(combos) > rn:
    combos = random.sample(combos, rn)
profile_corr_sum += [[spearmanr(profile_gps_sum[a], profile_other_sum[b])[0], \
                      'Diff. cmpd. aggr. GPS\nvs. other cells'] for a, b in combos]

profile_corr_sum = pd.DataFrame(profile_corr_sum, columns=['LINCS profile correlation', 'Group'])

render_order = ['Same cmpd. aggr. vs. MCF7', 'Same cmpd. aggr. vs. PC3', \
                'Same cmpd. aggr. vs. HEPG2', 'Same cmpd. aggr. vs. VCAP', \
                'Same cmpd. aggr. GPS\nvs. other cells', 'Diff. cmpd. aggr. GPS\nvs. other cells']
clr = dict(zip(render_order, ['royalblue' for c in gps_cells] + ['forestgreen', 'orange']))
nums = [profile_corr_sum.loc[profile_corr_sum['Group'] == grp].shape[0] for grp in render_order]
FIG = plt.figure(figsize=(6, 6))
xx = plt.subplot()
sns.boxenplot(x='LINCS profile correlation', y='Group', order=render_order, \
              palette=clr, data=profile_corr_sum, ax=xx)
xx.set_xlim(-1, 1)
xx.set_ylabel(None)
xx.set_yticklabels(['%s\n(%s pairs)'%(a, b) for a, b in zip(render_order, nums)])
plt.tight_layout()
#plt.savefig('./data_for_code/boxplot_cmpd_profiles_correlation_B.pdf', transparent=True)








