# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 15:37:35 2023

@author: Jing
"""

import pandas as pd
from itertools import combinations
from scipy.stats import zscore
from scipy.spatial.distance import cosine, jaccard, correlation
from umap import UMAP
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["font.family"] = "Arial"
from rdkit import Chem
from rdkit.Chem import AllChem


wrkdir = '../data/example_pathway_targets/'
meta = pd.read_csv(wrkdir + 'summary20_target_pathway_examples.csv', index_col='ID')
gfp = pd.read_csv(wrkdir + 'examplePwTar_MEDIAN_GeneExpressionChange.csv', index_col=0)
gfp = gfp.T

tmp = gfp.where(gfp != 0).count()
sele = tmp.loc[tmp > round(gfp.shape[0] * 0.05)].index
gfp = gfp[sele]


cfp = []
nb = 1024
for smi in meta['SMILES']:
    mol = Chem.MolFromSmiles(smi)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nb, useChirality=True, useFeatures=True)
    cfp.append([int(i) for i in fp.ToBitString()])
cfp = pd.DataFrame(data=cfp, index=meta.index, columns=['Feature%s'%i for i in range(nb)])

dist_mtx_g = gfp.T.corr(method=lambda a, b: cosine(a, b))
dist_mtx_c = cfp.T.corr(method=lambda a, b: jaccard(a, b))

nb = 15

tops = []
for cpd, row in dist_mtx_g.iterrows():
    top = row.sample(nb, random_state=1)
    qtar, qpw = meta.loc[cpd, 'Target'], meta.loc[cpd, 'Pathway']
    lst1, lst2 = meta.loc[top.index, 'Target'].to_list(), meta.loc[top.index, 'Pathway'].to_list()
    num_same_tar = sum([tar == qtar for tar in lst1])
    num_same_pw = sum([pw == qpw and tar != qtar for tar, pw in zip(lst1, lst2)])
    tops.append([cpd, num_same_tar, 'Same target', 'Random'])
    tops.append([cpd, num_same_pw, 'Same pathway diff. target', 'Random'])

for cpd, row in dist_mtx_c.iterrows():
    top = row.sort_values().iloc[:nb]
    qtar, qpw = meta.loc[cpd, 'Target'], meta.loc[cpd, 'Pathway']
    lst1, lst2 = meta.loc[top.index, 'Target'].to_list(), meta.loc[top.index, 'Pathway'].to_list()
    num_same_tar = sum([tar == qtar for tar in lst1])
    num_same_pw = sum([pw == qpw and tar != qtar for tar, pw in zip(lst1, lst2)])
    tops.append([cpd, num_same_tar, 'Same target', 'Chemistry'])
    tops.append([cpd, num_same_pw, 'Same pathway diff. target', 'Chemistry'])

for cpd, row in dist_mtx_g.iterrows():
    top = row.sort_values().iloc[:nb]
    qtar, qpw = meta.loc[cpd, 'Target'], meta.loc[cpd, 'Pathway']
    lst1, lst2 = meta.loc[top.index, 'Target'].to_list(), meta.loc[top.index, 'Pathway'].to_list()
    num_same_tar = sum([tar == qtar for tar in lst1])
    num_same_pw = sum([pw == qpw and tar != qtar for tar, pw in zip(lst1, lst2)])
    tops.append([cpd, num_same_tar, 'Same target', 'Gene expression'])
    tops.append([cpd, num_same_pw, 'Same pathway diff. target', 'Gene expression'])

tops = pd.DataFrame(data=tops, columns=['Cmpd', 'Value', 'Class', 'Feature type'])

FIG = plt.figure(figsize=(4, 4))
sns.violinplot(x='Feature type', y='Value', hue='Class', \
               inner='quart', split=True, data=tops)
a, b = plt.ylim()
plt.ylim(a, b + (b - a) * 0.25)
plt.ylabel('Num. in %s neighbors'%nb)
plt.tight_layout()
#plt.savefig('../..//manuscript_cpdProfilePred/doc_fig_table/chem_GE_example_tar_pw.pdf', transparent=True)



model = UMAP(random_state=1, n_neighbors=gfp.shape[0]-1, min_dist=0.5, metric='cosine')
Xe = model.fit_transform(gfp)
Xe = pd.DataFrame(data=Xe, columns=['Feature 0', 'Feature 1'], index=gfp.index)
Xe['Target'] = meta.reindex(Xe.index)['Target']
Xe['Pathway'] = meta.reindex(Xe.index)['Pathway']
FIG = plt.figure(figsize=(6, 6))
G = sns.scatterplot(x='Feature 0', y='Feature 1', data=Xe, hue='Target', \
                    style='Pathway', s=100)
a, b = plt.xlim()
plt.xlim(a, b + (b - a) * 0.4)
plt.title('Gene Expression Features')
plt.tight_layout()
#plt.savefig('../doc_fig_table/umap_GE_example_targets_pathways.pdf', transparent=True)


model = UMAP(random_state=1, n_neighbors=cfp.shape[0]-1, min_dist=0.5, metric='jaccard')
Xe = model.fit_transform(cfp.to_numpy())
Xe = pd.DataFrame(data=Xe, columns=['Feature0', 'Feature1'], index=gfp.index)
Xe['Target'] = meta.reindex(Xe.index)['Target']
Xe['Pathway'] = meta.reindex(Xe.index)['Pathway']
FIG = plt.figure(figsize=(6, 6))
G = sns.scatterplot(x='Feature0', y='Feature1', data=Xe, hue='Target', \
                    style='Pathway', s=100)
a, b = plt.xlim()
plt.xlim(a, b + (b - a) * 0.4)
plt.title('Chemical Features')
plt.tight_layout()
#plt.savefig('../doc_fig_table/umap_chem_example_targets_pathways.pdf', transparent=True)

