import os
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
import numpy as np
np.random.seed(1)
from sklearn.metrics import mean_squared_error
#from mpi4py import MPI
from scipy.stats import pearsonr, ttest_rel
from rdkit import Chem
from rdkit.Chem import AllChem


def get_morgan_fingerprint(mol, radius, nBits, FCFP=False):
    m = Chem.MolFromSmiles(mol)
    fp = AllChem.GetMorganFingerprintAsBitVect(m, radius=radius, nBits=nBits, useFeatures=FCFP)
    fp_bits = fp.ToBitString()
    finger_print = np.fromstring(fp_bits, 'u1')-ord('0')
    return finger_print


def get_drug_fp_batch(drug_smiles, radius=3, length=1024, FCFP=False):
    fp = []
    for mol in drug_smiles:
        fp.append(get_morgan_fingerprint(mol, radius, length, FCFP))
    fp = np.array(fp)
    return fp


def pred_true_vs_random(task_gene, run_times=1, cv_fold=10):
    X = dfX.to_numpy()
    y = dfY.loc[:, task_gene].to_numpy()
    mse_true, mse_random = [], []
    r_true, r_random = [], []
    for rt in range(run_times):
        tmp = np.arange(len(y))
        np.random.shuffle(tmp)
        y_random = y[tmp]
        kf = KFold(n_splits=cv_fold, shuffle=True, random_state=rt)
        for trnidx, tstidx in kf.split(X):
            X_trn, y_trn = X[trnidx], y[trnidx]
            X_tst, y_tst = X[tstidx], y[tstidx]
            model = RandomForestRegressor(random_state=1)
            model.fit(X_trn, y_trn)
            y_pred = model.predict(X_tst)
            mse = mean_squared_error(y_tst, y_pred)
            r = pearsonr(y_tst, y_pred.flat)[0]
            mse_true.append(mse)
            r_true.append(r)
            y_trn, y_tst = y_random[trnidx], y_random[tstidx]
            model = RandomForestRegressor(random_state=1)
            model.fit(X_trn, y_trn)
            y_pred = model.predict(X_tst)
            mse = mean_squared_error(y_tst, y_pred)
            r = pearsonr(y_tst, y_pred.flat)[0]
            mse_random.append(mse)
            r_random.append(r)     
    _, mse_p = ttest_rel(mse_true, mse_random)
    mse_mean = np.mean(mse_true)
    mse_diff = mse_mean - np.mean(mse_random)
    _, r_p = ttest_rel(r_true, r_random)
    r_mean = np.mean(r_true)
    r_diff = r_mean - np.mean(r_random)
    return mse_p, mse_mean, mse_diff, r_p, r_mean, r_diff


dfX = pd.read_csv('unique_drug_profile.csv')
dfY = pd.read_csv('unique_drug_value.csv')

train_smiles_feature = get_drug_fp_batch(dfX['SMILES']).astype(np.float32)
dfX = pd.DataFrame(train_smiles_feature, index = dfX['SMILES'])
dfY.index = dfX.index


gene_lst = dfY.columns.tolist()
for gene in gene_lst:
    outfile = gene
    #if os.path.exists(outfile) and os.path.getsize(outfile) > 0: continue
    res = pred_true_vs_random(gene)
    with open(outfile, 'w') as f:
        f.write('\t'.join(map(str, res)))


collect = []
for gene in gene_lst:
    f = open(gene)
    res = f.read().split('\t')
    f.close()
    collect.append([gene] + list(map(float, res)))

result_out = []
for r in zip(collect):
    result_out.append(list(r)[0])

pd.DataFrame(result_out, columns=['Gene', 'P_MSE', 'Average_MSE', 'Delta_MSE', 'P_PearsonR', 'Average_PearsonR', 'Delta_PearsonR']).to_csv('Predictabilities.csv')

#------------------------------------------------------
#------------------------------------------------------
#draft:



            
            
            









