import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import random
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score, roc_auc_score, accuracy_score, confusion_matrix
import pickle as pkl

def data():
    df = pd.read_csv('drug_induced_GEC_cholesterol_homeostasis.csv')
    print(len(df))
    props = ['NPC1', 'INSIG1', 'HMGCS1']
    ys = df[props].to_numpy()
    print(ys.shape)
    # ys[((ys > -1.5) & (ys < 1.5))] = 1
    # ys[ys <= -1.5] = 0
    # ys[ys >= 1.5] = 2

    ys[ys < 1.0] = 0
    ys[ys >= 1.0] = 1

    smiles = df['SMILES'].tolist()
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    fps = [get_fp(mol) for mol in mols]
    fps = np.array(fps)
    print(fps.shape)
    X=fps

    y1=ys[:, 0]
    y2=ys[:, 1]
    y3=ys[:, 2]
    print(y1.shape)
    print(y2.shape)
    print(y3.shape)

    v1, c1 = np.unique(y1, return_counts=True)
    v2, c2 = np.unique(y2, return_counts=True)
    v3, c3 = np.unique(y3, return_counts=True)
    print(v1, c1)
    print(v2, c2)
    print(v3, c3)
    return X, y1, y2, y3


def RF_single(X, y, n_est=100, md=10):
    random.seed(0)
    np.random.seed(0)
    skf = StratifiedKFold(n_splits=5, random_state=0, shuffle=True)

    RES_AUC = []
    RES_ACC = []
    RES_F1 = []
    for train_index, test_index in skf.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        print(y_train)
        print(y_train.shape)

        idx_balance = down_sampling(y_train)
        X_train_b, y_train_b = X_train[idx_balance], y_train[idx_balance]
        #X_train_b, y_train_b = X_train, y_train

        clf = RandomForestClassifier(max_depth=md, n_estimators=n_est, random_state=0)
        clf.fit(X_train_b, y_train_b)

        y_pred = clf.predict(X_test)
        y_prob = clf.predict_proba(X_test)
        auc, acc, f1, cm = eval(y_test, y_pred, y_prob)
        print("testing", auc, acc, f1)
        print(cm)
        RES_AUC.append(auc)
        RES_ACC.append(acc)
        RES_F1.append(f1)

        y_pred_train = clf.predict(X_train_b)
        y_prob_train = clf.predict_proba(X_train_b)
        auc, acc, f1, cm = eval(y_train_b, y_pred_train, y_prob_train)
        print("training", auc, acc, f1)
    RES_AUC = np.array(RES_AUC)
    RES_ACC = np.array(RES_ACC)
    RES_F1 = np.array(RES_F1)

    auc_avg = np.mean(RES_AUC, axis=0)
    acc_avg = np.mean(RES_ACC, axis=0)
    f1_avg = np.mean(RES_F1, axis=0)
    print(auc_avg, acc_avg, f1_avg)
    return auc_avg, acc_avg, f1_avg


def RF_CV(X, y, name):
    n_est_list= [200]
    md_list = [20]
    RES = []
    for n_est in n_est_list:
        for md in md_list:
            print("process n_est {:d} md {:d}".format(n_est, md))
            auc_avg, acc_avg, f1_avg = RF_single(X, y, n_est=n_est, md=md)
            RES.append(np.array([auc_avg, acc_avg, f1_avg]))
    RES = np.array(RES)
    print(RES)
    np.savetxt('CV_'+name+'.txt', RES)
    return


def eval(y_true, y_pred, y_prob):
    auc = roc_auc_score(y_true, y_prob[:, 1])#multi_class='ovr' for multi class
    acc = accuracy_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred) # , average='macro' for multi class
    cm = confusion_matrix(y_true, y_pred)
    return auc, acc, f1, cm


# def down_sampling(y):
#     unique, counts = np.unique(y, return_counts=True)
#     max_idx = np.argmax(counts)
#     max_value = unique[max_idx]
#     max_counts = counts[max_idx]
#     n_select = int((np.sum(counts) - max_counts) * 0.5)
#     print('max_value, max_counts, n_select')
#     print(max_value, max_counts, n_select)
#
#     random.seed(0)
#     #print(np.where(y == max_value)[0])
#     #print(list(np.where(y == max_value)[0]))
#     idx_select = random.sample(list(np.where(y == max_value)[0]), k=n_select)
#     idx_final = np.concatenate([np.where(y == 0)[0], idx_select, np.where(y == 2)[0]])
#     return idx_final


def down_sampling(y):
    unique, counts = np.unique(y, return_counts=True)
    max_idx = np.argmax(counts)
    max_value = unique[max_idx]
    max_counts = counts[max_idx]
    n_select = int((np.sum(counts) - max_counts))
    print('max_value, max_counts, n_select')
    print(max_value, max_counts, n_select)

    random.seed(0)
    #print(np.where(y == max_value)[0])
    #print(list(np.where(y == max_value)[0]))
    idx_select = random.sample(list(np.where(y == max_value)[0]), k=n_select)
    idx_final = np.concatenate([idx_select, np.where(y == 1)[0]])
    return idx_final



def get_fp(mol, radius=3, nBits=1024):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nBits, useFeatures=False)
    fp_bits = fp.ToBitString()
    finger_print = np.array(list(map(int, fp_bits)))
    return finger_print


def CV_all_tasks():
    X, y1, y2, y3 = data()
    name_list = ['NPC1', 'INSIG1', 'HMGCS1']
    print("====== Process NPC1 =========")
    RF_CV(X, y1, 'NPC1')
    print("====== Process INSIG1 =========")
    RF_CV(X, y2, 'INSIG1')
    print("====== Process HMGCS1 =========")
    RF_CV(X, y3, 'HMGCS1')


def final_model(X, y, name, balance=False):
    random.seed(0)
    np.random.seed(0)
    clf = RandomForestClassifier(max_depth=20, n_estimators=200, random_state=0)

    if balance:
        idx_balance = down_sampling(y)
        X_train, y_train = X[idx_balance], y[idx_balance]
    else:
        X_train, y_train = X, y

    clf.fit(X_train, y_train)
    fn = 'models/'+name+'.pkl'
    with open(fn, 'wb') as f:
        pkl.dump(clf, f)
    return clf


def final_models(balance=True):
    X, y1, y2, y3 = data()
    final_model(X, y1, 'npc1', balance=balance)
    final_model(X, y2, 'insig1', balance=balance)
    final_model(X, y3, 'hmgcs1', balance=balance)


if __name__=='__main__':
    #data()
    CV_all_tasks()
    final_models()
