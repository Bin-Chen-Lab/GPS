import numpy as np
from numpy import linalg as la
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem.Fingerprints import FingerprintMols
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split, StratifiedKFold, KFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
from sklearn.neural_network import MLPRegressor
import copy


def normalized_mean_squared_error(y_true, y_pred):
    return np.power(la.norm(y_true-y_pred), 2)/np.power(la.norm(y_true), 2)


def sigmoid(x):
    sigm = 1. / (1. + np.exp(-x))
    return sigm


def load_data_batch(start_idx, batch_size, max_length, x, y):
    if start_idx + batch_size <= max_length:
        x_b = x[start_idx:(start_idx + batch_size), :]
        y_b = y[start_idx:(start_idx + batch_size), :]
    else:
        x_b = np.concatenate([x[start_idx:max_length, :], x[:start_idx + batch_size - max_length, :]], axis=0)
        y_b = np.concatenate([y[start_idx:max_length, :], y[:start_idx + batch_size - max_length, :]], axis=0)
    return x_b, y_b


def generate_data(input_dim=100, n1=1000, n2=5000, eps1=0.01, eps2=0.05, hidden_dim=32, output_dim=1, plot_flag=False):
    """
    X -> W1 -> ReLU -> W2 -> Sigmoid -> Y_TRUE -> ERR -> Y_OBS
    :param d:
    :param n1:
    :param n2:
    :param eps1:
    :param eps2:
    :param hidden_dim:
    :param output_dim:
    :return:
    """

    # seed control each time before generation
    print('Start generating data...')
    print('Configurations: input_dim={:d} n1={:d} n2={:d} eps1={:.2f} eps2={:.2f} hidden_dim={:d}'.format(input_dim, n1, n2, eps1, eps2, hidden_dim))

    # x generate from Gaussian
    np.random.seed(1)
    X1 = np.random.randn(n1, input_dim)

    np.random.seed(2)
    X2 = np.random.randn(n2, input_dim)

    # weight from Gaussian
    np.random.seed(3)
    W1 = 0.1 * np.random.randn(input_dim, hidden_dim)

    np.random.seed(4)
    W2 = 0.1 * np.random.randn(hidden_dim, output_dim)

    # error from Gaussian
    np.random.seed(5)
    E1 = eps1 * np.random.randn(n1, 1)

    np.random.seed(6)
    E2 = eps2 * np.random.randn(n2, 1)

    if plot_flag:
        ax = sns.distplot(E1, label='eps {:.2f}'.format(eps1))
        ax.legend()
        ax = sns.distplot(E2, label='eps {:.2f}'.format(eps2))
        ax.legend()
        plt.show()

    Y1 = np.matmul(X1, W1)
    Y1[Y1 < 0] = 0
    Y1 = np.matmul(Y1, W2)
    Y1 = sigmoid(Y1)
    print('High Quality True: [{:.2f}, {:.2f}]'.format(np.min(Y1), np.max(Y1)))
    print('High Quality Error Range: [{:.2f}, {:.2f}]'.format(np.min(E1), np.max(E1)))

    Y1_true = copy.deepcopy(Y1)
    Y1 = Y1 + E1
    err_sig_ratio_h = normalized_mean_squared_error(Y1_true, Y1)
    print('Error-Signal-Ratio {:.4f}'.format(err_sig_ratio_h))

    Y2 = np.matmul(X2, W1)
    Y2[Y2 < 0] = 0
    Y2 = np.matmul(Y2, W2)
    Y2 = sigmoid(Y2)
    print('Low Quality True: [{:.2f}, {:.2f}]'.format(np.min(Y2), np.max(Y2)))
    print('Low Quality Error Range: [{:.2f}, {:.2f}]'.format(np.min(E2), np.max(E2)))

    Y2_true = copy.deepcopy(Y2)
    Y2 = Y2 + E2
    err_sig_ratio_l = normalized_mean_squared_error(Y2_true, Y2)
    print('Error-Signal-Ratio {:.4f}'.format(err_sig_ratio_l))

    # sort according to err
    idx1 = np.argsort(abs(E1.ravel()))
    idx2 = np.argsort(abs(E2.ravel()))

    E1, E2 = E1[idx1], E2[idx2]
    X1, Y1, Y1_true = X1[idx1], Y1[idx1], Y1_true[idx1]
    X2, Y2, Y2_true = X2[idx2], Y2[idx2], Y2_true[idx2]

    low_in_low_idx = np.where(abs(E2) > np.max(abs(E1)))[0]
    print('Number of High: {:d} Number of Low in Low: {:d}'.format(n1, low_in_low_idx.shape[0]))

    high_in_low_idx = np.where(abs(E2) <= np.max(abs(E1)))[0]
    print('Number of High: {:d} Number of High in Low: {:d}'.format(n1, high_in_low_idx.shape[0]))

    return X1, Y1_true, Y1, err_sig_ratio_h, X2, Y2_true, Y2, err_sig_ratio_l, low_in_low_idx, high_in_low_idx, E1, E2



# def generate_data(input_dim=100, n1=1000, n2=5000, eps1=0.01, eps2=0.05, hidden_dim=32, output_dim=1, plot_flag=False):
#     """
#     X -> W1 -> ReLU -> W2 -> Sigmoid -> Y_TRUE -> ERR -> Y_OBS
#     :param d:
#     :param n1:
#     :param n2:
#     :param eps1:
#     :param eps2:
#     :param hidden_dim:
#     :param output_dim:
#     :return:
#     """
#
#     # seed control each time before generation
#     print('Start generating data...')
#     print('Configurations: input_dim={:d} n1={:d} n2={:d} eps1={:.2f} eps2={:.2f} hidden_dim={:d}'.format(input_dim, n1, n2, eps1, eps2, hidden_dim))
#
#     # x generate from Gaussian
#     np.random.seed(1)
#     X1 = np.random.randn(n1, input_dim)
#
#     np.random.seed(2)
#     X2 = np.random.randn(n2, input_dim)
#
#     # weight from Gaussian
#     np.random.seed(3)
#     W1 = 0.1 * np.random.randn(input_dim, hidden_dim)
#
#     np.random.seed(4)
#     W2 = 0.1 * np.random.randn(hidden_dim, output_dim)
#
#     # error from Gaussian
#     np.random.seed(5)
#     E1 = eps1 * np.random.randn(n1, 1)
#
#     np.random.seed(6)
#     E2 = eps2 * np.random.randn(n2, 1)
#
#     if plot_flag:
#         ax = sns.distplot(E1, label='eps {:.2f}'.format(eps1))
#         ax.legend()
#         ax = sns.distplot(E2, label='eps {:.2f}'.format(eps2))
#         ax.legend()
#         plt.show()
#
#     Y1 = np.matmul(X1, W1)
#     Y1[Y1 < 0] = 0
#     Y1 = np.matmul(Y1, W2)
#     Y1 = sigmoid(Y1)
#     print('High Quality True: [{:.2f}, {:.2f}]'.format(np.min(Y1), np.max(Y1)))
#     print('High Quality Error Range: [{:.2f}, {:.2f}]'.format(np.min(E1), np.max(E1)))
#
#     Y1_true = Y1
#     Y1 = Y1 + E1
#     err_sig_ratio_h = normalized_mean_squared_error(Y1_true, Y1)
#     print('Error-Signal-Ratio {:.4f}'.format(err_sig_ratio_h))
#
#     Y2 = np.matmul(X2, W1)
#     Y2[Y2 < 0] = 0
#     Y2 = np.matmul(Y2, W2)
#     Y2 = sigmoid(Y2)
#     print('Low Quality True: [{:.2f}, {:.2f}]'.format(np.min(Y2), np.max(Y2)))
#     print('Low Quality Error Range: [{:.2f}, {:.2f}]'.format(np.min(E2), np.max(E2)))
#
#     Y2_true = Y2
#     Y2 = Y2 + E2
#     err_sig_ratio_l = normalized_mean_squared_error(Y2_true, Y2)
#     print('Error-Signal-Ratio {:.4f}'.format(err_sig_ratio_l))
#
#     low_in_low_idx = np.where(abs(E2) > np.max(abs(E1)))[0]
#     print('Number of High: {:d} Number of Low in Low: {:d}'.format(n1, low_in_low_idx.shape[0]))
#
#     high_in_low_idx = np.where(abs(E2) <= np.max(abs(E1)))[0]
#     print('Number of High: {:d} Number of High in Low: {:d}'.format(n1, high_in_low_idx.shape[0]))
#
#     return X1, Y1_true, Y1, err_sig_ratio_h, X2, Y2_true, Y2, err_sig_ratio_l, low_in_low_idx, high_in_low_idx


def trn_and_val(x_trn, y_trn, random_seed=1, val_size=0.1):
    x_trn, x_val, y_trn, y_val = train_test_split(x_trn, y_trn, test_size=val_size, random_state=random_seed)
    return x_trn, x_val, y_trn, y_val


def eval_metrcis(y_true, y_pred):
    """ y_true & y_pred shape: (n, 1) """
    y_true, y_pred = y_true.ravel(), y_pred.ravel()
    # print(y_pred)

    mse = mean_squared_error(y_true, y_pred)
    cor = pearsonr(y_true, y_pred)[0]
    return mse, cor


def pick_best_performance(p):
    """
    :param p:
    :return:
    Note: nanargmax will ignore "nan" correlation caused by constant array.
    """

    p_avg = np.mean(p, axis=0)
    p_best = p_avg[np.nanargmax(p_avg[:, 1]), :]
    print('===== 10-fold performance =====')
    print(p)
    print('===== 10-fold average for different reg =====')
    print(p_avg)
    print('===== 10-fold average best =====')
    print(p_best)
    return p_avg, p_best


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


def top_10_hit(y_true, y_pred):
    down_idx = np.argsort(y_pred)[:10]
    down_hit = np.sum(y_true[down_idx] < -1.5)
    down_n = np.sum(y_true < -1.5)

    up_idx = np.argsort(y_pred)[::-1][:10]
    up_hit = np.sum(y_true[up_idx] > 1.5)
    up_n = np.sum(y_true > 1.5)

    return down_hit, down_n, up_hit, up_n


def smiles_to_graph(mol):
    image = Draw.MolToImage(mol, size=(1000, 1000), kekulize=True)
    image.show()


def smiles_to_number_of_atoms(mol):
    mol = Chem.MolFromSmiles(mol)
    return len(mol.GetAtoms())




