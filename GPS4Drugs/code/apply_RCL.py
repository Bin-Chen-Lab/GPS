# -*- coding:utf-8 -*-
import os
import torch 
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
import torchvision.transforms as transforms
from model import MLP
import argparse, sys
import numpy as np
import datetime
import shutil
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

pathname = os.path.dirname(sys.argv[0])
path = os.path.abspath(pathname)
GATE = path.split('code')[0]

parser = argparse.ArgumentParser()
parser.add_argument('--cmpd_input', type=str, help='Input csv for compound ID and SMILES', default='../data/input_cmpd_gene/cmpd__TestJob0.csv')
parser.add_argument('--gene_input', type=str, help='Input csv for gene symbols and features', default='../data/input_gene_features/go_fingerprints_geneTestJob0.csv')
parser.add_argument('--result_dir', type = str, help = 'dir to save result txt files', default = GATE + 'data/profile_pred/')
parser.add_argument('--top_bn', action='store_true')
parser.add_argument('--cl', type = str, help = 'cell line', default='MCF7')
parser.add_argument('--num_workers', type=int, default=4, help='how many subprocesses to use for data loading')
parser.add_argument('--num_classes', type=int, default=3, help='number of classes')
parser.add_argument('--gene_start_idx', type=int, default=0)
parser.add_argument('--gene_end_idx', type=int, default=100000)
parser.add_argument('--seed', type=int, default=1)

args = parser.parse_args()

# Seed
def seed_everything(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)

    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.backends.cudnn.deterministic = True


# Hyper Parameters
batch_size = 128
input_size = 2131
num_classes = args.num_classes
args.top_bn = False

save_dir = args.result_dir + args.cl +'/'
if not os.path.exists(save_dir): os.mkdir(save_dir)
model_str = GATE + 'code/results/%s/multi/model.pkl'%(args.cl)
all_models = torch.load(model_str)
prefix = os.path.basename(args.cmpd_input).replace('.csv', '').replace('cmpd__', '')
save_dir_pred = args.result_dir + args.cl + '/prob_' + prefix + '/'  # change for new drug list
#if not os.path.exists(save_dir_pred): os.mkdir(save_dir_pred)


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


def main():
    # Seed
    seed_everything(args.seed)
    
    # load gene features and names
    g_df = pd.read_csv(args.gene_input, index_col=0)
    g_list = g_df.index.to_list()
    g_list = g_list[args.gene_start_idx:args.gene_end_idx]
    print('gene type: {:s} length {:d}'.format(os.path.basename(args.gene_input), len(g_list)))

    # load compounds
    c_df = pd.read_csv(args.cmpd_input, index_col='ID')
    drug_smiles = c_df['SMILES'].to_list()
    
    if len(drug_smiles) % batch_size == 0:
        n_batch = len(drug_smiles) // batch_size
    else:
        n_batch = len(drug_smiles) // batch_size + 1
    print('total drugs {:d} batch size {:d} n_batch {:d}'.format(len(drug_smiles), batch_size, n_batch))

    print('load model...')
    model = all_models['model0']
    model.eval()

    for gene_idx in range(len(g_list)):    #len(g_list)
        gene = g_list[gene_idx]
        print('Gene name: ' + gene)

        gene_feature = np.array(g_df.loc[gene].to_list()).astype(np.float32).reshape((1, -1))

        pred = np.empty((0, 3), dtype=np.float32)
        for i in range(n_batch):
            start = i * batch_size
            end = min((i+1)*batch_size, len(drug_smiles))

            drug_smiles_batch = drug_smiles[start:end]
            drug_features = get_drug_fp_batch(drug_smiles_batch).astype(np.float32)

            gene_features = np.repeat(gene_feature, end-start, axis=0)
            data = np.concatenate([drug_features, gene_features], axis=1)
            data = torch.from_numpy(data)
            data = Variable(data).cuda()

            logits = model(data)
            output = F.softmax(logits, dim=1)
            pred = np.concatenate([pred, output.data.cpu().numpy()], axis=0)

            if i%1000 == 0:
                print('finish {:d} batches'.format(i))
                fn = save_dir_pred + gene + '_probs.txt'
                np.savetxt(fn, pred, delimiter='\t')

        print(pred.shape)
        fn = save_dir_pred + gene + '_probs.txt'
        np.savetxt(fn, pred, delimiter='\t')


if __name__=='__main__':
    main()
    print('WORK DONE: %s %s %s %s'%(os.path.basename(args.cmpd_input), args.cl, \
                                    args.gene_start_idx, args.gene_end_idx))
