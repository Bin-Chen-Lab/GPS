import warnings
warnings.simplefilter("ignore", category=UserWarning)
warnings.simplefilter("ignore", category=DeprecationWarning)
warnings.simplefilter("ignore", category=FutureWarning)
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.warning')

import torch
device = "cuda" if torch.cuda.is_available() else "cpu"
from torch.autograd import Variable
import torch.nn.functional as F

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import QED

import numpy as np
import pickle as pkl
import pandas as pd
import os
from pathlib import Path

from scipy.stats import zscore

class RGESCalculator():

    def __init__(self, gen_df_path, **kwargs):
        super(RGESCalculator, self).__init__(**kwargs)
        self.gen_df_path = gen_df_path
        self.gene_names, self.gene_features, self.n_genes = self._get_gene_df() # gene_features are pre_loaded
        self.gene_na_idx = self._get_gene_na_idx()
        self.current_fp = Path(__file__)
        self.models = self._load_HCC_models()   # models are pre_loaded, fix seed and evaluation mode

    def rges_score(self, mol):
        fp = self._get_morgan_fingerprint(mol).astype(np.float32).reshape(1, -1)
        regulate_vector = self._merge_predicted_profiles(fp)
        reg = pd.DataFrame({"gene_names":self.gene_names, 'reg_effect': regulate_vector})
        reg = reg.set_index("gene_names")
        reg = reg.T
        # print("reg df: ", reg)
        # print("reg size: ", reg.shape)

        ## load disease signature and background
        dzde = pd.read_csv(f'libs/rges_input/DZSIG__{self.gen_df_path}.csv', index_col='GeneSymbol')
        dzde['RANK'] = dzde['Value'].rank(ascending=False)
        
        # permt_bgrd: the background distribution of raw rges values for a specific dysregulated gene number of a drug.
        permt_bgrd = pd.read_csv(f'libs/rges_input/BGRD__{self.gen_df_path}.csv', index_col=0)

        # get profile
        profile = reg.iloc[0]
        up_genes = set(profile.loc[profile == 1].index)
        down_genes = set(profile.loc[profile == -1].index)

        # Find the background of a given number of up/down genes
        bins_up = [int(b.split('_')[1]) for b in permt_bgrd.index if b.startswith('Up')]
        bins_down = [int(b.split('_')[1]) for b in permt_bgrd.index if b.startswith('Down')]
        mark_up = min(bins_up, key=lambda x: abs(x - len(up_genes)))
        mark_down = min(bins_down, key=lambda x: abs(x - len(down_genes)))
        bgrd_up = permt_bgrd.loc['Up_%s'%mark_up].to_list()
        bgrd_down = permt_bgrd.loc['Down_%s'%mark_down].to_list()
        rges_bgrd = (bgrd_up, bgrd_down)

        # calculate score
        s = self._score_norm_rges(up_genes, down_genes, dzde, rges_bgrd)
        return s

    def _score_norm_rges(self, gene_up, gene_down, dz_de, bgrd):
        # Raw RGES calculation
        gene_up = list(gene_up)
        gene_down = list(gene_down)
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
        
        # Normalization
        bgrd_up, bgrd_down = bgrd
        z_up = zscore([ks_up] + bgrd_up)[0]
        z_down = zscore([ks_down] + bgrd_down)[0]
        score = z_up - z_down
        
        return score

    def _merge_predicted_profiles(self, fp, threshold=0.95):

        HEPG2_NA_IDX, MCF7_NA_IDX, PC3_NA_IDX, VCAP_NA_IDX = self.gene_na_idx

        probs_HEPG2 = self._predicted_profiles(fp, "HEPG2", HEPG2_NA_IDX)
        probs_MCF7 = self._predicted_profiles(fp, "MCF7", MCF7_NA_IDX)
        probs_PC3 = self._predicted_profiles(fp, "PC3", PC3_NA_IDX)
        probs_VCAP = self._predicted_profiles(fp, "VCAP", VCAP_NA_IDX)

        probs = np.stack([probs_HEPG2, probs_MCF7, probs_PC3, probs_VCAP], axis=0)  # 4 x n_genes x 3
        merged_probs = np.nanmedian(probs, axis=0)  # n_genes x 3

        # up/down regulate/no change accoding to 3 colums
        # merged_probs[:, 0] > threshold, -1
        # merged_probs[:, 2] > threshold, 1
        # else, 0
        ans = np.zeros(merged_probs.shape[0], dtype=np.int8)
        ans[merged_probs[:, 0] > threshold] = -1
        ans[merged_probs[:, 2] > threshold] = 1

        return ans

    def _predicted_profiles(self, fp, cl, na_idx):
        # get predicted profiles for a given mol
        model = self.models['model_'+cl]

        drug_features = np.repeat(fp, self.n_genes, axis=0)
        # print("drug: ", drug_features.shape)
        # print("gene: ", self.gene_features.shape)

        data = np.concatenate([drug_features, self.gene_features], axis=1)
        data = torch.from_numpy(data)
        data = Variable(data).to(device)

        logits = model(data)
        output = F.softmax(logits, dim=1)
        probs = output.data.cpu().numpy()   # n_genes x 3

        # set to na if gene is not predicatable
        probs[na_idx, :] = np.nan
        return probs

    def _seed_everything(self, seed):
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        np.random.seed(seed)

        os.environ['PYTHONHASHSEED'] = str(seed)
        torch.backends.cudnn.deterministic = True

    def _get_morgan_fingerprint(self, mol, radius=3, nBits=1024):
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nBits, useFeatures=False)
        fp_bits = fp.ToBitString()
        #finger_print = np.fromstring(fp_bits, 'u1')-ord('0')
        finger_print = np.array(list(map(int, fp_bits)))
        return finger_print

    def _get_gene_na_idx(self):
        with open(f'libs/rges_input/NA_IDX_{self.gen_df_path}.pkl', 'rb') as f:
            gene_na_idx = pkl.load(f)
        return gene_na_idx

    def _get_gene_df(self):
        gene_df = pd.read_csv(f'libs/rges_input/{self.gen_df_path}_gene_feature.csv', index_col=0)
        gene_names = gene_df.index
        gene_features = gene_df.to_numpy().astype(np.float32)
        nrow, ncol = gene_features.shape
        # print("# genes: ", len(gene_names))
        # print("gene dim: ", gene_features.shape)
        return gene_names, gene_features, nrow

    def _load_HCC_models(self):
        self._seed_everything(1)

        cell_lines = ['HEPG2', 'MCF7', 'PC3', 'VCAP']
        models = dict()
        for cl in cell_lines:
            model_str = self.current_fp.parent / "models/{:s}_model.pkl".format(cl)
            all_models = torch.load(model_str, map_location=device)

            model = all_models['model0']
            model.eval()

            models['model_'+cl] = model
            print("model {:s} successfully loaded!".format(str(model_str)))
        return models


if __name__ == '__main__':
    rges_calculator=RGESCalculator()

 

