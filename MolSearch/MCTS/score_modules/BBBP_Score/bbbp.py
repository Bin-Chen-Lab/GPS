import os
from pathlib import Path
import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
import loader
from loader import mol_to_graph_data_obj_simple
from model import GNN_graphpred
import pandas as pd


class BBBPCalculator():
    def __init__(self, **kwargs):
        super(BBBPCalculator, self).__init__(**kwargs)
        self.current_fp = Path(__file__)
        self.model = self._load_BBBP_model()  # models are pre_loaded, fix seed and evaluation mode

    def bbbp_score(self, mol):
        data = mol_to_graph_data_obj_simple(mol)
        # print(data.x.shape[0])
        list = [0] * data.x.shape[0]
        batch = torch.tensor(np.array(list), dtype=torch.long)
        with torch.no_grad():
            pred = self.model(data.x, data.edge_index, data.edge_attr, batch)
        # print(type(pred.cpu().numpy()[0][0]))
        # print(pred.cpu().numpy()[0][0])
        pred = pred.cpu().numpy()[0][0]
        score = 1 / (1 + np.exp(-pred))
        # print(type(score))
        return score

    def _seed_everything(self, seed):
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        np.random.seed(seed)
        os.environ['PYTHONHASHSEED'] = str(seed)
        torch.backends.cudnn.deterministic = True
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)

    def _load_BBBP_model(self):
        self._seed_everything(1)
        model_str = self.current_fp.parent/"model/bbbp_semi_1.0_protocol_nonlinear_aug1_DK1_0.2_aug2_mask_edge_0.2_lamb_0.0_do_0.5_seed_0_1.pth"
        # print(model_str)
        model = GNN_graphpred(3, 512, 1, JK="last", drop_ratio=0.5, graph_pooling="mean", gnn_type="gin")
        model.load_state_dict(torch.load(model_str))
        model.eval()
        print("model {:s} successfully loaded!".format(str(model_str)))
        return model


if __name__ == '__main__':
    bbbp_calculator = BBBPCalculator()
    s = 'CC1=C2COC(=O)C2=C(C(=C1OC)CC=C(C)CCC(=O)O)O'
    # s = 'C=CC(=O)Nc1cccc(Oc2nc(Nc3ccc(N4CCN(C)CC4)cc3)ncc2Cl)c1'
    # s = 'Oc1c(Cl)cc(Cl)cc1/C=N/Nc1cccc2cccnc12'
    print(s)
    mol = AllChem.MolFromSmiles(s)
    score = bbbp_calculator.bbbp_score(mol)
    print(score)
    df = pd.read_csv("dipg_ZINC_top2k_gnn.csv")
    smiles = df.SMILES.tolist()
    scores = []
    for s in smiles:
        mol = AllChem.MolFromSmiles(s)
        score = bbbp_calculator.bbbp_score(mol)
        scores.append(score)
    df['score_mcts'] = scores
    df.to_csv("dipg_ZINC.csv")







