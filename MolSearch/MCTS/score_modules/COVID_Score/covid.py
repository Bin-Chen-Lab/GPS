import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

class COVIDCalculator():

    def __init__(self, **kwargs):
        super(COVIDCalculator, self).__init__(**kwargs)
        self.models = self._load_model()

    def _load_model(self):
        genes = ['npc1', 'insig1', 'hmgcs1']
        models = dict()
        for gene in genes:
            fn = 'score_modules/COVID_Score/models/'+gene+'.pkl' #score_modules/COVID_Score/
            with open(fn, 'rb') as f:
                model =  pickle.load(f)
                models[gene] = model
        return models

    def _get_morgan_fingerprint(self, mol, radius=3, nBits=1024): # original 2
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nBits, useFeatures=False)
        fp_bits = fp.ToBitString()
        #finger_print = np.fromstring(fp_bits, 'u1')-ord('0')
        finger_print = np.array(list(map(int, fp_bits)))
        return finger_print

    def get_score(self, mol, gene):
        fp = self._get_morgan_fingerprint(mol).astype(np.float32).reshape(1, -1)
        model = self.models[gene]
        score = model.predict_proba(fp)[0, 1]
        return score


if __name__ == '__main__':
    s = 'C1=CC(=C(C=C1[N+](=O)[O-])Cl)NC(=O)C2=C(C=CC(=C2)Cl)O'
    mol = Chem.MolFromSmiles(s)
    covid_calculator=COVIDCalculator()
    print(covid_calculator.get_score(mol, 'npc1'))
    print(covid_calculator.get_score(mol, 'insig1'))
    print(covid_calculator.get_score(mol, 'hmgcs1'))
