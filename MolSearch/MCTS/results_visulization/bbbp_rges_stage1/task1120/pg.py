import pandas as pd
df = pd.read_csv('bbbp_rges_maxchild_5_sim_20_scalar_0.7_idx_1_seed_0.txt', sep=' ')
# smiles = df.smiles.to_lsit()
print(df.shape)
df = df.drop_duplicates(subset=['smiles'])
smiles = df.smiles.to_list()
print(len(smiles))