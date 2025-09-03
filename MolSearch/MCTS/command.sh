### Instruction from Mengying

## Molsearch code location
/home/ubuntu/chenlab_deeplearning/chenlab_deeplearning_V2/sunmeng2/MCTS


## MCTS environment setup
conda create --name mcts python=3.7 pip
conda install -c conda-forge rdkit
pip install pandas seaborn pickle-mixin
pip install -U scikit-learn==0.21.3 (RF scorer of GSK3B and JNK3 requires this version of sklearn)
put chemblDB3.sqlitdb (1.5GB) under MCTS root folder

# outside of MCTS folder
git clone https://github.com/mahendra-awale/medchem_moves
cd medchem_moves
python setup.py install


## commnand
# stage 1
python frag_mcts_mo.py --goal gsk3b_jnk3 --start_mols task1 --max_child 5 --num_sims 20 --mol_idx 0 --seed 0 --scalar 0.7

# stage2
python frag_mcts_mo_stage2.py --goal gsk3b_jnk3 --constraint gsk3b_jnk3 --start_mols task1 --max_child 3 --num_sims 10 --scalar 0.7 --group_idx 0 --seed 0

# BBBP optimization
python frag_mcts_mo.py --goal bbbp --start_mols bbbp --max_child 5 --num_sims 20 --mol_idx 0 --seed 0 --scalar 0.7 > bbbp_0.txt 2>&1&
python frag_mcts_mo.py --goal bbbp --start_mols bbbp --max_child 5 --num_sims 20 --mol_idx 0 --seed 1 --scalar 0.7 > bbbp_1.txt 2>&1&
python frag_mcts_mo.py --goal bbbp --start_mols bbbp --max_child 5 --num_sims 20 --mol_idx 0 --seed 2 --scalar 0.7 > bbbp_2.txt 2>&1&
python frag_mcts_mo.py --goal bbbp --start_mols bbbp --max_child 5 --num_sims 20 --mol_idx 0 --seed 3 --scalar 0.7 > bbbp_3.txt 2>&1&
python frag_mcts_mo.py --goal bbbp --start_mols bbbp --max_child 5 --num_sims 20 --mol_idx 0 --seed 4 --scalar 0.7 > bbbp_4.txt 2>&1&

python frag_mcts_mo.py --goal bbbp_rges --start_mols bbbp --max_child 5 --num_sims 20 --mol_idx 0 --seed 0 --scalar 0.7 > bbbp_rges_0.txt 2>&1&
python frag_mcts_mo.py --goal bbbp_rges --start_mols bbbp --max_child 5 --num_sims 20 --mol_idx 0 --seed 1 --scalar 0.7 > bbbp_rges_1.txt 2>&1&
python frag_mcts_mo.py --goal bbbp_rges --start_mols bbbp --max_child 5 --num_sims 20 --mol_idx 0 --seed 2 --scalar 0.7 > bbbp_rges_2.txt 2>&1&
python frag_mcts_mo.py --goal bbbp_rges --start_mols bbbp --max_child 5 --num_sims 20 --mol_idx 0 --seed 3 --scalar 0.7 > bbbp_rges_3.txt 2>&1&
python frag_mcts_mo.py --goal bbbp_rges --start_mols bbbp --max_child 5 --num_sims 20 --mol_idx 0 --seed 4 --scalar 0.7 > bbbp_rges_4.txt 2>&1&
