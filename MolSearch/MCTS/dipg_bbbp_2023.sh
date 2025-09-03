### use the original get_mo_actions function to select moveable actions
python frag_mcts_mo_stage1.py --goal bbbp_rges_sa_qed --start_mols dipg_bbbp --max_child 30 --mol_idx 0
python frag_mcts_mo_stage1.py --goal bbbp_rges_sa_qed --start_mols dipg_bbbp --max_child 30 --mol_idx 1
### changed to use get_constraint_actions
python frag_mcts_mo_stage1.py --goal bbbp_rges_sa_qed --start_mols dipg_bbbp --max_child 30 --mol_idx 0
python frag_mcts_mo_stage1.py --goal bbbp_rges_sa_qed --start_mols dipg_bbbp --max_child 30 --mol_idx 1
### as it takes too long to run the simulation, I modified the restrictions
python frag_mcts_mo_stage1.py --goal bbbp_rges_sa_qed --start_mols dipg_bbbp --max_child 30 --mol_idx 0 --t 0.99
python frag_mcts_mo_stage1.py --goal bbbp_rges_sa_qed --start_mols dipg_bbbp --max_child 30 --mol_idx 1 --t 0.99

### as it takes too long to run the simulation, I modified the restrictions
python frag_mcts_mo_stage1.py --goal bbbp_rges_sa_qed --start_mols dipg_bbbp --max_child 10 --mol_idx 0 --t 0.99
python frag_mcts_mo_stage1.py --goal bbbp_rges_sa_qed --start_mols dipg_bbbp --max_child 30 --mol_idx 1 --t 0.99


### molecule 1 was not predicted to have good bbbp and RGES properties.
python frag_mcts_mo_stage1.py --goal bbbp_rges_sa_qed --start_mols dipg_bbbp --max_child 50 --mol_idx 1 --t 0.99