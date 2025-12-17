import os
import sys
import argparse
import subprocess
import multiprocessing as mt
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--num_drugs', type=int, help='number of drugs in csv file', default=1)
parser.add_argument('--sample_name', type=str, help='csv filename no .csv', default='taskdipg')
parser.add_argument('--pool_cores', type=int, help='number of cores to use', default=1)
parser.add_argument('--goals', type=str, help='plogp,qed,sa,rges,bbbp use _ to connect', default='bbbp_rges')
parser.add_argument('--sig_name', type=str, help='name of signature for rges_module', default='None')
args = parser.parse_args()

def parra(argi):
    sample_name, num, goals = argi
    print(num, 'started')
    cmd = ['python', 'frag_mcts_mo_stage1.py', '--goal', goals, '--start_mols', sample_name, '--max_child', '5', '--num_sims', '20', '--mol_idx', str(num), '--seed', '0', '--scalar', '0.7', '--sig_name', args.sig_name]
    try:
        output = subprocess.check_call(cmd)
        print(num, 'complete')
    except Exception as err:
        print('Error check nohup.out', err)

if __name__ == '__main__':
    data = []
    for num in range(args.num_drugs):
        data.append((args.sample_name, num, args.goals))
    with mt.Pool(args.pool_cores) as p:
        print(p.map(parra, data))