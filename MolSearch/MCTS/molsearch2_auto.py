import os
import sys
import argparse
import subprocess
import multiprocessing as mt
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--num_drugs', type=int, help='number of drugs in csv file', default=1)
parser.add_argument('--sample_name', type=str, help='csv filename no .csv', default='taskdipg')
parser.add_argument('--previous_goals', type=str, help='stage 1 goals', default='bbbp_rges')
parser.add_argument('--pool_cores', type=int, help='number of cores to use', default=1)
parser.add_argument('--goals', type=str, help='plogp,qed,sa,rges,bbbp use _ to connect', default='plogp_qed_sa_rges')
args = parser.parse_args()

currpathname = os.path.dirname(sys.argv[0])        
curr_path = os.path.abspath(currpathname)

def parra(argi):
    previous_goals, sample_name, num, goals, curr_path = argi
    print(num, 'started')
    stage1_path = f'{previous_goals}_stage1/{sample_name}/{previous_goals}_maxchild_5_sim_20_scalar_0.7_t_0.98_idx_{num}_seed_0.txt'
    stage1_df = pd.read_csv(curr_path + '/results_visulization/' + stage1_path, sep=' ')

    cmd = ['python', 'frag_mcts_mo_stage2.py', '--goal', goals, '--num', str(num),'--start_mols', sample_name, '--max_child', '5', '--num_sims', '20', '--group_idx', '0', '--group_length', str(len(stage1_df.index)), '--seed', '0', '--scalar', '0.7', '--stage1_result', stage1_path]
    try:
        output = subprocess.check_call(cmd)
        print(num, 'complete')
    except Exception as err:
        print('Error check nohup.out', err)

if __name__ == '__main__':
    data = []
    for num in range(args.num_drugs):
        data.append((args.previous_goals, args.sample_name, num, args.goals, curr_path))
    with mt.Pool(args.pool_cores) as p:
        print(p.map(parra, data))