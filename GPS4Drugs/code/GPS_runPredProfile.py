# -*- coding: utf-8 -*-
"""
Created on Tue May 25 16:15:21 2021

@author: jingxing
"""

import os
import sys
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--cmpd_input', type=str, help='Input csv for compound ID and SMILES', default='../data/input_cmpd_gene/cmpd__TestJob0.csv')
parser.add_argument('--gene_input', type=str, help='Input list of gene symbols', default='preselect')
parser.add_argument('--cpu_num', type=int, help='Number of cpu cores to use', default=10)
# parser.add_argument('--email', type=str, help='The user\'s email address for receiving the result', default='xingjin1@msu.edu')
args = parser.parse_args()

prefix = os.path.basename(args.cmpd_input).replace('.csv', '').replace('cmpd__', '')

pathname = os.path.dirname(sys.argv[0])        
path = os.path.abspath(pathname)
GATE = path.split('code')[0]

path_code = GATE + 'code'
sys.path.append(path_code)


logfilename = GATE + 'logs/logfile_profilePred_%s.txt'%prefix
flog = open(logfilename, 'w')
print('Created log file: %s'%logfilename)


try:
    import pandas as pd
    from subprocess import Popen, PIPE
    import time
    from input_error_check import CheckError
    # import smtplib
    # from email.header import Header
    # from email.mime.text import MIMEText
    # from email.mime.application import MIMEApplication
    # from email.mime.multipart import MIMEMultipart
except Exception as e:
    flog.write('ERROR: environment setting failed.\n')
    flog.write(str(e) + '\n')
    print(e)
    flog.close()
    raise(e)


# Clean folders
CLS = ['HEPG2_t0', 'MCF7_t1', 'PC3_t1', 'VCAP_t1']
for cl in CLS:
    dst = GATE + 'data/profile_pred/' + cl + '/prob_' + prefix
    if os.path.exists(dst): os.system('rm -r ' + dst)
    os.mkdir(dst)


VERYBGN = time.time()


print('Checking input variables...')
try:
    warnmsg = ''
    warnmsg += CheckError('cmpd_input', args.cmpd_input)
    warnmsg += CheckError('gene_input', args.gene_input)
    warnmsg += CheckError('cpu_num', args.cpu_num)
    if warnmsg != '': flog.write(warnmsg)
except Exception as e:
    flog.write(str(e))
    if warnmsg != '': flog.write(warnmsg)
    flog.close()
    raise Exception('Got error during input check, please check log file: %s'%logfilename)


print('Assigning tasks...')
task_num, task_sep = 0, []
if args.gene_input != 'preselect':
    print('Generating gene features...')
    gl = [l.strip() for l in open(args.gene_input)]
    df = pd.read_csv(GATE + 'data/input_gene_features/go_fingerprints_allGenesExt.csv', index_col=0)
    df = df.loc[df.index.isin(gl)]
    task_num = df.shape[0]
    df.to_csv(GATE + 'data/input_gene_features/go_fingerprints_%s.csv'%prefix)
    desc = 'Generated features of %s gene(s): '%task_num + '|'.join(df.index) + '\n'
    flog.write(desc)
    print(desc)
    del df
    t = int(round(task_num / 3))
    task_sep = [(i * t, (i + 1) * t) if i + 1 < 3 else (i * t, task_num) for i in range(3)]
    gene_inp = GATE + 'data/input_gene_features/go_fingerprints_%s.csv'%prefix

prs = []
for CL in CLS:
    if args.gene_input == 'preselect':
        df = pd.read_csv(GATE + 'data/input_gene_features/go_fingerprints_2k_%s.csv'%CL)
        task_num = df.shape[0]
        del df
        t = int(round(task_num / 3))
        task_sep = [(i * t, (i + 1) * t) if i + 1 < 3 else (i * t, task_num) for i in range(3)]
        gene_inp = GATE + 'data/input_gene_features/go_fingerprints_2k_%s.csv'%CL
    for st, ed in task_sep:
        cmdl = 'python %scode/apply_RCL.py --cmpd_input %s --gene_input %s --cl %s --gene_start_idx %s --gene_end_idx %s'%(GATE, args.cmpd_input, gene_inp, CL, st, ed)
        p = Popen(cmdl.split(' '), stdout=PIPE, stderr=PIPE)
        prs.append(('%s__%s__%s'%(CL, st, ed), p))
print('All tasks are running... This might take ~10 min per thousand compounds per thousand genes.')


time_start = time.time()
while True:
    proc_status = [p.poll() for pname, p in prs]
    err = any([s != 0 and s is not None for s in proc_status])
    timeout = (time.time() - time_start >= 14400)
    done = all([s == 0 for s in proc_status])
    desc = ''
    
    if err:
        desc = 'ERROR: one or more tasks error, all killed.\n'
        for pname, p in prs:
            if p.poll() is None: p.kill()
    elif done:
        desc = 'COMPUTING DONE\n'
    elif timeout:
        desc = 'TIME OUT\n'
        for pname, p in prs:
            if p.poll() is None: p.kill()
    
    if any([err, done, timeout]):
        print(desc)
        for pname, p in prs:
            flog.write('###### LOG %s ######\n'%pname)
            flog.write(str(p.stdout.read()))
            flog.write('\n')
            flog.write(str(p.stderr.read()))
            flog.write('\n')
        flog.write(desc)
        flog.flush()
    
    if err or timeout:
        flog.close()
        raise Exception('Got error or time out. Please check the logfile: %s'%logfilename)
    elif done:
        break
        
    time.sleep(10)


print('Start merging... This might take several minutes.')
flog.write('###### START MERGING... ######\n')
cmdl = 'python %scode/prob2pred.py --cmpd_input %s --cpu_num %s'%(GATE, args.cmpd_input, args.cpu_num)
p = Popen(cmdl.split(' '), stdout=PIPE, stderr=PIPE)

time_start = time.time()
while True:
    status = p.poll()
    err, done, timeout = (status != 0) and (status is not None), status == 0, time.time() - time_start >= 3600
    desc = ''
    
    if err:
        desc = 'ERROR: merging error\n'
    elif timeout:
        desc = 'MERGE TIME OUT\n'
        p.kill()
    elif done:
        desc = 'MERGING DONE\n'
    
    if any([err, done, timeout]):
        flog.write(p.stdout.read())
        flog.write('\n')
        flog.write(p.stderr.read())
        flog.write('\n')
        flog.write(desc)
    
    if err or timeout:
        flog.close()
        raise Exception('Got error or time out. Please check the logfile: %s'%logfilename)
    elif done:
        break
    
    time.sleep(10)
    
flog.write('TOTAL TIME: %.2f seconds.\n'%(time.time() - VERYBGN))

attm = GATE + 'data/profile_pred/MEDIAN/preds_%s/%s_MEDIAN_GeneExpressionChange.csv'%(prefix, prefix)
print('Job done. Results saved to file %s'%attm)


flog.close()



    

