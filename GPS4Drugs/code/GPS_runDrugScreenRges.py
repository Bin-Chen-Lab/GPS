#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 18:33:27 2021

@author: jingxing
"""

import os
import sys
import argparse
import pandas as pd
from subprocess import Popen, PIPE
import time
from input_error_check import CheckError
# import smtplib
# from email.header import Header
# from email.mime.text import MIMEText
# from email.mime.application import MIMEApplication
# from email.mime.multipart import MIMEMultipart

# GATE defines the location of the GPS4Drugs folder
pathname = os.path.dirname(sys.argv[0])        
path = os.path.abspath(pathname)
GATE = path.split('code')[0]

parser = argparse.ArgumentParser()
parser.add_argument('--dzSigFile', type = str, default=GATE + 'data/dzsig/DZSIG__TestJobNotExists.csv')
parser.add_argument('--cmpdLibID', type = str, default='ZINC')
parser.add_argument('--cpu_num', type = int, help = 'Number of cpu cores to use', default=10)
parser.add_argument('--RGES_bgrd_ID', type=str, default='NONE')
# parser.add_argument('--email', type=str, help='The user\'s email address for receiving the result', default='xingjin1@msu.edu')
args = parser.parse_args()

# defines the path of input and output files
# potential issue here and throughout GPS is that it's implied that these folders already exist and there are no checks or code to create them
job_prefix = os.path.basename(args.dzSigFile).replace('DZSIG__', '').replace('.csv', '')
bgrd_file = None if args.RGES_bgrd_ID == 'NONE' else GATE + 'data/dzsig/BGRD__%s.pkl'%args.RGES_bgrd_ID
# cmpdLibFile = 'ZINC' if args.cmpdLibID == 'ZINC' else GATE + 'data/profile_pred/MEDIAN/preds_%s/%s_MEDIAN_GeneExpressionChange.csv'%(args.cmpdLibID, args.cmpdLibID)

# added this for other libs
cmpdLibFile = args.cmpdLibID

# starts time countdown
VERYBGN = time.time()

# creates log file
logfilename = GATE + 'logs/logfile_drugScreen_%s.txt'%job_prefix
flog = open(logfilename, 'w')
print('Created log file: %s'%logfilename)

# # checks the inputs
# print('Checking input variables...')
# try:
#     warnmsg = ''
#     warnmsg += CheckError('dzSigFile', args.dzSigFile)
#     warnmsg += CheckError('cmpdLibID', args.cmpdLibID)
#     warnmsg += CheckError('cpu_num', args.cpu_num)
#     warnmsg += CheckError('RGES_bgrd_ID', args.RGES_bgrd_ID)
#     if warnmsg != '': flog.write(warnmsg)
# except Exception as e:
#     flog.write(str(e))
#     if warnmsg != '': flog.write(warnmsg)
#     flog.close()
#     raise Exception('Got error during input check, please check log file: %s'%logfilename)

if args.RGES_bgrd_ID == 'NONE':
    # Get compound profile included gene list
    trscrptm_file = 'preselect'
    # uploads predictions from GPS_runPredProfile.py for user defined compounds (those uploaded with SMILES formulas)
    # if args.cmpdLibID != 'ZINC':
    #     cmpd_lib = pd.read_csv(cmpdLibFile, index_col=0)
    #     with open(GATE + 'data/input_cmpd_gene/gene__%s.txt'%job_prefix, 'w') as f:
    #         f.write('\n'.join(cmpd_lib.index) + '\n')
    #     trscrptm_file = GATE + 'data/input_cmpd_gene/gene__%s.txt'%job_prefix
    
    # Permute RGES
    # this is the first external script 
    cmdl = 'python %scode/permute_rges.py --dzSigFile %s --gene_input %s --cpu_num %s'%(GATE, args.dzSigFile, trscrptm_file, args.cpu_num)
    p = Popen(cmdl.split(' '), stdout=PIPE, stderr=PIPE)
    desc = 'Start RGES permutation... This might take 0.5 to 6h.\n'
    print(desc)
    flog.write(desc)
    
    # Job monitoring
    # this part picks up the completion of permute_rges
    time_start = time.time()
    while True:
        time.sleep(30)
        
        status = p.poll()
        err, done, timeout = (status != 0) and (status is not None), status == 0, time.time() - time_start >= 21600
        desc = ''
    
        if err:
            desc = 'ERROR: RGES permutation error\n'
        elif timeout:
            desc = 'PERMUTATION TIME OUT\n'
            p.kill()
        elif done:
            desc = 'PERMUTATION DONE\nTo save the permutation time for the same disease signature, please specify --RGES_bgrd_ID as %s in the future.\n'%job_prefix
    
        if any([err, done, timeout]):
            print(desc)
            flog.write(p.stdout.read())
            flog.write('\n')
            flog.write(p.stderr.read())
            flog.write('\n')
            flog.write(desc)
            flog.flush()
    
        if err or timeout:
            flog.close()
            raise Exception('Got error or time out. Please check the logfile: %s'%logfilename)
        elif done:
            bgrd_file = GATE + 'data/dzsig/BGRD__%s.pkl'%job_prefix
            break

# Calculate normalized RGES for each compound profile
cmdl = 'python %scode/Run_reversal_score.py --dzSigFile %s --cmpdLibFile %s --RGES_bgrd_file %s --cpu_num %s'%\
       (GATE, args.dzSigFile, cmpdLibFile, bgrd_file, args.cpu_num)
p = Popen(cmdl.split(' '), stdout=PIPE, stderr=PIPE)
desc = 'Start screening... This might take ~40 seconds per hundred compounds.\n'
print(desc)
flog.write(desc)

# Job monitoring
# picks up the completion of Run_reversal_score
time_start = time.time()
while True:
    time.sleep(10)
    
    status = p.poll()
    err, done, timeout = (status != 0) and (status is not None), status == 0, time.time() - time_start >= 21600
    desc = ''
    
    if err:
        desc = 'ERROR: drug screening error\n'
    elif timeout:
        desc = 'SCREEN TIME OUT\n'
        p.kill()
    elif done:
        desc = 'SCREEN DONE\n'
    
    if any([err, done, timeout]):
        print(desc)
        flog.write(p.stderr.read())
        flog.write('\n')
        flog.write(p.stdout.read())
        flog.write('\n')
        flog.write(desc)
        flog.flush()
    
    if err or timeout:
        flog.close()
        raise Exception('Got error or time out. Please check the logfile: %s'%logfilename)
    elif done:
        break

flog.write('TOTAL TIME: %.2f seconds.\n'%(time.time() - VERYBGN))


# replaced the email sending with path to output file since not all servers allow connection to SMTP
attm = GATE + 'data/reversal_score/%s_RGES_norm.csv'%job_prefix
print('Job done. Results saved to file %s'%attm)


# print('Job done. Sending an email to the user...')
# smtp_server, from_addr, passwd = 'smtp.office365.com', 'contact@octad.org', 'Chenlab2018'
# txt = 'Your submission \'%s\' on GPS4Drugs is done. Here\'s the result file.\n\nBest,\nJing Xing, PhD\nGPS4Drugs Team\n'%job_prefix
# attm = GATE + 'data/reversal_score/%s_RGES_norm.csv'%job_prefix
# try:
#     msg = MIMEMultipart('mixed')
#     text_load = MIMEText(txt,'plain','utf-8')
#     msg.attach(text_load)
#     msg['Subject'] = Header('GPS4Drugs ' + job_prefix, charset='utf-8')
#     with open(attm, 'rb') as fat:
#         attm_load = MIMEApplication(fat.read())
#     attm_load.add_header('Content-Disposition', 'attachment', filename=os.path.basename(attm))
#     msg.attach(attm_load)
#     server = smtplib.SMTP(host=smtp_server, port=587)
#     server.starttls()
#     server.login(user=from_addr, password=passwd)
#     server.sendmail(from_addr=from_addr, to_addrs=args.email, msg=msg.as_string())
#     server.quit()
# except Exception as e:
#     flog.write(str(e) + '\n')
#     print('Got error sending the email. Please check the logfile: %s'%logfilename)
#     flog.close()
#     raise e
# print('Sent an email notification.')


flog.close()




