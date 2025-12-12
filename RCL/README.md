The RCL model trains on single cell line's drug-induced gene expression data, and predicts on the same cell line for new drug-induced gene expression values. The conda environment setting is provided [here](https://github.com/Bin-Chen-Lab/GPS/blob/main/GPS4Drugs/py27.yml).

# Step 1: Prepare your drug-induced gene expression data in your interested cell line. 
[Here](https://chenlab-data-public.s3.amazonaws.com/GPS4Drugs_RCL/LINCS_NEW_LOW_DOSE_L4.csv) is the raw data downloaded from LINCS L4 for drug-induced gene expression in various cell lines under the lose-dose condition. 

To select your interested cell line for training the RCL model and prediction, you can use the following codes in R:

```
cell = 'MCF7' #This is the cell line you use for training and prediction

LINCS_L4_low_dose <- read.csv('LINCS_NEW_LOW_DOSE_L4.csv', stringsAsFactors = F)
tmp = filter(LINCS_L4_low_dose, cell_id == cell & SMILES != 'restricted')
meta = tmp[, 1:6]
zscore = tmp[, 7:985]
save(meta, zscore, file = 'LINCS_NEW_LOW_DOSE_L4_MCF7.RData')
```
The demo generates MCF7 cell line's data. You can prepare data for other cell lines by changing the cell line name. Place the downloaded data in the RCL/code folder, in the same folder, run [this code file](https://github.com/Bin-Chen-Lab/GPS/blob/main/RCL/code/pre_process.R) for data pre-processing. This will generate two new files "unique_drug_profile.csv" and "unique_drug_value.csv" in the RCL/code folder.

It is suggested to do gene predictability estimation before training. To do this, in the RCL/code folder, run [this code file](https://github.com/Bin-Chen-Lab/GPS/blob/main/RCL/code/cal_predictability.py) for gene predictability estimation. For tutorial purpose, we provide a result file which already contains the gene predictability results, this file can be downloaded [here](https://chenlab-data-public.s3.amazonaws.com/GPS4Drugs_RCL/Predictabilities.csv).  

Then, in the RCL/code folder, run [this code file](https://github.com/Bin-Chen-Lab/GPS/blob/main/RCL/code/train_test_split.R) for training and test set splitting.

# Step 2: Training and evaluating RCL.
In the RCL/code folder, use the following bash command to train and evaluate the RCL model:
python main_multi.py --cl MCF7 --num_networks 4 --forget_rate 0.2 --stop_fuzzy 2 --fuzzy_exponent 3 --seed 1 > MCF7_net_4_fr_0.2_sf_2_fe_3_seed_1.out 2>&1

The result file "RCL/code/results/MCF7/multi/MCF7_multi_net_4_forget_rate_0.2_stop_fuzzy_2_fuzzy_exp_3_seed_1.txt" contains the performance on the test set after each training epoch.







