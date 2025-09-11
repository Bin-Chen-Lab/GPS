The RCL model trains on single cell line's drug-induced gene expression data, and predicts on the same cell line for new drug-induced gene expression values.

# Step 1: Prepare your drug-induced gene expression data in your interested cell line. 
Here we provide demo data of drug-induced gene expression in the MCF7 cell line (the data only serves for tutorial purpose), which can be downloaded: [https://chenlab-data-public.s3.amazonaws.com/GPS4Drugs_RCL/LINCS_NEW_LOW_DOSE_L4_A549.RData](https://chenlab-data-public.s3.amazonaws.com/GPS4Drugs_RCL/LINCS_NEW_LOW_DOSE_L4_MCF7.RData). 

Place the downloaded data in the RCL/code folder, in the same folder, run [this code file](https://github.com/Bin-Chen-Lab/GPS/blob/main/RCL/code/pre_process.R) for data pre-processing. This will generate two new files "unique_drug_profile.csv" and "unique_drug_value.csv" in the RCL/code folder.





