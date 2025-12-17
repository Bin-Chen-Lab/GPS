# HCC Novel Compound Discovery Tutorial

This tutorial provides a **simplified, end‑to‑end workflow** for transcriptomics‑guided novel compound discovery in hepatocellular carcinoma (HCC). It is intended to help users quickly understand the pipeline and adapt it to **other diseases or conditions of interest**.

> **Note**  
> This is a streamlined version of the full workflow described in the paper. Because of simplifications and stochastic components, results may differ slightly from those reported in the manuscript.

---

## Overview of the Workflow

1. **Predict GPS transcriptomic profiles** for novel compounds from the ZINC library
2. **Stage 1 optimization**: refine initial hits based on RGES and explore basic structure and RGES relationships.
3. **Stage 2 optimization**: further optimize compounds for drug‑like and ADME‑related properties

---

## GPS Profile Prediction on Novel Compounds from the ZINC Library

### Step 1. Set up the working directory

Create a working directory and initialize the following subfolders:

- `input/`
- `output/`
- `library/`

### Step 2. Prepare the disease signature

Place the file `DZSIG__HCC.csv` into the `input/` folder.

- `DZSIG__HCC` refers to the HCC disease signature
- This file can be replaced with signatures for other diseases
- Please ensure that **column names strictly follow the required format**

### Step 3. Set up the GPS4Drugs Docker environment

Follow the official instructions to set up the GPS4Drugs Docker image:

https://github.com/Bin-Chen-Lab/GPS/tree/main/GPS4Drugs

### Step 4. Download the ZINC compound library

Download the preprocessed ZINC library from the link below and place it into the `library/` folder:

https://chenlab-data-public.s3.amazonaws.com/ZINC_strong.npz

### Step 5. Run GPS‑based virtual screening

From the working directory, execute the following command:

```bash
sudo docker run --rm \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/library:/app/data/profile_pred/MEDIAN \
    -v $(pwd)/output:/app/data/reversal_score \
    leshchi4/gpsimage:latest \
    python code/GPS_runDrugScreenRges.py --dzSigFile input/DZSIG__<sample_name>.csv --cmpdLibID ZINC
```

### Output

- The output file `HCC_RGES_norm.csv` will be generated in the `output/` folder
- The file contains **ZINC compound IDs** and corresponding **Z‑RGES scores**
- Runtime is typically **~1 hour**, depending on hardware

> **Expected result**  
> Due to randomization, results may vary slightly. However, the compound **ZINC000000086363** should appear within the **top 100 ranked hits**.

---

## Stage 1: Optimization of Compounds

### Objective

Stage 1 focuses on optimizing the initial hit compound (**ZINC000000086363**) while maintaining comparable RGES values and exploring basic structure and RGES patterns.

### Step 1. Prepare input files

Create the following folders in the working directory:

- `input/`
- `output/`
- `rges_input/`

Then:

- Place the file `taskhcc` (containing the SMILES structure of the hit compound) into `input/`
- Place all required `*.csv` files into `rges_input/`

> **Note**  
> To apply this workflow to other disease signatures, you must create a corresponding **RGES module** following the MolSearch instructions:
>
> https://github.com/Bin-Chen-Lab/GPS/tree/main/MolSearch
>
> The current module is configured specifically for HCC drug discovery.

### Step 2. Run Stage 1 optimization

```bash
sudo docker run --rm \
    -v $(pwd)/input:/app/MCTS/libs/start_mols \
    -v $(pwd)/output:/app/MCTS/results_visulization \
    -v $(pwd)/rges_input:/app/MCTS/libs/rges_input \
    leshchi4/molsearch:latest \
    python MCTS/molsearch1_auto.py --num_drugs 1 --sample_name taskhcc --pool_cores 1 --goals rges --sig_name HCC
```

### Output

- The output folder `rges_stage1` will be generated under `taskhcc/`

---

## Stage 2: Multi‑objective Optimization of Compounds

### Objective

Stage 2 further optimizes compounds generated in Stage 1 by incorporating additional objectives, such as:

- LogP
- Drug‑likeness (QED)
- Synthetic accessibility
- RGES

### Run Stage 2 optimization

```bash
sudo docker run --rm \
    -v $(pwd)/output:/app/MCTS/results_visulization \
    -v $(pwd)/rges_input:/app/MCTS/libs/rges_input \
    leshchi4/molsearch:latest \
    python MCTS/molsearch2_auto.py --num_drugs 1 --sample_name taskhcc --previous_goals rges --pool_cores 1 --goals plogp_qed_sa_rges --sig_name HCC
```

### Output and interpretation

- The output folder `plogp_qed_sa_rges_stage2` will be created under `taskhcc/`
- This folder contains multiple files corresponding to different optimized compounds
- Files can be concatenated to generate a final candidate list

You may visualize and further evaluate compound structures using external tools such as:

https://www.swissadme.ch/

> **Observation**  
> You should frequently observe substitutions near the central nitrogen‑containing functional group. Runtime can vary substantially depending on starting compounds and parameter settings, and may take **several hours**.

---

## Extending the Workflow

To optimize additional ADME or physicochemical properties, corresponding modules must be implemented and added to the **MolSearch** framework.

---

## Summary

This tutorial demonstrates how to:

- Perform transcriptomics‑guided virtual screening using GPS
- Optimize hit compounds through multi‑stage molecular generation
- Extend the framework for other diseases and optimization objectives


