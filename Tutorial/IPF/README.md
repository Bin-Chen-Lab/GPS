# IPF Drug Repurposing and Novel Compound Screening Tutorial

This tutorial provides a **simplified, end-to-end workflow** for transcriptomics-guided **drug repurposing** and **novel compound screening** in idiopathic pulmonary fibrosis (IPF). It is designed to help users quickly understand the pipeline and adapt it to **other diseases or biological conditions of interest**.

> **Note**  
> This tutorial presents a streamlined version of the key workflow described in the paper. Due to simplifications and stochastic components, results may vary from those reported in the manuscript.

---

## Overview of the Workflow

1. **Identify drug repurposing candidates** using OCTAD and LINCS data
2. **Construct a virtual screening library** from the ENAMINE HTS compounds
3. **Perform transcriptomics-guided virtual screening** using GPS

---

## Step 1: Identification of Drug Repurposing Candidates Using OCTAD

This step uses **OCTAD** to identify candidate drugs that can reverse a **cell-type–specific IPF transcriptional signature**, based on the LINCS perturbation library.

### Step 1.1. Install OCTAD

Install the drug repurposing pipeline OCTAD following the official repository instructions:

https://github.com/Bin-Chen-Lab/octad

### Step 1.2. Run IPF-specific repurposing analysis

Execute the script `OCTAD_IPF.R` to generate repurposing results.

### Expected output

- Runtime is typically **< 1 minute**
- The compound **pyrithyldione** should appear among the **top-ranked drugs** predicted to reverse the **MUC5B⁺ transcriptional signature**

---

## Step 2: Construction of a Virtual Screening Library from the ENAMINE HTS Collection

This step demonstrates how to construct a **virtual compound library** for downstream screening. In this example, we use active compounds from the **ENAMINE HTS library**, defined as compounds that induce a significant number of dysregulated genes precomputed by GPS.

> **Note**  
> The example compounds provided here can be replaced with compounds from any other chemical library of interest.

### Step 2.1. Set up GPS profile prediction

Refer to the following documentation to set up transcriptomic profile prediction using GPS:

https://github.com/Bin-Chen-Lab/GPS/tree/main/GPS4Drugs#gps-prediction-profile-docker-image

### Step 2.2. Prepare input files

- Place the file `cmpd__HTS.csv` into the `input/` folder
- Ensure Docker is properly configured

### Step 2.3. Run GPS profile prediction

From the working directory, execute the following command:

```bash
sudo docker run --rm \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/output:/app/data/profile_pred/MEDIAN/preds_all \
    leshchi4/gpsimage:latest \
    python code/GPS_runPredProfile.py --cmpd_input input/cmpd__HTS.csv
```

### Output

- The file `HTS_MEDIAN_GeneExpressionChange.csv` will be generated in the `output/` folder
- Runtime depends strongly on the size of the compound library

---

## Step 3: Transcriptomics-Guided Virtual Screening Using the Virtual Library

This step uses the newly generated **virtual compound library** to identify compounds that reverse the **MUC5B⁺ IPF signature**.

### Step 3.1. Prepare the disease signature

- Place the file `DZSIG__MUC5B+.csv` into the `input/` folder
- Ensure that column names strictly follow the required format

### Step 3.2. Run GPS-based virtual screening

From the working directory, execute the following command:

```bash
sudo docker run --rm \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/output:/app/data/reversal_score \
    leshchi4/gpsimage:latest \
    python code/GPS_runDrugScreenRges.py --dzSigFile input/DZSIG__MUC5B+.csv --cmpdLibID HTS
```

### Output

- The file `MUC5B+_RGES_norm.csv` will be generated in the `output/` folder
- The file contains compound identifiers and corresponding **Z-RGES scores** for prioritization

---

## Summary

This tutorial demonstrates how to:

- Identify IPF drug repurposing candidates using OCTAD
- Generate a virtual transcriptomic library from HTS compounds
- Perform GPS-based virtual screening to identify compounds that reverse disease-specific signatures


