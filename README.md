# Deep learning-based screening and design of novel therapeutics that reverse disease-associated transcriptional phenotype<img width="468" height="61" alt="image" src="https://github.com/user-attachments/assets/11ec2749-9dd2-4a9f-a120-b7642b099110" />

## Table of Contents
- [Environment Setup](#environment-setup)
- [GPS Prediction](#gps-prediction)
  - [Input Format](#input-format)
  - [Command](#command)
  - [Flags](#flags)
- [GPS Drug Screening (Reversal Prediction)](#gps-drug-screening-reversal-prediction)
  - [Input Format](#input-format-1)
  - [Commands](#commands)
  - [Flags](#flags-1)
- [Molecular Search](#molecular-search)
  - [Input Format](#input-format-2)
  - [Stage 1: molsearch1_auto.py](#stage-1-molsearch1_autopy)
    - [Flags](#flags-2)
  - [Stage 2: molsearch2_auto.py](#stage-2-molsearch2_autopy)
    - [Flags](#flags-3)

---

## Environment Setup

To run both pipelines, install the Python 2.7 environment:

```bash
conda env create --file py27.yml
```

The `py27.yml` file can be found on GitHub in the **GPS4Drugs** folder.

---

## GPS Prediction

### Input Format

Input file for `GPS_runPredProfile.py`:

- Filename: `cmpd__<sample_name>.csv`  
- Required columns:  

| Column | Description |
|--------|-------------|
| ID     | Compound name |
| SMILES | Structure of compound in SMILES format |

---

### Command

```bash
python GPS_runPredProfile.py --cmpd_input cmpd__<sample_name>.csv
```

### Flags

| Flag          | Type | Description                           | Default                                               |
|---------------|------|---------------------------------------|-------------------------------------------------------|
| `--cmpd_input` | str  | Input csv for compound ID and SMILES | `../data/input_cmpd_gene/cmpd__TestJob0.csv` |
| `--gene_input` | str  | Input list of gene symbols           | `preselect`                                           |
| `--cpu_num`    | int  | Number of cpu cores to use           | `10`                                                  |

---

## GPS Drug Screening (Reversal Prediction)

### Input Format

Input file for `GPS_runDrugScreenRges.py`:

- Filename: `DZSIG__<sample_name>.csv`  
- Required columns:  

| Column     | Description |
|------------|-------------|
| GeneSymbol | Gene symbol ID |
| Value      | log2FC value or other directional expression value |

---

### Commands

Using included library:

```bash
python GPS_runDrugScreenRges.py --dzSigFile DZSIG__<sample_name>.csv --cmpdLibID HTS
```

Using output of `GPS_runPredProfile.py`:

```bash
python GPS_runDrugScreenRges.py --dzSigFile DZSIG__<sample_name>.csv --cmpdLibID input/<sample_name>_MEDIAN_GeneExpressionChange.csv
```

### Flags

| Flag         | Type | Description              | Default                                     |
|--------------|------|--------------------------|---------------------------------------------|
| `--dzSigFile` | str  | Disease signature file   | `data/dzsig/DZSIG__TestJobNotExists.csv`    |
| `--cmpdLibID` | str  | Library to use           | `ZINC`                                      |
| `--cpu_num`   | int  | Number of cpu cores      | `10`                                        |

---

## Molecular Search

### Input Format

Input file for `molsearch1_auto.py`:

- Filename: `<sample_name>.csv`  
- Required columns:  

| Column  | Description |
|---------|-------------|
| idx     | Numeric index for compound (starting at 0) |
| smiles  | Structure of compound in SMILES format |
| molname | Compound name (string) |

> ⚠️ If using the **RGES module**, first run `GPS_runDrugScreenRges.py` and prepare the output using `molsearch_rges.ipynb` from [MolSearch](https://github.com/Bin-Chen-Lab/GPS/tree/main/MolSearch).

---

### Stage 1: `molsearch1_auto.py`

```bash
python molsearch1_auto.py --num_drugs 1 --sample_name <sample_name> --pool_cores 1 --goals bbbp_rges
```

#### Flags

| Flag          | Type | Description                                      | Default    |
|---------------|------|--------------------------------------------------|------------|
| `--num_drugs`  | int  | Number of drugs in CSV file                      | `1`        |
| `--sample_name`| str  | CSV filename (no `.csv`)                         | `taskdipg` |
| `--pool_cores` | int  | Number of cores to use                           | `1`        |
| `--goals`      | str  | Objectives: `plogp`, `qed`, `sa`, `rges`, `bbbp` (use `_` to connect) | `bbbp_rges` |
| `--sig_name`   | int  | Signature name                                   | `None`     |

---

### Stage 2: `molsearch2_auto.py`

Input file: Output of Stage 1 (stored in `output/`).

```bash
python molsearch2_auto.py --num_drugs 1 --sample_name <sample_name> --previous_goals bbbp_rges --pool_cores 1 --goals plogp_qed_sa_rges
```

#### Flags

| Flag            | Type | Description              | Default             |
|-----------------|------|--------------------------|---------------------|
| `--num_drugs`    | int  | Number of drugs in CSV   | `1`                 |
| `--sample_name`  | str  | CSV filename (no `.csv`) | `mmfdipg`           |
| `--previous_goals` | str | Stage 1 goals            | `bbbp_rges`         |
| `--pool_cores`   | int  | Number of cores to use   | `1`                 |
| `--goals`        | str  | Objectives: `plogp`, `qed`, `sa`, `rges`, `bbbp` (use `_` to connect) | `plogp_qed_sa_rges` |
| `--sig_name`     | int  | Signature name           | `None`              |


The code is used to generate key figures in the following paper:
Jing Xing, et. al., Deep learning-based screening and design of novel therapeutics that reverse disease-associated transcriptional phenotypes, submitted.

Novel compound screening portal http://apps.octad.org/GPS/.

Drug repurposing web portal is available http://octad.org/ and the R package octad is available in Bioconductor. 

RNA-seq data are deposited in the GEO under accession number (GSE291867, GSE291190, and GSE291833). 
