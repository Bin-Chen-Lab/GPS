# Molecular Search

# Table of Contents

- [Molecular Search](#molecular-search)
  - [Stage 1: molsearch1_autopy](#stage-1-molsearch1_autopy)
  - [Stage 2: molsearch2_autopy](#stage-2-molsearch2_autopy)

- [Molecular Search Docker Image](#molecular-search-docker-image)
  - [MolSearch with Docker](#molsearch-with-docker)
    - [Installation](#installation)
    - [Input Files](#input-files)
    - [RGES Module (Optional)](#rges-module-optional)
    - [Running the Container](#running-the-container)

### MCTS Environment Setup

```bash
conda create --name mcts python=3.7 pip 
conda install -c conda-forge rdkit
pip install pandas seaborn pickle-mixin
pip install -U scikit-learn==0.21.3 (RF scorer of GSK3B and JNK3 requires this version of sklearn)
```
### In medchem_moves folder

```bash
python setup.py install
```

### Input Format

Input file for `molsearch1_auto.py`:

- Filename: `<sample_name>.csv`  
- Required columns:  

| Column  | Description |
|---------|-------------|
| idx     | Numeric index for compound (starting at 0) |
| smiles  | Structure of compound in SMILES format |
| molname | Compound name (string) |

Inputs should be located in MCTS/libs/start_mols.

> ⚠️ If using the **RGES module**, first run `GPS_runDrugScreenRges.py` and prepare the output using `molsearch_rges.ipynb` from [MolSearch](https://github.com/Bin-Chen-Lab/GPS/tree/main/MolSearch). Further details are in [GPS4Drugs](https://github.com/Bin-Chen-Lab/GPS/tree/main/GPS4Drugs#gps-drug-screening-reversal-prediction) under 'Using RGES Module in MolSearch'.


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
| `--sig_name`   | int  | Signature name for RGES module                   | `None`     |

---

### Stage 2: `molsearch2_auto.py`

Input file: Output of Stage 1 (stored in MCTS/results_visulization).

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
| `--sig_name`     | int  | Signature name for RGES module       | `None`              |

Output file: Output of Stage 2 is stored in MCTS/results_visulization.


# Molecular Search Docker Image

This Docker image runs Molecular Search for lead optimization using GPU acceleration.  
It comes with all dependencies pre-installed and organized project directories.

---
# MolSearch with Docker

## Installation

To install **NVIDIA Container Toolkit**, follow the official guide:  
[Installation Guide](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)

To pull the MolSearch Docker image from DockerHub:

```bash
sudo docker pull leshchi4/molsearch:latest
```

---

## Input Files

### Stage 1 Input (`molsearch1_auto.py`)
The input file should be named:

```
<sample_name>.csv
```

It must contain the following columns:

- **idx** → numeric index for compounds (starting with 0)  
- **smiles** → compound structure in SMILES format  
- **molname** → compound name (string)  

---

### Stage 2 Input (`molsearch2_auto.py`)
The input file will be the **output of Stage 1**, located in the `output` folder.

---

### RGES Module (Optional)

To run demo of the data preparation add contents of 'demo' folder to 'rges_input' folder.

If using the **RGES module**, first run `GPS_runDrugScreenRges.py` with the appropriate directory mounted.  
Then, prepare output using `molsearch_rges.ipynb`, which can be found in GitHub:  
[MolSearch RGES Notebook](https://github.com/Bin-Chen-Lab/GPS/tree/main/MolSearch)

After completing that step, the output of `molsearch_rges.ipynb` should be located in the `rges_input`, as well as data from the GPS run and add this command to your MolSearch Docker run:

```
-v $(pwd)/rges_input:/app/MCTS/libs/rges_input
```

---

## Running the Container

First, create the required folders:

```bash
mkdir -p input output
```

---

### Stage 1: Optimization (`molsearch1_auto.py`)

Run the container with your input file mounted:

```bash
sudo docker run --rm --gpus all \
    -v $(pwd)/input:/app/MCTS/libs/start_mols \
    -v $(pwd)/output:/app/MCTS/results_visulization \
    leshchi4/molsearch:latest \
    python MCTS/molsearch1_auto.py --num_drugs 1 --sample_name <sample_name> --pool_cores 1 --goals bbbp_rges
```
---

### Stage 2: Optimization (`molsearch2_auto.py`)

Run the container with your output file mounted:

```bash
sudo docker run --rm --gpus all \
    -v $(pwd)/output:/app/MCTS/results_visulization \
    leshchi4/molsearch:latest \
    python MCTS/molsearch2_auto.py --num_drugs 1 --sample_name <sample_name> --previous_goals bbbp_rges --pool_cores 1 --goals plogp_qed_sa_rges
```
