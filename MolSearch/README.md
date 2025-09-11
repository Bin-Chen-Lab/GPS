# Molecular Search

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
If using the **RGES module**, first run `GPS_runDrugScreenRges.py` with the appropriate directory mounted.  
Then, prepare output using `molsearch_rges.ipynb`, which can be found in GitHub:  
[MolSearch RGES Notebook](https://github.com/Bin-Chen-Lab/GPS/tree/main/MolSearch)

After completing that step, place the output of `molsearch_rges.ipynb` into a separate folder called `rges_input` and add this command to your MolSearch Docker run:

```
-v $(pwd)/rges_input:/app/MCTS/score_modules/RGES_Score
```

---

## Running the Container

First, create the required folders:

```bash
mkdir -p input output rges_input
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

#### Flags for `molsearch1_auto.py`
- `--num_drugs`, type=int, help='number of drugs in csv file', default=1  
- `--sample_name`, type=str, help='csv filename no .csv', default='taskdipg'  
- `--pool_cores`, type=int, help='number of cores to use', default=1  
- `--goals`, type=str, help='plogp,qed,sa,rges,bbbp use _ to connect', default='bbbp_rges'  
- `--sig_name`, type=int, help='number of cores to use', default='None'  

---

### Stage 2: Optimization (`molsearch2_auto.py`)

Run the container with your output file mounted:

```bash
sudo docker run --rm --gpus all \
    -v $(pwd)/output:/app/MCTS/results_visulization \
    leshchi4/molsearch:latest \
    python MCTS/molsearch2_auto.py --num_drugs 1 --sample_name <sample_name> --previous_goals bbbp_rges --pool_cores 1 --goals plogp_qed_sa_rges
```

#### Flags for `molsearch2_auto.py`
- `--num_drugs`, type=int, help='number of drugs in csv file', default=1  
- `--sample_name`, type=str, help='csv filename no .csv', default='mmfdipg'  
- `--previous_goals`, type=str, help='stage 1 goals', default='bbbp_rges'  
- `--pool_cores`, type=int, help='number of cores to use', default=1  
- `--goals`, type=str, help='plogp,qed,sa,rges,bbbp use _ to connect', default='plogp_qed_sa_rges'  
- `--sig_name`, type=int, help='number of cores to use', default='None'  

The code is used to generate key figures in the following paper:
Jing Xing, et. al., Deep learning-based screening and design of novel therapeutics that reverse disease-associated transcriptional phenotypes, submitted.

Novel compound screening portal http://apps.octad.org/GPS/.

Drug repurposing web portal is available http://octad.org/ and the R package octad is available in Bioconductor. 

RNA-seq data are deposited in the GEO under accession number (GSE291867, GSE291190, and GSE291833). 
