# Deep learning-based screening and design of novel therapeutics that reverse disease-associated transcriptional phenotype

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
# GPS Prediction Profile Docker Image

This Docker image runs the **GPS_runPredProfile.py** script for compound profiling and **GPS_runDrugScreenRges.py** script for drug screening using GPU acceleration.  
It comes with all dependencies pre-installed and organized project directories.

---

## Quick Start

To pull the image from DockerHub run:

```bash
sudo docker pull leshchi4/gpsimage:latest
```

### Input Files

- **GPS_runPredProfile.py** input file should be named:  
  ```
  cmpd__<sample_name>.csv
  ```
  It must contain columns:
  - **ID** → compound names  
  - **SMILES** → compound structure in SMILES format  

- **GPS_runDrugScreenRges.py** input file should be named:  
  ```
  DZSIG__<sample_name>.csv
  ```
  It must contain columns:
  - **GeneSymbol** → gene symbol ID  
  - **Value** → log2FC or other directional expression value  

---

## Running the Container

First, create input and output folders in your working directory:

```bash
mkdir -p input output
```

### GPS_runPredProfile.py (predicts expression profile)

```bash
sudo docker run --rm \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/output:/app/data/profile_pred/MEDIAN/preds_all \
    leshchi4/gpsimage:latest \
    python code/GPS_runPredProfile.py --cmpd_input input/cmpd__<sample_name>.csv
```

#### Flags for `GPS_runPredProfile.py`
- `--cmpd_input`, type=str, help='Input csv for compound ID and SMILES', default='../data/input_cmpd_gene/cmpd__TestJob0.csv'  
- `--gene_input`, type=str, help='Input list of gene symbols', default='preselect'  
- `--cpu_num`, type=int, help='Number of cpu cores to use', default=10  

---

### Using RGES Module in MolSearch

If using the **RGES module** in MolSearch, create an additional folder and mount it when running `GPS_runDrugScreenRges.py`:

```bash
-v $(pwd)/bgrd_pkl:/app/data/dzsig
```

The `BGRD__<sample_name>.pkl` file will be an input for `molsearch_rges.ipynb` along with `DZSIG__<sample_name>.csv`.  
You can find `molsearch_rges.ipynb` here:  
[MolSearch Notebook](https://github.com/Bin-Chen-Lab/GPS/tree/main/MolSearch)

---

### GPS_runDrugScreenRges.py (reversal prediction)

#### To use included library:

```bash
sudo docker run --rm \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/output:/app/data/reversal_score \
    leshchi4/gpsimage:latest \
    python code/GPS_runDrugScreenRges.py --dzSigFile input/DZSIG__<sample_name>.csv --cmpdLibID HTS
```

#### To use output of `GPS_runPredProfile.py`:

If you want to run a screen against your prediction, copy the file  
`<sample_name>_MEDIAN_GeneExpressionChange.csv` from the output folder into the input folder, then run:

```bash
sudo docker run --rm \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/output:/app/data/reversal_score \
    leshchi4/gpsimage:latest \
    python code/GPS_runDrugScreenRges.py --dzSigFile input/DZSIG__<sample_name>.csv --cmpdLibID input/<sample_name>_MEDIAN_GeneExpressionChange.csv
```

#### Flags for `GPS_runDrugScreenRges.py`
- `--dzSigFile`, type=str, help='Disease signature file', default=GATE + 'data/dzsig/DZSIG__TestJobNotExists.csv'  
- `--cmpdLibID`, type=str, help='Library to use', default='ZINC'  
- `--cpu_num`, type=int, help='Number of cpu cores to use', default=10  

---

## GPU Support

To use `GPS_runPredProfile.py` with GPU acceleration, install **nvidia-docker**:  
[Installation Guide](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)

You will also need a GPU that supports **CUDA 10.1**.

Run with GPU support:

```bash
sudo docker run --rm --gpus all --runtime=nvidia \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/output:/app/data/profile_pred/MEDIAN/preds_all \
    leshchi4/gpsimage:latest \
    python code/GPS_runPredProfile.py --cmpd_input input/cmpd__<sample_name>.csv
```

---

## Debugging

For debugging, replace your output mount with logs:

```bash
-v $(pwd)/output:/app/logs
```

Send the logfile to the developers to help diagnose issues.

The code is used to generate key figures in the following paper:
Jing Xing, et. al., Deep learning-based screening and design of novel therapeutics that reverse disease-associated transcriptional phenotypes, submitted.

Novel compound screening portal http://apps.octad.org/GPS/.

Drug repurposing web portal is available http://octad.org/ and the R package octad is available in Bioconductor. 

RNA-seq data are deposited in the GEO under accession number (GSE291867, GSE291190, and GSE291833). 
