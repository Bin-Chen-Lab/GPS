# GPS Prediction Profile Docker Image

This Docker image runs the **GPS_runPredProfile.py** script for compound profiling and *GPS_runDrugScreenRges.py** script for drug screening using GPU acceleration.
It comes with all dependencies pre-installed and organized project directories.

---

## Quick Start

To pull the image from DockerHub run command:

sudo docker pull leshchi4/gpsimage:latest

Input file for GPS_runPredProfile.py should be named cmpd__<sample_name>.csv and contain columns ID, SMILES
ID column should contain the compound names and SMILES column the structure of compound in SMILES format

Input file for GPS_runDrugScreenRges.py should be named DZSIG__<sample_name>.csv contain columns GeneSymbol, Value
GeneSymbol column should contain the GeneSymbol ID and Value column the log2FC value or any other directional expression value for that gene


Run the container with your input directory mounted:

To create two folders in your working directory run command:
mkdir -p input output

Commands for GPS_runPredProfile.py (predicts the expression profile):

sudo docker run --rm \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/output:/app/data/profile_pred/MEDIAN/preds_all \
    leshchi4/gpsimage:latest \
    python code/sudo docker run --rm \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/output:/app/data/reversal_score \
    leshchi4/gpsimage:latest \
    python code/GPS_runDrugScreenRges.py --dzSigFile input/DZSIG__<sample_name>.csv --cmpdLibID HTS
 --cmpd_input input/cmpd__<sample_name>.csv

Flags for GPS_runPredProfile.py:

--cmpd_input, type=str, help='Input csv for compound ID and SMILES', default='../data/input_cmpd_gene/cmpd__TestJob0.csv'
--gene_input, type=str, help='Input list of gene symbols', default='preselect'
--cpu_num, type=int, help='Number of cpu cores to use', default=10


If you will be using the RGES module in the MolSearch pipeline we advise to also create another folder and mount this output folder when running GPS_runDrugScreenRges.py:

-v $(pwd)/bgrd_pkl:/app/data/dzsig

The BGRD__<sample_name>.pkl file will be one of the input for molsearch_rges.ipynb script along side the DZSIG__<sample_name>.csv.
molsearch_rges.ipynb can be found in the GitHub at https://github.com/Bin-Chen-Lab/GPS/tree/main/MolSearch

Commands for GPS_runDrugScreenRges.py (reversal prediction against the disease signature):

To use included library:

sudo docker run --rm \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/output:/app/data/reversal_score \
    leshchi4/gpsimage:latest \
    python code/GPS_runDrugScreenRges.py --dzSigFile input/DZSIG__<sample_name>.csv --cmpdLibID HTS

To use output of GPS_runPredProfile.py:
If you want to run a screen against your prediction copy the <sample_name>_MEDIAN_GeneExpressionChange.csv from the output folder to the input folder.

sudo docker run --rm \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/output:/app/data/reversal_score \
    leshchi4/gpsimage:latest \
    python code/GPS_runDrugScreenRges.py --dzSigFile input/DZSIG__<sample_name>.csv --cmpdLibID input/<sample_name>_MEDIAN_GeneExpressionChange.csv


Flags for GPS_runDrugScreenRges.py:

--dzSigFile, type = str, help = 'Disease signature file', default=GATE + 'data/dzsig/DZSIG__TestJobNotExists.csv'
--cmpdLibID, type = str, help = 'Library to use', default='ZINC'
--cpu_num, type = int, help = 'Number of cpu cores to use', default=10

## GPU support (will run faster)

To use GPS_runPredProfile.py with GPU users will need to install nvidia-docker:

https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html

Users will also need a GPU that supports CUDA 10.1.

Run command with following flags:

sudo docker run --rm --gpus all --runtime=nvidia \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/output:/app/data/profile_pred/MEDIAN/preds_all \
    leshchi4/gpsimage:latest \
    python code/GPS_runPredProfile.py --cmpd_input input/cmpd__<sample_name>.csv

For debugging replace your output mounting with app/logs mounting:

-v $(pwd)/output:/app/logs \

Send the logfile to us to help us fix any issues.



