# Molecular Search Docker Image

This Docker image runs Molecular Search for lead optimization using GPU acceleration.  
It comes with all dependencies pre-installed and organized project directories.

---

## Quick Start

To install nvidia-docker follow up:

https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html

To pull the image from DockerHub run command:

sudo docker pull leshchi4/molsearch:latest

Input file for molsearch1_auto.py should be named <sample_name>.csv and contain columns idx,smiles,molname
isx column should contain the numeric index for compound starting with 0
smiles column the structure of compound in SMILES format
molname should contain the name of compound in string format

Input file for molsearch2_auto.py will be the output of stage one which will be located in the output folder.

If using the RGES module you must first run GPS_runDrugScreenRges.py script with appropriate directory mounted and prepare output using molsearch_rges.ipynb which can be found in the GitHub at https://github.com/Bin-Chen-Lab/GPS/tree/main/MolSearch

After completing that step place the output of molsearch_rges.ipynb into separate folder called rges_input and add this command to your MolSearch Docker run:
 -v $(pwd)/rges_input:/app/MCTS/score_modules/RGES_Score

Run the container with your directories mounted:

To create folders in your working directory run command:
mkdir -p input output rges_input

Run the container with your input file mounted:

Commands for molsearch1_auto.py (Stage 1 of optimization):

sudo docker run --rm --gpus all \
    -v $(pwd)/input:/app/MCTS/libs/start_mols \
    -v $(pwd)/output:/app/MCTS/results_visulization \
    leshchi4/molsearch:latest \
    python MCTS/molsearch1_auto.py --num_drugs 1 --sample_name <sample_name> --pool_cores 1 --goals bbbp_rges

Flags for molsearch1_auto.py:

--num_drugs, type=int, help='number of drugs in csv file', default=1
--sample_name, type=str, help='csv filename no .csv', default='taskdipg'
--pool_cores, type=int, help='number of cores to use', default=1
--goals, type=str, help='plogp,qed,sa,rges,bbbp use _ to connect', default='bbbp_rges'
--sig_name, type=int, help='number of cores to use', default='None'


Commands for molsearch2_auto.py (Stage 2 of optimization):

sudo docker run --rm --gpus all \
    -v $(pwd)/output:/app/MCTS/results_visulization \
    leshchi4/molsearch:latest \
    python MCTS/molsearch2_auto.py --num_drugs 1 --sample_name <sample_name> --previous_goals bbbp_rges --pool_cores 1 --goals plogp_qed_sa_rges

Flags for molsearch2_auto.py:

--num_drugs, type=int, help='number of drugs in csv file', default=1
--sample_name, type=str, help='csv filename no .csv', default='mmfdipg'
--previous_goals, type=str, help='stage 1 goals', default='bbbp_rges'
--pool_cores, type=int, help='number of cores to use', default=1
--goals, type=str, help='plogp,qed,sa,rges,bbbp use _ to connect', default='plogp_qed_sa_rges'
--sig_name, type=int, help='number of cores to use', default='None'




