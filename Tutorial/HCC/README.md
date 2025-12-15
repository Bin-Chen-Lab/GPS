# HCC Drug discovery tutorial

## GPS profile prediction on novel compounds from the ZINC library:

To run the prediction place the DZSIG__HCC.csv into input folder and run the following command from the working directory:

After setting up Docker download the ZINC preprocessed library follow instructions in 'To use ZINC preprocessed library' in the [GitHub](https://github.com/Bin-Chen-Lab/GPS/tree/main/GPS4Drugs#gps-prediction-profile-docker-image) page and run the following command from the pwd:

'''bash
sudo docker run --rm \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/library:/app/data/profile_pred/MEDIAN \
    -v $(pwd)/output:/app/data/reversal_score \
    leshchi4/gpsimage:latest \
    python code/GPS_runDrugScreenRges.py --dzSigFile input/DZSIG__<sample_name>.csv --cmpdLibID ZINC
'''

The file HCC_RGES_norm.csv should be located in the output folder.

## Stage 1 optimization of compounds 

Create input, output and rges_input folder in working directory. Place taskhcc into input and the rest of files into rges_input.

Run following command:

'''bash
sudo docker run --rm --gpus all \
    -v $(pwd)/input:/app/MCTS/libs/start_mols \
    -v $(pwd)/output:/app/MCTS/results_visulization \
    -v $(pwd)/rges_input:/app/MCTS/libs/rges_input \
    leshchi4/molsearch:latest \
    python MCTS/molsearch1_auto.py --num_drugs 1 --sample_name taskhcc --pool_cores 1 --goals rges --sig_name HCC
'''

## Stage 2 optimization of compounds 


Run following command:

'''bash
sudo docker run --rm --gpus all \
    -v $(pwd)/output:/app/MCTS/results_visulization \
    -v $(pwd)/rges_input:/app/MCTS/libs/rges_input \
    leshchi4/molsearch:latest \
    python MCTS/molsearch2_auto.py --num_drugs 1 --sample_name taskhcc --previous_goals rges --pool_cores 1 --goals plogp_qed_sa_rges --sig_name HCC
'''
