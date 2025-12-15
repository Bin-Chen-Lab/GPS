# IPF Drug discovery tutorial

## OCTAD repurposing candidate prediction:

Install drug repurposing pipeline [OCTAD](https://github.com/Bin-Chen-Lab/octad).

Run OCTAD_IPF.R to get result.

## GPS profile prediction on novel active compounds from the ENAMINE HTS library:

To run prediction refer to the [GPS_runPredProfile.py (predicts expression profile)]
(https://github.com/Bin-Chen-Lab/GPS/tree/main/GPS4Drugs#gps-prediction-profile-docker-image)

After setting up Docker place the cmpd__HTS.csv into the input folder and run the following command from the pwd:

```bash
sudo docker run --rm \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/output:/app/data/profile_pred/MEDIAN/preds_all \
    leshchi4/gpsimage:latest \
    python code/GPS_runPredProfile.py --cmpd_input input/cmpd__HTS.csv
```

The file HTS_MEDIAN_GeneExpressionChange.csv should be located in the output folder.

---

## Run GPS profile prediction on active compounds from the ENAMINE HTS library:

To run the prediction place the DZSIG__MUC5B+.csv into input folder and run the following command from the pwd:

```bash
sudo docker run --rm \
    -v $(pwd)/input:/app/input \
    -v $(pwd)/output:/app/data/reversal_score \
    leshchi4/gpsimage:latest \
    python code/GPS_runDrugScreenRges.py --dzSigFile input/DZSIG__MUC5B+.csv --cmpdLibID HTS
```

The file MUC5B+_RGES_norm.csv should be located in the output folder.

---