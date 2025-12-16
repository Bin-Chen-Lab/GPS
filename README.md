# GPS: an AI-platform for novel compound discovery using transcriptomics

![alt text](technical/GPS_GitHub.png)

## Background

Identifying drugs that reverse expression of disease-associated transcriptomic features has been widely explored as a strategy for discovering drug repurposing candidates, but its potential for novel compound discovery and optimization remains largely underexplored. We developed a deep learning-based platform that predicts gene expression changes from chemical structures, enabling high-throughput screening of large compound libraries. The platform refines compound scoring and employs a Monte Carlo Tree Search (MCTS) for multi-objective optimization and by incorporating Structure-Gene-Activity Relationships (SGAR), it can be used to uncover potential drug mechanisms directly from transcriptomic data, supporting both compound prioritization and lead optimization.

---

## About GPS and other Core Components

The **Gene expression profile Predictor on chemical Structures (GPS)** is an open-source platform under the Apache 2.0 Licence that allows researchers to predict the effects of chemical structures on gene expression, screen large-scale compound libraries, and optimize lead compounds for specific medicinal chemistry properties. Users can retrain models with their own data, and the platform is designed for community-driven improvements and adaptations.  

GPS consists of three core components:

- **RCL for training** – allows the user to retrain the model using custom drug/cell line data  
- **GPS4Drugs** – predict the effect of any structure on gene expression and screen drugs in custom or prescreened libraries  
- **MolSearch** – optimize chemical structures to improve multiple medicinal chemistry properties

---

## Repository


Clone the repository with:

```bash
git clone https://github.com/Bin-Chen-Lab/GPS/
```
Don't forget to download large files for each folder. See specific documentation in each folder.

---

## DockerHub Images

For ease of use, we provide DockerHub repositories:  

- [GPS4Drugs DockerHub](https://hub.docker.com/repository/docker/leshchi4/gpsimage/general)  
- [MolSearch DockerHub](https://hub.docker.com/repository/docker/leshchi4/molsearch)  

Detailed instructions on how to run Docker available inside respective folders.

---

## Online resources

Novel compound screening portal http://apps.octad.org/GPS/.

Drug repurposing web portal is available http://octad.org/ and the R package octad is available in Bioconductor. 

---

## GPS Documentation and Demo

Detailed documentation is available for each component within its respective folder.

A demo file is also available in the **`demo`** folders:

- for [GPS4Drugs](https://github.com/Bin-Chen-Lab/GPS/tree/main/GPS4Drugs/demo)
- for [MolSearch](https://github.com/Bin-Chen-Lab/GPS/tree/main/MolSearch/demo)

This file is the recommended starting point, and together with the documentation, can be used as input to run each pipeline starting from the first step.

---

## Figure generation code

The code in the [figure_code](https://github.com/Bin-Chen-Lab/GPS/tree/main/figure_code) folder is used to generate key figures in the paper

---

## Sequencing data

RNA-seq data are deposited in the GEO under accession number (GSE291867, GSE291190, and GSE291833). 

---

## Authors

Jing Xing, et. al., Deep learning-based screening and design of novel therapeutics that reverse disease-associated transcriptional phenotypes, submitted.

---

## Contact

If you have any questions please contact us at: contact@octad.org
