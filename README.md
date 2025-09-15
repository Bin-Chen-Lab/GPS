# Deep learning-based screening and design of novel therapeutics that reverse disease-associated transcriptional phenotype

![alt text](technical/GPS_GitHub.png)

Identifying drugs that reverse expression of disease-associated transcriptomic features has been widely explored as a strategy for discovering drug repurposing candidates, but its potential for novel compound discovery and optimization remains largely underexplored. 

Here, we present a **deep learning-based drug discovery platform**, guided by transcriptomic features, that screens large compound libraries and optimizes lead compounds. We first develop a model that predicts gene expression changes solely from chemical structures and deploy it to infer the expression changes of compounds in large screening libraries. We then refine compound scoring and employ a Monte Carlo Tree Search method for multi-objective optimization. By incorporating Structure-Gene-Activity Relationships, we uncover drug mechanisms directly from transcriptomic data.

---

## About GPS

The **Gene expression profile Predictor on chemical Structures (GPS)** is a set of software algorithms that allow for:

- Prediction of effects of chemical structures on gene expression  
- Screening of large-scale chemical libraries  
- Optimization of lead compounds for specific medicinal chemistry properties  

We published this code on GitHub to enable researchers to perform all of these tasks, as well as retrain the model using their own data. This code is open source to encourage the community to improve, adapt, and tailor this approach for their own needs, as well as to address wider scientific questions.

---

## Core Components

The GPS core software consists of 3 major components:

- **RCL for training** – retrain the model using your drug/cell line data  
- **GPS4Drugs** – predict the effect of any structure on gene expression and screen drugs in custom or prescreened libraries  
- **MolSearch** – optimize chemical structures to improve multiple medicinal chemistry properties  

---

## Repository

Clone the repository with:

```bash
git clone https://github.com/Bin-Chen-Lab/GPS/
```

Detailed documentation is available for each component within its respective folder.

---

## DockerHub Images

For ease of use, we provide DockerHub repositories:  

- [GPS4Drugs DockerHub](https://hub.docker.com/repository/docker/leshchi4/gpsimage/general)  
- [MolSearch DockerHub](https://hub.docker.com/repository/docker/leshchi4/molsearch)  

---

## GPS Documentation and Demo

The main documentation for each component is provided inside its respective folder.  

A demo file is also available in the **`demo`** folder. This file is the recommended starting point, and together with the documentation, can be used as input to run each pipeline starting from the first step.

The code is used to generate key figures in the following paper:
Jing Xing, et. al., Deep learning-based screening and design of novel therapeutics that reverse disease-associated transcriptional phenotypes, submitted.

Novel compound screening portal http://apps.octad.org/GPS/.

Drug repurposing web portal is available http://octad.org/ and the R package octad is available in Bioconductor. 

RNA-seq data are deposited in the GEO under accession number (GSE291867, GSE291190, and GSE291833). 
