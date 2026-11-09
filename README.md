# SMA-Numerical-Analyses
This repository contains the code, data, and results for all simulation studies and empirical analyses in the paper “Model Averaging Approach for Analyzing Complex Survey Data”.

## Organization
### code
This directory contains the R code for all simulation studies and empirical analyses in the paper. Each script performs 1,000 replications. The results can be reproduced by adjusting the settings as described below.

Sim1.R reproduces Tables 1–2 in the main manuscript and Tables S1–S6 in the Supplementary Material. Setting the number of candidate models to $M=[n/3]$ reproduces Tables S7–S14. Using the true variance $\sigma^2=1$ instead of its estimate reproduces Tables S15–S22.

Sim2.R reproduces Tables S23–S24 in the Supplementary Material.

Sim3.R reproduces Figure S1 in the Supplementary Material.

real1.R reproduces the results for the first empirical application in the main manuscript.

real2.R reproduces the results for the second empirical application in the main manuscript.
### data
This directory contains pm2.5.csv, the dataset used in the first empirical application in the main manuscript. The dataset used in the second empirical application, ethylene_methane.txt, can be downloaded from the website: https://archive.ics.uci.edu/dataset/322.
### manuscript
The results produced by the code are reported in Tables 1–6 in Section 5 of the main manuscript and in Tables S1–S24 and Figure S1 in Section S5 of the Supplementary Material.
