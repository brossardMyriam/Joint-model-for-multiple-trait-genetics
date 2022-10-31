# Joint-model-for-multiple-trait-genetics
From paper (under revision) :  "Characterization of direct and/or indirect genetic associations for multiple traits in longitudinal studies of disease progression"
by Myriam Brossard, Andrew D. Paterson, Osvaldo Espin-Garcia, Radu V. Craiu, Shelley B. Bull

Example of R scripts to simulate datasets genotypes and time-to-event traits based on DCCT dataset & to fit a joint model for multiple longitudinal and multiple time-to-event traits using  the two-stage approach with empirical estimation of the variance-covariance matrix by bootstrap (see details in Brossard et al, manucript under revision)

This directories contains:
- 3 R scripts:
  1. for data-based simulations 
- Example of R script and for application of a joint model of two longitudinal and two time-to-event traits, based on two-stage fitting, with stage 1 based on fitting of a multivariate mixed model, and stage 2 based on Cox PH frailty time-to-event models using fitted trajectories of longitudinal ttaits
- an artificial dataset with longitudinal HbA1c and SBP values  & baseline covariates (SEX, T1D duration_years) simulated in N=xx individuals. HbA1c and SBP values have been simulated based on the original DCCT dataset. The purpose of this dataset is to illustrate the DCCT-based simulation script presented in the paper.
- a dataset with 5 replicates of DCCT-based simulated datasets, including genotypes data (at 5 causal SNPs: SNP1A, SNP2A, SNP3, SNP4A, SNP5A & 5 Null SNPs with same MAFs as the causal SNPs but generated under the global null hypothesis, and named: SNP1R-SNP5R) & K=2 simulated time-to-event traits
- an example of output produced by the script 

