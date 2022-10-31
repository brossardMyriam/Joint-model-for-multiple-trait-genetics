# Joint-model-for-multiple-trait-genetics

This Github repository includes Example of R script for data simulation under multiple trait genetics model & to fit a joint model for multiple longitudinal and multiple time-to-event traits using  the two-stage approach with empirical estimation of the variance-covariance matrix by bootstrap, see details in the paper (under revision):
"Characterization of direct and/or indirect genetic associations for multiple traits in longitudinal studies of disease progression"
by Myriam Brossard, Andrew D. Paterson, Osvaldo Espin-Garcia, Radu V. Craiu, Shelley B. Bull

./simulation includes:
- the R script used to simulate replicates of datasets with  K=2 time-to-event traits (retinopathy, nephropathy), based on L=3 longitudinal risk factors, 2 observed in DCCT (hba1c, SBP) and one simulated (U) and 5 causal SNPs with direct and/or indirect effects on these time-to-event traits (see Fig. 3 from Brossard et al). The script also generates genotypes for 5 SNPs under the global null hypothesis (referred as SNP1R-SNP5R with same MAFS as the causal SNPs)
- an ariticial DCCT dataset with simulated longitudunal Hba1c and SBP values & 2 baseline covariates (SEX, T1D_diagnosis) for N=xxx DCCT individuals to illustrate the procedure (input dataset_
- an example of dataset produced by the script for 5 replicates (based on the same simulation parameters as used in Brossard et al,)

./JM_analysis includes:
- a main R script that fit a joint model of L=2 longitudinal QT (HbA1c, SBP) and K=2 time-to-event traits (DR, DN) for 1 simulation replicate & SNP5A
- a R script that includes some sub-functions to fit the joint model & the R script that mvlme.r
- an example of results produced by the ZZZZZ

Commants in the scripts help to interpret the different outputs
