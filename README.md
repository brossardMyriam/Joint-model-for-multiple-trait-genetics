# Joint-model-for-multiple-trait-genetics

This repository includes some Examples of R scripts to generate the DCCT-based simulated datasets under a multiple-trait genetics scenario, as described in the paper (under revision): "Characterization of direct and/or indirect genetic associations for multiple traits in longitudinal studies of disease progression" by Myriam Brossard, Andrew D. Paterson, Osvaldo Espin-Garcia, Radu V. Craiu, Shelley B. Bull, & to fit a Joint Model for multiple longitudinal and multiple time-to-event traits using the two-stage approach with estimation of the full variance-covariance matrix by the bootstrap. 

This repository includes two main directories:

## ./data_simulation, includes:
- the script "data_simulation.r" used to generate R=1000 data replicates of N=667 DCCT individuals with K=2 simulated time-to-event traits (retinopathy, nephropathy) & simulated genotypes at M=5 causal SNPs(referred as "SNP1A to SNP5A" with direct and/or indirect effects on the time-to-event traits via L=3 longitudinal risk factors: 2 observed in DCCT (HbA1c, SBP) and one simulated (U), see Fig. 3 in Brossard et al. This script also simulates genotypes at 5 SNPs with the same MAFs as SNP1A-SNP5A under the global null genetic hypothesis (referred as "SNP1R to SNP5R").
- "DCCT_ARTIFICIAL_longQT.Rdata", an aritificial DCCT dataset with simulated longitudinal values & baseline covariates (SEX, T1D_diagnosis) in N=667 DCCT individuals, used as a replacement of the confidential DCCT dataset to illustrate the simulation procedure (see details in "generate_artificial_DCCT_longdata.r")
- "DCCT_based_simulated_data_reps5.Rdata", an example of output produced by  "data_simulation.r" for R=5 data replicates.

## ./JM_analysis, includes:
- the script "example_JMfit.r", which illustates an example of application of a Joint Model for L=2 longitudinal quantitative traits (HbA1c, SBP), and K=2 time-to-event traits (DR, DN) for causal SNP5A (with direct and indirect effects on time-to-DN via SBP) to one simulation replicate from "DCCT_based_simulated_data_reps5.Rdata"
- the script "example_JMfit_subfunctions.r" which includes the sub-functions used to fit the Joint Model 
- the script "mvlme.r" that includes the function mvlme() to fit a multivariate longitudinal mixed model, from the JoineRML R package (Version 0.4.2 on CRAN)  
- the output file "JM_results_SNP5A_replicate1.Rdata" produced by "example_JMfit.r", that includes the joint model results for SNP5 with comments in "example_JMfit.r for the interpretation of the output
