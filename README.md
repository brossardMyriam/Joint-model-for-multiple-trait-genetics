# Joint-model-for-multiple-trait-genetics

Reference (paper under revision): "Characterization of direct and/or indirect genetic associations for multiple traits in longitudinal studies of disease progression" by Brossard M, Paterson AD, Espin-Garcia O, Craiu RV and Bull SB. 

This repository includes some Examples of R scripts to generate the DCCT-based simulated datasets under a multiple-trait genetics scenario (see Figs 3 and 4 in Brossard et al); & to fit the proposed Joint Model for multiple longitudinal and multiple time-to-event traits using a two-stage approach with estimation of the full variance-covariance matrix by the bootstrap. 

This repository includes two main directories:

## ./data_simulation:
- "data_simulation.r" is the the R script used to generate R=1000 data replicates of N=667 DCCT individuals with K=2 simulated time-to-T1D complications (retinopathy (referred as DR), nephropathy (referred as DN)) & simulated genotypes at M=5 causal SNPs (referred as "SNP1A to SNP5A") with direct and/or indirect effects on the time-to-event traits via L=3 longitudinal risk factors: 2 observed in DCCT (HbA1c, systolic blood pressure (SBP)) and one simulated (U), under scenario from Fig. 3 in Brossard et al. This script also simulates genotypes at 5 SNPs with the same MAFs as SNP1A-SNP5A but under the global null genetic hypothesis (referred as "SNP1R to SNP5R").
- "DCCT_ARTIFICIAL_longQT.Rdata", is an aritificial DCCT dataset with simulated longitudinal values & baseline covariates (SEX, T1D_diagnosis) in N=667 DCCT individuals, provided as a replacement of the confidential DCCT dataset to illustrate the simulation procedure (this dataset was generated using "generate_artificial_DCCT_longdata.r").
- "DCCT_based_simulated_data_reps5.Rdata", shows an example of output produced by "data_simulation.r" for R=5 data replicates.

## ./JM_analysis, includes:
- "example_JMfit.r", illustates an example of application of a Joint Model for L=2 longitudinal quantitative traits (HbA1c, SBP), and K=2 time-to-event traits (DR, DN) for causal SNP5 (with direct and indirect effects on time-to-DN via SBP) to one simulation replicate from "DCCT_based_simulated_data_reps5.Rdata".
- "example_JMfit_subfunctions.r", includes the sub-functions used to fit the proposed Joint Model.
- "mvlme.r", includes the function mvlme() to fit a multivariate longitudinal mixed model, from the JoineRML R package (Version 0.4.2 on CRAN)  
- "JM_results_SNP5A_replicate1.Rdata", is an example of output file produced by "example_JMfit.r". It includes the joint model results for SNP5 with comments in "example_JMfit.r" for the interpretation of the output.
