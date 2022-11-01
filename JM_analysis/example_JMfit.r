# Reference: Brossard  et al (manucript under revision) -  “Characterization of direct and/or indirect genetic associations for multiple traits in longitudinal studies of disease progression”. 
# Illustration of an application of the joint model for L=2 longitudinal quantitative (HbA1c, SBP) & K=2 time-to-event traits in the simulated dataset with data_simulation.r
# with Sex as a covariate for SBP & T1D_duration as a covariate for both time-to-event traits 
# As described in Brossard et al, we used a two-stage approach with estimation of the variance-covariance matrix by bootstrap 
# For this illustration, we focus on SNP5 simulated with a direct effect on time-to-DN (k=2) and an indirect effect via SBP longitudinal risk factor (l=2) - see Fig 3
# Sub-models specified in this application of joint model 
# l1: hbac ~ bi01 + b01 +(b11 + bi11)*obstime + bg1*SNP5A
# l2: sbp ~ bi02 + b02 +(b12 + bi12)*obstime + b22*sex_2 +bg2*SNP5A
# k1 : DR ~ gammag1*SNP5A + alpha_l1k1*hbac_fitted_l1k1 + alpha_Durk1*DURATION_Years
# k2 : DN ~ gammag2*SNP5A + alpha_l1k2*hbac_fitted_l1k2 + alpha_l2k2*sbpc_fitted_l2k2+ alpha_Durk2*DURATION_Years
# Note: The following functions need to be modified for other model specifications than the one presented here

set.seed(20221028)
ipak <- function(pkg){ sapply(pkg, require, character.only = TRUE) }
packages <- c("survival", "nlme", "aod", "numDeriv", "parallel")
ipak(packages)
source("mvlme.r") 
source("example_JMfit_subfunctions.r")

# Illustration for replicate 1 
load("DCCT_based_simulated_data_reps5.Rdata") # dataset generated from DCCT data using the script simulation.r
datalong_rep<-subset(datalong,datalong$rep==5) # selection of data replicate 1 (for illustration)
datasurv_rep<-subset(datasurv,datasurv$rep==5) 
datalong_rep$SNP<-datalong_rep$SNP5A #tested SNP
datasurv_rep$SNP<-datasurv_rep$SNP5A 

###############################################
# 	Stage 1 of the JM : Bivariate linear mixed model for hba1c (l=1), SBP (l=2)
# 	formLongFixed = Expresion for fixed effects for each QT
# 	formLongRandom = Expresion for random effects for each QT

stage1fit <- mvlme(
	formLongFixed = list("l1" = hbac ~ obstime + SNP, 
		"l2" = sbpc ~ obstime + SNP + SEX), 
	formLongRandom = list("l1" = hbac ~ 1+ obstime | PATIENT , 
		"l2" = sbpc ~ 1+ obstime | PATIENT ),
	data=datalong_rep, timeVar="obstime")

###############################################    
# Stage 2 of the JM : Bivariate Cox PH frailty model for time-to-DR (k=1) and time-to-DN (k=2)
# using fitted trajectory values for HbA1c & SBP 
# SNP_k1, Fitted_l1k1, DURATION_Years_k1: effects of SNP5, fitted  HbA1c trajectory & T1D duration on DR
# SNP_k2, Fitted_l1k2,Fitted_l2k2, DURATION_Years_k2: effects of SNP5, fitted  HbA1c & SBP trajectories & T1D duration on DN

stage2data<-prepare_for_stage2(stage1fit,datasurv_rep) # transform the survival data in long format with risk intervals
stage2fit <- coxph(formula = Surv(start, end, event_all) ~ SNP_k1 +  Fitted_l1k1 + DURATION_Years_k1 +
	SNP_k2 + Fitted_l1k2 +  Fitted_l2k2 +  DURATION_Years_k2 + 
	strata(event_type) + frailty(PATIENT,"gamma"), data = stage2data, 
	control = coxph.control(timefix = FALSE))
	
###############################################    
# Extraction of the results & estimation of the variance-covariance matrix 

result<-JM_res(stage1fit, stage2fit, datalong_rep, datasurv_rep, nboots=500)
save(result, file="JM_results_SNP5A_replicate1.Rdata")

# result
# $stage1 - estimates from the bivariate mixed model, bootstrap se & 1df Wald test P-values (based on the bootstrap se)
# $stage1$coef
                     # beta    lbootse    lwald1f.p
# (Intercept)_1 -0.04951817 0.05941478 4.046013e-01
# obstime_1      0.02489811 0.01055763 1.835860e-02
# SNP_1         -0.11138793 0.10344558 2.815794e-01
# (Intercept)_2 -7.07683964 0.51752423 1.443689e-42
# obstime_2      0.36899204 0.05392508 7.772544e-12
# SNP_2          6.47403000 0.63912584 4.087024e-24
# SEX_2          6.88572441 0.57545186 5.369738e-33

# $stage1$D - variance-covariance matrix of the random-effects coefficients in the bivariate mixed modell
              # (Intercept)_1   obstime_1 (Intercept)_2  obstime_2
# (Intercept)_1     1.9577957 -0.18192947     0.1833524 -0.0807241
# obstime_1        -0.1819295  0.06870721     0.0353017 -0.0155427
# (Intercept)_2     0.1833524  0.03530170    60.4719341 -3.3857045
# obstime_2        -0.0807241 -0.01554270    -3.3857045  1.2105402


# $stage2 - effect estimates in the Cox PH frailty model, bootstrap standard errors and Wald 1df test (based on the bootstrap standard errors)
                    # survcoef    sbootse    swald1f.p
# SNP_k1            0.06096215 0.12824686 6.345375e-01
# Fitted_l1k1       0.19225259 0.05096604 1.618336e-04
# DURATION_Years_k1 0.19084441 0.01775585 6.039171e-27
# SNP_k2            0.48660322 0.19147844 1.104426e-02
# Fitted_l1k2       0.17223678 0.07771523 2.667403e-02
# Fitted_l2k2       0.16914690 0.01663513 2.753537e-24
# DURATION_Years_k2 0.20837112 0.02876245 4.338668e-13

# $vcov - bootstrap variance-covariance matrix 
                  # (Intercept)_1     obstime_1         SNP_1 (Intercept)_2
# (Intercept)_1      3.530116e-03 -2.345077e-04 -2.933860e-03 -8.712009e-04
# obstime_1         -2.345077e-04  1.114636e-04 -1.017034e-04 -1.259363e-04
# SNP_1             -2.933860e-03 -1.017034e-04  1.070099e-02  2.587941e-03
# (Intercept)_2     -8.712009e-04 -1.259363e-04  2.587941e-03  2.678313e-01
# obstime_2         -1.450532e-04 -1.839108e-05  9.569905e-05 -7.787932e-03
# SNP_2              2.165702e-04  2.221937e-04  7.526049e-04 -1.291881e-01
# SEX_2              4.155480e-04  2.551929e-04 -2.520435e-03 -2.035682e-01
# SNP_k1             5.485556e-05 -2.367739e-05  6.511674e-04  6.261538e-03
# Fitted_l1k1        5.668633e-05  2.515461e-05 -7.454485e-04 -4.675503e-03
# DURATION_Years_k1  1.377878e-05 -8.471391e-06 -9.379659e-06  2.537439e-04
# SNP_k2            -2.777731e-04  1.137379e-04  1.301725e-03 -2.724497e-06
# Fitted_l1k2       -2.494058e-04  3.208341e-05  2.430283e-04 -2.728875e-03
# Fitted_l2k2        6.484854e-06  1.200301e-05 -2.029715e-04 -4.313562e-04
# DURATION_Years_k2  6.514494e-05  1.670181e-06  1.502479e-04  4.624762e-04
                      # obstime_2         SNP_2         SEX_2        SNP_k1
# (Intercept)_1     -1.450532e-04  0.0002165702  0.0004155480  5.485556e-05
# obstime_1         -1.839108e-05  0.0002221937  0.0002551929 -2.367739e-05
# SNP_1              9.569905e-05  0.0007526049 -0.0025204352  6.511674e-04
# (Intercept)_2     -7.787932e-03 -0.1291880696 -0.2035682402  6.261538e-03
# obstime_2          2.907914e-03  0.0007216405  0.0008210299 -1.264004e-04
# SNP_2              7.216405e-04  0.4084818345  0.0249103459 -2.460664e-03
# SEX_2              8.210299e-04  0.0249103459  0.3311448454 -9.247993e-03
# SNP_k1            -1.264004e-04 -0.0024606640 -0.0092479930  1.644726e-02
# Fitted_l1k1        7.940242e-06  0.0026215172  0.0029384573 -1.346323e-04
# DURATION_Years_k1  1.178299e-06 -0.0001461241 -0.0003162944  9.661791e-05
# SNP_k2             4.620263e-05 -0.0039527246  0.0091805152  3.958142e-03
# Fitted_l1k2       -8.370090e-05 -0.0026657026  0.0027324650 -3.246458e-04
# Fitted_l2k2        3.788119e-05  0.0006045565 -0.0005878866 -1.140450e-04
# DURATION_Years_k2 -1.110315e-04  0.0017897329 -0.0007104484  2.385072e-04
                    # Fitted_l1k1 DURATION_Years_k1        SNP_k2   Fitted_l1k2
# (Intercept)_1      5.668633e-05      1.377878e-05 -2.777731e-04 -2.494058e-04
# obstime_1          2.515461e-05     -8.471391e-06  1.137379e-04  3.208341e-05
# SNP_1             -7.454485e-04     -9.379659e-06  1.301725e-03  2.430283e-04
# (Intercept)_2     -4.675503e-03      2.537439e-04 -2.724497e-06 -2.728875e-03
# obstime_2          7.940242e-06      1.178299e-06  4.620263e-05 -8.370090e-05
# SNP_2              2.621517e-03     -1.461241e-04 -3.952725e-03 -2.665703e-03
# SEX_2              2.938457e-03     -3.162944e-04  9.180515e-03  2.732465e-03
# SNP_k1            -1.346323e-04      9.661791e-05  3.958142e-03 -3.246458e-04
# Fitted_l1k1        2.597537e-03      4.445676e-05 -4.912237e-04  2.661445e-04
# DURATION_Years_k1  4.445676e-05      3.152701e-04  1.823158e-05  3.228268e-05
# SNP_k2            -4.912237e-04      1.823158e-05  3.666399e-02  2.023052e-03
# Fitted_l1k2        2.661445e-04      3.228268e-05  2.023052e-03  6.039657e-03
# Fitted_l2k2        3.412358e-05      2.530808e-05 -8.747115e-04  1.244823e-04
# DURATION_Years_k2 -9.193203e-05      7.649288e-05  7.264560e-04  2.230159e-04
                    # Fitted_l2k2 DURATION_Years_k2
# (Intercept)_1      6.484854e-06      6.514494e-05
# obstime_1          1.200301e-05      1.670181e-06
# SNP_1             -2.029715e-04      1.502479e-04
# (Intercept)_2     -4.313562e-04      4.624762e-04
# obstime_2          3.788119e-05     -1.110315e-04
# SNP_2              6.045565e-04      1.789733e-03
# SEX_2             -5.878866e-04     -7.104484e-04
# SNP_k1            -1.140450e-04      2.385072e-04
# Fitted_l1k1        3.412358e-05     -9.193203e-05
# DURATION_Years_k1  2.530808e-05      7.649288e-05
# SNP_k2            -8.747115e-04      7.264560e-04
# Fitted_l1k2        1.244823e-04      2.230159e-04
# Fitted_l2k2        2.767275e-04      1.557240e-04
# DURATION_Years_k2  1.557240e-04      8.272787e-04
