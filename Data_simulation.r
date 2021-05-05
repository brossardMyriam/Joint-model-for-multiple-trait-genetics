# Brossard  et al
# Example of R script to generate the DCCT-data based simulation study under the causal (Figure 2) and under the global null genetic scenarios
# M=5 causal SNPs with direct associations on K=2 simulated time-to-event traits &/or indirect associations via L=3 longitudinal risk factors: 2 observed in DCCT subjects (HbA1c, SBP) & one simulated (U)

args=(commandArgs(TRUE))
print(args)
for(k in 1:length(args)){
	eval(parse(text=args[[k]]))
}

set.seed(13857)
dir_wk <- c("/home/bulllab/mbrossard/Joint_models/Example_Rscripts_paper/simulation/")
setwd(dir_wk)
ipak <- function(pkg){ sapply(pkg, require, character.only = TRUE) }
packages <- c( "survival","eha","mvtnorm", "nlme","psychTools","MASS" )
ipak( packages )

# DCCT longitudinal data processed 
long_data <- read.delim("/home/bulllab/mbrossard/Joint_models/DCCT_data/2.Data_Processed/DCCT_longdata_for_simu_nomiss.tsv",sep=" ")
long_data$SEX <- long_data$SEX-1
n <- length(unique(long_data$PATIENT))
long_data$hbac <- long_data$hba-mean(long_data$hba) # centered
long_data$sbpc <- long_data$SBP-mean(long_data$SBP) # centered 

## Parameters 
R <- reps

# Indirect via HbA1c
maf_SNP1 <- 0.3  
beta_SNP1 <- 0.7
alpha_l1k1 <- alpha_l1k2 <- 0.2

# Direct
maf_SNP2 <- 0.1
gamma_SNP2 <- 0.8

# Indirect via U
maf_SNP3 <- 0.4
beta_SNP3 <- 0.8
alpha_Uk1 <- alpha_Uk2 <- 0.4

# Direct
maf_SNP4 <- 0.3
gamma_SNP4 <- 0.7

# Direct & Indirect induced via SBP
maf_SNP5 <- 0.2
beta_SNP5 <- 7
gamma_SNP5 <- 0.1
alpha_l2k2 <- 0.2

## Fit the LMM in DCCT individuals(for parameters sepcification for the simulation study)
lmefitl1 <- lme( hbac ~ obstime , random = ~ 1+obstime|PATIENT, 
                  data = long_data[,c("hbac","obstime","PATIENT")],
                  control = lmeControl(opt = "optim"))
lmefitl2 <- lme( sbpc ~ obstime + SEX , random = ~ 1+obstime|PATIENT,
                  data = long_data[,c("sbpc","obstime","PATIENT","SEX")],
                  control = lmeControl(opt = "optim"))

## Set up parameters values in LMM1 for SNP1
beta_0<-summary(lmefitl1)$coefficients$fixed[1] 
beta_t<- summary(lmefitl1)$coefficients$fixed[2]  
parms_SNP1A <-list(
	beta = c(beta_0=beta_0,beta_t=beta_t,beta_g=beta_SNP1),
	Sigma = matrix(c(getVarCov( lmefitl1)[1,1], getVarCov(lmefitl1)[1,2], 
	getVarCov( lmefitl1)[1,2],getVarCov( lmefitl1)[2,2]),nrow=2), 
	sigma2 = as.numeric(VarCorr(lmefitl1)[3,1]))

## Set up parameters values in LMM2 for SNP5
beta_0 <-summary( lmefitl2 )$coefficients$fixed[1] 
beta_t <- summary( lmefitl2 )$coefficients$fixed[2]
beta_sex <- summary( lmefitl2 )$coefficients$fixed[3] 
parms_SNP5A <-list( 
	beta = c(beta_0, beta_t, beta_sex, beta_g=beta_SNP5), 
	Sigma = matrix(c(getVarCov( lmefitl2 )[1,1], getVarCov(lmefitl2)[1,2], 
	getVarCov(lmefitl2)[1,2], getVarCov(lmefitl2)[2,2]),nrow=2), 
	sigma2 = as.numeric(VarCorr(lmefitl2)[3,1]))

### Step 1 - Simulation of the shared longitudinal risk factor U (unobserved longitudinal trait)
######

parms_U <- list(beta = c(beta0=12, beta1=0.09), tau0=1.52, tau1=0.3, tau01=-0.5, sigma2=0.5^2)
simdata_QT <- do.call(rbind,by(long_data,INDICES=list(long_data$PATIENT),
  FUN=function(U){
  	U<-U[order(U$obstime),]
  	n<-nrow(U)
  	X<-cbind(rep(1,n),sort(U$obstime))
  	Z<-cbind(rep(1,n),sort(U$obstime))
  	parms_U$Sigma=matrix(c(parms_U$tau0^2,parms_U$tau01*parms_U$tau0*parms_U$tau1,parms_U$tau01*parms_U$tau0*parms_U$tau1, parms_U$tau1^2),nrow=2)
  	vcov<-Z%*%parms_U$Sigma%*%t(Z)+parms_U$sigma2*diag(n)
  	U$sim_U <- mvrnorm(1, mu= X%*%parms_U$beta,  Sigma = vcov)
  	return(U)
  }
  ))
simdata_QT$sim_Uc <- simdata_QT$sim_U-mean(simdata_QT$sim_U) 
parms_U <- list( beta = c(beta0=12, beta1=0.09, beta_g=beta_SNP3 ), tau0=1.52, tau1=0.3, tau01=-0.5, sigma2=0.5^2 )

lmefitU <- lme( sim_Uc ~ obstime, random = ~ 1+obstime|PATIENT, 
                data = simdata_QT[,c("sim_Uc","obstime","PATIENT")],
                control = lmeControl(opt = "optim") )

parms_SNP3A <-list(
	beta = c( beta_0, beta_t, beta_g=beta_SNP3 ), 
	Sigma = matrix( c(getVarCov( lmefitU)[1,1],getVarCov( lmefitU)[1,2],
	getVarCov( lmefitU)[1,2],getVarCov( lmefitU)[2,2]),nrow=2 ), 
	sigma2 = as.numeric( VarCorr(lmefitU)[3,1]) ) #parameters values in LMM3 for SNP3

				
### Step 2 - Simulation of the R replicates of the genotype data at M "causal" SNPs and M "null" SNPs
#######

pmf_g<-function( g, data, parms, maf, log=F, QT ){
  if( g %in% 0:2 ){
    pG <- ifelse(g==0,(1-maf)^2,ifelse(g==1,2*maf*(1-maf),maf^2)) 
    n_i <- nrow(data)  
    X <- cbind(rep(1,n_i), data$obstime, rep(g,n_i)) 
	Z <- cbind(rep(1,n_i),data$obstime) 
    vcov <- Z%*%parms$Sigma%*%t(Z)+parms$sigma2*diag(n_i) 
    f_G_Y <- (log==F)*(dmvnorm(data[,QT], mean = X%*%parms$beta, sigma = vcov)*pG) + 
		(log==T)*(dmvnorm(data[,QT], mean = X%*%parms$beta, sigma = vcov, log=T)+ log(pG))
  } else stop( "G not in the right range" )
  return(f_G_Y)
} # P(Y,G) = P(Y|G) U P(G)

pmf_l2<-function( g, data, parms, maf, log=F, QT ){
  if( g %in% 0:2 ){
    pG <- ifelse(g==0,(1-maf)^2,ifelse(g==1,2*maf*(1-maf),maf^2)) 
    n_i <- nrow(data)  
    X <- cbind(rep(1,n_i), data$obstime, rep(g,n_i), rep(unique(data$SEX),n_i)) 
	Z <- cbind(rep(1,n_i),data$obstime) 
    vcov <- Z%*%parms$Sigma%*%t(Z)+parms$sigma2*diag(n_i) 
	f_G_Y <- (log==F)*(dmvnorm(data[,QT], mean = X%*%parms$beta, sigma = vcov)*pG) + 
	    (log==T)*(dmvnorm(data[,QT], mean = X%*%parms$beta, sigma = vcov, log=T)+ log(pG))
  }else stop( "G not in the right range" )
  return( f_G_Y )
} # P(Y,G) = P(Y|G) U P(G) # specific to l=2 (includes effect of SEX)

Vpmf_g <- Vectorize( pmf_g, vectorize.args = c("g") ) #  P(Y|G=g) for each g 
Vpmf_sbp <- Vectorize( pmf_l2, vectorize.args = c("g") ) #  P(Y|G=g) for each g 

Generate_SNP1 <- do.call(rbind, by(simdata_QT, INDICES = list(simdata_QT$PATIENT),
	function(data, reps, parms,  maf){
		pA <- Vpmf_g( 0:2, data, parms, maf, QT="hbac" )
		if(any(pA==0)){ 
			p0A <- Vpmf_g( 0:2, data, parms, maf, log=T ) 
			p0bA <- p0A-maU( p0A ) 
			p1A <-eUp ( p0bA )/sum( eUp(p0bA) )
		} else p1A<-pA/sum(pA)
		data.frame(PATIENT=data$PATIENT[1],
			rep=1:reps,SNP1A=sample(0:2, reps, prob=p1A, rep=T)) 
	  },
reps=R, parms=parms_SNP1A, maf=maf_SNP1)) # Simulation of P(G|Y)=P(Y|G)/P(Y)

Generate_SNP3 <- do.call(rbind,by(simdata_QT, INDICES = list(simdata_QT$PATIENT),
  function(data, reps, parms, maf){
    pA <- Vpmf_g( 0:2, data, parms,  maf, QT="sim_Uc" )
    if(any(pA==0)){ 
		p0A <- Vpmf_g( 0:2, data, parms, maf, log=T ) 
		p0bA <- p0A-maU( p0A ) 
		p1A <- eUp( p0bA )/sum(eUp( p0bA ))
    } else p1A <- pA/sum(pA)
    data.frame(PATIENT=data$PATIENT[1], rep=1:reps,
        SNP3A=sample(0:2, reps, prob=p1A, rep=T)) 
  },reps=R, parms=parms_SNP3A, maf=maf_SNP3))
 
Generate_SNP5<-do.call(rbind,by(simdata_QT,INDICES = list(simdata_QT$PATIENT),
    function(data, reps, parms, maf){
    pA<-Vpmf_sbp(0:2, data, parms, maf, QT="sbpc")
    if(any(pA==0)){ 
        p0A <- Vpmf_sbp(0:2, data, parms, maf, log=T) 
        p0bA <- p0A-maU(p0A) 
        p1A <- eUp(p0bA)/sum(eUp(p0bA))
    } else p1A <- pA/sum(pA)
    data.frame(PATIENT=data$PATIENT[1], rep=1:reps,
              SNP5A=sample(0:2, reps, prob=p1A, rep=T))
   },reps=R, parms=parms_SNP5A, maf=maf_SNP5))

## Simulation of the M=5 SNPs under the null & SNP4 with a direct effect on k=2
N <- length(unique(simdata_QT$PATIENT))
dat_SNP<-do.call(rbind,lapply(1:R,function(i){data.frame(rep=i,PATIENT=unique(simdata_QT$PATIENT),
	SNP2A = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP2)^2, 2*maf_SNP2*(1-maf_SNP4), maf_SNP2^2)) ),
	SNP4A = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP4)^2, 2*maf_SNP4*(1-maf_SNP4), maf_SNP4^2)) ),
	SNP1R = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP1)^2, 2*maf_SNP1*(1-maf_SNP1), maf_SNP1^2)) ),
	SNP2R = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP1)^2, 2*maf_SNP1*(1-maf_SNP1), maf_SNP1^2)) ),
	SNP3R = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP3)^2, 2*maf_SNP3*(1-maf_SNP3), maf_SNP3^2)) ), 
	SNP4R = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP4)^2, 2*maf_SNP4*(1-maf_SNP4), maf_SNP4^2)) ),
	SNP5R = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP5)^2, 2*maf_SNP5*(1-maf_SNP5), maf_SNP5^2))))} )
	)

tmp1<-merge(Generate_SNP1, Generate_SNP3, by=c("PATIENT", "rep"))
tmp2<-merge(tmp1, Generate_SNP5, by=c("PATIENT", "rep"))
tmp3<-merge(tmp2, dat_SNP, by=c("PATIENT", "rep"))
tmp4<-merge(simdata_QT, tmp3, by=c("PATIENT"))
datalong<- tmp4[,c("rep","PATIENT", "obstime", "qv", "SEX", "GROUP","AGE", "DURATION_Years", 
	"hbac", "sbpc", "sim_Uc", "SNP1A", "SNP2A", "SNP3A", "SNP4A", "SNP5A", 
	"SNP1R", "SNP2R", "SNP3R", "SNP4R", "SNP5R")]
datalong<-dfOrder(datalong,c("rep","PATIENT","obstime"))
rm(tmp1, tmp2, tmp3, tmp4)

### Step 3: Simulation of R replicates of K=2 non-independent time-to-T1DC traits
# using separate Weibull  time-to-event models linked via trajectory of U
# and contemporaneous HbA1c & SBP values as association structures on T1DC traits

## Parameters for k=1 (time-to-DR) and k=2 (time-to-DN)
parms_K <-list(
  gamma_SNP2=gamma_SNP2,
  gamma_SNP4=gamma_SNP4, 
  gamma_SNP5=gamma_SNP5,
  phi_k1=1.01,
  xi_k1=0.1,
  alpha_l1k1=alpha_l1k1,
  alpha_Uk1=alpha_Uk1,
  alpha_t1dur_k1=0.2,
  phi_k2=1.01, 
  xi_k2=0.01, 
  alpha_l1k2=alpha_l1k2, 
  alpha_l2k2=alpha_l1k2,
  alpha_Uk2=alpha_Uk2, 
  alpha_t1dur_k2=0.2
)

datasurv <-  survival_sim <- NULL

## Simulation for each patient i & for each replicate r
# foreach ?
for ( r in 1:R ) {  
  
  datalong.byrep <- datalong[ which(datalong$rep==r) ,]
  
  lmefitl1 <- lme( hbac  ~ obstime + SNP1A , random = ~ 1+obstime|PATIENT, 
                 control = lmeControl(opt = "optim"), data = datalong.byrep) 
  lmefitl1.RE <- summary(lmefitl1)$coefficients$random$PATIENT
  lmefitl1.FE <- summary(lmefitl1)$coefficients$fixed
  
  lmefitl2 <- lme( sbpc  ~ obstime + SNP5A + SEX, random = ~ 1+obstime|PATIENT, 
                 control = lmeControl(opt = "optim"), data = datalong.byrep) 
  lmefitl2.RE <- summary(lmefitl2)$coefficients$random$PATIENT
  lmefitl2.FE <- summary(lmefitl2)$coefficients$fixed
  
  lmefitU <- lme( sim_Uc ~ obstime + SNP3A , random = ~ 1+obstime|PATIENT, 
                  control = lmeControl(opt = "optim"), data =datalong.byrep) 
  lmefitU.RE <- summary(lmefitU)$coefficients$random$PATIENT
  lmefitU.FE <- summary(lmefitU)$coefficients$fixed

 data.byrep<-unique(datalong.byrep[,c("PATIENT", "SEX", "DURATION_Years", colnames(datalong.byrep)[grep("SNP", colnames(datalong.byrep))])])
 
  survival_sim <- do.call(rbind,by(data.byrep, INDICES = list(data.byrep$PATIENT),
        function(data, parms, tmin, tmax){
			all_iter<-NULL
                             
            lmefitl1.REPAT<-lmefitl1.RE[which(rownames(lmefitl1.RE)==data$PATIENT),]	
            traj_l1<-function(t) { lmefitl1.FE[1] + lmefitl1.REPAT[1] + (lmefitl1.FE[2] + lmefitl1.REPAT[2])*t + lmefitl1.FE[3]*data$SNP1A }
                                     
            lmefitl2.REPAT<-lmefitl2.RE[which(rownames(lmefitl2.RE)==data$PATIENT),]
            traj_l2 <- function(t) { lmefitl2.FE[1]+ lmefitl2.REPAT[1] + (lmefitl2.FE[2] + lmefitl2.REPAT[2])*t + lmefitl2.FE[3]*data$SNP5A + lmefitl2.FE[4]*data$SEX	}
                                     
            lmefitU.REPAT<-lmefitU.RE[which(rownames(lmefitU.RE)==data$PATIENT),]
            traj_U<-function(t) { 
                    lmefitU.FE[1] + lmefitU.REPAT[1] + (lmefitU.FE[2] + lmefitU.REPAT[2])*t +
                    lmefitU.FE[3]*data$SNP3A }
                                     
            ## Simulation of the Time-to-retinopathy outcome (k=1)
            baseline_k1<-function(t) {	parms_K$xi_k1*parms_K$phi_k1*t^(parms_K$phi_k1-1) }
            hazard_k1 <- function(t) { 
                  baseline_k1(t)*exp( parms_K$alpha_l1k1*traj_l1(t) +  
      					  parms_K$alpha_t1dur_k1*data$DURATION_Years + 
      					  parms_K$gamma_SNP2*data$SNP2A + 
      					  parms_K$alpha_Uk1*traj_U(t) )} 
            HazCum_k1 <- function(t) { integrate(hazard_k1, lower=0, upper=t, subdivisions=10000, rel.tol=1.e-05)$value  }
            InvHazard_k1 <- function(HazCum_k1, hazard_k1, tmin, tmax) { 
                                U <- runif(1, 0, 1)
                                rootfct <- function(t) {   U - exp(-HazCum_k1(t)) }
                                root <- try(uniroot(rootfct, interval = c(tmin, tmax))$root, silent = TRUE)
                                root <- if (inherits(root, "try-error")) { tmax + 0.01} else { root }
                                return(root)								
                                     }                     
            time_k1 <- InvHazard_k1(HazCum_k1, hazard_k1, tmin, tmax) 
            cumhaz_k1 <- HazCum_k1(time_k1)
            survprob_k1 <- exp(-cumhaz_k1 )
                                     
            ## Simulation of the Time-to-nephropathy outcome (k=2) 
            baseline_k2 <- function(t) { parms_K$xi_k2*parms_K$phi_k2*t^(parms_K$phi_k2-1) }
            hazard_k2 <- function(t) { baseline_k2(t)*exp( parms_K$alpha_l1k2*traj_l1(t) + parms_K$alpha_l2k2*traj_l2(t) + parms_K$alpha_Uk2*traj_U(t) + 
                                parms_K$alpha_t1dur_k2*data$DURATION_Years + parms_K$gamma_SNP4*data$SNP4A + parms_K$gamma_SNP5*data$SNP5A)} 
            HazCum_k2 <- function(t) { integrate(hazard_k2, lower=0, upper=t, subdivisions=10000, rel.tol=1.e-05)$value  }
            InvHazard_k2 <- function(HazCum_k2, hazard_k2, tmin, tmax) { 
                                U <- runif(1, 0, 1)
                                rootfct <- function(t) {   U - exp(-HazCum_k2(t)) }
                                root <- try(uniroot(rootfct, interval = c(0, tmax))$root, silent = TRUE) 
                                root <- if (inherits(root, "try-error")) { tmax + 0.01} else { root }
                                return(root)								}                      
             time_k2<- InvHazard_k2(HazCum_k2, hazard_k2, tmin, tmax)
            cumhaz_k2 <- HazCum_k2(time_k2)
            survprob_k2 <- exp(-cumhaz_k2 )
                                     
            ## censoring time
            maxobst <- max(datalong.byrep[which(datalong.byrep$PATIENT==data$PATIENT),]$obstime)
            cens <- runif(1,min=0,max=maxobst)
            event_k1 <- as.numeric(time_k1 <= cens )
            t_obs_k1 <- min(time_k1,cens)
            event_k2 <- as.numeric(time_k2 <= cens )
			t_obs_k2 <- min(time_k2,cens)
                                     
            subj_iter<-data.frame(data$PATIENT, r, event_k1, t_obs_k1, event_k2, t_obs_k2 )
            all_iter<-rbind(all_iter, subj_iter)
                                     
        }, parms=parms_K, tmin=1e-5, tmax=100))
  
  datasurv <- rbind(datasurv, survival_sim)
}

colnames(datasurv) <- c("PATIENT", "rep", "event_k1", "tobs_k1", "event_k2", "tobs_k2") #merge with the SNPs, SNPR & covariates
dim(datasurv)
datasurv <- merge(unique(datalong[,c("rep", "PATIENT", "SEX", "DURATION_Years", colnames(datalong)[grep("SNP", colnames(datalong))])]), datasurv, by=c("PATIENT", "rep"))
dim(datasurv)
datasurv<-dfOrder(datasurv,c("rep","PATIENT"))
dim(datasurv)
dim(datalong)
save(datalong, datasurv, file=paste("simulated_data_reps", R, ".Rdata", sep=""))





