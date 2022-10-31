# Brossard  et al (manucript under revision) -  “Characterization of direct and/or indirect genetic associations for multiple traits in longitudinal studies of disease progression”. 
# Example of R script to generate the R=1000 replicates of DCCT-based simulated datasets under the scenario from Fig. 3 using the simulation algorithm described in Fig.4 (with details in Supplementary Information 1)
# Note: For this illustration, we have replaced the original DCCT data by an artitificial dataset simulated based on DCCT

set.seed(20221028)
ipak <- function(pkg){ sapply(pkg, require, character.only = TRUE) }
packages <- c( "iterators","foreach","doSNOW","snow","survival","eha","mvtnorm", "nlme","psychTools","MASS" )
ipak( packages )

load("DCCT_ARTIFICIAL_longQT.Rdata") # Using artificial DCCT data with longitudinal trait values (for illustration purposes)  
long_data <- ARTIFICIAL_long
N <- length(unique(long_data$PATIENT))
long_data$hbac <- long_data$hba-mean(long_data$hba) # centered
long_data$sbpc <- long_data$sbp-mean(long_data$sbp) # centered 

## Parameters 
R <- 5 # no. of data replicates to generate, set to 5 for illustration
maf_SNP1 <- 0.3  
beta_SNP1 <- 0.7
alpha_l1k1 <- alpha_l1k2 <- 0.2

maf_SNP2 <- 0.1
gamma_SNP2 <- 0.8

maf_SNP3 <- 0.4
beta_SNP3 <- 0.8
alpha_Uk1 <- alpha_Uk2 <- 0.4

maf_SNP4 <- 0.3
gamma_SNP4 <- 0.7

maf_SNP5 <- 0.2
beta_SNP5 <- 7
gamma_SNP5 <- 0.1
alpha_l2k2 <- 0.2

## Fit the LMM to determine the DCCT-based parameters values to specify for the simulation study
lmefitl1 <- lme( hbac ~ obstime, random = ~ 1+obstime|PATIENT, 
                  data = long_data, control = lmeControl(opt = "optim"))
lmefitl2 <- lme( sbpc ~ obstime + SEX , random = ~ 1+obstime|PATIENT,
                  data = long_data, control = lmeControl(opt = "optim"))

## Specification of parameters values for SNP1 effect on QT1 (HbA1c)
parms_SNP1A <-list(
	beta = c(beta_0=fixef(lmefitl1)[1],beta_t=fixef(lmefitl1)[2],beta_g=beta_SNP1),
	Sigma = getVarCov( lmefitl1)[1:2,1:2], 
	sigma2 = as.numeric(VarCorr(lmefitl1)[3,1]))

## Specification of parameters values for SNP5 effect on QT2 (SBP)
parms_SNP5A <-list( 
	beta = c(beta_0=fixef(lmefitl2)[1], beta_t=fixef(lmefitl2)[2],
		beta_sex=fixef(lmefitl2)[3], beta_g=beta_SNP5), 
	Sigma = getVarCov( lmefitl2 )[1:2,1:2], 
	sigma2 = as.numeric(VarCorr(lmefitl2)[3,1]))


### Step 1 - Simulation of the shared longitudinal risk factor U (unobserved longitudinal trait)
######

parms_U <- list(beta = c(beta0=12, beta1=0.09), tau0=1.52, tau1=0.3, tau01=-0.5, sigma2=0.5^2)
long_data <- do.call(rbind,by(long_data,INDICES=list(long_data$PATIENT),
  FUN=function(U){
  	U<-U[order(U$obstime),]
  	n<-nrow(U)
  	X<-cbind(rep(1,n),sort(U$obstime))
  	Z<-cbind(rep(1,n),sort(U$obstime))
  	parms_U$Sigma=matrix(c(parms_U$tau0^2,parms_U$tau01*parms_U$tau0*parms_U$tau1,
		parms_U$tau01*parms_U$tau0*parms_U$tau1, parms_U$tau1^2),nrow=2)
  	vcov<-Z%*%parms_U$Sigma%*%t(Z)+parms_U$sigma2*diag(n)
  	U$sim_U <- mvrnorm(1, mu= X%*%parms_U$beta, Sigma = vcov)
  	return(U)
}))

long_data$sim_Uc <- long_data$sim_U-mean(long_data$sim_U) 
lmefitU <- lme( sim_Uc ~ obstime, random = ~ 1+obstime|PATIENT, 
                data = long_data[,c("sim_Uc","obstime","PATIENT")],
                control = lmeControl(opt = "optim") )

### Set up parameters values for SNP3 effect on unmeasured longitudinal risk factor (U)
parms_SNP3A <-list(
	beta = c( beta_0=fixef(lmefitU)[1], beta_t=fixef(lmefitU)[2], beta_g=beta_SNP3 ), 
	Sigma = getVarCov( lmefitU)[1:2,1:2], 
	sigma2 = as.numeric(VarCorr(lmefitU)[3,1]))

				
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

#############################################
# Simulation of the R replicates of genotypes for SNPs with indirect effects (SNP1, SNP3 & SNP5) from observed HbA1c, SBP and simulated U values
SNP1_replicates <- do.call(rbind, by(long_data, INDICES = list(long_data$PATIENT),
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

SNP3_replicates <- do.call(rbind,by(long_data, INDICES = list(long_data$PATIENT),
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
 
SNP5_replicates<-do.call(rbind,by(long_data,INDICES = list(long_data$PATIENT),
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
   
# checking
# lmefitl1 <- lme( hbac ~ obstime + SNP2A , random = ~ 1+obstime|PATIENT, 
                  # data = datalong, control = lmeControl(opt = "optim"))
# lmefitl2 <- lme( sbpc ~ obstime + SEX + SNP5A , random = ~ 1+obstime|PATIENT,
                  # data = datalong, control = lmeControl(opt = "optim"))

#############################################
# Simulation of replicates of genotypes for SNPs with direct effects (SNP2, SNP4) & for the SNPs under the global null independently from longitudinal trait values

N <- length(unique(long_data$PATIENT))
otherSNPs_replicates<-do.call(rbind,lapply(1:R,function(i){data.frame(rep=i,PATIENT=unique(long_data$PATIENT),
	SNP2A = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP2)^2, 2*maf_SNP2*(1-maf_SNP2), maf_SNP2^2)) ),
	SNP4A = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP4)^2, 2*maf_SNP4*(1-maf_SNP4), maf_SNP4^2)) ),
	
	SNP1R = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP1)^2, 2*maf_SNP1*(1-maf_SNP1), maf_SNP1^2)) ),
	SNP2R = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP2)^2, 2*maf_SNP2*(1-maf_SNP2), maf_SNP2^2)) ),
	SNP3R = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP3)^2, 2*maf_SNP3*(1-maf_SNP3), maf_SNP3^2)) ), 
	SNP4R = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP4)^2, 2*maf_SNP4*(1-maf_SNP4), maf_SNP4^2)) ),
	SNP5R = t(c(0,1,2)%*%rmultinom( N, 1, c((1-maf_SNP5)^2, 2*maf_SNP5*(1-maf_SNP5), maf_SNP5^2))))} )
	)

tmp1<-merge(SNP1_replicates, SNP3_replicates, by=c("PATIENT", "rep"))
tmp2<-merge(tmp1, SNP5_replicates, by=c("PATIENT", "rep"))
tmp3<-merge(tmp2, otherSNPs_replicates, by=c("PATIENT", "rep"))
datalong<-merge(long_data, tmp3, by=c("PATIENT"))
datalong<-subset(datalong, select=c("PATIENT", "obstime", "SEX", "DURATION_Years", "hbac","sbpc", "sim_Uc","rep",paste("SNP",1:5,"A",sep=""),paste("SNP",1:5,"R",sep="")))
datalong<-dfOrder(datalong,c("rep","PATIENT","obstime"))
rm(tmp1, tmp2, tmp3)

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
  alpha_l1k1=alpha_l1k1, # effect of QT l=1 (hba1c) on TTE trait k=1 (DR)
  alpha_Uk1=alpha_Uk1, # effect of QT l=u on TTE trait k=1 (DR)
  alpha_t1dur_k1=0.2, # effect of Duration_Years on TTE trait k=1 (DR)
  phi_k2=1.01, 
  xi_k2=0.01, 
  alpha_l1k2=alpha_l1k2, # effect of l=1 (hba1c) on TTE trait k=2 (DN)
  alpha_l2k2=alpha_l1k2, # effect of l=2 (sbp) on TTE trait k=2 (DN)
  alpha_Uk2=alpha_Uk2, # effect of u on TTE trait k=2 (DN)
  alpha_t1dur_k2=0.2
)

datasurv <-  survival_sim <- NULL

## Simulation for each patient i & for each replicate r
for ( r in 1:R ) {  
  
  datalong.byrep <- datalong[ which(datalong$rep==r) ,]
  
  lmefitl1 <- lme( hbac  ~ obstime + SNP1A , random = ~ 1+obstime|PATIENT, 
                 control = lmeControl(opt = "optim"), data = datalong.byrep) 
     
  lmefitl2 <- lme( sbpc  ~ obstime + SNP5A + SEX, random = ~ 1+obstime|PATIENT, 
                 control = lmeControl(opt = "optim"), data = datalong.byrep) 
    
  lmefitU <- lme( sim_Uc ~ obstime + SNP3A , random = ~ 1+obstime|PATIENT, 
                  control = lmeControl(opt = "optim"), data =datalong.byrep) 
  
 data.byrep<-unique(datalong.byrep[,c("PATIENT", "SEX", "DURATION_Years", 
	paste("SNP",1:5,"A",sep=""),paste("SNP",1:5,"R",sep=""))])
 
  survival_sim <- do.call(rbind,by(data.byrep, INDICES = list(data.byrep$PATIENT),
        function(data, parms, tmin, tmax){
			all_iter<-NULL
                             
            lmefitl1.coefsPAT<-unname(as.numeric(as.character(coef(lmefitl1)[data$PATIENT,])))
			traj_l1<-function(t) { lmefitl1.coefsPAT[1] + lmefitl1.coefsPAT[2]*t + lmefitl1.coefsPAT[3]*data$SNP1A }
                                     
            lmefitl2.coefsPAT<-unname(as.numeric(as.character(coef(lmefitl2)[data$PATIENT,])))
            traj_l2 <- function(t) { lmefitl2.coefsPAT[1]+ lmefitl2.coefsPAT[2]*t + lmefitl2.coefsPAT[3]*data$SNP5A + lmefitl2.coefsPAT[4]*data$SEX	}
                                     
            lmefitU.coefsPAT<-unname(as.numeric(as.character(coef(lmefitU)[data$PATIENT,])))
            traj_U<-function(t) { lmefitU.coefsPAT[1] + lmefitU.coefsPAT[2]*t + lmefitU.coefsPAT[3]*data$SNP3A }
                                     
            ## Simulation of the Time-to-retinopathy outcome (k=1)
            baseline_k1<-function(t) {	parms_K$xi_k1*parms_K$phi_k1*t^(parms_K$phi_k1-1) }
            hazard_k1 <- function(t) { baseline_k1(t)*exp( parms_K$alpha_l1k1*traj_l1(t) + parms_K$alpha_t1dur_k1*data$DURATION_Years + parms_K$gamma_SNP2*data$SNP2A + parms_K$alpha_Uk1*traj_U(t) )} 
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
            hazard_k2 <- function(t) { baseline_k2(t)*exp( parms_K$alpha_l1k1*traj_l1(t) + parms_K$alpha_l2k2*traj_l2(t) + parms_K$alpha_Uk2*traj_U(t) +  parms_K$alpha_t1dur_k2*data$DURATION_Years + parms_K$gamma_SNP2*data$SNP2A + parms_K$gamma_SNP5*data$SNP5A)} 
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
            maxobst <- max(subset(datalong.byrep,datalong.byrep$PATIENT==data$PATIENT)$obstime)
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

colnames(datasurv) <- c("PATIENT", "rep", "event_k1", "tobs_k1", "event_k2", "tobs_k2") 
datasurv <- merge(unique(datalong[,c("PATIENT", "rep", "SEX", "DURATION_Years",
	paste("SNP",1:5,"A",sep=""),paste("SNP",1:5,"R",sep=""))]), datasurv, by=c("rep","PATIENT"))
datasurv<-dfOrder(datasurv,c("rep","PATIENT"))
dim(datasurv)
save(datalong, datasurv, file=paste("DCCT_based_simulated_data_reps", R, ".Rdata", sep=""))

# checking
# summary(lme( hbac  ~ obstime + SNP1A , random = ~ 1+obstime|PATIENT, 
                 # control = lmeControl(opt = "optim"), data = datalong) )
     
# summary(lme( sbpc  ~ obstime + SNP5A + SEX, random = ~ 1+obstime|PATIENT, 
                 # control = lmeControl(opt = "optim"), data = datalong)) 
    
 # summary(lme( sim_Uc ~ obstime + SNP3A , random = ~ 1+obstime|PATIENT, 
                  # control = lmeControl(opt = "optim"), data =datalong))
				  
# coxph(formula = Surv(tobs_k1, event_k1) ~  SNP2A+ DURATION_Years, data = datasurv, 
	# control = coxph.control(timefix = FALSE))				  

# coxph(formula = Surv(tobs_k2, event_k2) ~  SNP5A+ DURATION_Years, data = datasurv, 
	# control = coxph.control(timefix = FALSE))				  