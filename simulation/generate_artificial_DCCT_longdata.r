# This script generates an artificial DCCT longitudinal dataset to illustrate our simulation procedure
# N=667 individuals, L=2 QTs (hba1c, and SBP) and 2 baseline factors (sex and T1D duration)

set.seed(13857)
ipak <- function(pkg){ sapply(pkg, require, character.only = TRUE) }
packages <- c( "iterators","foreach","doSNOW","snow","survival","eha","mvtnorm", "nlme","psychTools","MASS" )
ipak( packages )

# Original DCCT data - with longitudinal HbA1c, SBP values & baseline covariates 
long_data <- read.delim("./../../DCCT_data/2.Data_Processed/DCCT_longdata_for_simu_nomiss.tsv",sep=" ")[,-c(12,13)]
long_data$SEX <- long_data$SEX-1
DCCT_baseline<-unique(long_data[,c("PATIENT","SEX","DURATION_Years")])
DCCT_baseline$PATIENT_OLD<-DCCT_baseline$PATIENT
DCCT_long<-long_data[,c("PATIENT","obstime")]

# Generation of the aritificial longitudinal dataset that mimicks DCCT
n<-nrow(DCCT_baseline) ; i_n<-1
ARTIFICIAL_baseline<-as.data.frame(cbind(PATIENTOLD=DCCT_baseline$PATIENT,PATIENT=1:n,
	SEX=rbinom(n,1,mean(DCCT_baseline$SEX)),
	DURATION_Years=rnorm(n,mean(DCCT_baseline$DURATION_Years),sd(DCCT_baseline$DURATION_Years))))

ARTIFICIAL_long<- do.call(rbind,by(DCCT_long,
	INDICES = list(DCCT_long$PATIENT),function(x) {
	y<-ARTIFICIAL_baseline[which(ARTIFICIAL_baseline$PATIENTOLD==unique(x$PATIENT)),]
	x$PATIENT<-y$PATIENT
	x$SEX<-y$SEX
	x$DURATION_Years<-y$DURATION_Years
	return(x) }))
	
## Simulation of hba values
lmefitl1 <- lme( hba ~ obstime , random = ~ 1+obstime|PATIENT, 
                  data = long_data[,c("hba","obstime","PATIENT")],
                  control = lmeControl(opt = "optim"))
	  
parms_hba <- lapply(list(beta = fixef(lmefitl1), tau0=VarCorr(lmefitl1)[1,2],
	tau1=VarCorr(lmefitl1)[2,2], 	tau01=VarCorr(lmefitl1)[2,3], 
	sigma2=VarCorr(lmefitl1)[3,1]), function(x) { as.numeric(as.character(x)) })
	
ARTIFICIAL_long <- do.call(rbind,by(ARTIFICIAL_long,INDICES=list(ARTIFICIAL_long$PATIENT),
  FUN=function(U){
  	U<-U[order(U$obstime),]
  	n<-nrow(U)
  	X<-cbind(rep(1,n),sort(U$obstime))
  	Z<-cbind(rep(1,n),sort(U$obstime))
  	parms_hba$Sigma=matrix(c(parms_hba$tau0^2,parms_hba$tau01*parms_hba$tau0*parms_hba$tau1,
		parms_hba$tau01*parms_hba$tau0*parms_hba$tau1, parms_hba$tau1^2),nrow=2)
  	vcov<-Z%*%parms_hba$Sigma%*%t(Z)+parms_hba$sigma2*diag(n)
  	U$hba <- mvrnorm(1, mu= X%*%parms_hba$beta,  Sigma = vcov)
  	return(U)
  }
))

  
## Simulation of SBP values 
lmefitl2 <- lme( SBP ~ obstime + SEX , random = ~ 1+obstime|PATIENT,
                  data = long_data[,c("SBP","obstime","PATIENT","SEX")],
                  control = lmeControl(opt = "optim"))
				  
parms_sbp <- lapply(list(beta = fixef(lmefitl2), tau0=VarCorr(lmefitl2)[1,2], 
	tau1=VarCorr(lmefitl2)[2,2], tau01=VarCorr(lmefitl2)[2,3], 
	sigma2=VarCorr(lmefitl2)[3,1]), function(x) { as.numeric(as.character(x)) })
	
ARTIFICIAL_long <- do.call(rbind,by(ARTIFICIAL_long,INDICES=list(ARTIFICIAL_long$PATIENT),
  FUN=function(U){
  	U<-U[order(U$obstime),]
  	n<-nrow(U)
  	X<-cbind(rep(1,n),sort(U$obstime),U$SEX)
  	Z<-cbind(rep(1,n),sort(U$obstime))
  	parms_sbp$Sigma=matrix(c(parms_sbp$tau0^2,parms_sbp$tau01*parms_sbp$tau0*parms_sbp$tau1,
		parms_sbp$tau01*parms_sbp$tau0*parms_sbp$tau1, parms_sbp$tau1^2),nrow=2)
  	vcov<-Z%*%parms_sbp$Sigma%*%t(Z)+parms_sbp$sigma2*diag(n)
  	U$sbp <- mvrnorm(1, mu= X%*%parms_sbp$beta, Sigma = vcov)
  	return(U)
  }
))
 
rownames(ARTIFICIAL_long)<-NULL
# checking
ARTIFICIAL_lmefitl1 <- lme( hba ~ obstime , random = ~ 1+obstime|PATIENT, 
                  data = ARTIFICIAL_long[,c("hba","obstime","PATIENT")],
                  control = lmeControl(opt = "optim")) # checking for simulated hba

summary(ARTIFICIAL_lmefitl1) #simulated
summary(lmefitl1) #original

ARTIFICIAL_lmefitl2 <- lme( sbp ~ obstime + SEX , random = ~ 1+obstime|PATIENT,
                  data = ARTIFICIAL_long[,c("sbp","obstime","PATIENT","SEX")],
                  control = lmeControl(opt = "optim")) # checking for simulated sbp

summary(ARTIFICIAL_lmefitl2) # simulated 
summary(lmefitl2) # original

save(ARTIFICIAL_long, file="DCCT_ARTIFICIAL_longQT.Rdata")