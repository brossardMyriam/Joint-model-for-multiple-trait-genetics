# Brossard  et al (manucript under revision) -  “Characterization of direct and/or indirect genetic associations for multiple traits in longitudinal studies of disease progression”. 

# Sub-functions() used in example_JMfit.r

# prepare_for_stage2() 
# inputs: output from stage1, survival data (datasurv_rep)
# calculate fitted trajectories for HbA1c, SBP at the beginning of each risk interval
# prepare a stacked dataset in long format for the K=2 time-to-event traits & define specific variables 
# 	For k1: SNP_k1, Fitted_l1k1, DURATION_Years_k1; those variables are all set to 0 for k=2
# 	For k2: SNP_k2, Fitted_l1k2, Fitted_l2k2, DURATION_Years_k2 ; those variables are all set to 0 for k=1
# output: a dataset in longformat for stage2 stacked for the K=2 time-to-event traits

prepare_for_stage2<- function(stage1fit,datasurv_rep){			
	beta <-stage1fit$coefs.beta[c("SNP_1","SNP_2"),c(1,2,4)]
	n<-length(unique(datasurv_rep$PATIENT))
	tmp <- rep(list(stage1fit$coefs.beta[,1]), n)
	names(tmp) <- unique(datasurv_rep$PATIENT)
	lmmfr.i <- Map(c,stage1fit$b[[1]],tmp)
	
	data_k1 <- datasurv_rep[,c("PATIENT", "rep", "event_k1", "tobs_k1", "SEX", "DURATION_Years","SNP")]
	data_k1$tobs_k1[data_k1$event_k1==1 & data_k1$tobs_k1==0] <- 1e-5
	data_k2 <- datasurv_rep[,c("PATIENT","rep","event_k2","tobs_k2","SEX", "DURATION_Years","SNP")]
	data_k2$tobs_k2[data_k2$event_k2==1 & data_k2$tobs_k2==0] <- 1e-5
	cut.pointsK <- unique(c(data_k1$tobs_k1[data_k1$event_k1==1],data_k2$tobs_k2[data_k2$event_k2==1]))
	data_k1.long <- survSplit(data=data_k1,cut=cut.pointsK, end = "tobs_k1", start = "start", event = "event_k1")
	colnames(data_k1.long)[which(colnames(data_k1.long)=="tobs_k1")] <- "end"
	colnames(data_k1.long)[which(colnames(data_k1.long)=="event_k1")] <- "event_all"
	data_k1.long$event_type <- 0
	
	data_k2.long <- survSplit(data=data_k2,cut=cut.pointsK, 
			end = "tobs_k2", start = "start", event = "event_k2")
	colnames(data_k2.long)[which(colnames(data_k2.long)=="tobs_k2")] <- "end"
	colnames(data_k2.long)[which(colnames(data_k2.long)=="event_k2")] <- "event_all"
	data_k2.long$event_type <- 1

	# k1,k2 stacked			
	data_K.long <- rbind(data_k1.long,data_k2.long)
	rm(data_k2.long) ; rm(data_k1.long)
	rm(data_k2); rm(data_k1)
			 
	# fitted QT values at start of each interval 		 
	stage2data <- do.call(rbind, by(data_K.long, INDICES=list(data_K.long$PATIENT), function (x) {
		lme.coef <- lmmfr.i
		lme.coef <- unlist(lmmfr.i[which(names(lmmfr.i)==unique(x$PATIENT))])
							  
		Traj_l1 <- function(x,t) { lme.coef[5] + lme.coef[1] + ( lme.coef[6]+lme.coef[2] )*t + lme.coef[7]*x$SNP  }
		Traj_l2 <- function(x,t) { lme.coef[8] + lme.coef[3] + ( lme.coef[4]+lme.coef[9] )*t + lme.coef[10]*x$SNP + lme.coef[11]*x$SEX }

		return( cbind(x,Fitted_l1=Traj_l1(x, x$start),Fitted_l2=Traj_l2(x, x$start)))
	}))

	rm(data_K.long)
				  
	# Coding of variables for event of type k=1 (retinopathy)
	stage2data$SNP_k1 <- stage2data$SNP*(1-stage2data$event_type)
	stage2data$DURATION_Years_k1 <- stage2data$DURATION_Years*(1-stage2data$event_type)
	stage2data$Fitted_l1k1 <- stage2data$Fitted_l1*(1-stage2data$event_type)
	stage2data$Fitted_l2k1 <- stage2data$Fitted_l2*(1-stage2data$event_type)
							   
	# Coding of variables for event of type k=2 (nephropathy)
	stage2data$SNP_k2 <- stage2data$SNP*stage2data$event_type
	stage2data$DURATION_Years_k2 <- stage2data$DURATION_Years*stage2data$event_type
	stage2data$Fitted_l1k2 <- stage2data$Fitted_l1*stage2data$event_type
	stage2data$Fitted_l2k2 <- stage2data$Fitted_l2*stage2data$event_type
	return( stage2data ) 
}	

# JM_bcov()
# inputs: longitudinal (datalong) & survival data (datasurv ) matrices used to & nboots=no. of bootstraps 
# output: bootstrap covariate matrix of  the stage1 & stage2 parameters

JM_bcov<-function(datalong_rep, datasurv_rep, nboots) {
	return(cov(
	do.call('rbind',(mclapply(1:nboots, mc.cores = parallel::detectCores()-2, function(boots) {
		subjID <- unique(datalong_rep$PATIENT)
		bootcur <- sample(subjID,length(subjID),replace=TRUE)
		boot.id <- as.data.frame(cbind(PATIENT=rep(1:length(bootcur)),oldid=bootcur))
		datalong.boot <- do.call(rbind, by(boot.id, INDICES=list(boot.id$PATIENT), function (x) {
			return(cbind(x$PATIENT, datalong_rep[which(datalong_rep$PATIENT==unique(x$oldid)),]))		}))
		datasurv.boot <-  do.call(rbind, by(boot.id, INDICES=list(boot.id$PATIENT), function (x) {
			return(cbind(x$PATIENT, datasurv_rep[which(datasurv_rep$PATIENT==unique(x$oldid)),]))		}))
		names(datalong.boot)[c(1,2)] <- names(datasurv.boot)[c(1,3)] <- c("PATIENT","PATIENT.OLD")
		stage1fitb <- mvlme(
			formLongFixed = list("l1" = hbac ~ obstime + SNP, 
				"l2" = sbpc ~ obstime + SNP + SEX), 
			formLongRandom = list("l1" = hbac ~ 1+ obstime | PATIENT , 
				"l2" = sbpc ~ 1+ obstime | PATIENT ),
			data=datalong.boot, timeVar="obstime")

		stage2datab<-prepare_for_stage2(stage1fitb,datasurv.boot)
		stage2fitb <- coxph(formula = Surv(start, end, event_all) ~ SNP_k1 +  Fitted_l1k1 +DURATION_Years_k1 + SNP_k2 + Fitted_l1k2 +  Fitted_l2k2 +  DURATION_Years_k2 + strata(event_type) + frailty(PATIENT,"gamma"), 
			data = stage2datab, control = coxph.control(timefix = FALSE))
		return(c(stage1fitb$coefs.beta[,"Value"], coef(stage2fitb)))
	})))))}


# JM_res()
# input: outputs from stage1 & stage2, datalong, datasurv, no. of bootstraps 
# extract & format the results, calculate the bootsrap variance-covariance matrix & Wald.1df test of single-parameter tests using bootstrap se
# output: list of coefficients from the longitudinal (stage 1) & time-to-event sub-models (stage 2), bootstrap SE, and Wald 1df P-values, bootstrap variance-covariance matrix

JM_res<-function(stage1fit, stage2fit, datalong_rep, datasurv_rep, nboots){ 
	JMbcov<-JM_bcov(datalong_rep, datasurv_rep, nboots)
	beta<-stage1fit$thetaLong.new$beta
	lbootse<-sqrt(diag(JMbcov)[1:length(beta)])
	lwald1f.p=pchisq((beta/lbootse)^2,1,lower.tail=F)
	D<-stage1fit$thetaLong.new$D
	sigma2<-stage1fit$thetaLong.new$sigma2
	survcoef<-coef(stage2fit)
	sbootse<-sqrt(diag(JMbcov))[names(survcoef)]
	swald1f.p<-pchisq((survcoef/sbootse)^2,1,lower.tail=F)
	return(list(stage1=list(coef=cbind(beta=beta,lbootse=lbootse,lwald1f.p=lwald1f.p), D=D, sigma2=sigma2),
		stage2=cbind(survcoef=survcoef,sbootse=sbootse,swald1f.p=swald1f.p), vcov=JMbcov))
}

