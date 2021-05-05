## R functions for MV LME
## Function from JoineRML package

mvlme <- function (formLongFixed, formLongRandom, data, timeVar,inits=NULL, 
tol.em = 1e-4 , verbose = FALSE )  {
time.start <- Sys.time()
balanced <- FALSE
  if (!is.list(formLongFixed)) {
        balanced <- TRUE
        formLongFixed <- list(formLongFixed)
        formLongRandom <- list(formLongRandom)
        L <- 1
    } else {
        L <- length(formLongFixed)
    }
    if (!("list" %in% class(data))) {
        balanced <- TRUE
        data <- list(data)
        if (L > 1) {
            for (l in 2:L) {
                data[[l]] <- data[[1]]
            }
        }
    }     else {
        balanced <- (length(unique(data)) == 1)
    }
    
	if (length(data) != L) {
        stop(paste("The number of datasets expected is L =", L))
    }
	
    data <- lapply(data, function(d) {
        if (any(c("tbl_df", "tbl") %in% class(d))) {
            return(as.data.frame(d))
        }
        else {
            return(d)
        }
    })
    id <- as.character(nlme::splitFormula(formLongRandom[[1]], "|")[[2]])[2]
    n <- length(unique(data[[1]][, id]))
    if (length(timeVar) == 1 & (L > 1)) {
        timeVar <- rep(timeVar, L)
    }
	
	
    if (length(timeVar) != L) {
        stop(paste("The length of timeVar must equal", L))
    }
    for (l in 1:L) {
        data[[l]] <- data[[l]][order(xtfrm(data[[l]][, id]), 
            data[[l]][, timeVar[l]]), ]
        data[[l]][, id] <- as.factor(data[[l]][, id])
        data[[l]] <- droplevels(data[[l]])
    }
    if (L > 1) {
        uniq.ids <- list(sort(unique(data[[1]][, id])))
        for (l in 2:L) {
            uniq.ids[[l]] <- sort(unique(data[[l]][, id]))
            if (!identical(uniq.ids[[l - 1]], uniq.ids[[l]])) {
                stop("Every subject must have at least one measurement per each outcome")
            }
        }
    }
	
  	lfit <- list()
    mf.fixed <- list()
    yik <- list()
    Xik <- list()
    nk <- vector(length = L)
    Xik.list <- list()
    nik.list <- list()
    Zik <- list()
    Zik.list <- list()
	
	# Preparation of design matrix 
    for (l in 1:L) {
        lfit[[l]] <- nlme::lme(fixed = formLongFixed[[l]], random = formLongRandom[[l]], 
            data = data[[l]], method = "ML", control = nlme::lmeControl(opt = "optim"))
        lfit[[l]]$call$fixed <- eval(lfit[[l]]$call$fixed)
        mf.fixed[[l]] <- model.frame(lfit[[l]]$terms, data[[l]][,all.vars(formLongFixed[[l]])])
        yik[[l]] <- by(model.response(mf.fixed[[l]], "numeric"), data[[l]][, id], as.vector)
        Xik[[l]] <- data.frame(id = data[[l]][, id], model.matrix(formLongFixed[[l]], data[[l]]))
        nk[l] <- nrow(Xik[[l]])
        Xik.list[[l]] <- by(Xik[[l]], Xik[[l]]$id, function(u) { as.matrix(u[, -1])})
        nik.list[[l]] <- by(Xik[[l]], Xik[[l]]$id, nrow)
        ffk <- nlme::splitFormula(formLongRandom[[l]], "|")[[1]]
        Zik[[l]] <- data.frame(id = data[[l]][, id], model.matrix(ffk, data[[l]])) # Z design matrix 
        Zik.list[[l]] <- by(Zik[[l]], Zik[[l]]$id, function(u) { as.matrix(u[, -1]) }) # Z design matrix split in list of matrix per individual
    }
	
	# Yi = Per individual list of phenotypes values concatenated in one vector
    yi <- sapply(names(yik[[1]]), function(i) {
        unlist(lapply(yik, "[[", i))
    }, USE.NAMES = TRUE, simplify = FALSE)
    
	# Xi = list of X design matrix of fixed effects for all L traits & per individual
	Xi <- sapply(names(Xik.list[[1]]), function(i) {
        as.matrix(Matrix::bdiag(lapply(Xik.list, "[[", i)))
    }, USE.NAMES = TRUE, simplify = FALSE)
    
	# Zi = list of X design matrix of random effects for all L traits & per individual
	Zi <- sapply(names(Zik.list[[1]]), function(i) {
        as.matrix(Matrix::bdiag(lapply(Zik.list, "[[", i)))
    }, USE.NAMES = TRUE, simplify = FALSE)
	# Zi = list of transposed Zi  design matrix & per individual
	Zit <- lapply(Zi, t)
    # XtXi = equivalent to (but usually slightly faster than) the call t(x) %*% y (crossprod) or x %*% t(y) (tcrossprod). 
	# XtXi = list of XtX per sample
	XtXi <- lapply(Xi, crossprod)
	# Reduce sum all XtXi matrixes
	XtX <- Reduce("+", XtXi)
    XtX.inv <- solve(XtX)
    
	# List of Xtyi per individual
	Xtyi <- mapply(function(x, y) {
        crossprod(x, y)
    }, x = Xi, y = yi, SIMPLIFY = FALSE)
	
	# List of XtZi per individual
    XtZi <- mapply(function(x, z) {
        crossprod(x, z)
    }, x = Xi, z = Zi, SIMPLIFY = FALSE)
    
	nik <- sapply(names(nik.list[[1]]), function(i) {
        unlist(lapply(nik.list, "[[", i))
    }, USE.NAMES = TRUE, simplify = FALSE)
    
	p <- sapply(1:L, function(i) { ncol(Xik[[i]]) - 1 })
	r <- sapply(1:L, function(i) { ncol(Zik[[i]]) - 1 })
 
	list1 <- list(yi = yi, Xi = Xi, Zi = Zi, Zit = Zit, nik = nik, 
        yik = yik, Xik.list = Xik.list, Zik.list = Zik.list, 
        XtX.inv = XtX.inv, Xtyi = Xtyi, XtZi = XtZi, p = p, r = r, 
        L = L, n = n, nk = nk)

	out <- mvlmefit(lfit = lfit, inits = inits, list1 = list1, 
        z = z, L = L, p = p, r = r, tol.em = 1e-4, verbose = verbose)
	time.end <- Sys.time()
	timeprint <- paste('total time MVLM:', time.end-time.start)
	#message(timeprint)
 return(out)
}	
	
	
mvlmefit <- function (lfit, inits, list1, z, L, p, r, tol.em, verbose) {
	# Extract covariance for random effects for each trait & put the together 
  D <- Matrix::bdiag(lapply(lfit, function(u) matrix(nlme::getVarCov(u), dim(nlme::getVarCov(u)))))
  D <- as.matrix(D)

# Add colnames to D matrix  
  D.names <- c()
    for (l in 1:L) {
        D.names.l <- paste0(rownames(nlme::getVarCov(lfit[[l]])), "_", l)
        D.names <- c(D.names, D.names.l)
    }
 rownames(D) <- colnames(D) <- D.names
 # Extract fixed effects from lfit for each trait & put them together
 beta <- do.call("c", lapply(lfit, fixef))
 names(beta) <- paste0(names(beta), "_", rep(1:L, p))
 # Idem for variance of error terms
 sigma2 <- unlist(lapply(lfit, function(u) u$sigma))^2

# EM iterations for multivariate LMM (EM MV LMMM)
   if ((L > 1) && !all(c("beta", "D", "sigma2") %in% names(inits))) {
        message("Running multivariate LMM EM algorithm...")
        out <- emmvlme(thetaLong = list(beta = beta, D = D, sigma2 = sigma2), 
            list1 = list1, z = z, tol.em = tol.em, verbose = verbose)
		message("Finished multivariate LMM EM algorithm...")
    }   else {
        out <- list(D = D, beta = beta, sigma2 = sigma2)
    }


    if ("beta" %in% names(inits)) {
        if (length(inits$beta) != sum(p)) {
            stop("Dimension of beta inits does not match model.")
        }
        beta <- inits$beta
        names(beta) <- names(out[["beta"]])
        out[["beta"]] <- beta
    }
	
    if ("D" %in% names(inits)) {
        if (nrow(inits$D) != sum(r)) {
            stop("Dimension of D inits does not match model.")
        }
        is.posdef <- all(eigen(inits$D)$values > 0)
        if (is.posdef) {
            D <- inits$D
            rownames(D) <- colnames(D) <- rownames(out[["D"]])
            out[["D"]] <- D
        }
        else {
            warning("Initial parameter matrix D is non positive definite: falling back to automated value")
        }
    }
    if ("sigma2" %in% names(inits)) {
        sigma2 <- inits$sigma2
        if (length(sigma2) != L) {
            stop("Dimension of sigma2 inits does not match model")
        }
        names(sigma2) <- paste0("sigma2_", 1:L)
        out[["sigma2"]] <- sigma2
    }
    return(out)
} 

		
emmvlme <- function (thetaLong, list1, z, tol.em, verbose) {
    yi <- list1$yi
    Xi <- list1$Xi
    XtX.inv <- list1$XtX.inv
    Xtyi <- list1$Xtyi
    XtZi <- list1$XtZi
    Zi <- list1$Zi
    Zit <- list1$Zit
    nik <- list1$nik
    yik <- list1$yik
    Xik.list <- list1$Xik.list
    Zik.list <- list1$Zik.list
    n <- list1$n
    p <- list1$p
    r <- list1$r
    L <- list1$L
    nk <- list1$nk
    delta <- 1
    while (delta > tol.em) {
        D <- thetaLong$D
        beta <- thetaLong$beta
        sigma2 <- thetaLong$sigma2
        # Construction of inverse sigma2 matrix for each individual for L traits 
		# dimension for each sample = sum of no. of measures for L traits
		Sigmai.inv <- lapply(nik, function(i) {
            diag(x = rep(1/sigma2, i), ncol = sum(i)) })
			
		Dinv <- solve(D)
		# Ai = in Eq. 6 (computed for each sample i) variance-covariance matrix 
		Ai <- mapply(FUN = function(zt, s, z) {
            solve((zt %*% s %*% z) + Dinv)
        }, z = Zi, zt = Zit, s = Sigmai.inv, SIMPLIFY = FALSE)
    
		# Mean of bi MVN from Eq. 6
		# Eb = List of random effects per individual
		Eb <- mapply(function(a, z, s, y, X) {
            as.vector(a %*% (z %*% s %*% (y - X %*% beta)))
        }, a = Ai, z = Zit, s = Sigmai.inv, y = yi, X = Xi, SIMPLIFY = FALSE)
		# EbbT = Ai*t(b)
		EbbT <- mapply(function(v, e) {
            v + tcrossprod(e)
        }, v = Ai, e = Eb, SIMPLIFY = FALSE)
		
		# New covariance matrix of RE (for the 2 traits)
		D.new <- Reduce("+", EbbT)/n
		rownames(D.new) <- colnames(D.new) <- rownames(D)
		# rr = Xtyi - (XtZi %*% Eb) = new residuals computed / individual (dim = nrow= n fixed variables, ncol=n samples)
		rr <- mapply(function(x1, x2, b) {
				x1 - (x2 %*% b)
			}, x1 = Xtyi, x2 = XtZi, b = Eb)
   
		rr.sum <- rowSums(rr)
		# new fixed effects , FORMULA FROM SUPPLEMEMTS - CLOSED FORM UPDATE EQUATION
		beta.new <- as.vector(XtX.inv %*% rr.sum)
		names(beta.new) <- names(beta)
		beta.inds <- cumsum(c(0, p))
		b.inds <- cumsum(c(0, r))
		sigma2.new <- vector(length = L)
		
		# Computation of new residuals (for each trait)
		for (l in 1:L) {
			beta.l <- beta.new[(beta.inds[l] + 1):(beta.inds[l + 1])]
			SSq <- mapply(function(y, x, z, b, b2) {
                b.l <- b[(b.inds[l] + 1):(b.inds[l + 1])]
                bbT.l <- b2[(b.inds[l] + 1):(b.inds[l + 1]), (b.inds[l] + 1):(b.inds[l + 1])]
                residFixed <- (y - x %*% beta.l)
                t(residFixed) %*% (residFixed - 2 * (z %*% b.l)) + sum(diag(crossprod(z) %*% bbT.l))
		         }, y = yik[[l]], x = Xik.list[[l]], z = Zik.list[[l]], b = Eb, b2 = EbbT)
				 
            sigma2.new[l] <- sum(SSq)/nk[[l]] # nk =  total number of observation for all samples by trait l
        }
        names(sigma2.new) <- paste0("sigma2_", 1:L)
		#var.beta.new <- sqrt(as.vector(c(rep(sigma2.new[1],3),rep(sigma2.new[2],3))*XtX.inv))

		#thetaLong.new <- list(D = D.new, beta = beta.new, sigma2 = sigma2.new,b=Eb)
		thetaLong.new <- list(D = D.new, beta = beta.new, sigma2 = sigma2.new)
        delta <- sapply(c("D", "beta", "sigma2"), function(i) {
            abs(thetaLong[[i]] - thetaLong.new[[i]])/(abs(thetaLong[[i]]) + 0.001) })
        delta <- max(unlist(delta))
        thetaLong <- thetaLong.new
        if (verbose) {
            print(thetaLong.new)
			#print(list(Eb))
        }
    }
	   Var.Covar.Mat <- sefixed(estim = thetaLong.new, list1 = list1)
     indSNP=NULL
     indSNP=grep("SNP",names(beta)) 
     Std.Err <- sqrt(diag(Var.Covar.Mat))
       
     coefs.beta <- cbind(Value = beta, Std.Err = Std.Err ,
  	 `z-value` = beta/Std.Err,`p-value` = pchisq((beta/Std.Err)^2,1, lower.tail = FALSE))
         
     if (length(indSNP)>0) {
	     wald2df=pchisq(wald.test(Var.Covar.Mat,beta,Terms=indSNP)$result$chi2[1],
			wald.test(Var.Covar.Mat,beta,Terms=indSNP)$result$chi2[2],lower.tail=FALSE)
     
     	   return(list(thetaLong.new=thetaLong.new, coefs.beta = coefs.beta, b = list(Eb),wald2df=wald2df))
      } else {  return(list(thetaLong.new=thetaLong.new, coefs.beta = coefs.beta, b = list(Eb))) }
}

sefixed = function ( estim, list1 ) {	
    yi <- list1$yi
    Xi <- list1$Xi
    Zi <- list1$Zi
    Zit <- list1$Zit
    nik <- list1$nik
    yik <- list1$yik
    Xik.list <- list1$Xik.list
    Zik.list <- list1$Zik.list
    p <- list1$p
    r <- list1$r
    L <- list1$L
    n <- list1$n
    nk <- list1$nk
	
	D <- estim$D
    beta <- estim$beta
    sigma2 <- estim$sigma2
    Sigmai <- lapply(nik, function(i) {
        diag(x = rep(sigma2, i), ncol = sum(i))
    })


	######################################################
	# log density of y
	loglikf=function(betavar=c(b1,b2,b3,b4,b5,b6)) {
	Reduce('+', fy <- mapply(function(y, x, z, zt, s, nik) {
		r <- y - (x %*% betavar)
        v <- s + z %*% D %*% zt
        vinv <- solve(v)
        -0.5 * (sum(nik) * log(2 * pi) + as.numeric(determinant(v,logarithm = TRUE)$modulus) + colSums(t(r) %*% vinv %*% r)) 
    }, y = yi, x = Xi, z = Zi, zt = Zit, s = Sigmai, nik = nik))
	}
	
	return(solve(-hessian(loglikf,beta)))
}

modfit_datasim <- function( longfile, survfile ) {
	stage1fit <- NULL
	res <- interm.res <- NULL
  
	list.SNP<-c( "SNP1A", "SNP1R", "SNP2A", "SNP2R", "SNP3A", 
		"SNP3R", "SNP4A", "SNP4R", "SNP5A", "SNP5R" )

	for(snp in 1:length(list.SNP)) {
		SNP <- list.SNP[snp]
		longfile0 <- longfile
    survfile0 <- survfile
    colnames(longfile0)[grep(SNP,colnames(longfile0))]<-colnames(survfile0)[grep(SNP,colnames(survfile0))]<-"SNP"
	  
	  	tryCatch({
				# Stage 1 of the JM - multivariate mixed model
				stage1fit <- mvlme(formLongFixed = 
					list("l1" = hbac ~ obstime + SNP, "l2" = sbpc ~ obstime + SNP + SEX), 
					formLongRandom = list("l1" = hbac ~ 1+ obstime | PATIENT , "l2" = sbpc ~ 1+ obstime | PATIENT ),
					data=longfile0, timeVar="obstime")
				beta <-stage1fit$coefs.beta[grep("SNP", rownames(stage1fit$coefs.beta)),c(1,2,4)]
				sigma2_l1 <-stage1fit$thetaLong.new$sigma2[1]  
				sigma2_l2<-stage1fit$thetaLong.new$sigma2[2] 
			  tmp <- rep(list(stage1fit$coefs.beta[,1]), length(unique(longfile0$PATIENT)))
				names(tmp) <- unique(longfile0$PATIENT)
				lmmfr.i <- merge(stage1fit$b[[1]],tmp)
        
       	lmfit_U <- lme( sim_Uc  ~ obstime + SNP , random = ~ 1 + obstime | PATIENT, 
					data = longfile0, control = lmeControl(opt = "optim"))
				sigma2_U <-summary(lmfit_U)$sigma^2

				data_k1 <- survfile0[,c("PATIENT", "rep", "event_k1", "tobs_k1", "SEX", "DURATION_Years",
					colnames(survfile0)[grep("SNP",colnames(survfile0))])]
				data_k1$tobs_k1[data_k1$event_k1==1 & data_k1$tobs_k1==0] <- 1e-5
				data_k2 <- survfile0[,c("PATIENT","rep","event_k2","tobs_k2","SEX", "DURATION_Years",
					colnames(survfile0)[grep("SNP",colnames(survfile0))])]
				data_k2$tobs_k2[data_k2$event_k2==1 & data_k2$tobs_k2==0] <- 1e-5
				cut.pointsK <- unique(c(data_k1$tobs_k1[data_k1$event_k1==1],
					data_k2$tobs_k2[data_k2$event_k2==1]))
			  
				data_k1.long <- survSplit(data=data_k1,cut=cut.pointsK, 
          end = "tobs_k1", start = "start", event = "event_k1")
				colnames(data_k1.long)[which(colnames(data_k1.long)=="tobs_k1")] <- "end"
				colnames(data_k1.long)[which(colnames(data_k1.long)=="event_k1")] <- "event_all"
				data_k1.long$event_type <- 0
				
				data_k2.long <- survSplit(data=data_k2,cut=cut.pointsK, 
          end = "tobs_k2", start = "start", event = "event_k2")
				colnames(data_k2.long)[which(colnames(data_k2.long)=="tobs_k2")] <- "end"
				colnames(data_k2.long)[which(colnames(data_k2.long)=="event_k2")] <- "event_all"
				data_k2.long$event_type <- 1
				
				data_K.long <- rbind(data_k1.long,data_k2.long)
				rm(data_k2.long) ; rm(data_k1.long)
				rm(data_k2); rm(data_k1)
			  
				stage2data <- do.call(rbind, by(data_K.long, INDICES=list(data_K.long$PATIENT), 
					function (x) {
						  lme.coef <- lmmfr.i
						  lme.coef <- unlist(lmmfr.i[which(names(lmmfr.i)==unique(x$PATIENT))])
						  lat.rcoef <- lmfit_U$coefficients$random$PATIENT[which(rownames(lmfit_U$coefficients$random$PATIENT)==unique(x$PATIENT)),]
						  Traj_l1 <- function(x,t) { lme.coef[5] + lme.coef[1] + ( lme.coef[6]+lme.coef[2] )*t + lme.coef[7]*x$SNP  }
						  Traj_l2 <- function(x,t) { lme.coef[8] + lme.coef[3] + ( lme.coef[4]+lme.coef[9] )*t + lme.coef[10]*x$SNP + lme.coef[11]*x$SEX }
						  Traj_U <- function(x,t) { lmfit_U$coefficients$fixed[1] + lat.rcoef[1] + ( lmfit_U$coefficients$fixed[2]+lat.rcoef[2] )*t + lmfit_U$coefficients$fixed[3]*x$SNP }
						  e1 <- rnorm(nrow(x),0,sd=sqrt(sigma2_l1))
						  e2 <- rnorm(nrow(x),0,sd=sqrt(sigma2_l2))
						  return( cbind(x,Fitted_l1=Traj_l1(x, x$start),Fitted_l2=Traj_l2(x, x$start),
							Fitted_U=Traj_U(x, x$start),Obs_l1=Traj_l1(x, x$start)+e1,Obs_l2=Traj_l2(x, x$start)+e2) )
					})
			  )
			  rm(longfile0); rm(survfile0)
			  rm(data_K.long)
			  
			  # Coding of variables for time-to-event sub-models for event of type k=1 (retinopathy)
			  stage2data$SNP_k1 <- stage2data$SNP*(1-stage2data$event_type)
			  stage2data$DURATION_Years_k1 <- stage2data$DURATION_Years*(1-stage2data$event_type)
			  stage2data$Fitted_l1k1 <- stage2data$Fitted_l1*(1-stage2data$event_type)
			  stage2data$Fitted_l2k1 <- stage2data$Fitted_l2*(1-stage2data$event_type)
			  stage2data$Fitted_Uk1 <- stage2data$Fitted_U*(1-stage2data$event_type)
			  stage2data$Obs_l1k1 <- stage2data$Obs_l1*(1-stage2data$event_type)
			   
			 # Coding of variables for time-to-event sub-models for event of type k=2  (nephropathy)
			  stage2data$SNP_k2 <- stage2data$SNP*stage2data$event_type
			  stage2data$DURATION_Years_k2 <- stage2data$DURATION_Years*stage2data$event_type
			  stage2data$Fitted_l1k2 <- stage2data$Fitted_l1*stage2data$event_type
			  stage2data$Fitted_l2k2 <- stage2data$Fitted_l2*stage2data$event_type
			  stage2data$Fitted_Uk2 <- stage2data$Fitted_U*stage2data$event_type
			  stage2data$Obs_l1k2 <- stage2data$Obs_l1*stage2data$event_type
			  stage2data$Obs_l2k2 <- stage2data$Obs_l2*stage2data$event_type
			  
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
			 
	  res0 <- res1 <- NULL
	  
		gamma <- indirect_k1 <-indirect_l1k2 <- indirect_l2k2 <- indirect_k2  <- NULL
		total_k1 <- total_k2 <- total_K <- NULL          

	   tryCatch( 
	   { 
		MMfit <- coxph(formula = Surv(start, end, event_all) ~ SNP_k1 +  SNP_k2  + 
				DURATION_Years_k1+ DURATION_Years_k2 +
				strata(event_type) , data = stage2data, control = coxph.control(timefix = FALSE))
		MM <- summary(MMfit)$coefficients[grep("SNP", rownames(summary(MMfit)$coefficients)) ,c(1,3,5)]       
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
		)

		## Stage 2 of JM-cmp 
		gamma <- indirect_k1 <-indirect_l1k2 <- indirect_l2k2 <- indirect_k2  <- NULL
		total_k1 <- total_k2 <- total_K <- NULL  
		tryCatch({ JM_cmpfit <- coxph(formula = Surv(start, end, event_all) 
			~ SNP_k1 +  SNP_k2 +  Fitted_l1k1 + Fitted_l1k2 +  Fitted_l2k2 +
			   DURATION_Years_k1 + Fitted_Uk1 + DURATION_Years_k2 + Fitted_Uk2 +
			   strata(event_type) + frailty(PATIENT,"gamma"), data = stage2data, 
         control = coxph.control(timefix = FALSE)) 
			  JM_cmp_sum <- summary(JM_cmpfit)$coefficients[1:5,c(1,3,6)]       
			  gamma <- JM_cmp_sum[1:2,1]
			  alpha <- JM_cmp_sum[3:5,1]
			  indirect_k1 <- alpha[1]*beta[1,1]
			  indirect_l1k2 <- alpha[2]*beta[1,1]
			  indirect_l2k2 <- alpha[3]*beta[2,1]
			  indirect_k2 <- indirect_l1k2 + indirect_l2k2 
			  indirect_K <- indirect_k1 +indirect_k2
			  total_k1 <- indirect_k1 + gamma[1]
			  total_k2 <- indirect_k2 + gamma[2]
			  totalK <- total_k1+total_k2
			  JM_cmp <- list(JM_cmp_sum,beta=beta,gamma=gamma,alpha=alpha,
					indirect_k1=indirect_k1, indirect_l1k2=indirect_l1k2, indirect_l2k2=indirect_l2k2, 
					indirect_k2=indirect_k2,indirect_K=indirect_K,
					total_k1=total_k1,
					total_k2=total_k2,total_K=total_K)
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

		## Stage 2 of JM-mis 
		gamma <- indirect_k1 <-indirect_l1k2 <- indirect_l2k2 <- indirect_k2  <- NULL
		total_k1 <- total_k2 <- total_K <- NULL 		
		tryCatch({ JM_misfit <- coxph(formula = Surv(start, end, event_all) 
			~  SNP_k1 +  SNP_k2 +  Fitted_l1k1 + Fitted_l1k2 +  Fitted_l2k2 
			  + DURATION_Years_k1  + DURATION_Years_k2  
			  + strata(event_type) + frailty(PATIENT,"gamma"), control = coxph.control(timefix = FALSE), data = stage2data) 
			JM_mis_sum <- summary(JM_misfit)$coefficients[1:5,c(1,3,6)]       
			gamma <- JM_mis_sum[1:2,1]
			alpha <- JM_mis_sum[3:5,1]
			indirect_k1 <- alpha[1]*beta[1,1]
			indirect_l1k2 <- alpha[2]*beta[1,1]
			indirect_l2k2 <- alpha[3]*beta[2,1]
			indirect_k2 <- indirect_l1k2 + indirect_l2k2 
			indirect_K <- indirect_k1 +indirect_k2
			total_k1 <- indirect_k1 + gamma[1]
			total_k2 <- indirect_k2 + gamma[2]
			total_K <- total_k1+total_k2
			JM_mis <- list(JM_mis_sum,beta=beta,gamma=gamma,alpha=alpha,
				indirect_k1=indirect_k1, indirect_l1k2=indirect_l1k2, indirect_l2k2=indirect_l2k2, 
				indirect_k2=indirect_k2,indirect_K=indirect_K,
				total_k1=total_k1,
				total_k2=total_k2,total_K=total_K)
			  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

		## Stage 2 of JM-sep (Cox PH frailty model)
		gamma <- indirect_k1 <-indirect_l1k2 <- indirect_l2k2 <- indirect_k2  <- NULL
		total_k1 <- total_k2 <- total_K <- NULL  
		tryCatch({ JM_sepfit <- coxph(formula = Surv(start, end, event_all) 
			~ SNP_k1 +  SNP_k2 +  Fitted_l1k1 + Fitted_l1k2 +  Fitted_l2k2 
			  + DURATION_Years_k1  + DURATION_Years_k2 , data = stage2data, control = coxph.control(timefix = FALSE))
			JM_sep_sum <- summary(JM_sepfit)$coefficients[1:5,c(1,3,5)]       
			gamma <- JM_sep_sum[1:2,1] 
			 alpha <- JM_sep_sum[3:5,1]
			indirect_k1 <- alpha[1]*beta[1,1]
			indirect_l1k2 <- alpha[2]*beta[1,1]
			indirect_l2k2 <- alpha[3]*beta[2,1]
			indirect_k2 <- indirect_l1k2 + indirect_l2k2 
			indirect_K <- indirect_k1 +indirect_k2
			total_k1 <- indirect_k1 + gamma[1]
			total_k2 <- indirect_k2 + gamma[2]
			totalK <- total_k1+total_k2
			JM_sep <- list(JM_sep_sum,beta=beta,gamma=gamma,alpha=alpha,
				indirect_k1=indirect_k1, indirect_l1k2=indirect_l1k2, indirect_l2k2=indirect_l2k2, 
				indirect_k2=indirect_k2,indirect_K=indirect_K,
				total_k1=total_k1,
				total_k2=total_k2,total_K=total_K)
			  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
			  
		
		## Stage 2 of CM-obs (Cox PH frailty model)
		gamma <- indirect_k1 <-indirect_l1k2 <- indirect_l2k2 <- indirect_k2  <- NULL
		total_k1 <- total_k2 <- total_K <- NULL  
		# Corresponds to the "naive" Cox PH frailty model (CM-obs) fitted in the paper	  
		tryCatch({ CMobs_fit <- coxph(formula = Surv(start, end, event_all) 
			~ SNP_k1 +  SNP_k2 +  Obs_l1k1 + Obs_l1k2 + Obs_l2k2 
			  + DURATION_Years_k1  + DURATION_Years_k2 + strata(event_type) + 
			  frailty(PATIENT,"gamma"), data = stage2data, control = coxph.control(timefix = FALSE))
			CMobs_sum <- summary(CMobs_fit)$coefficients[1:5,c(1,3,6)]       
			gamma <- CMobs_sum[1:2,1]
			alpha <- CMobs_sum[3:5,1]
			indirect_k1 <- alpha[1]*beta[1,1]
			indirect_l1k2 <- alpha[2]*beta[1,1]
			indirect_l2k2 <- alpha[3]*beta[2,1]
			indirect_k2 <- indirect_l1k2 + indirect_l2k2 
			indirect_K <- indirect_k1 +indirect_k2
			total_k1 <- indirect_k1 + gamma[1]
			total_k2 <- indirect_k2 + gamma[2]
			total_K <- total_k1+total_k2
		   CM_obs <- list( CMobs_sum ,beta=beta,gamma=gamma,alpha=alpha,
				indirect_k1=indirect_k1, indirect_l1k2=indirect_l1k2, indirect_l2k2=indirect_l2k2, 
				indirect_k2=indirect_k2,indirect_K=indirect_K,
				total_k1=total_k1,
				total_k2=total_k2,total_K=total_K)
			  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	  
		  
		
		interm.res <- list(
			naive=list(MM=MM, CM_obs=CM_obs), 
			stage1=list(beta),
			stage2=list(JM_cmp=JM_cmp,JM_mis=JM_mis,JM_sep=JM_sep)
			)
		
		res <- cbind(res,list(interm.res))
	  }
		rm(list=c("longfile", "survfile"))
		return(cbind(rep=list(unique(stage2data$rep)),res))
  	rm(stage2data)
}

