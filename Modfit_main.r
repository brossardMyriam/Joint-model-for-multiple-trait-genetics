## Without bootstrap (Parameters estimation only)
args=(commandArgs(TRUE))
print(args)
for(k in 1:length(args)){
  eval(parse(text=args[[k]]))
}

message(args)
ipak <- function(pkg){ sapply(pkg, require,character.only = TRUE) }
# "Rmpi","doMPI"
packages <- c("snow", "iterators", "foreach", "doSNOW", "rlecuyer", 
  "survival", "nlme", "aod", "taRifx", "numDeriv")
ipak(packages)

# for parallelization
#cl <- makeCluster(ncores, type="MPI", outfile="JM")
#registerDoSNOW(cl) 
#clusterSetupRNG(cl, seed=rep(12345))
#print(ncores)

R <- 1000 ; nboost<- 500
nseed<- R
nseed <- nboost*R
vectseed <- seq(1:nseed)

system.time( { results <- foreach( rep=srep:erep, .combine='rbind',.verbose=T, .noexport=c('rep', 'datalong', 
	'datasurv', 'subjID', 'bootcur' , 'boot.id', 'datalong.boot','datasurv.boot',	'test','stage2data','outstage2',
	'lmmfr.i','data_k1','data_k2','data_K','data_K.long','JX','data_k1.long',  'data_k2.long','datasurv','datalong',
	'snplist','seed.init','coef.surv','coef.long','tmp.res'),
	.packages = packages, .multicombine=T, .errorhandling = "remove") %dopar% {
   
	  results <- NULL
    ipak <- function(pkg){ sapply(pkg, require, character.only = TRUE) }
    packages <- c("survival","nlme", "aod","taRifx", "numDeriv")
    ipak(packages)
      
    load("/home/bulllab/mbrossard/Joint_models/Example_Rscripts_paper/simulation/simulated_data_reps5.Rdata")    
    source("/home/bulllab/mbrossard/Joint_models/Example_Rscripts_paper/analysis/modfit_subfunctions.r") # 
    longfile <- datalong[which(datalong$rep==rep),]
    survfile <- datasurv[which(datasurv$rep==rep),]
    result <- modfit_datasim(longfile,survfile)
    rm(list=c("data_k2.long","data_k1.long","longfile","survfile",
        "longQT.SNP","survalllongSNP","datalong","stage2data", "data_k1","data_k1"))
    return(result)
  }
} )

colnames( results ) <- c("rep", "SNP1A", "SNP1R", "SNP2A", "SNP2R", "SNP3A", "SNP3R", "SNP4A", "SNP4R", "SNP5A", "SNP5R")
filename <- paste("simulation.result_allSNPs",srep,"_",erep,".RData",sep="")
save( results, file=filename )
warnings()
#stopCluster(cl)

  
