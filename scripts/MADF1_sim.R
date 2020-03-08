##############################################################
#
# Script to run simulation with MADF1 model
#
# Optimized to parallel run (8 cores)
##############################################################


source(file="scenario_generation.R")



#### model:

model= "MADF1"

require("rstan")
require(parallel)

# non centered parametrization, it is good
sm_meta <- stan_model(file="MADF1.stan", model_name='sm-meta1',verbose=FALSE)

SED4=1478963

# for stan
nchains=4

estimfunc <- function(tr, datagen){
  data_ma <- data.frame(datagen[[tr]]$res)
  
  
  # using all data
  K= max(data_ma$trial)
  N = dim(data_ma)[1]
  studyIdx = data_ma$trial
  doselevel_tot = sort(unique(data_ma$p_true))
  I = length(doselevel_tot)
  doselevel=match(data_ma$p_true, doselevel_tot)
  events = data_ma$nDLT
  total = data_ma$nPT
  
  index_res <- match(doselevel_tot, p_true)
  
  stanData <- list(K        = K,
                   N       = N,
                   I       = I,
                   studyIdx = studyIdx,
                   doselevel     = doselevel,
                   events  = events,
                   total   = total,
                   doses_modif = doses_modif[index_res])
  
  set.seed(SED4+tr)
  
  fit_mod <- sampling(sm_meta, data=stanData ,iter=4000, chains=nchains, control = list(adapt_delta = 0.95))
  results <- extract(fit_mod, pars="ptox")$ptox
  
  # checking if we used or not all doses
  ptox_res <- ptox_res2 <- q1_res <- q3_res <- rep(NA, length(p_true))
  
  
  ptox_res[index_res] <- apply(results,2,median)
  ptox_res2[index_res] <- apply(results,2,mean)
  q1_res[index_res]<- apply(results,2,quantile, probs=0.25)
  q3_res[index_res] <- apply(results,2,quantile, probs=0.75)
  
  return(list(ptox_res=ptox_res, ptox_res2=ptox_res2,
              q1_res=q1_res, q3_res= q3_res))
}


results_parallel <- mclapply(1:TR, estimfunc, datagen=datagen,mc.cores=8)

#### combining results

median_sim <- NULL    # median
mean_sim <- NULL   # mean
CI_lower <- NULL
CI_upper <- NULL

for (tr in 1:TR){
  median_sim <- rbind(median_sim, results_parallel[[tr]]$ptox_res)
  mean_sim <- rbind(mean_sim, results_parallel[[tr]]$ptox_res2)
  CI_lower <- rbind(CI_lower, results_parallel[[tr]]$q1_res)
  CI_upper <- rbind(CI_upper, results_parallel[[tr]]$q3_res)
}

write.table(median_sim, file=paste("res_median_",model,"scen",scen,"_",NTcrm+NT3p3,".txt", sep="") )
write.table(mean_sim, file=paste("res_mean_",model,"scen",scen,"_",NTcrm+NT3p3,".txt", sep="") )
write.table(CI_lower, file=paste("res_q1_",model,"scen",scen,"_",NTcrm+NT3p3,".txt", sep="") )
write.table(CI_upper, file=paste("res_q3_",model,"scen",scen,"_",NTcrm+NT3p3,".txt", sep="") )






