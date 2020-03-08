##############################################################
#
# Script to run simulation with MADF model
#
# Optimized to parallel run (8 cores)
##############################################################


source(file="scenario_generation.R")



#### model:

model= "MADF"

require("rstan")
require(parallel)
require(Iso)

# non centered parametrization, it is good
sm_meta <- stan_model(file="MADF.stan", model_name='sm-meta1',verbose=FALSE)

priors_matrix = data.frame(rbind(
c(-2.0,5.0,0.667,0.5),
c(-4.0,3.5,0.642,0.5)
))

colnames(priors_matrix) <-c("mu_norm", "sigma_norm", "slope_gamma", "var_gamma")

unit_dose= doses_modif[2]-doses_modif[1] # defyining the measure of the jump
theta_MTD = 0.33

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
  
  # dose jumps
  index_res <- match(doselevel_tot, p_true)
  doses_used <- doses_modif[index_res]
  delta_dose <- (doses_used[2:length(doses_used)] - doses_used[1:(length(doses_used)-1)]) / unit_dose
  
  
  nDLT_meta <- function(p_true, data_ma) {
      sum(data_ma[data_ma$p_true == p_true, "nDLT"])
  }
  nPT_meta <- function(p_true, data_ma) {
      sum(data_ma[data_ma$p_true == p_true, "nPT"])
  }
  ti <- round(sapply(doselevel_tot, nDLT_meta, data_ma = data_ma), 0)
  ni <- round(sapply(doselevel_tot, nPT_meta, data_ma = data_ma), 0)
  
  EB_prior= pava(ti/ni)
  pos_EB = order(abs(EB_prior-theta_MTD))[1]
  pos_prior = ( ((doses_used[pos_EB] - doses_used[1])/unit_dose) >=2 )*1
  
  
  var_coeff= priors_matrix[pos_prior+1, "var_gamma"]
  slope_gamma = priors_matrix[pos_prior+1, "slope_gamma"]
  mu_first = priors_matrix[pos_prior+1, "mu_norm"]
  sd_first = priors_matrix[pos_prior+1, "sigma_norm"]
  
  shape_gamma = 1/var_coeff^2
  rate_gamma = (var_coeff^2*slope_gamma)^-1
  
  stanData <- list(K        = K,
                   N       = N,
                   I       = I,
                   studyIdx = studyIdx,
                   doselevel     = doselevel,
                   events  = events,
                   total   = total,
                   doses_modif = doses_modif[index_res],
                   shape_gamma  = shape_gamma,
                   rate_gamma  = rate_gamma,
                   delta_dose = delta_dose,
                   mu_first = mu_first,
                   sd_first = sd_first)
  
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






