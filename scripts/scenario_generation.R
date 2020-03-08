##############################################################################
#
# Script to generate scenario datasets
#
# In order to change scenario specifications:
#
# scen          integer number to identify the scenario
# p_true        scenario specific true probability of toxicity
# NTcrm         how many crm dose-finding trials
# NT3p3         how many 3+3 trials
# Sigma         variance-covariance matrix for random effects
#
###############################################################################


library(UBCRM)
library(MASS)
library(dfcrm)

scen= 100
p_true <<- c(.15, .20, .33, .45, .55, .60, .65)
doses_true <<- c(100, 200, 300, 400, 600, 800, 1000)    #in mg for example
doses_modif <<- doses_true/mean(doses_true)


prior_all <- c(0.01,0.05,0.1,0.2,0.33,0.40,0.45)
ndoses <- 3:7
ndoses3p3 <- 5:7
ncohort <- 2:3
npat <- 18:24

SED= 5467
SED2=34
SED3=963299

TR=1000
NTcrm=5
NT3p3=5
target=0.33

sigma =0.3
l=1

sigmas = rep(sigma,length(p_true))
Sigma = matrix(0, ncol=length(p_true), nrow=length(p_true))
for (j in 1:(length(p_true)-1)){
  for (k in (j+1):length(p_true)){
    Sigma[j,k] = Sigma[k,j] = exp(-abs(doses_modif[k]-doses_modif[j])/l)*sigmas[j]*sigmas[k]      #OU formula
  }
}
Sigma = Sigma + diag(sigmas^2)

datagen <- NULL

for (tr in 1:TR){
  set.seed(SED+tr)
  res=NULL
  MTD=NULL
  
  # CRM trials
  for (i in 1:NTcrm){
    cohort = sample(ncohort,1)         # sampling number of patients per cohort
    np = sample(npat,1)                # sampling total number of patients
    np = np - np%%cohort               # adjusting for the cohort selected
    
    ndos = sample(ndoses,1)
    doses = c(sample(p_true[-3],ndos-1), p_true[3])  # deciding the doses for the trial - including the mtd
    
    ptox = ptox_true = doses[order(doses)]         # ordering the ptox

    prior = prior_all[match(ptox_true, p_true)]
    
    # between trial heterogeniety computation for all doses and then selected the trial ones
    mu = qnorm(p_true)

    ptox_trial = mvrnorm(n = 1, mu, Sigma, tol = 1e-6)
    ptox_trial= pnorm(ptox_trial)
    ptox_trial = ptox_trial[match(ptox_true, p_true)]
    
    trial = crmsim(PI=ptox_trial, prior, target, n=np, x0=1, nsim=1, mcohort=cohort, seed=SED2+tr+i)
    
    
    tab= table(trial$level, trial$tox)
    index= as.numeric(dimnames(tab)[[1]])
    
    # check if we have toxicities
    if (dim(tab)[2]==2) {
      tox_print = tab[,2]
    } else tox_print = rep(0,max(index))
    
    res= rbind(res, cbind(ptox[index],ptox_trial[index], tox_print, 
                          apply(tab,1,sum), rep(i,max(index)))  )
    
    MTD = c(MTD, ptox[trial$MTD])
  }
  
  # 3+3 trials
  for (i in 1:NT3p3){
    
    ndos = sample(ndoses,1)
    doses = c(sample(p_true[-3],ndos-1), p_true[3])  # deciding the doses for the trial - including the mtd
    
    ptox = ptox_true = doses[order(doses)]         # ordering the ptox

    prior = prior_all[match(ptox_true, p_true)]
    
    # between trial heterogeniety computation for all doses and then selected the trial ones
    mu = qnorm(p_true)
    
    ptox_trial = mvrnorm(n = 1, mu, Sigma, tol = 1e-6)
    ptox_trial= pnorm(ptox_trial)
    ptox_trial = ptox_trial[match(ptox_true, p_true)]
    
    trial = sim3p3(ptox_trial, seed=SED3+tr+i)
    
    
    index= which(trial$data$npt>0)
    
    res= rbind(res, cbind(ptox[index],ptox_trial[index], trial$data$ndlt[index], 
                          trial$data$npt[index], rep(i+NTcrm,max(index)))  )
    
    MTD = c(MTD, ptox[trial$MTD])
  }
  
 dimnames(res)[[2]] <- c("p_true", "p_trial", "nDLT", "nPT", "trial" )
 datagen[[tr]] <- list(res=res, MTD=MTD) 
}

save.image()




