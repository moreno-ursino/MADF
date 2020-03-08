// model for efficacy
data {
  int<lower=0> N; // number of data points
  int<lower=0> K; // number of studies
  int<lower=0> I; // number of doses
  int studyIdx[N]; // index of studies 1 < ... < k
  int doselevel[N]; // log of dose assigned to each data point
  int events[N]; // number of tox per each data point (y)
  int total[N]; // total number of patients accrued at each dose (data point)
  real<lower=0> shape_gamma; // mean prior gamma
  real<lower=0> rate_gamma;
  real delta_dose[I-1]; // deltas for the gamma process
  real mu_first;
  real<lower=0> sd_first;
}
parameters {
  real mu0b;
  real<lower=0> mub[I-1];
  real<lower=0> sigma;
  real beta_norm[N]; // study random effects, no shared random effects
}
transformed parameters{
  real mu[I];
  real ptox[I];
  real beta[N];
  
  mu[1] = mu0b;
  ptox[1] = inv_logit(mu[1]);
  for (i in 2:I) {
    mu[i] = mu[i-1] + mub[i-1]; // add the sum here if you want to add constraints
    ptox[i] = inv_logit(mu[i]);
  }
  // random effects
  for (i in 1:N){
    beta[i] = sigma*beta_norm[i];
  }
}
model {
  real logitp[N]; // logit of probs
  
  for (j in 1:N){
    logitp[j] = mu[doselevel[j]] + beta[j];
    //logitp[j] = mu[doselevel[j]] + beta[studyIdx[j]];
    events[j] ~ binomial_logit(total[j], logitp[j]);     // in logit parametrization  
  }

// priors
  mu0b ~ normal(mu_first, sd_first);
  // other jumps:
    for (s in 1:(I-1)){
  mub[s] ~ gamma(shape_gamma*delta_dose[s], rate_gamma);
  }
  sigma ~ normal(0, 1); // actually it is a truncated one - see the lower bound
  beta_norm ~ normal(0, 1);
}
