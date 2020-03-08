// model for efficacy
data {
  int<lower=0> N; // number of data points
  int<lower=0> K; // number of studies
  int<lower=0> I; // number of doses
  int studyIdx[N]; // index of studies 1 < ... < k
  int doselevel[N]; // log of dose assigned to each data point
  int events[N]; // number of tox per each data point (y)
  int total[N]; // total number of patients accrued at each dose (data point)
  real doses_modif[I];
}
parameters {
  real mu0b;
  real<lower=0> mub[I-1];
  real<lower=0> sigma;
  real<lower=0.001> l;
  vector[I] beta_norm[K]; // study random effects, one for each dose
}
transformed parameters{
  real mu[I];
  real ptox[I];
  vector[I] beta[K];
  matrix[I,I] Sigma_matrix;
  matrix[I,I] L;
    
  mu[1] = mu0b;
  ptox[1] = inv_logit(mu[1]);
  for (i in 2:I) {
    mu[i] = mu[i-1] + mub[i-1]; // add the sum here if you want to add constraints
    ptox[i] = inv_logit(mu[i]);
  }
  
  // random effects
  // Sigma
    
  for (j in 1:(I-1)){
    for (k in (j+1):I){
      Sigma_matrix[j,k] =  exp(-fabs(doses_modif[k]-doses_modif[j])/l)*sigma^2;
      Sigma_matrix[k,j] = Sigma_matrix[j,k];
    }
  }
  for (j in 1:I){
      Sigma_matrix[j,j] =  sigma^2;
  }
  
  // cholesky 
  L = cholesky_decompose(Sigma_matrix);
  
  for (i in 1:K){
    beta[i] = L * beta_norm[i];
  }
}
model {
  real logitp[N]; // logit of probs
  
  
  for (j in 1:N){
    //logitp[j] = mu[doselevel[j]] + beta[j];
    logitp[j] = mu[doselevel[j]] + beta[ studyIdx[j],doselevel[j] ];
    events[j] ~ binomial_logit(total[j], logitp[j]);     // in logit parametrization  
  }

// priors
  mu0b ~ normal(-2, 10);
  mub ~ gamma(3, 0.5);
  sigma ~ normal(0, 1); // actually it is a truncated one - see the lower bound
  l ~ inv_gamma(1,1);
  for (i in 1:K) beta_norm[i] ~ normal(0, 1);
}
