
data {
  int<lower=0> N;
  vector[N] percept;
  vector[N] stim;
  int predict[N];
  vector[N] outcome;
  
  
  
}

parameters {
  real<lower=0> kappa;
  real b0a;
  real b1a;
  real b0p;
  real b1p;
  real b2p;
  
}

transformed parameters{
  
  vector[N] expectation;
  vector[N] alpha;
  vector[N] alphap;
  vector[N] mupercept;
  
  
  expectation[1] = 0.5;
  mupercept[1] = 0.01;
  
  
  for (i in 2:N){
    
    alpha[i] = inv_logit(b0a+b1a*(percept[i-1]*(1-percept[i-1])));
    expectation[i] = expectation[i-1]+alpha[i]*(outcome[i-1]-expectation[i-1]);
    
    alphap[i] = inv_logit(b0p+b1p*(expectation[i]*(1-expectation[i])));
    
    mupercept[i] = mupercept[i-1]+alphap[i]*(b2p*stim[i]-mupercept[i-1]);
    
  }
  
  
  
  
  
  
}

model {
    target += normal_lpdf(b0a | 0, 2);
    target += normal_lpdf(b1a | 0, 2);
    target += normal_lpdf(b0p | 0, 2);
    target += normal_lpdf(b1p | 0, 2);
    target += normal_lpdf(b2p | 0, 2);
    target += lognormal_lpdf(kappa |3.5,0.5);
    
    
    
  
    target += beta_proportion_lpdf(percept | mupercept,kappa);
    
    for(i in 2:N){
      target += bernoulli_lpmf(predict[i] | expectation[i]);
    }
}


generated quantities{
  
    real prior_b0a = normal_rng(0, 2);
    real prior_b1a = normal_rng(0, 2);
    
    real prior_b0p = normal_rng(0, 2);
    real prior_b1p = normal_rng(0, 2);
    
    real prior_b2p = normal_rng(0, 2);
    real prior_b3p = normal_rng(0, 2);
    
    real prior_kappa = lognormal_rng(3.5,0.5);
    
}

