
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
  vector[N] mupercept;
  
  
  expectation[1] = 0.5;
  
  
  
  for (i in 2:N){
  
    alpha[i] = inv_logit(b0a+b1a+percept[i-1]);
    expectation[i] = expectation[i-1]+alpha[i]*(outcome[i]-expectation[i-1]);
    mupercept[i] = inv_logit(b0p+b1p*expectation[i]+b2p*stim[i]);
    
  }
  
  
  
  
  
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  
    target += beta_proportion_lpdf(percept | mupercept,kappa);
    
    for(i in 2:N){
      target += bernoulli_lpmf(predict[i] | expectation[i]);
    }
}

