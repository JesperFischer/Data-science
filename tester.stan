data {

  int<lower=1> nSubject; // number of subjects
  int<lower=1> nTrials; // number of subjects
  matrix[nTrials, nSubject] painRating; // observations
  matrix[nTrials, nSubject] expectPain; // prediciton
  matrix[nTrials, nSubject] pain; // prediciton
  
  matrix[nTrials, nSubject] noxInput; // observations
  matrix[nTrials, nSubject] cues; // observations
  
}

parameters {
  vector<lower=0,upper=1>[nSubject] alphaCoefI;
  vector<lower=0>[nSubject] painError;
  vector<lower=0>[nSubject] beta;
  
  vector<lower=0, upper=1>[nSubject] yParam;

  // Group-level parameters
  real <lower=0> alphaI;
  real <lower=0, upper  = 1> betaI;
  real <lower=0, upper = 1> betaY;
  real <lower=0> alphaY;
 // Group-level parameters
  real <lower=0> betaScale;
  real <lower=0> painScale;
}

transformed parameters {
  matrix <lower=0, upper  = 1> [nTrials, nSubject] painMu; 
  matrix <lower=0, upper  = 1> [nTrials+1, nSubject] association; 
  matrix <lower=0, upper  = 1> [nTrials+1, nSubject] expectMu;
  matrix[nTrials, nSubject] predErr;
  

  for (s in 1:nSubject) {
    association[1, s] = 0.5;
      
    for (t in 1:nTrials){
      
      
      
      
      
      
      
     if(yParam[s]*noxInput[t,s]+(1-yParam[s])*expectMu[t,s] == 0){
        painMu[t,s] = 0.01;
      }else if (yParam[s]*noxInput[t,s]+(1-yParam[s])*expectMu[t,s] == 1){
        painMu[t,s] = 0.99;
      }else{
        painMu[t,s] = yParam[s]*noxInput[t,s]+(1-yParam[s])*expectMu[t,s];
        }
        
        
        
      
      
      if(cues[t,s] == 1){
        predErr[t,s] = (painMu[t,s] - expectMu[t,s]);
      }else{
        predErr[t,s] = -(painMu[t,s] - expectMu[t,s]));
      }
      
      
      association[t+1,s] = association[t,s] + alphaCoefI[s] * predErr[t,s];
    
    
    
    }
  }
  
  

}
model {
  for (s in 1:nSubject){

    for (t in 1:nTrials){
      target += beta_proportion_lpdf(painRating[t,s] | painMu[t,s], 1/painError[s]);
      
      target += bernoulli_lpmf(pain[t,s] | (painMu[t,s]^beta[s])/((painMu[t,s]^beta[s])+(1-painMu[t,s])^(beta[s])));

      target += bernoulli_lpmf(expectPain[t,s] |  (expectMu[t,s]^beta[s])/((expectMu[t,s]^beta[s])+(1-expectMu[t,s])^(beta[s])));
      
    }
    
    target += beta_proportion_lpdf(alphaCoefI[s] | betaI , alphaI);
    target += beta_proportion_lpdf(yParam[s] | betaY , alphaY);

    target += lognormal_lpdf(painError[s] | 0, painScale);
    target += lognormal_lpdf(beta[s] | 0, betaScale);
    
  }
  
  // Hierarchical Priors
  target += normal_lpdf(alphaI | 10 , 10)-normal_lccdf(0 | 10, 10); 
  target += beta_proportion_lpdf(betaI | 0.1 , 10) ; 
  target += beta_proportion_lpdf(betaY | 0.1 , 10) ; 
  target += normal_lpdf(alphaY | 10 , 10)-normal_lccdf(0 | 10, 10);
  
  target += normal_lpdf(painScale | 0 , 5)-normal_lccdf(0 | 0, 5);
  target += normal_lpdf(betaScale | 0 , 5)-normal_lccdf(0 | 0, 5);
}

generated quantities{
  real prior_painScale;
  real prior_expectScale;
  real prior_betaI;
  real prior_betaY;
  real prior_alphaCoefI[nSubject];
  real prior_painError[nSubject];
  real prior_expectError[nSubject];
  real prior_yParam[nSubject];
  
  for (s in 1:nSubject){
    prior_alphaCoefI[s] = beta_proportion_rng(betaI , alphaI);
    prior_yParam[s] = beta_proportion_rng(betaY , alphaY);
    
    
    prior_painError[s] = lognormal_rng(-2, painScale);
    prior_expectError[s] = lognormal_rng(-2, expectScale);
  }
  
  prior_betaI = beta_proportion_rng(0.1,10);
  prior_betaY = beta_proportion_rng(0.1,10);
  
  prior_painScale = lognormal_rng(0,0.3);
  prior_expectScale = lognormal_rng(0,0.3);
  
  
  
}

