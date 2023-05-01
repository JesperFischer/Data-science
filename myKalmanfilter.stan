data {
  
  int<lower=1> nSubject; // number of subjects
  int<lower=1> nTrials; // number of subjects
  int<lower=0> sigmaEpsilon;
  matrix[nTrials, nSubject] painRating; // observations
  matrix[nTrials, nSubject] expectPain; // observations
  matrix[nTrials, nSubject] noxInput; // observations
  
}

parameters {
  vector<lower=0>[nSubject] sigmaEta;
  vector<lower=0>[nSubject] sigmaPsi;
  vector<lower=0>[nSubject] painError;
  vector<lower=0>[nSubject] expectError;
  
  // Group-level parameters
  real <lower=0> expectScale;
  real <lower=0> painScale;
  real <lower=0> etaScale;
  real <lower=0> psiScale;
}

transformed parameters {
  matrix[nTrials, nSubject] painMu; 
  matrix[nTrials, nSubject] painVar; 
  matrix[nTrials, nSubject] expectMu; 
  matrix[nTrials, nSubject] expectVar; 


  for (s in 1:nSubject) {
      expectMu[1,s] = 0.5;
      expectVar[1,s] = sigmaEta[s];
    for (t in 1:(nTrials)){
      painMu[t,s] =  ( sigmaEpsilon * expectMu[t,s] + (sigmaPsi[s] + expectVar[t,s]) * noxInput[t,s] ) / 
                       (sigmaEpsilon + sigmaPsi[s] + expectVar[t,s]);
                       
      painVar[t,s] = ( sigmaEpsilon * (sigmaPsi[s] + expectVar[t,s]) ) / (sigmaEpsilon + sigmaPsi[s] + expectVar[t,s]) ;    
      
      expectMu[t,s] = ((sigmaEpsilon + sigmaPsi[s]) * expectMu[t,s] + (expectVar[t,s] * noxInput[t,s])) / 
                                       (sigmaEpsilon + sigmaPsi[s] + expectVar[t,s]) ;
                                       
      expectVar[t,s] = ((sigmaEpsilon + sigmaPsi[s]) * expectVar[t,s] / (sigmaEpsilon + sigmaPsi[s] + expectVar[t,s])) + sigmaEta[s];
    }
  }
}

model {
   for (s in 1:nSubject){
    
    for (t in 1:nTrials){
      target += normal_lpdf(painRating[t,s] | painMu[t,s], painError[s]);
      target += normal_lpdf(expectPain[t,s] | expectMu[t,s], expectError[s]);
    }
    
    target += lognormal_lpdf(sigmaEta[s] | 0, etaScale);
    target += lognormal_lpdf(sigmaPsi[s] | 0, psiScale);
  
    
    
    target += lognormal_lpdf(painError[s] | 0, painScale);
    target += lognormal_lpdf(expectError[s] | 0, expectScale);
  }
  
  // Hierarchical Priors
  
  target += normal_lpdf(etaScale | 0 , 5)-normal_lccdf(0 | 0, 5);
  target += normal_lpdf(psiScale | 0 , 5)-normal_lccdf(0 | 0, 5);
  
  target += normal_lpdf(painScale | 0 , 5)-normal_lccdf(0 | 0, 5);
  target += normal_lpdf(expectScale | 0 , 5)-normal_lccdf(0 | 0, 5);
}
