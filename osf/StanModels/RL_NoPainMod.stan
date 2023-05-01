data {
  int<lower=1> nCues; // number of cues
  int<lower=1> nSubject; // number of subjects
  int<lower=1> nTrials; // number of subjects
  int<lower=1> uniformPriorParameter;
  matrix[nTrials, nSubject] painRating; // observations
  matrix[nTrials, nSubject] expectPain; // observations
  matrix[nTrials, nSubject] noxInput; // observations
  int<lower=1> myMulti;
  int cueIndexVec[(nTrials-nCues), nSubject]; // indexing next trial for spec cue
  int firstTimeCue[nCues, nSubject];
  matrix[nTrials, nSubject] cueType; // 
  matrix[nCues,  nSubject] initExpect; //  initial expectation mu = expectation rating
}

parameters {
  vector<lower=0,upper=1>[nSubject] alphaCoefI;
  vector<lower=0,upper=1>[nSubject] alphaCoefC;
  vector<lower=0>[nSubject] painError;
  vector<lower=0>[nSubject] expectError;

  // Group-level parameters
  real <lower=0> alphaI;
  real <lower=0> betaI;
  real <lower=0> alphaC;
  real <lower=0> betaC;
 // Group-level parameters
  real <lower=0> expectScale;
  real <lower=0> painScale;
}

transformed parameters {
  matrix[nTrials, nSubject] painMu; 
  matrix[nTrials, nSubject] expectMu; 
  matrix[nTrials, nSubject] predErr;
  vector<lower=0>[nSubject] stretchAlphaCoefI;
  vector<lower=0>[nSubject] stretchAlphaCoefC;
  
  for (s in 1:nSubject) {
    stretchAlphaCoefC[s] = myMulti * alphaCoefC[s];
    stretchAlphaCoefI[s] = myMulti * alphaCoefI[s];
  }

  for (s in 1:nSubject) {
    for (m in 1:nCues) {
      expectMu[firstTimeCue[m,s], s] = initExpect[m, s];
    }
    for (t in 1:(nTrials-nCues)){
      painMu[t,s] = noxInput[t,s];
      predErr[t,s] = painMu[t,s] - expectMu[t,s];
      if (((predErr[t,s] > 0) && (cueType[t,s] == -1)) || ((predErr[t,s] < 0) && (cueType[t,s] == 1))) {
          expectMu[cueIndexVec[t,s],s] = expectMu[t,s] + stretchAlphaCoefI[s] * predErr[t,s] ;
      } else {
          expectMu[cueIndexVec[t,s],s] = expectMu[t,s] + stretchAlphaCoefC[s] * predErr[t,s] ;
      }
    }
    for (t in (nTrials - nCues):nTrials){
      painMu[t,s] =  noxInput[t,s];
      predErr[t,s] = painMu[t,s] - expectMu[t,s];
    }
  }
}

model {
  for (s in 1:nSubject){

    for (t in 1:nTrials){
      target += normal_lpdf(painRating[t,s] | painMu[t,s], painError[s]);
      target += normal_lpdf(expectPain[t,s] | expectMu[t,s], expectError[s]);
    }
    
    target += beta_lpdf(alphaCoefI[s] | alphaI, betaI);
    target += beta_lpdf(alphaCoefC[s] | alphaC, betaC);

    target += student_t_lpdf(painError[s] | 1, 0, painScale) - log(0.5);
    target += student_t_lpdf(expectError[s] | 1, 0, expectScale) - log(0.5);

  }
  
  // Hierarchical Priors
  target += uniform_lpdf(alphaI | 0 , uniformPriorParameter) ; 
  target += uniform_lpdf(alphaC | 0 , uniformPriorParameter) ; 
  target += uniform_lpdf(betaI | 0 , uniformPriorParameter) ; 
  target += uniform_lpdf(betaC | 0 , uniformPriorParameter) ; 

  target += uniform_lpdf(painScale | 0 , uniformPriorParameter) ; 
  target += uniform_lpdf(expectScale | 0 , uniformPriorParameter) ; 
}

