data {
  int<lower=1> nCues; // number of cues
  int<lower=1> nSubject; // number of subjects
  int<lower=1> nTrials; // number of subjects
  int<lower=0> sigmaEpsilon;
  int<lower=1> sigmaPsi;
  int<lower=1> myMultiBeta;
  int<lower=1> uniformPriorParameter;
  matrix[nTrials, nSubject] painRating; // observations
  matrix[nTrials, nSubject] expectPain; // observations
  matrix[nTrials, nSubject] noxInput; // observations
  matrix[nCues, nSubject] sortedInitExp;
  int cueIndexVec[(nTrials - nCues), nSubject]; // indexing next trial for spec cue
  int firstTimeCue[nCues, nSubject];
  int allCues[nTrials, nSubject];
  matrix[nCues,  nSubject] initExpect; //  initial expectation mu = expectation rating
}

parameters {
  vector<lower=0, upper=1>[nSubject] myBeta;
  vector<lower=0>[nSubject] sigmaEta;
  vector<lower=0>[nSubject] painError;
  vector<lower=0>[nSubject] expectError;
  
  // Group-level parameters
  real <lower=0> alphaBeta;
  real <lower=0> betaBeta;
  real <lower=0> expectScale;
  real <lower=0> painScale;
  real <lower=0> etaScale;
}

transformed parameters {
  matrix[nTrials, nSubject] painMu; 
  matrix[nTrials, nSubject] painVar; 
  matrix[nTrials, nSubject] expectMu; 
  matrix[nTrials, nSubject] expectVar; 
  vector<lower=0, upper=2>[nSubject] stretchMyBeta;


  for (s in 1:nSubject) {
    stretchMyBeta[s] = myMultiBeta * myBeta[s];

    for (m in 1:nCues) {
      expectMu[firstTimeCue[m, s],s] = initExpect[m, s];
      expectVar[firstTimeCue[m,s], s] = sigmaEta[s]    ;
    }
    for (t in 1:(nTrials - nCues)){
      painMu[t,s] =  ( sigmaEpsilon * expectMu[t,s] + (sigmaPsi + expectVar[t,s]) * noxInput[t,s] ) / 
                       (sigmaEpsilon + sigmaPsi + expectVar[t,s]);
      painVar[t,s] = ( sigmaEpsilon * (sigmaPsi + expectVar[t,s]) ) / (sigmaEpsilon + sigmaPsi + expectVar[t,s]) ;       
      expectMu[cueIndexVec[t,s],s] = (stretchMyBeta[s] * (sigmaEpsilon + sigmaPsi) * expectMu[t,s] + (stretchMyBeta[s] * expectVar[t,s] * noxInput[t,s]) + 
                                      ((1 - stretchMyBeta[s]) * (sigmaEpsilon + sigmaPsi + expectVar[t,s]) * sortedInitExp[allCues[t,s],s] )) / 
                                      (sigmaEpsilon + sigmaPsi + expectVar[t,s]) ;
      expectVar[cueIndexVec[t,s],s] = (stretchMyBeta[s] * (sigmaEpsilon + sigmaPsi) * expectVar[t,s] / 
                                       (sigmaEpsilon + sigmaPsi + expectVar[t,s])) + sigmaEta[s];
    }
    for (t in (nTrials - nCues):nTrials){
      painMu[t,s] =  ( sigmaEpsilon * expectMu[t,s] + (sigmaPsi + expectVar[t,s]) * noxInput[t,s] ) / 
                       (sigmaEpsilon + sigmaPsi + expectVar[t,s]);
      painVar[t,s] = ( sigmaEpsilon * (sigmaPsi + expectVar[t,s]) ) / (sigmaEpsilon + sigmaPsi + expectVar[t,s]);   
    }
  }
}

model {
   for (s in 1:nSubject){
    
    for (t in 1:nTrials){
      target += normal_lpdf(painRating[t,s] | painMu[t,s], painError[s]);
      target += normal_lpdf(expectPain[t,s] | expectMu[t,s], expectError[s]);
    }
    
    target += beta_lpdf(myBeta[s] | alphaBeta, betaBeta);
    
    target += student_t_lpdf(sigmaEta[s] | 1, 0, etaScale) - log(0.5);
    target += student_t_lpdf(painError[s] | 1, 0, painScale) - log(0.5);
    target += student_t_lpdf(expectError[s] | 1, 0, expectScale) - log(0.5);
  
  }
  
  // Hierarchical Priors
  target += uniform_lpdf(alphaBeta | 0 , uniformPriorParameter) ; 
  target += uniform_lpdf(betaBeta | 0 , uniformPriorParameter) ; 

  target += uniform_lpdf(etaScale | 0 , uniformPriorParameter) ; 
  target += uniform_lpdf(painScale | 0 , uniformPriorParameter) ; 
  target += uniform_lpdf(expectScale | 0 , uniformPriorParameter) ; 
}
