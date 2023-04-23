

agent = function(b0a,b1a,b0p,b1p,b2p,b3p,kappa,beta,outcome,stim, trials){
  

  
  expectation = array(NA, trials)
  percept = array(NA, trials)  
  alpha = array(NA, trials)  
  mupercept = array(NA, trials)    
  prediction = array(NA, trials)   
  prob = array(NA, trials)   
  
  
  
  expectation[1] = 0.5
  percept[1] = 0
  mupercept[1] = 0
  
  for (i in 2:trials){
  
    alpha[i] = inv_logit_scaled(b0a+b1a*mupercept[i-1])
    
    expectation[i] = expectation[i-1]+alpha[i]*(outcome[i-1]-expectation[i-1])
    
    mupercept[i] = inv_logit_scaled(b0p+b1p*(expectation[i]*(1-expectation[i]))+b2p*stim[i])
    
    
    prob[i] = expectation[i]^beta/((expectation[i]^beta)+(1-expectation[i])^(beta))
    
    percept[i] = extraDistr::rprop(1,kappa,mupercept[i])
    
    prediction[i] = rbinom(1,1,prob[i])
   
  }
  
 
  return(data.frame(expectation, alpha, percept, prediction, trials = 1:trials, outcome = outcome, stim = stim, prob = prob)) 
}




agent2 = function(b0a,b1a,b0p,b1p,b2p,kappa,outcome,stim, trials){
  
  
  
  expectation = array(NA, trials)
  percept = array(NA, trials)  
  alpha = array(NA, trials)  
  mupercept = array(NA, trials)    
  prediction = array(NA, trials)
  alphap = array(NA, trials)
  
  
  expectation[1] = 0.5
  percept[1] = 0
  mupercept[1] = 0.5
  
  for (i in 2:trials){
    
    alpha[i] = inv_logit_scaled(b0a+b1a*(percept[i-1]*(1-percept[i-1])))
    expectation[i] = expectation[i-1]+alpha[i]*(outcome[i-1]-expectation[i-1])
    
    alphap[i] = inv_logit_scaled(b0p+b1p*(expectation[i]*(1-expectation[i])))
    
    mupercept[i] = mupercept[i-1]+alphap[i]*(b2p*stim[i]-mupercept[i-1])
      
    
    percept[i] = extraDistr::rprop(1,kappa,mupercept[i])
    
    prediction[i] = rbinom(1,1,expectation[i])
    
  }
  
  
  return(data.frame(expectation, alpha, percept, prediction, trials = 1:trials, outcome = outcome, stim = stim, alphap = alphap, mupercept = mupercept)) 
}