

agent = function(b0a,b1a,b0p,b1p,b2p,kappa,outcome,stim, trials){
  

  
  expectation = array(NA, trials)
  percept = array(NA, trials)  
  alpha = array(NA, trials)  
  mupercept = array(NA, trials)    
  prediction = array(NA, trials)    
  
  
  expectation[1] = 0.5
  percept[1] = 0
  
  
  for (i in 2:trials){
  
    alpha[i] = inv_logit_scaled(b0a+b1a+percept[i-1])
    expectation[i] = expectation[i-1]+alpha[i]*(outcome[i]-expectation[i-1])
    mupercept[i] = inv_logit_scaled(b0p+b1p*expectation[i]+b2p*stim[i])
    
    percept[i] = extraDistr::rprop(1,kappa,mupercept[i])
    
    prediction[i] = rbinom(1,1,expectation[i])
   
  }
  
 
  return(data.frame(expectation, alpha, percept, prediction, trials = 1:trials, outcome = outcome, stim = stim)) 
}