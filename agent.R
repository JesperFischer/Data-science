

rm_agent = function(bias,trials){
  u = c()
  for (i in 1:length(trials)){
    u1 = rbinom(n = trials[i], size = 1, prob = bias[i])
    u = c(u,u1)
  }
  return(u)
}




our_rw_agent = function(parameters){
  

  trials = c(40,40,40,40)
  bias = rep(c(0.2,0.8,0.3,0.7),trials)
  
  
  cue = rbinom(sum(trials),1,0.5)
  
  stim = array(sum(trials, NA))
  for(i in 1:sum(trials)){
    stim[i] = ifelse(cue[i] == 1, rbinom(1,1,bias[i]), rbinom(1,1,(1-bias[i])))
    
  }
  
  u = ifelse(cue == stim, 1,0)
  dd = data.frame(stim = stim, cue = cue, u = u)
  
  ntrials = sum(trials)
  
  percept = array(NA, ntrials)
  perceptmu = array(NA, ntrials)
  exp = array(NA, ntrials)
  
  association = array(NA, ntrials)
  pred = array(NA,ntrials)
  pe = array(NA,ntrials)
  percept_bin = array(NA,ntrials)
  
  
  w1 = parameters$w1
  alpha = parameters$alpha
  precision_percept = parameters$precision_percept
  beta = parameters$beta
  
  
  association[1] = 0.5
  

  for (i in seq(ntrials)){
    
    exp[i] = ifelse(cue[i] == 1, association[i], 1-association[i])
    
    
    
    pred[i] = rbinom(1,1,(exp[i]^beta)/((exp[i]^beta)+(1-exp[i])^(beta)))
    
    
    perceptmu[i] = ifelse(w1*stim[i]+(1-w1)*exp[i] == 0, 0.01,
                          ifelse(w1*stim[i]+(1-w1)*exp[i] == 1, 0.99, w1*stim[i]+(1-w1)*exp[i]))
    
    percept[i] = extraDistr::rprop(1,precision_percept,perceptmu[i])
    
    percept_bin[i] = rbinom(1,1,(perceptmu[i]^beta)/((perceptmu[i]^beta)+(1-perceptmu[i])^(beta)))
    
    pe[i] = ifelse(cue[i] == 1, (perceptmu[i] -exp[i]),(-(perceptmu[i] -exp[i])))  
    
    association[i+1] = association[i]+alpha*pe[i]
    
  }
  
  length = sum(trials)
  
  df = data.frame(u = u[1:length], stim = stim[1:length], percept = percept[1:length], pred = pred[1:length],association = association[1:length],
                  exp = exp[1:length],pe = pe[1:length],cue = cue[1:length], alpha = alpha,
                  percept_bin = percept_bin[1:length],
                  w1 = w1,precision_percept = precision_percept,
                  x = 1:length, perceptmu = perceptmu[1:length], desired = rep(bias,1), id = rnorm(1,0,1))
  
  return(df)

  
  
}








