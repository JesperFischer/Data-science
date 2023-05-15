
get_experiment = function(){
  trials = c(40,40,40,40)
  bias = rep(c(0.2,0.8,0.3,0.7),trials)
  
  
  cue = rbinom(sum(trials),1,0.5)
  
  stim = array(sum(trials, NA))
  for(i in 1:sum(trials)){
    stim[i] = ifelse(cue[i] == 1, rbinom(1,1,bias[i]), rbinom(1,1,(1-bias[i])))
    
  }
  
  u = ifelse(cue == stim, 1,0)
  dd = data.frame(stim = stim, cue = cue, u = u, bias = bias)
  
  return(dd)
  
  
}




our_rw_agent = function(parameters){
  

  dd = get_experiment()
  
  cue = dd$cue
  stim = dd$stim
  u = dd$u
  bias = dd$bias
  
  
  
  ntrials = nrow(dd)
  
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
    
    #get the expectation of the first trial i.e. (how likeli is hot?)
    exp[i] = ifelse(cue[i] == 1, association[i], 1-association[i])
    #the prediction
    pred[i] = rbinom(1,1,(exp[i]^beta)/((exp[i]^beta)+(1-exp[i])^(beta)))
    
    #the mean percept given from the stimulus received and the expectation
    perceptmu[i] = ifelse(w1*stim[i]+(1-w1)*exp[i] == 0, 0.01,
                          ifelse(w1*stim[i]+(1-w1)*exp[i] == 1, 0.99, w1*stim[i]+(1-w1)*exp[i]))
    
    #generating the percept
    percept[i] = extraDistr::rprop(1,precision_percept,perceptmu[i])
    
    #now we get the binary percept (did you feel hot or warm?)
    percept_bin[i] = rbinom(1,1,(perceptmu[i]^beta)/((perceptmu[i]^beta)+(1-perceptmu[i])^(beta)))
    
    #prediction error for the association
    pe[i] = ifelse(cue[i] == 1, (perceptmu[i] -exp[i]),(-(perceptmu[i] -exp[i])))  
    
    #association update
    association[i+1] = association[i]+alpha*pe[i]
    
  }
  
  df = data.frame(u = u[1:ntrials], stim = stim[1:ntrials], percept = percept[1:ntrials], pred = pred[1:ntrials],association = association[1:ntrials],
                  exp = exp[1:ntrials],pe = pe[1:ntrials],cue = cue[1:ntrials], alpha = alpha,
                  percept_bin = percept_bin[1:ntrials],
                  w1 = w1,precision_percept = precision_percept,
                  x = 1:ntrials, perceptmu = perceptmu[1:ntrials], desired = rep(bias,1), id = rnorm(1,0,1))
  
  return(df)

  
  
}

our_rw_agent_v2 = function(w1,alpha,precision_percept,beta){
  
  
  dd = get_experiment()
  
  cue = dd$cue
  stim = dd$stim
  u = dd$u
  bias = dd$bias
  
  
  
  ntrials = nrow(dd)
  
  percept = array(NA, ntrials)
  perceptmu = array(NA, ntrials)
  exp = array(NA, ntrials)
  
  association = array(NA, ntrials)
  pred = array(NA,ntrials)
  pe = array(NA,ntrials)
  percept_bin = array(NA,ntrials)
  
  
  association[1] = 0.5
  
  
  for (i in seq(ntrials)){
    
    #get the expectation of the first trial i.e. (how likeli is hot?)
    exp[i] = ifelse(cue[i] == 1, association[i], 1-association[i])
    #the prediction
    pred[i] = rbinom(1,1,(exp[i]^beta)/((exp[i]^beta)+(1-exp[i])^(beta)))
    
    #the mean percept given from the stimulus received and the expectation
    perceptmu[i] = ifelse(w1*stim[i]+(1-w1)*exp[i] == 0, 0.01,
                          ifelse(w1*stim[i]+(1-w1)*exp[i] == 1, 0.99, w1*stim[i]+(1-w1)*exp[i]))
    
    #generating the percept
    percept[i] = extraDistr::rprop(1,precision_percept,perceptmu[i])
    
    #now we get the binary percept (did you feel hot or warm?)
    percept_bin[i] = rbinom(1,1,(perceptmu[i]^beta)/((perceptmu[i]^beta)+(1-perceptmu[i])^(beta)))
    
    #prediction error for the association
    pe[i] = ifelse(cue[i] == 1, (perceptmu[i] -exp[i]),(-(perceptmu[i] -exp[i])))  
    
    #association update
    association[i+1] = association[i]+alpha*pe[i]
    
  }
  
  df = data.frame(u = u[1:ntrials], stim = stim[1:ntrials], percept = percept[1:ntrials], pred = pred[1:ntrials],association = association[1:ntrials],
                  exp = exp[1:ntrials],pe = pe[1:ntrials],cue = cue[1:ntrials], alpha = alpha,
                  percept_bin = percept_bin[1:ntrials],
                  w1 = w1,precision_percept = precision_percept,
                  x = 1:ntrials, perceptmu = perceptmu[1:ntrials], desired = rep(bias,1), id = rnorm(1,0,1))
  
  return(df)
  
  
  
}


our_hier_rw_agent = function(parameters){
  nsubs = parameters$nsubs
  
  df = data.frame()
  for (s in 1:nsubs){
    df1 = our_rw_agent_v2(w1 = extraDistr::rprop(1,parameters$kappa_w1,parameters$mu_w1),
                 alpha = extraDistr::rprop(1,parameters$kappa_alpha,parameters$mu_alpha),
                 precision_percept = rlnorm(1,parameters$mu_precision_percept,parameters$sd_precision_percept),
                 beta = rlnorm(1,parameters$mu_beta,parameters$sd_beta)
                 )
    df = rbind(df,df1)
    
  }
  
  
  return(df)
  
  
  
}







