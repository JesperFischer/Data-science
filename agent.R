

rm_agent = function(bias,trials){
  u = c()
  for (i in 1:length(trials)){
    u1 = rbinom(n = trials[i], size = 1, prob = bias[i])
    u = c(u,u1)
  }
  return(u)
}




our_rw_agent = function(parameters){
  
  library(here)
  
  source((here("~","Advanced-cognitive-modeling","assignment2","hgf_agent.R")))
  
  trials = c(50,50,50,50)
  bias = rep(c(0.2,0.8,0.3,0.7),trials)
  
  
  cue = rbinom(sum(trials),1,0.5)
  
  stim = array(sum(trials, NA))
  for(i in 1:200){
    stim[i] = ifelse(cue[i] == 1, rbinom(1,1,bias[i]), rbinom(1,1,(1-bias[i])))
    
  }
  
  u = ifelse(cue == stim, 1,0)
  dd = data.frame(stim = stim, cue = cue, u = u)
  
  dd %>% ggplot(aes(x = 1:nrow(.), y = u, col = as.factor(stim)))+geom_point()+theme_classic()
  
  ntrials = sum(trials)
  
  percept = array(NA, ntrials)
  perceptmu = array(NA, ntrials)
  exp = array(NA, ntrials)
  
  expectation = array(NA, ntrials)
  per_con = array(NA, ntrials)
  pred = array(NA,ntrials)
  pe = array(NA,ntrials)
  percept_bin = array(NA,ntrials)
  
  
  w1 = parameters$w1
  alpha = parameters$alpha
  
  
  expectation[1] = 0.5
  

  for (i in seq(ntrials)){
    pred[i] = rbinom(1,1,expectation[i])
    exp[i] = ifelse(cue[i] == 1, expectation[i], 1-expectation[i])
    
    perceptmu[i] = ifelse(w1*stim[i]+(1-w1)*exp[i] == 0, 0.01, ifelse(w1*stim[i]+(1-w1)*exp[i] == 1, 0.99,w1*stim[i]+(1-w1)*exp[i]))
    
    percept[i] = extraDistr::rprop(1,100,perceptmu[i])
    
    percept_bin[i] = rbernoulli(1,perceptmu[i])
    
    pe[i] = ifelse(cue[i] == 1, (perceptmu[i] -exp[i]),(-(perceptmu[i] -exp[i])))  
    
    expectation[i+1] = expectation[i]+alpha*pe[i]
    
  }
  
  
  df = data.frame(u = u[1:200], per_con = per_con[1:200],stim = stim[1:200], percept = percept[1:200], pred = pred[1:200],association = expectation[1:200],
                  exp = exp[1:200],pe = pe[1:200],cue = cue[1:200], alpha = alpha, w1 = w1, x = 1:200, perceptmu = perceptmu[1:200], desired = rep(bias,1))
  
  return(df)

  
  
}








