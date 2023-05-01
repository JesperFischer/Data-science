#osf replica full RW


RW_agent = function(alphaI, betaI, alphaY, betaY, painScale, expectScale, nTrials, nSubject){
  
  
  painRating = array(NA, c(nTrials,nSubject))
  expectpain = array(NA, c(nTrials,nSubject))
  
  noxInput = matrix(rnorm(79, 0,1), nrow = nTrials-1, ncol=nSubject)
  noxInput = ifelse(noxInput > 0, 1, 0)
  
  
  #int<lower=1> myMulti;
  #int<lower=1> uniformPriorParameter;
  
  painMu = array(NA, c(nTrials, nSubject))
  expectMu = array(NA, c(nTrials, nSubject))
  predErr = array(NA, c(nTrials, nSubject))
  
  alphaCoefI = array(NA, nSubject)
  yParam = array(NA, nSubject)
  painError = array(NA, nSubject)
  expectError = array(NA, nSubject)
  
  
  
  
  
  for (s in 1:nSubject) {
    expectMu[1,s] = 0.5
    
    alphaCoefI[s] = extraDistr::rprop(1, alphaI, betaI)
    yParam[s] = extraDistr::rprop(1, alphaY, betaY)
    
    painError[s] = rlnorm(1,-2,painScale)
    expectError[s] = rlnorm(1,-2,expectScale)
    
    
    
    for (t in 1:(nTrials-1)){
      painMu[t,s] = (1-yParam[s]) * noxInput[t,s] + yParam[s] * expectMu[t,s]
      predErr[t,s] = painMu[t,s] - expectMu[t,s]
      
      expectMu[t+1,s] = expectMu[t,s] + alphaCoefI[s] * predErr[t,s]
      
      
      painRating[t,s]  = extraDistr::rprop(1, 1/painError[s],painMu[t,s])
      expectpain[t,s]  = extraDistr::rprop(1, 1/expectError[s],expectMu[t,s])
    }
  }
  
  qq = data.frame(painRating) %>% mutate(trial = 1:nrow(.)) %>% filter(trial != nTrials) %>% pivot_longer(cols = starts_with("X"), values_to = "PainRating") %>% mutate(name = as.factor(name))
  qq1 = data.frame(expectpain) %>% mutate(trial = 1:nrow(.)) %>% filter(trial != nTrials) %>% pivot_longer(cols = starts_with("X"), values_to = "expectRating")%>% mutate(name = as.factor(name))
  qq2 = data.frame(noxInput) %>% mutate(trial = 1:nrow(.)) %>% pivot_longer(cols = starts_with("X"), values_to = "nox")%>% mutate(name = as.factor(name))
  
  
  df = inner_join(qq,qq1)
  df = inner_join(df, qq2)
  
  df$name = as.factor(df$name)
  
  subparm = data.frame(yParam =yParam, alphaCoefI = alphaCoefI, painError = painError, expectError = expectError)
  
  return(list(df,subparm))
  
}



RW_agent2 = function(alphaI, betaI, alphaY, betaY, painScale, expectScale, nTrials, nSubject, noxInput){
  
  
  painRating = array(NA, c(nTrials,nSubject))
  expectpain = array(NA, c(nTrials,nSubject))
  
  
  
  #int<lower=1> myMulti;
  #int<lower=1> uniformPriorParameter;
  
  painMu = array(NA, c(nTrials, nSubject))
  expectMu = array(NA, c(nTrials, nSubject))
  predErr = array(NA, c(nTrials, nSubject))
  
  alphaCoefI = array(NA, nSubject)
  yParam = array(NA, nSubject)
  painError = array(NA, nSubject)
  expectError = array(NA, nSubject)
  
  
  
  

    expectMu[1,1] = 0.5
    
    alphaCoefI[1] = betaI
    yParam[1] = betaY
    
    painError[1] = painScale
    expectError[1] = expectScale
    
    
    
    for (t in 1:(nTrials-1)){
      painMu[t,1] = (1-yParam[1]) * noxInput[t,1] + yParam[1] * expectMu[t,1]
      predErr[t,1] = painMu[t,1] - expectMu[t,1]
      
      expectMu[t+1,1] = expectMu[t,1] + alphaCoefI[1] * predErr[t,1]
      
      
      painRating[t,1]  = extraDistr::rprop(1, 1/painError[1],painMu[t,1])
      expectpain[t,1]  = extraDistr::rprop(1, 1/expectError[1],expectMu[t,1])
    }

  
  qq = data.frame(painRating) %>% mutate(trial = 1:nrow(.)) %>% filter(trial != nTrials) %>% pivot_longer(cols = starts_with("X"), values_to = "PainRating") %>% mutate(name = as.factor(name))
  qq1 = data.frame(expectpain) %>% mutate(trial = 1:nrow(.)) %>% filter(trial != nTrials) %>% pivot_longer(cols = starts_with("X"), values_to = "expectRating")%>% mutate(name = as.factor(name))
  qq2 = data.frame(noxInput) %>% mutate(trial = 1:nrow(.)) %>% pivot_longer(cols = starts_with("X"), values_to = "nox")%>% mutate(name = as.factor(name))
  
  
  df = inner_join(qq,qq1)
  df = inner_join(df, qq2)
  
  df$name = as.factor(df$name)
  
  subparm = data.frame(yParam =yParam, alphaCoefI = alphaCoefI, painError = painError, expectError = expectError)
  
  return(list(df,subparm))
  
}



plot_prior_posterior_update_hier = function(fit,parameters,reals){
  
  
  post = as_draws_df(fit$draws()) %>% 
    select(all_of(parameters)) %>% 
    mutate(posterior = T)
  
  
  priors = as_draws_df(fit$draws()) %>% 
    select(all_of(paste0("prior_",parameters)))%>% 
    mutate(posterior = F) %>%
    rename_with(~sub("^prior_", "", .), starts_with("prior_"))
  df = rbind(post,priors)
  
  df %>% pivot_longer(cols = parameters) %>% ggplot(aes(x = value, fill = posterior))+
    geom_histogram(alpha = 0.5, position="identity")+
    geom_vline(data = data.frame(reals, name = names(((reals)))), aes(xintercept = reals))+
    facet_wrap(~name, scales = "free")+
    theme_classic()
  
}




plot_prior_posterior_update_sub = function(fit,parameters,reals){


  
  post = as_draws_df(fit$draws()) %>% 
    select(starts_with(parameters)) %>% 
    mutate(posterior = T)
  
  priors = as_draws_df(fit$draws()) %>% 
    select(starts_with(paste0("prior_",parameters)))%>% 
    mutate(posterior = F) %>%
    rename_with(~sub("^prior_", "", .), starts_with("prior_"))
  
  
  
  
  df = rbind(post,priors)
  
  df %>% pivot_longer(cols = -posterior) %>% ggplot(aes(x = value, fill = posterior))+
    geom_histogram(alpha = 0.5, position="identity")+facet_wrap(~name, nrow = length(parameters), ncol = 10, scales = "free")+
    geom_vline(data = reals, aes(xintercept = value))+
    theme_classic()
  
}
