



model_comparison = function(tgi, model1, model2, brms = FALSE, formula = FALSE){

  index = runif(1,-1000,1000)
  
  if(brms == TRUE){
    if(str_detect(model1, "brms")){
      model11 = cmdstan_model(model1)
      model22 = cmdstan_model(model2)
    }else{
      print("first model has to be the brms thanks,...")
    }
    
    standata = make_standata(formula, tgi)
    
    fitbm <- model11$sample(
      data = standata, 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    
    
    data = list(Tw1 = tgi$Tw1,Tw2 = (tgi$Tw2),Tc1 = (tgi$Tc1),Tc2 = (tgi$Tc2), N = nrow(tgi), w = (tgi$warm1), c = (tgi$cold1), p = tgi$pain1)
    
    fit <- model22$sample(
      data = data, 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    
    
    loobm = loo(fitbm$draws(c("log_lik")), r_eff = relative_eff(exp(fitbm$draws(c("log_lik")))))
    
    loofit = loo(fit$draws(c("log_lik")), r_eff = relative_eff(exp(fit$draws(c("log_lik")))))
    
    loocompar = loo_compare(loobm, loofit)
    
    loocompar = data.frame(loocompar, id = index)
    loocompar$models = c(basename(get(rownames(loocompar)[1])),basename(get(rownames(loocompar)[2])))
    
    weights = loo::loo_model_weights(list(loobm,loofit))
    
    weights = data.frame(weights = weights[1:2], id = index)
    weights$model = c(basename(get(rownames(weights)[1])),basename(get(rownames(weights)[2])))


  return(list(loocompar = loocompar,weights = weights))
    
  }
  else{
    
  data = list(Tw1 = tgi$Tw1,Tw2 = (tgi$Tw2),Tc1 = (tgi$Tc1),Tc2 = (tgi$Tc2), N = nrow(tgi), w = (tgi$warm1), c = (tgi$cold1), p = tgi$pain1)
  
  model11 = cmdstan_model(model1)
  
  fit1 <- model11$sample(
    data = data, 
    seed = 123, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500
  )
  
  model22 = cmdstan_model(model2)
  
  fit2 <- model22$sample(
    data = data, 
    seed = 123, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500
  )
  
  loofit1 = loo(fit1$draws(c("log_lik")), r_eff = relative_eff(exp(fit1$draws(c("log_lik")))))
  
  loofit2 = loo(fit2$draws(c("log_lik")), r_eff = relative_eff(exp(fit2$draws(c("log_lik")))))
  
  loocompar = loo_compare(loofit1, loofit2)
  loocompar = data.frame(loocompar, id = index)
  
  
  loocompar$models = c(basename(get(rownames(loocompar)[1])),basename(get(rownames(loocompar)[2])))
  
  
  weights = loo::loo_model_weights(list(loofit1,loofit2))
  
  weights = data.frame(weights = weights[1:2], id = index)
  weights$model = c(basename(model1),basename(model2))
  
  
  return(list(loocompar = loocompar,weights = weights))
  
  }
  
}

