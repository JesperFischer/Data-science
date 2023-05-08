


do_parameterrecover = function(data){
  index = rnorm(1,-1000,1000)
  
  df = poweragent_gen(m = data$m,
                      n = data$n,
                      q1 = data$q1,
                      q2 = data$q2,
                      kw = data$kw,
                      kc = data$kc,
                      kappa = data$kappa,
                      autocor = data$autocor,
                      bw = data$bw,
                      bc = data$bc,
                      bp = data$bp,
                      dist = data$dist)
  
  
  data1 = list(Tw1 = df$Tw1,
               Tw2 = df$Tw2,
               Tc1 = df$Tc1,
               Tc2 = df$Tc2, 
               N = nrow(df), 
               w = df$warm, 
               c = df$cold, 
               p = df$pain,
               dist = dist)
  
    model = here::here("concise","realstuf","stan","diamond_simplest.stan")
    mod = cmdstan_model(model)
    
    fit <- mod$sample(
      data = data1, 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    
    return(data.frame(fit$summary(c("q1","kappa")), reals = c(data$q1, data$kappa), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
}

do_parameterrecover2 = function(data){
  index = rnorm(1,-1000,1000)
  
  df = poweragent_gen(m = data$m,
                      n = data$n,
                      q1 = data$q1,
                      q2 = data$q2,
                      kw = data$kw,
                      kc = data$kc,
                      kappa = data$kappa,
                      autocor = data$autocor,
                      bw = data$bw,
                      bc = data$bc,
                      bp = data$bp,
                      dist = data$dist)
  
  
  data1 = list(Tw1 = df$Tw1,
               Tw2 = df$Tw2,
               Tc1 = df$Tc1,
               Tc2 = df$Tc2, 
               N = nrow(df), 
               w = df$warm, 
               c = df$cold, 
               p = df$pain,
               dist = dist)
  
  model = here::here("concise","realstuf","stan","real_diamond.stan")
  mod = cmdstan_model(model)
  
  fit <- mod$sample(
    data = data1, 
    seed = 123, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500
  )
  
  return(data.frame(fit$summary(c("q1","kw","kc","kappa")), reals = c(data$q1,data$kw, data$kc, data$kappa), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
}

do_parameterrecover3 = function(data){
  index = rnorm(1,-1000,1000)
  
  df = poweragent_gen(m = data$m,
                      n = data$n,
                      q1 = data$q1,
                      q2 = data$q2,
                      kw = data$kw,
                      kc = data$kc,
                      kappa = data$kappa,
                      autocor = data$autocor,
                      bw = data$bw,
                      bc = data$bc,
                      bp = data$bp,
                      dist = data$dist)
  
  
  data1 = list(Tw1 = df$Tw1,
               Tw2 = df$Tw2,
               Tc1 = df$Tc1,
               Tc2 = df$Tc2, 
               N = nrow(df), 
               w = df$warm, 
               c = df$cold, 
               p = df$pain,
               dist = dist)
  
  model = here::here("concise","realstuf","stan","real_diamond_with_auto_on_thermo.stan")
  mod = cmdstan_model(model)
  
  
  fitter = function(){
    fit <- mod$sample(
      data = data1,
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    return(fit)
  }
  
  
  fit <- tryCatch({
    fit = withTimeout(fitter(), timeout = 600, cpu=600)
  }, error = function(e) {
    print("An error occurred!")
    return(NA)
  })
  
  if(is.environment(fit)){
    parameters = c("q1","kw","kc","bw","bc","kappa")
    np = length(parameters)
    return(data.frame(fit$summary(parameters), reals = c(data$q1,data$kw, data$kc,data$bw,data$bc, data$kappa), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
  }else{
    parameters = c("q1","kw","kc","bw","bc","kappa")
    np = length(parameters)
    return(data.frame(variable = parameters, 
               mean = rep(NA,np),
               median = rep(NA,np), 
               sd = rep(NA,np), 
               q5 = rep(NA,np), 
               q95 = rep(NA,np), rhat = rep(NA,np), ess_bulk= rep(NA,np), ess_tail= rep(NA,np),
               reals = c(data$q1,data$kw, data$kc,data$bw,data$bc, data$kappa),
               div = rep(NA,np), index = index))
  }
  
  
  
}


do_parameterrecover4 = function(data){
  index = rnorm(1,-1000,1000)
  
  df = poweragent_gen(m = data$m,
                      n = data$n,
                      q1 = data$q1,
                      q2 = data$q2,
                      kw = data$kw,
                      kc = data$kc,
                      kappa = data$kappa,
                      autocor = data$autocor,
                      bw = data$bw,
                      bc = data$bc,
                      bp = data$bp,
                      dist = data$dist)
  
  
  data1 = list(Tw1 = df$Tw1,
               Tw2 = df$Tw2,
               Tc1 = df$Tc1,
               Tc2 = df$Tc2, 
               N = nrow(df), 
               w = df$warm, 
               c = df$cold, 
               p = df$pain,
               dist = dist)
  
  model = here::here("concise","realstuf","stan","powerlaw_fixed_k.stan")
  mod = cmdstan_model(model)
  
  fit <- mod$sample(
    data = data1, 
    seed = 123, 
    chains = 4, 
    parallel_chains = 1,
    refresh = 500
  )
  
  return(data.frame(fit$summary(c("q1","m","n","kappa")), reals = c(data$q1,data$m, data$n, data$kappa), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
}


do_parameterrecover5 = function(data){
  index = rnorm(1,-1000,1000)
  
  df = poweragent_gen(m = data$m,
                      n = data$n,
                      q1 = data$q1,
                      q2 = data$q2,
                      kw = data$kw,
                      kc = data$kc,
                      kappa = data$kappa,
                      autocor = data$autocor,
                      bw = data$bw,
                      bc = data$bc,
                      bp = data$bp,
                      dist = data$dist)
  
  
  data1 = list(Tw1 = df$Tw1,
               Tw2 = df$Tw2,
               Tc1 = df$Tc1,
               Tc2 = df$Tc2, 
               N = nrow(df), 
               w = df$warm, 
               c = df$cold, 
               p = df$pain,
               dist = dist)
  
  model = here::here("concise","realstuf","stan","powerlaw_fixed_k_with_auto_on_thermo.stan")
  mod = cmdstan_model(model)
  
  fit <- mod$sample(
    data = data1, 
    seed = 123, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500
  )
  
  return(data.frame(fit$summary(c("q1","m","n","bw","bc","kappa")), reals = c(data$q1,data$m, data$n,data$bw,data$bc, data$kappa), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
}


do_parameterrecover6 = function(data){
  index = rnorm(1,-1000,1000)
  
  df = poweragent_gen(m = data$m,
                      n = data$n,
                      q1 = data$q1,
                      q2 = data$q2,
                      kw = data$kw,
                      kc = data$kc,
                      kappa = data$kappa,
                      autocor = data$autocor,
                      bw = data$bw,
                      bc = data$bc,
                      bp = data$bp,
                      dist = data$dist)
  
  
  data1 = list(Tw1 = df$Tw1,
               Tw2 = df$Tw2,
               Tc1 = df$Tc1,
               Tc2 = df$Tc2, 
               N = nrow(df), 
               w = df$warm, 
               c = df$cold, 
               p = df$pain,
               dist = dist)
  
  model = here::here("concise","realstuf","stan","powerlaw_fixed_k_with_auto_general.stan")
  mod = cmdstan_model(model)
  
  fit <- mod$sample(
    data = data1, 
    seed = 123, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500
  )
  
  return(data.frame(fit$summary(c("q1","m","n","bw","bc","bp","kappa")), reals = c(data$q1,data$m, data$n,data$bw,data$bc,data$bp, data$kappa), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
}


do_parameterrecover7 = function(data){
  index = rnorm(1,-1000,1000)
  
  df = poweragent_gen(m = data$m,
                      n = data$n,
                      q1 = data$q1,
                      q2 = data$q2,
                      kw = data$kw,
                      kc = data$kc,
                      kappa = data$kappa,
                      autocor = data$autocor,
                      bw = data$bw,
                      bc = data$bc,
                      bp = data$bp,
                      dist = data$dist)
  
  
  data1 = list(Tw1 = df$Tw1,
               Tw2 = df$Tw2,
               Tc1 = df$Tc1,
               Tc2 = df$Tc2, 
               N = nrow(df), 
               w = df$warm, 
               c = df$cold, 
               p = df$pain,
               dist = dist)
  
  model = here::here("concise","realstuf","stan","full_power(sameq).stan")
  mod = cmdstan_model(model)
  
  fit <- mod$sample(
    data = data1, 
    seed = 123, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500
  )
  
  return(data.frame(fit$summary(c("q1","m","n","bw","bc","bp","kc","kw","kappa")), reals = c(data$q1,data$m, data$n,data$bw,data$bc,data$bp,data$kc,data$kw, data$kappa), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
}






plot_prior_posterior_update = function(fit,parameters,reals){
  
  
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
    geom_vline(data = data.frame(reals, name = as.factor(rownames(data.frame(reals)))), aes(xintercept = reals))+
    facet_wrap(~name, scales = "free")+
    theme_classic()
  
}





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
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    
    
    data = list(Tw1 = tgi$Tw1,Tw2 = (tgi$Tw2),Tc1 = (tgi$Tc1),Tc2 = (tgi$Tc2), N = nrow(tgi), w = (tgi$warm1), c = (tgi$cold1), p = tgi$pain1, dist = tgi$dist[1])
    
    fit <- model22$sample(
      data = data, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    
    
    loobm = loo(fitbm$draws(c("log_lik")), r_eff = relative_eff(exp(fitbm$draws(c("log_lik")))))
    
    loofit = loo(fit$draws(c("log_lik")), r_eff = relative_eff(exp(fit$draws(c("log_lik")))))
    
    loocompar = loo_compare(loobm, loofit)
    
    loocompar = data.frame(loocompar, id = index)
    loocompar$models = c(basename(get(rownames(loocompar)[1])),basename(get(rownames(loocompar)[2])))
    
    div = data.frame(diamond = sum(fit$diagnostic_summary()$num_divergent), brms = sum(fitbm$diagnostic_summary()$num_divergent))
    
    
    weights = loo::loo_model_weights(list(loobm,loofit))
    
    weights = data.frame(weights = weights[1:2], id = index)
    weights$model = c(basename(get(rownames(weights)[1])),basename(get(rownames(weights)[2])))
    
    
    return(list(loocompar = loocompar,weights = weights, div = div))
    
  }
  else{
    
    data = list(Tw1 = tgi$Tw1,Tw2 = (tgi$Tw2),Tc1 = (tgi$Tc1),Tc2 = (tgi$Tc2), N = nrow(tgi), w = (tgi$warm1), c = (tgi$cold1), p = tgi$pain1, dist = tgi$dist[1])
    
    model11 = cmdstan_model(model1)
    
    fit1 <- model11$sample(
      data = data,
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    
    model22 = cmdstan_model(model2)
    
    fit2 <- model22$sample(
      data = data,
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    
    loofit1 = loo(fit1$draws(c("log_lik")), r_eff = relative_eff(exp(fit1$draws(c("log_lik")))))
    
    loofit2 = loo(fit2$draws(c("log_lik")), r_eff = relative_eff(exp(fit2$draws(c("log_lik")))))
    
    loocompar = loo_compare(loofit1, loofit2)
    loocompar = data.frame(loocompar, id = index)
    
    
    loocompar$models = c(basename(get(rownames(loocompar)[1])),basename(get(rownames(loocompar)[2])))
    
    div = data.frame(model1= sum(fit1$diagnostic_summary()$num_divergent), model2 = sum(fit2$diagnostic_summary()$num_divergent))
    
    
    weights = loo::loo_model_weights(list(loofit1,loofit2))
    
    weights = data.frame(weights = weights[1:2], id = index)
    weights$model = c(basename(model1),basename(model2))
    
    
    return(list(loocompar = loocompar,weights = weights, div = div))
    
  }
  
}

