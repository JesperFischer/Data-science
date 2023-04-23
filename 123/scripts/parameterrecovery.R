


do_parameterrecover = function(m,n,q1,q2,kappa,autocor,bw,bc,bp, dist){
  index = rnorm(1,-1000,1000)
  
  df = poweragent_generalized(m = m,
                              n = n,
                              q1 = q1,
                              q2 = q2,
                              kappa = kappa,
                              autocor = autocor,
                              bw = bw,
                              bc = bc,
                              bp = bp,
                              dist = dist)
  

  data1 = list(Tw1 = df$Tw1,
               Tw2 = df$Tw2,
               Tc1 = df$Tc1,
               Tc2 = df$Tc2, 
               N = nrow(df), 
               w = df$warm, 
               c = df$cold, 
               p = df$pain,
               dist = dist)
  
  if(q1 == q2 && autocor == T){
    print("same q and autocor")
    model = here::here("concise","stan","generalizedpower_sameq_with_auto.stan")
    mod = cmdstan_model(model)
  
    fit <- mod$sample(
      data = data1, 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    
    return(data.frame(fit$summary(c("q","m","n","kappa","bw","bc","bp")), reals = c(q1,m,n,kappa,bw,bc,bp), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
  }else if(q1 == q2 && autocor == F){
    print("same q and no autocor")
    model = here::here("concise","stan","generalizedpower_sameq_without_auto.stan")
    mod = cmdstan_model(model)
    
    fit <- mod$sample(
      data = data1, 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    return(data.frame(fit$summary(c("q","m","n","kappa")), reals = c(q1,m,n,kappa), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
  
  }else if(q1 != q2 && autocor == F){
    print("dif q and no autocor")
    model = here::here("concise","stan","generalizedpower_difq_without_auto.stan")
    mod = cmdstan_model(model)
    
    fit <- mod$sample(
      data = data1, 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    return(data.frame(fit$summary(c("q1","q2","m","n","kappa")), reals = c(q1,q2,m,n,kappa), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
    
  }else if(q1 != q2 && autocor == T){
    print("dif q and autocor")
    model = here::here("concise","stan","generalizedpower_difq_with_auto.stan")
    mod = cmdstan_model(model)
    
    fit <- mod$sample(
      data = data1, 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    return(data.frame(fit$summary(c("q1","q2","m","n","kappa","bw","bc","bp")), reals = c(q1,q2,m,n,kappa,bw,bc,bp), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
    
    
  }else{
      print("idk")
    }

  
  
}



do_parameterrecover_v2 = function(m,n,q1,q2,kappa,autocor,bw,bc,bp, dist){
  index = rnorm(1,-1000,1000)
  
  df = poweragent_generalized_v2(m = m,
                              n = n,
                              q1 = q1,
                              q2 = q2,
                              kappa = kappa,
                              autocor = autocor,
                              bw = bw,
                              bc = bc,
                              bp = bp,
                              dist = dist)
  
  
  data1 = list(Tw1 = df$Tw1,
               Tw2 = df$Tw2,
               Tc1 = df$Tc1,
               Tc2 = df$Tc2, 
               N = nrow(df), 
               w = df$warm, 
               c = df$cold, 
               p = df$pain,
               dist = dist)
  
  if(q1 == q2 && autocor == T){
    print("same q and autocor")
    model = here::here("concise","stan","generalizedpower_sameq_with_newauto.stan")
    mod = cmdstan_model(model)
    
    fit <- mod$sample(
      data = data1, 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    
    return(data.frame(fit$summary(c("q","m","n","kappa","bw","bc","bp")), reals = c(q1,m,n,kappa,bw,bc,bp), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
  }else if(q1 == q2 && autocor == F){
    print("same q and no autocor")
    model = here::here("concise","stan","generalizedpower_sameq_without_auto.stan")
    mod = cmdstan_model(model)
    
    fit <- mod$sample(
      data = data1, 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    return(data.frame(fit$summary(c("q","m","n","kappa")), reals = c(q1,m,n,kappa), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
    
  }else if(q1 != q2 && autocor == F){
    print("dif q and no autocor")
    model = here::here("concise","stan","generalizedpower_difq_without_auto.stan")
    mod = cmdstan_model(model)
    
    fit <- mod$sample(
      data = data1, 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    return(data.frame(fit$summary(c("q1","q2","m","n","kappa")), reals = c(q1,q2,m,n,kappa), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
    
  }else if(q1 != q2 && autocor == T){
    print("dif q and autocor")
    model = here::here("concise","stan","generalizedpower_difq_with_auto.stan")
    mod = cmdstan_model(model)
    
    fit <- mod$sample(
      data = data1, 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    return(data.frame(fit$summary(c("q1","q2","m","n","kappa","bw","bc","bp")), reals = c(q1,q2,m,n,kappa,bw,bc,bp), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
    
    
  }else{
    print("idk")
  }
  
  
  
}



do_parameterrecover_v3 = function(m,n,q1,q2,kappa,autocor,b, dist){
  index = rnorm(1,-1000,1000)
  
  df = poweragent_generalized_v2(m = m,
                              n = n,
                              q1 = q1,
                              q2 = q2,
                              kappa = kappa,
                              autocor = autocor,
                              bw = b,
                              bc = b,
                              bp = 0,
                              dist = dist)
  
  
  data1 = list(Tw1 = df$Tw1,
               Tw2 = df$Tw2,
               Tc1 = df$Tc1,
               Tc2 = df$Tc2, 
               N = nrow(df), 
               w = df$warm, 
               c = df$cold, 
               p = df$pain,
               dist = dist)
    print("same q and autocor")
    model = here::here("concise","stan","generalizedpower_sameq_with_newauto_simp.stan")
    mod = cmdstan_model(model)
    
    fit <- mod$sample(
      data = data1, 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500
    )
    
    return(data.frame(fit$summary(c("q","m","n","kappa","b")), reals = c(q1,m,n,kappa,b), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
  
  
}



do_parameterrecover_doubled = function(m,n,q1,q2,kappa,autocor,b, dist){
  index = rnorm(1,-1000,1000)
  
  df = poweragent_generalized_v2_doubled(m = m,
                              n = n,
                              q1 = q1,
                              q2 = q2,
                              kappa = kappa,
                              autocor = autocor,
                              bw = b,
                              bc = b,
                              bp = 0,
                              dist = dist)
  
  
  data1 = list(Tw1 = df$Tw1,
               Tw2 = df$Tw2,
               Tc1 = df$Tc1,
               Tc2 = df$Tc2, 
               N = nrow(df), 
               w = df$warm, 
               c = df$cold, 
               p = df$pain,
               dist = dist)
  print("same q and autocor")
  model = here::here("concise","stan","generalizedpower_sameq_with_newauto_simp.stan")
  mod = cmdstan_model(model)
  
  fit <- mod$sample(
    data = data1, 
    seed = 123, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500
  )
  
  return(data.frame(fit$summary(c("q","m","n","kappa","b")), reals = c(q1,m,n,kappa,b), div = sum(fit$diagnostic_summary()$num_divergent),index = index))
  
  
}


