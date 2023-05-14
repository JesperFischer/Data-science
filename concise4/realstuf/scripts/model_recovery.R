

model_fitter = function(data){
  

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
  
  brmsdata = data.frame(temp1 = df$T1, 
                        temp2 = df$T2, 
                        cold = df$cold, 
                        warm = df$warm, 
                        pain = df$pain)
  
  
  loo1 = get_loo(here::here("concise","realstuf","stan","diamond_simplest.stan"), data1)
  loo2 = get_loo(here::here("concise","realstuf","stan","real_diamond.stan"), data1)
  loo3 = get_loo(here::here("concise","realstuf","stan","real_diamond_with_auto_on_thermo.stan"), data1)
  loo4 = get_loo(here::here("concise","realstuf","stan","powerlaw_fixed_k.stan"), data1)
  loo5 = get_loo(here::here("concise","realstuf","stan","powerlaw_fixed_k_with_auto_on_thermo.stan"), data1)
  loo6 = get_loo(here::here("concise","realstuf","stan","powerlaw_fixed_k_with_auto_general.stan"), data1)
  loo7 = get_loo(here::here("concise","realstuf","stan","full_power(sameq).stan"), data1)
  loo8 = get_loo(here::here("concise","realstuf","stan","full_power(sameq)_fix_bp.stan"), data1)
  loo9 = get_brms_loo(here::here("concise","realstuf","stan","brms_ind.stan"),brmsdata)
  
  comparison = loo_compare(list(loo1[[3]],loo2[[3]],loo3[[3]],loo4[[3]],loo5[[3]],loo6[[3]],loo7[[3]], loo8[[3]], loo9[[3]]))
  
  weights = loo_model_weights(list(loo1[[3]],loo2[[3]],loo3[[3]],loo4[[3]],loo5[[3]],loo6[[3]],loo7[[3]], loo8[[3]], loo9[[3]]))
  
  diag = rbind(loo1[[4]],loo2[[4]],loo3[[4]],loo4[[4]],loo5[[4]],loo6[[4]],loo7[[4]],loo8[[4]],loo9[[4]])
  
  return(list(data.frame(names = c("simp_diamond",
             "real_diamond",
             "real_diamond_auto_on_thermo",
             "powerlaw_fixed_k",
             "powerlaw_fixed_k_auto_on_thermo",
             "powerlaw_fixed_k_full_auto",
             "fullpower_sameq",
             "full_power(sameq)_fix_bp" ,
             "brms"),
             means = c(loo1[[1]],loo2[[1]],loo3[[1]],loo4[[1]],loo5[[1]],loo6[[1]],loo7[[1]], loo8[[1]],loo9[[1]]),
             sds = c(loo1[[2]],loo2[[2]],loo3[[2]],loo4[[2]],loo5[[2]],loo6[[2]],loo7[[2]], loo8[[2]],loo9[[2]])),
             comparison,
             weights,
             diag))
}

get_loo = function(file, data1){
  mod = cmdstan_model(file)
  
    fit <- tryCatch({
      fit = withTimeout(fitter(mod, data1), timeout = 200, cpu=200)
    }, error = function(e) {
      print("An error occurred!")
      return(NA)
    })
    
  if(is.environment(fit)){
    q = fit$loo()
    loo1 = q$estimates[3,1]
    loo1sd = q$estimates[3,2]
    diag = as.data.frame(fit$diagnostic_summary()) %>% mutate(id = file)
    
    return(list(loo1,loo1sd,q, diag))
  }else{
    
    diag = data.frame(array(NA,c(4,4)))
    names(diag) = c("num_divergent","num_max_treedepth","ebfmi","id")
    diag$id = file
    
    return(list(NA,NA,diag))
  }
  
  }




fitter = function(mod, data1){
  fit <- mod$sample(
    data = data1,
    chains = 4, 
    parallel_chains = 4,
    refresh = 500
  )
  return(fit)
}


get_brms_loo = function(file, data1){
  
  mod = cmdstan_model(file)
  
  formula = brms::bf(mvbind(cold,warm,pain) ~ temp1+temp2,
                     zoi ~ temp1+temp2,
                     family = zero_one_inflated_beta())
  
  
  standata = make_standata(formula, data1)
  
  fit <- mod$sample(
    data = standata, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500
  )
  
  q = fit$loo()
  loo1 = q$estimates[3,1]
  loo1sd = q$estimates[3,2]
  diag = as.data.frame(fit$diagnostic_summary()) %>% mutate(id = file)
  
  
  return(list(loo1,loo1sd,q,diag))
  
}
