modelrecovery = function(m,n,q1,q2,kappa,b, dist){
  index = rnorm(1,-1000,1000)
  qq1 = data.frame()
  df_auto = poweragent_generalized_v2_simp(m = m,
                                 n = n,
                                 q1 = q1,
                                 q2 = q2,
                                 kappa = kappa,
                                 autocor = TRUE,
                                 b = b,
                                 dist = dist)
  
  
  
  df_noauto = poweragent_generalized_v2_simp(m = m,
                                      n = n,
                                      q1 = q1,
                                      q2 = q2,
                                      kappa = kappa,
                                      autocor = FALSE,
                                      b = 0,
                                      dist = dist)
  
  
  
  
  data_auto = list(Tw1 = df_auto$Tw1,
               Tw2 = df_auto$Tw2,
               Tc1 = df_auto$Tc1,
               Tc2 = df_auto$Tc2, 
               N = nrow(df_auto), 
               w = df_auto$warm, 
               c = df_auto$cold, 
               p = df_auto$pain,
               dist = dist)
  
  
  
  data_noauto = list(Tw1 = df_noauto$Tw1,
                   Tw2 = df_noauto$Tw2,
                   Tc1 = df_noauto$Tc1,
                   Tc2 = df_noauto$Tc2, 
                   N = nrow(df_noauto), 
                   w = df_noauto$warm, 
                   c = df_noauto$cold, 
                   p = df_noauto$pain,
                   dist = dist)
  
  
  
  model = here::here("concise","stan","generalizedpower_sameq_with_newauto_simp.stan")
  mod_auto = cmdstan_model(model)
  
  fit_auto_auto <- mod_auto$sample(
    data = data_auto, 
    seed = 123, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500
  )
  
  fit_auto_noauto <- mod_auto$sample(
    data = data_noauto, 
    seed = 123, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500
  )
  loopower_auto = loo(fit_auto_auto$draws(c("log_lik")), r_eff = relative_eff(exp(fit_auto_auto$draws(c("log_lik")))))
  loopower_noauto = loo(fit_auto_noauto$draws(c("log_lik")), r_eff = relative_eff(exp(fit_auto_noauto$draws(c("log_lik")))))
  
  loo = loo::loo_compare(list(loopower_auto, loopower_noauto))
  
  
  qq = data.frame(loo = loo[1:2], id = index)
  qq$models = rownames(qq)
  qq1 = rbind(qq1,qq)
  
  
  model = here::here("concise","stan","generalizedpower_sameq_with_newnoauto_simp.stan")
  mod_noauto = cmdstan_model(model)
  
  fit_noauto_noauto <- mod_noauto$sample(
    data = data_noauto, 
    seed = 123, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500
  )
  
  fit_noauto_auto <- mod_noauto$sample(
    data = data_auto, 
    seed = 123, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500
  )
  loopower_noauto = loo(fit_noauto_noauto$draws(c("log_lik")), r_eff = relative_eff(exp(fit_noauto_noauto$draws(c("log_lik")))))
  
  loopower_auto = loo(fit_noauto_auto$draws(c("log_lik")), r_eff = relative_eff(exp(fit_noauto_auto$draws(c("log_lik")))))

  weights = loo::loo_model_weights(list(loopower_noauto, loopower_auto))
  
  
  qq = data.frame(weightss = weights[1:2], id = index)
  qq$models = rownames(qq)
  
  qq1 = cbind(qq1,qq)
  
  
  return(qq1)
  
  
}