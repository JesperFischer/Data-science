#Stan functions


get_diag = function(fit, parameter, prior = FALSE){
  
  parameters = grep(parameter, dimnames(fit$draws())$variable, value = TRUE)
  
  priors = grep("prior", parameters, value = TRUE)
  
  posteriors = parameters[!parameters %in% priors]
  
  if(!prior){
    
    nsubs = length(posteriors)
    nplots = ceiling(nsubs/10)
    
    trace_plots = list()
    pairs_plots = list()
    
    for(i in 1:nplots){
      
      if(i != nplots){
        
        indicies = c(((i - 1) * 10 + 1):(i * 10),nsubs-1,nsubs)
        
        trace = mcmc_trace(fit$draws(variables = posteriors[indicies]))+
          theme(strip.text = element_text(size = 18), axis.text = element_text(size = 18),axis.title = element_text(size = 18))
        
        
        pairs = mcmc_pairs(fit$draws(variables = posteriors[indicies]), np = nuts_params(fit), pars = posteriors[indicies],
                   off_diag_args = list(size = 0.75))
        
        trace_plots[[i]] = trace
        pairs_plots[[i]] = pairs
        
      }else{
        
        
        trace = mcmc_trace(fit$draws(variables = posteriors[((i-1) * 10):nsubs]))+
          theme(strip.text = element_text(size = 18), axis.text = element_text(size = 18),axis.title = element_text(size = 18))
        
        
        pairs = mcmc_pairs(fit$draws(variables = posteriors[((i-1) * 10):nsubs]), np = nuts_params(fit), pars = posteriors[((i-1) * 10):nsubs],
                           off_diag_args = list(size = 0.75))
        
        trace_plots[[i]] = trace
        pairs_plots[[i]] = pairs
      }
    }
    
  
    return(list(trace_plots = trace_plots, pairs_plots = pairs_plots))
  }
}

ppu_hier = function(fit,parameters,reals){
  
  #things to search for 
  search = c("mu","kappa","sd")
  
  post_parameters = paste0(search,"_",parameters)
  
  post_parameters = post_parameters[post_parameters %in% dimnames(fit$draws())$variable]
  
  prior_parameters = paste0("prior_",search,"_",parameters)
  
  prior_parameters = prior_parameters[prior_parameters %in% dimnames(fit$draws())$variable]
  
  
  post = as_draws_df(fit$draws(variables = post_parameters)) %>% 
    mutate(posterior = T)
  
  priors = as_draws_df(fit$draws(variables = prior_parameters)) %>% 
    mutate(posterior = F) %>%
    rename_with(~sub("^prior_", "", .), starts_with("prior_"))
  
  
  df = rbind(post,priors)
  
  df %>% pivot_longer(cols = post_parameters)%>% 
    ggplot(aes(x = value, fill = posterior))+
    geom_histogram(alpha = 0.5, position="identity")+
    facet_wrap(~name, scales = "free")+
    theme_classic()+
    geom_vline(data = reals %>% pivot_longer(cols = post_parameters), aes(xintercept = value))
  
  
}




ppu_sub = function(fit,parameter,reals, lim = 10){
  
  reals = reals %>% select(parameter,"id") %>% group_by(id) %>% summarize(mean = mean(!!sym(parameter)))
  
  #get rid of these:
  search = c("mu","kappa","sd")
  #get the names of the parameters
  post_parameters = grep(parameter, dimnames(fit$draws())$variable, value = TRUE)
  #find where the search term is.
  result <- sapply(search, function(x) grepl(x, post_parameters, fixed = TRUE))
  #find where any of them are
  result <- apply(result, 1, any)
  #get rid of them
  parameters = post_parameters[!result]
  
  
  prior_parameters = grep("prior",parameters, value = TRUE)
  
  post_parameters = parameters[!parameters %in% prior_parameters]
  
  post = as_draws_df(fit$draws(variables = post_parameters)) %>%
    mutate(posterior = T)
  
  priors = as_draws_df(fit$draws(variables = prior_parameters)) %>%
    mutate(posterior = F) %>%
    rename_with(~sub("^prior_", "", .), starts_with("prior_"))
  
  
  df = rbind(post,priors)
  
  df %>% pivot_longer(cols = post_parameters) %>% 
    ggplot(aes(x = value, fill = posterior))+
    geom_histogram(alpha = 0.5, position="identity")+
    facet_wrap(~name, nrow = length(parameters)/2, ncol = 5, scales = "free")+
    geom_vline(data = reals %>% mutate(name = paste0(parameter,"[", row_number(), "]")), aes(xintercept = mean))+
    theme_classic()+
    scale_x_continuous(limits = c(0,lim), breaks = scales::pretty_breaks(n = 5))
  
}




