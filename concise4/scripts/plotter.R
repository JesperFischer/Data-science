#plotter


get_pr_plot = function(parameterrecover){
  kappas = parameterrecover %>% filter(str_detect(variable, "kappa"))
  values = unique(kappas$reals)
  q = 0
  row_plots = list()
  for(i in values){
    q = q+1
    print(i)
    pr = parameterrecover %>% filter(index %in% parameterrecover[parameterrecover[,1] == "kappa" & parameterrecover[,11] == i,13]) %>% 
      filter(!str_detect(variable, "kappa")) %>% 
      ggplot(aes(x = mean, y = reals))+
      geom_point()+
      facet_wrap(~variable, scales = "free", nrow = 1)+theme_classic()+geom_abline(slope = 1, intercept = 0)
    
    
    kappa = ggplot(data.frame(),aes())+geom_histogram(aes(y = rprop(10000, i, 0.5)))+theme_classic()+ylab(" ")+xlab(" ")+
      scale_y_continuous(breaks = seq(0,1,by = 0.1), labels  = seq(0,1,by = 0.1))
    
    layout <- c(
      area(1, 1,1,3),
      area(1, 4)
    )
    
    
    p = pr + kappa + plot_layout(design = layout)
    
    row_plots[[q]] = p
  }
  
  plot = wrap_plots(
    row_plots[[1]],
    row_plots[[2]],
    row_plots[[3]],
    row_plots[[4]],
    row_plots[[5]],
    row_plots[[6]],
    row_plots[[7]],
    nrow = length(row_plots),
    byrow = TRUE
  )
  
  plot
  
  
  
}