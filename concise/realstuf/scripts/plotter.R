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
      facet_wrap(~variable, scales = "free", nrow = 1)+theme_classic()+
      geom_abline(slope = 1, intercept = 0)+
      scale_y_continuous(breaks=pretty_breaks(n = 4), labels = comma)+scale_x_continuous(breaks=pretty_breaks(n = 4), labels = comma)
    
    
    kappa = ggplot(data.frame(y = rprop(10000, i, 0.5)),aes(y = y))+geom_histogram()+theme_classic()+ylab(" ")+xlab(" ")+
      scale_y_continuous(breaks=pretty_breaks(n = 4), labels = comma)+scale_x_continuous(breaks=pretty_breaks(n = 4), labels = comma)
    
    layout <- c(
      area(1, 1,1,3),
      area(1, 4)
    )
    
    
    p = pr + kappa + plot_layout(design = layout)
    
    row_plots[[q]] = p
    
    }
  
  plots_list <- lapply(row_plots, function(p) p)
  
  
  plot = wrap_plots(
    plots_list,
    nrow = length(row_plots),
    byrow = TRUE
  )
  
  plot
  
  
  
}