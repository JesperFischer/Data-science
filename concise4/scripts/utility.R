get_temps = function(ntrials){#temperatures on thermode 1 - to match experiment
  T1 <- rep(c(seq(18,42,4),
              seq(18,38,4),
              seq(18,34,4),
              seq(18,30,4),
              seq(18,26,4),
              seq(18,22,4),
              seq(18,18,4)),ntrials)
  
  #temperatures on thermode 2
  T2 <- rep(c(rep(42,7),
              rep(38,6),
              rep(34,5),
              rep(30,4),
              rep(26,3),
              rep(22,2),
              rep(18,1)),ntrials)
  
  ntrials <- length(T1)
  
  
  Tc2 = ifelse(T2-30<0,abs(T2-30),0) 
  Tw2 = ifelse(T2-30>0,T2-30,0) 
  Tc1 = ifelse(T1-30<0,abs(T1-30),0) 
  Tw1 = ifelse(T1-30>0,T1-30,0) 
  
  #temps = list(Tc1 = Tc1,Tw1 = Tw1, Tw2 = Tw2, Tc2 = Tc2, temp1 = T1, temp2 = T2)
  temps = list(T1 = T1, T2 = T2)
  
  return(temps)
  
}
