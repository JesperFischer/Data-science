

poweragent_gen = function(m,n,q1,q2,kw,kc,kappa, autocor = FALSE,bw,bc,bp, dist = TRUE){
  
  powerlaw_w1 = function(x, q, n, kw){
    return(kw * (x+q)^n)
  }
  powerlaw_w1_int= function(x, q, n, kw){
    return(kw * (((x+q)^(n+1))/(n+1)))
  }
  
  powerlaw_w2= function(x,  q,  n, kw){
    return(kw * (-x+q)^n)
  }
  
  powerlaw_w2_int= function(x, q, n, kw){
    return(-kw * (((-x+q)^(n+1))/(n+1)))
  }
  
  
  powerlaw_c1 = function(x,  q, m, kc){
    return(-kc * (x+q)^m)
  }
  
  powerlaw_c1_int= function( x, q, m, kc){
    return(-kc * (((x+q)^(m+1))/(m+1)))
  }
  
  powerlaw_c2 = function(x, q, m, kc){
    return(-kc * (-x+q)^m)
  }
  
  powerlaw_c2_int= function(x, q, m, kc){
    return(kc * (((-x+q)^(m+1))/(m+1)))
  }
  
  linear = function(xs, x1, x2, y1, y2){
    a = ((y2-y1)/(x2-x1))
    return(a*xs+y1-(a*x1))
    
  }
  
  linear_int = function(xs, x1, x2, y1, y2){
    a = (y2-y1)/(x2-x1)
    b = y1-a*x1
    return(((a*xs^2)/2)+(b*xs))
  }
  
  
  get_points = function(Tw1 , Tw2 , Tc1 , Tc2 , q1, q2 , n , m, kw, kc){
    #both warm
    if(Tw1>Tc1){
      x1 = Tw1-q1
      x2 = q2-Tw2
      y1 = powerlaw_w1(x1,q1,n, kw)
      y2 = powerlaw_w2(x2,q2,n, kw)
      #both cold
    }else if(Tc1>0 && Tc2 >0){
      x1 = Tc1-q1
      x2 = q2-Tc2
      y1 = powerlaw_c1(x1,q1,m, kc)
      y2 = powerlaw_c2(x2,q2,m, kc)
      #1 cold 2 warm
    }else if(Tc1>0 && Tw2 >0){
      x1 = Tc1-q1
      x2 = q2-Tw2
      y1 = powerlaw_c1(x1,q1,m, kc)
      y2 = powerlaw_w2(x2,q2,n, kw)
      #second neutral and first cold
    }else if(Tw2 == 0 && Tc2 == 0 && Tc1>0){
      x1 = Tc1-q1
      x2 = q2
      y1 = powerlaw_c1(x1,q1,m, kc)
      y2 = 0
      #first neutral and therefore other warm (given the above takes care if its neutral)
    }else{
      x1 = -q1
      x2 = q2-Tw2
      y1 = 0
      y2 = powerlaw_w2(x2,q2,n, kw)
    }
    
    points = c("x1" = x1,"x2" = x2,"y1" = y1, "y2" = y2)
    
    return(data.frame(t(points)))
  }
  
  get_end_areas = function(Tw1,  Tw2, Tc1, Tc2, q1, q2, n, m, x1,x2, kw, kc){
    
    #being warm for the first thermode: we integrate from -q to x1
    if(Tw1>Tc1){
      a1_w = powerlaw_w1_int(x1, q1, n, kw)-powerlaw_w1_int(-q1, q1, n, kw)
      a1_c = 0
      #cold for the first thermode: we integrate from -q to x1 and
    }else if(Tw1<Tc1){
      a1_c = -(powerlaw_c1_int(x1, q1, m, kc)-powerlaw_c1_int(-q1, q1, m, kc))
      a1_w = 0
      #if they are equal then it must mean that they are both baseline
    }else{
      a1_w = 0
      a1_c = 0
    }
    #second thermode warm: integrate from x2 to q
    if(Tw2>Tc2){
      a2_w = powerlaw_w2_int(q2,q2,n, kw)-powerlaw_w2_int(x2,q2,n, kw)
      a2_c = 0
      #second thermode cold integrate from x2 to q and -
    }else if(Tw2<Tc2){
      a2_c = -(powerlaw_c2_int(q2,q2,m, kc)-powerlaw_c2_int(x2,q2,m, kc))
      a2_w = 0
      #if they are equal then it must mean that they are both baseline
    }else{
      a2_w = 0
      a2_c = 0
    }
    areas = c("a1_w" = a1_w, "a2_w" = a2_w, "a1_c" = a1_c, "a2_c" = a2_c)
    
    return(data.frame(t(areas)))
    
  }
  
  get_middle_areas = function( Tw1,  Tw2,  Tc1,  Tc2, q1, q2,  n,  m, x1 , x2 , y1 , y2){
    
    if(Tw2== 0 && Tw1 == 0 && Tc2 == 0 && Tc1 == 0){
      a3_w = 0
      a3_c = 0
      #if both are warm
    }else if(Tw2>0 && Tw1 > 0){
      #if they are not equal i.e. a vertical line:
      if(Tw2 != Tw1){
        #integrate linear from x1 to x2
        a3_w = linear_int(x2,x1,x2,y1,y2)-linear_int(x1,x1,x2,y1,y2)
        a3_c = 0
      }else{
        #if they are equal then the slope is = 0 and its just a square
        a3_w = (-x1+x2)*(y1)
        a3_c = 0
      }
      # if both are cold
    }else if(Tc2>0 && Tc1 > 0){
      if(Tc2 != Tc1){
        a3_c = -(linear_int(x2,x1,x2,y1,y2)-linear_int(x1,x1,x2,y1,y2))
        a3_w = 0
      }else{
        a3_c = -(-x1+x2)*y1
        a3_w = 0
      }
      # if 1 is cold and 2 is warm
    }else if(Tc1>=0 && Tw2>=0){
      
      a = (y2-y1)/(x2-x1)
      xint = ((a*x1)-y1)/a
      
      a3_w = linear_int(x2,x1,x2,y1,y2)-linear_int(xint,x1,x2,y1,y2)
      a3_c = -(linear_int(xint,x1,x2,y1,y2)-linear_int(x1,x1,x2,y1,y2))
    }else{
      a3_w = 0
      a3_c = 0
    }
    
    areas3 = c("a3_w" = a3_w, "a3_c" = a3_c)
    
    return(data.frame(t(areas3)))
    
  }
  
  betaprop <- function(kappa, mean) {
    value = rprop(1, kappa, mean)
    return(data.frame(value))
  }
  
  
  gaus = function(kappa, mean) {
    value = rnorm(1, mean, 1/kappa^2)
    return(data.frame(value))
    
  }
  
  source(here("concise","scripts","utility.R"))
  
  temp = get_temps(ntrials = 3)
  
  Tw1 = ifelse(temp$T1 <= 30, 0 ,ifelse(temp$T1 > 30, temp$T1-30, NA))
  Tc1 = ifelse(temp$T1 >= 30, 0 ,ifelse(temp$T1 < 30, 30-temp$T1, NA))
  Tw2 = ifelse(temp$T2 <= 30, 0 ,ifelse(temp$T2 > 30, temp$T2-30, NA))
  Tc2 = ifelse(temp$T2 >= 30, 0 ,ifelse(temp$T2 < 30, 30-temp$T2, NA))
  
  perm <- sample(length(Tw1))
  
  Tw1 = Tw1[perm]
  Tw2 = Tw2[perm]
  Tc1 = Tc1[perm]
  Tc2 = Tc2[perm]
  
  
  
  Atotw = (powerlaw_w1_int(0,q1,n, kw)-powerlaw_w1_int(-q1,q1,n, kw))+(powerlaw_w2_int(q2,q2,n, kw)-powerlaw_w2_int(0,q2,n, kw))
  Atotc = -((powerlaw_c1_int(0,q1,m, kc)-powerlaw_c1_int(-q1,q1,m, kc))   +   (powerlaw_c2_int(q2,q2,m, kc)-powerlaw_c2_int(0,q2,m, kc)))
  
  data = data.frame(Tw1 = Tw1, Tw2 = Tw2, Tc1 = Tc1, Tc2 = Tc2, q1 = q1,q2 = q2, n = n, m = m, kw = kw, kc = kc)
  
  
  data = pmap(data, get_points) %>% bind_rows() %>% bind_cols(data, .)
  
  data = data %>% dplyr::select(-c("y1","y2")) %>% pmap(get_end_areas) %>% bind_rows() %>% bind_cols(data, .)
  
  data = data %>% dplyr::select(-c("a1_w","a2_w","a1_c","a2_c","kw","kc")) %>% pmap(get_middle_areas) %>% bind_rows() %>% bind_cols(data, .)
  
  
  if (dist == TRUE){
    data = data %>% mutate(Atotw = Atotw,
                           Atotc = Atotc,
                           uw = ((a1_w+a2_w+a3_w)/Atotw)+0.001,
                           uc = ((a1_c+a2_c+a3_c)/Atotc)+0.001,
                           up = uc*uw+0.001)
    
    if(autocor == TRUE){
      data$uw = ifelse((lag(data$uw, 1, default = 0)+(1/(1+bw))*(data$uw-lag(data$uw, 1, default = 0))) <= 0 , 0.01,
                ifelse((lag(data$uw, 1, default = 0)+(1/(1+bw))*(data$uw-lag(data$uw, 1, default = 0))) >= 1 , 0.99,
                       lag(data$uw, 1, default = 0)+(1/(1+bw))*(data$uw-lag(data$uw, 1, default = 0))))
      
      
      data$uc = ifelse((lag(data$uc, 1, default = 0)+(1/(1+bc))*(data$uc-lag(data$uc, 1, default = 0))) <= 0, 0.01,
                ifelse((lag(data$uc, 1, default = 0)+(1/(1+bc))*(data$uc-lag(data$uc, 1, default = 0))) >= 1, 0.99,       
                       lag(data$uc, 1, default = 0)+(1/(1+bc))*(data$uc-lag(data$uc, 1, default = 0))))
      
      
      data$up = ifelse((lag(data$up, 1, default = 0)+(1/(1+bp))*(data$up-lag(data$up, 1, default = 0))) <= 0, 0.01,
                ifelse((lag(data$up, 1, default = 0)+(1/(1+bp))*(data$up-lag(data$up, 1, default = 0))) >= 1, 0.99,
                       lag(data$up, 1, default = 0)+(1/(1+bp))*(data$up-lag(data$up, 1, default = 0))))
      
    }
    
    data = data %>% mutate(kappa = kappa) %>% dplyr::select("kappa","uw") %>% rename(mean = uw) %>% pmap(betaprop) %>% bind_rows() %>% rename(warm = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% dplyr::select("kappa","uc") %>% rename(mean = uc) %>% pmap(betaprop) %>% bind_rows() %>% rename(cold = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% dplyr::select("kappa","up") %>% rename(mean = up) %>% pmap(betaprop) %>% bind_rows() %>% rename(pain = value) %>% bind_cols(data, .)
    
    data$T1 = temp$T2
    data$T2 = temp$T1
    
    
    return(data)
    
  }
  
  
  if (dist == FALSE){
    data = data %>% mutate(Atotw = Atotw,
                           Atotc = Atotc,
                           uw = ((a1_w+a2_w+a3_w)/Atotw)+0.01,
                           uc = ((a1_c+a2_c+a3_c)/Atotc)+0.01,
                           up = uc*uw)
    
    if(autocor == TRUE){
      data$uw = lag(data$uw, 1, default = 0)+(1/(1+bw))*(data$uw-lag(data$uw, 1, default = 0))
      
      data$uc = lag(data$uc, 1, default = 0)+(1/(1+bc))*(data$uc-lag(data$uc, 1, default = 0))
      
      data$up = lag(data$up, 1, default = 0)+(1/(1+bp))*(data$up-lag(data$up, 1, default = 0))
      
    }
    
    data = data %>% mutate(kappa = kappa) %>% dplyr::select("kappa","uw") %>% rename(mean = uw) %>% pmap(gaus) %>% bind_rows() %>% rename(warm = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% dplyr::select("kappa","uc") %>% rename(mean = uc) %>% pmap(gaus) %>% bind_rows() %>% rename(cold = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% dplyr::select("kappa","up") %>% rename(mean = up) %>% pmap(gaus) %>% bind_rows() %>% rename(pain = value) %>% bind_cols(data, .)
    
    data$T1 = temp$T2
    data$T2 = temp$T1
    
    return(data)
    
  }
  
  
  
}









model1 = function(h1,h2,w1,w2, sigma){ 

  source(here("concise","scripts","utility.R"))
  
  temp = get_temps(ntrials = 3)
  
  ntrials = length(temp$T1)
  
  
  Tw1 = ifelse(temp$T1 <= 30, 0 ,ifelse(temp$T1 > 30, temp$T1-30, NA))
  Tc1 = ifelse(temp$T1 >= 30, 0 ,ifelse(temp$T1 < 30, 30-temp$T1, NA))
  Tw2 = ifelse(temp$T2 <= 30, 0 ,ifelse(temp$T2 > 30, temp$T2-30, NA))
  Tc2 = ifelse(temp$T2 >= 30, 0 ,ifelse(temp$T2 < 30, 30-temp$T2, NA))
  
  
  Tw1_t = sqrt(w1^2+h1^2)
  Tw2_t = sqrt(w2^2+h1^2)
  Tc1_t = sqrt(w2^2+h1^2)
  Tc2_t = sqrt(w2^2+h2^2)
  
  
  arc_length_1 = array(NA, ntrials)
  arc_length_2 = array(NA, ntrials)
  arc_length_2_fake= array(NA, ntrials)
  x1 = array(NA, ntrials)
  x2 = array(NA, ntrials)
  y1 = array(NA, ntrials)
  y2= array(NA, ntrials)
  a= array(NA, ntrials)
  b= array(NA, ntrials)
  Atot = array(NA, ntrials)
  Awarm = array(NA, ntrials)
  Apain = array(NA, ntrials)
  area_w = array(NA, ntrials)
  Acold = array(NA, ntrials)
  area_c = array(NA, ntrials)
  xinter = array(NA,ntrials)
  alpha_w= array(NA,ntrials)
  alpha_c= array(NA,ntrials)
  beta_w= array(NA,ntrials)
  beta_c= array(NA,ntrials)
  alpha_p= array(NA,ntrials)
  beta_p= array(NA,ntrials)
  
  c= array(NA,ntrials)
  w= array(NA,ntrials)
  p = array(NA,ntrials)
  area_c_1 = array(NA,ntrials)
  area_c_2= array(NA,ntrials)
  area_c_3= array(NA,ntrials)
  area_w_1 = array(NA,ntrials)
  area_w_2= array(NA,ntrials)
  area_w_3= array(NA,ntrials)
  
  
  y1_w = function(xw1){
    return((h1/(w1))*xw1+h1)
  }
  y1_c = function(xc1){
    return((-h2/w1)*xc1-h2)
  }
  
  y2_w = function(xw2){
    return((-h1/w2)*xw2+h1)
  }
  y2_c = function(xc2){
    return((h2/w2)*xc2-h2)
  }
  
  
  for (i in 1:ntrials){
    
    arc_length_1[i] = max(Tc1[i],Tw1[i])
    arc_length_2[i] = max(Tc2[i],Tw2[i])
    #need to redefine the arclengths for the right side of the triangle to account for the fact that we go from the x-axis and up or downwards and not from the y-axis
    
    arc_length_2[i] = ifelse(Tw2[i] > Tc2[i],Tw2_t-arc_length_2[i], Tc2_t-arc_length_2[i])
    
    
    
    x1[i] = ifelse(Tw1[i] > Tc1[i], (arc_length_1[i]/sqrt(1+(h1/w1)^2))-w1, (arc_length_1[i]/sqrt(1+(h2/w1)^2))-w1)
    
    
    
    x2[i] = ifelse(Tw2[i] > Tc2[i],(arc_length_2[i]/sqrt(1+(-h1/w2)^2)), (arc_length_2[i]/sqrt(1+(h2/w2)^2)))
    
    
    
    
    #line connecting the two points:
    
    
    y1[i] <- ifelse(Tc1[i] > Tw1[i], y1_c(x1[i]), y1_w(x1[i]))
    y2[i] <- ifelse(Tc2[i] > Tw2[i], y2_c(x2[i]), y2_w(x2[i]))
    
    
    a[i] <- ifelse(y2[i] != y1[i], (y2[i]-y1[i])/(x2[i]-x1[i]),-100000)
    xinter[i] = ifelse(y2[i] != y1[i], -(y1[i]-a[i]*x1[i])/a[i],-100000)
    b[i] <- ifelse(y2[i] != y1[i], y1[i]-a[i]*x1[i],y1[i])
    
    
    ycon1 = function(x,a,b){
      return(a*x+b)   
    }
    
    
    
    if(Tw2[i]>Tc2[i] & Tc1[i] > Tw1[i]){
      
      area_c[i] = integrate(y1_c, lower = -w1, upper = x1[i])[[1]]+
        integrate(ycon1, lower = x1[i], upper = xinter[i], a = a[i], b = b[i])[[1]]
      
      area_w[i] = integrate(ycon1, lower = xinter[i], upper = x2[i], a = a[i], b = b[i])[[1]]+
        integrate(y2_w, lower = x2[i], upper = w2)[[1]]
    }
    if(Tc2[i] > Tw2[i] | Tc2[i] == Tw2[i]){
      area_c[i] = integrate(y1_c, lower = -w1, upper = x1[i])[[1]]+
        integrate(ycon1, lower = x1[i], upper = x2[i], a = a[i], b = b[i])[[1]]+
        integrate(y2_c, lower = x2[i], upper = w2)[[1]]
      area_w[i] = 0
    }
    if(Tw1[i] > Tc1[i] | Tw1[i] == Tc1[i]){
      area_w[i] = integrate(y1_w, lower = -w1, upper = x1[i])[[1]]+
        integrate(ycon1, lower = x1[i], upper = x2[i], a = a[i], b = b[i])[[1]]+
        integrate(y2_w, lower = x2[i], upper = w2)[[1]]
      
      area_c[i] = 0
    }
    
    if(Tw1[i] == Tc1[i] & Tc2[i] == Tw2[i]){
      area_w[i] = 0
      area_c[i] = 0
    }
    
    
    Atot_w = ((w1+w2)*h1)/2
    Atot_c = ((w1+w2)*h2)/2
    
    
    Awarm[i] = area_w[i]/Atot_w
    Acold[i] = -area_c[i]/Atot_c
    
    Apain[i] = (-area_c[i]*area_w[i])/(Atot_w*Atot_c)
    
    kappa = 1/sigma
    
    Awarm[i] = Awarm[i]+0.01
    Acold[i] = Acold[i]+0.01
    Apain[i] = Apain[i]+0.01
    
    
    alpha_w[i] = Awarm[i]*(kappa-1)
    beta_w[i] = (kappa-1)*(1-Awarm[i])
    
    w[i] <- rbeta(1,alpha_w[i],beta_w[i])
    
    alpha_c[i] = Acold[i]*(kappa-1)
    beta_c[i] = (kappa-1)*(1-Acold[i])
    
    c[i] <- rbeta(1,alpha_c[i],beta_c[i])
    
    alpha_p[i] = Apain[i]*(kappa-1)
    beta_p[i] = (kappa-1)*(1-Apain[i])
    
    p[i] <- rbeta(1,alpha_p[i],beta_p[i])
    
    
    
    
    
    
  }
  return(data.frame(Atot_w = Atot_w, Atot_c = Atot_c, w = w,c =  c, p = p, Acold = Acold, Awarm = Awarm, Apain = Apain,temp1 = temp$T1,temp2 = temp$T2, Tc1 = Tc1, Tw1 = Tw1,Tc2 = Tc2,Tw2 = Tw2, area_w = area_w, area_c = area_c))
}


