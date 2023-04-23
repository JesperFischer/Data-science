#agent

poweragent_generalized = function(m,n,q1,q2,kappa, autocor = FALSE,bw = 0,bc = 0,bp = 0, dist = TRUE){
  
  powerlaw_w1 = function(x, q, n){
    return((x+q)^n)
  }
  powerlaw_w1_int= function(x, q, n){
    return((((x+q)^(n+1))/(n+1)))
  }
  
  powerlaw_w2= function(x,  q,  n){
    return((-x+q)^n)
  }
  
  powerlaw_w2_int= function(x, q, n){
    return(-(((-x+q)^(n+1))/(n+1)))
  }
  
  
  powerlaw_c1 = function(x,  q, m){
    return(-(x+q)^m)
  }
  
  powerlaw_c1_int= function( x, q, m){
    return(-(((x+q)^(m+1))/(m+1)))
  }

  powerlaw_c2 = function(x, q, m){
    return(-(-x+q)^m)
  }
  
  powerlaw_c2_int= function(x, q, m){
    return((((-x+q)^(m+1))/(m+1)))
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
  
  
  get_points = function(Tw1 , Tw2 , Tc1 , Tc2 , q1, q2 , n , m){
    #both warm
    if(Tw1>Tc1){
      x1 = Tw1-q1
      x2 = q2-Tw2
      y1 = powerlaw_w1(x1,q1,n)
      y2 = powerlaw_w2(x2,q2,n)
      #both cold
    }else if(Tc1>0 && Tc2 >0){
      x1 = Tc1-q1
      x2 = q2-Tc2
      y1 = powerlaw_c1(x1,q1,m)
      y2 = powerlaw_c2(x2,q2,m)
      #1 cold 2 warm
    }else if(Tc1>0 && Tw2 >0){
      x1 = Tc1-q1
      x2 = q2-Tw2
      y1 = powerlaw_c1(x1,q1,m)
      y2 = powerlaw_w2(x2,q2,n)
      #second neutral and first cold
    }else if(Tw2 == 0 && Tc2 == 0 && Tc1>0){
      x1 = Tc1-q1
      x2 = q2
      y1 = powerlaw_c1(x1,q1,m)
      y2 = 0
      #first neutral and therefore other warm (given the above takes care if its neutral)
    }else{
      x1 = -q1
      x2 = q2-Tw2
      y1 = 0
      y2 = powerlaw_w2(x2,q2,n)
    }
    
    points = c("x1" = x1,"x2" = x2,"y1" = y1, "y2" = y2)
    
    return(data.frame(t(points)))
  }
  
  get_end_areas = function(Tw1,  Tw2, Tc1, Tc2, q1, q2, n, m, x1,x2){
    
    #being warm for the first thermode: we integrate from -q to x1
    if(Tw1>Tc1){
      a1_w = powerlaw_w1_int(x1, q1, n)-powerlaw_w1_int(-q1, q1, n)
      a1_c = 0
      #cold for the first thermode: we integrate from -q to x1 and
    }else if(Tw1<Tc1){
      a1_c = -(powerlaw_c1_int(x1, q1, m)-powerlaw_c1_int(-q1, q1, m))
      a1_w = 0
      #if they are equal then it must mean that they are both baseline
    }else{
      a1_w = 0
      a1_c = 0
    }
    #second thermode warm: integrate from x2 to q
    if(Tw2>Tc2){
      a2_w = powerlaw_w2_int(q2,q2,n)-powerlaw_w2_int(x2,q2,n)
      a2_c = 0
      #second thermode cold integrate from x2 to q and -
    }else if(Tw2<Tc2){
      a2_c = -(powerlaw_c2_int(q2,q2,m)-powerlaw_c2_int(x2,q2,m))
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
  
  
  Atotw = (powerlaw_w1_int(0,q1,n)-powerlaw_w1_int(-q1,q1,n))+(powerlaw_w2_int(q2,q2,n)-powerlaw_w2_int(0,q2,n))
  Atotc = -((powerlaw_c1_int(0,q1,m)-powerlaw_c1_int(-q1,q1,m))   +   (powerlaw_c2_int(q2,q2,m)-powerlaw_c2_int(0,q2,m)))
  
  data = data.frame(Tw1 = Tw1, Tw2 = Tw2, Tc1 = Tc1, Tc2 = Tc2, q1 = q1,q2 = q2, n = n, m = m)
  
  data = pmap(data, get_points) %>% bind_rows() %>% bind_cols(data, .)
  data = data %>% select(-c("y1","y2")) %>% pmap(get_end_areas) %>% bind_rows() %>% bind_cols(data, .)
  data = data %>% select(-c("a1_w","a2_w","a1_c","a2_c")) %>% pmap(get_middle_areas) %>% bind_rows() %>% bind_cols(data, .)
  
  if (dist == TRUE){
    data = data %>% mutate(Atotw = Atotw,
                           Atotc = Atotc,
                           uw = ((a1_w+a2_w+a3_w)/Atotw)+0.001,
                           uc = ((a1_c+a2_c+a3_c)/Atotc)+0.001,
                           up = uc*uw+0.001)
    
    if(autocor == TRUE){
      data$uw = ifelse((data$uw+bw*lag(data$uw, 1, default = 0)) <= 0 , 0.01, ifelse((data$uw+bw*lag(data$uw, 1, default = 0)+0.01) >= 1, 0.99, data$uw+bw*lag(data$uw, 1, default = 0)))
      
      data$uc = ifelse((data$uc+bc*lag(data$uc, 1, default = 0)) <= 0, 0.01, ifelse((data$uc+bc*lag(data$uc, 1, default = 0)+0.01) >= 1, 0.99, data$uc+bc*lag(data$uc, 1, default = 0)))
      
      data$up = ifelse((data$up+bp*lag(data$up, 1, default = 0)) <= 0, 0.01, ifelse((data$up+bp*lag(data$up, 1, default = 0)+0.01) >= 1, 0.99, data$up+bp*lag(data$up, 1, default = 0)))
    }
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uw") %>% rename(mean = uw) %>% pmap(betaprop) %>% bind_rows() %>% rename(warm = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uc") %>% rename(mean = uc) %>% pmap(betaprop) %>% bind_rows() %>% rename(cold = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","up") %>% rename(mean = up) %>% pmap(betaprop) %>% bind_rows() %>% rename(pain = value) %>% bind_cols(data, .)
    
  
  return(data)
  
  }
  
  
  if (dist == FALSE){
    data = data %>% mutate(Atotw = Atotw,
                           Atotc = Atotc,
                           uw = ((a1_w+a2_w+a3_w)/Atotw),
                           uc = ((a1_c+a2_c+a3_c)/Atotc),
                           up = uc*uw)
    
    if(autocor == TRUE){
      data$uw = data$uw+bw*lag(data$uw, 1, default = 0)
      
      data$uc = data$uc+bc*lag(data$uc, 1, default = 0)
      
      data$up = data$up+bp*lag(data$up, 1, default = 0)
    }
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uw") %>% rename(mean = uw) %>% pmap(gaus) %>% bind_rows() %>% rename(warm = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uc") %>% rename(mean = uc) %>% pmap(gaus) %>% bind_rows() %>% rename(cold = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","up") %>% rename(mean = up) %>% pmap(gaus) %>% bind_rows() %>% rename(pain = value) %>% bind_cols(data, .)
    
    
    return(data)
    
  }
  
  
  
}




poweragent_generalized_v2 = function(m,n,q1,q2,kappa, autocor = FALSE,bw = 0,bc = 0,bp = 0, dist = TRUE){
  
  powerlaw_w1 = function(x, q, n){
    return((x+q)^n)
  }
  powerlaw_w1_int= function(x, q, n){
    return((((x+q)^(n+1))/(n+1)))
  }
  
  powerlaw_w2= function(x,  q,  n){
    return((-x+q)^n)
  }
  
  powerlaw_w2_int= function(x, q, n){
    return(-(((-x+q)^(n+1))/(n+1)))
  }
  
  
  powerlaw_c1 = function(x,  q, m){
    return(-(x+q)^m)
  }
  
  powerlaw_c1_int= function( x, q, m){
    return(-(((x+q)^(m+1))/(m+1)))
  }
  
  powerlaw_c2 = function(x, q, m){
    return(-(-x+q)^m)
  }
  
  powerlaw_c2_int= function(x, q, m){
    return((((-x+q)^(m+1))/(m+1)))
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
  
  
  get_points = function(Tw1 , Tw2 , Tc1 , Tc2 , q1, q2 , n , m){
    #both warm
    if(Tw1>Tc1){
      x1 = Tw1-q1
      x2 = q2-Tw2
      y1 = powerlaw_w1(x1,q1,n)
      y2 = powerlaw_w2(x2,q2,n)
      #both cold
    }else if(Tc1>0 && Tc2 >0){
      x1 = Tc1-q1
      x2 = q2-Tc2
      y1 = powerlaw_c1(x1,q1,m)
      y2 = powerlaw_c2(x2,q2,m)
      #1 cold 2 warm
    }else if(Tc1>0 && Tw2 >0){
      x1 = Tc1-q1
      x2 = q2-Tw2
      y1 = powerlaw_c1(x1,q1,m)
      y2 = powerlaw_w2(x2,q2,n)
      #second neutral and first cold
    }else if(Tw2 == 0 && Tc2 == 0 && Tc1>0){
      x1 = Tc1-q1
      x2 = q2
      y1 = powerlaw_c1(x1,q1,m)
      y2 = 0
      #first neutral and therefore other warm (given the above takes care if its neutral)
    }else{
      x1 = -q1
      x2 = q2-Tw2
      y1 = 0
      y2 = powerlaw_w2(x2,q2,n)
    }
    
    points = c("x1" = x1,"x2" = x2,"y1" = y1, "y2" = y2)
    
    return(data.frame(t(points)))
  }
  
  get_end_areas = function(Tw1,  Tw2, Tc1, Tc2, q1, q2, n, m, x1,x2){
    
    #being warm for the first thermode: we integrate from -q to x1
    if(Tw1>Tc1){
      a1_w = powerlaw_w1_int(x1, q1, n)-powerlaw_w1_int(-q1, q1, n)
      a1_c = 0
      #cold for the first thermode: we integrate from -q to x1 and
    }else if(Tw1<Tc1){
      a1_c = -(powerlaw_c1_int(x1, q1, m)-powerlaw_c1_int(-q1, q1, m))
      a1_w = 0
      #if they are equal then it must mean that they are both baseline
    }else{
      a1_w = 0
      a1_c = 0
    }
    #second thermode warm: integrate from x2 to q
    if(Tw2>Tc2){
      a2_w = powerlaw_w2_int(q2,q2,n)-powerlaw_w2_int(x2,q2,n)
      a2_c = 0
      #second thermode cold integrate from x2 to q and -
    }else if(Tw2<Tc2){
      a2_c = -(powerlaw_c2_int(q2,q2,m)-powerlaw_c2_int(x2,q2,m))
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
  
  
  Atotw = (powerlaw_w1_int(0,q1,n)-powerlaw_w1_int(-q1,q1,n))+(powerlaw_w2_int(q2,q2,n)-powerlaw_w2_int(0,q2,n))
  Atotc = -((powerlaw_c1_int(0,q1,m)-powerlaw_c1_int(-q1,q1,m))   +   (powerlaw_c2_int(q2,q2,m)-powerlaw_c2_int(0,q2,m)))
  
  data = data.frame(Tw1 = Tw1, Tw2 = Tw2, Tc1 = Tc1, Tc2 = Tc2, q1 = q1,q2 = q2, n = n, m = m)
  
  data = pmap(data, get_points) %>% bind_rows() %>% bind_cols(data, .)
  data = data %>% select(-c("y1","y2")) %>% pmap(get_end_areas) %>% bind_rows() %>% bind_cols(data, .)
  data = data %>% select(-c("a1_w","a2_w","a1_c","a2_c")) %>% pmap(get_middle_areas) %>% bind_rows() %>% bind_cols(data, .)
  
  if (dist == TRUE){
    data = data %>% mutate(Atotw = Atotw,
                           Atotc = Atotc,
                           uw = ((a1_w+a2_w+a3_w)/Atotw)+0.001,
                           uc = ((a1_c+a2_c+a3_c)/Atotc)+0.001,
                           up = uc*uw+0.001)
    
    if(autocor == TRUE){
      data$uw = ifelse((lag(data$uw, 1, default = 0)+(1/(1+bw))*(data$uw-lag(data$uw, 1, default = 0))) <= 0 , 0.01, lag(data$uw, 1, default = 0)+(1/(1+bw))*(data$uw-lag(data$uw, 1, default = 0)))
      
      data$uc = ifelse((lag(data$uc, 1, default = 0)+(1/(1+bc))*(data$uc-lag(data$uc, 1, default = 0))) <= 0, 0.01, lag(data$uc, 1, default = 0)+(1/(1+bc))*(data$uc-lag(data$uc, 1, default = 0)))
      
      data$up = ifelse((lag(data$up, 1, default = 0)+(1/(1+bp))*(data$up-lag(data$up, 1, default = 0))) <= 0, 0.01, lag(data$up, 1, default = 0)+(1/(1+bp))*(data$up-lag(data$up, 1, default = 0)))
    }
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uw") %>% rename(mean = uw) %>% pmap(betaprop) %>% bind_rows() %>% rename(warm = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uc") %>% rename(mean = uc) %>% pmap(betaprop) %>% bind_rows() %>% rename(cold = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","up") %>% rename(mean = up) %>% pmap(betaprop) %>% bind_rows() %>% rename(pain = value) %>% bind_cols(data, .)
    
    
    return(data)
    
  }
  
  
  if (dist == FALSE){
    data = data %>% mutate(Atotw = Atotw,
                           Atotc = Atotc,
                           uw = ((a1_w+a2_w+a3_w)/Atotw),
                           uc = ((a1_c+a2_c+a3_c)/Atotc),
                           up = uc*uw)
    
    if(autocor == TRUE){
      data$uw = lag(data$uw, 1, default = 0)+(1/(1+bw))*(data$uw-lag(data$uw, 1, default = 0))
      
      data$uc = lag(data$uc, 1, default = 0)+(1/(1+bc))*(data$uc-lag(data$uc, 1, default = 0))
      
      data$up = lag(data$up, 1, default = 0)+(1/(1+bp))*(data$up-lag(data$up, 1, default = 0))
    }
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uw") %>% rename(mean = uw) %>% pmap(gaus) %>% bind_rows() %>% rename(warm = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uc") %>% rename(mean = uc) %>% pmap(gaus) %>% bind_rows() %>% rename(cold = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","up") %>% rename(mean = up) %>% pmap(gaus) %>% bind_rows() %>% rename(pain = value) %>% bind_cols(data, .)
    
    
    return(data)
    
  }
  
  
  
}



poweragent_generalized_v2_simp = function(m,n,q1,q2,kappa, autocor = FALSE,b, dist = TRUE){
  
  powerlaw_w1 = function(x, q, n){
    return((x+q)^n)
  }
  powerlaw_w1_int= function(x, q, n){
    return((((x+q)^(n+1))/(n+1)))
  }
  
  powerlaw_w2= function(x,  q,  n){
    return((-x+q)^n)
  }
  
  powerlaw_w2_int= function(x, q, n){
    return(-(((-x+q)^(n+1))/(n+1)))
  }
  
  
  powerlaw_c1 = function(x,  q, m){
    return(-(x+q)^m)
  }
  
  powerlaw_c1_int= function( x, q, m){
    return(-(((x+q)^(m+1))/(m+1)))
  }
  
  powerlaw_c2 = function(x, q, m){
    return(-(-x+q)^m)
  }
  
  powerlaw_c2_int= function(x, q, m){
    return((((-x+q)^(m+1))/(m+1)))
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
  
  
  get_points = function(Tw1 , Tw2 , Tc1 , Tc2 , q1, q2 , n , m){
    #both warm
    if(Tw1>Tc1){
      x1 = Tw1-q1
      x2 = q2-Tw2
      y1 = powerlaw_w1(x1,q1,n)
      y2 = powerlaw_w2(x2,q2,n)
      #both cold
    }else if(Tc1>0 && Tc2 >0){
      x1 = Tc1-q1
      x2 = q2-Tc2
      y1 = powerlaw_c1(x1,q1,m)
      y2 = powerlaw_c2(x2,q2,m)
      #1 cold 2 warm
    }else if(Tc1>0 && Tw2 >0){
      x1 = Tc1-q1
      x2 = q2-Tw2
      y1 = powerlaw_c1(x1,q1,m)
      y2 = powerlaw_w2(x2,q2,n)
      #second neutral and first cold
    }else if(Tw2 == 0 && Tc2 == 0 && Tc1>0){
      x1 = Tc1-q1
      x2 = q2
      y1 = powerlaw_c1(x1,q1,m)
      y2 = 0
      #first neutral and therefore other warm (given the above takes care if its neutral)
    }else{
      x1 = -q1
      x2 = q2-Tw2
      y1 = 0
      y2 = powerlaw_w2(x2,q2,n)
    }
    
    points = c("x1" = x1,"x2" = x2,"y1" = y1, "y2" = y2)
    
    return(data.frame(t(points)))
  }
  
  get_end_areas = function(Tw1,  Tw2, Tc1, Tc2, q1, q2, n, m, x1,x2){
    
    #being warm for the first thermode: we integrate from -q to x1
    if(Tw1>Tc1){
      a1_w = powerlaw_w1_int(x1, q1, n)-powerlaw_w1_int(-q1, q1, n)
      a1_c = 0
      #cold for the first thermode: we integrate from -q to x1 and
    }else if(Tw1<Tc1){
      a1_c = -(powerlaw_c1_int(x1, q1, m)-powerlaw_c1_int(-q1, q1, m))
      a1_w = 0
      #if they are equal then it must mean that they are both baseline
    }else{
      a1_w = 0
      a1_c = 0
    }
    #second thermode warm: integrate from x2 to q
    if(Tw2>Tc2){
      a2_w = powerlaw_w2_int(q2,q2,n)-powerlaw_w2_int(x2,q2,n)
      a2_c = 0
      #second thermode cold integrate from x2 to q and -
    }else if(Tw2<Tc2){
      a2_c = -(powerlaw_c2_int(q2,q2,m)-powerlaw_c2_int(x2,q2,m))
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
  
  
  Atotw = (powerlaw_w1_int(0,q1,n)-powerlaw_w1_int(-q1,q1,n))+(powerlaw_w2_int(q2,q2,n)-powerlaw_w2_int(0,q2,n))
  Atotc = -((powerlaw_c1_int(0,q1,m)-powerlaw_c1_int(-q1,q1,m))   +   (powerlaw_c2_int(q2,q2,m)-powerlaw_c2_int(0,q2,m)))
  
  data = data.frame(Tw1 = Tw1, Tw2 = Tw2, Tc1 = Tc1, Tc2 = Tc2, q1 = q1,q2 = q2, n = n, m = m)
  
  data = pmap(data, get_points) %>% bind_rows() %>% bind_cols(data, .)
  data = data %>% select(-c("y1","y2")) %>% pmap(get_end_areas) %>% bind_rows() %>% bind_cols(data, .)
  data = data %>% select(-c("a1_w","a2_w","a1_c","a2_c")) %>% pmap(get_middle_areas) %>% bind_rows() %>% bind_cols(data, .)
  
  if (dist == TRUE){
    data = data %>% mutate(Atotw = Atotw,
                           Atotc = Atotc,
                           uw = ((a1_w+a2_w+a3_w)/Atotw)+0.001,
                           uc = ((a1_c+a2_c+a3_c)/Atotc)+0.001,
                           up = uc*uw+0.001)
    
    if(autocor == TRUE){
      data$uw = ifelse((lag(data$uw, 1, default = 0)+(1/(1+b))*(data$uw-lag(data$uw, 1, default = 0))) <= 0 , 0.001, lag(data$uw, 1, default = 0)+(1/(1+b))*(data$uw-lag(data$uw, 1, default = 0)))
      
      data$uc = ifelse((lag(data$uc, 1, default = 0)+(1/(1+b))*(data$uc-lag(data$uc, 1, default = 0))) <= 0, 0.001, lag(data$uc, 1, default = 0)+(1/(1+b))*(data$uc-lag(data$uc, 1, default = 0)))
      
      }
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uw") %>% rename(mean = uw) %>% pmap(betaprop) %>% bind_rows() %>% rename(warm = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uc") %>% rename(mean = uc) %>% pmap(betaprop) %>% bind_rows() %>% rename(cold = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","up") %>% rename(mean = up) %>% pmap(betaprop) %>% bind_rows() %>% rename(pain = value) %>% bind_cols(data, .)
    
    
    return(data)
    
  }
  
  
  if (dist == FALSE){
    data = data %>% mutate(Atotw = Atotw,
                           Atotc = Atotc,
                           uw = ((a1_w+a2_w+a3_w)/Atotw),
                           uc = ((a1_c+a2_c+a3_c)/Atotc),
                           up = uc*uw)
    
    if(autocor == TRUE){
      data$uw = lag(data$uw, 1, default = 0)+(1/(1+b))*(data$uw-lag(data$uw, 1, default = 0))
      
      data$uc = lag(data$uc, 1, default = 0)+(1/(1+b))*(data$uc-lag(data$uc, 1, default = 0))
      
    }
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uw") %>% rename(mean = uw) %>% pmap(gaus) %>% bind_rows() %>% rename(warm = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uc") %>% rename(mean = uc) %>% pmap(gaus) %>% bind_rows() %>% rename(cold = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","up") %>% rename(mean = up) %>% pmap(gaus) %>% bind_rows() %>% rename(pain = value) %>% bind_cols(data, .)
    
    data$T1 = temp$T2
    data$T2 = temp$T1
    
    return(data)
    
  }
  
  
  
}


brms_agent = function(i,b1,b2,b3){
  
  
  
  
  
}


poweragent_generalized_v2_doubled = function(m,n,q1,q2,kappa, autocor = FALSE,bw = 0,bc = 0,bp = 0, dist = TRUE){
  
  powerlaw_w1 = function(x, q, n){
    return((x+q)^n)
  }
  powerlaw_w1_int= function(x, q, n){
    return((((x+q)^(n+1))/(n+1)))
  }
  
  powerlaw_w2= function(x,  q,  n){
    return((-x+q)^n)
  }
  
  powerlaw_w2_int= function(x, q, n){
    return(-(((-x+q)^(n+1))/(n+1)))
  }
  
  
  powerlaw_c1 = function(x,  q, m){
    return(-(x+q)^m)
  }
  
  powerlaw_c1_int= function( x, q, m){
    return(-(((x+q)^(m+1))/(m+1)))
  }
  
  powerlaw_c2 = function(x, q, m){
    return(-(-x+q)^m)
  }
  
  powerlaw_c2_int= function(x, q, m){
    return((((-x+q)^(m+1))/(m+1)))
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
  
  
  get_points = function(Tw1 , Tw2 , Tc1 , Tc2 , q1, q2 , n , m){
    #both warm
    if(Tw1>Tc1){
      x1 = Tw1-q1
      x2 = q2-Tw2
      y1 = powerlaw_w1(x1,q1,n)
      y2 = powerlaw_w2(x2,q2,n)
      #both cold
    }else if(Tc1>0 && Tc2 >0){
      x1 = Tc1-q1
      x2 = q2-Tc2
      y1 = powerlaw_c1(x1,q1,m)
      y2 = powerlaw_c2(x2,q2,m)
      #1 cold 2 warm
    }else if(Tc1>0 && Tw2 >0){
      x1 = Tc1-q1
      x2 = q2-Tw2
      y1 = powerlaw_c1(x1,q1,m)
      y2 = powerlaw_w2(x2,q2,n)
      #second neutral and first cold
    }else if(Tw2 == 0 && Tc2 == 0 && Tc1>0){
      x1 = Tc1-q1
      x2 = q2
      y1 = powerlaw_c1(x1,q1,m)
      y2 = 0
      #first neutral and therefore other warm (given the above takes care if its neutral)
    }else{
      x1 = -q1
      x2 = q2-Tw2
      y1 = 0
      y2 = powerlaw_w2(x2,q2,n)
    }
    
    points = c("x1" = x1,"x2" = x2,"y1" = y1, "y2" = y2)
    
    return(data.frame(t(points)))
  }
  
  get_end_areas = function(Tw1,  Tw2, Tc1, Tc2, q1, q2, n, m, x1,x2){
    
    #being warm for the first thermode: we integrate from -q to x1
    if(Tw1>Tc1){
      a1_w = powerlaw_w1_int(x1, q1, n)-powerlaw_w1_int(-q1, q1, n)
      a1_c = 0
      #cold for the first thermode: we integrate from -q to x1 and
    }else if(Tw1<Tc1){
      a1_c = -(powerlaw_c1_int(x1, q1, m)-powerlaw_c1_int(-q1, q1, m))
      a1_w = 0
      #if they are equal then it must mean that they are both baseline
    }else{
      a1_w = 0
      a1_c = 0
    }
    #second thermode warm: integrate from x2 to q
    if(Tw2>Tc2){
      a2_w = powerlaw_w2_int(q2,q2,n)-powerlaw_w2_int(x2,q2,n)
      a2_c = 0
      #second thermode cold integrate from x2 to q and -
    }else if(Tw2<Tc2){
      a2_c = -(powerlaw_c2_int(q2,q2,m)-powerlaw_c2_int(x2,q2,m))
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
  
  temp = get_temps(ntrials = 6)
  
  Tw1 = ifelse(temp$T1 <= 30, 0 ,ifelse(temp$T1 > 30, temp$T1-30, NA))
  Tc1 = ifelse(temp$T1 >= 30, 0 ,ifelse(temp$T1 < 30, 30-temp$T1, NA))
  Tw2 = ifelse(temp$T2 <= 30, 0 ,ifelse(temp$T2 > 30, temp$T2-30, NA))
  Tc2 = ifelse(temp$T2 >= 30, 0 ,ifelse(temp$T2 < 30, 30-temp$T2, NA))
  
  
  Atotw = (powerlaw_w1_int(0,q1,n)-powerlaw_w1_int(-q1,q1,n))+(powerlaw_w2_int(q2,q2,n)-powerlaw_w2_int(0,q2,n))
  Atotc = -((powerlaw_c1_int(0,q1,m)-powerlaw_c1_int(-q1,q1,m))   +   (powerlaw_c2_int(q2,q2,m)-powerlaw_c2_int(0,q2,m)))
  
  data = data.frame(Tw1 = Tw1, Tw2 = Tw2, Tc1 = Tc1, Tc2 = Tc2, q1 = q1,q2 = q2, n = n, m = m)
  
  data = pmap(data, get_points) %>% bind_rows() %>% bind_cols(data, .)
  data = data %>% select(-c("y1","y2")) %>% pmap(get_end_areas) %>% bind_rows() %>% bind_cols(data, .)
  data = data %>% select(-c("a1_w","a2_w","a1_c","a2_c")) %>% pmap(get_middle_areas) %>% bind_rows() %>% bind_cols(data, .)
  
  if (dist == TRUE){
    data = data %>% mutate(Atotw = Atotw,
                           Atotc = Atotc,
                           uw = ((a1_w+a2_w+a3_w)/Atotw)+0.001,
                           uc = ((a1_c+a2_c+a3_c)/Atotc)+0.001,
                           up = uc*uw+0.001)
    
    if(autocor == TRUE){
      data$uw = ifelse((lag(data$uw, 1, default = 0)+(1/(1+bw))*(data$uw-lag(data$uw, 1, default = 0))) <= 0 , 0.01, lag(data$uw, 1, default = 0)+(1/(1+bw))*(data$uw-lag(data$uw, 1, default = 0)))
      
      data$uc = ifelse((lag(data$uc, 1, default = 0)+(1/(1+bc))*(data$uc-lag(data$uc, 1, default = 0))) <= 0, 0.01, lag(data$uc, 1, default = 0)+(1/(1+bc))*(data$uc-lag(data$uc, 1, default = 0)))
      
      data$up = ifelse((lag(data$up, 1, default = 0)+(1/(1+bp))*(data$up-lag(data$up, 1, default = 0))) <= 0, 0.01, lag(data$up, 1, default = 0)+(1/(1+bp))*(data$up-lag(data$up, 1, default = 0)))
    }
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uw") %>% rename(mean = uw) %>% pmap(betaprop) %>% bind_rows() %>% rename(warm = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uc") %>% rename(mean = uc) %>% pmap(betaprop) %>% bind_rows() %>% rename(cold = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","up") %>% rename(mean = up) %>% pmap(betaprop) %>% bind_rows() %>% rename(pain = value) %>% bind_cols(data, .)
    
    
    return(data)
    
  }
  
  
  if (dist == FALSE){
    data = data %>% mutate(Atotw = Atotw,
                           Atotc = Atotc,
                           uw = ((a1_w+a2_w+a3_w)/Atotw),
                           uc = ((a1_c+a2_c+a3_c)/Atotc),
                           up = uc*uw)
    
    if(autocor == TRUE){
      data$uw = lag(data$up, 1, default = 0)+(1/(1+bp))*(data$up-lag(data$up, 1, default = 0))
      
      data$uc = lag(data$up, 1, default = 0)+(1/(1+bp))*(data$up-lag(data$up, 1, default = 0))
      
      data$up = lag(data$up, 1, default = 0)+(1/(1+bp))*(data$up-lag(data$up, 1, default = 0))
    }
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uw") %>% rename(mean = uw) %>% pmap(gaus) %>% bind_rows() %>% rename(warm = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","uc") %>% rename(mean = uc) %>% pmap(gaus) %>% bind_rows() %>% rename(cold = value) %>% bind_cols(data, .)
    
    data = data %>% mutate(kappa = kappa) %>% select("kappa","up") %>% rename(mean = up) %>% pmap(gaus) %>% bind_rows() %>% rename(pain = value) %>% bind_cols(data, .)
    
    
    return(data)
    
  }
  
  
  
}
