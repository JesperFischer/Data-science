functions {
  
  real powerlaw_w1(real x,real q,real n,real kw){
    return(kw*(x+q)^n);
  }
  real powerlaw_w1_int(real x,real q,real n,real kw){
    return(kw * (((x+q)^(n+1))/(n+1)) );
  }
  
  real powerlaw_w2 (real x, real q, real n, real kw){
    return(kw * (-x+q)^n);
  }
  
  real powerlaw_w2_int(real x,real q,real n,real kw){
    return(-kw* (((-x+q)^(n+1))/(n+1)));
  }
  
  real powerlaw_c1 (real x, real q,real m,real kc){
    return(-kc*(x+q)^m);
  }
  
  real powerlaw_c1_int(real x,real q,real m,real kc){
    return(-kc*(((x+q)^(m+1))/(m+1)));
  }
  
  real powerlaw_c2 (real x,real q,real m,real kc){
    return(-kc * (-x+q)^m);
  }
  
  real powerlaw_c2_int(real x,real q,real m,real kc){
    return(kc * (((-x+q)^(m+1))/(m+1)));
  }
  
  
  real linear (real xs,real x1,real x2,real y1,real y2){
    real a = ((y2-y1)/(x2-x1));
    return(a*xs+y1-(a*x1));
    
  }
  
  real linear_int (real xs,real x1,real x2,real y1,real y2){
    
    real a = (y2-y1)/(x2-x1);
    real b = y1-a*x1;
    return(((a*xs^2)/2.0)+(b*xs));
  }
 
  vector get_points(real Tw1, real Tw2, real Tc1, real Tc2,real q1, real q2, real n, real m, real kw, real kc){
    real  x1;
    real  x2;
    real  y1;
    real  y2;
    vector[4] points;
    
    //both warm
    if(Tw1>Tc1){
      x1 = Tw1-q1;
      x2 = q2-Tw2;
      y1 = powerlaw_w1(x1,q1,n,kw);
      y2 = powerlaw_w2(x2,q2,n,kw);
    //both cold
    }else if(Tc1>0 && Tc2 >0){
      x1 = Tc1-q1;
      x2 = q2-Tc2;
      y1 = powerlaw_c1(x1,q1,m,kc);
      y2 = powerlaw_c2(x2,q2,m,kc);
      //1 cold; 2 warm
    }else if(Tc1>0 && Tw2 >0){
      x1 = Tc1-q1;
      x2 = q2-Tw2;
      y1 = powerlaw_c1(x1,q1,m,kc);
      y2 = powerlaw_w2(x2,q2,n,kw);
      //second neutral and first cold
    }else if(Tw2 == 0 && Tc2 == 0 && Tc1>0){
      x1 = Tc1-q1;
      x2 = q2;
      y1 = powerlaw_c1(x1,q1,m,kc);
      y2 = 0;
      //first neutral and therefore other warm (given the above takes care if its neutral)
    }else{
      x1 = -q1;
      x2 = q2-Tw2;
      y1 = 0;
      y2 = powerlaw_w2(x2,q2,n,kw);
    }
    
    points[1] = x1;
    points[2] = x2;
    points[3] = y1;
    points[4] = y2;
    
    return(points);
  }
 
  vector get_end_areas(real Tw1, real Tw2, real Tc1, real Tc2,real q1, real q2, real n, real m, vector points, real kw, real kc){
    real a1_w;
    real a2_w;
    real a1_c;
    real a2_c;
    real x1;
    real x2;
    
    
    vector[4] areas;
    x1 = points[1];
    x2 = points[2];
    
      
    //being warm for the first thermode: we integrate from -q to x1
    if(Tw1>Tc1){
      a1_w = powerlaw_w1_int(x1, q1, n, kw)-powerlaw_w1_int(-q1, q1, n, kw);
      a1_c = 0;
    //cold for the first thermode: we integrate from -q to x1 and
    }else if(Tw1<Tc1){
      a1_c = -(powerlaw_c1_int(x1, q1, m, kc)-powerlaw_c1_int(-q1, q1, m, kc));
      a1_w = 0;
    //if they are equal then it must mean that they are both baseline
    }else{
      a1_w = 0;
      a1_c = 0;
    }
    //second thermode warm: integrate from x2 to q
    if(Tw2>Tc2){
      a2_w = powerlaw_w2_int(q2,q2,n, kw)-powerlaw_w2_int(x2,q2,n, kw);
      a2_c = 0;
    //second thermode cold integrate from x2 to q and -
    }else if(Tw2<Tc2){
      a2_c = -(powerlaw_c2_int(q2,q2,m, kc)-powerlaw_c2_int(x2,q2,m, kc));
      a2_w = 0;
    //if they are equal then it must mean that they are both baseline
    }else{
      a2_w = 0;
      a2_c = 0;
    }
     
    areas[1] = a1_w;
    areas[2] = a2_w;
    areas[3] = a1_c;
    areas[4] = a2_c;
    
    return(areas);
    
  }
  
  vector get_middle_areas(real Tw1, real Tw2, real Tc1, real Tc2,real q1,real q2, real n, real m, vector points){
    
    real a3_w;
    real a3_c;
    real a;
    real xint;
    
    real x1;
    real x2;
    real y1;
    real y2;
    
    vector[2] areas3;
    
    
    x1 = points[1];
    x2 = points[2];
    y1 = points[3];
    y2 = points[4];
    


    if(Tw2== 0 && Tw1 == 0 && Tc2 == 0 && Tc1 == 0){
      a3_w = 0;
      a3_c = 0;
      //if both are warm
    }else if(Tw2>0 && Tw1 > 0){
      //if they are not equal i.e. a vertical line:
      if(Tw2 != Tw1){
        //integrate linear from x1 to x2
        a3_w = linear_int(x2,x1,x2,y1,y2)-linear_int(x1,x1,x2,y1,y2);
        a3_c = 0;
        }else{
          //if they are equal then the slope is = 0 and its just a square
          a3_w = (-x1+x2)*(y1);
          a3_c = 0;
          }
    // if both are cold
    }else if(Tc2>0 && Tc1 > 0){
      if(Tc2 != Tc1){
        a3_c = -(linear_int(x2,x1,x2,y1,y2)-linear_int(x1,x1,x2,y1,y2));
        a3_w = 0;
      }else{
        a3_c = -(-x1+x2)*y1;
        a3_w = 0;
      }
    // if 1 is cold and 2 is warm
    }else if(Tc1>=0 && Tw2>=0){
    
      a = (y2-y1)/(x2-x1);
      xint = ((a*x1)-y1)/a;
      
      a3_w = linear_int(x2,x1,x2,y1,y2)-linear_int(xint,x1,x2,y1,y2);
      a3_c = -(linear_int(xint,x1,x2,y1,y2)-linear_int(x1,x1,x2,y1,y2));
    }else{
      a3_w = 0;
      a3_c = 0;
    }

    areas3[1] = a3_w;
    areas3[2] = a3_c;
    
    return(areas3);
    
  }

}


data {
  int<lower=0, upper = 1> dist;
  int<lower=0> N;
  vector[N] Tw1;
  vector[N] Tw2;
  vector[N] Tc1;
  vector[N] Tc2;
  vector[N] w;
  vector[N] c;
  vector[N] p;
  
  
  
}

parameters {
  real<lower = 0, upper =100> q1;
  real<lower=0> kappa;
  
}



transformed parameters{
  vector[N] aw1;
  vector[N] aw2;
  vector[N] aw3;
  vector[N] ac1;
  vector[N] ac2;
  vector[N] ac3;
  vector[N] wa = rep_vector(0, N);
  vector[N] ca = rep_vector(0, N);
  vector[N] pa = rep_vector(0, N);
  vector[N] wa1;
  vector[N] ca1;
  vector[N] pa1;
  vector<lower =0> [N] uw;
  vector<lower =0> [N] uc;
  vector<lower =0> [N] up;
  real Atotw;
  real Atotc;
  
  real q2 = q1;
  real n = 1;
  real m = 1;
  real kw = 1;
  real kc = 1;
  real bw = 0;
  real bc = 0;
  real bp = 0;
  
  
  
  
  Atotw = (powerlaw_w1_int(0,q1,n,kw)-powerlaw_w1_int(-q1,q1,n,kw))+(powerlaw_w2_int(q2,q2,n,kw)-powerlaw_w2_int(0,q2,n,kw));
  Atotc = -((powerlaw_c1_int(0,q1,m,kc)-powerlaw_c1_int(-q1,q1,m,kc)) + ((powerlaw_c2_int(q2,q2,m,kc)-powerlaw_c2_int(0,q2,m,kc))));

  for (i in 1:N){
    // [aw1,aw2,ac1,ac2]
    aw1[i] = get_end_areas(Tw1[i],Tw2[i],Tc1[i],Tc2[i],q1,q2,n,m,get_points(Tw1[i],Tw2[i],Tc1[i],Tc2[i],q1,q2,n,m,kw,kc),kw,kc)[1];
    
    aw2[i] = get_end_areas(Tw1[i],Tw2[i],Tc1[i],Tc2[i],q1,q2,n,m,get_points(Tw1[i],Tw2[i],Tc1[i],Tc2[i],q1,q2,n,m,kw,kc),kw,kc)[2];
    
    ac1[i] = get_end_areas(Tw1[i],Tw2[i],Tc1[i],Tc2[i],q1,q2,n,m,get_points(Tw1[i],Tw2[i],Tc1[i],Tc2[i],q1,q2,n,m,kw,kc),kw,kc)[3];
    
    ac2[i] = get_end_areas(Tw1[i],Tw2[i],Tc1[i],Tc2[i],q1,q2,n,m,get_points(Tw1[i],Tw2[i],Tc1[i],Tc2[i],q1,q2,n,m,kw,kc),kw,kc)[4];
    // [aw3, ac3]
    
    aw3[i] = get_middle_areas(Tw1[i],Tw2[i],Tc1[i],Tc2[i],q1,q2,n,m,get_points(Tw1[i],Tw2[i],Tc1[i],Tc2[i],q1,q2,n,m,kw,kc))[1];
    
    ac3[i] = get_middle_areas(Tw1[i],Tw2[i],Tc1[i],Tc2[i],q1,q2,n,m,get_points(Tw1[i],Tw2[i],Tc1[i],Tc2[i],q1,q2,n,m,kw,kc))[2];
    
    wa1[i] = (aw1[i]+aw2[i]+aw3[i])/Atotw;
    ca1[i] = (ac1[i]+ac2[i]+ac3[i])/Atotc;
    pa1[i] = wa1[i]*ca1[i];
    

    if(i == 1){
      //prior + weighting (likelihood - prior)
      uw[i] = wa[i]+(1/(1+bw))*(wa1[i]-wa[i]);
      uc[i] = ca[i]+(1/(1+bc))*(ca1[i]-ca[i]);
      up[i] = pa[i]+(1/(1+bp))*(pa1[i]-pa[i]);
      
      wa[i] = wa1[i];
      ca[i] = ca1[i];
      pa[i] = pa1[i];
      
    }else{
      if(wa[i-1]+(1/(1+bw))*(wa1[i]-wa[i-1]) < 0){uw[i] = 0;}
      else{
        uw[i] = wa[i-1]+(1/(1+bw))*(wa1[i]-wa[i-1]);
        }
      
      
      if(ca[i-1]+(1/(1+bc))*(ca1[i]-ca[i-1]) < 0){uc[i] = 0;}
      else{
        uc[i] = ca[i-1]+(1/(1+bc))*(ca1[i]-ca[i-1]);
        }
      
      
      if(pa[i-1]+(1/(1+bp))*(pa1[i]-pa[i-1]) < 0){up[i] = 0;}
      else{
        up[i] = pa[i-1]+(1/(1+bp))*(pa1[i]-pa[i-1]);
        }
      
      
      wa[i] = wa1[i];
      ca[i] = ca1[i];
      pa[i] = pa1[i];
    }
    
  
  }
  

}


model{
  
  target += lognormal_lpdf(q1|3,0.5);
  
  target += lognormal_lpdf(kappa|3.5,0.5);
  
  
  for(i in 1:N){
    if(dist){
      target += beta_proportion_lpdf(w[i] | uw[i]+0.001, kappa);
      target += beta_proportion_lpdf(c[i] | uc[i]+0.001, kappa);
      target += beta_proportion_lpdf(p[i] | up[i]+0.001, kappa);
      }else{
        target += normal_lpdf(w[i] | uw[i], 1/kappa^2);
        target += normal_lpdf(c[i] | uc[i], 1/kappa^2);
        target += normal_lpdf(p[i] | up[i], 1/kappa^2);
      }
    }
}

generated quantities{
 vector[N] log_lik;

 real prior_q1 = lognormal_rng(3,0.5);


 
    for(i in 1:N){
      if(dist){
        log_lik[i] = beta_proportion_lpdf(w[i] | uw[i]+0.001, kappa)+beta_proportion_lpdf(c[i] | uc[i]+0.001, kappa)+beta_proportion_lpdf(p[i] | up[i]+0.01, kappa);
        }else{
         log_lik[i] = normal_lpdf(w[i] | uw[i], kappa)+normal_lpdf(c[i] | uc[i], kappa)+normal_lpdf(p[i] | up[i], kappa);
        }
    }
 
 

}

