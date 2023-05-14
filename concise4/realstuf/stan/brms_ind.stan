// generated with brms 2.19.0
functions {
  /* zero-one-inflated beta log-PDF of a single response
   * Args:
   *   y: response value
   *   mu: mean parameter of the beta part
   *   phi: precision parameter of the beta part
   *   zoi: zero-one-inflation probability
   *   coi: conditional one-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real zero_one_inflated_beta_lpdf(real y, real mu, real phi,
                                    real zoi, real coi) {
     row_vector[2] shape = [mu * phi, (1 - mu) * phi];
     if (y == 0) {
       return bernoulli_lpmf(1 | zoi) + bernoulli_lpmf(0 | coi);
     } else if (y == 1) {
       return bernoulli_lpmf(1 | zoi) + bernoulli_lpmf(1 | coi);
     } else {
       return bernoulli_lpmf(0 | zoi) + beta_lpdf(y | shape[1], shape[2]);
     }
   }
   
   
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_cold;  // number of observations
  vector[N_cold] Y_cold;  // response variable
  int<lower=1> K_cold;  // number of population-level effects
  matrix[N_cold, K_cold] X_cold;  // population-level design matrix
  int<lower=1> K_zoi_cold;  // number of population-level effects
  matrix[N_cold, K_zoi_cold] X_zoi_cold;  // population-level design matrix
  int<lower=1> N_warm;  // number of observations
  vector[N_warm] Y_warm;  // response variable
  int<lower=1> K_warm;  // number of population-level effects
  matrix[N_warm, K_warm] X_warm;  // population-level design matrix
  int<lower=1> K_zoi_warm;  // number of population-level effects
  matrix[N_warm, K_zoi_warm] X_zoi_warm;  // population-level design matrix
  int<lower=1> N_pain;  // number of observations
  vector[N_pain] Y_pain;  // response variable
  int<lower=1> K_pain;  // number of population-level effects
  matrix[N_pain, K_pain] X_pain;  // population-level design matrix
  int<lower=1> K_zoi_pain;  // number of population-level effects
  matrix[N_pain, K_zoi_pain] X_zoi_pain;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc_cold = K_cold - 1;
  matrix[N_cold, Kc_cold] Xc_cold;  // centered version of X_cold without an intercept
  vector[Kc_cold] means_X_cold;  // column means of X_cold before centering
  int Kc_zoi_cold = K_zoi_cold - 1;
  matrix[N_cold, Kc_zoi_cold] Xc_zoi_cold;  // centered version of X_zoi_cold without an intercept
  vector[Kc_zoi_cold] means_X_zoi_cold;  // column means of X_zoi_cold before centering
  int Kc_warm = K_warm - 1;
  matrix[N_warm, Kc_warm] Xc_warm;  // centered version of X_warm without an intercept
  vector[Kc_warm] means_X_warm;  // column means of X_warm before centering
  int Kc_zoi_warm = K_zoi_warm - 1;
  matrix[N_warm, Kc_zoi_warm] Xc_zoi_warm;  // centered version of X_zoi_warm without an intercept
  vector[Kc_zoi_warm] means_X_zoi_warm;  // column means of X_zoi_warm before centering
  int Kc_pain = K_pain - 1;
  matrix[N_pain, Kc_pain] Xc_pain;  // centered version of X_pain without an intercept
  vector[Kc_pain] means_X_pain;  // column means of X_pain before centering
  int Kc_zoi_pain = K_zoi_pain - 1;
  matrix[N_pain, Kc_zoi_pain] Xc_zoi_pain;  // centered version of X_zoi_pain without an intercept
  vector[Kc_zoi_pain] means_X_zoi_pain;  // column means of X_zoi_pain before centering
  
  for (i in 2:K_cold) {
    means_X_cold[i - 1] = mean(X_cold[, i]);
    Xc_cold[, i - 1] = X_cold[, i] - means_X_cold[i - 1];
  }
  for (i in 2:K_zoi_cold) {
    means_X_zoi_cold[i - 1] = mean(X_zoi_cold[, i]);
    Xc_zoi_cold[, i - 1] = X_zoi_cold[, i] - means_X_zoi_cold[i - 1];
  }
  for (i in 2:K_warm) {
    means_X_warm[i - 1] = mean(X_warm[, i]);
    Xc_warm[, i - 1] = X_warm[, i] - means_X_warm[i - 1];
  }
  for (i in 2:K_zoi_warm) {
    means_X_zoi_warm[i - 1] = mean(X_zoi_warm[, i]);
    Xc_zoi_warm[, i - 1] = X_zoi_warm[, i] - means_X_zoi_warm[i - 1];
  }
  for (i in 2:K_pain) {
    means_X_pain[i - 1] = mean(X_pain[, i]);
    Xc_pain[, i - 1] = X_pain[, i] - means_X_pain[i - 1];
  }
  for (i in 2:K_zoi_pain) {
    means_X_zoi_pain[i - 1] = mean(X_zoi_pain[, i]);
    Xc_zoi_pain[, i - 1] = X_zoi_pain[, i] - means_X_zoi_pain[i - 1];
  }
}
parameters {
  vector[Kc_cold] b_cold;  // population-level effects
  real Intercept_cold;  // temporary intercept for centered predictors
  real<lower=0> phi_cold;  // precision parameter
  vector[Kc_zoi_cold] b_zoi_cold;  // population-level effects
  real Intercept_zoi_cold;  // temporary intercept for centered predictors
  real<lower=0,upper=1> coi_cold;  // conditional one-inflation probability
  vector[Kc_warm] b_warm;  // population-level effects
  real Intercept_warm;  // temporary intercept for centered predictors
  real<lower=0> phi_warm;  // precision parameter
  vector[Kc_zoi_warm] b_zoi_warm;  // population-level effects
  real Intercept_zoi_warm;  // temporary intercept for centered predictors
  real<lower=0,upper=1> coi_warm;  // conditional one-inflation probability
  vector[Kc_pain] b_pain;  // population-level effects
  real Intercept_pain;  // temporary intercept for centered predictors
  real<lower=0> phi_pain;  // precision parameter
  vector[Kc_zoi_pain] b_zoi_pain;  // population-level effects
  real Intercept_zoi_pain;  // temporary intercept for centered predictors
  real<lower=0,upper=1> coi_pain;  // conditional one-inflation probability
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += student_t_lpdf(Intercept_cold | 3, 0, 2.5);
  lprior += gamma_lpdf(phi_cold | 0.01, 0.01);
  lprior += logistic_lpdf(Intercept_zoi_cold | 0, 1);
  lprior += beta_lpdf(coi_cold | 1, 1);
  lprior += student_t_lpdf(Intercept_warm | 3, 0, 2.5);
  lprior += gamma_lpdf(phi_warm | 0.01, 0.01);
  lprior += logistic_lpdf(Intercept_zoi_warm | 0, 1);
  lprior += beta_lpdf(coi_warm | 1, 1);
  lprior += student_t_lpdf(Intercept_pain | 3, 0, 2.5);
  lprior += gamma_lpdf(phi_pain | 0.01, 0.01);
  lprior += logistic_lpdf(Intercept_zoi_pain | 0, 1);
  lprior += beta_lpdf(coi_pain | 1, 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N_cold] mu_cold = rep_vector(0.0, N_cold);
    // initialize linear predictor term
    vector[N_cold] zoi_cold = rep_vector(0.0, N_cold);
    // initialize linear predictor term
    vector[N_warm] mu_warm = rep_vector(0.0, N_warm);
    // initialize linear predictor term
    vector[N_warm] zoi_warm = rep_vector(0.0, N_warm);
    // initialize linear predictor term
    vector[N_pain] mu_pain = rep_vector(0.0, N_pain);
    // initialize linear predictor term
    vector[N_pain] zoi_pain = rep_vector(0.0, N_pain);
    mu_cold += Intercept_cold + Xc_cold * b_cold;
    zoi_cold += Intercept_zoi_cold + Xc_zoi_cold * b_zoi_cold;
    mu_warm += Intercept_warm + Xc_warm * b_warm;
    zoi_warm += Intercept_zoi_warm + Xc_zoi_warm * b_zoi_warm;
    mu_pain += Intercept_pain + Xc_pain * b_pain;
    zoi_pain += Intercept_zoi_pain + Xc_zoi_pain * b_zoi_pain;
    mu_cold = inv_logit(mu_cold);
    zoi_cold = inv_logit(zoi_cold);
    mu_warm = inv_logit(mu_warm);
    zoi_warm = inv_logit(zoi_warm);
    mu_pain = inv_logit(mu_pain);
    zoi_pain = inv_logit(zoi_pain);
    for (n in 1:N_cold) {
      target += zero_one_inflated_beta_lpdf(Y_cold[n] | mu_cold[n], phi_cold, zoi_cold[n], coi_cold);
    }
    for (n in 1:N_warm) {
      target += zero_one_inflated_beta_lpdf(Y_warm[n] | mu_warm[n], phi_warm, zoi_warm[n], coi_warm);
    }
    for (n in 1:N_pain) {
      target += zero_one_inflated_beta_lpdf(Y_pain[n] | mu_pain[n], phi_pain, zoi_pain[n], coi_pain);
    }
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  
  vector[N] log_lik;
  // actual population-level intercept
  real b_cold_Intercept = Intercept_cold - dot_product(means_X_cold, b_cold);
  // actual population-level intercept
  real b_zoi_cold_Intercept = Intercept_zoi_cold - dot_product(means_X_zoi_cold, b_zoi_cold);
  // actual population-level intercept
  real b_warm_Intercept = Intercept_warm - dot_product(means_X_warm, b_warm);
  // actual population-level intercept
  real b_zoi_warm_Intercept = Intercept_zoi_warm - dot_product(means_X_zoi_warm, b_zoi_warm);
  // actual population-level intercept
  real b_pain_Intercept = Intercept_pain - dot_product(means_X_pain, b_pain);
  // actual population-level intercept
  real b_zoi_pain_Intercept = Intercept_zoi_pain - dot_product(means_X_zoi_pain, b_zoi_pain);
  

  // initialize linear predictor term
  vector[N_cold] mu_cold = rep_vector(0.0, N_cold);
  // initialize linear predictor term
  vector[N_cold] zoi_cold = rep_vector(0.0, N_cold);
  // initialize linear predictor term
  vector[N_warm] mu_warm = rep_vector(0.0, N_warm);
  // initialize linear predictor term
  vector[N_warm] zoi_warm = rep_vector(0.0, N_warm);
  // initialize linear predictor term
  vector[N_pain] mu_pain = rep_vector(0.0, N_pain);
  // initialize linear predictor term
  vector[N_pain] zoi_pain = rep_vector(0.0, N_pain);
  mu_cold += Intercept_cold + Xc_cold * b_cold;
  zoi_cold += Intercept_zoi_cold + Xc_zoi_cold * b_zoi_cold;
  mu_warm += Intercept_warm + Xc_warm * b_warm;
  zoi_warm += Intercept_zoi_warm + Xc_zoi_warm * b_zoi_warm;
  mu_pain += Intercept_pain + Xc_pain * b_pain;
  zoi_pain += Intercept_zoi_pain + Xc_zoi_pain * b_zoi_pain;
  mu_cold = inv_logit(mu_cold);
  zoi_cold = inv_logit(zoi_cold);
  mu_warm = inv_logit(mu_warm);
  zoi_warm = inv_logit(zoi_warm);
  mu_pain = inv_logit(mu_pain);
  zoi_pain = inv_logit(zoi_pain);
  
  
  for (n in 1:N_cold) {
    
    log_lik[n] = zero_one_inflated_beta_lpdf(Y_cold[n] | mu_cold[n], phi_cold, zoi_cold[n], coi_cold)+zero_one_inflated_beta_lpdf(Y_warm[n] | mu_warm[n], phi_warm, zoi_warm[n], coi_warm)+zero_one_inflated_beta_lpdf(Y_pain[n] | mu_pain[n], phi_pain, zoi_pain[n], coi_pain);
  }




  
  
}

