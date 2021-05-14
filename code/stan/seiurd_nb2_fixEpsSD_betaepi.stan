// Single Region
// Revised Model Structure
//
// Constant Params:
// tau, alpha, nu, kappa,
// gamma0, delta0, rho0, eta0
//
// Parameters with regression model:
// beta, omega (and thus gamma and delta)
//
// Psi = a - b *exp(-c*t) 
//
// Fixed value of v (sd of epsilon's)
//
// Note: no epsilonD[1] -- modify this if there are any known UD on day 0
//

#include functions/functions_SEIURD.stan

data {
  int<lower=2> nT;  // number of time points observed
  int yc[nT]; // number of new cases at time t
  int yd[nT]; // number of new deaths at time t
  int<lower=1> betaU_d; // number of splines/predictors for betaU
  matrix[nT, betaU_d] betaU_Mx; // spline/predictor values for betaU
  int<lower=1> omega_d; // number of splines/predictors in nu
  matrix[nT, omega_d] omega_Mx; // spline/predictor values for nu

  // initial values for latent states (S, E, IU, RU, ID, UD, RD, DD)
  vector<lower=0>[8] state_init;

  // Prior information
  // betaU = exp(betaU_b0 + betaU_Mx %*% betaU_b )
  // betaU_b0 ~ N(betaU_b0_mean, betaU_b0_sd)
  real betaU_b0_mean;
  real<lower=0> betaU_b0_sd;
  // betaU_b ~ N(betaU_b_mu, betaU_b_sigma)
  // betaU_b_sigma ~ N+(betaU_b_sigma_mean, betaU_b_sigma_sd)
  real betaU_b_mu;
  real<lower=0> betaU_b_sigma_mean;
  real<lower=0> betaU_b_sigma_sd;
  // Prior for psi elements 
  real<lower=0,upper=1> psi_a_limit_low; //
  real<lower=0,upper=1> psi_a_limit_high; //
  real<lower=0> psi_a_a; // a ~ Beta(a,b)
  real<lower=0> psi_a_b;
  real<lower=0> psi_b0_a; // b0 ~ Beta(a, b)
  real<lower=0> psi_b0_b;
  real<lower=0> psi_c_a; // c ~ Gamma(a, b)
  real<lower=0> psi_c_b;
  // Omega prior, logit scale regression model
  real omega_b0_mean;
  real<lower=0> omega_b0_sd;
  real omega_b_mu;
  real<lower=0> omega_b_sigma_mean;
  real<lower=0> omega_b_sigma_sd;
  // others
  real<lower=0> alpha_a; // alpha prior is Beta(a,b)
  real<lower=0> alpha_b; //
  real<lower=0,upper=1> alpha_limit_low; //
  real<lower=0,upper=1> alpha_limit_high; //
  real<lower=0> eta0_a; // eta0 is Beta(a,b)
  real<lower=0> eta0_b;
  real<lower=0,upper=1> eta0_limit_low; //
  real<lower=0,upper=1> eta0_limit_high; //
  real<lower=0> gamma0_a; // gamma0 is Beta(a,b)
  real<lower=0> gamma0_b;
  real<lower=0,upper=1> gamma0_limit_low; //
  real<lower=0,upper=1> gamma0_limit_high; //
  real<lower=0> delta0_a; // delta0 is Beta(a,b)
  real<lower=0> delta0_b; //
  real<lower=0,upper=1> delta0_limit_low; //
  real<lower=0,upper=1> delta0_limit_high; //
  real<lower=0> nu_a; // nu is Beta(a,b)
  real<lower=0> nu_b; //
  real<lower=0,upper=1> nu_limit_low; //
  real<lower=0,upper=1> nu_limit_high; //
  real<lower=0> tau_a; // tau prior is beta(a, b)
  real<lower=0> tau_b;
  real<lower=0> kappa_a; // kappa prior is gamma(a, b)
  real<lower=0> kappa_b;
  real<lower=0> EEIU0_a; // EEIU0 prior is beta(a, b)
  real<lower=0> EEIU0_b;
  real<lower=0> epsilonI_sd; // FIXED quantity
  real<lower=0> epsilonD_sd; // FIXED quantity
  real<lower=0> phiC_a; // phiC has a gamma(a,b) prior
  real<lower=0> phiC_b;
  real<lower=0> phiD_a; // phiD has a gamma(a,b) prior
  real<lower=0> phiD_b;
}

parameters {
  // ODE model parameters
  real betaU_b0_raw;
  vector[betaU_d] betaU_b_raw;
  real<lower=0> betaU_b_sigma;
  real<lower=psi_a_limit_low,upper=psi_a_limit_high> psi_a;
  real<lower=0,upper=1> psi_b0;
  real<lower=0> psi_c;
  real omega_b0_raw;
  vector[omega_d] omega_b_raw;
  real<lower=0> omega_b_sigma;
  real<lower=alpha_limit_low,upper=alpha_limit_high> alpha;
  real<lower=eta0_limit_low,upper=eta0_limit_high> eta0;
  real<lower=nu_limit_low,upper=nu_limit_high> nu;
  real<lower=gamma0_limit_low,upper=gamma0_limit_high> gamma0;
  real<lower=delta0_limit_low,upper=delta0_limit_high> delta0;
  real<lower=0,upper=1> tau;
  real<lower=0> kappa;
  real<lower=0, upper=1> EEIU0; // E(0)/[E(0) + IN(0)]

  // Process model parameters
  vector<lower=0>[nT] epsilonI; // 1-dimensional version
  vector<lower=0>[nT-1] epsilonD; // 1-dimensional version
  real<lower=0> phiC;
  real<lower=0> phiD;
}

transformed parameters {
  // ODE model parameters
  vector<lower=0>[8] state[nT]; // process states at each time point
  vector<lower=0>[8] state_init_mod;
  // vector[8] state_change; // temporary holder for state updates
  real<lower=0> EIN0; // E(0) + IN(0)
  real betaU_b0;
  vector[betaU_d] betaU_b;
  vector<lower=0>[nT] betaU;
  real<lower=0, upper=psi_a> psi_b;
  vector<lower=0,upper=1>[nT] psi;
  vector<lower=0,upper=1>[nT] eta;
  real omega_b0;
  vector[omega_d] omega_b;
  vector<lower=0,upper=1>[nT] omega;
  real<lower=0,upper=1> rho0;
  vector<lower=0,upper=1>[nT] rho;
  vector<lower=0,upper=1>[nT] gamma;
  vector<lower=0,upper=1>[nT] delta;
  real<lower=0> theta[8]; // (betaU, betaD, alpha, rho, eta, nu, gamma, delta)
  // Process model parameters
  real extra_count;
  // Data model parameters
  vector<lower=0>[nT] muC;
  vector<lower=0>[nT] muD;
  vector<lower=0>[nT] muCorig;
  vector<lower=0>[nT] muDorig;
  // vector<lower=0>[nT] phiC;
  // vector<lower=0>[nT] phiD;

  // Compute betaU_b as transformation of
  // the "raw" values which are N(0, 1)
  betaU_b0=betaU_b0_mean + betaU_b0_raw*betaU_b0_sd;
  betaU_b=betaU_b_mu + betaU_b_raw*betaU_b_sigma;
  betaU=exp(betaU_b0 + betaU_Mx*betaU_b);
  // eta, rho
  psi_b=psi_b0*psi_a;
  psi[1] = psi_a - psi_b*exp(-psi_c*1);
  eta[1]=psi[1]*eta0;
  rho0 = eta0*nu / (eta0 + nu);
  rho[1]=(1-psi[1])*rho0;
    // omega, gamma, delta
  omega_b0=omega_b0_mean + omega_b0_raw*omega_b0_sd;
  omega_b=omega_b_mu + omega_b_raw*omega_b_sigma;
  omega=inv_logit(omega_b0 + omega_Mx*omega_b);
  gamma = (1-omega)*gamma0;
  delta = omega*delta0;

  theta[1] = betaU[1];
  theta[2] = tau*betaU[1];
  theta[3] = alpha;
  theta[4] = rho[1];
  theta[5] = eta[1];
  theta[6] = nu;
  theta[7] = gamma[1];
  theta[8] = delta[1];

  state_init_mod = state_init;
  EIN0 = kappa*state_init[5];
  state_init_mod[1] -= EIN0;
  state_init_mod[2] = EEIU0*EIN0;
  state_init_mod[3] = (1 - EEIU0)*EIN0;

  // update states with ODE solution
  state[1] = state_init_mod + SEIURD(state_init_mod, theta);
  muCorig[1] = compute_new_cases(state_init_mod,state[1,]);
  muDorig[1] = compute_new_deaths(state_init_mod,state[1,]);

  // add multiplicative process on IP
  extra_count = muCorig[1]*(epsilonI[1] - 1);
  muC[1] = muCorig[1] + extra_count;
  state[1,5] += extra_count;
  state[1,3] -= extra_count;
  
  muD[1] = muDorig[1];

  // repeat for each time point
  for (t in 2:nT){
    psi[t] = psi_a - psi_b*exp(-psi_c*t);
    eta[t]=psi[t]*eta0;
    rho[t]=(1-psi[t])*rho0;
    
    theta[1] = betaU[t];
    theta[2] = tau*betaU[t];
    // theta[3] = alpha;
    theta[4] = rho[t];
    theta[5] = eta[t];
    // theta[6] = nu;
    theta[7] = gamma[t];
    theta[8] = delta[t];
    state[t] = state[t-1] + SEIURD(state[t-1], theta);
    muCorig[t] = compute_new_cases(state[t-1],state[t]);
    muDorig[t] = compute_new_deaths(state[t-1],state[t]);
    // add multiplicative process on IP
    extra_count = muCorig[t]*(epsilonI[t] - 1);
    muC[t] = muCorig[t] + extra_count;
    state[t,5] += extra_count;
    state[t,3] -= extra_count;
    // add multiplicative process on DD
    extra_count = muDorig[t]*(epsilonD[t-1] - 1);
    muD[t] = muDorig[t] + extra_count;
    state[t,8] += extra_count;
    state[t,6] -= extra_count;
    
    
  }
}

model {
  // Data model
  target += neg_binomial_2_lpmf(yc | muC, phiC);
  // target += neg_binomial_2_lpmf(yd[2:] | muD[2:], phiD[2:]);
  target += neg_binomial_2_lpmf(yd[2:] | muD[2:], phiD);
  
  target += gamma_lpdf(phiC | phiC_a, phiC_b);
  target += gamma_lpdf(phiD | phiD_a, phiD_b);

  // Process model
  target += gamma_lpdf(epsilonI | 1/epsilonI_sd^2, 1/epsilonI_sd^2);
  target += gamma_lpdf(epsilonD | 1/epsilonD_sd^2, 1/epsilonD_sd^2);

  // ODE model parameters
  target += normal_lpdf(betaU_b0_raw | 0, 1);
  target += normal_lpdf(betaU_b_raw | 0, 1);
  target += normal_lpdf(betaU_b_sigma | betaU_b_sigma_mean, betaU_b_sigma_sd);
  
  target += beta_lpdf(psi_a | psi_a_a, psi_a_b);
  target += beta_lpdf(psi_b0 | psi_b0_a, psi_b0_b);
  target += gamma_lpdf(psi_c | psi_c_a, psi_c_b);
  
  target += normal_lpdf(omega_b0_raw | 0, 1);
  target += normal_lpdf(omega_b_raw | 0, 1);
  target += normal_lpdf(omega_b_sigma | omega_b_sigma_mean, omega_b_sigma_sd);
  target += beta_lpdf(tau | tau_a, tau_b);
  target += beta_lpdf(alpha | alpha_a, alpha_b);
  target += beta_lpdf(eta0 | eta0_a, eta0_b);
  target += beta_lpdf(nu | nu_a, nu_b);
  target += beta_lpdf(gamma0 | gamma0_a, gamma0_b);
  target += beta_lpdf(delta0 | delta0_a, delta0_b);

  // Initial state parameters
  target +=  gamma_lpdf(kappa | kappa_a, kappa_b);
  target +=  beta_lpdf(EEIU0 | EEIU0_a, EEIU0_b);
}

generated quantities {
  vector[nT] R0;
  vector[nT] R0_noepsilon;
  int yc_pred[nT];
  int yd_pred[nT-1];

  yc_pred = neg_binomial_2_rng(muC, phiC);
  // yd_pred = neg_binomial_2_rng(muD[2:], phiD[2:]);
  yd_pred = neg_binomial_2_rng(muD[2:], phiD);

  R0 = betaU ./ (eta .*epsilonI + rho) .* (1 + eta .* epsilonI * tau / nu);
  R0_noepsilon = betaU ./ (eta + rho) .* (1 + eta * tau / nu);
}

