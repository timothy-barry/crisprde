data {
  // COUNTS
  int<lower=1> m; // number of mutation types
  int<lower=1> q; // total number of alleles
  int<lower=1> p; // dimension of the regression coefficient
  int<lower=1> r_cntrl; // number of control replicates
  int<lower=1> r_trt; // number of treated replicates
  array[r_cntrl, q + 1] int<lower=0> Y_cntrl; // control count matrix
  array[r_trt, q + 1] int<lower=0> Y_trt; // treated count matrix

  // FIXED DM LIKELIHOOD OVERDISPERSION RHO
  real<lower=1e-9, upper=1 - 1e-9> rho; // DM overdispersion parameter

  // TYPE STRUCTURE FOR THE FLATTENED MUTATED-ALLELE VECTOR
  array[m] int<lower=1> q_t; // number of distinct alleles per type
  array[m] int<lower=1, upper=q> type_start; // within mutated subvector, vector of type starting positions (first entry should be 1)
  array[m] int<lower=1, upper=q> type_end; // within mutated subvector, vector of type ending positions (final entry should be q)

  // DESIGN MATRIX FOR MUTATED ALLELES
  matrix[q, p] X; // stacked design matrix; alleles of type 1 form the first q_1 rows; alleles of type 2 form the second q_2 rows, etc.

  // PRIORS
  // BACKGROUND SPECTRUM
  simplex[m + 1] mu_pi_block; // mutation rate in control condition across types; first element is background mutation rate, second element is type-1 mutation rate, etc.
  real<lower=1e-9> kappa_pi_block; // Dirichlet concentration parameter kappa

  // PRIOR FOR MARGINAL EDITING RATE
  real<lower=1e-9> alpha_theta; // marginal editing rate alpha
  real<lower=1e-9> beta_theta; // marginal editing rate beta

  // PRIORS FOR BLOCK-LEVEL EDITING SPECTRUM
  simplex[m] mu_phi_block; // editing rate across types; first element is type-1 mutation contribution, second element is type-2 mutation contribution, etc.
  real<lower=1e-9> kappa_phi_block; // Dirichlet concentration parameter kappa

  // PER-TYPE PRIORS FOR CONTROL-SIDE REGRESSION
  matrix[m, p] mu_gamma; // control condition regression parameters
  matrix<lower=0>[m, p] sigma_gamma; // control condition regression parameter variances

  // PER-TYPE PRIORS FOR TREATED-SIDE REGRESSION
  matrix[m, p] mu_delta; // treated conition editing rate parameters
  matrix<lower=0>[m, p] sigma_delta; // treated condition editing rate parameter variances
}

transformed data {
  real<lower=0> kappa_lik;
  kappa_lik = (1 - rho) / rho; // Kappa parameter for the DM likelihood
}

parameters {
  simplex[m + 1] pi_block; // control mutation rate over mutation types (first entry is unmutated control rate)
  real<lower=1e-9, upper=1 - 1e-9> theta; // marginal editing rate
  simplex[m] phi_block; // editing rate over types
  matrix[m, p] gamma; // control condition regression coefficients
  matrix[m, p] delta; // treated condition regression coefficients
}

transformed parameters {
  vector<lower=0, upper=1>[q] psi; // flattened vector of type-wise background mutation spectra (does not sum to 1 for m >= 2)
  vector<lower=0, upper=1>[q] phi; // flattened vector of type-wise editing spectra (does not sum to 1 m >= 2)
  simplex[q + 1] pi; // flattened vector of background mutation rates (sums to 1); first entry unmutated
  simplex[q + 1] tau; // flattened vector of treated mutation rates (sums to 1); first entry unmutated
  real log_pi0;
  real log_theta;

  log_pi0 = log(pi_block[1]);
  log_theta = log(theta);

  pi[1] = pi_block[1]; //  background unmutated rate
  tau[1] = exp(log_pi0 + log1m(theta)); // treated unmutated rate

  // for each mutation type
  for (t in 1:m) {
    // initialize key integers and vectors
    int edited_start;
    int edited_stop;
    int full_start;
    int full_stop;
    vector[q_t[t]] linpred_gamma;
    vector[q_t[t]] linpred_delta;
    vector[q_t[t]] log_psi_t;
    vector[q_t[t]] log_phi_t;
    vector[q_t[t]] log_pi_t;
    vector[q_t[t]] log_tau_t;
    vector[q_t[t]] psi_t;
    vector[q_t[t]] phi_t;
    real log_pi_block_t;
    real log_phi_block_t;

    edited_start = type_start[t]; // get the subvector start index
    edited_stop = type_end[t]; // get the subvector end index
    full_start = edited_start + 1; // get start relative to full vector
    full_stop = edited_stop + 1; // get stop relative to full vector

    linpred_gamma = X[edited_start:edited_stop, ] * gamma[t]'; // gamma linear predictor across alleles
    linpred_delta = X[edited_start:edited_stop, ] * delta[t]'; // delta linear predictor across alleles

    log_psi_t = log_softmax(linpred_gamma); // log background mutation spectrum
    log_phi_t = log_softmax(linpred_delta); // log editing spectrum
    psi_t = exp(log_psi_t); // background mutation spectrum
    phi_t = exp(log_phi_t); // editing spectrum
    log_pi_block_t = log(pi_block[t + 1]);
    log_phi_block_t = log(phi_block[t]);
    log_pi_t = log_pi_block_t + log_psi_t;

    for (a in 1:q_t[t]) {
      log_tau_t[a] = log_sum_exp(log_pi_t[a], log_pi0 + log_theta + log_phi_block_t + log_phi_t[a]);
    }

    psi[edited_start:edited_stop] = psi_t; // fill component of psi
    phi[edited_start:edited_stop] = phi_t; // fill component of phi

    pi[full_start:full_stop] = exp(log_pi_t); // fill component of pi with background editing rate
    tau[full_start:full_stop] = exp(log_tau_t); // fill component of tau with treated mutation rate
  }
}

model {
  // initialize vectors needed for stan DM parameterization
  vector[q + 1] alpha_dm_cntrl;
  vector[q + 1] alpha_dm_trt;
  real eps_alpha;

  // PRIORS
  pi_block ~ dirichlet(kappa_pi_block * mu_pi_block);
  phi_block ~ dirichlet(kappa_phi_block * mu_phi_block);
  theta ~ beta(alpha_theta, beta_theta);
  for (t in 1:m) {
    gamma[t] ~ normal(mu_gamma[t], sigma_gamma[t]);
    delta[t] ~ normal(mu_delta[t], sigma_delta[t]);
  }

  // LIKELIHOOD
  eps_alpha = 1e-12;
  alpha_dm_cntrl = kappa_lik * pi + rep_vector(eps_alpha, q + 1);
  alpha_dm_trt = kappa_lik * tau + rep_vector(eps_alpha, q + 1);
  for (i in 1:r_cntrl) {
    Y_cntrl[i] ~ dirichlet_multinomial(alpha_dm_cntrl);
  }
  for (i in 1:r_trt) {
    Y_trt[i] ~ dirichlet_multinomial(alpha_dm_trt);
  }
}

generated quantities {
  vector<lower=0, upper=1>[m] theta_tilde_block;
  vector<lower=0, upper=1>[q] theta_tilde;
  array[r_cntrl] real log_lik_cntrl;
  array[r_trt] real log_lik_trt;
  vector[q + 1] alpha_dm_cntrl;
  vector[q + 1] alpha_dm_trt;
  real eps_alpha;

  // blocked theta tilde (i.e., overall editing rate for each mutation type; length m)
  theta_tilde_block = theta * to_vector(phi_block);

  // allele-level theta tilde (i.e., allele-specific editing rate for each allele; length q)
  for (t in 1:m) {
    int edited_start;
    int edited_stop;
    edited_start = type_start[t];
    edited_stop = type_end[t];
    theta_tilde[edited_start:edited_stop] = theta * phi_block[t] * phi[edited_start:edited_stop];
  }

  // observation-wise log likelihood
  eps_alpha = 1e-12;
  alpha_dm_cntrl = kappa_lik * pi + rep_vector(eps_alpha, q + 1);
  alpha_dm_trt = kappa_lik * tau + rep_vector(eps_alpha, q + 1);

  for (i in 1:r_cntrl) {
    log_lik_cntrl[i] = dirichlet_multinomial_lpmf(Y_cntrl[i] | alpha_dm_cntrl);
  }
  for (i in 1:r_trt) {
    log_lik_trt[i] = dirichlet_multinomial_lpmf(Y_trt[i] | alpha_dm_trt);
  }
}
