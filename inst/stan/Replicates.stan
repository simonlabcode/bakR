data {
  int NE; //num entries in data frame.
  int NF; //num features
  int TP[NE]; // type for each entry (+s4U = 1, -s4U = 0)
  int MT[NE]; // mutant type for each entry (WT = 1, treatment A = 2, treatment B = 3, treatment C = 4)
  int nMT; // number of different mutation types
  //int N_obs[NF, nMT]; // number of total reads for each feature in each mutation type
  int FE[NE]; //  feature  for each entry
  int num_mut[NE]; //  mutations in each entry
  int num_obs[NE]; // number of times its observed
  int R[NE]; //Replicate ID
  int nrep;
  real tl;
}

parameters {
  real log_lambda_o_m; // mean background mutation rate
  real log_lambda_o_sd; // sd for background mutations
  vector[NF] z_o; // z score for old mutation rate lvector
  real log_mu_TL;
  real log_sig_TL;
  //vector<lower=0>[nMT] TL_lambda_eff[NF];
  vector[nMT] z_TL [NF];
  vector[NF] alpha;
  real mu_fn;
  real log_sig_fn;
  vector[nMT-1] log_sig_e;
  vector[nMT-1] mu_e;
  vector[nMT-1] z_e [NF];
  //real mu_rep_logit_fn[NF, nMT, nrep]; // Inferred fraction new of obs reads on native scale.
  real z_fn[NF, nMT, nrep];
  real log_sd_r_mu;
  //real<lower=0> sd_rep[NF, nMT];
}

transformed parameters {
  real frac_new[NF, nMT, nrep];
  real mu_rep_logit_fn[NF, nMT, nrep]; // Inferred fraction new of obs reads on native scale.
  vector[NF] log_lambda_o; // Inferred background muatation rate per read
  vector[nMT] log_lambda_n[NF]; // Inferred s4U-induced muatation rate per read
  vector[nMT-1] eff[NF]; //  Parameter for fraction new of observed reads
  real lambda_o_sd = exp(log_lambda_o_sd);
  real sig_fn = exp(log_sig_fn);
  real sd_r_mu = exp(log_sd_r_mu);
  real mu_TL = exp(log_mu_TL);
  vector[nMT-1] sig_e = exp(log_sig_e);
  real sig_TL = exp(log_sig_TL);
  vector[nMT] TL_lambda_eff [NF];

  // Left here if needed for non-centered parameterization:
    for (i in 1:NF) {
      log_lambda_o[i] = log_lambda_o_m + (z_o[i] * lambda_o_sd);
      for(j in 1:nMT){
        TL_lambda_eff[i, j] = exp(log_mu_TL + z_TL[i,j]*sig_TL);
        log_lambda_n[i,j] = log(exp(log_lambda_o[i]) + TL_lambda_eff[i,j]);
        if(j > 1){
          eff[i, j-1] = mu_e[j-1] + z_e[i,j-1]*sig_e[j-1];
        }
        for(k in 1:nrep){
          if(j == 1){
            mu_rep_logit_fn[i, j, k] = alpha[i] + z_fn[i,j,k]*sd_r_mu;
            frac_new[i, j, k] = inv_logit(mu_rep_logit_fn[i,j, k]);
          }else{
            mu_rep_logit_fn[i, j, k] = (alpha[i] + eff[i, j - 1]) + z_fn[i,j,k]*sd_r_mu;
            frac_new[i, j, k] = inv_logit(mu_rep_logit_fn[i,j, k]);
          }

        }
      }
    }
}

model {
  // Priors:

    log_lambda_o_m ~ normal(-3.5, 0.5);
  log_lambda_o_sd ~ normal(-1, 1);

  z_o ~ normal(0, 1);

  log_mu_TL ~ normal(-1, 0.3);
  log_sig_TL ~ normal(-1, 1);

  mu_fn ~ normal(0, 1.5);
  log_sig_fn ~ normal(-1, 1);
  alpha ~ normal(mu_fn, sig_fn);
  log_sig_e ~ normal(-1, 1);
  mu_e ~ normal(0, 1);
  log_sd_r_mu ~ normal(-1,1);

  for (i in 1:NF) {
    //eff[i] ~ normal(mu_e, sig_e);
    z_e[i] ~ normal(0, 1);
    z_TL[i] ~ normal(0, 1);
    //sd_rep[i] ~ normal(sd_r_mu, 1);
    for(j in 1:nMT){

      for(k in 1:nrep){
        z_fn[i, j, k] ~ normal(0, 1);
      }
    }
  }

  // Model fraction of reads that are new.
  for (i in 1:NE) {
    target += num_obs[i] * log_mix(TP[i] * frac_new[FE[i], MT[i], R[i]],
                                   poisson_log_lpmf(num_mut[i] | log_lambda_n[FE[i], MT[i]]),
                                   poisson_log_lpmf(num_mut[i] | log_lambda_o[FE[i]]));
  } //reads
}

generated quantities {

  vector[nMT] kd[NF]; // decay rate
  vector[nMT-1] L2FC_kd[NF]; // change in decay rate

  for (i in 1:NF) {
    for (j in 1:nMT) {
      if(j == 1){
        kd[i,j] = -log(1-inv_logit(alpha[i]))/tl;  // divide by time if not 1
      }else{
        kd[i,j] = -log(1-inv_logit(alpha[i] + eff[i, j-1]))/tl;
      }
    } // treatment type
    for(k in 1:nMT-1){
      L2FC_kd[i,k] = log2(kd[i,k+1]/kd[i,1]);
    }
  }

}
