data {
  int NE; //num entries in data frame.
  int NF; //num features
  int TP[NE]; // type for each entry (+s4U = 1, -s4U = 0)
  int MT[NE]; // mutant type for each entry (WT = 1, treatment A = 2, treatment B = 3, treatment C = 4)
  int nMT; // number of different mutation types
  int FE[NE]; //  feature  for each entry
  int num_mut[NE]; //  mutations in each entry
  int num_obs[NE]; // number of times its observed
  int R[NE]; //Replicate ID
  real U_cont[NE]; // Feature specific log-fold difference in U content
  int nrep;  // Number of replicates
  real tl;   // Label time
}

parameters {
  vector[NF] alpha;  // Reference sample mean logit(fn)
  real mu_fn;        // Global mean fn
  real log_sig_fn;   // Global variation in fn
  vector[nMT-1] log_sig_e;   // Global variation in effect sizes (change in logit(fn))
  vector[nMT-1] mu_e;        // Global effect size mean
  vector[nMT-1] z_e [NF];    // Effect size non-centered parameterization z-score
  real z_fn[NF, nMT, nrep];  // Fn non-centered parameterization z-score
  real log_sd_r_mu;          // Log replicate variability
  vector<lower=0>[nrep] TL_lambda_eff[nMT];   // s4U induced increase in mutation rate
  vector[nrep] log_lambda_o[nMT];  // Background mutation rate
}

transformed parameters {
  real frac_new[NF, nMT, nrep];   // Fraction new (fn) of obs reads on native scale.
  real mu_rep_logit_fn[NF, nMT, nrep]; // Fn on logit scale
  vector[nrep] log_lambda_n[nMT]; // Mutation rate of s4U labelled reads
  vector[nMT-1] eff[NF];           // Effect size (change in logit(fn) w.r.t. reference sample)
  real sig_fn = exp(log_sig_fn);   // Replicate variability
  real sd_r_mu = exp(log_sd_r_mu); //
  vector[nMT-1] sig_e = exp(log_sig_e);  // Global variability in effect size

    for (i in 1:NF) {
      for(j in 1:nMT){
        if(j > 1){
          eff[i, j-1] = mu_e[j-1] + z_e[i,j-1]*sig_e[j-1];
        }
        for(k in 1:nrep){
          log_lambda_n[j,k] = log_lambda_o[j,k] + TL_lambda_eff[j, k];
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

  for(m in 1:nMT){
    TL_lambda_eff[m] ~ lognormal(-1, 0.5);
    log_lambda_o[m] ~ normal(-3.5, 0.5);
  }

  mu_fn ~ normal(0, 1.5);
  log_sig_fn ~ normal(-1, 1);
  alpha ~ normal(mu_fn, sig_fn);
  log_sig_e ~ normal(-1, 1);
  mu_e ~ normal(0, 1);
  log_sd_r_mu ~ normal(-1,1);

  for (i in 1:NF) {
    //eff[i] ~ normal(mu_e, sig_e);
    z_e[i] ~ normal(0, 1);
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
                                   poisson_log_lpmf(num_mut[i] | (U_cont[i] + log_lambda_n[MT[i], R[i]])),
                                   poisson_log_lpmf(num_mut[i] | (U_cont[i] + log_lambda_o[MT[i], R[i]])));
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
