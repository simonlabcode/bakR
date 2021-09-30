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
  real tl[nMT];   // Label time
  real Avg_Reads[NF, nMT]; // Average read counts in transcript i and sample j (avg. across replicates)
}

parameters {
  vector<lower=0>[nrep] TL_lambda_eff[nMT];   // s4U induced increase in mutation rate
  vector[nrep] log_lambda_o[nMT];  // Background mutation rate
  vector[NF] alpha;  // Reference sample mean logit(fn)
  real mu_fn;        // Global mean fn
  vector[nMT-1] mu_e;        // Global effect size mean
  real<lower=0> sig_fn;
  vector<lower=0>[nMT - 1] sig_e;
  real z_fn[NF, nMT, nrep];
  //vector[nMT-1] z_e[NF];
  vector<lower=0>[nMT] a;
  vector[nMT] b;
  //vector[nMT] z_sd_r [NF];
  //vector<lower=0>[nMT] sig_rep;
  vector[nMT-1] eff[NF];
  //vector<lower=0>[nMT] sd_r_mu [NF];
}

transformed parameters {
  real<lower=0,upper=1> frac_new[NF, nMT, nrep];   // Fraction new (fn) of obs reads on native scale.
  real mu_rep_logit_fn[NF, nMT, nrep]; // Fn on logit scale
  vector[nrep] log_lambda_n[nMT]; // Mutation rate of s4U labelled reads
  //vector<lower=0>[nMT] sd_r_mu[NF];
  vector<lower=0>[nMT] sd_rep[NF];

      for(i in 1:NF){
      for(j in 1:nMT){
        sd_rep[i,j] = exp(-a[j]*Avg_Reads[i,j] + b[j]);
        //sd_r_mu[i, j] = exp(sd_rep[i,j] + z_sd_r[i,j]*sig_rep[j]);
        // if(j > 1){
        //   eff[i, j-1] = mu_e[j-1] + sig_e[j-1]*z_e[i, j-1];
        // }
        for(k in 1:nrep){
          log_lambda_n[j,k] = log_lambda_o[j,k] + TL_lambda_eff[j, k];
          if(j == 1){
             mu_rep_logit_fn[i, j, k] = alpha[i] + sd_rep[i,j]*z_fn[i,j,k];
             frac_new[i, j, k] = inv_logit(mu_rep_logit_fn[i,j,k]);
          }else{
             mu_rep_logit_fn[i, j, k] = (alpha[i] + eff[i, j-1]) + sd_rep[i,j]*z_fn[i, j, k];
             frac_new[i, j, k] = inv_logit(mu_rep_logit_fn[i,j,k]);
          }
        }
      }
      }
}

model {
  // Priors:

  mu_fn ~ normal(0, 1.5);
  sig_fn ~ lognormal(-1, 0.75);
  alpha ~ normal(mu_fn, sig_fn);

  sig_e ~ lognormal(-1, 0.7);
  mu_e ~ normal(0, 1);

  a ~ normal(0.5, 0.25);
  b ~ normal(-1, 0.5);


  for(m in 1:nMT){
      TL_lambda_eff[m] ~ lognormal(-1, 0.5);
      log_lambda_o[m] ~ normal(-3.5, 0.5);
  }

  //sig_rep ~ lognormal(-1,0.5);

  for (i in 1:NF) {

    // Remove if non-centered
    //eff[i] ~ normal(mu_e, sig_e);

    //z_sd_r[i] ~ normal(0,1);
    for(j in 1:nMT){
        //sd_r_mu[i,j] ~ exponential(1/exp(-a[j]*log10(Avg_Reads[i,j]) + b[j]));
        z_fn[i,j] ~ normal(0, 1);

        //sd_r_mu[i,j] ~ lognormal(sd_rep[i,j], sig_rep[j]);

        if(j > 1){
          eff[i, j-1] ~ normal(mu_e[j-1], sig_e[j-1]);
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
        kd[i,j] = -log(1-inv_logit(alpha[i]))/tl[j];  // divide by time if not 1
      }else{
        kd[i,j] = -log(1-inv_logit(alpha[i] + eff[i, j-1]))/tl[j];
      }
    } // treatment type
    for(k in 1:nMT-1){
      L2FC_kd[i,k] = log2(kd[i,k+1]/kd[i,1]);
    }
  }

}


