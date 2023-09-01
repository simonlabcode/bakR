// Used to be called Full_Model.stan
data {
  int NE; //num entries in data frame.
  int NF; //num features
  array[NE] int TP; // type for each entry (+s4U = 1, -s4U = 0)
  array[NE] int MT; // mutant type for each entry (WT = 1, treatment A = 2, treatment B = 3, treatment C = 4)
  int nMT; // number of different mutation types
  array[NE] int FE; //  feature  for each entry
  array[NE] int num_mut; //  mutations in each entry
  array[NE] int num_obs; // number of times its observed
  array[NE] int R; //Replicate ID
  array[NE] real U_cont; // Feature specific log-fold difference in U content
  int nrep;  // Number of replicates
  array[nMT] real tl;   // Label time
  array[NF, nMT] real Avg_Reads; // Average read counts in transcript i and sample j (avg. across replicates)
  int<lower = 0, upper = 1> Chase;
  real nU;
}

parameters {
  array[nMT] vector<lower=0>[nrep] TL_lambda_eff;   // s4U induced increase in mutation rate
  array[nMT] vector[nrep] log_lambda_o;  // Background mutation rate
  array[NF] vector[nMT] alpha;
  array[nMT] real mu_fn;        // Global mean fn
  array[nMT] real<lower=0> sig_fn;
  array[NF, nMT, nrep] real z_fn;
  //vector[nMT-1] z_e[NF];
  vector<lower=0>[nMT] a;
  vector[nMT] b;
  array[NF] vector[nMT] z_sd_r ;
  vector<lower=0>[nMT] sig_rep;
  //vector<lower=0>[nMT] sd_r_mu [NF];
}

transformed parameters {
  array[NF, nMT, nrep] real<lower=0,upper=1> frac_new;   // Fraction new (fn) of obs reads on native scale.
  array[NF, nMT, nrep] real mu_rep_logit_fn; // Fn on logit scale
  array[nMT] vector[nrep] log_lambda_n; // Mutation rate of s4U labelled reads
  array[NF] vector<lower=0>[nMT] sd_r_mu;
  array[NF] vector[nMT] sd_rep;

      for(m in 1:nMT){
        for(r in 1:nrep){
          log_lambda_n[m,r] = log_lambda_o[m,r] + TL_lambda_eff[m, r];

        }
      }

      for(i in 1:NF){
      for(j in 1:nMT){
        sd_rep[i,j] = -a[j]*Avg_Reads[i,j] + b[j];
        sd_r_mu[i, j] = exp(sd_rep[i,j] + z_sd_r[i,j]*sig_rep[j]);
        // if(j > 1){
        //   eff[i, j-1] = mu_e[j-1] + sig_e[j-1]*z_e[i, j-1];
        // }
        for(k in 1:nrep){
          mu_rep_logit_fn[i, j, k] = alpha[i,j] + sd_r_mu[i,j]*z_fn[i,j,k];
          frac_new[i, j, k] = inv_logit(mu_rep_logit_fn[i,j,k]);
        }
      }
      }
}

model {
  // Priors:

  mu_fn ~ normal(0, 1.25);
  sig_fn ~ lognormal(-1, 0.5);

  a ~ normal(0.3, 0.2);
  b ~ normal(-1.5, 0.35);


  for(m in 1:nMT){
      TL_lambda_eff[m] ~ lognormal(1.4, 0.15);
      log_lambda_o[m] ~ normal(log(0.001*nU), 0.5);
  }

  sig_rep ~ lognormal(-2,0.25);

  for (i in 1:NF) {

    // Remove if non-centered
    //eff[i] ~ normal(mu_e, sig_e);

    z_sd_r[i] ~ normal(0,1);
    for(j in 1:nMT){
        alpha[i,j] ~ normal(mu_fn[j], sig_fn[j]);

        //sd_r_mu[i,j] ~ exponential(1/exp(-a[j]*log10(Avg_Reads[i,j]) + b[j]));
        z_fn[i,j] ~ normal(0, 1);

        //sd_r_mu[i,j] ~ lognormal(sd_rep[i,j], sig_rep[j]);


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

  array[NF] vector[nMT] kd; // decay rate
  array[NF] vector[nMT-1] L2FC_kd; // change in decay rate
  array[NF] vector[nMT] log_kd; //log(kdeg)

  for (i in 1:NF) {
    for (j in 1:nMT) {

      if(Chase == 1){
        kd[i,j] = -log(inv_logit(alpha[i,j]))/tl[j];  // divide by time if not 1

        log_kd[i,j] = log(kd[i,j]);

      }else{
        kd[i,j] = -log(1-inv_logit(alpha[i,j]))/tl[j];  // divide by time if not 1

        log_kd[i,j] = log(kd[i,j]);

      }


    } // treatment type
    for(k in 1:nMT-1){
      L2FC_kd[i,k] = log2(kd[i,k+1]/kd[i,1]);
    }
  }

}
