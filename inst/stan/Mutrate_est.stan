// Completely centered parameterization with no hierarchy
// Used to be called Mut_rate_Model.stan

// Data block corresponds to what gets passed to Stan in data_list
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
  int nrep;
  real tl[nMT];
}

// Parameter block tells model all the things it will estimate
parameters {
  vector[nrep] log_lambda_o[nMT]; // Background mutation rate
  vector<lower=0>[nrep] TL_lambda_eff[nMT]; // Increase to mutation rate due to s4U
  real mu_rep_logit_fn[NF, nMT, nrep]; // Inferred fraction new of obs reads on native scale.
}

// Transformed parameters are parameters that are defined in terms of data or other parameters
transformed parameters {
  real frac_new[NF, nMT, nrep] = inv_logit(mu_rep_logit_fn); // fraction new
  vector[nrep] log_lambda_n[nMT]; // Inferred s4U-induced muatation rate per read

  for(j in 1:nMT){
    for(k in 1:nrep){
      log_lambda_n[j, k] = log_lambda_o[j, k] + TL_lambda_eff[j, k];
    }
  }

}

// The actual statistical model is defined in the model block
model {


  for(i in 1:nMT){
    for(k in 1:nrep){
      log_lambda_o[i,k] ~ normal(-3.5, 0.5);
      TL_lambda_eff[i,k] ~ lognormal(-1, 0.5);
    }
  }

  // For loop because can't totally vectorize in Stan
  for (i in 1:NF) {


    for(j in 1:nMT){

      mu_rep_logit_fn[i, j] ~ normal(0, 1);


    }
  }

  // Model fraction of reads that are new with mixture of Poisson distributions
  for (i in 1:NE) {
    target += num_obs[i] * log_mix(TP[i] * frac_new[FE[i], MT[i], R[i]],
                                   poisson_log_lpmf(num_mut[i] | (U_cont[i] + log_lambda_n[MT[i], R[i]])),
                                   poisson_log_lpmf(num_mut[i] | (U_cont[i] + log_lambda_o[MT[i], R[i]])));
  } //reads
}

// Things you can estimate once you have finished estimating parameters and transformed parameters
generated quantities {

  vector[nMT] kd[NF]; // decay rate
  vector[nMT-1] L2FC_kd[NF]; // change in decay rate

  for (i in 1:NF) {
    for (j in 1:nMT) {
      kd[i, j] = -log(1 - mean(frac_new[i, j]))/tl[j];
    } // treatment type
    for(k in 1:nMT-1){
      L2FC_kd[i,k] = log2(kd[i,k+1]/kd[i,1]);
    }
  }

}



