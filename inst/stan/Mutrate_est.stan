// Completely centered parameterization with no hierarchy
// Used to be called Mut_rate_Model.stan

// Data block corresponds to what gets passed to Stan in data_list
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
  int nrep;
  array[nMT] real tl;
  real nU;
}

// Parameter block tells model all the things it will estimate
parameters {
  array[nMT] vector[nrep] log_lambda_o; // Background mutation rate
  array[nMT] vector<lower=0>[nrep] TL_lambda_eff; // Increase to mutation rate due to s4U
  array[NF, nMT, nrep] real mu_rep_logit_fn; // Inferred fraction new of obs reads on native scale.
}

// Transformed parameters are parameters that are defined in terms of data or other parameters
transformed parameters {
  array[NF, nMT, nrep] real frac_new = inv_logit(mu_rep_logit_fn); // fraction new
  array[nMT] vector[nrep] log_lambda_n; // Inferred s4U-induced muatation rate per read

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
      log_lambda_o[i,k] ~ normal(log(nU*0.001), 0.5);
      TL_lambda_eff[i,k] ~ lognormal(1.4, 0.15);
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

  array[NF] vector[nMT] kd; // decay rate
  array[NF] vector[nMT-1] L2FC_kd; // change in decay rate

  for (i in 1:NF) {
    for (j in 1:nMT) {
      kd[i, j] = -log(1 - mean(frac_new[i, j]))/tl[j];
    } // treatment type
    for(k in 1:nMT-1){
      L2FC_kd[i,k] = log2(kd[i,k+1]/kd[i,1]);
    }
  }

}

