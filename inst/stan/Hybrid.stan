//STAN MODEL DETAILS:
//Takes as input fraction news estimated for
//each replicate, as well as the standard error
//of the fn estimates. Uses a model of the error
//in fn inference to estimate the true fn as well
//as to pool replicate information to infer
//the mean fn and L2FC(kdeg)
data {
  int NE; //num entries in data frame.
  int NF; //num features
  int MT[NE]; // mutant type for each entry (WT = 1, treatment A = 2, treatment B = 3, treatment C = 4)
  int nMT; // number of different mutation types
  int FE[NE]; //  feature  for each entry
  int R[NE]; //Replicate ID
  int nrep;
  real tl;
  real logit_fn_rep[NE]; //Replicate fn estimate
  real fn_se[NE] ; //Standard error of replicate fn estimate
  real Avg_Reads[NF, nMT]; // Average read counts in transcript i and sample j (avg. across replicates)
}

parameters {
  vector[NF] alpha;
  real mu_fn;
  real log_sig_fn;
  vector[nMT-1] log_sig_e;
  vector[nMT-1] mu_e;
  vector[nMT-1] z_e [NF];
  //real mu_rep_logit_fn[NF, nMT, nrep]; // Inferred fraction new of obs reads on native scale.
  real z_fn[NF, nMT, nrep];
  vector<lower=0>[nMT] a;
  vector[nMT] b;
  vector<lower=0>[nMT] sd_rep;
  vector[nMT] z_rep [NF];
  //vector<lower=0>[nMT] sd_r_mu[NF];
}

transformed parameters {
  real mu_rep_logit_fn[NF, nMT, nrep]; // Inferred fraction new of obs reads on native scale.
  vector[nMT-1] eff[NF]; //  Parameter for fraction new of observed reads
  real<lower=0> sig_fn = exp(log_sig_fn);
  vector<lower=0>[nMT-1] sig_e = exp(log_sig_e);
  vector<lower=0>[nMT] sd_r_mu[NF];
  //vector[nMT] sd_mean[NF];

    // Left here if needed for non-centered parameterization:
    for (i in 1:NF) {
       for(j in 1:nMT){
          //sd_mean[i,j] = -a[j]*Avg_Reads[i,j] + b[j];
          sd_r_mu[i,j] = exp(-a[j]*Avg_Reads[i,j] + b[j] + sd_rep[j]*z_rep[i,j]);

          if(j > 1){
             eff[i, j-1] = mu_e[j-1] + z_e[i,j-1]*sig_e[j-1];
          }
          for(k in 1:nrep){
            if(j == 1){
               mu_rep_logit_fn[i, j, k] = alpha[i] + z_fn[i,j,k]*sd_r_mu[i,j];
            }else{
               mu_rep_logit_fn[i, j, k] = (alpha[i] + eff[i, j - 1]) + z_fn[i,j,k]*sd_r_mu[i,j];
            }

          }
       }
    }
}

model {
  // Priors:
  mu_fn ~ normal(0, 1.5);
  log_sig_fn ~ normal(-1, 1);
  alpha ~ normal(mu_fn, sig_fn);

  log_sig_e ~ normal(-1, 1);
  mu_e ~ normal(0, 1);

  a ~ normal(0.6, 0.25);
  b ~ normal(0, 0.5);
  sd_rep ~ lognormal(-1,1);

  for (i in 1:NF) {
    //eff[i] ~ normal(mu_e, sig_e);
    z_e[i] ~ normal(0, 1);
    z_rep[i] ~ normal(0,1);
    for(j in 1:nMT){
     //sd_r_mu[i] ~ lognormal(sd_mean[i,j], sd_rep[j]);

     for(k in 1:nrep){
      z_fn[i, j, k] ~ normal(0, 1);
     }
    }
  }

    // Model fraction of reads that are new.
    for (i in 1:NE) {
        logit_fn_rep[i] ~ normal( mu_rep_logit_fn[FE[i], MT[i], R[i]], fn_se[i]);

        // if(MT[i] > 1){
        //   logit_fn_rep[i] ~ normal(alpha[FE[i]] + eff[FE[i], MT[i] - 1], sd_r_mu[FE[i], MT[i]] );
        // }else{
        //   logit_fn_rep[i] ~ normal(alpha[FE[i]], sd_r_mu[FE[i], MT[i]] );
        // }
    }
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
   for(j in 1:nMT-1){
     L2FC_kd[i,j] = log2(kd[i,j+1]/kd[i,1]);
   }
 }

 }

