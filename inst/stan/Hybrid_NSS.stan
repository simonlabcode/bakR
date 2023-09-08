// Used to be called Hybrid_Model.stan
data {
  int NE; //num entries in data frame.
  int NF; //num features
  array[NE] int MT; // mutant type for each entry (WT = 1, treatment A = 2, treatment B = 3, treatment C = 4)
  int nMT; // number of different mutation types
  array[NE] int FE; //  feature  for each entry
  array[NE] int R; //Replicate ID
  int nrep;
  array[nMT] real tl;
  array[NE] real logit_fn_rep; //Replicate fn estimate
  array[NE] real fn_se ; //Standard error of replicate fn estimate
  array[NF, nMT] real Avg_Reads; // Average read counts in transcript i and sample j (avg. across replicates)
  int<lower=0, upper=1> Chase;
}

parameters {
  array[NF] vector[nMT] alpha;
  vector[nMT] mu_fn;
  vector[nMT] log_sig_fn;
  //vector[nMT-1] z_e [NF];
  //real mu_rep_logit_fn[NF, nMT, nrep]; // Inferred fraction new of obs reads on native scale.
  array[NF, nMT, nrep] real z_fn;
  vector<lower=0>[nMT] a;
  vector[nMT] b;
  vector<lower=0>[nMT] sd_rep;
  array[NF] vector[nMT] z_rep;
  //vector<lower=0>[nMT] sd_r_mu[NF];
}

transformed parameters {
  array[NF, nMT, nrep] real mu_rep_logit_fn; // Inferred fraction new of obs reads on native scale.
  //vector[nMT-1] eff[NF]; //  Parameter for fraction new of observed reads
  vector<lower=0>[nMT] sig_fn = exp(log_sig_fn);
  array[NF] vector<lower=0>[nMT] sd_r_mu;
  //vector[nMT] sd_mean[NF];

    // Left here if needed for non-centered parameterization:
    for (i in 1:NF) {
       for(j in 1:nMT){
          //sd_mean[i,j] = -a[j]*Avg_Reads[i,j] + b[j];
          sd_r_mu[i,j] = exp(-a[j]*Avg_Reads[i,j] + b[j] + sd_rep[j]*z_rep[i,j]);

          // if(j > 1){
          //    eff[i, j-1] = mu_e[j-1] + z_e[i,j-1]*sig_e[j-1];
          // }
          for(k in 1:nrep){
              mu_rep_logit_fn[i, j, k] = alpha[i,j] + z_fn[i,j,k]*sd_r_mu[i,j];


          }
       }
    }
}

model {
  // Priors:
  mu_fn ~ normal(0, 1.25);
  log_sig_fn ~ normal(-1, 0.5);



  a ~ normal(0.3, 0.2);
  b ~ normal(-1.5, 0.35);
  sd_rep ~ lognormal(-2, 0.25);

  for (i in 1:NF) {
    //eff[i] ~ normal(mu_e, sig_e);
    //z_e[i] ~ normal(0, 1);
    z_rep[i] ~ normal(0,1);
    for(j in 1:nMT){
      alpha[i,j] ~ normal(mu_fn[j], sig_fn[j]);
     //sd_r_mu[i,j] ~ lognormal(sd_mean[i,j], sd_rep[j]);


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

   array[NF] vector[nMT] kd; // decay rate
   array[NF] vector[nMT-1] L2FC_kd; // change in decay rate
   array[NF] vector[nMT] log_kd; // log(kdeg)

  for (i in 1:NF) {
    for (j in 1:nMT) {

      if(Chase == 1){

        kd[i,j] = -log(inv_logit(alpha[i,j]))/tl[j];  // divide by time if not 1

      }else{

        kd[i,j] = -log(1 - inv_logit(alpha[i,j]))/tl[j];  // divide by time if not 1

      }
      
      log_kd[i,j] = log(kd[i,j]);


    } // treatment type
    for(k in 1:nMT-1){
      if(Chase == 1){

        L2FC_kd[i,k] = -alpha[i,k+1] + alpha[i,1];

      }else{

        L2FC_kd[i,k] = alpha[i,k+1] - alpha[i,1];

      }
    }
  }

}
