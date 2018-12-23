// -*- mode: C -*-

data {

  int<lower=1> num_experiments;
  int<lower=1> num_peptides;
  int<lower=1> num_total_observations;
  int<lower=1> num_pep_exp_pairs;

  int<lower=1> muij_map[num_total_observations];
  int<lower=1> muij_to_exp[num_pep_exp_pairs];
  int<lower=1> muij_to_pep[num_pep_exp_pairs];
  
  int<lower=1, upper=num_experiments> experiment_id[num_total_observations];
  int<lower=1, upper=num_peptides> peptide_id[num_total_observations];
  real<lower=0> retention_times[num_total_observations];
}
parameters {

  // canonical retention time for given peptide
  real<lower=0> mu[num_peptides];


  real<lower=0> beta_1[num_experiments];
  real<lower=0> beta_2[num_experiments];
  real<lower= -1.0*min(beta_1)*min(mu)> beta_0[num_experiments];
  real<lower=0> split_point[num_experiments];

  
  real<lower=0> sigma_slope_global;
  real<lower=0> sigma_slope[num_experiments];
  real<lower=0> sigma_intercept[num_experiments];

}
transformed parameters {

  real<lower=0> sigma_ij[num_pep_exp_pairs];
  
  // alignment is based on a segemented linear regression

  // muij is the mean retention time for peptide i in experiment j 
  real<lower=0> muij[num_pep_exp_pairs];
  for (i  in 1:num_pep_exp_pairs) {
    if(mu[muij_to_pep[i]] < split_point[muij_to_exp[i]]) {
      muij[i] = beta_0[muij_to_exp[i]] + 
        beta_1[muij_to_exp[i]] * mu[muij_to_pep[i]];
    } else if( mu[muij_to_pep[i]] >=  split_point[muij_to_exp[i]] ) {
      muij[i] = beta_0[muij_to_exp[i]] + 
        beta_1[muij_to_exp[i]] * split_point[muij_to_exp[i]] + 
        beta_2[muij_to_exp[i]] * (mu[muij_to_pep[i]] - split_point[muij_to_exp[i]]);
      //muij[i] = beta_0[muij_to_exp[i]] + beta_1[muij_to_exp[i]] * mu[muij_to_pep[i]];
    
    }
  }

  // standard deviation grows linearly with time
  for(i in 1:num_pep_exp_pairs) { 
    sigma_ij[i] = sigma_intercept[muij_to_exp[i]] + 
      sigma_slope[muij_to_exp[i]] / 100 * mu[muij_to_pep[i]];
  }

}
model {
  mu ~ lognormal(mean(log(retention_times)), sd(log(retention_times)));

  sigma_slope_global ~ lognormal(0.1, 0.5);
  //sigma_slope ~ lognormal(log(sigma_slope_global), 0.5);
  sigma_slope ~ lognormal(log(sigma_slope_global), 1);
  sigma_intercept ~ lognormal(0, 2);
  
  beta_0 ~ normal(0, 1);
  beta_1 ~ lognormal(0, 0.5);
  beta_2 ~ lognormal(0, 0.5);

  for(i in 1:num_experiments){
    split_point ~ uniform(min(muij), max(muij));
  }
  
  for (i in 1:num_total_observations) {
    retention_times[i] ~ double_exponential(muij[muij_map[i]], sigma_ij[peptide_id[i]]);
  }


}
