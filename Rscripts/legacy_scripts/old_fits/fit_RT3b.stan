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
  real<lower=0> mean_log_rt;
  real<lower=0> sd_log_rt;
  real<lower=0, upper=1> pep[num_total_observations];

  real<lower=0> max_retention_time;
}
parameters {
  // canonical retention time for given peptide
  real<lower=0, upper=max_retention_time> mu[num_peptides];

  // segmented linear regression parameters for the 
  // retention time of each pep-exp pair
  //
  // beta_0 = y-intercept
  // beta_1 = slope of first segment
  // beta_2 = slope of second segment
  // split_point = x-value of the split b/n the two segments
  real<lower=0> beta_1[num_experiments];
  real<lower=0> beta_2[num_experiments];
  real<lower= -1.0 * min(beta_1) * min(mu) > beta_0[num_experiments];
  real<lower=0> split_point[num_experiments];

  // linear regression parameters for the standard deviation/spread 
  // for each pep-exp pair
  //
  // sigma_slope_global = average sigma slope for the set of pep-exp pairs.
  //                      "moving prior" for sigma_slopes of pairs moving forwards
  // sigma_slope = slope of line
  // sigma_intercept = y-intercept of line
  real<lower=0> sigma_slope_global;
  real<lower=0> sigma_slope[num_experiments];
  real<lower=0> sigma_intercept[num_experiments];

}
transformed parameters {
  // alignment is based on a segemented linear regression
  
  // sigma_ij is the spread of the laplace distribution of 
  // peptide i in experiment j
  real<lower=0> sigma_ij[num_pep_exp_pairs];
  
  // muij is the mean retention time for peptide i in experiment j 
  real<lower=0, upper=max_retention_time> muij[num_pep_exp_pairs];
  
  for (i in 1:num_pep_exp_pairs) {
    if(mu[muij_to_pep[i]] < split_point[muij_to_exp[i]]) {
      // first segment of the linear regression (before split_point)
      muij[i] = beta_0[muij_to_exp[i]] + 
        beta_1[muij_to_exp[i]] * mu[muij_to_pep[i]];
    } else if( mu[muij_to_pep[i]] >=  split_point[muij_to_exp[i]]) {
      // second segment of the linear regression (after split_point)
      muij[i] = beta_0[muij_to_exp[i]] + 
        beta_1[muij_to_exp[i]] * split_point[muij_to_exp[i]] + 
        beta_2[muij_to_exp[i]] * (mu[muij_to_pep[i]] - split_point[muij_to_exp[i]]);
    }
    
    // simple linear regression
    //muij[i] = beta_0[muij_to_exp[i]] + beta_1[muij_to_exp[i]] * mu[muij_to_pep[i]];

    // make sure that the retention time is within the physical limits
    if( muij[i] > max_retention_time ) {
      muij[i] = max_retention_time * 0.99;
    }
  }

  // standard deviation grows linearly with time - simple linear regression
  for(i in 1:num_pep_exp_pairs) { 
    sigma_ij[i] = sigma_intercept[muij_to_exp[i]] + 
      sigma_slope[muij_to_exp[i]] / 100 * mu[muij_to_pep[i]];
  }

}
model {
  //mu ~ lognormal(mean(log(retention_times)), sd(log(retention_times)));
  mu ~ lognormal(mean_log_rt, sd_log_rt);

  sigma_slope_global ~ lognormal(0.1, 0.5);
  //sigma_slope ~ lognormal(log(sigma_slope_global), 0.5);
  sigma_slope ~ lognormal(log(sigma_slope_global), 1);
  sigma_intercept ~ lognormal(0, 2);
  
  beta_0 ~ normal(0, 10);
  beta_1 ~ lognormal(0, 0.5);
  beta_2 ~ lognormal(0, 0.5);

  for(i in 1:num_experiments){
    split_point ~ uniform(0, max_retention_time);
  }
  
  for (i in 1:num_total_observations) {
    // mixture between laplace and a uniform prior (distribution over all peptides)
    real comp1 = log(1-pep[i]) + double_exponential_lpdf(retention_times[i] | muij[muij_map[i]], sigma_ij[peptide_id[i]]);
    // real comp2 = log(pep[i]) + uniform_lpdf(retention_times[i] | 0, max_retention_time);
    real comp2 = log(pep[i]) + lognormal_lpdf(retention_times[i] | mean_log_rt, sd_log_rt);

    real lse = log_sum_exp(comp1, comp2); 

    target += lse;

  }


}
