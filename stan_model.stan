data{
    // data things
    int<lower=1> n_subjects; //Number of subjects
    int<lower=1> n_trials; //Number of trials per stimulus
    int<lower=0> hits[n_subjects];
    int<lower=0> fa[n_subjects];
    
    // individual d' is normal distributed
    // individual c is normal distributed
    // hyperprior for d' mean is cauchy
    real group_dprime_mu_hyper_loc;
    real<lower=0> group_dprime_mu_hyper_scale;
    // hyperprior for d' sd is half-cauchy
    real<lower=0> group_dprime_sd_hyper_loc;
    real<lower=0> group_dprime_sd_hyper_scale;
    
    // hyperprior for c mean is cauchy
    real group_bias_mu_hyper_loc;
    real<lower=0> group_bias_mu_hyper_scale;
    // hyperprior for c sd is half cauchy
    real<lower=0> group_bias_sd_hyper_loc;
    real<lower=0> group_bias_sd_hyper_scale;
}

parameters{
    //group dprime comes from a normal distribution
    real group_dprime_mean;
    real<lower=0> group_dprime_sd;
    //group bias comes from a normal distribution
    real group_bias_mean;
    real<lower=0> group_bias_sd;
    
    // dprime for each subject
    real d_prime_central[n_subjects];
    // use a non-central parameterisation
    real bias_central[n_subjects];
   
}
transformed parameters{
    real bias[n_subjects];
    real d_prime[n_subjects];
    real prob_fa[ n_subjects];
    real prob_hit[n_subjects];
    
    // Transform the centered parameters
    for ( n in 1:n_subjects ){
        bias[n] =  group_bias_mean + bias_central[n] * group_bias_sd;
        d_prime[n] = group_dprime_mean + d_prime_central[n] * group_dprime_sd;
        prob_fa[n] = Phi( (-2*bias[n] - d_prime[n])/2 );
        prob_hit[n] = Phi( (-2*bias[n] + d_prime[n] )/2 );
    }

}

model {
    // 4 hyper priors
    group_dprime_mean ~ cauchy( group_dprime_mu_hyper_loc, group_dprime_mu_hyper_scale );
    group_dprime_sd ~ cauchy( group_dprime_sd_hyper_loc, group_dprime_sd_hyper_scale );
    group_bias_mean ~ cauchy( group_bias_mu_hyper_loc, group_bias_mu_hyper_scale );
    group_bias_sd ~ cauchy( group_bias_sd_hyper_loc, group_bias_sd_hyper_scale );
    // 2 priors
    // use a centered parameterisation to reduce transition problems
    // See transformed parameters
    bias_central ~ normal( 0, 1 );
    d_prime_central ~ normal( 0, 1 );

    hits ~ binomial( n_trials, prob_hit );
    fa ~ binomial( n_trials, prob_fa );
}
