library( data.table )
library( rstan )

# Setting the seed will set the rng for R functions but not for stan, but this should be okay
set.seed(1)

make_data = function( n_subjects, n_trials, group_dprime_mu, group_dprime_sd, group_bias_mu, group_bias_sd ){
    
    # Convert mean and sd to alpha and beta for the gamma distribution
    # These equations are obtained using the method of moments
    alpha = group_dprime_mu^2/group_dprime_sd^(2)
    beta = group_dprime_mu/group_dprime_sd^2
    
    subject_dprime = rgamma( n_subjects, alpha, beta )
    subject_bias = rnorm( n_subjects, group_bias_mu, group_bias_sd )
    
    fa_rates = pnorm( (-2*subject_bias - subject_dprime)/2 )
    hit_rates = pnorm( (-2*subject_bias + subject_dprime)/2 )

    hits = rbinom( n_subjects, n_trials/2, hit_rates )
    fas = rbinom( n_subjects, n_trials/2, fa_rates )
    
    as.data.table( list(
        n_subjects = n_subjects,
        n_trials = n_trials,
        subject = rep( 1:n_subjects ),
        group_dprime_mu = group_dprime_mu, group_dprime_sd = group_dprime_sd,
        group_bias_mu = group_bias_mu, group_bias_sd = group_bias_sd,
        subject_dprime = subject_dprime, subject_bias = subject_bias,
        hit_rates = hit_rates, fa_rates = fa_rates,
        hits = hits, fas = fas
    ) )
}

##### Make some data #####
true_group_dprime_mean = 1.5
true_group_dprime_sd = 1.5/5
true_group_c_sd = 0.2
true_group_c_mean = 0
number_of_trials = 32
subjects_pooled = 16
dataset = make_data( subjects_pooled, number_of_trials, true_group_dprime_mean, true_group_dprime_sd,
    true_group_c_mean, true_group_dprime_sd )
##########################

##### Set up the model #####
# Compile the model
hierarchical_model = stan_model( "stan_model.stan" )

# Calculate the d' for each subject
# the constants 0.5 and 1 are for the log-linear correction
# n_trials divided by 2 because it is the sum of stimulus-present and noise trials
dataset[ , dprime_log_linear := qnorm( (hits+0.5)/(n_trials/2+1) ) - qnorm( (fas+0.5)/(n_trials/2+1)) ]

# Averaged d'
averaged_dprime = dataset[ , mean(dprime_log_linear) ]
# Pooled d'
pooled_dprime = dataset[ , {
    qnorm( (sum(hits)+0.5)/(sum(n_trials)/2+1) ) - qnorm( (sum(fas)+0.5)/(sum(n_trials)/2+1))
    } ]

# Set up data and hyperpriors for the bayesian model
# All list elements are defined in the data block in stan_model.stan
# To fit any custom data, just substitute the n_subjects, n_trials, fa, and hits with your data.
data_list = list(
    n_subjects = dataset$n_subjects[1],
    n_trials = dataset$n_trials[1]/2,
    fa = dataset$fas,
    hits = dataset$hits,
    # The following list elements are prior parameters
    # They are outside in the R code so we can change them without recompiling the stan model
    group_dprime_mu_hyper_loc = 1,
    group_dprime_mu_hyper_scale = 1,
    group_dprime_sd_hyper_loc = 0.5,
    group_dprime_sd_hyper_scale = 2,
    group_bias_mu_hyper_loc = 0,
    group_bias_mu_hyper_scale = 1,
    group_bias_sd_hyper_loc = 0.5,
    group_bias_sd_hyper_scale = 2
)

# Run the MCMC to generate a posterior
# adapt_delta and max_treedepth chosen based on trial and error
posterior = sampling( hierarchical_model, data = data_list,
    control = list( adapt_delta = 0.995, max_treedepth = 15), iter = 5000, refresh = 1000, warmup = 2000 )

# We have a lot of parameters
# 4 hyperpriors and n_subjects*2 prior parameters, so it is quite difficult to do comprehensive 
# diganostics for all of them
# We also have many transformed parameters because we used a centralised parameterisation, see the 
# stan model code for details.
# We will look at the hyperprior parameters for the diagnostics
traceplot(posterior, c("group_dprime_mean", "group_dprime_sd", "group_bias_mean", "group_bias_sd"))
# We want the trace plots to look like white noise and the trace plots for each chain to be
# indistinguishable from one another

# Diagnostics seem okay, so we can look at the marginal posteriors
hypers = extract( posterior, c("group_dprime_mean", "group_dprime_sd", "group_bias_mean", "group_bias_sd") )

dprime_range = seq(-1,4,by=0.01)
group_dprime_density = dgamma( dprime_range, true_group_dprime_mean^2/true_group_dprime_sd^2, 
    true_group_dprime_mean/true_group_dprime_sd^2 )

# Plot the true distribution of the group d'
plot( dprime_range, group_dprime_density, type = "l", ylab = "density", xlab = "dprime" )
# Plot the distribution of group d' from the hierarhical model
lines( dprime_range, dnorm( dprime_range, mean(hypers$group_dprime_mean), mean(hypers$group_dprime_sd) ),
    col = "red")
# Also plot a point estimate from the hierarhical model so that we can compare it with the
# pooled and averaged estimators
abline( v = mean(hypers$group_dprime_mean), col = "red" )
abline( v = averaged_dprime, col = "green" )
abline( v = pooled_dprime, col = "blue" )
legend( "topleft", legend = c("True", "Bayes model", "Averaged", "Pooled"), 
    col = c("black", "red", "green", "blue"), lty = 1 )

# Have a look at each subject
par( mfrow = c(4,4))
individual_dprime_posterior = extract( posterior, "d_prime" )[[1]]

# Compare the posterior individual d' with the sample estimated d'
for ( s_ in 1:subjects_pooled ){
    hist( individual_dprime_posterior[,s_], main = paste0("Subject ", s_),
        xlab = "individual d'", xlim = c(0, 4))
    abline( v = mean( individual_dprime_posterior[,s_] ) )
    abline( v = dataset[ s_, dprime_log_linear ], col = "red" )
    # plot the legend if you want, but the graphs are probably pretty small
    # if ( s_ == 16 ){
    #     legend( "topleft", legend = c("Sample", "Bayes model"), 
    #         col = c("black", "red"), lty = 1 )
    # }
}
# The posterior for the group d' is about 1.36
# The posteriors for each individual d' are pulled to around 1.36 also. This is a phenomenon
# called shrinkage, in that individual estimates are 'shrunk' towards the group mean.
# If you compare the individual sample estimates with the bayesian model, you will find that the
# sample estimates tend to be more extreme than the bayesian estimates on whatever side of 1.36



