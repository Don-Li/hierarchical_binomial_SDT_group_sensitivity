make_data = function( n_subjects, n_signal_trials, n_noise_trials, group_dprime_mu, group_dprime_sd, group_c_mu, group_c_sd ){

    subject_dprime = rnorm( n_subjects, group_dprime_mu, group_dprime_sd )
    subject_c = rnorm( n_subjects, group_c_mu, group_c_sd )
    
    fa_rates = pnorm( (-2*subject_c - subject_dprime)/2 )
    hit_rates = pnorm( (-2*subject_c + subject_dprime)/2 )

    hits = rbinom( n_subjects, n_signal_trials, hit_rates )
    fas = rbinom( n_subjects, n_noise_trials, fa_rates )
    
    true_params = list(
        n_subjects = n_subjects,
        n_signal_trials = n_signal_trials, n_noise_trials = n_noise_trials,
        subject = rep(1:n_subjects),
        group_dprime_mu = group_dprime_mu, group_dprime_sd = group_dprime_sd,
        group_c_mu = group_c_mu, group_c_sd = group_c_sd,
        subject_dprime = subject_dprime, subject_c = subject_c,
        hit_rates = hit_rates, fa_rates = fa_rates,
        hits = hits, fas = fas
    )
    
    dataset = data.frame(
        subject = 1:n_subjects,
        n_signal_trials = n_signal_trials, n_noise_trials = n_noise_trials,
        hit = hits, fa = fas
    )
    list( true_params = true_params, dataset = dataset )
}

restructure_import_df = function(import_dataframe){
    if ( class(import_dataframe) == "error" ){
        if( grepl( "rows are empty", import_dataframe ) ){
            error_message = "Error: It appears that the data file is empty."
        }
        return( list( error = TRUE, error_message = error_message, df = NULL ) )
    }
    if ( ncol( import_dataframe ) != 5 ){
        error_message = "Error: The data file must have five columns. See the import help for more details."
        return( list( error = TRUE, error_message = error_message, df = NULL ) )
    }
    
    df_names = names(import_dataframe)
    subject_col = grep( "ubject", df_names )
    n_signal_trials_col = grep( "signal_trials", df_names )
    n_noise_trials = grep( "noise_trials", df_names )
    hit_col = grep( "hit", df_names )
    fa_col = grep( "fa", df_names )
    
    match_lengths = lengths(list(subject_col, n_signal_trials_col,n_noise_trials, hit_col, fa_col))
    names(match_lengths) = c("Subject", "Signal trials", "Noise trials", "Hit", "FA")
    if ( any(match_lengths>1) ){
        error_message = paste0( "Error: Multiple matching names for the ",
            paste( names(match_lengths)[match_lengths>1], collapse = ", " ),
            " columns. See import help for more details." )
        return( list( error = TRUE, error_message = error_message, df = NULL ) )
    }
    
    import_dataframe = import_dataframe[ ,c(subject_col, n_signal_trials_col, n_noise_trials, hit_col, fa_col) ]
    names(import_dataframe) = c("subject", "n_signal_trials", "n_noise_trials", "hit", "fa" )
    list( error = FALSE, error_message = "", df = import_dataframe )
}

plot_hyperpriors = function( 
    dprime_mu_hyperloc, dprime_mu_hyperscale, dprime_sd_hyperloc, dprime_sd_hyperscale,
    c_mu_hyperloc, c_mu_hyperscale, c_sd_hyperloc, c_sd_hyperscale, dataset ){
    
    dprime_mu_x_space = seq( -10, 10, by = 0.01 )
    dprime_mu_density = dcauchy( dprime_mu_x_space, dprime_mu_hyperloc, dprime_mu_hyperscale )
    
    dprime_sd_x_space = seq( 0, 10, by = 0.01 )
    dprime_sd_density = dcauchy( dprime_sd_x_space, dprime_sd_hyperloc, dprime_sd_hyperscale )
    zero_density = pcauchy( 0, dprime_sd_hyperloc, dprime_sd_hyperscale )
    dprime_sd_density = dprime_sd_density/zero_density
    
    c_mu_x_space = seq( -10, 10, by = 0.01 )
    c_mu_density = dcauchy( c_mu_x_space, c_mu_hyperloc, c_mu_hyperscale )
    
    c_sd_x_space = seq( 0, 10, by = 0.01 )
    c_sd_density = dcauchy( c_sd_x_space, c_sd_hyperloc, c_sd_hyperscale )
    zero_density = pcauchy( 0, c_sd_hyperloc, c_sd_hyperscale )
    c_sd_density = c_sd_density/zero_density
    
    point_estimates = c(NULL,NULL)
    if ( !is.null(dataset) ){
        point_estimates = compute_avg_dprime_estimates(dataset)
    }
    dprime_mu_vline = c(dprime_mu_hyperloc, point_estimates$avg_dprime )
    pt_estimate_cols = rgb(c(0,0), c(0,0), c(0,1), c(1,0.5))
    
    c_mu_vline = c(c_mu_hyperloc, point_estimates$avg_c )
    dprime_sd_vline = c(dprime_sd_hyperloc, point_estimates$sd_dprime)
    c_sd_vline = c(c_sd_hyperloc, point_estimates$sd_c )

    par(mfrow = c(2,2))
    plot( dprime_mu_x_space, dprime_mu_density, type = "l",
        main = "Group d' mu hyperprior", ylab = "Density", xlab = "Group d' mu")
    abline( v = dprime_mu_vline, col = pt_estimate_cols )
    
    plot( dprime_sd_x_space, dprime_sd_density, type = "l",
        main = "Group d' sd hyperprior", ylab = "Density", xlab = "Group d' sd")
    abline( v = dprime_sd_vline, col = pt_estimate_cols )
    
    plot( c_mu_x_space, c_mu_density, type = "l",
        main = "Group c mu hyperprior", ylab = "Density", xlab = "Group c mu")
    abline( v = c_mu_vline, col = pt_estimate_cols )
    
    plot( c_sd_x_space, c_sd_density, type = "l",
        main = "Group c sd hyperprior", ylab = "Density", xlab = "Group c sd")
    abline( v = c_sd_vline, col = pt_estimate_cols )
}

compute_avg_dprime_estimates = function( dataset ){
    avg_hit_q = qnorm( (dataset$hit + 0.5)/(dataset$n_signal_trials+1) )
    avg_fa_q = qnorm( (dataset$fa + 0.5)/(dataset$n_noise_trials+1) )
    averaged_dprime = mean( avg_hit_q - avg_fa_q )
    averaged_c = mean( 0.5*(avg_hit_q+avg_fa_q) )
    sd_dprime = sd(avg_hit_q - avg_fa_q)
    sd_c = sd(-0.5*(avg_hit_q+avg_fa_q))
    
    list( avg_dprime = averaged_dprime, avg_c = averaged_c, sd_dprime = sd_dprime,
        sd_c = sd_c)
}
# 
run_mcmc = function(dataset, dprime_mu_hyperloc, dprime_mu_hyperscale, dprime_sd_hyperloc, dprime_sd_hyperscale,
    c_mu_hyperloc, c_mu_hyperscale, c_sd_hyperloc, c_sd_hyperscale, adapt_delta, max_treedepth, chains,
    iter, warmup_iter){
    
    data_list = list(
        n_subjects = max(dataset$subject),
        n_signal_trials = dataset$n_signal_trials, n_noise_trials = dataset$n_noise_trials,
        hits = dataset$hit, fa = dataset$fa,
        group_dprime_mu_hyper_loc = dprime_mu_hyperloc, group_dprime_mu_hyper_scale = dprime_mu_hyperscale,
        group_dprime_sd_hyper_loc = dprime_sd_hyperloc, group_dprime_sd_hyper_scale = dprime_sd_hyperscale,
        group_bias_mu_hyper_loc = c_mu_hyperloc, group_bias_mu_hyper_scale = c_mu_hyperscale,
        group_bias_sd_hyper_loc = c_sd_hyperloc, group_bias_sd_hyper_scale = c_sd_hyperscale
    )
    control_list= list(
        adapt_delta = adapt_delta, max_treedepth = max_treedepth
    )
    
    hierarhical_model = sampling( model_obj, data = data_list, control = control_list,
        iter = iter, warmup = warmup_iter, refresh = 1000, chains = chains )
    cat("Sampling finished\n")
    hierarhical_model
}

plot_group_posteriors = function(model_fits, prior_things){
    extract_ = extract(model_fits)
    
    dprime_mu_x_space = seq( min(extract_$group_dprime_mean), max(extract_$group_dprime_mean), by = 0.01 )
    dprime_mu_density = dcauchy( dprime_mu_x_space, prior_things$dprime_mu_hyperloc, prior_things$dprime_mu_hyperscale )
    
    dprime_sd_x_space = seq( min( extract_$group_dprime_sd), max( extract_$group_dprime_sd), by = 0.01 )
    dprime_sd_density = dcauchy( dprime_sd_x_space, prior_things$dprime_sd_hyperloc, prior_things$dprime_sd_hyperscale )
    zero_density = pcauchy( 0, prior_things$dprime_sd_hyperloc, prior_things$dprime_sd_hyperscale )
    dprime_sd_density = dprime_sd_density/zero_density
    
    c_mu_x_space = seq( min(extract_$group_bias_mean), max(extract_$group_bias_mean), by = 0.01 )
    c_mu_density = dcauchy( c_mu_x_space, prior_things$c_mu_hyperloc, prior_things$c_mu_hyperscale )
    
    c_sd_x_space = seq( min(extract_$group_bias_sd), max(extract_$group_bias_sd), by = 0.01 )
    c_sd_density = dcauchy( c_sd_x_space, prior_things$c_sd_hyperloc, prior_things$c_sd_hyperscale )
    zero_density = pcauchy( 0, prior_things$c_sd_hyperloc, prior_things$c_sd_hyperscale )
    c_sd_density = c_sd_density/zero_density
    
    dataset = prior_things$dataset
    dprimes = qnorm((dataset$hit + 0.5)/(dataset$n_signal_trials+1)) - qnorm((dataset$fa + 0.5)/(dataset$n_noise_trials+1))
    avg_dprime = mean(dprimes)
    
    criterion = -0.5*(qnorm((dataset$hit + 0.5)/(dataset$n_signal_trials+1)) + qnorm((dataset$fa + 0.5)/(dataset$n_noise_trials+1)))
    avg_criterion = mean(criterion)
    
    pooed_dprime = qnorm((sum(dataset$hit) + 0.5)/(sum(dataset$n_signal_trials)+1)) - 
        qnorm((sum(dataset$fa) + 0.5)/(sum(dataset$n_noise_trials)+1))
    pooled_criterion = -0.5* (qnorm((sum(dataset$hit) + 0.5)/(sum(dataset$n_signal_trials)+1)) +
        qnorm((sum(dataset$fa) + 0.5)/(sum(dataset$n_noise_trials)+1)))
    
    col_ = rgb( c(0,0,1), c(0,0,0), c(0,1,0), c(0.5,0.5,0.5) )
    
    par( mfrow = c(2,2) )
    hist( extract_$group_dprime_mean, xlab = "Group d' mean", main = NULL, prob = T )
    lines( dprime_mu_x_space, dprime_mu_density )
    abline( v = c(mean(extract_$group_dprime_mean), avg_dprime, pooed_dprime), col = col_ )
    
    hist( extract_$group_dprime_sd, xlab = "Group d' sd", main = NULL, prob = T )
    lines( dprime_sd_x_space, dprime_sd_density )
    
    hist( extract_$group_bias_mean, xlab = "Group c mean", main = NULL, prob = T )
    lines( c_mu_x_space, c_mu_density )
    abline( v = c(mean(extract_$group_bias_mean), avg_criterion, pooled_criterion), col = col_ )
    
    hist( extract_$group_bias_sd, xlab = "Group c sd", main = NULL, prob = T )
    lines( c_sd_x_space, c_sd_density )
}

text_group_posteriors = function(model_fits, prior_things){
    extract_ = extract(model_fits)
    group_params = c("group_dprime_mean", "group_dprime_sd", "group_bias_mean", "group_bias_sd" )
    marginal_posterior_summaries = lapply( group_params, function(x){
        summary( extract_[[x]] )
    } )
    names(marginal_posterior_summaries) = c("Group d' mu", "Group d' SD", "Group c mu", "Group c SD")
    
    dataset = prior_things$dataset
    dprimes = qnorm((dataset$hit + 0.5)/(dataset$n_signal_trials+1)) - qnorm((dataset$fa + 0.5)/(dataset$n_noise_trials+1))
    avg_dprime = mean(dprimes)
    
    criterion = -0.5*(qnorm((dataset$hit + 0.5)/(dataset$n_signal_trials+1)) + qnorm((dataset$fa + 0.5)/(dataset$n_noise_trials+1)))
    avg_criterion = mean(criterion)
    
    pooed_dprime = qnorm((sum(dataset$hit) + 0.5)/(sum(dataset$n_signal_trials)+1)) - 
        qnorm((sum(dataset$fa) + 0.5)/(sum(dataset$n_noise_trials)+1))
    pooled_criterion = -0.5* (qnorm((sum(dataset$hit) + 0.5)/(sum(dataset$n_signal_trials)+1)) +
        qnorm((sum(dataset$fa) + 0.5)/(sum(dataset$n_noise_trials)+1)))

    
    marginal_posterior_summaries[["Averaged estimate"]] = avg_dprime
    marginal_posterior_summaries[["Pooled dprime"]] = pooed_dprime
    
    marginal_posterior_summaries
}

plot_individual_posteriors_dprime = function(model_fits, prior_things){
    extract_ = extract(model_fits)
    n_subjects = ncol( extract_$d_prime )
    dataset = prior_things$dataset
    dprimes = qnorm((dataset$hit + 0.5)/(dataset$n_signal_trials+1)) - qnorm((dataset$fa + 0.5)/(dataset$n_noise_trials+1))
    criterion = -0.5*(qnorm((dataset$hit + 0.5)/(dataset$n_signal_trials+1)) + qnorm((dataset$fa + 0.5)/(dataset$n_noise_trials+1)))

    col_ = rgb( c(0,0,1), c(0,0,0), c(0,1,0), c(0.5,0.5,0.5) )
    xlim_range = c( min(extract_$d_prime)-0.5, max(extract_$d_prime)+0.5)
    
    par( mfrow = plot_grid_dims(n_subjects) )
    for ( s in 1:n_subjects ){
        hist( extract_$d_prime[,s], xlab = "d'", main = s, prob = T, xlim = xlim_range )
        abline( v = c(mean(extract_$d_prime[,s]), dprimes[s] ), col = col_ )
    }
}

text_individual_posteriors_dprime = function( model_fits, prior_things ){
    extract_ = extract(model_fits)
    dataset = prior_things$dataset

    dprimes = extract_$d_prime
    marginal_posterior_summaries = lapply( 1:ncol(dprimes), function(x){
        summary( dprimes[,x] )
    } )
    names(marginal_posterior_summaries) = paste0("Subject ", 1:ncol(dprimes))
    marginal_posterior_summaries
    
    dprimes = qnorm((dataset$hit + 0.5)/(dataset$n_signal_trials+1)) - qnorm((dataset$fa + 0.5)/(dataset$n_noise_trials+1))
    names(dprimes) = names(marginal_posterior_summaries)
    marginal_posterior_summaries[["Point estimate"]] = dprimes
    marginal_posterior_summaries
}

plot_grid_dims = function(x){
    # possible_rows = x/1:x
    # possible_cols = 1:x
    # possible_cols = possible_cols[ possible_rows %% 1 == 0 ]
    # possible_rows = possible_rows[ possible_rows %% 1 == 0 ]
    # 
    # get_this_one = which.max( 1/(possible_rows-possible_cols) )
    # 
    # rows = possible_rows[get_this_one]
    # cols = possible_cols[get_this_one]
    # c(rows,cols)
    cols = 3
    rows = ceiling(x/3)
    c(rows,cols)
}

null_plot = function(message){
    plot(0,0,xlab = "", ylab = "", xaxt = "n", yaxt = "n", pch = NA)
    text( 0, 0, "Posterior results will be shown here." )
}

plot_individual_posteriors_c = function(model_fits, prior_things){
    extract_ = extract(model_fits)
    n_subjects = ncol( extract_$bias )
    dataset = prior_things$dataset
    criterion = -0.5*(qnorm((dataset$hit + 0.5)/(dataset$n_signal_trials+1)) + qnorm((dataset$fa + 0.5)/(dataset$n_noise_trials+1)))

    col_ = rgb( c(0,0,1), c(0,0,0), c(0,1,0), c(0.5,0.5,0.5) )
    
    xlim_range = c( min(extract_$bias)-0.5, max(extract_$bias)+0.5 )
    
    par( mfrow = plot_grid_dims(n_subjects) )
    for ( s in 1:n_subjects ){
        hist( extract_$bias[,s], xlab = "c", main = s, prob = T, xlim = xlim_range )
        abline( v = c(mean(extract_$bias[,s]), criterion[s] ), col = col_ )
    }
}

text_individual_posteriors_c = function( model_fits, prior_things ){
    extract_ = extract(model_fits)
    dataset = prior_things$dataset

    criterion = extract_$bias
    marginal_posterior_summaries = lapply( 1:ncol(criterion), function(x){
        summary( criterion[,x] )
    } )
    names(marginal_posterior_summaries) = paste0("Subject ", 1:ncol(criterion))
    marginal_posterior_summaries
    
    criterion = -0.5*(qnorm((dataset$hit + 0.5)/(dataset$n_signal_trials+1)) + qnorm((dataset$fa + 0.5)/(dataset$n_noise_trials+1)))
    names(criterion) = names(marginal_posterior_summaries)
    marginal_posterior_summaries[["Point estimate"]] = criterion
    marginal_posterior_summaries
}



