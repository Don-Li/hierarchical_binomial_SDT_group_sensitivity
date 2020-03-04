cat( "Loading libraries\n" )

if (!require(shiny)) install.packages('shiny')
library(shiny)
if (!require(rstan)) install.packages('rstan')
library(rstan)

cat( "Sourcing files\n" )

source( "app_make_data.R" )

cat( "Application starting up\n" )
cat( "Compiling stan model. This may take around a minute.\n")
model_obj = stan_model( "app_stan_model.stan" )
cat( "Stan compilation finished. Application starting.\n" )

session_storage = new.env()

##### Import panel #####
data_choice_panel = radioButtons( inputId = "data_or_simulation", label = "Data input", 
    choices = c("Real data", "Simulated data"), inline = TRUE )
import_file_title = tagList( "Select file", actionLink( "import_help", "(Help)" ) )
import_file_panel = fileInput( inputId = "file_import", label = import_file_title, multiple = FALSE )

data_input_panel = conditionalPanel(
    condition = "input.data_or_simulation == 'Real data'",
    import_file_panel
)

import_help_dialog = tagList(
    tags$strong("Help for importing files."),
    hr(),
    tags$strong("File types:"),
    tags$p("Only .csv files are accepted."),
    hr(),
    tags$strong("File structure:"),
    tags$p("If the data file contains no data, then you should receive an error."),
    tags$p("Your data file should contain five columns containing the data for the subject ID, number of signal trials,
        the number of noise trials, the number of hits, and the number of false alarms."),
    tags$p("The columns are identified based on regular expression (essentially a character search)."),
    tags$p("The subject ID column is searched for using 'ubject', so any variation on this (e.g., Subject, subject, subject_id)
        will tell the program that this column is the subject ID column" ),
    tags$p("The number of signal trials and number of noise trials are searched for using 'signal_trials' and 'noise_trials', respectively" ),
    tags$p("The hits and false alarms are searched for using 'hit' and 'fa'."),
    tags$p("If there are more than six columns and/or there are multiple columns that match each of these criteria, 
        then an error will be thrown."),
    tags$p("The columns will be rearranged in the order described above."),
    tags$p("If in doubt, you may simulate some data, save the dataset, and substitute your own data into that file." )
)
########################

##### Simulation parameters #####
sim_param_panel = conditionalPanel(
    condition = "input.data_or_simulation == 'Simulated data'",
    tags$strong("True group parameters"),
    actionLink( "sim_param_help", "(Help)" ),
    
        splitLayout(
            numericInput( inputId = "Group_dprime_mean", 
                label = div( style = "font-weight:normal", "Group d' mean" ), value = 1),
            numericInput( inputId = "Group_dprime_sd", 
                label = div( style = "font-weight:normal", "Group d' sd"), value = 0.2, min = 0)
        ),
        splitLayout(
            numericInput( inputId = "Group_c_mean", 
                label = div( style = "font-weight:normal", "Group c mean"), value = 0),
            numericInput( inputId = "Group_c_sd", 
                label = div( style = "font-weight:normal", "Group c sd"), value = 0.1, min = 0)
        ),
        splitLayout(
            numericInput( inputId = "n_noise_trials", 
                label = div( style = "font-weight:normal", "Noise trials per subject"), value = 8, min = 0 ),
            numericInput( inputId = "n_signal_trials", 
                label = div( style = "font-weight:normal", "Signal trials per subject"), value = 8, min = 0 )
        ),
        numericInput( inputId = "n_subjects", 
            label = div( style = "font-weight:normal", "Number of subjects"), value = 8, min = 0 ),
        splitLayout(
            actionButton( inputId = "simulate_data", 
                label = "Simulate data"),
            downloadButton( outputId = "download_data", 
                label = "Save data" )
        )
    
)
all_input_panel = wellPanel( data_choice_panel, 
    hr(), 
    # Following two are conditional on input.data_or_simulation
    data_input_panel,
    sim_param_panel
)

parameter_help_dialog = tagList(
    tags$p("Select the values for the true group parameters for a simulated dataset."),
    tags$p("For the simulation, individual values of d' and criterion (c) are sampled from hierarhical normal distributions."),
    tags$p("For simplicity, there will be an equal number of trials for each subject and half the trials will be stimulus-present trials." ),
    tags$p("If you choose to import your own data, there are no restrictions on the allocation of trials." )
)
##################################

##### Hyperparameters #####
hyperprior_panel = wellPanel(
    tags$strong("Hyperprior parameters"),
    actionLink( "hyper_prior_help", "(Help)" ),
    hr(),
    splitLayout(
        numericInput( inputId = "Group_dprime_mu_hyperloc", 
            label = div( style = "font-weight:normal", "Group d' mu hyperlocation"), value = 1),
        numericInput( inputId = "Group_dprime_mu_hyperscale", 
            label = div( style = "font-weight:normal", "Group d' mu hyperscale"), value = 1)
    ),
    splitLayout(
        numericInput( inputId = "Group_dprime_sd_hyperloc", 
            label = div( style = "font-weight:normal", "Group d' sd hyperlocation"), value = 0.5),
        numericInput( inputId = "Group_dprime_sd_hyperscale", 
            label = div( style = "font-weight:normal", "Group d' sd hyperscale"), value = 2)
    ),
    splitLayout(
        numericInput( inputId = "Group_c_mu_hyperloc", 
            label = div( style = "font-weight:normal", "Group c mu hyperlocation"), value = 0),
        numericInput( inputId = "Group_c_mu_hyperscale", 
            label = div( style = "font-weight:normal", "Group c mu hyperscale"), value = 1)
    ),
    splitLayout(
        numericInput( inputId = "Group_c_sd_hyperloc", 
            label = div( style = "font-weight:normal", "Group c sd hyperlocation"), value = 0.5),
        numericInput( inputId = "Group_c_sd_hyperscale", 
            label = div( style = "font-weight:normal", "Group c sd hyperscale"), value = 2)
    )
)

hyper_prior_help_dialog = tagList(
    tags$p("In our hierarhical SDT model, the d' and c values for each individual are sampled from group-level distributions."),
    tags$p("The group-level d' and c distributions are normal. We have prior distributions for each of these parameters."),
    tags$p("The group d' mu and group c mu are Cauchy distributed, with corresponding location and scale parameters."),
    tags$p("The group d' sd and group c sd are Half-Cauchy distributed with corresponding location and scale, parameterised as 
        if they were full Cauchy distributions")
)
##########################

##### Other options #####
options_panel = wellPanel(
    tags$strong("Options"),
    actionLink( "options_help", "(Help)" ),
    hr(),
    tags$p("MCMC options:"),
    splitLayout(
        numericInput( "mcmc_warmup_iter", "Warmup iterations", value = 1000, min = 0 ),
        numericInput( "mcmc_max_iter", "MCMC iterations", value = 5000, min = 0 )
    ),
    numericInput( "mcmc_chains", "Number of chains", value = 4, min = 1 ),
    splitLayout(
        numericInput( "mcmc_adapt_delta", "Adapt delta", value = 0.95, min = 0, max = 1 ),
        numericInput( "mcmc_max_treedepth", "Max treedepth", value = 15, min = 0, max = 30 )
    ),
    hr(),
    actionButton( "clear_dataset", "Clear data" ),
    actionButton( "run_mcmc", "Run MCMC" ),
    downloadButton( "save_Rdata", "Save .RData")
)

options_help_dialog = tagList(
    tags$strong("MCMC options help"),
    tags$p("To fit the model, we use Stan. Stan uses a variation of Markov Chain Monte Carlo called Hamiltonian Monte Carlo"),
    tags$p("Like most MCMC schemes, we need to specify some number of iterions for the Markov chain. 
        Usually longer chains are better because they give the chain a better chance of exploring more of the
        posterior distribution. But there are diminishing returns. For our SDT model, 4000-5000 seeems to be a good number."),
    tags$p("In addition to the maximum number of iterions, we also need to specify some number of 'warmup' iterations. 
        In other MCMC schemes, this is analogous to 'burn-in' or 'adaptation' iterations. For our SDT model, 1000 seems okay."),
    tags$p("For HMC, we also need to specify 'adapt delta' and 'max treedepth'"),
    tags$p("Adapt delta is somewhat like the resolution of the sampler. Occasionally, we will get warnings from the sampler about
        some number of divergent transitions. Setting 'adapt delta' closer to 1.0 will sometimes help with this problem. Another way
        to resolve divergent transitions is to parameterise the model so that transitions are easier. On inspection of the stan code for
        this model, you can see we have used non-central parameterisations for d' and c."),
    tags$p("However, sometimes divergent transitions cannot be helped. This is often due to a combination between the complexity of the
        space that the model has to sample (e.g., if you want to sample covariance matrices) and the data. For example, fitting our
        hierarhical model to a dataset with two subjects and 4 trials will often lead to divergent transitions that cannot
        be resolved."),
    tags$p("Max treedepth is related to the efficiency of the sampler. It is usually not as severe of a problem as divergent transitions."),
    tags$p("It is also usually a good idea to run multiple Markov chains instead of just one. In one sense, this allows us to check 
        the validity of the posterior. We want the chains to 'overlap' to convince us that the chains are adequately exploring the posterior. 
        More chains also give us more posterior samples. In general, we prefer a few long chains rather than many short chains."),
    tags$p("Refer to the Stan documentation if anything is confusing:", 
        tags$a( "https://mc-stan.org/users/documentation/", href = "https://mc-stan.org/users/documentation/" )
    )
)
##########################

##### MCMC results #####
mcmc_results_panel = wellPanel(
    tags$strong("MCMC results"),
    hr(),
    tags$strong("Group posterior"),
    actionLink("group_posterior_help", "(Help)"),
    actionLink("group_posterior_plot_popout", "(Popout)"),
    HTML('<button data-toggle="collapse" data-target="#group_posterior_plots" style = "float: right">-</button>'),
    tags$div( class = "collapse in", id = "group_posterior_plots",
        plotOutput("group_posterior_plots"),
        verbatimTextOutput("posterior_group_text")
    ),
    hr(),
    tags$strong("Individual posteriors for d'"),
    actionLink("individual_posterior_plots_dprime_popout", "(Popout)"),
    HTML('<button data-toggle="collapse" data-target="#individual_posterior_plots" style = "float: right">-</button>'),
    tags$div( class = "collapse in", id = "individual_posterior_plots",
        plotOutput("individual_posterior_plots_dprime"),
        verbatimTextOutput("posterior_individual_text_dprime")
    ),
    hr(),
    tags$strong("Individual posteriors for c"),
    actionLink("individual_posterior_plots_c_popout", "(Popout)"),
    HTML('<button data-toggle="collapse" data-target="#individual_posterior_plots_c" style = "float: right">-</button>'),
    tags$div( class = "collapse in", id = "individual_posterior_plots_c",
        plotOutput("individual_posterior_plots_c"),
        verbatimTextOutput("posterior_individual_text_c")
    )
)

group_posterior_help_dialog = tagList(
    tags$p("These four plots show the marginal posterior for the group parameters. The lines show the prior density."),
    tags$p("The black vertical line is the mean of the marginal posterior."),
    tags$p("The blue vertical line is the point estimate obtained by averaging (either d' or c)."),
    tags$p("The red vertical line is the point estimate obtained by pooling trials across subjects."),
    tags$p("For the averaged and pooled estimators, a loglinear correction (+0.5 to cells, +1 to trials) was used.")
)


########################

##### Data summary panel #####
data_summary_panel = wellPanel(
    tags$strong( "Data summaries" ),
    HTML('<button data-toggle="collapse" data-target="#data_summary_panel" style = "float: right">-</button>'),
    hr(),
    tags$div( class = "collapse in", id = "data_summary_panel",
        dataTableOutput("data_summary")
    )
)
##############################

##### Hyperprior #####
hyperprior_plot_panel = wellPanel(
    tags$strong( "Hyperprior visualisations" ),
    actionLink("hyperprior_plot_help", "(Help)"),
    actionLink("hyperprior_plot_popout", "(Popout)"),
    HTML('<button data-toggle="collapse" data-target="#hyperprior_plot_panel" style = "float: right">-</button>'),
    hr(),
    tags$div( class = "collapse in", id = "hyperprior_plot_panel",
        plotOutput("hyperprior_plots")
    )
)

hyperprior_plot_info = tagList(
    tags$p("In our hierarhical SDT model, individual values of d' and c come from normal distributions"),
    tags$p("The mu for the group d' and c have Cauchy hyperpriors, while the sd parameters have half-Cauchy hyperpriors"),
    tags$p("The parameters for these plots are controlled by the 'Hyperprior parameters' panel."),
    tags$p("These hyperprior densities are shown in the plots."),
    tags$p("The vertical black lines indicate the location/mode for these distributions"),
    tags$p("If there is data loaded into the application, there will be vertical blue lines on the plot"),
    tags$p("The blue vertical lines indicate sample values for these parameters. 
        The sample values are computed by the average or standard deviation of d' or c.
        A log-linear correction of 0.5 is used."),
    tags$p("The blue vertical line can be used as a guide to help you choose your hyperprior parameters.")
)
######################

main_title = "Bayesian hierarhical signal detection"

ui = fluidPage(
    titlePanel(main_title),
    hr(),
    fluidRow(
        column( 4, 
            all_input_panel,
            hyperprior_panel,
            options_panel
        ),
        column( 8, 
            data_summary_panel,
            hyperprior_plot_panel,
            mcmc_results_panel
        )
    )
)

datasetupload = reactiveVal(0)
posterior_results = reactiveVal(0)

server = function( input, output ){
    output$data_summary = renderDataTable( data.frame( "No data loaded" = "No data loaded" ) )
    
    ##### Help button observers #####
    observeEvent( input$sim_param_help,
        showModal(modalDialog(parameter_help_dialog))
    )
    observeEvent( input$hyper_prior_help,
        showModal(modalDialog(hyper_prior_help_dialog))
    )
    observeEvent( input$import_help,
        showModal(modalDialog(import_help_dialog))
    )
    observeEvent( input$hyperprior_plot_help,
        showModal(modalDialog(hyperprior_plot_info))
    )
    observeEvent( input$options_help,
        showModal(modalDialog(options_help_dialog))
    )
    observeEvent( input$group_posterior_help,
        showModal(modalDialog(group_posterior_help_dialog))
    )
    ################################
    
    ##### Hyperprior plots #####
    output$hyperprior_plots = renderPlot({
        datasetupload()
        arglist = list(
            dprime_mu_hyperloc  = input$Group_dprime_mu_hyperloc, 
            dprime_mu_hyperscale = input$Group_dprime_mu_hyperscale, 
            dprime_sd_hyperloc = input$Group_dprime_sd_hyperloc, 
            dprime_sd_hyperscale = input$Group_dprime_sd_hyperscale,
            c_mu_hyperloc = input$Group_c_mu_hyperloc, 
            c_mu_hyperscale = input$Group_c_mu_hyperscale, 
            c_sd_hyperloc = input$Group_c_sd_hyperloc, 
            c_sd_hyperscale = input$Group_c_sd_hyperscale,
            dataset = session_storage$dataset
        )
        session_storage$hyperprior_plot = arglist
        do.call( "plot_hyperpriors", session_storage$hyperprior_plot )
    })
    
    observeEvent( input$hyperprior_plot_popout, {
        dev.new()
        do.call( "plot_hyperpriors", session_storage$hyperprior_plot )
    } )
    ################################
    
    ##### Group posterior plots #####
    output$group_posterior_plots = renderPlot({
        posterior_results()
        if ( !is.null(session_storage$model_fits) ){
            arglist = list( model_fits = session_storage$model_fits,
                prior_things = session_storage$hyperprior_plot )
            do.call( "plot_group_posteriors", arglist )
        } else{
            null_plot("Posterior results will be shown here")
        }
    })
    
    observeEvent( input$group_posterior_plot_popout, {
        dev.new()
        if ( !is.null(session_storage$model_fits) ){
            arglist = list( model_fits = session_storage$model_fits,
                prior_things = session_storage$hyperprior_plot )
            do.call( "plot_group_posteriors", arglist )
        } else{
            null_plot("Posterior results will be shown here")
        }
    } )
    
    output$posterior_group_text = renderPrint({
        posterior_results()
        if ( !is.null(session_storage$group_posterior_summaries) ){
            print(session_storage$group_posterior_summaries)
        } else{
            print("Text summaries will be shown here")
        }
    })
    #################################
    
    ##### Individual posterior plots #####
    output$individual_posterior_plots_dprime = renderPlot({
        posterior_results()
        if ( !is.null(session_storage$model_fits) ){
            arglist = list(
                model_fits=session_storage$model_fits,
                prior_things = session_storage$hyperprior_plot
                )
            do.call( "plot_individual_posteriors_dprime", arglist )
        } else{
            null_plot("Posterior results will be shown here")
        }
    })
    observeEvent( input$individual_posterior_plots_dprime_popout, {
        dev.new()
        if ( !is.null(session_storage$model_fits) ){
            arglist = list(
                model_fits=session_storage$model_fits,
                prior_things = session_storage$hyperprior_plot
                )
            do.call( "plot_individual_posteriors_dprime", arglist )
        } else{
            null_plot("Posterior results will be shown here")
        }
    })
    output$posterior_individual_text_dprime = renderPrint({
        posterior_results()
        if ( !is.null(session_storage$individual_posterior_summaries_dprime) ){
            print(session_storage$individual_posterior_summaries_dprime)
        } else{
            print("Text summaries will be shown here")
        }
    })
    
    output$individual_posterior_plots_c = renderPlot({
        posterior_results()
        if ( !is.null(session_storage$model_fits) ){
            arglist = list(
                model_fits=session_storage$model_fits,
                prior_things = session_storage$hyperprior_plot
                )
            do.call( "plot_individual_posteriors_c", arglist )
        } else{
            null_plot("Posterior results will be shown here")
        }
    })
    observeEvent( input$individual_posterior_plots_c_popout, {
        dev.new()
        if ( !is.null(session_storage$model_fits) ){
            arglist = list(
                model_fits=session_storage$model_fits,
                prior_things = session_storage$hyperprior_plot
                )
            do.call( "plot_individual_posteriors_c", arglist )
        } else{
            null_plot("Posterior results will be shown here")
        }
    })
    output$posterior_individual_text_c = renderPrint({
        posterior_results()
        if ( !is.null(session_storage$individual_posterior_summaries_c) ){
            print(session_storage$individual_posterior_summaries_c)
        } else{
            print("Text summaries will be shown here")
        }
    })
    #####################################
    
    
    observeEvent( input$simulate_data, {
        param_list = list(
            group_dprime_mu = input$Group_dprime_mean,
            group_dprime_sd = input$Group_dprime_sd,
            group_c_mu = input$Group_c_sd,
            group_c_sd = input$Group_c_sd,
            n_subjects = input$n_subjects,
            n_signal_trials = input$n_signal_trials,
            n_noise_trials = input$n_noise_trials
        )
        datasetupload( datasetupload() + 1 )
        dataset = do.call( "make_data", param_list )
        session_storage$true_params = dataset$true_params
        session_storage$dataset = dataset$dataset
        
        output$data_summary = renderDataTable( dataset$dataset )
    } )
    
    output$download_data = downloadHandler(
        filename = function(){
            filename = paste0("sim_SDT_data_", format(Sys.time(), "%H_%M_%d_%m_%y"), ".csv" )
            filename
        },
        content = function(file){
            if ( is.null(session_storage$dataset) ){
                showModal(modalDialog("File contains no data."))
            }
            write.csv( session_storage$dataset, file, row.names = FALSE )
        }
    )
    
    output$save_Rdata = downloadHandler(
        filename = function(){
            filename = paste0("sim_SDT_data_", format(Sys.time(), "%H_%M_%d_%m_%y"), ".RData" )
            filename
        },
        content = function(file){
            ls(session_storage)
            save_vars = c("model_fits", "dataset", "individual_posterior_summaries_c",
                "individual_posterior_summaries_dprime", "group_posterior_summaries" )
            save( list = save_vars, envir = session_storage, file = "x.RData" )
        }
    )
    
    observeEvent( input$file_import, {
        import_filepath = input$file_import$datapath
        session_storage$import_filepath = import_filepath
        import_dataframe = try(read.csv( import_filepath ))
        
        import_dataframe_check = restructure_import_df(import_dataframe)
        if ( import_dataframe_check$error ){
            showModal(modalDialog(import_dataframe_check$error_message))
        } else{
            session_storage$dataset = import_dataframe_check$df
            datasetupload( datasetupload() + 1 )
            output$data_summary = renderDataTable( session_storage$dataset )
        }
    } )
    
    observeEvent( input$clear_dataset, {
        datasetupload( datasetupload() + 1 )
        session_storage$dataset = NULL
        output$data_summary = renderDataTable( data.frame( "No data loaded" = "No data loaded" ) )
    } )
    
    observeEvent(input$run_mcmc, {
        if ( is.null(session_storage$hyperprior_plot$dataset) ){
            showModal(modalDialog("Model fitting failed because there is no dataset."))
        } else{
            showModal(modalDialog("Model is now being fit in stan. Check your R console for progress outputs."))
            mcmc_options = list(
                adapt_delta  = input$mcmc_adapt_delta,
                max_treedepth = input$mcmc_max_treedepth,
                iter = input$mcmc_max_iter,
                warmup_iter = input$mcmc_warmup_iter,
                chains = input$mcmc_chains
            )
            session_storage$mcmc_options = mcmc_options
            session_storage$model_fits = do.call( "run_mcmc", c( mcmc_options, session_storage$hyperprior_plot ) )
            session_storage$group_posterior_summaries = text_group_posteriors(session_storage$model_fits, 
                session_storage$hyperprior_plot)
            session_storage$individual_posterior_summaries_dprime = text_individual_posteriors_dprime(
                session_storage$model_fits, session_storage$hyperprior_plot
            )
            session_storage$individual_posterior_summaries_c = text_individual_posteriors_c(
                session_storage$model_fits, session_storage$hyperprior_plot
            )
            posterior_results( posterior_results() + 1)
        }
    } )
}

shinyApp( ui = ui, server = server )
