#' bayesNMF_sampler R6 class
#' 
#' @details
#' This class implements the Gibbs sampler for Bayesian Nonnegative Matrix Factorization with automatically learned rank and optional Metropolis-Hastings updates for computational efficiency.
#' It holds the data, model specifications, current sampler state, sampler history, inference results, and various utility functions for visualization and postprocessing.
#' 
#' @export
bayesNMF_sampler <- R6::R6Class(
  "bayesNMF_sampler",
  public = list(
    #' @field data matrix, observed data (rows are variables, columns are samples)
    data = NULL,

    #' @field dims list, dimensions of the data (K is variables, N is latent dimension, G is samples)
    dims = list(),

    #' @field specs list, model specifications (rank, likelihood, prior, MH, learning_rank, convergence_control, output_dir, overwrite, verbosity, periodic_save, save_all_samples)
    specs = list(),

    #' @field temperature_schedule vector, temperature schedule (temperature_schedule[i] is the temperature for iteration i)
    temperature_schedule = NULL,

    #' @field hyperprior_params list, fixed hyperprior parameters
    hyperprior_params = list(),

    #' @field prior_params list, current values of prior parameters
    prior_params = list(),

    #' @field params list, current values of parameters
    params = list(),

    #' @field samples list, sampler history (samples[[name]][[iter]] is the value of parameter name at iteration iter)
    samples = list(),

    #' @field acceptance_rates list, acceptance rates for Metropolis-Hastings updates (all 1 if MH is FALSE)
    acceptance_rates = list(P_acceptance_rate = NULL, E_acceptance_rate = NULL),

    #' @field state list, current sampler state (iter is the current iteration, indent is the current indentation level, converged is whether the sampler has converged)
    state = list(
      iter = 1,
      indent = 0,
      converged = FALSE
    ),

    #' @field MAP list, current maximum a posteriori (MAP) estimate from most recent convergence check
    MAP = list(assignment_res = NULL),

    #' @field credible_intervals list, 95% credible intervals for the current MAP estimates (credible_intervals[[name]][["lower"]] and credible_intervals[[name]][["upper"]] are the bounds for parameter name)
    credible_intervals = list(),

    #' @field reference_comparison list, includes assignments, votes, and plots for a given reference set
    reference_comparison = list(
      reference_P = NULL,
      assignments = NULL,
      keep_sigs = NULL,
      idxs = NULL,
      votes = NULL,
      summary = NULL,
      plots = list(),
      label_switching_df = NULL
    ),

    #' @field time list, time taken for various operations (total is the total time taken, per_iter is the average time per iteration)
    time = list(),

    #' @field log_con file, connection to log file
    log_con = NULL,

    ########################################
    ########## KEY UTILITIES ###############
    ########################################
    
    #' @description
    #' Initialize a bayesNMF_sampler object
    #' @param data matrix, observed data (rows are variables, columns are samples)
    #' @param rank integer, number of signatures or range of ranks to consider (default 1:10)
    #' @param likelihood string, likelihood distribution for M, one of "poisson", "normal" (default "poisson")
    #' @param prior string, prior distribution for P and E, one of "truncnormal", "exponential", "gamma" (default "truncnormal")
    #' @param rank_method string, method to determine rank, one of "SBFI", "BFI", "BIC" (heuristic) (default "SBFI")
    #' @param MH boolean, if likelihood is Poisson and prior is "truncnormal" or "exponential", whether to use Metropolis-Hastings updates described in the paper (default TRUE when possible)
    #' @param convergence_control list, convergence control parameters (see `new_convergence_control()` for default values)
    #' @param prop_temp float, proportion of iterations to use for tempering if learning rank (default 0.2)
    #' @param post_warmup integer, number of iterations post-warmup to sample from (if MH = TRUE) (default 2 * convergence_control$MAP_over)
    #' @param output_dir string, output directory, will be created if it does not exist
    #' @param overwrite boolean, whether to overwrite existing output, will delete existing output directory if it exists (default FALSE)
    #' @param hyperprior_params list, fixed hyperprior parameters (will be filled with default values if not provided)
    #' @param init_prior_params list, initial prior parameters (will be initialized from hyperpriors if not provided)
    #' @param init_params list, initial parameters (will be initialized from priors if not provided)
    #' @param verbosity integer, verbosity level (default 1)
    #' @param periodic_save boolean, whether to save the sampler object periodically (default TRUE, set to FALSE to reduce I/O and computational time)
    #' @param save_all_samples boolean, whether to save all samples (default TRUE, set to FALSE to reduce memory usage)
    initialize = function(
      data, rank,
      likelihood = "poisson",
      prior = "truncnormal",
      rank_method = "SBFI",
      MH = likelihood == "poisson" & prior %in% c('truncnormal','exponential'),
      convergence_control = new_convergence_control(),
      prop_temp = 0.2,
      post_warmup = 2 * convergence_control$MAP_over,
      output_dir = paste0('nmf_', likelihood, '_', prior),
      overwrite = FALSE,
      hyperprior_params = list(),
      init_prior_params = list(),
      init_params = list(),
      verbosity = 1,
      periodic_save = TRUE,
      save_all_samples = FALSE
    ) {
      # add extension to output dir to avoid overwrite
      final_dir = output_dir; tail = 0
      while (!overwrite & dir.exists(final_dir)) {
        tail = tail + 1
        final_dir = paste0(output_dir, '_', tail)
      }
      # empty outputdir if it exsits and overwrite is TRUE
      if (overwrite & dir.exists(final_dir)) {
        unlink(final_dir, recursive = TRUE)
      }
      # create output dir if it doesn't exist
      if (!dir.exists(final_dir)) { dir.create(final_dir, recursive = TRUE) }
  
      # check if rank is fixed or will be learned
      learning_rank <- length(rank) > 1
      if (learning_rank & min(rank) != 0) { rank = 0:max(rank) }
      display_rank = ifelse(learning_rank, paste0(min(rank), ':', max(rank)), rank)
      
      # set tempering schedule
      n_iters = convergence_control$maxiters + ifelse(MH, post_warmup, 0)
      if (learning_rank) {
        self$temperature_schedule <- private$get_temp_sched(
          len = n_iters,
          n_temp = round(prop_temp * convergence_control$maxiters)
        )
      } else {
        self$temperature_schedule <- rep(1, n_iters)
      }

      # set provided attributes
      self$data <- data
      self$dims <- list(
        K = nrow(data),
        N = max(rank),
        G = ncol(data)
      )
      self$specs <- list(
        rank = rank,
        likelihood = likelihood,
        prior = prior,
        MH = MH,
        learning_rank = learning_rank,
        convergence_control = convergence_control,
        output_dir = final_dir,
        overwrite = overwrite,
        verbosity = verbosity,
        periodic_save = periodic_save,
        save_all_samples = save_all_samples
      )
      if (learning_rank) {
        self$specs$prop_temp <- prop_temp
        self$specs$rank_method <- rank_method
      }
      if (MH) {
        self$specs$post_warmup <- post_warmup
      }

      # set up logging
      self$log_con <- file(file.path(final_dir, "log.txt"), open = "wt")
      self$log("Initialized sampler", verbosity = 1)
      self$state$indent <- 1
      self$log(glue::glue("likelihood = {self$specs$likelihood}, prior = {self$specs$prior}, MH = {self$specs$MH}"), verbosity = 1)
      self$log(glue::glue("learning_rank = {learning_rank}, rank = {display_rank}"), verbosity = 1)
      self$log(glue::glue("maxiters = {self$specs$convergence_control$maxiters}"), verbosity = 1)
      self$log(glue::glue("MAP_over = {self$specs$convergence_control$MAP_over}"), verbosity = 1)
      self$log(glue::glue("MAP_every = {self$specs$convergence_control$MAP_every}"), verbosity = 1)
      
      # set up state
      self$state[['MAP_metrics']] <- data.frame(
        iter = numeric(),
        RMSE = numeric(),
        KL = numeric(),
        loglikelihood = numeric(),
        logposterior = numeric(),
        n_params = numeric(),
        BIC = numeric(),
        rank = numeric(),
        MAP_A_counts = numeric(),
        mean_temp = numeric()
      )
      self$state[['sample_metrics']] <- data.frame(
        iter = numeric(),
        RMSE = numeric(),
        KL = numeric(),
        loglikelihood = numeric(),
        logposterior = numeric(),
        n_params = numeric(),
        BIC = numeric(),
        rank = numeric(),
        temp = numeric()
      )
      self$state[['MAP_idx']] <- 1:self$specs$convergence_control$MAP_over
      if (self$specs$MH) {
        self$state[['sample_metrics']]$P_mean_acceptance_rate <- numeric()
        self$state[['sample_metrics']]$E_mean_acceptance_rate <- numeric()
        self$state[['MAP_metrics']]$P_mean_acceptance_rate <- numeric()
        self$state[['MAP_metrics']]$E_mean_acceptance_rate <- numeric()
      }

      self$state$indent <- 0
      self$log(glue::glue("Setup"), verbosity = 1)
      self$state$indent <- 1
      self$log("Setting hyperprior parameters", verbosity = 1)
      self$hyperprior_params <- hyperprior_params
      self$hyperprior_params <- private$fill_hyperprior_params()

      # check model definition is available
      private$check_model()
      self$log("Model check passed", verbosity = 1)

      # set initial values and sample other prior parameters and parameters
      self$log("Initializing prior parameters and parameters", verbosity = 1)
      if (self$specs$likelihood == "normal") {
        # set default alpha and beta priors for sigmasq if not provided
        if (!("alpha" %in% names(init_prior_params))) {
          init_prior_params$alpha <- 3
        }
        if (!("beta" %in% names(init_prior_params))) {
          init_prior_params$beta <- 3
        }
      }

      private$init_prior_params(init_prior_params)
      private$init_params(init_params)

      if (self$specs$MH) {
        self$acceptance_rates$P_acceptance_rate <- matrix(nrow = self$dims$K, ncol = self$dims$N)
        self$acceptance_rates$E_acceptance_rate <- matrix(nrow = self$dims$N, ncol = self$dims$G)
      }

      self$log("Sampling parameters from priors", verbosity = 1)
      private$sample_params(skip = names(init_params), from_prior = TRUE)

      # record initial sample to stored samples
      self$log("Logging initial sample", verbosity = 1)
      self$samples <- list()
      for (name in names(self$params)) { self$samples[[name]] <- list() }
      for (name in names(self$prior_params)) { self$samples[[name]] <- list() }
      if (self$specs$MH) {
        for (name in names(self$acceptance_rates)) {
          self$samples[[name]] <- list()
        }
      }
      private$record_sample()

      # compute metrics for initial sample
      self$log("Logging initial sample metrics", verbosity = 1)
      private$update_sample_metrics(update_trace = FALSE)

      self$state$indent <- 0
    },

    #' @description
    #' Perform Gibbs sampling
    #' Samples until convergence or maxiters is reached, plus additional post_warmup samples after convergence if MH is TRUE
    run_gibbs_sampler = function() {
      self$log("Starting Gibbs sampler", verbosity = 1)
      start_time <- Sys.time()
      while (
        !self$state$converged &
        self$state$iter < self$specs$convergence_control$maxiters
      ) {
        self$state$indent <- 1
        self$state$iter <- self$state$iter + 1
        self$log("sampling prior parameters", verbosity = 3)
        private$sample_prior_params()
        self$log("sampling parameters", verbosity = 3)
        private$sample_params()
        self$log("recording sample", verbosity = 3)
        private$record_sample()
        private$update_sample_metrics(
          update_trace = (
            self$state$iter %% self$specs$convergence_control$MAP_every == 0 &
            self$specs$periodic_save
          )
        )

        # check convergence
        if (
          (self$state$iter %% self$specs$convergence_control$MAP_every == 0 &
            self$state$iter >= max(
              self$specs$convergence_control$MAP_over,
              self$specs$convergence_control$MAP_every
            )
          )
          | self$state$iter >= self$specs$convergence_control$maxiters
        ) {
          if (self$specs$save_all_samples) {
            # update MAP_idx to last MAP_over samples
            self$state$MAP_idx <- seq(
              from = self$state$iter -
                self$specs$convergence_control$MAP_over + 1,
              to = self$state$iter
            )
          }
          
          self$log(glue::glue("iter = {self$state$iter}"), verbosity = 1)
          self$state$indent <- 2
          
          self$log("Computing MAP", verbosity = 1)
          self$get_MAP()
          if (self$specs$learning_rank) {
            self$log(log_table(t(self$MAP$A_counts)), verbosity = 1)
          }

          self$log("Checking convergence", verbosity = 1)
          msg = private$check_convergence()
          self$log(msg, verbosity = 1)

          self$state$indent <- 1
          if (self$state$converged) {
            self$state$converged_iter <- self$state$iter
            self$log(glue::glue("Converged at {self$state$iter} due to {self$state$why}"), verbosity = 1)
          }

          if (self$specs$periodic_save) {
            self$log("Saving object", verbosity = 1)
            self$save_object()
          }
        }
      }
      self$state$indent <- 1
      if (self$specs$MH) {
        start_MH <- Sys.time()
        self$time$warmup <- difftime(start_MH, start_time, units = 'mins')
        self$log(glue::glue("Warmup done, sampling {self$specs$post_warmup} with MH for inference"), verbosity = 1)

        for (i in 1:self$specs$post_warmup) {
          self$state$indent <- 1
          self$state$iter <- self$state$iter + 1
          private$sample_prior_params()
          private$sample_params()
          private$record_sample()
          private$update_sample_metrics(
            update_trace = (
              i %% self$specs$convergence_control$MAP_every == 0 &
              self$specs$periodic_save
            ) | i == self$specs$post_warmup
          )

          if (
            self$state$iter %% self$specs$convergence_control$MAP_every == 0 |
            i == self$specs$post_warmup
          ) {
            if (self$specs$save_all_samples) {
              # update MAP_idx to last MAP_over samples
              self$state$MAP_idx <- seq(
                from = self$state$iter -
                  self$specs$convergence_control$MAP_over + 1,
                to = self$state$iter
              )
            }

            final = i == self$specs$post_warmup
            self$log(glue::glue("iter = {self$state$iter}"), verbosity = 1)
            self$state$indent <- 2
            self$log("Computing MAP", verbosity = 1)
            # final MAP filters to included signatures only
            self$get_MAP(final = final)
            if (self$specs$learning_rank) {
              self$log(log_table(t(self$MAP$A_counts)), verbosity = 1)
            }
            
            # just to update MAP metrics (already converged)
            msg = private$check_convergence(final = final)

            if (self$specs$periodic_save) {
              self$log("Saving object", verbosity = 1)
              self$save_object()
            }
          }
        }
        self$state$indent <- 1
        self$log(glue::glue("Additional {self$specs$post_warmup} MH samples done"), verbosity = 1)
        self$time$MH <- difftime(Sys.time(), start_MH, units = 'mins')
      } else {
        # final MAP filters to included signatures only
        self$get_MAP(final = TRUE)
        if (self$specs$learning_rank) {
          self$log(log_table(t(self$MAP$A_counts)), verbosity = 1)
        }
        self$log("Final MAP computed", verbosity = 1)

        # final updates of trace plots
        trace_plot(self, MAP_means = TRUE, save = TRUE)
        trace_plot(self, save = TRUE)

      }
      self$log("Sampler done", verbosity = 1)
      self$time$total <- difftime(Sys.time(), start_time, units = 'mins')
      self$time$per_iter <- self$time$total / self$state$iter
      self$log(glue::glue(
        "Total time: {round(self$time$total, 2)} minutes"
      ), verbosity = 1)

      # save final object
      self$log("Saving final object", verbosity = 1)
      self$save_object()
    },

    #' @description
    #' Save the sampler object to the output directory
    #' 
    #' @return None, saves the sampler object to the output directory
    save_object = function() {
      saveRDS(self, file.path(self$specs$output_dir, "sampler.rds"))
    },

    #' @description
    #' Log a message to the log file
    #' 
    #' @param msg string, message to log
    #' @param verbosity integer, verbosity level (default 5, only logs messages with verbosity <= self$specs$verbosity)
    log = function(msg, verbosity = 5) {
      if (verbosity > self$specs$verbosity) return()
      if (length(msg) == 0) return()   # nothing to log
      if (length(msg) == 1 && grepl("\n", msg)) {
        # split single string with newlines
        msg_lines <- strsplit(msg, "\n", fixed = TRUE)[[1]]
      } else {
        msg_lines <- msg
      }

      # format message with timestamp and indent
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      indent <- paste(rep("\t", self$state$indent), collapse = "")
      indent_vec <- sapply(1:length(msg_lines), function(i) {
        if (i > 1) {
          indent_i <- paste(c(indent, rep(" ", nchar(timestamp) + 1)), collapse = "")
        } else {
          indent_i <- indent
        }
        return(indent_i)
      })

      # combine indented lines with timestamp
      msg_lines <- paste0(indent_vec, msg_lines)
      msg_lines <- msg_lines[trimws(msg_lines) != ""]
      msg_lines <- paste(msg_lines, collapse = "\n")

      # write to log file
      if (!is.null(self$log_con)) {
        writeLines(sprintf("[%s] %s", timestamp, msg_lines), con = self$log_con)
        flush(self$log_con)
      }
    },

    ########################################
    ######## FUNCTIONS FROM UTILS.R ########
    ########################################

    #' @description
    #' Compute the expected data matrix
    #' 
    #' @param P matrix, optional, dimensions K x N, uses self$params$P if not provided
    #' @param A matrix, optional, dimensions 1 x N, uses self$params$A if not provided
    #' @param E matrix, optional, dimensions N x G, uses self$params$E if not provided
    #' 
    #' @return matrix, dimensions K x G
    get_Mhat = function(P = NULL, A = NULL, E = NULL) {
      get_Mhat_(self, P, A, E)
    },

    #' @description
    #' Compute the log likelihood
    #' 
    #' @param P matrix, optional, dimensions K x N, uses self$params$P if not provided
    #' @param A matrix, optional, dimensions 1 x N, uses self$params$A if not provided
    #' @param E matrix, optional, dimensions N x G, uses self$params$E if not provided
    #' @param sigmasq vector, optional, dimensions G, uses self$params$sigmasq if not provided, only used if likelihood is "normal"
    #' @param likelihood string, likelihood distribution for M, one of "poisson", "normal", uses self$specs$likelihood if not provided
    #' @param return_matrix boolean, whether to return the log likelihood matrix, used to compute acceptance ratios in MH updates (default FALSE)
    #' 
    #' @return scalar log likelihood if return_matrix is FALSE, otherwise matrix of log likelihoods per observed value (K x G)
    get_loglik = function(P = NULL, A = NULL, E = NULL, sigmasq = NULL, likelihood = self$specs$likelihood, return_matrix = FALSE) {
      get_loglik_(self, P, A, E, sigmasq, likelihood, return_matrix)
    },

    #' @description
    #' Compute the log posterior using model specifications
    #' 
    #' @param P matrix, optional, dimensions K x N, uses self$params$P if not provided
    #' @param A matrix, optional, dimensions 1 x N, uses self$params$A if not provided
    #' @param E matrix, optional, dimensions N x G, uses self$params$E if not provided
    #' @param sigmasq vector, optional, dimensions G, uses self$params$sigmasq if not provided, only used if likelihood is "normal"
    #' 
    #' @return scalar log posterior
    get_logpost = function(P = NULL, A = NULL, E = NULL, sigmasq = NULL) {
      get_logpost_(self, P, A, E, sigmasq)
    },

    #' @description
    #' Compute the maximum a posteriori (MAP) estimate
    #' 
    #' @param end_iter integer, last iteration to consider for inference, defaults to current iteration
    #' @param n_samples integer, number of samples to consider for inference, defaults to MAP_over specified in convergence control
    #' @param final boolean, if TRUE, subset to only included signatures
    #' @param credible_interval float, credible interval width (default 0.95)
    #' 
    #' @return None, updates self$MAP and self$credible_intervals
    get_MAP = function(
      end_iter = self$state$iter,
      n_samples = self$specs$convergence_control$MAP_over,
      final = FALSE, credible_interval = 0.95
    ) {
      get_MAP_(self, end_iter, n_samples, final, credible_interval)
    },

    ########################################
    ### FUNCTIONS FROM POSTPROCESSING.R ####
    ########################################

    #' @description
    #' Assign signatures for each posterior sample based on cosine similarity and ensemble with majority vote
    #'
    #' @param reference_P matrix, "cosmic", or NULL, reference signatures to align to
    #' @param idxs vector of indices to consider (default "MAP_idx" indicates indices used to compute most recent MAP)
    #' @param credible_interval numeric, credible interval to compute for cosine similarities (default 0.95)
    #'
    #' @return list of two data frames:
    #' \itemize{
    #'   \item assignments: holds final assignments for each signature, cosine similarities between MAP estimates and reference signatures, and 95% credible intervals for the cosine similarities. 
    #'   \item votes: holds the proportion of votes for each signature assignment for each posterior sample. This includes signatures that did not receive majority vote, allowing users to understand posterior uncertainty in assignment.
    #' }
    #' and stores the result in self$reference_comparison
    assign_signatures_ensemble = function(
      reference_P = "cosmic",
      idxs = self$MAP$idx,
      credible_interval = 0.95
    ) {
      assign_signatures_ensemble_(self, reference_P, idxs, credible_interval)
    }
  ),

  private = list(
    # ROxygen documentation not allowed for private methods, so using #"
    # Keeping documentation for clarity of code for future contributors

    ########################################
    ### FUNCTIONS SETUP.R ##################
    ########################################

    #" @description
    #" Fill hyperprior parameters
    #" Based on user-provided hyperprior parameters and default hyperprior parameters, create full hyperprior parameters matrices
    #" 
    #" @return None, updates self$hyperprior_params
    fill_hyperprior_params = function() { fill_hyperprior_params_(self) },

    ########################################
    ### FUNCTIONS SAMPLE_PRIORS.R #########
    ########################################

    #" @description
    #" Initialize prior parameters based on specified initial values or sampled from hyperprior distributions if not provided
    #" 
    #" @return None, updates self$prior_params
    init_prior_params = function(init_prior_params) { init_prior_params_(self, init_prior_params) },

    #" @description
    #" Sample prior parameters from full conditionals
    #" 
    #" @param skip list, names of prior parameters to skip
    #" 
    #" @return None, updates self$prior_params
    sample_prior_params = function(skip = c()) {
      sample_prior_params_(self, skip)
    },

    ########################################
    ### FUNCTIONS SAMPLE_PARAMS.R #########
    ########################################

    #" @description
    #" Initialize parameters based on specified initial values or sampled from prior distributions if not provided
    #" 
    #" @param init_params list, initial values for parameters
    #" 
    #" @return None, updates self$params
    init_params = function(init_params) { init_params_(self, init_params) },

    #" @description
    #" Sample parameters from full conditionals (one iteration of Gibbs sampling)
    #" 
    #" @param skip list, names of parameters to skip
    #" @param from_prior boolean, whether to sample directlyfrom prior distributions
    #" 
    #" @return None, updates self$params
    sample_params = function(skip = c(), from_prior = FALSE) {
      sample_params_(self, skip, from_prior)
    },

    #" @description
    #" Custom error logger
    #" 
    #" @param msg string, message to log
    #" 
    #" @return None, stops execution with error message
    error = function(msg) {
      msg = glue::glue(msg)
      msg = paste("ERROR:", msg)
      if (!is.null(self$log_con)) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        writeLines(sprintf("[%s] %s", timestamp, msg), con = self$log_con)
        flush(self$log_con)
      }
      stop(msg)
    },

    #" @description
    #" Check model specifications
    #" 
    #" @return None, stops execution if model specifications are invalid
    check_model = function() {
      available_likelihoods = c('normal','poisson')

      if (!(self$specs$likelihood %in% available_likelihoods)) {
        private$error("likelihood must be one of {paste(available_likelihoods, collapse = ', ')}")
      }

      if (self$specs$likelihood == 'normal') {
        if (!(self$specs$prior %in% c('truncnormal','exponential'))) {
          private$error("prior must be one of c('truncnormal','exponential') with `likelihood = 'normal'`")
        }
      } else if (self$specs$likelihood == 'poisson') {
        if (!(self$specs$prior %in% c('gamma','exponential','truncnormal'))) {
          private$error("prior must be one of c('gamma','exponential','truncnormal') with `likelihood = 'poisson'`")
        }
        if (self$specs$prior == 'gamma' & self$specs$MH) {
          private$error('gamma prior cannot be used in a MH-within-gibbs sampler')
        }
        if (self$specs$prior == 'truncnormal' & !self$specs$MH) {
          private$error('truncnormal prior can only be used in a MH-within-gibbs sampler')
        }
      }
    },

    #" @description
    #" Record current iteration of parameters and prior parameters
    #" 
    #" @return None, updates self$samples
    record_sample = function() {
      lists <- list('params' = self$params, 'prior_params' = self$prior_params)
      if (self$specs$MH) {
        lists$acceptance_rates <- self$acceptance_rates
      }
      for (li in lists) {
        for (name in names(li)) {
          if (!self$specs$save_all_samples) {
            # update_list keeps only the last MAP_over samples
            # if not saving all samples
            self$samples[[name]] <- update_list(
              list = self$samples[[name]],
              new_value = li[[name]],
              index = self$state$iter,
              max_length = self$specs$convergence_control$MAP_over
            )
          } else {
            self$samples[[name]][[self$state$iter]] <- li[[name]]
          }
        }
      }
    },


    ########################################
    ######## FUNCTIONS FROM UTILS.R ########
    ########################################

    #" @description
    #" Get tempering schedule
    #" Temperature parameter is used to control the balance between exploration and exploitation in the sampler.
    #" The temperature starts at 0 and slowly increases to 1 over the course of the first `n_temp` iterations.
    #" 
    #" @param len integer, total number of iterations
    #" @param n_temp integer, number of temperature levels
    #" 
    #" @return vector of temperatures, length `len`
    get_temp_sched = function(len, n_temp) {
      get_temp_sched_(len, n_temp)
    },

    #" @description
    #" Update sample metrics
    #" 
    #" @param self bayesNMF_sampler object
    #" @param update_trace boolean, whether to update the trace plot (default FALSE)
    #" 
    #" @return None, updates self$state$sample_metrics and saves trace plot if update_trace is TRUE
    update_sample_metrics = function(update_trace = FALSE) {
      update_sample_metrics_(
        self,
        update_trace = update_trace
      )
    },

    #" @description
    #" Update MAP metrics
    #" 
    #" @param self bayesNMF_sampler object
    #" @param final boolean, if TRUE, subset to only included signatures
    #" 
    #" @return None, updates self$state$MAP_metrics and saves trace plot if periodic save or final is TRUE
    update_MAP_metrics = function(final = FALSE) {
      update_MAP_metrics_(self, final = final)
    },

    ########################################
    ##### FUNCTIONS FROM CONVERGENCE.R #####
    ########################################

    #" @description
    #" Check for convergence
    #" 
    #" @param final boolean, if TRUE, subset to only included signatures
    #" 
    #" @return string, message indicating convergence status
    check_convergence = function(final = FALSE) {
      check_convergence_(self, private, final = final)
    },


    ########################################
    ############### CLEAN UP ###############
    ########################################

    #" @description
    #" Finalize the sampler by closing the log file
    #" Automatically called when the sampler is destroyed
    #" @return None, closes the log file
    finalize = function() {
      if (!is.null(self$log_con)) {
        writeLines("[INFO] Finalizing and closing log file", self$log_con)
        close(self$log_con)
      }
    }
  )
)