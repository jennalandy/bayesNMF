#' Bayesian Non-Negative Matrix Factorization
#' @description Perform single-study Bayesian NMF with the provided likelihood and prior
#' combination. Exact rank `N` or maximum rank `max_N` must be provided.
#'
#' @param M mutational catalog matrix, K x G
#' @param N fixed number of latent factors
#' @param max_N maximum number of latent factors if learning rank
#' @param likelihood string, one of c('normal','poisson')
#' @param prior string, one of c('truncnormal','exponential')
#' @param inits (optional) list of initial values for P and E as well as sigmasq
#' if `likelihood = "normal"`
#' @param fixed (ptional) list of parameters to fix and not include in Gibbs
#' updates.
#' @param prior_parameters list, optional specification of prior parameters
#' @param file file name without extension of log, save, and plot files
#' @param true_P (optional) true latent factors matrix P to compare to in a heatmap
#' @param convergence_control list, specification of convergence parameters.
#' See documentation for `new_convergence_control`.
#' @param store_logs boolean, whether to store each sample in resulting .rds file
#' @param overwrite if `overwrite = TRUE`, the log, safe, and plot files of
#' previous runs with the same `file` will be overwritten
#'
#' @return list
#' @export
bayesNMF <- function(
        M,
        N = NULL,
        max_N = NULL,
        likelihood = "normal",
        prior = "truncnormal",
        inits = NULL,
        fixed = NULL,
        prior_parameters = NULL,
        recovery = FALSE,
        recovery_priors = "cosmic",
        file = paste0('nmf_', likelihood, '_', prior),
        true_P = NULL,
        convergence_control = new_convergence_control(),
        store_logs = TRUE,
        overwrite = FALSE
) {
    START = Sys.time()

    # rescale input data to fit scale of priors
    rescale_by = mean(M)/100
    M_truescale = M
    M = M/rescale_by
    if (!is.null(fixed$sigmasq)) {
        fixed$sigmasq = fixed$sigmasq/(rescale_by**2)
    }
    if (!is.null(fixed$E)) {
        fixed$E <- fixed$E/rescale_by
    }

    # check recovery/discovery
    if (recovery) {
        if (is.character(recovery_priors)) {
            if (recovery_priors == "cosmic") {
                if (likelihood == 'normal' & prior == 'truncnormal') {
                    recovery_priors <- normal_truncnormal_recovery_priors
                } else if (likelihood == 'normal' & prior == 'exponential') {
                    recovery_priors <- normal_exponential_recovery_priors
                } else {
                    stop("Recovery priors not defined for this likelihood/prior combination")
                }
            }
        }
    } else {
        recovery_priors = list()
    }

    # check N/max_N combination is valid
    N <- validate_N(N, max_N, recovery_priors)

    # check prior and likelihood are valid
    validate_model(likelihood, prior)

    # set up tempering schedule
    learn_A <- !is.null(max_N) & is.null(fixed$A)
    if (learn_A) {
        gamma_sched <- get_gamma_sched(len = convergence_control$maxiters)
    } else {
        gamma_sched <- rep(1, convergence_control$maxiters)
    }

    # precompute log factorials for Poisson likelihood
    if (likelihood == 'poisson') {
        logfac = vector(length = max(M))
        logfac[1] = 0
        for (i in 2:length(logfac)) {
            logfac[i] = log(i) + logfac[i-1]
        }
    } else {
        logfac = NULL
    }

    # set up dimensions
    dims = list(
        K = dim(M)[1],
        G = dim(M)[2],
        N = N,
        S = 1
    )

    # set up file names
    savefile = paste0(file, '.rds')
    logfile = paste0(file, '.log')
    plotfile = paste0(file, '.pdf')
    tail = 0
    while (!overwrite & (file.exists(savefile) | file.exists(logfile))) {
        tail = tail + 1
        savefile = paste0(file, '_', tail, '.rds')
        logfile = paste0(file, '_', tail, '.log')
        plotfile = paste0(file, '_', tail, '.pdf')
    }

    # set up Theta
    Theta <- initialize_Theta(
        M = M,
        likelihood = likelihood,
        prior = prior,
        learn_A = learn_A,
        dims = dims,
        inits = inits, fixed = fixed,
        prior_parameters = prior_parameters,
        recovery = recovery,
        recovery_priors = recovery_priors
    )
    Theta$prob_inclusion <- Theta$q

    # set up metrics and logs
    metrics <- list(
        sample_idx = list(),
        RMSE = list(),
        KL = list(),
        loglik = list(),
        logpost = list(),
        N = list(),
        n_params = list(),
        MAP_A_counts = list(),
        BIC = list()
    )

    logs <- list(
        P = list(),
        E = list(),
        A = list(),
        q = list(),
        prob_inclusion = list()
    )
    if (likelihood == "normal") {
        logs$sigmasq <- list()
    } else if (likelihood == "poisson") {
        logs$Z <- list()
    }

    # initialize convergence status
    convergence_status <- list(
        prev_MAP_metric = Inf,
        best_MAP_metric = Inf,
        inarow_no_change = 0,
        inarow_no_best = 0,
        converged = FALSE,
        best_iter = NULL
    )

    # start logging
    sink(file = logfile)
    print(START)
    print(paste("maxiters =", convergence_control$maxiters))
    PREV = Sys.time()
    print(paste("starting iterations,", PREV))
    avg_time = 0

    # Gibbs sampler
    iter = 1
    logiter = 1
    done = FALSE
    first_MAP = TRUE
    stop = NULL
    START_ITER <- Sys.time()
    while (iter <= convergence_control$maxiters & !done) {

        # update P
        if (!Theta$is_fixed$P) {
            for (n in sample(1:dims$N)) {
                Theta$P[, n] <- sample_Pn(n, M, Theta, dims, likelihood = likelihood, prior = prior, gamma = gamma_sched[iter])
            }
        }

        # update E
        if (!Theta$is_fixed$E) {
            for (n in sample(1:dims$N)) {
                Theta$E[n, ] <- sample_En(n, M, Theta, dims, likelihood = likelihood, prior = prior, gamma = gamma_sched[iter])
            }
        }

        # if Normal likelihood, update sigmasq
        if (likelihood == 'normal') {
            if (!Theta$is_fixed$sigmasq) {
                Theta$sigmasq <- sample_sigmasq_normal(M, Theta, dims, gamma = gamma_sched[iter])
            }
        }

        # update A and q
        if (!Theta$is_fixed$A) {
            for (n in sample(1:dims$N)) {
                sample_An_out <- sample_An(n, M, Theta, dims, logfac, likelihood = likelihood,  prior = prior, gamma = gamma_sched[iter])
                Theta$A[1, n] <- sample_An_out$sampled
                Theta$prob_inclusion[1,n] <- sample_An_out$prob_inclusion
            }
        }
        if (!Theta$is_fixed$q) {
            for (n in sample(1:dims$N)) {
                Theta$q[1, n] <- sample_qn(n, Theta, gamma = gamma_sched[iter])
            }
        }

        # if Poisson likelihood, update latent counts Z
        if (likelihood == 'poisson') {
            for (k in sample(1:dims$K)) {
                for (g in sample(1:dims$G)) {
                    Theta$Z[k,,g] <- sample_Zkg_poisson(k, g, M, Theta, dims, gamma = gamma_sched[iter])
                }
            }
        }

        # update priors
        for (n in 1:dims$N) {
            if (prior == "truncnormal") {
                if (!Theta$is_fixed$prior_P[n]) {
                    Theta$Mu_p[,n] <- sample_Mu_Pn(n, Theta, dims, gamma = gamma_sched[iter])
                    Theta$Sigmasq_p[,n] <- sample_Sigmasq_Pn(n, Theta, dims, gamma = gamma_sched[iter])
                }
                Theta$Mu_e[n,] <- sample_Mu_En(n, Theta, dims, gamma = gamma_sched[iter])
                Theta$Sigmasq_e[n,] <- sample_Sigmasq_En(n, Theta, dims, gamma = gamma_sched[iter])
            } else if (prior == "exponential") {
                if (!Theta$is_fixed$prior_P[n]) {
                    Theta$Lambda_p[,n] <- sample_Lambda_Pn(n, Theta, dims, gamma = gamma_sched[iter])
                }
                Theta$Lambda_e[n,] <- sample_Lambda_En(n, Theta, dims, gamma = gamma_sched[iter])
            } else if (prior == "gamma") {
                if (!Theta$is_fixed$prior_P[n]) {
                    Theta$Beta_p[,n] <- sample_Beta_Pn(n, Theta, dims, gamma = gamma_sched[iter])
                    for (k in 1:dims$K) {
                        Theta$Alpha_p[k,n] <- sample_Alpha_Pkn(k, n, Theta, dims, gamma = gamma_sched[iter])
                    }
                }
                Theta$Beta_e[n,] <- sample_Beta_En(n, Theta, dims, gamma = gamma_sched[iter])
                for (g in 1:dims$G) {
                    Theta$Alpha_e[n,g] <- sample_Alpha_Eng(n, g, Theta, dims, gamma = gamma_sched[iter])
                }
            }
        }

        # log on original scale
        # only if storing logs or if we will use it for MAP
        if (store_logs | iter >= convergence_control$MAP_every + 1) {
            logs$P[[logiter]] <- Theta$P
            logs$E[[logiter]] <- Theta$E * rescale_by
            logs$A[[logiter]] <- Theta$A
            logs$q[[logiter]] <- Theta$q
            logs$prob_inclusion[[logiter]] <- Theta$prob_inclusion
            if (likelihood == "normal") {
                logs$sigmasq[[logiter]] <- Theta$sigmasq * (rescale_by**2)
            } else if (likelihood == "poisson") {
                logs$Z[[logiter]] <- Theta$Z * rescale_by
            }
            logiter = logiter + 1
        }

        # periodically check convergence and log progress
        if (
            (iter %% convergence_control$MAP_every == 0 &
             iter >= convergence_control$MAP_over + convergence_control$MAP_every)
            | iter == convergence_control$maxiters
        ) {
            # get MAP over past convergence_control$MAP_over iterations
            burn_in <- logiter - convergence_control$MAP_over
            keep <- burn_in:(logiter - 1)
            MAP <- get_MAP(logs, keep)

            # log metrics
            out <- update_metrics(
                metrics, MAP, iter, Theta, M_truescale, M,
                likelihood, prior, dims, logfac, rescale_by
            )
            metrics <- out$metrics
            Theta_MAP_rescaled <- out$Theta_MAP_rescaled

            # check convergence
            convergence_status <- check_converged(
                iter, gamma_sched[iter],
                Mhat_MAP_rescaled, M,
                convergence_status,
                convergence_control,
                first_MAP,
                Theta = Theta_MAP_rescaled,
                likelihood = likelihood,
                prior = prior,
                dims = dims,
                logfac = logfac
            )

            # check whether convergence is allowed (i.e., whether tempering is over)
            if (gamma_sched[iter] == 1 & first_MAP & MAP$top_counts[1] >= convergence_control$minA) {
                first_MAP = FALSE
                # forces convergence after gamma == 1
                convergence_status$best_MAP_metric = Inf
            }

            # if not storing logs, store current best MAP to return if needed
            # and remove logs to save memory
            if (
                !store_logs &
                !is.null(convergence_status$best_iter)
            ) {
                # update stored MAP if this is the best_iter
                if (convergence_status$best_iter == iter) {
                    store_MAP <- MAP
                    store_MAP_iter <- iter
                    store_credible_intervals <- get_credible_intervals(logs, store_MAP$idx)
                }

                # remove logs to save memory
                drop <- 1:convergence_control$MAP_every
                for (name in names(logs)) {
                    logs[[name]][drop] <- NULL
                }
                logiter = logiter - length(drop)
            }

            # log progress
            NOW = Sys.time()
            diff = as.numeric(difftime(NOW, PREV, units = "secs"))
            PREV = NOW
            log_MAP(
                iter, done, diff, convergence_control,
                convergence_status, gamma_sched, MAP, learn_A
            )

            # if converged, stop
            if (convergence_status$converged){
                stop = convergence_status$best_iter
                log_converged(convergence_control, convergence_status)
                done = TRUE

                # re-compute MAP at stop, compute 95% credible intervals
                if (store_logs) {
                    burn_in <- stop - convergence_control$MAP_over
                    keep <- burn_in:stop
                    MAP <- get_MAP(logs, keep, final = TRUE)
                    credible_intervals <- get_credible_intervals(logs, MAP$idx)
                } else {
                    MAP <- store_MAP
                    keep_sigs <- which(MAP$A[1,] == 1)
                    MAP$P <- MAP$P[, keep_sigs]
                    MAP$E <- MAP$E[keep_sigs, ]
                    credible_intervals <- store_credible_intervals
                }
            }

            # plot metrics
            plot_metrics(metrics, plotfile, stop, learn_A, gamma_sched, iter, true_P)

            # save results
            metrics_df = data.frame(matrix(nrow = length(unlist(metrics$BIC)), ncol = 0))
            for (name in names(metrics)) {
                metrics_df[,name] <- unlist(metrics[[name]])
            }
            keep_sigs = as.logical(MAP$A[1, ])
            res <- list(
                MAP = MAP,
                metrics = metrics_df,
                model = list(
                    M = M_truescale,
                    true_P = true_P,
                    likelihood = likelihood,
                    prior = prior,
                    prior_parameters = prior_parameters,
                    fixed = fixed,
                    inits = inits,
                    convergence_control = convergence_control,
                    dims = dims,
                    gamma_sched = gamma_sched
                ),
                totaliters = iter,
                converged_at = stop,
                final_Theta = Theta,
                time = list(
                    "avg_secs_per_iter" = as.numeric(difftime(Sys.time(), START_ITER, units = "secs"))/iter,
                    "total_secs" = as.numeric(difftime(Sys.time(), START, units = "secs"))
                )
            )
            if (done) {
                res$credible_intervals <- credible_intervals
            }
            if (store_logs) {
                res$logs = logs
            }
            saveRDS(res, file = savefile)
        }

        iter = iter + 1
    }

    # similarity matrix with true_P, if provided
    if (!is.null(true_P)) {
        sim_mat <- pairwise_sim(
            res$MAP$P, true_P,
            name1 = "estimated", name2 = "true",
            which = "cols"
        )
        heatmap <- get_heatmap(res$MAP$P, true_P)

        res$sim_mat <- sim_mat
        res$heatmap <- heatmap
        saveRDS(res, file = savefile)
    }

    # end log
    sink()
    return(res)
}
