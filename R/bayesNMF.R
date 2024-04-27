#' Perform single-study Bayesian NMF with the provided likelihood and prior
#' combination. Exact rank `N` or maximum rank `max_N` must be provided.
#'
#' @param M mutational catalog matrix, K x G
#' @param N fixed number of signatures
#' @param max_N maximum number of signatures if learning rank
#' @param inits (optional) list of initial values for P and E as well as sigmasq
#' if `likelihood = "normal"`
#' @param fixed (ptional) list of parameters to fix and not include in Gibbs
#' updates.
#' @param likelihood string, one of c('normal','poisson')
#' @param prior string, one of c('truncnormal','exponential')
#' @param sigmasq_type string, one of c('eq_mu','invgamma','noninformative')
#' @param prior_parameters list, optional specification of prior parameters
#' @param file file name without extension of log, save, and plot files
#' @param overwrite if `overwrite = TRUE`, the log, safe, and plot files of
#' previous runs with the same `file` will be overwritten
#' @param true_P (optional) true signatures matrix P to compare to in a heatmap
#' @param convergence_control list, specification of convergence parameters.
#' See documentation for `new_convergence_control`.
#' @param store_logs boolean, whether to store each sample in resulting .RData file
#'
#' @return list
#' @export
bayesNMF <- function(
        M,
        N = NULL,
        max_N = NULL,
        inits = NULL,
        fixed = NULL,
        likelihood = "normal",
        prior = "truncnormal",
        sigmasq_type = "noninformative",
        prior_parameters = NULL,
        file = paste0('nmf_', likelihood, '_', prior),
        overwrite = FALSE,
        true_P = NULL,
        convergence_control = new_convergence_control(),
        store_logs = FALSE,
        temper = FALSE
) {
    START = Sys.time()

    rescale_by = mean(M)/100
    M_truescale = M
    M = M/rescale_by
    if (!is.null(fixed$sigmasq)) {
        fixed$sigmasq = fixed$sigmasq/(rescale_by**2)
    }
    if (!is.null(fixed$E)) {
        fixed$E <- fixed$E/rescale_by
    }

    # check N/max_N combination is valid
    N <- validate_N(N, max_N)

    # check prior and likelihood are valid
    validate_model(likelihood, prior)

    # set up tempering schedule
    learn_A <- !is.null(max_N) & is.null(fixed$A)
    if (learn_A | temper) {
        gamma_sched <- get_gamma_sched(len = convergence_control$maxiters)
    } else {
        gamma_sched <- rep(1, convergence_control$maxiters)
    }

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
    savefile = paste0(file, '.RData')
    logfile = paste0(file, '.log')
    plotfile = paste0(file, '.pdf')
    tail = 0
    while (!overwrite & (file.exists(savefile) | file.exists(logfile))) {
        tail = tail + 1
        savefile = paste0(file, '_', tail, '.RData')
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
        sigmasq_type = sigmasq_type,
        inits = inits, fixed = fixed,
        prior_parameters = prior_parameters
    )

    # set up logs
    metrics <- list(
        sample_idx = list(),
        RMSE = list(),
        KL = list(),
        loglik = list(),
        logpost = list(),
        N = list(),
        n_params = list(),
        BIC = list()
    )

    logs <- list(
        P = list(),
        E = list(),
        A = list(),
        q = list()
    )
    if (likelihood == "normal") {
        logs$sigmasq <- list()
    } else if (likelihood == "poisson") {
        logs$Z <- list()
    }

    convergence_status <- list(
        prev_MAP_RMSE = 1000,
        best_MAP_RMSE = 1000,
        inarow_no_change = 0,
        inarow_no_best = 0,
        converged = FALSE
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
    done = FALSE
    first_MAP = TRUE
    stop = NULL
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
                if (sigmasq_type == 'eq_mu') {
                    Theta$sigmasq <- rowMeans(get_Mhat(Theta))
                } else {
                    Theta$sigmasq <- sample_sigmasq_normal(M, Theta, dims, sigmasq_type, gamma = gamma_sched[iter])
                }
            }
        }

        # if Poisson likelihood, update latent counts Z
        if (likelihood == 'poisson') {
            for (k in sample(1:dims$K)) {
                for (g in sample(1:dims$G)) {
                    Theta$Z[k,,g] <- sample_Zkg_poisson(k, g, M, Theta, dims)
                }
            }
        }

        # update A and q if learn_A = TRUE
        if (!Theta$is_fixed$A) {
            for (n in sample(1:dims$N)) {
                Theta$A[1, n] <- sample_An(n, M, Theta, dims, logfac, likelihood = likelihood,  prior = prior, gamma = gamma_sched[iter])
            }
        }

        if (!Theta$is_fixed$q) {
            for (n in sample(1:dims$N)) {
                Theta$q[1, n] <- sample_qn(n, Theta, gamma = gamma_sched[iter])
            }
        }

        # log on original scale
        logs$P[[iter]] <- Theta$P
        logs$E[[iter]] <- Theta$E * rescale_by
        logs$A[[iter]] <- Theta$A
        logs$q[[iter]] <- Theta$q

        if (likelihood == "normal") {
            logs$sigmasq[[iter]] <- Theta$sigmasq * (rescale_by**2)
        } else if (likelihood == "poisson") {
            logs$Z[[iter]] <- Theta$Z * rescale_by
        }

        # periodically check convergence and log progress
        if (
            (iter %% convergence_control$MAP_every == 0 &
             iter >= convergence_control$MAP_over + convergence_control$MAP_every)
            | iter == convergence_control$maxiters
        ) {
            # get MAP over past convergence_control$MAP_over iterations
            burn_in <- iter - convergence_control$MAP_over
            keep <- burn_in:iter
            MAP <- get_MAP(logs, keep)

            # log metrics
            out <- update_metrics(
                metrics, MAP, iter, Theta, M_truescale, M,
                likelihood, prior, dims, logfac, rescale_by,
                sigmasq_type
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
                logfac = logfac,
                sigmasq_eq_mu = sigmasq_type == 'eq_mu'
            )
            if (gamma_sched[iter] == 1 & first_MAP & MAP$top_counts[1] >= convergence_control$minA) {
                first_MAP = FALSE
                # forces convergence after gamma == 1
                convergence_status$best_MAP_metric = Inf
            }

            NOW = Sys.time()
            diff = as.numeric(difftime(NOW, PREV, units = "secs"))
            PREV = NOW
            log_MAP(iter, done, diff, convergence_control, convergence_status, gamma_sched)

            if (learn_A) {
                print(MAP$top_counts)
                cat("\n")
            }

            if (convergence_status$converged){
                stop = convergence_status$best_iter
                log_converged(convergence_control, convergence_status)
                done = TRUE

                # re-compute MAP at stop, compute 95% credible intervals
                burn_in <- stop - convergence_control$MAP_over
                keep <- burn_in:stop
                MAP <- get_MAP(logs, keep)
                credible_intervals <- get_credible_intervals(logs, MAP$idx)
            }

            # plot metrics
            plot_metrics(metrics, plotfile, stop, learn_A, gamma_sched, iter, true_P)

            # save results
            keep_sigs = as.logical(MAP$A[1, ])
            res <- list(
                M = M_truescale,
                true_P = true_P,
                MAP = MAP,
                metrics = metrics,
                convergence_control = convergence_control,
                totaliters = iter,
                converged_at = stop,
                final_Theta = Theta,
                dims = dims,
                time = list(
                    "avg_secs_per_iter" = as.numeric(difftime(Sys.time(), START, units = "secs"))/iter,
                    "total_secs" = as.numeric(difftime(Sys.time(), START, units = "secs"))
                )
            )
            if (done) {
                res$Credible_Intervals <- credible_intervals
            }
            if (store_logs) {
                res$logs = logs
            }
            save(res, file = savefile)
        }

        iter = iter + 1
    }

    # similarity matrix with true_P, if provided
    # only at end because it slows things down
    if (!is.null(true_P)) {
        sim_mat <- pairwise_sim(
            res$MAP$P, true_P,
            name1 = "estimated", name2 = "true",
            which = "cols"
        )
        heatmap <- get_heatmap(res$MAP$P, true_P)

        res$sim_mat <- sim_mat
        res$heatmap <- heatmap
        save(res, file = savefile)
    }

    # end log
    sink()
    return(res)
}

plot_one <- function(x, y, vblue = NULL, vgreen = NULL, xlab = "", ylab = "") {
    if (sum(!is.na(y)) > 0){
        if (sum(y != -Inf, na.rm = TRUE) > 0) {
            plot(x, y, ylab = ylab, xlab = xlab)
            if (!is.null(vblue)) {
                abline(v = vblue, col = 'blue')
            }
            if (!is.null(vgreen)) {
                abline(v = vgreen, col = 'green')
            }
        }
    }
}

plot_metrics <- function(metrics, plotfile, stop, learn_A, gamma_sched, iter, true_P) {

    if (!is.null(stop)) {
        vblue = stop
    } else { vblue = NULL }
    if (learn_A & gamma_sched[iter] == 1) {
        vgreen = which(gamma_sched == 1)[1]
    } else { vgreen = NULL }


    grDevices::pdf(plotfile)
    graphics::par(mfrow = c(3,1))

    x = unlist(metrics$sample_idx)
    plot_one(x, unlist(metrics$loglik), vblue, vgreen, xlab = "Iteration", ylab = "Log Likelihood")
    plot_one(x, unlist(metrics$logpost), vblue, vgreen, xlab = "Iteration", ylab = "Log Posterior")
    plot_one(x, unlist(metrics$BIC), vblue, vgreen, xlab = "Iteration", ylab = "BIC")
    plot_one(x, unlist(metrics$RMSE), vblue, vgreen, xlab = "Iteration", ylab = "RMSE")
    plot_one(x, unlist(metrics$KL), vblue, vgreen, xlab = "Iteration", ylab = "KL Divergence")

    plot_one(x, unlist(metrics$N), vblue, vgreen, xlab = "Iteration", ylab = "Number of Signatures")
    if (!is.null(true_P)) {
        abline(h = ncol(true_P))
    }

    grDevices::dev.off()
}

get_MAP <- function(logs, keep) {

    # get MAP of A matrix (fine to do even if learn_A = FALSE)
    A_MAP = get_mode(logs$A[keep])
    map.idx = keep[A_MAP$idx]

    # get MAP of P, E conditional on MAP of A
    MAP <- list(
        A = A_MAP$matrix,
        P = get_mean(logs$P[map.idx]),
        E = get_mean(logs$E[map.idx]),
        q = get_mean(logs$q[map.idx]),
        idx = map.idx,
        top_counts = A_MAP$top_counts
    )
    if ("sigmasq" %in% names(logs)) {
        MAP$sigmasq <- get_mean(logs$sigmasq[map.idx])
    }

    return(MAP)
}

update_metrics <- function(
    metrics, MAP, iter, Theta, M_truescale, M,
    likelihood, prior, dims, logfac, rescale_by,
    sigmasq_type
) {
    Theta_MAP <- Theta
    Theta_MAP$P = MAP$P
    Theta_MAP$E = MAP$E
    Theta_MAP$A = MAP$A
    Theta_MAP$q = MAP$q
    if (likelihood == 'normal') {
        Theta_MAP$sigmasq = MAP$sigmasq
    }

    Theta_MAP_rescaled <- Theta_MAP
    Theta_MAP_rescaled$E = MAP$E/rescale_by
    if (likelihood == 'normal') {
        Theta_MAP_rescaled$sigmasq = MAP$sigmasq/(rescale_by**2)
    }

    Mhat_MAP <- get_Mhat(Theta_MAP)
    Mhat_MAP_rescaled <- get_Mhat(Theta_MAP_rescaled)
    metrics$sample_idx[[iter]] <- iter
    metrics$RMSE[[iter]] <- get_RMSE(M_truescale, Mhat_MAP)
    metrics$KL[[iter]] <- get_KLDiv(M_truescale, Mhat_MAP)
    if (likelihood == 'normal') {
        metrics$loglik[[iter]] <- get_loglik_normal(M, Theta_MAP_rescaled, dims)
    } else if (likelihood == 'poisson') {
        metrics$loglik[[iter]] <- get_loglik_poisson(M, Theta_MAP_rescaled, dims, logfac)
    }
    metrics$N[[iter]] <- sum(Theta_MAP$A[1,])
    metrics$n_params[[iter]] <- metrics$N[[iter]] * (dims$G + dims$K + 2)
    metrics$BIC[[iter]] <- metrics$n_params[[iter]] * log(dims$G) - 2 * metrics$loglik[[iter]]
    metrics$logpost[[iter]] <- metrics$loglik[[iter]] + get_logprior(
        Theta_MAP_rescaled, likelihood, prior, sigmasq_type
    )

    return(list(
        metrics = metrics,
        Theta_MAP_rescaled = Theta_MAP_rescaled
    ))
}

log_MAP <- function(iter, done, diff, convergence_control, convergence_status, gamma_sched) {
    cat(paste(
        iter, "/", ifelse(done, "", "(up to)"), convergence_control$maxiters,
        "-", round(diff, 4), "seconds",
        ifelse(
            gamma_sched[iter] == 1,
            paste("-", paste0(round(convergence_status$prev_percent_change * 100, 4), "% change"),
                  "-", convergence_status$inarow_no_best, "no best",
                  "-", convergence_status$inarow_no_change, "no change"),
            ""
        ),
        "\n"
    ))
}

log_converged <- function(convergence_control, convergence_status) {
    cat(paste("\n\nCONVERGED at", convergence_status$best_iter))
    if (convergence_status$why %in% c("no best", "max iters")) {
        cat(paste(
            "\nNo best MAP since sample",
            convergence_status$best_iter, "\n\n"
        ))
    } else {
        cat(paste(
            "\nNo change in MAP over past",
            convergence_control$Ninarow_nochange, "samples\n\n"
        ))
    }
}

get_credible_intervals <- function(logs, map.idx) {
    credible_intervals <- list()
    credible_intervals[["P"]] <- get_quantile(logs$P[map.idx])
    credible_intervals[["E"]] <- get_quantile(logs$E[map.idx])
    credible_intervals[["q"]] <- get_quantile(logs$q[map.idx])
    if ("sigmasq" %in% names(logs)) {
        credible_intervals[["sigmasq"]] <- get_quantile(logs$sigmasq[map.idx])
    }
    return(credible_intervals)
}

validate_N <- function(N, max_N) {
    if (is.null(N) & is.null(max_N)) {
        stop("Either `N` or `max_N` must be provided.")
    } else if (!is.null(N) & !is.null(max_N)) {
        message("Both `N` and `max_N` provided, using `N`.")
        max_N = NULL
    } else if (is.null(N)) {
        N = max_N
    }
    return(N)
}

validate_model <- function(likelihood, prior) {
    if (!(likelihood %in% c('normal', 'poisson'))) {
        stop("likelihood must be one of c('normal')")
    } else if (likelihood == 'normal') {
        if (!(prior %in% c('truncnormal','exponential'))) {
            stop("prior must be one of c('truncnormal','exponential') with `likelihood = 'normal'`")
        }
    } else if (likelihood == 'poisson') {
        if (!(prior %in% c('gamma','exponential'))) {
            stop("prior must be one of c('gamma','exponential') with `likelihood = 'poisson'`")
        }
    }
}

