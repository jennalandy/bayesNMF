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
    if (is.null(N) & is.null(max_N)) {
        stop("Either `N` or `max_N` must be provided.")
    } else if (!is.null(N) & !is.null(max_N)) {
        message("Both `N` and `max_N` provided, using `N`.")
        max_N = NULL
    } else if (is.null(N)) {
        N = max_N
    }

    # check prior and likelihood are valid
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
    sample_idx <- c()
    RMSE <- c()
    KL <- c()
    loglik <- c()
    logpost <- c()

    P.log <- list()
    E.log <- list()
    sigmasq.log <- list()
    A.log <- list()
    q.log <- list()

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
            for (n in 1:dims$N) {
                Theta$P[, n] <- sample_Pn(n, M, Theta, dims, likelihood = likelihood, prior = prior, gamma = gamma_sched[iter])
            }
        }

        # update E
        if (!Theta$is_fixed$E) {
            for (n in 1:dims$N) {
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
            for (k in 1:dims$K) {
                for (g in 1:dims$G) {
                    Theta$Z[k,,g] <- sample_Zkg_poisson(k, g, M, Theta, dims)
                }
            }
        }

        # update A and q if learn_A = TRUE
        if (!Theta$is_fixed$A) {
            for (n in 1:dims$N) {
                Theta$A[1, n] <- sample_An(n, M, Theta, dims, logfac, likelihood = likelihood, gamma = gamma_sched[iter])
            }
        }

        if (!Theta$is_fixed$q) {
            for (n in 1:dims$N) {
                Theta$q[1, n] <- sample_qn(n, Theta, gamma = gamma_sched[iter])
            }
        }

        # log on original scale
        P.log[[iter]] <- Theta$P
        E.log[[iter]] <- Theta$E * rescale_by
        sigmasq.log[[iter]] <- Theta$sigmasq * (rescale_by**2)
        A.log[[iter]] <- Theta$A
        q.log[[iter]] <- Theta$q

        # periodically check convergence and log progress
        if (
            (iter %% convergence_control$MAP_every == 0 &
             iter >= convergence_control$MAP_over + convergence_control$MAP_every)
            | iter == convergence_control$maxiters
        ) {
            # get MAP over past convergence_control$MAP_over iterations
            burn_in <- iter - convergence_control$MAP_over
            keep <- burn_in:iter
            # get MAP of A matrix (fine to do even if learn_A = FALSE)
            A_MAP = get_mode(A.log[keep])
            map.idx = keep[A_MAP$idx]

            # only keep signatures present in MAP A matrix
            keep_sigs = as.logical(A_MAP$matrix[1, ])
            # get MAP of P, E conditional on MAP of A
            P_MAP <- get_mean(P.log[map.idx])
            E_MAP <- get_mean(E.log[map.idx])
            q_MAP <- get_mean(q.log[map.idx])
            sigmasq_MAP <- get_mean(sigmasq.log[map.idx])

            # log metrics
            Theta_MAP <- Theta
            Theta_MAP$P = P_MAP
            Theta_MAP$E = E_MAP
            Theta_MAP$A = A_MAP$matrix
            Theta_MAP$q = q_MAP
            Theta_MAP$sigmasq = sigmasq_MAP

            Theta_MAP_rescaled <- Theta_MAP
            Theta_MAP_rescaled$E = E_MAP/rescale_by
            Theta_MAP_rescaled$sigmasq = sigmasq_MAP/(rescale_by**2)

            Mhat_MAP <- get_Mhat(Theta_MAP)
            Mhat_MAP_rescaled <- get_Mhat(Theta_MAP_rescaled)
            sample_idx <- c(sample_idx, iter)
            RMSE <- c(RMSE, get_RMSE(M_truescale, Mhat_MAP))
            KL <- c(KL, get_KLDiv(M_truescale, Mhat_MAP))
            if (likelihood == 'normal') {
                this_loglik <- get_loglik_normal(M, Theta_MAP_rescaled, dims)
                loglik <- c(loglik, this_loglik)
            } else if (likelihood == 'poisson') {
                this_loglik <- get_loglik_poisson(M, Theta_MAP_rescaled, dims, logfac)
                loglik <- c(loglik, this_loglik)
            }
            logpost <- c(logpost, this_loglik + get_logprior(
                Theta_MAP_rescaled, likelihood, prior, sigmasq_type == 'eq_mu'
            ))
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
            first_MAP = FALSE

            NOW = Sys.time()
            diff = as.numeric(difftime(NOW, PREV, units = "secs"))
            PREV = NOW
            cat(paste(
                iter, "/", ifelse(done, "", "(up to)"), convergence_control$maxiters,
                "-", round(diff, 4), "seconds",
                "-", paste0(round(convergence_status$prev_percent_change * 100, 4), "% change"),
                ifelse(
                    gamma_sched[iter] == 1,
                    paste("-", convergence_status$inarow_no_best, "no best",
                         "-", convergence_status$inarow_no_change, "no change"),
                    ""
                ),
                "\n"
            ))
            if (learn_A) {
                print(A_MAP$top_counts)
                cat("\n")
            }

            if (convergence_status$converged){
                cat(paste("\n\nCONVERGED at", convergence_status$best_iter))
                if (convergence_status$why %in% c("no best", "max iters")) {
                    cat(paste("\nNo best MAP since sample", convergence_status$best_iter, "\n\n"))
                    stop = convergence_status$best_iter
                } else {
                    cat(paste(
                        "\nNo change in MAP over past",
                        convergence_control$Ninarow_nochange,
                        "samples\n\n"
                    ))
                    stop = convergence_status$best_iter
                }
                done = TRUE

                #re-compute MAP at stop
                # get MAP over past convergence_control$MAP_over iterations
                burn_in <- stop - convergence_control$MAP_over
                keep <- burn_in:stop
                # get MAP of A matrix (fine to do even if learn_A = FALSE)
                A_MAP = get_mode(A.log[keep])
                map.idx = keep[A_MAP$idx]

                # only keep signatures present in MAP A matrix
                keep_sigs = as.logical(A_MAP$matrix[1, ])
                # get MAP of P, E conditional on MAP of A
                P_MAP <- get_mean(P.log[map.idx])
                E_MAP <- get_mean(E.log[map.idx])
                q_MAP <- get_mean(q.log[map.idx])
                sigmasq_MAP <- get_mean(sigmasq.log[map.idx])
            }

            # plot metrics
            grDevices::pdf(plotfile)
            graphics::par(mfrow = c(3,1))
            plot(sample_idx, RMSE)
            if (!is.null(stop)) {
                abline(v = stop, col = 'blue')
            }
            if (learn_A & gamma_sched[iter] == 1) {
                abline(v = which(gamma_sched == 1)[1], col = 'green')
            }
            if (sum(!is.na(KL)) > 0) {
                if (sum(KL != -Inf, na.rm = TRUE) > 0) {
                    plot(sample_idx, KL)
                    if (!is.null(stop)) {
                        abline(v = stop, col = 'blue')
                    }
                    if (learn_A & gamma_sched[iter] == 1) {
                        abline(v = which(gamma_sched == 1)[1], col = 'green')
                    }
                }
            }
            if (sum(!is.na(loglik)) > 0){
                if (sum(loglik != -Inf, na.rm = TRUE) > 0) {
                    plot(sample_idx, loglik)
                    if (!is.null(stop)) {
                        abline(v = stop, col = 'blue')
                    }
                    if (learn_A & gamma_sched[iter] == 1) {
                        abline(v = which(gamma_sched == 1)[1], col = 'green')
                    }
                }
            }
            if (sum(!is.na(logpost)) > 0){
                if (sum(logpost != -Inf, na.rm = TRUE) > 0) {
                    plot(sample_idx, logpost)
                    if (!is.null(stop)) {
                        abline(v = stop, col = 'blue')
                    }
                    if (learn_A & gamma_sched[iter] == 1) {
                        abline(v = which(gamma_sched == 1)[1], col = 'green')
                    }
                }
            }
            grDevices::dev.off()

            # save results
            res <- list(
                M = M_truescale,
                true_P = true_P,
                MAP = list(
                    A = A_MAP$matrix,
                    q = q_MAP,
                    P = P_MAP[,keep_sigs],
                    E = E_MAP[keep_sigs,],
                    sigmasq = sigmasq_MAP
                ),
                metrics = list(
                    sample_idx = sample_idx,
                    loglik = loglik,
                    logpost = logpost,
                    RMSE = RMSE,
                    KL = KL
                ),
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
            if (store_logs) {
                res$logs = list(
                    P = P.log,
                    E = E.log,
                    sigmasq = sigmasq.log,
                    q = q.log,
                    A = A.log
                )
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
