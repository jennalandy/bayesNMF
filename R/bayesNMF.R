#' Perform single-study Bayesian NMF with the provided likelihood and prior
#' combination. Exact rank `N` or maximum rank `max_N` must be provided.
#'
#' @param M mutational catalog matrix, K x G
#' @param N number of signatures
#' @param P (optional) initial signatures matrix, K x N
#' @param E (optional) initial exposure matrix, N x G
#' @param sigmasq (optinal) initial variance vector, length K
#' @param prior string, one of c('truncnormal','exponential')
#' @param niters how many iterations the Gibbs sampler is run
#' @param burn_in the first `burn_in` iterations will be ignored when computing
#' MAP estimate
#' @param logevery the log, save, and plot files will be updated every
#' `logevery` iterations
#' @param file file name without extension of log, save, and plot files
#' @param overwrite if `overwrite = TRUE`, the log, safe, and plot files of
#' previous runs with the same `file` will be overwritten
#' @param true_P (optional) true signatures matrix P to compare to in a heatmap
#'
#' @return list
#' @export
bayesNMF <- function(
        M,
        N = NULL,
        max_N = NULL,
        inits = NULL,
        likelihood = "normal",
        prior = "truncnormal",
        prior_parameters = NULL,
        logevery = 100,
        file = paste0('nmf_', likelihood, '_', prior),
        overwrite = FALSE,
        true_P = NULL,
        niters = NULL,
        burn_in = NULL
) {
    START = Sys.time()

    # check N/max_N combination is valid
    if (is.null(N) & is.null(max_N)) {
        stop("Either `N` or `max_N` must be provided.")
    } else if (!is.null(N) & !is.null(max_N)) {
        message("Both `N` and `max_N` provided, using `N`.")
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

    # check whether A is given or learned
    learn_A <- !is.null(max_N) & is.null(inits$A)

    # number of iterations depends on whether we're learning A and likelihood
    if (is.null(niters)) {
        if (!learn_A) {
            niters = 1500
            burn_in = 1000
        } else if (likelihood == 'normal') {
            niters = 5000
            burn_in = 3500
        } else if (likelihood == 'poisson') {
            niters = 10000
            burn_in = 7500
        }
    }

    # set up tempering schedule
    if (learn_A) {
        gamma_sched <- get_gamma_sched(len = niters)
    } else {
        gamma_sched <- rep(1, niters)
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

    # check burn_in is valid
    if (burn_in > niters) {
        message(paste0(
            "Burn in ", burn_in, " is greater than niters ",
            niters, ", setting burn_in = 0"
        ))
        burn_in = 0
    }

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
        likelihood, prior,
        learn_A, dims,
        inits = inits,
        prior_parameters = prior_parameters
    )

    # set up logs
    RMSE <- c()
    KL <- c()
    loglik <- c()
    logpost <- c()

    P.log <- list()
    E.log <- list()
    sigmasq.log <- list()
    A.log <- list()
    q.log <- list()

    # start logging
    sink(file = logfile)
    print(START)
    print(paste("niters =", niters, "| burn_in =", burn_in))
    PREV = Sys.time()
    print(paste("starting iterations,", PREV))
    avg_time = 0

    # Gibbs sampler: sample niters times
    for (iter in 1:niters) {

        # update P
        for (n in 1:dims$N) {
            Theta$P[, n] <- sample_Pn(n, M, Theta, dims, likelihood = likelihood, prior = prior, gamma = gamma_sched[iter])
        }

        # update E
        for (n in 1:dims$N) {
            Theta$E[n, ] <- sample_En(n, M, Theta, dims, likelihood = likelihood, prior = prior, gamma = gamma_sched[iter])
        }

        # if Normal likelihood, update sigmasq
        if (likelihood == 'normal') {
            Theta$sigmasq <- sample_sigmasq_normal(M, Theta, dims, gamma = gamma_sched[iter])
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
        if (learn_A) {
            for (n in 1:dims$N) {
                Theta$A[1, n] <- sample_An(n, M, Theta, dims, logfac, likelihood = likelihood, gamma = gamma_sched[iter])
                Theta$q[1, n] <- sample_qn(n, Theta, gamma = gamma_sched[iter])
            }
        }

        # record in logs
        Mhat <- Theta$P %*% diag(Theta$A[1, ]) %*% Theta$E
        RMSE <- c(RMSE, get_RMSE(M, Mhat))
        KL <- c(KL, get_KLDiv(M, Mhat))
        if (likelihood == 'normal') {
            loglik <- c(loglik, get_loglik_normal(M, Theta, dims))
            logpost <- c(logpost, get_proportional_log_posterior(
                Theta, M, Theta$P, Theta$E,
                sigmasq = Theta$sigmasq,
                likelihood = likelihood, prior = prior
            ))
        } else if (likelihood == 'poisson') {
            loglik <- c(loglik, get_loglik_poisson(M, Theta, dims, logfac))
            logpost <- c(logpost, get_proportional_log_posterior(
                Theta, M, Theta$P, Theta$E,
                likelihood = likelihood, prior = prior
            ))
        }

        P.log[[iter]] <- Theta$P
        E.log[[iter]] <- Theta$E
        sigmasq.log[[iter]] <- Theta$sigmasq
        A.log[[iter]] <- Theta$A
        q.log[[iter]] <- Theta$q

        # periodically save every logevery iters and at the end
        if (iter %% logevery == 0 | iter == niters) {
            NOW = Sys.time()
            diff = as.numeric(difftime(NOW, PREV, units = "secs"))
            avg_time = (avg_time * (iter - logevery) + diff)/iter
            PREV = NOW
            cat(paste(iter, "/", niters, "-", diff, "seconds", "\n"))

            # plot metrics
            grDevices::pdf(plotfile)
            graphics::par(mfrow = c(3,1))
            plot(RMSE)
            if (sum(!is.na(KL)) > 0) {
                if (sum(KL != -Inf, na.rm = TRUE) > 0) {
                    plot(KL)
                }
            }
            if (sum(!is.na(loglik)) > 0){
                if (sum(loglik != -Inf, na.rm = TRUE) > 0) {
                    plot(loglik)
                }
            }
            if (sum(!is.na(logpost)) > 0) {
                if (sum(logpost != -Inf, na.rm = TRUE) > 0) {
                    plot(logpost)
                }
            }
            grDevices::dev.off()

            # only use burn_in if iter > burn_in
            if (iter > burn_in) {
                keep = burn_in:iter
            } else {
                keep = 1:iter
            }

            # get MAP of A matrix (fine to do even if learn_A = FALSE)
            A.map = get_mode(A.log[keep])
            map.idx = keep[A.map$idx]
            if (learn_A) {
                print(A.map$top_counts)
            }

            # only keep signatures present in MAP A matrix
            keep_sigs = as.logical(A.map$matrix[1, ])

            # save results
            res <- list(
                M = M,
                true_P = true_P,
                logs = list(
                    P = P.log,
                    E = E.log,
                    sigmasq = sigmasq.log,
                    q = q.log,
                    A = A.log
                ),
                MAP = list(
                    A = A.map$matrix,
                    q = get_mean(q.log[map.idx]),
                    P = get_mean(P.log[map.idx])[, keep_sigs],
                    E = get_mean(E.log[map.idx])[keep_sigs, ],
                    sigmasq = get_mean(sigmasq.log[map.idx])
                ),
                metrics = list(
                    loglik = loglik,
                    logpost = logpost,
                    RMSE = RMSE,
                    KL = KL
                ),
                burn_in = burn_in,
                niters = niters,
                final_Theta = Theta,
                dims = dims,
                time = list(
                    "avg_secs_per_iter" = avg_time,
                    "total_secs" = as.numeric(difftime(NOW, START, units = "secs"))
                )
            )
            save(res, file = savefile)
        }
    }

    # similarity matrix with true_P, if provided
    # only at end because it slows things down
    if (!is.null(true_P) & dims$N > 1) {
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
