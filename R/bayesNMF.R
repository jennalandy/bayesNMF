#' Bayesian Non-Negative Matrix Factorization
#' @description Bayesian Non-Negative Matrix Factorization.
#'
#' @param M matrix, data with samples as columns and features as rows
#' @param rank integer or vector, number of latent factors, integers
#' (e.g., \code{rank = 5}) or vector (e.g., \code{rank = 1:5}) are accepted.
#' If a vector is provided, values must be sequential and start at 0 or 1
#' (e.g., \code{1:5} or \code{0:5}) for \code{learn_rank_method = "BFI"}, but can be
#' nonsequential (e.g., \code{c(2,4,5)}) for \code{learn_rank_method = "heuristic"}.
#' @param learn_rank_method string, method used to learn latent rank, used if
#' a vector is provided for \code{rank}, one of c("BFI","heuristic") is used.
#' Bayesian Factor Inclusion (BFI) learns rank automatically as a part of the
#' Bayesian model with a binary Factor inclusion matrix. The heuristic approach
#' runs a fixed-rank Bayesian NMF for each value in \code{rank}, and selects
#' the model that minimizes BIC.
#' @param likelihood string, one of \code{c('normal','poisson')}, represents the
#' distribution used for likelihood f(M|P, E).
#' @param prior string, one of \code{c('truncnormal','exponential','gamma')}, represents
#' the distribution used for priors on P and E.
#' @param fast boolean, if \code{bast = TRUE}, \code{likelihood = 'poisson'},
#' and \code{fast = TRUE}, then updates from the corresponding
#' \code{likelihood = 'normal'} model are used as proposals
#' in an efficient Gibb's sampler. Only available for \code{likelihood = 'poisson'}
#' and \code{prior = c('truncnormal', 'exponential')}. Defaults \code{TRUE} when possible.
#' @param inits (optional) list, initial values for P and E (may also provide
#' sigmasq if \code{likelihood = "normal"} and/or A if \code{rank} is a vector).
#' @param fixed (optional) list, parameters values to fix at constant
#' values, rather than learning them through Gibbs updates.
#' @param clip numeric, if \code{rank} is a vector and \code{learn_rank_method = "BFI"},
#'  prior probabilities of factor inclusion will be clipped by
#' \code{clip}/N away from 0 and 1.
#' @param prior_parameters (optional) list, specification of prior parameters.
#' @param recovery boolean, if TRUE, additional factors are included with priors
#' at previously discovered factors. In this case, \code{rank} denotes the number of
#' additional latent factors on top of those with priors in \code{recovery_priors}.
#' @param recovery_priors "cosmic" or list, prior parameters for recovered
#' latent factors. If \code{recovery_priors = "cosmic"}, pre-computed priors
#' based on the 79 COSMIC signatures are used.
#' @param file string, file name (without extension) used for log, rds,
#' and pdf files created by this function.
#' @param true_P (optional) matrix, reference latent factors matrix P to
#' compare estimated factors to with a heatmap.
#' @param convergence_control list, specification of convergence parameters.
#' See documentation for \code{new_convergence_control}.
#' @param store_logs boolean, if \code{store_logs = TRUE}, each iteration of the
#' Gibb's sampler is stored in resulting \code{.rds} file. Otherwise, only MAP is
#' saved.
#' @param overwrite boolean, if \code{overwrite = TRUE}, the log, safe, and plot files of
#' previous runs with the same \code{file} will be overwritten
#'
#' @return list
#' @export
bayesNMF <- function(
    M, rank,
    learn_rank_method = "BFI",
    sparse_rank = FALSE,
    likelihood = "poisson",
    prior = "truncnormal",
    fast = (likelihood == "poisson" &
                prior %in% c('truncnormal','exponential')),
    inits = NULL,
    fixed = NULL,
    clip = 0.4,
    prior_parameters = NULL,
    recovery = FALSE,
    recovery_priors = "cosmic",
    file = paste0('nmf_', likelihood, '_', prior),
    true_P = NULL,
    convergence_control = new_convergence_control(),
    store_logs = TRUE,
    overwrite = FALSE
) {
    # set up file names
    final_file = file
    savefile = paste0(file, '.rds')
    logfile = paste0(file, '.log')
    tail = 0
    while (!overwrite & (file.exists(savefile) | file.exists(logfile))) {
        tail = tail + 1
        final_file = paste0(file, '_', tail)
        savefile = paste0(file, '_', tail, '.rds')
        logfile = paste0(file, '_', tail, '.log')
    }

    # if learning rank, consider "learn_rank_method"
    learn_A <- length(rank) > 1 & is.null(fixed$A)
    if (learn_A) {
        if (learn_rank_method == "BFI") {
            max_N = max(rank)
            min_N = min(rank)
            if (!setequal(min_N:max_N, rank)) {
                stop(paste(
                    "rank =", paste0("c(", paste(rank, collapse = ','), ")"),
                    "is not sequential. rank must be sequential",
                    "to use learn_rank_method = 'BFI'. Try a sequential rank",
                    "or learn_rank_method = 'heuristic'."
                ))
            }
            if (min_N > 1) {
                stop(paste(
                    "Minimum rank > 1 is not permitted with learn_rank_method =",
                    "'BFI'. Try rank =", paste0(1, ":", max_N),
                    "or use learn_rank_method = 'heuristic'."
                ))
            }
            tryCatch({
                sink(file = logfile)
                res <- inner_bayesNMF(
                    M = M,
                    N = NULL,
                    max_N = max_N,
                    sparse_rank = sparse_rank,
                    likelihood = likelihood,
                    prior = prior,
                    fast = fast,
                    inits = inits,
                    fixed = fixed,
                    clip = clip,
                    prior_parameters = prior_parameters,
                    recovery = recovery,
                    recovery_priors = recovery_priors,
                    file = final_file,
                    true_P = true_P,
                    convergence_control = convergence_control,
                    store_logs = store_logs,
                    overwrite = overwrite
                )
                sink()
            }, error = function(e) {
                sink()
                stop(e)
            }, interrupt = function(e) {
                sink()
            })
        } else if (learn_rank_method == "heuristic") {
            BICs <- data.frame(
                rank = rank,
                BIC = rep(NA, length(rank))
            )
            all_models <- list()
            for (r in rank) {
                tryCatch({
                    sink(file = paste0(final_file, '_rank', r, '.log'))
                    res_N <- inner_bayesNMF(
                        M = M,
                        N = r,
                        likelihood = likelihood,
                        prior = prior,
                        fast = fast,
                        inits = inits,
                        fixed = fixed,
                        clip = clip,
                        prior_parameters = prior_parameters,
                        recovery = recovery,
                        recovery_priors = recovery_priors,
                        file = paste0(final_file, "_rank", r),
                        true_P = true_P,
                        convergence_control = convergence_control,
                        store_logs = store_logs,
                        overwrite = overwrite
                    )
                    sink()
                }, error = function(e) {
                    sink()
                    stop(e)
                }, interrupt = function(e) {
                    sink()
                })
                BICs$BIC[BICs$rank == r] <- res_N$metrics$BIC[
                    res_N$metrics$sample_idx == res_N$converged_at
                ]
                all_models[[r]] <- res_N
            }
            best_rank <- BICs %>%
                dplyr::filter(BIC == min(BIC)) %>%
                dplyr::pull(rank)
            plot <- BICs %>%
                ggplot2::ggplot(ggplot2::aes(x = rank, y = BIC)) +
                ggplot2::geom_point() +
                ggplot2::geom_line()
            res <- list(
                best_rank = best_rank,
                best_model = all_models[[best_rank]],
                BICs = BICs,
                all_models = all_models,
                plot = plot
            )
            saveRDS(res, file = savefile)
        } else {
            stop(paste("learn_rank_method =", paste0("'", learn_rank_method, "'"),
                       "not defined,", "must be one of c('BFI','heuristic')."))
        }
    } else {
        tryCatch({
            sink(file = logfile)
            res <- inner_bayesNMF(
                M = M,
                N = rank,
                max_N = NULL,
                likelihood = likelihood,
                prior = prior,
                fast = fast,
                inits = inits,
                fixed = fixed,
                clip = clip,
                prior_parameters = prior_parameters,
                recovery = recovery,
                recovery_priors = recovery_priors,
                file = final_file,
                true_P = true_P,
                convergence_control = convergence_control,
                store_logs = store_logs,
                overwrite = overwrite
            )
            sink()
        }, error = function(e) {
            sink()
            stop(e)
        })
    }
    return(res)
}


#' Bayesian Non-Negative Matrix Factorization
#' @description Perform single-study Bayesian NMF with the provided likelihood and prior
#' combination. Exact rank `N` or maximum rank `max_N` must be provided.
#'
#' @param M mutational catalog matrix, K x G
#' @param N fixed number of latent factors
#' @param max_N maximum number of latent factors if learning rank
#' @param likelihood string, one of c('normal','poisson')
#' @param prior string, one of c('truncnormal','exponential','gamma')
#' @param fast boolean, if `likelihood == 'poisson'` and `fast = TRUE`, updates
#' from the corresponding `likelihood == 'normal'` model are used as proposals
#' in an efficient Gibb's sampler. Defaults TRUE when possible.
#' @param inits (optional) list of initial values for P and E (and sigmasq
#' if `likelihood = "normal"`)
#' @param fixed (optional) list of parameters to fix and not include in Gibbs
#' updates.
#' @param clip numeric, prior probabilities of inclusion will be clipped by
#' `clip`/N away from 0 and 1
#' @param prior_parameters list, optional specification of prior parameters
#' @param recovery boolean, whether to set priors of a subset of factors at
#' previously discovered factors
#' @param recovery_priors "cosmic" or list of prior parameters
#' @param file file name without extension of log, save, and plot files
#' @param true_P (optional) true latent factors matrix P to compare to in a heatmap
#' @param convergence_control list, specification of convergence parameters.
#' See documentation for `new_convergence_control`.
#' @param store_logs boolean, whether to store each sample in resulting .rds file
#' @param overwrite if `overwrite = TRUE`, the log, safe, and plot files of
#' previous runs with the same `file` will be overwritten
#'
#' @return list
#' @noRd
inner_bayesNMF <- function(
        M,
        N = NULL,
        max_N = NULL,
        sparse_rank = FALSE,
        likelihood = "poisson",
        prior = "truncnormal",
        fast = (likelihood == "poisson" &
                prior %in% c('truncnormal','exponential')),
        inits = NULL,
        fixed = NULL,
        clip = 0.4,
        prior_parameters = NULL,
        recovery = FALSE,
        recovery_priors = "cosmic",
        file = paste0('nmf_', likelihood, '_', prior),
        true_P = NULL,
        convergence_control = new_convergence_control(
            maxiters = ifelse(recovery, 5000, 2000)
        ),
        store_logs = TRUE,
        overwrite = FALSE
) {
    savefile = paste0(file, '.rds')
    plotfile = paste0(file, '.pdf')

    START = Sys.time()
    print(START)
    print(paste("maxiters =", convergence_control$maxiters))
    print(paste("fast", fast))

    # check recovery/discovery
    if (recovery) {
        if (is.character(recovery_priors)) {
            if (recovery_priors == "cosmic") {
                if (likelihood == 'normal' & prior == 'truncnormal') {
                    recovery_priors <- normal_truncnormal_recovery_priors
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
    if (recovery) {
        scale_by <- sqrt(mean(M)/N) / mean(recovery_priors$Mu_p)
        recovery_priors$Mu_p <- recovery_priors$Mu_p * scale_by
        recovery_priors$Sigmasq_p <- recovery_priors$Sigmasq_p * (scale_by**2)
    }

    # check prior and likelihood are valid
    validate_model(likelihood, prior, fast)

    # set up tempering schedule
    learn_A <- !is.null(max_N) & is.null(fixed$A)
    if (learn_A) {
        gamma_sched <- get_gamma_sched(len = convergence_control$maxiters)
    } else {
        gamma_sched <- rep(1, convergence_control$maxiters)
    }
    print(paste("learn_A",learn_A))

    # precompute log factorials for Poisson likelihood
    logfac = vector(length = max(M))
    logfac[1] = 0
    for (i in 2:length(logfac)) {
        logfac[i] = log(i) + logfac[i-1]
    }

    # set up dimensions
    dims = list(
        K = dim(M)[1],
        G = dim(M)[2],
        N = N,
        S = 1
    )

    # set up Theta
    Theta <- initialize_Theta(
        M = M,
        likelihood = likelihood,
        prior = prior,
        fast = fast,
        learn_A = learn_A,
        dims = dims,
        inits = inits, fixed = fixed,
        prior_parameters = prior_parameters,
        recovery = recovery,
        recovery_priors = recovery_priors,
        clip = clip
    )
    Theta$prob_inclusion <- Theta$A
    Theta$P_acceptance <- Theta$P
    Theta$E_acceptance <- Theta$E

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
        P_acceptance = list(),
        E_acceptance = list(),
        A = list(),
        prob_inclusion = list(),
        n = list()
    )
    if (likelihood == "normal" | (likelihood == "poisson" & fast)) {
        logs$sigmasq <- list()
    } else {
        # likelihood == "poisson" & !fast
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
                sample_Pn_out <- sample_Pn(
                    n, M, Theta, dims,
                    likelihood, prior, fast
                )
                Theta$P[, n] <- sample_Pn_out$sampled
                Theta$P_acceptance[, n] <- sample_Pn_out$acceptance
            }
        }

        # update E
        if (!Theta$is_fixed$E) {
            for (n in sample(1:dims$N)) {
                sample_En_out <- sample_En(
                    n, M, Theta, dims,
                    likelihood, prior, fast
                )
                Theta$E[n, ] <- sample_En_out$sampled
                Theta$E_acceptance[n, ] <- sample_En_out$acceptance
            }
        }

        # if Normal likelihood, update sigmasq
        if (likelihood == 'normal' | (likelihood == 'poisson' & fast)) {
            if (!Theta$is_fixed$sigmasq) {
                Theta$sigmasq <- sample_sigmasq_normal(
                    M, Theta, dims, gamma = 1
                )
            }
        }

        # update A and n
        if (!Theta$is_fixed$A) {
            Theta$n <- sample_n(Theta, dims, clip)
            Theta$q <- update_q(Theta, dims, clip)
            for (n in sample(1:dims$N)) {
                sample_An_out <- sample_An(
                    n, M, Theta, dims,
                    likelihood, prior, logfac,
                    sparse_rank = sparse_rank,
                    gamma = gamma_sched[iter]
                )
                Theta$A[1, n] <- sample_An_out$sampled
                Theta$prob_inclusion[1,n] <- sample_An_out$prob_inclusion
            }
        }

        # if Poisson likelihood, update latent counts Z
        if (likelihood == 'poisson' & !fast) {
            for (k in sample(1:dims$K)) {
                for (g in sample(1:dims$G)) {
                    Theta$Z[k,,g] <- sample_Zkg_poisson(
                        k, g, M, Theta, dims,
                        gamma = 1
                    )
                }
            }
        }

        # update priors
        for (n in 1:dims$N) {
            if (prior == "truncnormal") {
                if (!Theta$is_fixed$prior_P[n]) {
                    Theta$Mu_p[,n] <- sample_Mu_Pn(
                        n, Theta, dims, gamma = 1
                    )
                    Theta$Sigmasq_p[,n] <- sample_Sigmasq_Pn(
                        n, Theta, dims, gamma = 1
                    )
                }
                Theta$Mu_e[n,] <- sample_Mu_En(
                    n, Theta, dims, gamma = 1
                )
                Theta$Sigmasq_e[n,] <- sample_Sigmasq_En(
                    n, Theta, dims, gamma = 1
                )
            } else if (prior == "exponential") {
                if (!Theta$is_fixed$prior_P[n]) {
                    Theta$Lambda_p[,n] <- sample_Lambda_Pn(
                        n, Theta, dims, gamma = gamma_sched[iter]
                    )
                }
                Theta$Lambda_e[n,] <- sample_Lambda_En(
                    n, Theta, dims, gamma = gamma_sched[iter]
                )
            } else if (prior == "gamma") {
                if (!Theta$is_fixed$prior_P[n]) {
                    Theta$Beta_p[,n] <- sample_Beta_Pn(
                        n, Theta, dims, gamma = gamma_sched[iter]
                    )
                    for (k in 1:dims$K) {
                        Theta$Alpha_p[k,n] <- sample_Alpha_Pkn(
                            k, n, Theta, dims, gamma = gamma_sched[iter]
                        )
                    }
                }
                Theta$Beta_e[n,] <- sample_Beta_En(
                    n, Theta, dims, gamma = gamma_sched[iter]
                )
                for (g in 1:dims$G) {
                    Theta$Alpha_e[n,g] <- sample_Alpha_Eng(
                        n, g, Theta, dims, gamma = gamma_sched[iter]
                    )
                }
            }
        }

        # log on original scale
        # only if storing logs or if we will use it for MAP
        if (store_logs | iter >= convergence_control$MAP_every + 1) {
            logs$P[[logiter]] <- Theta$P
            logs$E[[logiter]] <- Theta$E
            logs$P_acceptance[[logiter]] <- Theta$P_acceptance
            logs$E_acceptance[[logiter]] <- Theta$E_acceptance
            logs$A[[logiter]] <- Theta$A
            logs$n[[logiter]] <- Theta$n
            logs$prob_inclusion[[logiter]] <- Theta$prob_inclusion
            if (likelihood == "normal" | (likelihood == "poisson" & fast)) {
                logs$sigmasq[[logiter]] <- Theta$sigmasq
            } else {
                # likelihood == "poisson" & !fast
                logs$Z[[logiter]] <- Theta$Z
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
                metrics, MAP, iter, Theta, M,
                likelihood, prior, dims, logfac
            )
            metrics <- out$metrics
            Theta_MAP <- out$Theta_MAP
            Mhat_MAP <- out$Mhat_MAP

            # check convergence
            convergence_status <- check_converged(
                iter, gamma_sched[iter],
                Mhat_MAP, M,
                convergence_status,
                convergence_control,
                first_MAP,
                Theta = Theta_MAP,
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
                # if (store_logs) {
                #     keep = (1:length(logs$A))[gamma_sched[1:length(logs$A)] == 1]
                #     running_posterior_counts <- get_posterior_counts_N(logs$A[keep])
                # } else {
                #     running_posterior_counts <- get_posterior_counts_N(logs$A)
                # }
            } else if (gamma_sched[iter] == 1) {
                # if (store_logs) {
                #     keep = (1:length(logs$A))[gamma_sched[1:length(logs$A)] == 1]
                #     running_posterior_counts <- running_posterior_counts + get_posterior_counts_N(logs$A[keep])
                # } else {
                #     running_posterior_counts <- running_posterior_counts + get_posterior_counts_N(logs$A)
                # }
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
                    if (dims$N > 1) {
                        MAP$P <- MAP$P[, keep_sigs]
                        MAP$E <- MAP$E[keep_sigs, ]
                    } else {
                        if (length(keep_sigs) == 0) {
                            MAP$P <- matrix(0, nrow = dims$K, ncol = 1)
                            MAP$E <- matrix(0, nrow = 1, ncol = dims$G)
                        }
                    }
                    credible_intervals <- store_credible_intervals
                }
                # posterior_pmf_N <- running_posterior_counts/sum(running_posterior_counts)
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
                    M = M,
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
                # res$posterior_N <- posterior_pmf_N
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

    return(res)
}
