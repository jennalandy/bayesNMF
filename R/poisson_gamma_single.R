#' sample P[k,n]
#'
#' @param k mutation type index
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#'
#' @return vector length K
sample_Pkn_poisson_gamma_single <- function(k, n, M, Theta) {
    rgamma(1, Theta$Alpha_p[k,n] + sum(Theta$Z[k,n,]), Theta$Beta_p[k,n] + sum(Theta$E[n,]))
}

#' Sample E[n,g]
#'
#' @param n signature index
#' @param g tumor genome index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#'
#' @return vector of length G
sample_Eng_poisson_gamma_single <- function(n, g, M, Theta) {
    rgamma(1, Theta$Alpha_e[n,g] + sum(Theta$Z[,n,g]), Theta$Beta_e[n,g] + sum(Theta$P[,n]))
}

#' sample Z[k,,g]
#'
#' @param k mutation type index
#' @param g tumor genome index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#'
#' @return vector length K
sample_Zkg_poisson_gamma_single <- function(k, g, M, Theta, dims){
    probs = sapply(1:dims$N, function(n) {
        Theta$P[k,n] * Theta$E[n,g]
    })
    probs = probs/sum(probs)
    rmultinom(1, size = M[k,g], prob = probs)
}

#' get log likelihood
#'
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#'
#' @return scalar
get_loglik_poisson_gamma_single <- function(M, Theta, dims, logfac) {
    - sum((Theta$P %*% Theta$E)) +
        sum(M * log(Theta$P %*% Theta$E)) -
        sum(logfac[M])
}

#' Run gibbs sampler for poisson-gamma single-study NMF
#'
#' @param M mutational catalog matrix, K x G
#' @param N number of signatures
#' @param P (optional) initial signatures matrix, K x N
#' @param E (optional) initial exposure matrix, N x G
#' @param niters how many iterations the Gibbs sampler is run
#' @param burn_in the first `burn_in` iterations will be ignored when computing
#' MAP estimate
#' @param logevery the log, save, and plot files will be updated every
#' `logevery` iterations
#' @param file file name without extension of log, save, and plot files
#' @param overwrite if `overwrite = TRUE`, the log, safe, and plot files of
#' previous runs with the same `file` will be overwritten
#' @param alpha_p see `Alpha_p`
#' @param Alpha_p shape parameter for the gamma prior on `P`, matrix
#' size K x N. Defaults to all same value `alpha_p`
#' @param beta_p see `Beta_p`
#' @param Beta_p rate parameter for the gamma prior on `P`, matrix
#' size K x N. Defaults to all same value `beta_p`
#' @param alpha_e see `Alpha_e`
#' @param Alpha_e shape parameter for the gamma prior on `E`, matrix
#' size N x G. Defaults to all same value `alpha_e`
#' @param beta_e see `Beta_e`
#' @param Beta_e rate parameter for the gamma prior on `E`, matrix
#' size N x G. Defaults to all same value `beta_e`
#' @param true_P (optional) true signatures matrix P to compare to in a heatmap
#'
#' @return list
#' @export
nmf_poisson_gamma <- function(
        M, N,
        P = NULL,
        E = NULL,
        niters = 10000,
        burn_in = round(2*niters/3),
        logevery = 100,
        file = 'nmf_poisson_gamma',
        overwrite = FALSE,
        alpha_p = 10,
        Alpha_p = matrix(alpha_p, nrow = dim(M)[1], ncol = N),
        beta_p = sqrt(N),
        Beta_p = matrix(beta_p, nrow = dim(M)[1], ncol = N),
        alpha_e = 10,
        Alpha_e = matrix(alpha_e, nrow = N, ncol = dim(M)[2]),
        beta_e = sqrt(N),
        Beta_e = matrix(beta_e, nrow = N, ncol = dim(M)[2]),
        true_P = NULL
) {
    if (burn_in > niters) {
        message(paste0(
            "Burn in ", burn_in, " is greater than niters ",
            niters, ", setting burn_in = 0"
        ))
        burn_in = 0
    }

    savefile = paste0(file, '.res')
    logfile = paste0(file, '.log')
    plotfile = paste0(file, '.pdf')
    tail = 0
    while (!overwrite & (file.exists(savefile) | file.exists(logfile))) {
        tail = tail + 1
        savefile = paste0(file, '_', tail, '.res')
        logfile = paste0(file, '_', tail, '.log')
        plotfile = paste0(file, '_', tail, '.pdf')
    }
    sink(file = logfile)
    print(Sys.time())

    dims = list(
        K = dim(M)[1],
        G = dim(M)[2],
        N = N
    )

    Theta <- list(
        P = P,
        E = E,
        Alpha_p = Alpha_p,
        Beta_p = Beta_p,
        Alpha_e = Alpha_e,
        Beta_e = Beta_e
    )

    logfac = vector(length = max(M))
    logfac[1] = 0
    for (i in 2:length(logfac)) {
        logfac[i] = log(i) + logfac[i-1]
    }

    if (is.null(P)) {
        Theta$P <- matrix(nrow = dims$K, ncol = dims$N)
        for (k in 1:dims$K) {
            for (n in 1:dims$N) {
                Theta$P[k,n] <- stats::rgamma(1, Theta$Alpha_p[k,n], Theta$Beta_p[k,n])
            }
        }
    }
    if (is.null(E)) {
        Theta$E <- matrix(nrow = dims$N, ncol = dims$G)
        for (n in 1:dims$N) {
            for (g in 1:dims$G) {
                Theta$E[n,g] <- stats::rgamma(1, Theta$Alpha_e[n,g], Theta$Beta_e[n,g])
            }
        }
    }
    Theta$Z <- array(dim = c(dims$K, dims$N, dims$G))
    for (k in 1:dims$K) {
        for (g in 1:dims$G) {
            Theta$Z[k,,g] <- sample_Zkg_poisson_gamma_single(k, g, M, Theta, dims)
        }
    }

    RMSE <- c()
    KL <- c()
    loglik <- c()

    P.log <- list()
    E.log <- list()
    Z.log <- list()

    for (iter in 1:niters) {
        for (k in 1:dims$K) {
            for (n in 1:dims$N) {
                Theta$P[k, n] <- sample_Pkn_poisson_gamma_single(k, n, M, Theta)
            }
        }

        for (n in 1:N) {
            for (g in 1:dims$G) {
                Theta$E[n, g] <- sample_Eng_poisson_gamma_single(n, g, M, Theta)
            }
        }

        for (k in 1:dims$K) {
            for (g in 1:dims$G) {
                Theta$Z[k,,g] <- sample_Zkg_poisson_gamma_single(k, g, M, Theta, dims)
            }
        }

        RMSE <- c(RMSE, get_RMSE(M, Theta))
        KL <- c(KL, get_KLDiv(M, Theta))
        loglik <- c(loglik, get_loglik_poisson_gamma_single(M, Theta, dims, logfac))

        P.log[[iter]] <- Theta$P
        E.log[[iter]] <- Theta$E
        Z.log[[iter]] <- Theta$Z

        if (iter %% logevery == 0 | iter == niters) {
            cat(paste(iter, "/", niters, "\n"))

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
            grDevices::dev.off()

            keep = burn_in:length(P.log)
            res <- list(
                M = M,
                true_P = true_P,
                P.log = P.log,
                E.log = E.log,
                Z.log = Z.log,
                P.mean = Reduce(`+`, P.log[keep])/length(keep),
                E.mean = Reduce(`+`, E.log[keep])/length(keep),
                Z.mean = Reduce(`+`, Z.log[keep])/length(keep),
                burn_in = burn_in,
                loglik.chain = loglik,
                RMSE.chain = RMSE,
                KLDiv.chain = KL,
                final_metrics = list(
		    Theta = Theta,
		    dims = dims,
                    loglik = loglik[iter],
                    RMSE = RMSE[iter],
                    KLDiv = KL[iter]
                )
            )
            save(res, file = savefile)
        }
    }
    if (!is.null(true_P)) {
        sim_mat <- pairwise_sim(res$P.mean, true_P, which = 'cols')
        heatmap <- get_heatmap(res$P.mean, true_P)

        res$sim_mat <- sim_mat
        res$heatmap <- heatmap
        save(res, file = savefile)
    }
    sink()
    return(res)
}
