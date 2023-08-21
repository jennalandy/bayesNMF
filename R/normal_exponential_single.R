#' get mu and sigmasq for P[,n]
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#'
#' @return list of two items, mu and sigmasq
get_mu_sigmasq_Pn <- function(n, M, Theta) {
    Mhat_no_n <- Theta$P[, -n] %*% Theta$E[-n, ]
    sum_E_sq <- sum(Theta$E[n, ] ** 2)

    # compute mean
    mu_num_term_1 <- sweep(
        (M - Mhat_no_n), # dim KxG
        2, # multiply each row by E[n,]
        Theta$E[n, ], # length G
        "*"
    ) %>% # dim KxG
        rowSums() # length K

    mu_num_term_2 <- Theta$Lambda_p[, n] * Theta$sigmasq # length K
    mu_P <- (mu_num_term_1 - mu_num_term_2) / sum_E_sq # length K
    sigmasq_P <- Theta$sigmasq / sum_E_sq # length K

    return(list(
        mu = mu_P,
        sigmasq = sigmasq_P
    ))
}

#' sample P[,n]
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#'
#' @return vector length K
sample_Pn <- function(n, M, Theta) {
    mu_sigmasq_P <- get_mu_sigmasq_Pn(n, M, Theta)
    mu_P = mu_sigmasq_P$mu
    sigmasq_P = mu_sigmasq_P$sigmasq

    # sample from truncated normal
    truncnorm::rtruncnorm(1, mean = mu_P, sd = sqrt(sigmasq_P), a = 0, b = Inf)
}

#' Get mu and sigmasq for E[n,]
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#'
#' @return list of two items, mu and sigmasq
get_mu_sigmasq_En <- function(n, M, Theta) {
    Mhat_no_n <- Theta$P[, -n] %*% Theta$E[-n, ]

    # compute mean
    mu_num_term_1 <- sweep(
        (M - Mhat_no_n), # dim KxG
        1, # multiply each column by P[,n]
        Theta$P[, n] / Theta$sigmasq, # length K
        "*"
    ) %>% # dim KxG
        colSums() # length G
    mu_num_term_2 <- Theta$Lambda_e[n, ] # length G
    denom <- sum(Theta$P[, n] ** 2 / Theta$sigmasq)

    mu_E <- (mu_num_term_1 - mu_num_term_2) / denom # length G
    sigmasq_E <- 1 / denom

    return(list(
        mu = mu_E,
        sigmasq = sigmasq_E
    ))
}

#' Sample E[n,]
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#'
#' @return vector of length G
sample_En <- function(n, M, Theta) {
    # compute mean
    mu_sigmasq_E <- get_mu_sigmasq_En(n, M, Theta)
    mu_E = mu_sigmasq_E$mu
    sigmasq_E = mu_sigmasq_E$sigmasq

    # sample from truncated normal
    truncnorm::rtruncnorm(1, mean = mu_E, sd = sqrt(sigmasq_E), a = 0, b = Inf)
}

#' sample sigmasq
#'
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#'
#' @return vector length K
sample_sigmasq <- function(M, Theta, dims){
    Mhat <- Theta$P %*% Theta$E
    sigmasq <- sapply(1:dims$K, function(k) {
        m_k <- max(M[k,])
        sigmasq_k <- invgamma::rinvgamma(
            n = 1,
            shape = Theta$Alpha[k] + dims$G / 2,
            scale = Theta$Beta[k] + sum(((M - Mhat)[k,])**2) / 2
        )
        sigmasq_k <- min(sigmasq_k, m_k)
    })

    return(sigmasq)
}

#' get RMSE
#'
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#'
#' @return scalar
get_RMSE <- function(M, Theta) {
    sqrt(mean((M - Theta$P %*% Theta$E)**2))
}

#' get log likelihood
#'
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#'
#' @return scalar
get_loglik <- function(M, Theta, dims) {
    - dims$G * sum(log(2 * pi * Theta$sigmasq)) / 2 -
    sum(sweep(
        (M - Theta$P %*% Theta$E)**2,
        1,
        1/(2 * Theta$sigmasq),
        '*'
    ))
}

#' Run gibbs sampler for normal-exponential single-study NMF
#'
#' @param M mutational catalog matrix, K x G
#' @param N number of signatures
#' @param P (optional) initial signatures matrix, K x N
#' @param E (optional) initial exposure matrix, N x G
#' @param sigmasq (optinal) initial variance vector, length K
#' @param niters how many iterations the Gibbs sampler is run
#' @param burn_in the first `burn_in` iterations will be ignored when computing
#' MAP estimate
#' @param logevery the log, save, and plot files will be updated every
#' `logevery` iterations
#' @param file file name without extension of log, save, and plot files
#' @param overwrite if `overwrite = TRUE`, the log, safe, and plot files of
#' previous runs with the same `file` will be overwritten
#' @param lambda_p see `Lambda_p`
#' @param Lambda_p rate parameter for the exponential prior on `P`, matrix
#' size K x N. Defaults to all same value `lambda_p`
#' @param lambda_e see `Lambda_e`
#' @param Lambda_e rate parameter for the exponential prior on `E`, matrix
#' size N x G. Defaults to all same value `lambda_e`
#' @param alpha see `Alpha`
#' @param Alpha shape parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `alpha`
#' @param beta see `Beta`
#' @param Beta rate parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `beta`
#' @param true_Theta (optional) true parameters to compare to in a heatmap
#'
#' @return list
#' @export
nmf_normal_exponential <- function(
    M, N,
    P = NULL,
    E = NULL,
    sigmasq = NULL,
    niters = 1000,
    burn_in = 500,
    logevery = 100,
    file = 'normal_exponential_gibbs_sampler',
    overwrite = FALSE,
    lambda_p = 5,
    Lambda_p = matrix(lambda_p, nrow = dim(M)[1], ncol = N),
    lambda_e = 3,
    Lambda_e = matrix(lambda_e, nrow = N, ncol = dim(M)[2]),
    alpha = 0.1,
    Alpha = rep(alpha, dim(M)[1]),
    beta = 2,
    Beta = rep(beta, dim(M)[1]),
    true_Theta = NULL
) {
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
        sigmasq = sigmasq,
        Lambda_e = Lambda_e,
        Lambda_p = Lambda_p,
        Alpha = Alpha,
        Beta = Beta
    )

    if (is.null(P)) {
        Theta$P <- matrix(nrow = dims$K, ncol = dims$N)
        for (k in 1:dims$K) {
            for (n in 1:dims$N) {
                Theta$P[k,n] <- stats::rexp(1, Theta$Lambda_p[k,n])
            }
        }
    }
    if (is.null(E)) {
        Theta$E <- matrix(nrow = dims$N, ncol = dims$G)
        for (n in 1:dims$N) {
            for (g in 1:dims$G) {
                Theta$E[n,g] <- stats::rexp(1, Theta$Lambda_e[n,g])
            }
        }
    }
    if (is.null(sigmasq)) {
        Theta$sigmasq <- sapply(1:dims$K, function(k) {
            invgamma::rinvgamma(n = 1, shape = Theta$Alpha[k], scale = Theta$Beta[k])
        })
    }

    RMSE <- c()
    RMSE <- c(RMSE, get_RMSE(M, Theta))
    if (!is.null(true_Theta)) {
        true_RMSE <- get_RMSE(M, true_Theta)
    }

    loglik <- c()
    loglik <- c(loglik, get_loglik(M, Theta, dims))
    if (!is.null(true_Theta)) {
        true_loglik <- get_loglik(M, true_Theta, dims)
    }

    P.log <- list()
    E.log <- list()
    sigmasq.log <- list()

    for (iter in 1:niters) {
        for (n in 1:dims$N) {
            Theta$P[, n] <- sample_Pn(n, M, Theta)
        }
        for (n in 1:N) {
            Theta$E[n, ] <- sample_En(n, M, Theta)
        }
        Theta$sigmasq <- sample_sigmasq(M, Theta, dims)

        RMSE <- c(RMSE, get_RMSE(M, Theta))
        loglik <- c(loglik, get_loglik(M, Theta, dims))

        P.log[[iter]] <- Theta$P
        E.log[[iter]] <- Theta$E
        sigmasq.log[[iter]] <- Theta$sigmasq

        if (iter %% logevery == 0 | iter == niters) {
            cat(paste(iter, "/", niters, "\n"))

            grDevices::pdf(plotfile)
            graphics::par(mfrow = c(1,2))
            plot(RMSE)
            if (!is.null(true_Theta)) {
                graphics::abline(h = true_RMSE)
            }
            plot(loglik)
            if (!is.null(true_Theta)) {
                graphics::abline(h = true_loglik)
            }
            grDevices::dev.off()

            keep = burn_in:length(P.log)
            res <- list(
                P.log = P.log,
                E.log = E.log,
                sigmasq.log = sigmasq.log,
                P.mean = Reduce(`+`, P.log[keep])/length(keep),
                E.mean = Reduce(`+`, E.log[keep])/length(keep),
                sigmasq.mean = Reduce(`+`, sigmasq.log[keep])/length(keep),
                burn_in = burn_in
            )
            save(res, file = savefile)
        }
    }
    if (!is.null(true_Theta)) {
        sim_mat <- get_sim_mat(res$P.mean, true_Theta$P)
        heatmap <- get_heatmap(res$P.mean, true_Theta$P)

        res$sim_mat <- sim_mat
        res$heatmap <- heatmap
    }
    return(res)
}
