#' sample sigmasq for Normal likelihood
#'
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_sigmasq_normal <- function(M, Theta, dims, gamma = 1){
    Mhat <- get_Mhat(Theta)
    sigmasq <- sapply(1:dims$G, function(g) {
        invgamma::rinvgamma(
            1,
            shape = Theta$Alpha[g] + gamma * dims$K/2,
            rate = Theta$Beta[g] + gamma * (1/2) * sum((M[,g] - Mhat[,g])**2)
        )
    })
    return(sigmasq)
}


#' sample Z[k,,g] for Poisson likelihood
#'
#' @param k mutation type index
#' @param g tumor genome index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_Zkg_poisson <- function(k, g, M, Theta, dims, gamma = 1){
    probs = sapply(1:dims$N, function(n) {
        Theta$P[k,n] * Theta$A[1,n] * Theta$E[n,g]
    })
    if (sum(probs) == 0) {
        # happens if all factors are excluded
        # setting Zs to 0 means P, E sampled from their priors
        return(rep(0, dims$N))
    }
    probs = probs/sum(probs)
    rmultinom(1, size = M[k,g], prob = probs)
}

#' Sample n
#'
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param clip numeric, prior probabilities of inclusion will be clipped by
#' `clip`/N away from 0 and 1
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_n <- function(Theta, dims, clip, gamma = 1) {
    probs = sapply(0:dims$N, function(n) {
        tmp <- list(n = n)
        prob_A <- update_q(tmp, dims, clip)
        prob <- 1/(dims$N + 1) * prod(
            prob_A ** Theta$A * (1-prob_A)**(1-Theta$A)
        ) ** gamma
        return(prob)
    })
    probs = probs/sum(probs)
    n <- sample(0:dims$N, size = 1, prob = probs)
    return(n)
}

#' Sample An
#'
#' @param n integer, signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimension values
#' @param likelihood string, one of c('normal','poisson')
#' @param prior string, one of c('gamma','exponential','truncnormal')
#' @param logfac vector, logfac[i] = log(i!), use NULL if
#' `likelihood == 'normal'`.
#' @param gamma double, tempering parameter
#'
#' @return integer
#' @noRd
sample_An <- function(n, M, Theta, dims, likelihood, prior, logfac, sparse_rank, gamma) {
    Theta_A0 <- Theta
    Theta_A0$A[1,n] <- 0

    Theta_A1 <- Theta
    Theta_A1$A[1,n] <- 1

    if (likelihood == 'poisson') {
        loglik_0 <- get_loglik_poisson(M, Theta_A0, dims, logfac)
        loglik_1 <- get_loglik_poisson(M, Theta_A1, dims, logfac)
    } else {
        # likelihood == 'normal'
        loglik_0 <- get_loglik_normal(M, Theta_A0, dims)
        loglik_1 <- get_loglik_normal(M, Theta_A1, dims)
    }

    if (sparse_rank) {
        n_params_0 <- sum(Theta_A0$A) * (dims$G + dims$K)
        n_params_1 <- sum(Theta_A1$A) * (dims$G + dims$K)

        neg_BIC_0 <- 2 * loglik_0 - n_params_0 * log(dims$G)
        neg_BIC_1 <- 2 * loglik_1 - n_params_1 * log(dims$G)

        log_p0 = log(1 - Theta$q[n]) + gamma * neg_BIC_0
        log_p1 = log(Theta$q[n]) + gamma * neg_BIC_1
    } else {
        log_p0 = log(1 - Theta$q[n]) + gamma * loglik_0
        log_p1 = log(Theta$q[n]) + gamma * loglik_1
    }


    log_p = log_p1 - sumLog(c(log_p0, log_p1))
    p = exp(log_p)
    if (is.na(p)) {
        if (log_p1 > log_p0) {
            p = 1
        } else if (log_p1 < logp0) {
            p = 0
        } else {
            p = 0.5
        }
    }

    sampled <- sample(c(0, 1), size = 1, prob = c(1-p, p))
    return(list(
        sampled = sampled,
        prob_inclusion = p
    ))
}

#' Update prior probability of inclusion, q
#'
#' @param Theta list of parameters
#' @param dims list of dimension values
#' @param clip numeric, prior probabilities of inclusion will be clipped by
#' `clip`/N away from 0 and 1
#'
#' @return scalar
#' @noRd
update_q <- function(Theta, dims, clip) {
    if (sum(Theta$recovery) > 0) {
        Theta$q <- c(
            rep(Theta$weight_r * Theta$n / sum(Theta$recovery), sum(Theta$recovery)),
            rep((1-Theta$weight_r) * Theta$n / sum(Theta$recovery == 0), sum(Theta$recovery == 0))
        )
    } else {
        Theta$q <- rep(Theta$n/dims$N, dims$N)
    }
    Theta$q[Theta$q <= 0] <- clip/dims$N
    Theta$q[Theta$q >= 1] <- 1 - clip/dims$N
    return(Theta$q)
}
