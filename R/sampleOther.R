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

sample_n <- function(Theta, dims, clip, gamma = 1) {
    probs_A <- (0:dims$N)/dims$N
    names(probs_A) <- 0:dims$N
    probs_A[probs_A == 0] <- clip/dims$N
    probs_A[probs_A == 1] <- 1 - clip/dims$N
    probs = sapply(0:dims$N, function(n) {
        prob_A <- probs_A[as.character(n)]
        (1/(dims$N + 1)) * ((prob_A ** sum(Theta$A)) * ((1 - prob_A) ** (dims$N - sum(Theta$A)))) ** gamma
    })
    probs = probs/sum(probs)
    sample(0:dims$N, size = 1, prob = probs)
}

#' Sample An
#'
#' @param n integer, signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimension values
#' @param logfac vector, logfac[i] = log(i!), only needed for `likelihood = 'poisson'`
#' @param clip numeric, prior probabilities of inclusion will be clipped by
#' `clip`/N away from 0 and 1
#' @param likelihood string, one of c('normal','poisson')
#' @param gamma double, tempering parameter
#'
#' @return integer
#' @noRd
sample_An <- function(n, M, Theta, dims, logfac, clip, likelihood = 'normal', prior = "truncnormal", gamma = 1) {
    Theta_A0 <- Theta
    Theta_A0$A[1,n] <- 0

    Theta_A1 <- Theta
    Theta_A1$A[1,n] <- 1

    loglik_0 <- get_loglik_poisson(M, Theta_A0, dims, logfac)
    loglik_1 <- get_loglik_poisson(M, Theta_A1, dims, logfac)

    n_params_0 <- sum(Theta_A0$A) * (dims$G + dims$K)
    n_params_1 <- sum(Theta_A1$A) * (dims$G + dims$K)

    neg_BIC_0 <- 2 * loglik_0 - n_params_0 * log(dims$G)
    neg_BIC_1 <- 2 * loglik_1 - n_params_1 * log(dims$G)

    prior_prob <- Theta$n/dims$N
    if (prior_prob == 0) {prior_prob = prior_prob + clip/dims$N}
    if (prior_prob == 1) {prior_prob = prior_prob - clip/dims$N}

    # log_p0 = log(1 - prior_prob) + gamma * loglik_0
    # log_p1 = log(prior_prob) + gamma * loglik_1
    log_p0 = log(1 - prior_prob) + gamma * neg_BIC_0
    log_p1 = log(prior_prob) + gamma * neg_BIC_1

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

#' Sample qn
#'
#' @param n integer, signature index
#' @param Theta list of parameters
#' @param gamma double, tempering parameter
#'
#' @return double
#' @noRd
sample_qn <- function(n, Theta, gamma = 1) {
    rbeta(1, Theta$a + gamma*Theta$A[1, n], Theta$b + gamma*(1 - Theta$A[1, n]))
}
