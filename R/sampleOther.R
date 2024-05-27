#' sample sigmasq for Normal likelihood
#'
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_sigmasq_normal <- function(M, Theta, dims, sigmasq_type, gamma = 1){
    Mhat <- get_Mhat(Theta)
    Sigmasq <- matrix(nrow = dims$K, ncol = 1)
    if (sigmasq_type == 'invgamma') {
        for (k in 1:dims$K) {
            Sigmasq[k,1] <- invgamma::rinvgamma(
                n = 1,
                shape = Theta$Alpha[k,1] + gamma * dims$G / 2,
                rate = Theta$Beta[k,1] + gamma * sum(((M - Mhat)[k,])**2 / (2 * colSums(M)))
            )
        }
    } else if (sigmasq_type == 'noninformative') {
        for (k in 1:dims$K) {
            Sigmasq[k,] <- invgamma::rinvgamma(
                n = dims$G,
                shape = dims$G / 2,
                rate = sum(((M - Mhat)[k,])**2 / (2 * colSums(M)))
            )
        }
        # Sigmasq[Sigmasq > max(M)] <- max(M)
    }

    return(Sigmasq)
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

#' Sample An
#'
#' @param n integer, signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimension values
#' @param logfac vector, logfac[i] = log(i!), only needed for `likelihood = 'poisson'`
#' @param likelihood string, one of c('normal','poisson')
#' @param gamma double, tempering parameter
#'
#' @return integer
#' @noRd
sample_An <- function(n, M, Theta, dims, logfac, likelihood = 'normal', prior = "truncnormal", gamma = 1) {
    Theta_A0 <- Theta
    Theta_A0$A[1,n] <- 0

    Theta_A1 <- Theta
    Theta_A1$A[1,n] <- 1

    if (likelihood == 'normal') {
        loglik_0 <- get_loglik_normal(M, Theta_A0, dims)
        loglik_1 <- get_loglik_normal(M, Theta_A1, dims)
    } else if (likelihood == 'poisson') {
        loglik_0 <- get_loglik_poisson(M, Theta_A0, dims, logfac)
        loglik_1 <- get_loglik_poisson(M, Theta_A1, dims, logfac)
    }

    log_p0 = log(1 - Theta$q[1,n]) + gamma * loglik_0
    log_p1 = log(Theta$q[1,n]) + gamma * loglik_1

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
