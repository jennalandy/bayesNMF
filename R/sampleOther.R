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
    Mhat <- Theta$P %*% diag(Theta$A[1, ]) %*% Theta$E
    sigmasq <- sapply(1:dims$K, function(k) {
        sigmasq_k <- invgamma::rinvgamma(
            n = 1,
            shape = Theta$Alpha[k] + gamma * dims$G / 2,
            scale = Theta$Beta[k] + gamma * sum(((M - Mhat)[k,])**2) / 2
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
sample_Zkg_poisson <- function(k, g, M, Theta, dims){
    probs = sapply(1:dims$N, function(n) {
        Theta$P[k,n] * Theta$E[n,g]
    })
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
sample_An <- function(n, M, Theta, dims, logfac, likelihood = 'normal', gamma = 1) {
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
    return(sample(c(0, 1), size = 1, prob = c(1-p, p)))
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
