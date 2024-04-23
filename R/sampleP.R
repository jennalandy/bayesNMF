#' get mu and sigmasq for P[,n] in Normal-Exponential model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_Pn_normal_exponential <- function(n, M, Theta, dims, gamma = 1) {
    Mhat_no_n <- get_Mhat_no_n(Theta, dims, n)
    sum_E_sq <- gamma * sum(Theta$A[1,n] * Theta$E[n, ] ** 2)

    # compute mean
    mu_num_term_1 <- gamma * sweep(
        (M - Mhat_no_n), # dim KxG
        2, # multiply each row by E[n,]
        Theta$A[1,n] * Theta$E[n, ], # length G
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

#' get mu and sigmasq for P[,n] in Normal-TruncNormal model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_Pn_normal_truncnormal <- function(n, M, Theta, dims, gamma = 1) {
    Mhat_no_n <- get_Mhat_no_n(Theta, dims, n)

    # compute mean
    mu_num_term_1 <- gamma * (1/Theta$sigmasq) * (sweep(
        (M - Mhat_no_n), # dim KxG
        2, # multiply each row by E[n,]
        Theta$A[1,n] * Theta$E[n, ], # length G
        "*"
    ) %>% # dim KxG
        rowSums()) # length K
    mu_num_term_2 <- Theta$Mu_p[, n] / Theta$Sigmasq_p[,n] # length K
    denom <- (1/Theta$Sigmasq_p[,n]) + gamma * sum(Theta$A[1,n] * Theta$E[n, ] ** 2) / Theta$sigmasq


    mu_P <- (mu_num_term_1 + mu_num_term_2) / denom # length K
    sigmasq_P <- 1 / denom # length K

    return(list(
        mu = mu_P,
        sigmasq = sigmasq_P
    ))
}

#' sample P[,n] for Normal likelihood
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param prior string, one of c('truncnormal','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_Pn_normal <- function(n, M, Theta, dims, prior = 'truncnormal', gamma = 1) {
    if (prior == 'truncnormal') {
        mu_sigmasq_P <- get_mu_sigmasq_Pn_normal_truncnormal(n, M, Theta, dims, gamma = gamma)
    } else if (prior == 'exponential') {
        mu_sigmasq_P <- get_mu_sigmasq_Pn_normal_exponential(n, M, Theta, dims, gamma = gamma)
    }


    mu_P = mu_sigmasq_P$mu
    sigmasq_P = mu_sigmasq_P$sigmasq

    # sample from truncated normal
    truncnorm::rtruncnorm(1, mean = mu_P, sd = sqrt(sigmasq_P), a = 0, b = Inf)
}

#' sample P[,n] for Poisson likelihood
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param prior string, one of c('gamma','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_Pn_poisson <- function(n, M, Theta, dims, prior = 'gamma', gamma = 1) {
    if (prior == 'gamma') {
        sampled <- sapply(1:dims$K, function(k) {
            rgamma(1, Theta$Alpha_p[k,n] + gamma * sum(Theta$Z[k,n,]), Theta$Beta_p[k,n] + gamma * Theta$A[1,n] * sum(Theta$E[n,]))
        })
    } else if (prior == 'exponential') {
        sampled <- sapply(1:dims$K, function(k) {
            rgamma(1, 1 + gamma * sum(Theta$Z[k,n,]), Theta$Lambda_p[k,n] + gamma * Theta$A[1,n] * sum(Theta$E[n,]))
        })
    }
    return(sampled)
}

#' Sample Pkn for Poisson-Exponential Model with adaptive rejection metropolis sampling
#'
#' @param k mutation type index
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Pkn_poisson_exp <- function(k, n, M, Theta, gamma) {
    log_pdf <- function(Pkn) {
        gamma * sum(Theta$Z[k,n,]) * log(Pkn) -
            Pkn * (Theta$Lambda_p[k,n] + gamma * Theta$A[1,n] * sum(Theta$E[n,]))
    }
    sample <- armspp::arms(
        n_samples = 1,
        log_pdf = log_pdf,
        lower = 0,
        upper = 100
    )
    return(sample)
}

#' Sample Pn for Poisson-Exponential Model with adaptive rejection metropolis sampling
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Pn_poisson_exp <- function(n, M, Theta, dims, gamma) {
    sapply(1:dims$K, function(k) {
        sample_Pkn_poisson_exp(k, n, M, Theta, gamma)
    })
}

#' Sample Pkn for Normal-Exponential Model with adaptive rejection metropolis sampling
#'
#' @param k mutation type index
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Pkn_norm_exp <- function(k, n, M, Theta, gamma) {
    log_pdf <- function(Pkn) {
        Theta$P[k,n] <- Pkn
        Mhat <- get_Mhat(Theta)

        -Theta$Lambda_p[k,n] * Pkn -
            gamma/(2*Theta$sigmasq[k]) * sum((M[k,] - Mhat[k,]) ** 2)
    }
    sample <- armspp::arms(
        n_samples = 1,
        log_pdf = log_pdf,
        lower = 0,
        upper = 100
    )
    return(sample)
}

#' Sample Pn for Normal-Exponential Model with adaptive rejection metropolis sampling
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Pn_norm_exp <- function(n, M, Theta, dims, gamma) {
    sapply(1:dims$K, function(k) {
        sample_Pkn_norm_exp(k, n, M, Theta, gamma)
    })
}

#' sample P[,n] wrapper function
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param likelihood string, one of c('normal','exponential')
#' @param prior string, one of c('gamma','exponential','truncnormal')
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_Pn <- function(n, M, Theta, dims, likelihood = 'normal', prior = 'truncnormal', gamma = 1) {
    if (likelihood == 'normal') {
        if (prior == 'truncnormal' | gamma > 0.5) {
            sample_Pn_normal(n, M, Theta, dims, prior, gamma)
        } else {
            sample_Pn_norm_exp(n, M, Theta, dims, gamma)
        }
    } else if (likelihood == 'poisson') {
        if (prior == 'gamma' | gamma > 0.5) {
            sample_Pn_poisson(n, M, Theta, dims, prior, gamma)
        } else {
            sample_Pn_poisson_exp(n, M, Theta, dims, gamma)
        }

    }
}
