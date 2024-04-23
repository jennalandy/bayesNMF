#' Get mu and sigmasq for E[n,] in Normal-Exponential model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_En_normal_exponential <- function(n, M, Theta, dims, gamma = 1) {
    Mhat_no_n <- get_Mhat_no_n(Theta, dims, n)

    # compute mean
    mu_num_term_1 <- gamma * Theta$A[1,n] * sweep(
        (M - Mhat_no_n), # dim KxG
        1, # multiply each column by P[,n]
        Theta$P[, n] / Theta$sigmasq, # length K
        "*"
    ) %>% # dim KxG
        colSums() # length G
    mu_num_term_2 <- Theta$Lambda_e[n, ] # length G
    denom <- sum(gamma * Theta$A[1,n] * Theta$P[, n] ** 2 / Theta$sigmasq)

    mu_E <- (mu_num_term_1 - mu_num_term_2) / denom # length G
    sigmasq_E <- 1 / denom

    return(list(
        mu = mu_E,
        sigmasq = sigmasq_E
    ))
}

#' Get mu and sigmasq for E[n,] in Normal-TruncNormal model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_En_normal_truncnormal <- function(n, M, Theta, dims, gamma = 1) {
    Mhat_no_n <- get_Mhat_no_n(Theta, dims, n)

    # compute mean
    mu_num_term_1 <- gamma * Theta$A[1,n] * sweep(
        (M - Mhat_no_n), # dim KxG
        1, # multiply each column by P[,n]
        Theta$P[, n] / Theta$sigmasq, # length K
        "*"
    ) %>% # dim KxG
        colSums() # length G
    mu_num_term_2 <- Theta$Mu_e[n, ] / Theta$Sigmasq_e[n,] # length G
    denom <- (1/Theta$Sigmasq_e[n,]) + gamma * sum(Theta$A[1,n] * Theta$P[, n] ** 2 / Theta$sigmasq)

    mu_E <- (mu_num_term_1 + mu_num_term_2) / denom # length G
    sigmasq_E <- 1 / denom

    return(list(
        mu = mu_E,
        sigmasq = sigmasq_E
    ))
}

#' Sample E[n,] for Normal likelihood
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param prior string, one of c('truncnormal','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_En_normal <- function(n, M, Theta, dims, prior = 'truncnormal', gamma = 1) {
    if (prior == 'truncnormal') {
        mu_sigmasq_E <- get_mu_sigmasq_En_normal_truncnormal(n, M, Theta, dims, gamma = gamma)
    } else if (prior == 'exponential') {
        mu_sigmasq_E <- get_mu_sigmasq_En_normal_exponential(n, M, Theta, dims, gamma = gamma)
    }
    mu_E = mu_sigmasq_E$mu
    sigmasq_E = mu_sigmasq_E$sigmasq

    # sample from truncated normal
    truncnorm::rtruncnorm(1, mean = mu_E, sd = sqrt(sigmasq_E), a = 0, b = Inf)
}

#' Sample Eng for Normal-Exponential Model with adaptive rejection metropolis sampling
#'
#' @param n signature index
#' @param g tumor genome index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Eng_norm_exp <- function(n, g, M, Theta, gamma) {
    log_pdf <- function(Eng) {
        Theta$E[n,g] <- Eng
        Mhat <- get_Mhat(Theta)
        -Theta$Lambda_e[n,g] * Eng -
            gamma * sum((M[,g] - Mhat[,g]) ** 2 / (2 * Theta$sigmasq))
    }
    armspp::arms(
        n_samples = 1,
        log_pdf = log_pdf,
        lower = 0,
        upper = 100
    )
}

#' Sample En for Normal-Exponential Model with adaptive rejection metropolis sampling
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_En_norm_exp <- function(n, M, Theta, dims, gamma) {
    sapply(1:dims$G, function(g) {
        sample_Eng_norm_exp(n, g, M, Theta, gamma)
    })
}

#' Sample E[n,g] for Poisson likelihood
#'
#' @param n signature index
#' @param g tumor genome index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param prior string, one of c('gamma','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_En_poisson <- function(n, M, Theta, dims, prior = 'gamma', gamma = 1) {
    if (prior == 'gamma') {
        sampled <- sapply(1:dims$G, function(g) {
            rgamma(1, Theta$Alpha_e[n,g] + gamma*sum(Theta$Z[,n,g]), Theta$Beta_e[n,g] + gamma * Theta$A[1,n] * sum(Theta$P[,n]))
        })
    } else if (prior == 'exponential') {
        sampled <- sapply(1:dims$G, function(g) {
            rgamma(1, 1 + gamma * sum(Theta$Z[,n,g]), Theta$Lambda_e[n,g] + gamma * Theta$A[1,n] * sum(Theta$P[,n]))
        })
    }
    return(sampled)
}

#' Sample Eng for Poisson-Exponential Model with adaptive rejection metropolis sampling
#'
#' @param n signature index
#' @param g tumor genome index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Eng_poisson_exp <- function(n, g, M, Theta, gamma) {
    log_pdf <- function(Eng) {
        gamma * sum(Theta$Z[,n,g]) * log(Eng) -
            Eng * (Theta$Lambda_e[n,g] + gamma * Theta$A[1,n] * sum(Theta$P[,n]))
    }
    sample <- armspp::arms(
        n_samples = 1,
        log_pdf = log_pdf,
        lower = 0,
        upper = 100
    )
    return(sample)
}

#' Sample En for Poisson-Exponential Model with adaptive rejection metropolis sampling
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_En_poisson_exp <- function(n, M, Theta, dims, gamma) {
    sapply(1:dims$G, function(g) {
        sample_Eng_poisson_exp(n, g, M, Theta, gamma)
    })
}


#' Sample E[n,] wrapper function
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param likelihood string, one of c('normal','poisson')
#' @param prior string, one of c('truncnormal','exponential','gamma')
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_En <- function(n, M, Theta, dims, likelihood = 'normal', prior = 'truncnormal', gamma = 1) {
    if (likelihood == 'normal') {
        if ((prior == 'truncnormal' | gamma > 0.5) & Theta$A[1,n] == 1) {
            sample_En_normal(n, M, Theta, dims, prior = prior, gamma = gamma)
        } else {
            # aarms sampling when gamma is small and prior is not conjugate
            sample_En_norm_exp(n, M, Theta, dims, gamma = gamma)
        }
    } else if (likelihood == 'poisson') {
        if (prior == 'gamma' | gamma > 0.5) {
            sample_En_poisson(n, M, Theta, dims, prior = prior, gamma = gamma)
        } else {
            # aarms sampling when gamma is small and prior is not conjugate
            sample_En_poisson_exp(n, M, Theta, dims, gamma = gamma)
        }
    }
}
