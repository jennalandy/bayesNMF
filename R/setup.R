#' Set prior parameters for Truncated Normal prior
#'
#' @param Theta list of parameters
#' @param dims named list of dimensions N, K, G
#' @param mean_p mean for the prior location parameter of `P`
#' @param m_p see `M_p`
#' @param M_p location parameter for the prior location parameter of `P`
#' size K x N. Defaults to all same value `m_p`
#' @param s_p see `S_p`
#' @param S_p scale prameter for the prior location parameter of `P`
#' size K x N. Defaults to all same value `s_p`
#' @param a_p see `A_p`
#' @param A_p shape parameter for the gamma prior on `P`, matrix
#' @param b_p see `B_p`
#' @param B_p rate parameter for the gamma prior on `P`, matrix
#' @param mean_e mean for the prior location parameter of `E`
#' @param m_e see `M_e`
#' @param M_e location parameter for the prior location parameter of `E`
#' size N x G. Defaults to all same value `m_e`
#' @param s_e see `S_e`
#' @param S_e scale prameter for the prior location parameter of `E`
#' size N x G. Defaults to all same value `s_e`
#' @param a_e see `A_e`
#' @param A_e shape parameter for the gamma prior on `E`, matrix
#' size N x G. Defaults to all same value `a_e`
#' @param b_e see `B_e`
#' @param B_e rate parameter for the gamma prior on `E`, matrix
#' size N x G. Defaults to all same value `b_e`
#' @param alpha see `Alpha`
#' @param Alpha shape parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `alpha`
#' @param beta see `Beta`
#' @param Beta rate parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `beta`
#' @param a shape1 parameter for beta prior on `q`
#' @param b shape2 parameter for beta prior on `q`
#'
#' @return named list of prior parameters
#' @noRd
set_truncnorm_hyperprior_parameters <- function(
        Theta,
        dims,
        mean_p = 10/sqrt(dims$N),
        m_p = uniroot(function(x) {
            x + sqrt(x) * dnorm(-sqrt(x))/(pnorm(-sqrt(x)) - 1) - mean_p
        }, interval = c(0, 100))$root,
        M_p = matrix(m_p, nrow = dims$K, ncol = dims$N),
        s_p = m_p, #mu_p/10,
        S_p = matrix(s_p, nrow = dims$K, ncol = dims$N),
        a_p = dims$N + 1,
        A_p = matrix(a_p, nrow = dims$K, ncol = dims$N),
        b_p = sqrt(dims$N) + 1,
        B_p = matrix(b_p, nrow = dims$K, ncol = dims$N),
        mean_e = 10/sqrt(dims$N),
        m_e = uniroot(function(x) {
            x + sqrt(x) * dnorm(-sqrt(x))/(pnorm(-sqrt(x)) - 1) - mean_e
        }, interval = c(0, 100))$root,
        M_e = matrix(m_e, nrow = dims$N, ncol = dims$G),
        s_e = m_e, #mu_e/10,
        S_e = matrix(s_e, nrow = dims$N, ncol = dims$G),
        a_e = dims$N + 1,
        A_e = matrix(a_e, nrow = dims$N, ncol = dims$G),
        b_e = sqrt(dims$N) + 1,
        B_e = matrix(b_e, nrow = dims$N, ncol = dims$G),
        alpha = 0.1,
        Alpha = rep(alpha, dims$K),
        beta = 0.1,
        Beta = rep(beta, dims$K),
        a = 0.8,
        b = 0.8
) {
    if ("m_p" %in% names(Theta) & !("M_p" %in% names(Theta))) {
        Theta$M_p = matrix(Theta$m_p, nrow = dims$K, ncol = dims$N)
    } else if ((
        "mean_p" %in% names(Theta) &
        !("m_p" %in% names(Theta)) &
        !("M_p" %in% names(Theta))
    )) {
        Theta$m_p = uniroot(function(x) {
            x + sqrt(x) * dnorm(-sqrt(x))/(pnorm(-sqrt(x)) - 1) - Theta$mean_p
        }, interval = c(0, 100))$root
        Theta$M_p = matrix(Theta$m_p, nrow = dims$N, ncol = dims$G)
    }

    if ("m_e" %in% names(Theta) & !("M_e" %in% names(Theta))) {
        Theta$M_e = matrix(Theta$m_e, nrow = dims$N, ncol = dims$G)
    } else if ((
        "mean_e" %in% names(Theta) &
        !("m_e" %in% names(Theta)) &
        !("M_e" %in% names(Theta))
    )) {
        Theta$m_e = uniroot(function(x) {
            x + sqrt(x) * dnorm(-sqrt(x))/(pnorm(-sqrt(x)) - 1) - Theta$mean_e
        }, interval = c(0, 100))$root
        Theta$M_e = matrix(Theta$m_e, nrow = dims$N, ncol = dims$G)
    }

    for (matrix in c("S_p", "A_p", "B_p")) {
        element = tolower(matrix)
        print(paste(element, matrix))
        Theta <- fill_matrix(
            Theta,
            element = element, matrix = matrix,
            nrow = dims$K, ncol = dims$N
        )
        print(Theta[[matrix]])
    }

    for (matrix in c("S_e", "A_e", "B_e")) {
        element = tolower(matrix)
        Theta <- fill_matrix(
            Theta,
            element = element, matrix = matrix,
            nrow = dims$N, ncol = dims$G
        )
    }

    if ("alpha" %in% names(Theta) & !("Alpha" %in% names(Theta))) {
        Theta$Alpha = rep(Theta$alpha, dims$K)
    }
    if ("beta" %in% names(Theta) & !("Beta" %in% names(Theta))) {
        Theta$Beta = rep(Theta$beta, dims$K)
    }

    fill_list(Theta, list(
        M_p = M_p,
        M_e = M_e,
        S_p = S_p,
        S_e = S_e,
        A_p = A_p,
        A_e = A_e,
        B_p = B_p,
        B_e = B_e,
        Alpha = Alpha,
        Beta = Beta,
        a = a,
        b = b
    ))
}

#' Sample Truncated Normal Prior Parameters from Hyperpriors
#'
#' @param Theta list of parameters
#' @param dims named list of dimensions N, K, G
#'
#' @return updated list of parameters
#'
#' @noRd
sample_truncnormal_prior_parameters <- function(Theta, dims, recovery, recovery_priors) {
    Theta$Mu_p <- matrix(nrow = dims$K, ncol = dims$N)
    Theta$Sigmasq_p <- matrix(nrow = dims$K, ncol = dims$N)
    for (k in 1:dims$K) {
        for (n in 1:dims$N) {
            Theta$Mu_p[k,n] <- rnorm(
                1, Theta$M_p[k,n], sqrt(Theta$S_p[k,n])
            )
            Theta$Sigmasq_p[k,n] <- invgamma::rinvgamma(
                n = 1, shape = Theta$A_p[k, n], rate = Theta$B_p[k, n]
            )
        }
    }

    Theta$Mu_e <- matrix(nrow = dims$N, ncol = dims$G)
    Theta$Sigmasq_e <- matrix(nrow = dims$N, ncol = dims$G)
    for (n in 1:dims$N) {
        for (g in 1:dims$G) {
            Theta$Mu_e[n,g] <- rnorm(
                1, Theta$M_e[n,g], sqrt(Theta$S_e[n,g])
            )
            Theta$Sigmasq_e[n,g] <- invgamma::rinvgamma(
                n = 1, shape = Theta$A_e[n, g], rate = Theta$B_e[n, g]
            )
        }
    }

    if (recovery) {
        Theta$Mu_p[,1:recovery_priors$N_r] <- recovery_priors$Mu_p
        Theta$Sigmasq_p[,1:recovery_priors$N_r] <- recovery_priors$Sigmasq_p
    }

    return(Theta)
}

#' Set prior parameters for Exponential prior
#'
#' @param Theta list of parameters
#' @param dims named list of dimensions N, K, G
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
#'
#' @return named list of prior parameters
#' @noRd
set_exponential_hyperprior_parameters <- function(
        Theta,
        dims,
        a_p = sqrt(dims$N),
        A_p = matrix(a_p, nrow = dims$K, ncol = dims$N),
        b_p = 10,
        B_p = matrix(b_p, nrow = dims$K, ncol = dims$N),
        a_e = sqrt(dims$N),
        A_e = matrix(a_e, nrow = dims$N, ncol = dims$G),
        b_e = 10,
        B_e = matrix(b_e, nrow = dims$N, ncol = dims$G),
        alpha = 0.1,
        Alpha = rep(alpha, dims$K),
        beta = 0.1,
        Beta = rep(beta, dims$K),
        a = 0.8,
        b = 0.8
) {
    for (matrix in c("A_p", "B_p")) {
        element = tolower(matrix)
        Theta <- fill_matrix(
            Theta,
            element = element, matrix = matrix,
            nrow = dims$K, ncol = dims$N
        )
    }
    for (matrix in c("A_e", "B_e")) {
        element = tolower(matrix)
        Theta <- fill_matrix(
            Theta,
            element = element, matrix = matrix,
            nrow = dims$N, ncol = dims$G
        )
    }

    if ("alpha" %in% names(Theta) & !("Alpha" %in% names(Theta))) {
        Theta$Alpha = rep(Theta$alpha, dims$K)
    }
    if ("beta" %in% names(Theta) & !("Beta" %in% names(Theta))) {
        Theta$Beta = rep(Theta$beta, dims$K)
    }

    fill_list(Theta, list(
        A_p = A_p,
        B_p = B_p,
        A_e = A_e,
        B_e = B_e,
        Alpha = Alpha,
        Beta = Beta,
        a = a,
        b = b
    ))
}

#' Sample Exponential Prior Parameters from Hyperprior Distributions
#'
#' @param Theta list of parameters
#' @param dims named list of dimensions N, K, G
#'
#' @return updated list of parameters
#' @noRd
sample_exponential_prior_parameters <- function(Theta, dims, recovery, recovery_priors) {
    Theta$Lambda_p <- matrix(nrow = dims$K, ncol = dims$N)
    for (k in 1:dims$K) {
        for (n in 1:dims$N) {
            Theta$Lambda_p[k,n] <- rgamma(1, Theta$A_p[k,n], Theta$B_p[k,n])
        }
    }

    Theta$Lambda_e <- matrix(nrow = dims$N, ncol = dims$G)
    for (n in 1:dims$N) {
        for (g in 1:dims$G) {
            Theta$Lambda_e[n,g] <- rgamma(1, Theta$A_e[n,g], Theta$B_e[n,g])
        }
    }

    if (recovery) {
        Theta$Lambda_p[,1:recovery_priors$N_r] <- recovery_priors$Lambda_p
    }

    return(Theta)
}

#' Set prior parameters for Gamma prior
#'
#' @param Theta list of parameters
#' @param dims named list of dimensions N, K, G
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
#' #'
#' @return named list of prior parameters
#' @noRd
set_gamma_hyperprior_parameters <- function(
        Theta,
        dims,
        a_p = 10 * sqrt(dims$N),
        A_p = matrix(a_p, nrow = dims$K, ncol = dims$N),
        b_p = 10,
        B_p = matrix(b_p, nrow = dims$K, ncol = dims$N),
        c_p = 100,
        C_p = matrix(c_p, nrow = dims$K, ncol = dims$N),
        d_p = 10,
        D_p = matrix(d_p, nrow = dims$K, ncol = dims$N),
        a_e = 10 * sqrt(dims$N),
        A_e = matrix(a_e, nrow = dims$N, ncol = dims$G),
        b_e = 10,
        B_e = matrix(b_e, nrow = dims$N, ncol = dims$G),
        c_e = 100,
        C_e = matrix(c_e, nrow = dims$N, ncol = dims$G),
        d_e = 10,
        D_e = matrix(d_e, nrow = dims$N, ncol = dims$G),
        a = 0.8,
        b = 0.8
) {

    for (matrix in c("A_p", "B_p", "C_p", "D_p")) {
        element = tolower(matrix)
        Theta <- fill_matrix(
            Theta,
            element = element, matrix = matrix,
            nrow = dims$K, ncol = dims$N
        )
    }

    for (matrix in c("A_e", "B_e", "C_e", "D_e")) {
        element = tolower(matrix)
        Theta <- fill_matrix(
            Theta,
            element = element, matrix = matrix,
            nrow = dims$N, ncol = dims$G
        )
    }

    fill_list(Theta, list(
        A_p = A_p,
        A_e = A_e,
        B_p = B_p,
        B_e = B_e,
        C_p = C_p,
        C_e = C_e,
        D_p = D_p,
        D_e = D_e,
        a = a,
        b = b
    ))
}

sample_gamma_prior_parameters <- function(Theta, dims, recovery, recovery_priors) {
    Theta$Alpha_p <- matrix(nrow = dims$K, ncol = dims$N)
    Theta$Beta_p <- matrix(nrow = dims$K, ncol = dims$N)
    for (k in 1:dims$K) {
        for (n in 1:dims$N) {
            Theta$Beta_p[k,n] <- rgamma(1, Theta$A_p[k,n], Theta$B_p[k,n])
            Theta$Alpha_p[k,n] <- rgamma(1, Theta$C_p[k,n], Theta$D_p[k,n])
        }
    }

    Theta$Alpha_e <- matrix(nrow = dims$N, ncol = dims$G)
    Theta$Beta_e <- matrix(nrow = dims$N, ncol = dims$G)
    for (n in 1:dims$N) {
        for (g in 1:dims$G) {
            Theta$Beta_e[n,g] <- rgamma(1, Theta$A_e[n,g], Theta$B_e[n,g])
            Theta$Alpha_e[n,g] <- rgamma(1, Theta$C_e[n,g], Theta$D_e[n,g])
        }
    }

    if (recovery) {
        Theta$Alpha_p[,1:recovery_priors$N_r] <- recovery_priors$Alpha_p
        Theta$Beta_p[,1:recovery_priors$N_r] <- recovery_priors$Beta_p
    }

    return(Theta)
}

#' Sample from prior distribution of P
#'
#' @param Theta list of parameters
#' @param dims named list of dimensions N, K, G
#' @param prior string, one of c('exponential','truncnormal','gamma')
#'
#' @return matrix, prior sample of P
#' @noRd
sample_prior_P <- function(Theta, dims, prior) {
    P <- matrix(nrow = dims$K, ncol = dims$N)
    for (k in 1:dims$K) {
        for (n in 1:dims$N) {
            if (prior == 'truncnormal') {
                P[k,n] <- truncnorm::rtruncnorm(
                    1, a = 0, b = Inf,
                    mean = Theta$Mu_p[k,n],
                    sd = sqrt(Theta$Sigmasq_p[k,n])
                )
            } else if (prior == 'exponential') {
                P[k,n] <- stats::rexp(1, Theta$Lambda_p[k,n])
            } else if (prior == 'gamma') {
                P[k,n] <- stats::rgamma(1, Theta$Alpha_p[k,n], Theta$Beta_p[k,n])
            }
        }
    }
    return(P)
}

#' Sample from prior distribution of E
#'
#' @param Theta list of parameters
#' @param dims named list of dimensions N, K, G
#' @param prior string, one of c('exponential','truncnormal','gamma')
#'
#' @return matrix, prior sample of E
#' @noRd
sample_prior_E <- function(Theta, dims, prior) {
    E <- matrix(nrow = dims$N, ncol = dims$G)
    for (n in 1:dims$N) {
        for (g in 1:dims$G) {
            if (prior == 'truncnormal') {
                E[n,g] <- truncnorm::rtruncnorm(
                    1, a = 0, b = Inf,
                    mean = Theta$Mu_e[n,g],
                    sd = sqrt(Theta$Sigmasq_e[n,g])
                )
            } else if (prior == 'exponential') {
                E[n,g] <- stats::rexp(1, Theta$Lambda_e[n,g])
            } else if (prior == 'gamma') {
                E[n,g] <- stats::rgamma(1, Theta$Alpha_e[n,g], Theta$Beta_e[n,g])
            }
        }
    }
    return(E)
}

#' Sample from prior distribution of sigmasq
#'
#' @param Theta list of parameters
#' @param dims named list of dimensions N, K, G
#'
#' @return matrix, prior sample of sigmasq
#' @noRd
sample_prior_sigmasq <- function(Theta, dims, sigmasq_type) {
    if (sigmasq_type == 'invgamma') {
        sapply(1:dims$K, function(k) {
            1/rgamma(n = 1, shape = Theta$Alpha[k], rate = Theta$Beta[k])
        })
    } else if (sigmasq_type == 'noninformative') {
        armspp::arms(n_samples = dims$K, log_pdf = function(x) {-1*log(x)}, lower = 0, upper = 1000)
    } else if (sigmasq_type == 'eq_mu') {
        rowMeans(get_Mhat(Theta))
    }

}

#' initialize Theta
#'
#' @param likelihood string, one of c('poisson','normal')
#' @param prior string, one of c('exponential','truncnormal','gamma')
#' @param learn_A boolean, whether A will be learned or will be fixed
#' @param dims named list of dimensions
#' @param inits list of initial values, optional
#' @param prior_parameters list of named prior parameters, optional
#'
#' @return named list of initialized unknowns
#' @noRd
initialize_Theta <- function(
        M, likelihood, prior, learn_A,
        dims, sigmasq_type,
        inits,
        fixed,
        prior_parameters,
        recovery,
        recovery_priors
) {
    is_fixed = list(
        A = !learn_A,
        prior_P = rep(FALSE, dims$N)
    )

    if (recovery) {
        is_fixed$prior_P[1:recovery_priors$N_r] <- TRUE
    }

    # hyperprior and prior parameters
    Theta = prior_parameters
    if (prior == 'truncnormal') {
        Theta <- set_truncnorm_hyperprior_parameters(Theta, dims)
        Theta <- sample_truncnormal_prior_parameters(Theta, dims, recovery, recovery_priors)
    } else if (prior == 'exponential') {
        Theta <- set_exponential_hyperprior_parameters(Theta, dims)
        Theta <- sample_exponential_prior_parameters(Theta, dims, recovery, recovery_priors)
    } else if (prior == 'gamma') {
        Theta <- set_gamma_hyperprior_parameters(Theta, dims)
        Theta <- sample_gamma_prior_parameters(Theta, dims, recovery, recovery_priors)
    }

    # signatures P
    if (!is.null(fixed$P)) {
        Theta$P <- fixed$P
        is_fixed$P <- TRUE
    } else if (!is.null(inits$P)) {
        Theta$P <- inits$P
        is_fixed$P <- FALSE
    } else {
        Theta$P <- sample_prior_P(Theta, dims, prior)
        is_fixed$P <- FALSE
    }

    # exposures E
    if (!is.null(fixed$E)) {
        Theta$E <- fixed$E
        is_fixed$E <- TRUE
    } else if (!is.null(inits$E)) {
        Theta$E <- inits$E
        is_fixed$E <- FALSE
    } else {
        Theta$E <- sample_prior_E(Theta, dims, prior)
        is_fixed$E <- FALSE
    }

    # signature assignment A
    if (!is.null(fixed$q)) {
        Theta$q <- fixed$q
        is_fixed$q <- TRUE
    } else if (!is.null(inits$q)) {
        Theta$q <- inits$q
        is_fixed$q <- FALSE
    } else {
        Theta$q <- matrix(
            rbeta(dims$S * dims$N, Theta$a, Theta$b),
            nrow = dims$S, ncol = dims$N
        )
        is_fixed$q <- FALSE
    }

    if (!is.null(fixed$A)) {
        is_fixed$A <- TRUE
    } else if (!is.null(inits$A)) {
        Theta$A <- inits$A
    } else if (!learn_A) {
        Theta$A <- matrix(
            as.numeric(rep(1, dims$S * dims$N)),
            nrow = dims$S, ncol = dims$N
        )
    } else {
        Theta$A <- matrix(
            as.numeric(runif(dims$S * dims$N) < c(Theta$q)),
            nrow = dims$S, ncol = dims$N
        )
    }

    if (likelihood == 'normal') {
        if (!is.null(fixed$sigmasq)) {
            Theta$sigmasq <- fixed$sigmasq
            is_fixed$sigmasq <- TRUE
        } else if (!is.null(inits$sigmasq)) {
            Theta$sigmasq <- inits$sigmasq
            is_fixed$sigmasq <- FALSE
        } else {
            Theta$sigmasq <- sample_prior_sigmasq(Theta, dims, sigmasq_type)
            is_fixed$sigmasq <- FALSE
        }
    } else if (likelihood == 'poisson') {
        Theta$Z <- array(dim = c(dims$K, dims$N, dims$G))
        for (k in 1:dims$K) {
            for (g in 1:dims$G) {
                Theta$Z[k,,g] <- sample_Zkg_poisson(k, g, M, Theta, dims)
            }
        }
    }

    Theta$is_fixed <- is_fixed

    return(Theta)
}
