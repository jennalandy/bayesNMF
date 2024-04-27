#' Set prior parameters for Truncated Normal prior
#'
#' @param Theta list of parameters
#' @param dims named list of dimensions N, K, G
#' @param mu_p see `Mu_p`
#' @param Mu_p mean for the truncated normal prior on `P`, matrix
#' size K x N. Defaults to all same value `mu_p`
#' @param sigmasq_p see `Sigmasq_p`
#' @param Sigmasq_p variance for the truncated normal prior on `P`, matrix
#' size K x N. Defaults to all same value `sigmasq_p`
#' @param mu_e see `Mu_e`
#' @param Mu_e mean for the truncated normal prior on `E`, matrix
#' size N x G. Defaults to all same value `mu_e`
#' @param sigmasq_e see `Sigmasq_e`
#' @param Sigmasq_e variance for the truncated normal prior on `P`, matrix
#' size N x G. Defaults to all same value `sigmasq_e`
#' @param alpha see `Alpha`
#' @param Alpha shape parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `alpha`
#' @param beta see `Beta`
#' @param Beta rate parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `beta`
#'
#' @return named list of prior parameters
#' @export
set_truncnorm_prior_parameters <- function(
        Theta,
        dims,
        mean_p = sqrt(100/dims$N),
        mu_p = uniroot(function(x) {
            x + sqrt(x) * dnorm(-sqrt(x))/(pnorm(-sqrt(x)) - 1) - mean_p
        }, interval = c(0, 100))$root,
        Mu_p = matrix(mu_p, nrow = dims$K, ncol = dims$N),
        sigmasq_p = mu_p, #mu_p/10,
        Sigmasq_p = matrix(sigmasq_p, nrow = dims$K, ncol = dims$N),
        mean_e = sqrt(100/dims$N),
        mu_e = uniroot(function(x) {
            x + sqrt(x) * dnorm(-sqrt(x))/(pnorm(-sqrt(x)) - 1) - mean_e
        }, interval = c(0, 100))$root,
        Mu_e = matrix(mu_e, nrow = dims$N, ncol = dims$G),
        sigmasq_e = mu_e, #mu_e/10,
        Sigmasq_e = matrix(sigmasq_e, nrow = dims$N, ncol = dims$G),
        alpha = 0.1,
        Alpha = rep(alpha, dims$K),
        beta = 0.1,
        Beta = rep(beta, dims$K),
        a = 0.8,
        b = 0.8
) {
    if ("mu_p" %in% names(Theta) & !("Mu_p" %in% names(Theta))) {
        Theta$Mu_p = matrix(Theta$mu_p, nrow = dims$K, ncol = dims$N)
    } else if ((
        "mean_p" %in% names(Theta) &
        !("mu_p" %in% names(Theta)) &
        !("Mu_p" %in% names(Theta))
    )) {
        Theta$mu_p = uniroot(function(x) {
            x + sqrt(x) * dnorm(-sqrt(x))/(pnorm(-sqrt(x)) - 1) - Theta$mean_p
        }, interval = c(0, 100))$root
        Theta$Mu_p = matrix(Theta$mu_p, nrow = dims$N, ncol = dims$G)
    }

    if ("mu_e" %in% names(Theta) & !("Mu_e" %in% names(Theta))) {
        Theta$Mu_e = matrix(Theta$mu_e, nrow = dims$N, ncol = dims$G)
    } else if ((
        "mean_e" %in% names(Theta) &
        !("mu_e" %in% names(Theta)) &
        !("Mu_e" %in% names(Theta))
    )) {
        Theta$mu_e = uniroot(function(x) {
            x + sqrt(x) * dnorm(-sqrt(x))/(pnorm(-sqrt(x)) - 1) - Theta$mean_e
        }, interval = c(0, 100))$root
        Theta$Mu_e = matrix(Theta$mu_e, nrow = dims$N, ncol = dims$G)
    }

    if ("sigmasq_p" %in% names(Theta) & !("Sigmasq_p" %in% names(Theta))) {
        Theta$Sigmasq_p = matrix(Theta$sigmasq_p, nrow = dims$K, ncol = dims$N)
    }
    if ("sigmasq_e" %in% names(Theta) & !("Sigmasq_e" %in% names(Theta))) {
        Theta$Sigmasq_e = matrix(Theta$sigmasq_e, nrow = dims$N, ncol = dims$G)
    }
    if ("alpha" %in% names(Theta) & !("Alpha" %in% names(Theta))) {
        Theta$Alpha = rep(Theta$alpha, dims$K)
    }
    if ("beta" %in% names(Theta) & !("Beta" %in% names(Theta))) {
        Theta$Beta = rep(Theta$beta, dims$K)
    }
    fill_list(Theta, list(
        Mu_p = Mu_p,
        Sigmasq_p = Sigmasq_p,
        Mu_e = Mu_e,
        Sigmasq_e = Sigmasq_e,
        Alpha = Alpha,
        Beta = Beta,
        a = a,
        b = b
    ))
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
#' @export
set_exponential_prior_parameters <- function(
        Theta,
        dims,
        lambda_p = sqrt(dims$N/100),
        Lambda_p = matrix(lambda_p, nrow = dims$K, ncol = dims$N),
        lambda_e = sqrt(dims$N/100),
        Lambda_e = matrix(lambda_e, nrow = dims$N, ncol = dims$G),
        alpha = 0.1,
        Alpha = rep(alpha, dims$K),
        beta = 0.1,
        Beta = rep(beta, dims$K),
        a = 0.8,
        b = 0.8
) {
    if ("lambda_p" %in% names(Theta) & !("Lambda_p" %in% names(Theta))) {
        Theta$Lambda_p = matrix(Theta$lambda_p, nrow = dims$K, ncol = dims$N)
    }
    if ("lambda_e" %in% names(Theta) & !("Lambda_e" %in% names(Theta))) {
        Theta$Lambda_e = matrix(Theta$lambda_e, nrow = dims$N, ncol = dims$G)
    }
    if ("alpha" %in% names(Theta) & !("Alpha" %in% names(Theta))) {
        Theta$Alpha = rep(Theta$alpha, dims$K)
    }
    if ("beta" %in% names(Theta) & !("Beta" %in% names(Theta))) {
        Theta$Beta = rep(Theta$beta, dims$K)
    }
    fill_list(Theta, list(
        Lambda_p = Lambda_p,
        Lambda_e = Lambda_e,
        Alpha = Alpha,
        Beta = Beta,
        a = a,
        b = b
    ))
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
#' @export
set_gamma_prior_parameters <- function(
        Theta,
        dims,
        alpha_p = 10,
        Alpha_p = matrix(alpha_p, nrow = dims$K, ncol = dims$N),
        beta_p = sqrt(dims$N),
        Beta_p = matrix(beta_p, nrow = dims$K, ncol = dims$N),
        alpha_e = 10,
        Alpha_e = matrix(alpha_e, nrow = dims$N, ncol = dims$G),
        beta_e = sqrt(dims$N),
        Beta_e = matrix(beta_e, nrow = dims$N, ncol = dims$G),
        a = 0.8,
        b = 0.8
) {
    if ("alpha_p" %in% names(Theta) & !("Alpha_p" %in% names(Theta))) {
        Theta$Alpha_p = matrix(Theta$alpha_p, nrow = dims$K, ncol = dims$N)
    }
    if ("beta_p" %in% names(Theta) & !("Beta_p" %in% names(Theta))) {
        Theta$Beta_p = matrix(Theta$beta_p, nrow = dims$K, ncol = dims$N)
    }
    if ("alpha_e" %in% names(Theta) & !("Alpha_e" %in% names(Theta))) {
        Theta$Alpha_e = matrix(Theta$alpha_e, nrow = dims$N, ncol = dims$G)
    }
    if ("beta_e" %in% names(Theta) & !("Beta_e" %in% names(Theta))) {
        Theta$Beta_e = matrix(Theta$beta_e, nrow = dims$N, ncol = dims$G)
    }
    fill_list(Theta, list(
        Alpha_p = Alpha_p,
        Beta_p = Beta_p,
        Alpha_e = Alpha_e,
        Beta_e = Beta_e,
        a = a,
        b = b
    ))
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
        inits = NULL,
        fixed = NULL,
        prior_parameters = NULL
) {
    is_fixed = list(A = !learn_A)

    # prior parameters
    Theta = prior_parameters
    if (prior == 'truncnormal') {
        Theta = set_truncnorm_prior_parameters(Theta, dims)
    } else if (prior == 'exponential') {
        Theta = set_exponential_prior_parameters(Theta, dims)
    } else if (prior == 'gamma') {
        Theta = set_gamma_prior_parameters(Theta, dims)
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
