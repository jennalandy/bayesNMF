#' Get Single Recovery Prior
#' @description Builds recover prior for a single signature
#'
#' @param Pn vector, K
#' @param likelihood string, one of c('normal','poisson')
#' @param prior string, one of c('truncnormal','exponential','gamma')
#'
#' @return list, prior parameters corresponding to given signature
#' @noRd
get_one_recovery_prior <- function(Pn, likelihood, prior) {
    # simulate data from signature Pn
    K <- length(Pn)
    E <- matrix(500, nrow = 1, ncol = 500)
    M <- matrix(nrow = K, ncol = 500)
    mean = Pn %*% E
    if (likelihood == "poisson") {
        for (k in 1:K) {
            for (g in 1:500) {
                M[k, g] <- rpois(1, mean[k, g])
            }
        }
    } else if (likelihood == "normal") {
        for (k in 1:K) {
            for (g in 1:500) {
                M[k, g] <- rnorm(1, mean[k, g], 2)
            }
        }
        M[M <= 0] <- 0
    }

    # run bayesNMF for 1 signature
    res <- bayesNMF(
        M, N = 1,
        likelihood = likelihood,
        prior = prior,
        file = tempfile(),
        store_logs = TRUE
    )
    # center P at 10
    res$MAP$P = 10 * K * res$MAP$P / sum(res$MAP$P)

    # set hyperpriors to center prior around MAP estimate
    start = max((res$converged_at - 1000), 1)
    keep = start:res$converged_at
    if (prior == 'truncnormal') {
        Mu_p = sapply(1:K, function(k) {
            uniroot(function(x) {
                x + sqrt(x) * dnorm(-sqrt(x))/(pnorm(-sqrt(x)) - 1) - res$MAP$P[k]
            }, interval = c(0, 1000))$root
        })
        prior_parameters = list(
            Mu_p = Mu_p,
            Sigmasq_p = Mu_p/10
        )
    } else if (prior == 'exponential') {
        prior_parameters = list(
            Lambda_p = 1/res$MAP$P
        )
    } else if (prior == 'gamma') {
        prior_parameters = list(
            Alpha_p = res$MAP$P * 10,
            Beta_p = rep(10, K)
        )
    }

    return(prior_parameters)
}


#' Get Recovery Prior Parameters
#' @description Builds recovery priors for a given signatures matrix `P`.
#' For each signature in `P`, data is simulated from that P alone and
#' `bayesNMF` is run with `N = 1`. Recovery priors are centered around the
#' estimated maximum a-posteriori learned by `bayesNMF`.
#'
#' @param P matrix, K x N
#' @param likelihood string, one of c('normal','poisson')
#' @param prior string, one of c('truncnormal','exponential','gamma')
#' @param verbose boolean
#' @param savefile string, path to save recovery priors to
#'
#' @return list, prior parameters corresponding to given signatures matrix
#' @export
get_recovery_priors <- function(
    P, likelihood, prior,
    verbose = FALSE,
    savefile = file.path(
        "recovery_priors", paste0(likelihood, "_", prior, ".rds")
    )
) {
    for (n in 1:ncol(P)) {
        prior_parameters = get_one_recovery_prior(P[,n], likelihood, prior)
        if (n == 1) {
            all_prior_parameters = prior_parameters
        } else {
            for (name in names(prior_parameters)) {
                all_prior_parameters[[name]] = cbind(
                    all_prior_parameters[[name]],
                    prior_parameters[[name]]
                )
            }
        }
        if (verbose) {
            print(paste(n, "/", ncol(P)))
        }
    }

    all_prior_parameters$N_r <- ncol(P)
    all_prior_parameters$P <- P

    saveRDS(all_prior_parameters, file = savefile)
    return(all_prior_parameters)
}
