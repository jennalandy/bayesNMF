#' Create new convergence control object
#' @description Specify convergence criteria for the bayesNMF Gibbs sampler.
#'
#' @param MAP_over integer, number of samples to average over for MAP
#' @param MAP_every integer, how often (in samples) to compute MAP for evaluation
#' @param tol numeric, tolerance for what is a meaninful "percent change"
#' in MAP metric
#' @param Ninarow_change integer, convergence may be determined by the number
#' of MAPs in a row with no change
#' @param Ninarow_nobest integer, convergence may be determined by the number
#' of MAPs in a row with no new "best"
#' @param maxiters integer, absolute maximum number of samples to explore
#' @param minA integer, minimum number of counts A matrix must have to consider it a valid MAP
#' @param metric string, one of c('loglikelihood','logposterior','RMSE','KL')
#'
#' @return list
#' @export
new_convergence_control <- function(
    MAP_over = 500,
    MAP_every = 100,
    tol = 0.001,
    Ninarow_nochange = 10,
    Ninarow_nobest = 20,
    miniters = 500,
    maxiters = 2000,
    minA = 0,
    metric = "logposterior"
) {
    if (miniters >= maxiters) {
        warning("miniters >= maxiters, setting miniters to 0.")
        miniters = 0
    }
    list(
        MAP_over = MAP_over,
        MAP_every = MAP_every,
        tol = tol,
        Ninarow_nochange = Ninarow_nochange,
        Ninarow_nobest = Ninarow_nobest,
        miniters = miniters,
        maxiters = maxiters,
        minA = minA,
        metric = metric
    )
}

#' Get metric
#'
#' @param metric string, one of c('loglikelihood','logposterior',RMSE','KL')
#' @param Mhat matrix, reconstructed mutational catalog
#' @param M matrix, true mutational catalog
#' @param Theta list, current values of all unkowns
#' @param likelihood string, one of c("normal", "poisson")
#' @param prior string, one of c("truncnormal","exponential","gamma")
#' @param dims list, named list of dimensions
#' @param logfac vector
#'
#' @return scalar
#' @noRd
get_metric <- function(
    metric, Mhat, M,
    Theta = NULL,
    likelihood = NULL,
    prior = NULL,
    dims = NULL,
    logfac = NULL
) {
    if (metric %in% c('loglikelihood', 'logposterior', 'BIC')) {
        if (likelihood == 'normal') {
            loglik = get_loglik_normal(M, Theta, dims)
        } else if (likelihood == 'poisson') {
            loglik = get_loglik_poisson(M, Theta, dims, logfac)
        }
        if (metric == 'loglikelihood') {
            return(-1 * loglik)
        } else if (metric == 'BIC') {
            return(get_BIC(loglik, Theta, dims, likelihood, prior))
        }
        logpost = loglik + get_logprior(Theta, likelihood, prior)
        return(-1 * logpost)
    } else if (metric == 'RMSE') {
        get_RMSE(M, Mhat)
    } else if (metric == 'KL') {
        get_KLDiv(M, Mhat)
    }
}

#' Check whether Gibb's sampler has converged
#'
#' @param iter integer, iteration/sample
#' @param gamma numeric, tempering parameter
#' @param Mhat matrix, reconstructed mutational catalog
#' @param M matrix, true mutational catalog
#' @param Theta list, current values of all unkowns
#' @param convergence_status list, current status of convergence
#' @param convergence_control list, control parameters
#' @param first_MAP boolean, whether this is the first MAP computed
#' @param metric string, one of c('loglikelihood','RMSE','KL')
#' @param likelihood string, one of c("normal", "poisson")
#' @param prior string, one of c("truncnormal","exponential","gamma")
#' @param dims list, named list of dimensions
#' @param logfac vector
#'
#' @return list, updated status of convergence
#' @noRd
check_converged <- function(
    iter, gamma, Mhat, M,
    convergence_status,
    convergence_control,
    first_MAP,
    Theta = NULL,
    likelihood = NULL,
    prior = NULL,
    dims = NULL,
    logfac = NULL
) {
    MAP_metric = get_metric(
        convergence_control$metric, Mhat, M,
        Theta, likelihood, prior, dims, logfac
    )

    # for first one, force % change < 0
    if (first_MAP) {
        convergence_status$prev_MAP_metric = MAP_metric + 1
        convergence_status$best_MAP_metric = MAP_metric + 1
        convergence_status$inarow_na = 0
        convergence_status$inarow_no_change = 0
        convergence_status$inarow_no_best = 0
    }

    percent_change = (
        MAP_metric - convergence_status$prev_MAP_metric
    )/convergence_status$prev_MAP_metric
    convergence_status$prev_percent_change = percent_change
    convergence_status$prev_MAP_metric = MAP_metric

    # if NA percent change, no way it's converged
    if (is.na(percent_change)) {
        convergence_status$inarow_no_change = 0
        convergence_status$inarow_no_best = 0
        convergence_status$inarow_na = convergence_status$inarow_na + 1
        return(convergence_status)
    }

    # if |change| < tol, then there was "no change"
    if (abs(percent_change) < convergence_control$tol) {
        convergence_status$inarow_no_change =
            convergence_status$inarow_no_change + 1
    } else {
        convergence_status$inarow_no_change = 0
    }

    # if this is new best
    if (MAP_metric < convergence_status$best_MAP_metric) {
        convergence_status$best_MAP_metric = MAP_metric
        convergence_status$best_iter = iter
        convergence_status$inarow_no_best = 0
    } else {
        convergence_status$inarow_no_best =
            convergence_status$inarow_no_best + 1
    }

    # stop if
    # gamma == 1 AND iter is at least miniters
    # AND
    #   (no change for Ninarow)
    #   OR (no best for Ninarow)
    #   OR (change is less than mintol)
    #   OR (iter hit maxiters)
    if (gamma == 1 & iter > convergence_control$miniters) {
        if (convergence_status$inarow_no_change >=
            convergence_control$Ninarow_nochange
        ) {
            convergence_status$converged = TRUE
            convergence_status$why = "no change"
        } else if (convergence_status$inarow_no_best >=
                   convergence_control$Ninarow_nobest
        ) {
            convergence_status$converged = TRUE
            convergence_status$why = "no best"
        } else if (iter >= convergence_control$maxiters) {
            convergence_status$converged = TRUE
            convergence_status$why = "max iters"
        }
    }


    return(convergence_status)
}
