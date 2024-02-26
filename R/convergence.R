#' Create new convergence control object
#'
#' @param MAP_over integer, number of samples to average over for MAP
#' @param MAP_every integer, how often (in samples) to compute MAP for evaluation
#' @param tol numeric, tolerance for what is a meaninful "percent change"
#' in MAP RMSE
#' @param Ninarow_change integer, convergence may be determined by the number
#' of MAPs in a row with no change
#' @param Ninarow_nobest integer, convergence may be determined by the number
#' of MAPs in a row with no new "best"
#' @param maxiters integer, absolute maximum number of samples to explore
#'
#' @return list
#' @export
new_convergence_control <- function(
    MAP_over = 1000,
    MAP_every = 100,
    tol = 0.001,
    Ninarow_nochange = 10,
    Ninarow_nobest = 20,
    maxiters = 10000
) {
    list(
        MAP_over = MAP_over,
        MAP_every = MAP_every,
        tol = tol,
        Ninarow_nochange = Ninarow_nochange,
        Ninarow_nobest = Ninarow_nobest,
        maxiters = maxiters
    )
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
#'
#' @return list, updated status of convergence
#' @noRd
check_converged <- function(
    iter, gamma, Mhat, M, Theta,
    convergence_status,
    convergence_control,
    first_MAP
) {
    MAP_RMSE = sqrt(mean((M - Mhat) ** 2))
    # for first one, force % change < 0
    if (first_MAP) {
        convergence_status$prev_MAP_RMSE = MAP_RMSE + 1
        convergence_status$best_MAP_RMSE = MAP_RMSE + 1
    }

    percent_change = (
        MAP_RMSE - convergence_status$prev_MAP_RMSE
    )/convergence_status$prev_MAP_RMSE
    convergence_status$prev_percent_change = percent_change
    convergence_status$prev_MAP_RMSE = MAP_RMSE

    # if |change| < tol, then there was "no change"
    if (abs(percent_change) < convergence_control$tol) {
        convergence_status$inarow_no_change =
            convergence_status$inarow_no_change + 1
    } else {
        convergence_status$inarow_no_change = 0
    }

    # if this is new best
    if (MAP_RMSE < convergence_status$best_MAP_RMSE) {
        convergence_status$best_MAP_RMSE = MAP_RMSE
        convergence_status$best_iter = iter
        convergence_status$best_Theta = Theta
        convergence_status$inarow_no_best = 0
    } else {
        convergence_status$inarow_no_best =
            convergence_status$inarow_no_best + 1
    }

    # stop if
    # (no change for Ninarow or no best for Ninarow)
    # OR (change is less than mintol AND iter is at least 1000)
    if (
        gamma == 1 &
        convergence_status$inarow_no_change >= convergence_control$Ninarow_nochange
    ) {
        convergence_status$converged = TRUE
        convergence_status$burn_in = iter
        convergence_status$why = "no change"
    } else if (
        gamma == 1 &
        convergence_status$inarow_no_best >= convergence_control$Ninarow_nobest &
        iter > 5000
    ) {
        convergence_status$converged = TRUE
        convergence_status$burn_in = convergence_status$best_iter
        convergence_status$why = "no best"
    } else if (
        iter >= convergence_control$maxiters
    ) {
        convergence_status$converged = TRUE
        convergence_status$burn_in = convergence_status$best_iter
        convergence_status$why = "max iters"
    }

    return(convergence_status)
}
