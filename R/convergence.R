#' Create new convergence control object
#' @description Specify convergence criteria for the bayesNMF Gibbs sampler.
#'
#' @param MAP_over integer, number of samples to average over for MAP
#' @param MAP_every integer, how often (in samples) to compute MAP for evaluation
#' @param tol numeric, tolerance for what is a meaningful "change" in MAP metric (default 0.001 means 0.1% change)
#' @param Ninarow_nochange integer, number of consecutive MAPs with no change before declaring convergence (default 5)
#' @param Ninarow_nobest integer, number of consecutive MAPs with no new "best" before declaring convergence (default 10)
#' @param miniters integer, minimum number of samples to consider for convergence (default 1000)
#' @param maxiters integer, absolute maximum number of samples to explore (default 5000)
#' @param minA integer, minimum number of counts A matrix must have to consider it a valid MAP (default 0)
#' @param metric string, one of c('loglikelihood','logposterior','RMSE','KL') (default "logposterior")
#'
#' @return list of convergence control parameters
#' @export
new_convergence_control <- function(
  MAP_over = 1000,
  MAP_every = 100,
  tol = 0.001,
  Ninarow_nochange = 5,
  Ninarow_nobest = 10,
  miniters = 1000,
  maxiters = 5000,
  minA = 0,
  metric = "logposterior"
) {
  # check for invalid inputs
  if (miniters >= maxiters) {
    warning("miniters >= maxiters, setting miniters to 0.")
    miniters <- 0
  }

  # return list of convergence control parameters
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

# These functions are copied over to be methods of the bayesNMF_sampler class

####################################
###### PRIVATE METHODS ##############
####################################

#' Check for convergence
#' @param self bayesNMF_sampler object
#' @param private list of private methods of bayesNMF_sampler object
#' @param final boolean, if TRUE, subset to only included signatures
#' 
#' @return string, message indicating convergence status
#' @noRd
check_convergence_ <- function(self, private, final = FALSE) {
  convergence_control <- self$specs$convergence_control

  # compute new MAP metrics
  private$update_MAP_metrics(final = final)

  # extract MAP metric of interest for current iteration
  MAP_metric <- self$state$MAP_metrics[
    self$state$MAP_metrics$iter == self$state$iter,
    self$specs$convergence_control$metric
  ]

  # log likelihood or log posterior should be maximized not minimized
  # so we need to invert the sign of the metric
  if (
    self$specs$convergence_control$metric %in% 
    c("loglikelihood", "logposterior")
  ) {
    MAP_metric <- -1 * MAP_metric
  }

  # for first one, force % change < 0
  if (!("prev_MAP_metric" %in% names(self$state))) {
    self$state$prev_MAP_metric <- MAP_metric + 1
    self$state$best_MAP_metric <- MAP_metric + 1
    self$state$inarow_na <- 0
    self$state$inarow_no_change <- 0
    self$state$inarow_no_best <- 0
  }
  percent_change <- (
    MAP_metric - self$state$prev_MAP_metric
  ) / self$state$prev_MAP_metric
  self$state$prev_percent_change <- percent_change
  self$state$prev_MAP_metric <- MAP_metric
  if (is.na(percent_change)) {
    # if NA percent change, reset convergence
    self$state$inarow_no_change <- 0
    self$state$inarow_no_best <- 0
    self$state$inarow_na <- self$state$inarow_na + 1
  } else if (abs(percent_change) < convergence_control$tol) {
    # if |change| < tol, then there was "no change"
    self$state$inarow_no_change <- self$state$inarow_no_change + 1
    self$state$inarow_na <- 0
  } else {
    # otherwise |change| > tol, so there was a change
    self$state$inarow_no_change <- 0
    self$state$inarow_na <- 0
  }

  # qualifies if 
  #     (temperature == 1 for all samples considered) 
  #     AND (iter is at least miniters)
  if (
    all(self$temperature_schedule[
        # only equivalent to self$MAP$idx if self$specs$save_all_samples is TRUE
        (self$state$iter - self$specs$convergence_control$MAP_over):(self$state$iter)
      ] == 1) &
      self$state$iter >= convergence_control$miniters
  ) {
    # if iter is the new best, store info
    if (MAP_metric < self$state$best_MAP_metric) {
      self$state$best_MAP_metric <- MAP_metric
      self$state$best_iter <- self$state$iter
      self$state$inarow_no_best <- 0
    } else {
      self$state$inarow_no_best <- self$state$inarow_no_best + 1
    }

    # converged if
    #   (no change for Ninarow_nochange)
    #   OR (no best for Ninarow_nobest)
    #   OR (iter hit maxiters)
    if (self$state$inarow_no_change >= convergence_control$Ninarow_nochange) {
      self$state$converged <- TRUE
      self$state$why <- "no change"
    } else if (self$state$inarow_no_best >= convergence_control$Ninarow_nobest) {
      self$state$converged <- TRUE
      self$state$why <- "no best"
    } else if (self$state$iter >= convergence_control$maxiters) {
      self$state$converged <- TRUE
      self$state$why <- "max iters"
    }
  }

  # construct message indicating convergence status
  flip <- if (self$specs$convergence_control$metric %in% c("loglikelihood", "logposterior")) -1 else 1
  msg <- glue::glue(paste0(
    "{self$specs$convergence_control$metric} = {round(MAP_metric, 2)} | ",
    "{flip * round(percent_change*100, 2)}% change | ",
    "{self$state$inarow_no_change} no change | ",
    "{self$state$inarow_no_best} no best | ",
    "{self$state$inarow_na} NA"
  ))
  return(msg)
}