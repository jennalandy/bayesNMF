#' Run bayesNMF
#' @description Gibbs sampling for Bayesian Nonnegative Matrix Factorization with automatically learned rank and optional Metropolis-Hastings updates for computational efficiency
#' 
#' @param data matrix, counts matrix (rows are variables, columns are samples)
#' @param rank integer, number of signatures or range of ranks to consider
#' @param likelihood string, likelihood distribution for M, one of "poisson", "normal" (default "poisson")
#' @param prior string, prior distribution for P and E, one of "truncnormal", "exponential", "gamma" (default "truncnormal")
#' @param rank_method string, method to determine rank, one of "SBFI", "BFI", "BIC" (heuristic) (default "SBFI")
#' @param MH boolean, if likelihood is Poisson, whether to use Metropolis-Hastings updates described in the paper (default TRUE if likelihood is Poisson and prior is "truncnormal" or "exponential")
#' @param convergence_control list, convergence control parameters (default new_convergence_control())
#' @param prop_temp float, proportion of iterations to use for tempering if learning rank (default 0.2)
#' @param post_warmup integer, number of iterations post-warmup to sample from (if MH = TRUE) (default 2 * convergence_control$MAP_over)
#' @param output_dir string, output directory, will be created if it does not exist
#' @param overwrite boolean, whether to overwrite existing output, will delete existing output directory if it exists (default FALSE)
#' @param hyperprior_params list, fixed hyperprior parameters
#' @param init_prior_params list, initial prior parameters
#' @param init_params list, initial parameters
#' @param periodic_save boolean, whether to save the sampler object periodically (default TRUE, set to FALSE to reduce I/O and computational time)
#' @param save_all_samples boolean, whether to save all samples (default TRUE, set to FALSE to reduce memory usage)
#' 
#' @return bayesNMF_sampler object if a fixed rank is provided or rank_method is not "BIC", list of results, best rank, and sampler if "BIC"
#' 
#' @export
bayesNMF <- function(
  data, rank,
  likelihood = "poisson",
  prior = "truncnormal",
  rank_method = "SBFI",
  MH = likelihood == "poisson" & prior %in% c('truncnormal','exponential'),
  convergence_control = new_convergence_control(),
  prop_temp = 0.2,
  post_warmup = 2 * convergence_control$MAP_over,
  output_dir = paste0('nmf_', likelihood, '_', prior),
  overwrite = FALSE,
  hyperprior_params = list(),
  init_prior_params = list(),
  init_params = list(),
  periodic_save = TRUE,
  save_all_samples = TRUE
) {
  verbosity <- 1
  
  # initialize bayesNMF_sampler object
  sampler <- bayesNMF_sampler$new(
    data = data,
    rank = rank,
    likelihood = likelihood,
    prior = prior,
    rank_method = rank_method,
    MH = MH,
    convergence_control = convergence_control,
    prop_temp = prop_temp,
    post_warmup = post_warmup,
    output_dir = output_dir,
    overwrite = overwrite,
    hyperprior_params = hyperprior_params,
    init_prior_params = init_prior_params,
    init_params = init_params,
    verbosity = verbosity,
    periodic_save = periodic_save,
    save_all_samples = save_all_samples
  )
  if (sampler$specs$learning_rank) {
    # sampler only has rank_method if learning rank
    # special case for BIC, separate sampler for each rank, compare final BICs
    if (sampler$specs$rank_method == "BIC") {
      results <- lapply(rank, function(k) {
        output_dir_k <- file.path(output_dir, paste0("rank_", k))

        # initialize sampler for each rank
        sampler_k <- bayesNMF_sampler$new(
          data = data,
          rank = k,
          likelihood = likelihood,
          prior = prior,
          rank_method = rank_method,
          MH = MH,
          convergence_control = convergence_control,
          prop_temp = prop_temp,
          post_warmup = post_warmup,
          output_dir = output_dir_k,
          overwrite = overwrite,
          hyperprior_params = hyperprior_params,
          init_prior_params = init_prior_params,
          init_params = init_params,
          verbosity = verbosity,
          save_all_samples = save_all_samples,
          periodic_save = periodic_save
        )

        # Gibbs sample for each rank
        sampler_k$run_gibbs_sampler()
        
        # extract final BIC
        BIC <- sampler_k$state$MAP_metrics %>%
          dplyr::filter(iter == sampler_k$state$iter) %>%
          dplyr::pull(BIC)

        # return list of output directory, BIC, and time
        return(list(
          dir = output_dir_k,
          BIC = BIC,
          time = as.numeric(sampler_k$time$total)
        ))
      })

      # combine results into a data frame and find best rank
      results <- results %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        dplyr::arrange(BIC)

      best_rank <- results$rank[which.min(results$BIC)]
      print(glue::glue("Best rank: {best_rank}"))

      best_rank_dir <- results$dir[which.min(results$BIC)]
      sampler <- readRDS(file.path(best_rank_dir, "sampler.rds"))
      
      # return list of results, best rank, and sampler
      out <- list(
        results = results,
        best_rank = best_rank,
        sampler = sampler
      )
      saveRDS(sampler, file.path(output_dir, "sampler.rds"))
    } else if (!(sampler$specs$rank_method %in% c("SBFI", "BFI"))) {
      stop("Rank method must be SBFI, BFI, or BIC")
    } 
  }
  
  # standard Gibbs sampling for fixed rank, SBFI, or BFI
  sampler$run_gibbs_sampler()
  out <- sampler

  # return output (sampler object or list of results, best rank, and sampler if BIC)
  return(out)
}