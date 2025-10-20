#' Summary method for bayesNMF_sampler object
#'
#' @param sampler bayesNMF_sampler object
#' @param reference_P matrix, "cosmic", or NULL, reference signatures to align to
#'
#' @return data frame with columns:
#' \itemize{
#'   \item G: number of signatures
#'   \item N: number of samples
#'   \item K: number of mutations
#'   \item Signature: assigned reference signature
#'   \item Med_Contribution: median contribution of the signature to the samples among samples with at least 1 mutation attributed to the signature
#'   \item Prop_atleast_1: proportion of samples with at least 1 mutation attributed to the signature
#'   \item Reference_Signature: assigned reference signature (if reference_P is provided)
#'   \item Cosine_Similarity: cosine similarity between the estimated and assigned reference signatures (if reference_P is provided)
#' }
#' @export
summary.bayesNMF_sampler <- function(sampler, reference_P = "cosmic") {
  if ("character" %in% class(reference_P)) {
    if (reference_P == "cosmic") {
      reference_P = get_cosmic()
    } else {
      stop("Parameter `reference_P` must be a matrix or 'cosmic'")
    }
  }
  if (!is.null(reference_P)) {
    if (sampler$dims$K != nrow(reference_P)) {
      warning(glue::glue("Reference matrix has {nrow(reference_P)} rows, but data has {sampler$dims$K} rows. Setting `reference_P` to NULL."))
      reference_P <- NULL
    }
  }

  if (!is.null(reference_P)) {
    # only reassigns signatures if reference set has changed
    signature_assignments <- sampler$assign_signatures_ensemble(
      reference_P = reference_P
    )
    cosine_sim <- pairwise_sim(
      # use MAP$sig_idx when accessing the MAP
      sampler$MAP$P[,sampler$MAP$sig_idx, drop = FALSE],
      reference_P[,signature_assignments$assignments$sig_ref, drop = FALSE]
    )
  } else {
    signature_assignments <- list(
      'assignments' = data.frame(
        'sig_est' = seq_len(ncol(sampler$MAP$P)),
        'sig_ref' = seq_len(ncol(sampler$MAP$P))
      ),
      'votes' = NULL
    )
  }

  # only recomputes summary if reference set has changed
  if (!is.null(sampler$reference_comparison$summary)) {
    sampler_summary <- sampler$reference_comparison$summary
    return(sampler_summary)
  }

  # otherwise, for each assigned signature,
  # compute median contributions across samples
  sampler_summary <- lapply(
    1:nrow(signature_assignments$assignments), function(i) {
      sig_est <- signature_assignments$assignments$sig_est[i]
      sig_ref <- signature_assignments$assignments$sig_ref[i]
      atleast_1 <- sampler$MAP$E[sampler$MAP$sig_idx[sig_est], ] >= 1 # use MAP$sig_idx when accessing the MAP
      contribs <- sampler$MAP$E[sampler$MAP$sig_idx[sig_est], atleast_1]
      
      out <- data.frame(
        G = sampler$dims$G,
        N = sampler$dims$N,
        K = sampler$dims$K,
        Signature = sig_est,
        Med_Contribution = median(contribs),
        Prop_atleast_1 = mean(atleast_1)
      )
      if (!is.null(reference_P)) {
        out$Reference_Signature <- sig_ref
        cos_sim <- cosine_sim[sig_est, sig_ref]
        out$Cosine_Similarity <- cos_sim
      }
      return(out)
    }
  ) %>%
    do.call(rbind, .)

  # store summary in sampler object
  sampler$reference_comparison$summary <- sampler_summary

  return(sampler_summary)
}


#' For a named list of results, gets signature assignments, cosine similarities,
#' and distribution of mutations attributed to signatures in each.
#'
#' @param sampler_list Named list, containing one or more bayesNMF sampler objects. 
#'        Names will become identifiers in the data frame.
#' @param reference_P matrix, "cosmic" , or NULL, reference signatures to align to
#'
#' @return data frame with columns:
#' \itemize{
#'   \item Name: name of the sampler
#'   \item G: number of signatures
#'   \item N: number of samples
#'   \item K: number of mutations
#'   \item Signature: assigned reference signature
#'   \item Med_Contribution: median contribution of the signature to the samples among samples with at least 1 mutation attributed to the signature
#'   \item Prop_atleast_1: proportion of samples with at least 1 mutation attributed to the signature
#'   \item Reference_Signature: assigned reference signature (if reference_P is provided)
#'   \item Cosine_Similarity: cosine similarity between the estimated and assigned reference signatures (if reference_P is provided)
#' }
#' @export
summarize_samplers <- function(sampler_list, reference_P = "cosmic") {
  if ("character" %in% class(reference_P)) {
    if (reference_P == "cosmic") {
      reference_P = get_cosmic()
    } else {
      stop("Parameter `reference_P` must be a matrix or 'cosmic'")
    }
  }

  all_rows <- sapply(sampler_list, function(sampler) {
    sampler$dims$K
  })
  if (length(unique(all_rows)) > 1) {
    stop("All samplers must have the same number of variables.")
  }

  if (!is.null(reference_P)) {
    if (all_rows[1] != nrow(reference_P)) {
      warning(glue::glue("Reference matrix has {nrow(reference_P)} rows, but data has {all_rows[1]} rows. Setting `reference_P` to NULL."))
      reference_P <- NULL
    }
  }

  sampler_list_summary <- lapply(names(sampler_list), function(sampler_name) {
    sampler <- sampler_list[[sampler_name]]
    if (is.null(sampler$state$converged)) {
      print(paste('not done:', sampler_name))
      return(NULL)
    }

    sampler_summary <- summary(sampler, reference_P)
    sampler_summary$Name <- paste0(sampler_name, ' (', sampler$dims$G, ')')
    return(sampler_summary)

  }) %>%
  do.call(rbind, .)

  return(sampler_list_summary)
}

# These functions are copied over to be methods of the bayesNMF_sampler class

####################################
###### PUBLIC METHODS ##############
####################################

#' Assign signatures for each posterior sample based on cosine similarity and ensemble with majority vote
#'
#' @param self bayesNMF sampler object
#' @param reference_P matrix, "cosmic", or NULL, reference signatures to align to
#' @param idxs vector of indices to consider (default "MAP_idx" indicates indices used to compute most recent MAP)
#' @param credible_interval numeric, credible interval to compute for cosine similarities (default 0.95)
#'
#' @return list of two data frames:
#' \itemize{
#'   \item assignments: holds final assignments for each signature, cosine similarities between MAP estimates and reference signatures, and 95% credible intervals for the cosine similarities. 
#'   \item votes: holds the proportion of votes for each signature assignment for each posterior sample. This includes signatures that did not receive majority vote, allowing users to understand posterior uncertainty in assignment.
#' }
#' and stores the result in self$reference_comparison
#' 
#' @noRd
assign_signatures_ensemble_ <- function(
  self,
  reference_P = "cosmic",
  idxs = "MAP_idx",
  credible_interval = 0.95
) {
  if ("character" %in% class(reference_P)) {
    if (reference_P == "cosmic") {
      reference_P = get_cosmic()
    } else {
      stop("Parameter `reference_P` must be a matrix or 'cosmic'")
    }
  }

  if (!is.null(reference_P)) {
    if (self$dims$K != nrow(reference_P)) {
      warning(glue::glue("Reference matrix has {nrow(reference_P)} rows, but data has {self$dims$K} rows. Setting `reference_P` to NULL."))
      reference_P <- NULL
    }
  }
  
  if ("character" %in% class(idxs)) {
    if (idxs == "MAP_idx") {
      idxs = self$MAP$idx
    } else {
      stop("Parameter `idxs` must be a vector of indices or 'MAP_idx'")
    }
  }

  # only assign included signatures, even if user ran get_MAP with final = FALSE
  # - keep_sigs is the set of all indices that are included (with respect to theoriginal rank vector)
  #   always use keep_sigs when accessing individual samples
  # - MAP$sig_idx is the set of MAP indices that are included
  #   always use MAP$sig_idx when accessing the MAP estimate
  #   - if get_MAP was run with final = FALSE, MAP$sig_idx will be the same as keep_sigs
  #   - if get_MAP was run with final = TRUE, MAP$sig_idx will just be 1:<final rank>

  if (length(self$MAP$keep_sigs) == self$dims$N & sum(self$MAP$A[1,] == 0) > 0) {
    # this is only true if final = FALSE was used in get_MAP (or all signatures are included)
    keep_sigs <- which(self$MAP$A[1,] == 1)
    self$MAP$sig_idx <- keep_sigs
  } else {
    # otherwise, final = True was used in get_MAP
    keep_sigs <- self$MAP$keep_sigs
    self$MAP$sig_idx <- 1:length(keep_sigs)
  }


  # check if reference has already been assigned
  # over the same indices and with the same MAP estimate
  if (
    check_matrix_equal(reference_P, self$reference_comparison$reference_P) &
    check_vector_equal(idxs, self$reference_comparison$idxs) &
    check_vector_equal(keep_sigs, self$reference_comparison$keep_sigs)
  ) {
    out <- list(
      'assignments' = self$reference_comparison$assignments,
      'votes' = self$reference_comparison$votes
    )
    return(out)
  } else {
    # reset all reference comparison data
    self$reference_comparison <- list(
      'reference_P' = reference_P,
      'MAP_idx' = idxs,
      'keep_sigs' = keep_sigs,
      'assignments' = NULL,
      'votes' = NULL,
      'summary' = NULL,
      'plots' = list(),
      'label_switching_df' = NULL
    )
  }
  # otherwise, assign signatures and store the result in self$reference_comparison

  # across idxs, assign signatures and use cosine similarity as voting weight
  votes <- lapply(idxs, function(idx) {
    assignment_res <- hungarian_assignment(
      self$samples$P[[idx]][, keep_sigs, drop = FALSE], # use keep_sigs when accessing individual samples
      reference_P
    )
    assignment_res$i = idx
    return(assignment_res)
  }) %>%
    do.call(rbind, .) %>%
    dplyr::group_by(sig_est, sig_ref) %>%
    dplyr::summarize(
      votes = sum(cos_sim),
      .groups = "keep"
    ) %>%
    dplyr::ungroup() 

  votes <- votes %>%
    merge(
      votes %>%
        dplyr::group_by(sig_est) %>%
        dplyr::summarize(total_votes = sum(votes))
    ) %>% 
    dplyr::mutate(
      prop_votes = votes / total_votes
    ) %>%
    dplyr::select(sig_est, sig_ref, prop_votes) %>%
    dplyr::arrange(sig_est, -prop_votes)

  assignments <- votes %>%
    dplyr::group_by(sig_est) %>%
    dplyr::summarize(
      sig_ref = sig_ref[which.max(prop_votes)],
      prop_votes = max(prop_votes),
      .groups = "keep"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(sig_est) %>%
    dplyr::select(-prop_votes)

  assignments$MAP_cosine <- diag(pairwise_sim(
    self$MAP$P[, self$MAP$sig_idx, drop = FALSE], # use MAP$sig_idx when accessing the MAP
    reference_P[, assignments$sig_ref, drop = FALSE]
  ))

  sample_cosines <- lapply(idxs, function(idx) {
    sim <- pairwise_sim(
      self$samples$P[[idx]][, keep_sigs, drop = FALSE], # use keep_sigs when accessing individual samples
      reference_P[,assignments$sig_ref, drop = FALSE]
    )
    sapply(1:nrow(assignments), function(i) {
      sim[assignments$sig_est[i], assignments$sig_ref[i]]
    })
  }) %>%
    do.call(rbind, .)
  
  sample_cosine_quantiles <- apply(sample_cosines, 2, function(col) {
    quantile(col, c((1 - credible_interval)/2, 1-(1 - credible_interval)/2))
  })

  assignments$lower_cosine <- sample_cosine_quantiles[1,]
  assignments$upper_cosine <- sample_cosine_quantiles[2,]

  out <- list(
    'assignments' = assignments,
    'votes' = votes
  )
  self$reference_comparison$reference_P <- reference_P
  self$reference_comparison$idxs <- idxs
  self$reference_comparison$assignments <- assignments
  self$reference_comparison$votes <- votes
  
  return(out)
}