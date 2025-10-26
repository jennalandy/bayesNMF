#' Plot a bayesNMF_sampler object
#'
#' @param sampler bayesNMF_sampler object
#' @param sigs boolean, whether to plot signatures
#' @param reference_P matrix or "cosmic", reference signatures to align to
#' @param combine_below numeric, threshold for combining signatures into "other"
#' @param verbose boolean, whether to print verbose output listing the saved plots and file locations
#' @param update_saved_object boolean, whether to update the bayesNMF_sampler object with stored plots and assignments
#'
#' @return named list of ggplot2 objects
#' @export
plot.bayesNMF_sampler <- function(
  sampler,
  sigs = FALSE,
  reference_P = "cosmic",
  combine_below = 0.8,
  verbose = TRUE,
  update_saved_object = TRUE
) {
  options(bitmapType = "cairo")
  if ('character' %in% class(reference_P)) {
    if (reference_P == 'cosmic') {
      reference_P = get_cosmic()
    } else {
      stop("Parameter `reference_P` must be a matrix, 'cosmic', or NULL")
    }
  }

  if (!is.null(reference_P)) {
    if (sampler$dims$K != nrow(reference_P)) {
      warning(glue::glue("Reference matrix has {nrow(reference_P)} rows, but data has {sampler$dims$K} rows. Setting `reference_P` to NULL."))
      reference_P <- NULL
    }
    # only reassigns signatures if reference set has changed
    signature_assignments <- sampler$assign_signatures_ensemble(reference_P = reference_P)
    reference_P <- sampler$reference_comparison$reference_P
  }

  # summary plot
  summary_plot <- plot_summary(sampler, reference_P = reference_P)
  summary_filename <- file.path(sampler$specs$output_dir, "summary.png")
  sampler$reference_comparison$plots$summary <- summary_plot
  ggplot2::ggsave(
    filename = summary_filename,
    plot = summary_plot,
    # type = "cairo",
    height = 20, width = 7
  )
  if (verbose) {
    print(glue::glue("{summary_filename}: Summary plot"))
  }

  if (!is.null(reference_P)) {
    # similarity heatmap plot
    # use assignments based on ensemble
    similarity_heatmap_plot <- plot_similarity_heatmap(
      sampler$MAP$P[,sampler$MAP$sig_idx, drop = FALSE], # use MAP$sig_idx when accessing the MAP
      reference_P = reference_P[,signature_assignments$assignments$sig_ref, drop = FALSE]
    )
    similarity_heatmap_filename <- file.path(
      sampler$specs$output_dir, "similarity_heatmap.png"
    )
    sampler$reference_comparison$plots$similarity_heatmap <- similarity_heatmap_plot
    ggplot2::ggsave(
      filename = similarity_heatmap_filename,
      plot = similarity_heatmap_plot,
      # type = "cairo",
      height = 10, width = 10
    )
    if (verbose) {
      print(glue::glue(
        "{similarity_heatmap_filename}: Similarity heatmap"
      ))
    }

    # label switching plot
    if (sampler$specs$learning_rank && sampler$specs$save_all_samples) {
      label_switching_plot <- plot_label_switching(
        sampler, reference_P = reference_P, combine_below = combine_below
      )
      label_switching_filename <- file.path(
        sampler$specs$output_dir, "label_switching.png"
      )
      sampler$reference_comparison$plots$label_switching <- label_switching_plot
      ggplot2::ggsave(
        filename = label_switching_filename,
        plot = label_switching_plot,
        # type = "cairo",
        height = min(sampler$dims$N*3, 20),
        width = min(sampler$state$iter/10 + 3, 40),
        limitsize = FALSE
      )
      if (verbose) {
        print(glue::glue(
          "{label_switching_filename}: Label switching dignostic plot"
        ))
      }
    }
  }

  # signature distribution plot
  signature_dist_plot <- plot_signature_dist(sampler, reference_P = reference_P)
  signature_dist_filename <- file.path(
    sampler$specs$output_dir, "signature_dist.png"
  )
  sampler$reference_comparison$plots$signature_dist <- signature_dist_plot
  ggplot2::ggsave(
    filename = signature_dist_filename,
    plot = signature_dist_plot,
    # type = "cairo",
    height = 5, width = 10
  )
  if (verbose) {
    print(glue::glue(
      "{signature_dist_filename}: Signature distribution plot"
    ))
  }

  # signature plots
  if (sigs) {
    for (sig in 1:sum(sampler$MAP$A)) {
      if (!is.null(reference_P)) {
        ref = signature_assignments$assignments$sig_ref[
          signature_assignments$assignments$sig_est == sig
        ]
        title = paste0("Signature ", sig, ", Reference ", ref, "")
      } else {
        ref = NULL
        title = paste0("Signature ", sig)
      }

      sig_plot <- plot_sig(
        sampler, sig, ref = ref,
        reference_P = reference_P,
        title = title
      )
      sig_filename <- file.path(
        sampler$specs$output_dir, paste0("sig_", sig, ".png")
      )
      sampler$reference_comparison$plots[[paste0("sig_", sig)]] <- sig_plot
      ggplot2::ggsave(
        filename = sig_filename,
        plot = sig_plot,
        # type = "cairo",
        height = 3, width = 10
      )
    }
    if (verbose) {
      print(glue::glue("{sig_filename}: Signature plots"))
    }
  }
  # update saved object after storing samples, label_switching_df, assignments, etc.
  if (update_saved_object) {
    sampler$save_object()
  }

  return(sampler$reference_comparison$plots)
}

#' Plot heatmap of cosine similarities
#' @description Plot a heatmap of cosine similarities between two matrices
#' with `ggplot2`. Can be similarity between rows or columns with the `which`
#' parameter.
#'
#' @param estimated_P matrix, estimated P (factors matrix)
#' @param reference_P matrix or "cosmic", reference factors to align to
#' @param est_names names of estimated factors
#' @param ref_names names of reference factors
#' @param which string, one of c("rows","cols")
#' @param keep_all_est boolean, whether to keep all estimated factors (versus only assigned estimated factors)
#' @param keep_all_ref boolean, whether to keep all reference factors (versus only assigned reference factors)
#'
#' @return ggplot object
#' @export
plot_similarity_heatmap <- function(
  estimated_P, reference_P = "cosmic",
  est_names = NULL,
  ref_names = NULL,
  which = "cols",
  keep_all_est = TRUE,
  keep_all_ref = FALSE
) {
  if ("character" %in% class(reference_P)) {
    if (reference_P == "cosmic") {
      reference_P = get_cosmic()
    } else {
      stop("Parameter `reference_P` must be a matrix or 'cosmic'")
    }
  }

  if (!is.null(reference_P)) {
    if (nrow(estimated_P) != nrow(reference_P)) {
      warning(glue::glue("Reference matrix has {nrow(reference_P)} rows, but data has {nrow(estimated_P)} rows. Setting `reference_P` to NULL."))
      reference_P <- NULL
    }
  }

  sim_mat <- hungarian_assignment(
    estimated_P =  estimated_P,
    reference_P = reference_P,
    which = which,
    keep_all_est = keep_all_est,
    keep_all_ref = keep_all_ref,
    return_mat = TRUE
  )

  sim_mat_melted <- reshape2::melt(sim_mat)

  # text needs to be smaller for more signatures, don't go over 10
  text_size = min(50/ncol(estimated_P), 10)
  text_size = max(text_size, 1) # or under 1
  heatmap <- sim_mat_melted %>%
    dplyr::mutate(
      Var1 = factor(
        sim_mat_melted$Var1,
        levels = unique(sim_mat_melted$Var1)
      ),
      Var2 = factor(
        sim_mat_melted$Var2,
        levels = unique(sim_mat_melted$Var2)
      )
    ) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = Var1, y = Var2,
      fill = value, label = round(value, 2)
    )) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(size = text_size) +
    ggplot2::labs(
      x = "Estimated Signatures",
      y = "Reference Signatures",
      fill = "Cosine\nSimilarity"
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 30)
    ) +
    ggplot2::scale_fill_gradient(
      limits = c(0, 1),
      low = "white", high = "steelblue"
    )

  return(heatmap)
}

#' Plot a single mutational signature.
#' 
#' @description
#' This function can be used to plot a reference mutational signature, an
#' estimated mutational signature with posterior uncertainty as points with
#' error bars, or both. The default reference is COSMIC, and unless otherwise
#' specified, the reference signature chosen to plot is the one assigned using
#' our posterior ensemble assignment procedure.
#'
#' One of the following must be true
#' - sampler and sig are specified
#' - reference_P is specified and ref is integer, a column name, or one of "assigned", "best"
#' - ref is a vector
#' 
#' Often ref = "assigned" is the same as ref = "best", but not necessarily as the assignment
#' procedure assumes only one estimated signature can be assigned to a reference signature, and
#' vise versa.
#'
#' @param sampler bayesNMF_sampler object, results list
#' @param sig integer, index of the estimated signature to plot
#' @param ref integer, column name string, or one of "assigned", "best" (default "assigned")
#' @param reference_P matrix or "cosmic", reference factors for comparison
#' @param title string
#' @param cosine boolean, whether to report cosine similarity between estimated
#' and reference. Ignored if either of sig, ref are NULL.
#'
#' @return ggplot2 object
#' @export
plot_sig <- function(
  sampler = NULL, sig = NULL,
  ref = "assigned", reference_P = "cosmic",
  title = "", cosine = TRUE
) {
  if (is.null(sampler) & is.null(ref)) {
    stop("Must provide either `sampler` (results of bayesNMF) or `ref` (vector for reference signature or column name of `reference_P`).")
  } else if (is.null(ref) & is.null(sampler) & is.null(sig)) {
    stop("To plot a signature from `sampler`, must specify signature index with `sig`.")
  }

  if ('character' %in% class(reference_P)) {
    if (reference_P == 'cosmic') {
      reference_P = get_cosmic()
    } else {
      stop("Parameter `reference_P` must be a matrix or the string 'cosmic'.")
    }
  }

  if (!is.null(sampler) & !is.null(reference_P)) {
    if (sampler$dims$K != nrow(reference_P)) {
      warning(glue::glue("Reference matrix has {nrow(reference_P)} rows, but data has {sampler$dims$K} rows. Setting `reference_P` to NULL."))
      reference_P <- NULL
    }
  }

  if ('character' %in% class(ref) & !is.null(reference_P)) {
    if (ref == 'assigned') {
      # only reassigns signatures if reference set has changed
      signature_assignments <- sampler$assign_signatures_ensemble(reference_P = reference_P)
      reference_P <- sampler$reference_comparison$reference

      ref <- signature_assignments$assignments$sig_ref[
        signature_assignments$assignments$sig_est == sig
      ]
      title <- glue::glue("{title}\nAssigned signature is {ref}")
    } else if (ref == 'best') {
      assignment_res <- hungarian_assignment(sampler$MAP$P[, sampler$MAP$sig_idx[sig], drop = FALSE], reference_P)
      ref <- assignment_res$sig_ref
      title <- glue::glue("{title}\nBest match in reference is {ref}")
    }
    if (!(ref %in% colnames(reference_P)) & !is.numeric(ref)) {
      stop(glue::glue("The name {ref} is not found as a column of `reference_P`."))
    }
    ref <- reference_P[, ref]
  }

  # mutation names are rownames of reference_P if available, otherwise rownames of cosmic
  plot_dat <- NULL
  cosmic <- get_cosmic()
  use_sig_theme <- FALSE
  if (!is.null(reference_P)) {
    if(setequal(rownames(reference_P), rownames(cosmic))) {
      use_sig_theme <- TRUE
    }
  }
  if (!is.null(sig)) {
    if(setequal(names(sig), rownames(cosmic))) {
      use_sig_theme <- TRUE
    }
  }
  if (use_sig_theme) {
    COSMIC_colors <- get_cosmic_colors()

    if (!is.null(reference_P)) {
      mutations <- rownames(reference_P)
    } else {
      mutations <- names(sig)
    }

    plot_dat <- data.frame(
      'mutation' = mutations
    ) %>%
      # decompose left[center]right into its parts
      # where center is reference>mutated
      dplyr::mutate(
        center = substr(mutation, 3, 5),
        left = substr(mutation, 1, 1),
        right = substr(mutation, 7, 7),
        reference = substr(mutation, 3, 3)
      ) %>%
      # sort by center first, then left, then right
      dplyr::arrange(
        center, left, right
      ) %>%
      # set factor levels to match this ordering
      dplyr::mutate(
        center = factor(center, levels = unique(center))
      ) %>%
      # create facet labels and x tick labels with COSMIC colors
      dplyr::mutate(
        x.label = paste0(
          left, "<span style = 'color: ", COSMIC_colors[center],
          ";'>", reference, "</span>", right
        ),
        facet.label = paste0(
          center, "<br><span style = 'color: ", COSMIC_colors[center],
          ";'>", paste0(rep('|', 50), collapse = '') , "</span>"
        )
      )
  } else {
    if (!is.null(reference_P)) {
      mutations <- 1:nrow(reference_P)
    } else {
      mutations <- 1:length(sig)
    }
    plot_dat <- data.frame(
      'mutation' = mutations
    ) %>%
      dplyr::mutate(
        x.label = mutations,
        facet.label = rep("", length(mutations)),
        center = rep("", length(mutations)),
        left = rep("", length(mutations)),
        right = rep("", length(mutations)),
        reference = rep("", length(mutations))
      )
  }

  plot_ref <- !is.null(ref) & !is.null(reference_P)
  plot_est <- !is.null(sampler)

  if (plot_ref & plot_est & cosine) {
    cosine_sim <- pairwise_sim(sampler$MAP$P[, sampler$MAP$sig_idx[sig], drop = FALSE], ref)[1,1]
    cosine_sim <- round(cosine_sim, 3)
    title <- glue::glue("{title}\nCosine similarity = {cosine_sim}")
  }

  if (plot_ref) {
    if (is.null(names(ref))) {
      names(ref) <- mutations
    }
  }  

  if (plot_ref) {
    ref_plot_dat <- data.frame(
      'ref' = ref,
      'mutation' = names(ref)
    )
    plot_dat <- merge(plot_dat, ref_plot_dat)
  }
  if (plot_est) {
    if (sum(sampler$MAP$A) == 1) {
      est_plot_dat <- data.frame(
        'mutation' = mutations,
        'mean' = sampler$MAP$P[, sampler$MAP$sig_idx[sig]],
        'lower' = sampler$credible_intervals$P[[1]], # already a single column
        'upper' = sampler$credible_intervals$P[[2]]
      )
    } else {
      est_plot_dat <- data.frame(
        'mutation' = mutations,
        'mean' = sampler$MAP$P[, sampler$MAP$sig_idx[sig]],
        'lower' = sampler$credible_intervals$P[[1]][, sampler$MAP$sig_idx[sig]],
        'upper' = sampler$credible_intervals$P[[2]][, sampler$MAP$sig_idx[sig]]
      )
    }
    plot_dat <- merge(plot_dat, est_plot_dat)
  }

  # build bar plot
  plot <- plot_dat %>%
    ggplot2::ggplot(ggplot2::aes(x = x.label, fill = center)) +
    ggplot2::facet_grid(cols = ggplot2::vars(facet.label), scales = 'free') +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::ggtitle(title) +
    ggthemes::theme_few()

  if (use_sig_theme) {
    plot <- plot +
      ggplot2::scale_fill_manual(values = COSMIC_colors) +
      get_sig_theme()
  } else {
    plot <- plot + 
      ggplot2::scale_fill_manual(values = c("grey")) +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_blank()
      )
  }

  if (plot_ref) {
    plot <- plot +
      ggplot2::geom_bar(ggplot2::aes(y = ref), stat = "identity", width = 0.5)
  }
  if (plot_est) {
    plot <- plot +
      ggplot2::geom_point(ggplot2::aes(y = mean)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper))
  }

  return(plot)
}


#' Get theme used for signature plots
#' @description
#' Helper function for plot_sig and plot_sig_distribution
#'
#' @return ggplot2 theme object
#' @noRd
get_sig_theme <- function() {
  ggplot2::theme(
    axis.text.x = ggtext::element_markdown(
      size = 6, angle = 90, vjust = 0.5, hjust=1
    ),
    axis.title = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(linewidth = 0.3),
    panel.spacing.x = ggplot2::unit(0, "lines"),
    panel.border = ggplot2::element_blank(),
    strip.text = ggtext::element_markdown(),
    legend.position = "none"
  )
}


#' Plot a summary of the sampler's results
#' @description
#' Plot a summary of the sampler's results, including the median contribution of each signature to the samples
#' and the cosine similarity between the estimated and assigned reference signatures. Each sampler is plotted as its own
#' x-axis tick, labeled by sampler name and sample size (G) of its dataset.
#'
#' @param sampler_list named list of bayesNMF_sampler objects
#' @param reference_P matrix or "cosmic", reference signatures to align to
#' @param title string, title of the plot
#' @param fontsize integer, font size
#' @param keep_all_ref boolean, whether to keep all reference signatures in the plot, or only the ones assigned to at least one sampler
#'
#' @return ggplot2 object
#' @export
plot_summary <- function(sampler_list, reference_P = "cosmic", title = "", fontsize = 20, keep_all_ref = FALSE) {
  if ("bayesNMF_sampler" %in% class(sampler_list)) {
    sampler_list <- list("Sampler" = sampler_list)
  }
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

  summary_df <- summarize_samplers(sampler_list, reference_P)
  color_range = 2**seq(
    floor(min(log(summary_df$Med_Contribution, base = 2))),
    ceiling(max(log(summary_df$Med_Contribution, base = 2))),
    by = 1
  )

  yaxis_limits = rev(colnames(reference_P))
  if (!keep_all_ref) {
    yaxis_limits = yaxis_limits[yaxis_limits %in% summary_df$Reference_Signature]
  }

  if (!is.null(reference_P)) {
    plot <- summary_df %>%
      ggplot2::ggplot(ggplot2::aes(
        x = Name, y = as.factor(Reference_Signature), color = Med_Contribution
      )) +
      ggplot2::geom_point(ggplot2::aes(size = Cosine_Similarity)) +
      ggplot2::scale_y_discrete(limits = yaxis_limits) +
      ggplot2::labs(size = "Cosine Similarity to\nReference Signature")
  } else {
    plot <- summary_df %>%
      ggplot2::ggplot(ggplot2::aes(
          x = Name, y = as.factor(Signature), color = Med_Contribution
      )) +
      ggplot2::geom_point(size = 6)
  }

  ylab = ifelse(is.null(reference_P), "Signature", "Reference Signature")
  plot <- plot +
    ggplot2::scale_x_discrete(
      position = "top"
    ) +
    ggplot2::scale_size_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      limits = c(0, 1)
    ) +
    ggplot2::scale_color_gradient(
      low = 'yellow', high = 'blue', trans = 'log2',
      breaks = color_range,
      limits = c(min(color_range), max(color_range))
    ) +
    ggplot2::labs(
      color = "Median Contribution\n(# Mutations)",
      title = title,
      y = ylab
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(size = fontsize),
      axis.title = ggplot2::element_text(size = fontsize),
      legend.title = ggplot2::element_text(size = fontsize),
      axis.text.x = ggplot2::element_text(size = fontsize, angle = 45, hjust = 0, vjust = 0.5),
      axis.title.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = fontsize, hjust = 0.5)
    )

  return(plot)
}

#' Plot label switching diagnostic
#' @description
#' Plot for each iteration, the assigned reference and inclusion status for each estimated signature.
#' If label switching occurs, colors will switch between estimated signatures.
#' This is a visual diagnostic for label switching, not a statistical test.
#' Recommended reference matrices include COSMIC reference signatures (default), a different literature reference, or the final MAP estimate of P if a reference does not exist.
#' 
#' @param sampler bayesNMF_sampler object
#' @param reference_P matrix or "cosmic", reference signatures to align to
#' @param idx vector of indices or "all", indices of iterations to plot
#' @param combine_below numeric, threshold for combining signatures into "other"
#'
#' @return ggplot2 object
#' @export
plot_label_switching <- function(
  sampler,
  reference_P = "cosmic",
  idx = "all",
  combine_below = 0.8 # combine signatures into "other" if cos is below threshold
) {
  save_df <- FALSE
  if ("character" %in% class(reference_P)) {
    if (reference_P == "cosmic") {
      reference_P <- get_cosmic()
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

  if ("character" %in% class(idx)) {
    if (idx == "all") {
      save_df <- TRUE
      idx <- 1:length(sampler$samples$P)
    } else {
      stop("Parameter `idx` must be a vector of indices or 'all'")
    }
  }
  if (min(idx) > sampler$state$iter) {
    stop(glue::glue("Parameter `idx` must be a vector of indices less than or equal to the total number of iterations: {sampler$state$iter}"))
  }

  other_name <- glue::glue("Other (similarity < {combine_below})")
  # only recomputes label switching df if reference set has changed
  if (
    is.null(sampler$reference_comparison$label_switching_df) |
      !check_matrix_equal(reference_P, sampler$reference_comparison$reference_P)
  ) {
    df <- lapply(idx, function(i) {
      P <- sampler$samples$P[[i]]
      A <- sampler$samples$A[[i]]
      cos <- hungarian_assignment(
        P, reference_P = reference_P,
        return_mat = TRUE, keep_all_est = TRUE
      )

      df_iter <- data.frame(
        iter = i,
        estimated = rownames(cos),
        assigned = colnames(cos),
        cosine_sim = diag(cos)
      ) %>%
        dplyr::mutate(
          k = as.numeric(stringr::str_remove(estimated, "Est"))
        ) %>%
        dplyr::arrange(k)

      # after arrange(k), latent dimensions are in same order as A
      df_iter$included <- ifelse(as.logical(A), "Included", "Excluded")
      return(df_iter)
    }) %>% 
      dplyr::bind_rows()

    if (save_df) {
      sampler$reference_comparison$label_switching_df <- df
    }
  } else {
    df <- sampler$reference_comparison$label_switching_df %>%
      dplyr::filter(iter %in% idx)
  }

  df <- df %>%
    dplyr::mutate(
      assigned = ifelse(
        cosine_sim < combine_below,
        other_name,
        assigned
      )
    )

  # order assigned with final assignments earliest in the factor
  # to give them the most unique colors
  df <- df %>%
    dplyr::arrange(-iter)
 
  levels <- unique(df$assigned)
  levels <- c(other_name, levels[levels != other_name])

  df <- df %>%
    dplyr::mutate(
      assigned = factor(assigned, levels = levels)
    )
    
  idx_annotations <- get_idx_annotations(sampler)
  p_main <- df %>%
    ggplot2::ggplot(ggplot2::aes(
      x = iter,
      y = factor(k, levels = sort(unique(k))),
      fill = assigned, alpha = cosine_sim
    )) +
    ggplot2::geom_tile() +
    ggplot2::geom_point(
      ggplot2::aes(shape = factor(included)),
      size = 2,
      stroke = 1.2,
      color = "black",
      na.rm = TRUE
    )

  # add annotations for tempering / convergence
  if (any(idx_annotations$line_df$x %in% idx)) {
    p_main <- p_main + ggplot2::geom_vline(
      data = idx_annotations$line_df %>% dplyr::filter(x %in% idx),
      ggplot2::aes(xintercept = x, color = what),
      linewidth = 1
    )
  }

  p_main <- p_main +
    ggplot2::scale_alpha(limits = c(0, 1), range = c(0.2, 1)) +
    ggplot2::scale_shape_manual(
      name = "Status",
      values = c("Excluded" = 4,  # "x" for excluded
                "Included" = NA) # nothing for included
    ) +
    ggplot2::labs(
      x = "Iteration",
      y = "Estimated Signature",
      fill = "Assigned Reference",
      alpha = "Cosine Similarity"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 20)
    )

  # add annotation for inference window
  if (
    sampler$state$iter >= sampler$specs$convergence_control$MAP_over &
    (idx_annotations$rect_df$xmin %in% idx | idx_annotations$rect_df$xmax %in% idx)
  ) {
    p_main <- p_main +
      ggplot2::geom_rect(
        data = idx_annotations$rect_df %>% dplyr::mutate(xmin = max(xmin, min(idx)), xmax = min(xmax, max(idx))),
        ggplot2::aes(
          xmin = xmin, xmax = xmax,
          ymin = -Inf, ymax = Inf,
          color = what
        ),
        inherit.aes = FALSE,
        alpha = 0.5, fill = "lightblue"
      )
  } 

  # update annotation colors and legend
  if (length(idx_annotations$used_levels) > 0) {
    p_main <- p_main +
      ggplot2::scale_color_manual(
        name = "Annotations",
        values = idx_annotations$colors,
        breaks = names(idx_annotations$colors)
      )
  }

  # update fill color
  # set to grey if assigned is "Other" and default colors otherwise
  cols <- setNames(
    c("grey", scales::hue_pal()(length(unique(df$assigned)) - 1)),
    c(other_name, setdiff(unique(df$assigned), other_name))
  )
  p_main <- p_main +
    ggplot2::scale_fill_manual(values = cols)

  # update legend order
  p_main <- p_main +
   ggplot2::guides(
    shape = ggplot2::guide_legend(order = 1),
    fill  = ggplot2::guide_legend(order = 2),
    alpha = ggplot2::guide_legend(order = 3),
    color = ggplot2::guide_legend(order = 4)
  )

  # add text annotations
  p <- suppressWarnings(add_annotations(p_main, sampler, idx_annotations, idx))

  # return plot
  return(p)
}


#' Plot estimated signature distribution
#' @description
#' Plot the distribution of mutational counts across signatures for each mutation
#' type across all, a subset, or one subject.
#'
#' @param sampler bayesNMF_sampler object
#' @param subjects integer or vector, indices of subjects to include (default is all subjects)
#' @param reference_P matrix, "cosmic" (default), or NULL, reference signatures to align to
#' @param title string, title of the produced plot
#'
#' @return ggplot2 object
#' @export
plot_signature_dist <- function(
  sampler, subjects = 1:ncol(sampler$data),
  reference_P = "cosmic",
  title = "Distribution of Signature Allocation"
) {
  if ('character' %in% class(reference_P)) {
    if (reference_P == 'cosmic') {
      reference_P = get_cosmic()
    } else {
      stop("Parameter `reference_P` must be a matrix, 'cosmic', or NULL")
    }
  }

  if (!is.null(reference_P)) {
    if (sampler$dims$K != nrow(reference_P)) {
      warning(glue::glue("Reference matrix has {nrow(reference_P)} rows, but data has {sampler$dims$K} rows. Setting `reference_P` to NULL."))
      reference_P <- NULL
    }
  }

  M = matrix(sampler$data[,subjects], ncol = length(subjects))

  all_counts <- do.call(cbind, args = lapply(1:sum(sampler$MAP$A), function(n) {
    if (length(subjects) == 1) {
      Mhat_n <- sampler$MAP$P[,sampler$MAP$sig_idx[n], drop = FALSE] * sampler$MAP$E[sampler$MAP$sig_idx[n], subjects]
    } else {
      Mhat_n <- sampler$MAP$P[,sampler$MAP$sig_idx[n], drop = FALSE] %*% sampler$MAP$E[sampler$MAP$sig_idx[n], subjects, drop = FALSE]
    }
    count_n <- rowSums(Mhat_n)
    return(count_n)
  }))
  rownames(all_counts) <- rownames(get_cosmic())

  Mhat <- sampler$get_Mhat()
  all_counts <- cbind(all_counts,  rowSums(M - Mhat[,subjects, drop = FALSE]))
  colnames(all_counts) <- c(paste0("Signature", 1:sum(sampler$MAP$A)), 'resid')

  if (!is.null(reference_P)) {
    # only reassigns signatures if reference set has changed
    signature_assignments <- sampler$assign_signatures_ensemble(reference_P = reference_P)
    reference_P <- sampler$reference_comparison$reference
    ref <- signature_assignments$assignments$sig_ref

    colnames(all_counts)[1:sum(sampler$MAP$A)] = ref
  }
    
  COSMIC_colors <- get_cosmic_colors()
  all_counts <- all_counts %>%
    data.frame() %>%
    dplyr::mutate(total = rowSums(M)) %>%
    tibble::rownames_to_column("mutation") %>%
    dplyr::mutate(
      center = substr(mutation, 3, 5),
      left = substr(mutation, 1, 1),
      right = substr(mutation, 7, 7),
      reference = substr(mutation, 3, 3),mutation,
      x.label = paste0(
        left, "<span style = 'color: ", COSMIC_colors[center],
        ";'>", reference, "</span>", right
      ),
      facet.label = paste0(
        center, "<br><span style = 'color: ", COSMIC_colors[center],
        ";'>", paste0(rep('|', 50), collapse = '') , "</span>"
      )
    ) 
    
  plot_dat <- all_counts %>%
    tidyr::pivot_longer(
      2:(sum(sampler$MAP$A) + 1),
      names_to = 'signature',
      values_to = 'count'
    ) %>%
    dplyr::mutate(
      signature = dplyr::case_when(
        signature == 'resid' & count < 0 ~ 'neg_resid',
        signature == 'resid' & count >=0 ~ 'pos_resid',
        TRUE ~ signature
      )
    ) 

  plot <- plot_dat %>%
    ggplot2::ggplot(ggplot2::aes(x = x.label, fill = signature, y = count)) +
    ggplot2::facet_grid(cols = ggplot2::vars(facet.label), scales = 'free') +
    ggplot2::geom_bar(
      ggplot2::aes(color = NULL),
      position = "stack", stat = 'identity'
    ) +
    get_sig_theme() +
    ggplot2::geom_point(
      data = all_counts %>% dplyr::mutate(observed = ' '),
      ggplot2::aes(y = total, color = observed), fill = "black"
    ) +
    ggplot2::scale_color_manual(values = c(' ' = 'black')) +
    ggplot2::theme(
      legend.position = 'top',
      axis.title.y = ggplot2::element_text()
    ) +
    ggplot2::labs(
      fill = ifelse(is.null(reference_P), "Assigned Signature", "Signature"),
      color = "Observed Count",
      y = "Count",
      title = title
    )

  return(plot)
}