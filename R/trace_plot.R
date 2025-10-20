#' Plot trace of metrics
#' @description
#' Plot trace of metrics for the sampler, including log posterior, log likelihood, BIC, RMSE, KL-Divergence.
#' If rank is being learned, also reports rank, number of functional parameters (used in BIC), and temperature.
#' If Metropolis-Hastings steps are used, also reports average acceptance rate for elements of $P$ and $E$ (initally all 1 to find high posterior region, then true MH samples for inference).
#' Note that for log likelihood and log posterior, values from the Poisson models are not comparable to those from Normal models.
#'
#' @param sampler bayesNMF_sampler object
#' @param MAP_means boolean, whether to plot means of MAP metrics
#' @param idx vector of indices or "all" (default), indices of iterations to plot
#' @param save boolean, whether to save the plot to the sampler's output directory (default FALSE)
#'
#' @return ggplot2 object
#' @export
trace_plot <- function(sampler, MAP_means = FALSE, idx = "all", save = FALSE) {
  if ("character" %in% class(idx)) {
    if (idx == "all") {
      idx <- 1:sampler$state$iter
    } else {
      stop("Parameter `idx` must be a vector of indices or 'all'")
    }
  }
  if (min(idx) > sampler$state$iter) {
    stop(glue::glue("Parameter `idx` must be a vector of indices less than or equal to the total number of iterations: {sampler$state$iter}"))
  }

  idx_annotations <- get_idx_annotations(sampler)

  if (MAP_means) {
    dat <- sampler$state$MAP_metrics
    metric_levels <- c(
      "logposterior", "loglikelihood", "BIC",
      "RMSE", "KL"
    )
    if (sampler$specs$learning_rank) {
      metric_levels <- c(metric_levels, "rank", "MAP_A_counts", "n_params", "mean_temp")
    }
  } else {
    dat <- sampler$state$sample_metrics
    metric_levels <- c(
      "logposterior", "loglikelihood", "BIC",
      "RMSE", "KL"
    )
    if (sampler$specs$learning_rank) {
      metric_levels <- c(metric_levels, "rank", "n_params", "temp")
    }
  }
  if (sampler$specs$MH) {
    metric_levels <- c(metric_levels, "P_mean_acceptance_rate", "E_mean_acceptance_rate")
  }

  metric_names <- list(
    "logposterior" = "Log\nPosterior",
    "loglikelihood" = "Log\nLikelihood",
    "BIC" = "BIC",
    "RMSE" = "RMSE",
    "KL" = "KL-\nDivergence",
    "rank" = "Rank",
    "MAP_A_counts" = "Rank\nCounts",
    "n_params" = "# Parameters",
    "P_mean_acceptance_rate" = "P Mean\nAcceptance",
    "E_mean_acceptance_rate" = "E Mean\nAcceptance",
    "mean_temp" = "Mean\nTemperature",
    "temp" = "Temperature"
  )

  dat <- dat %>%
    dplyr::filter(iter %in% idx)

  # plot
  p_main = dat %>%
    tidyr::pivot_longer(2:ncol(.), names_to = "metric", values_to = "value") %>%
    dplyr::filter(metric %in% metric_levels, !is.na(value), !is.infinite(value)) %>%
    dplyr::mutate(
      metric = factor(metric, levels = metric_levels, labels = metric_names[metric_levels])
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = iter, y = value)) +
    ggplot2::facet_grid(rows = ggplot2::vars(metric), scales = 'free', switch = 'y')

  # only add rectangle if iter > MAP_over AND inference is in idx
  if (
    sampler$state$iter >= sampler$specs$convergence_control$MAP_over &
    (idx_annotations$rect_df$xmin %in% idx | idx_annotations$rect_df$xmax %in% idx)
  ) {
    p_main <- p_main +
      ggplot2::geom_rect(
        data = idx_annotations$rect_df %>% dplyr::mutate(xmin = max(xmin, min(idx)), xmax = min(xmax, max(idx))),
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, color = what),
        inherit.aes = FALSE,
        alpha = 0.5, fill = "lightblue"
      )
  }

  if (any(idx_annotations$line_df$x %in% idx)) {
    p_main <- p_main +
      ggplot2::geom_vline(
        data = idx_annotations$line_df %>% dplyr::filter(x %in% idx),
        ggplot2::aes(xintercept = x, color = what),
        linewidth = 1
      ) 
  }
  p_main <- p_main +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::labs(x = "Iteration") +
    ggplot2::theme(
      text = ggplot2::element_text(size = 15),
      strip.text.y.left = ggplot2::element_text(angle = 0),
      strip.background = ggplot2::element_rect(fill = 'lightgrey'),
      strip.placement = 'outside',
      legend.position = 'bottom',
      legend.justification = 'left',
      axis.title.y = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(fill = NA, colour = NA)
    ) 

  if (length(idx_annotations$used_levels) > 0) {
    p_main <- p_main +
      ggplot2::scale_color_manual(
        name = "Annotations",
        values = idx_annotations$colors,
        breaks = names(idx_annotations$colors)
      )
  }

  if (sampler$state$converged) {
    p_main <- suppressWarnings(add_annotations(p_main, sampler, idx_annotations, idx))
  }
  if (save) {
    height <- ifelse(sampler$specs$learning_rank, 23, 18)
    name <- ifelse(MAP_means, "trace_plot_MAP.png", "trace_plot.png")
    suppressMessages(
      ggplot2::ggsave(
        file.path(sampler$specs$output_dir, name),
        plot = p_main,
        height = height, width = 12
      )
    )
  }
  return(p_main)
}

#' Add text annotations to the trace plot
#' @description
#' Helper function for trace_plot and plot_label_switching
#'
#' @param p_main ggplot2 object, the main trace plot
#' @param sampler bayesNMF_sampler object
#' @param idx_annotations list, the annotations to add
#' @param idx vector of indices or "all" (default), indices of iterations to plot
#' @param size integer, the size of the text labels
#' 
#' @return ggplot2 object
#' @noRd
add_annotations <- function(p_main, sampler, idx_annotations, idx, size = 10) {
  main_range <- ggplot2::ggplot_build(p_main)$layout$panel_params[[1]]$x.range

  if (sampler$specs$learning_rank) {
    seg <- data.frame(
      xmin = c(
        1,
        idx_annotations$line_df$x[idx_annotations$line_df$what == "Convergence"],
        idx_annotations$rect_df$xmin
      ),
      xmax = c(
        idx_annotations$line_df$x[idx_annotations$line_df$what == "Tempering done"],
        sampler$state$iter,
        idx_annotations$rect_df$xmax
      ),
      y = 1,
      eta = c(0.01, 0.01, 0.03),
      lab = c("Tempering", "MH Samples", "Inference"),
      color = c("Tempering done", "Convergence", "Inference label"),
      placement = c(0.5, 0.25, 0.5)
    )
  } else {
    seg <- data.frame(
      xmin = c(
        idx_annotations$line_df$x[idx_annotations$line_df$what == "Convergence"],
        idx_annotations$rect_df$xmin
      ),
      xmax = c(
        sampler$state$iter,
        idx_annotations$rect_df$xmax
      ),
      y = 1,
      eta = c(0.01, 0.03),
      lab = c("MH Samples", "Inference"),
      color = c("Convergence", "Inference label"),
      placement = c(0.25, 0.5)
    )
  }

  # remove MH Samples if !MH
  if (!sampler$specs$MH) {
    seg <- seg %>%
      dplyr::filter(lab != "MH Samples")
  }

  # only plot as endpoints if they are within idx
  seg <- seg %>%
    dplyr::mutate(
      end_min = xmin >= min(idx),
      end_max = xmax <= max(idx),
      xmin = pmax(min(idx), xmin),
      xmax = pmin(max(idx), xmax)
    )
  
  idx_annotations$colors['Inference label'] = 'turquoise3'

  p_anno <- ggplot2::ggplot(seg, ggplot2::aes(color = color)) +
    # horizontal
    ggplot2::geom_segment(ggplot2::aes(x = xmin, xend = xmax, y = 1 + eta, yend = 1 + eta)) +
    # tips
    ggplot2::geom_segment(ggplot2::aes(x = xmin, xend = xmin, y = ifelse(end_min, 1, 1 + eta), yend = 1 + eta)) +
    ggplot2::geom_segment(ggplot2::aes(x = xmax, xend = xmax, y = ifelse(end_max, 1, 1 + eta), yend = 1 + eta)) +
    # labels
    ggplot2::geom_text(ggplot2::aes(
      x = xmin + (xmax - xmin) * placement,
      y = 1.01 + eta,
      label = lab
    ), size = size, vjust = 0) +
    # match x-axis range of main plot
    ggplot2::coord_cartesian(xlim = c(main_range[1], main_range[2])) +
    ggplot2::scale_y_continuous(limits = c(1, 1.2), expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 0, l = 5),
      legend.position = 'none',
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
    ) +
    ggplot2::scale_color_manual(
      values = idx_annotations$colors[seg$color],
      breaks = names(idx_annotations$colors[seg$color])
    )

  # align limits and margins
  p_main2 <- p_main +
    ggplot2::scale_x_continuous(limits = main_range, expand = c(0,0)) +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 0, r = 5, b = 5, l = 5))

  p_anno2 <- p_anno +
    ggplot2::scale_x_continuous(limits = main_range, expand = c(0,0)) +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 5, r = 5, b = 0, l = 5))

  aligned <- cowplot::align_plots(p_anno2, p_main2, align = "v", axis = "lr")
  
  # combine plots
  p <- cowplot::plot_grid(
    aligned[[1]], aligned[[2]], ncol = 1,
    rel_heights = c(0.10, 1)
  )
  return(p)
}

#' Get annotations for the sampler's iterations
#' @description
#' Helper function for trace_plot and plot_label_switching
#'
#' @param sampler bayesNMF_sampler object
#'
#' @return list with line_df, rect_df, colors, and used_levels
#' @noRd
get_idx_annotations <- function(sampler, idx) {
  done_temp <- which(sampler$temperature_schedule == 1)[1]
  line_df <- data.frame(
    x = c(
      ifelse(sampler$state$converged, sampler$state$converged_iter, NA_integer_),
      ifelse(sampler$specs$learning_rank & sampler$state$iter >= done_temp, done_temp, NA_integer_)
    ),
    what = c("Convergence", "Tempering done")
  ) %>% 
    tidyr::drop_na()

  rect_df <- data.frame(
    xmin = sampler$state$iter - sampler$specs$convergence_control$MAP_over,
    xmax = sampler$state$iter,
    what = "Inference"
  )
  used_levels <- unique(c(as.character(line_df$what), if (exists("rect_df")) as.character(rect_df$what)))
  used_levels <- used_levels[!is.na(used_levels)]
  colors <- c(
    "Tempering done" = "orange",
    "Convergence" = "blue",
    "Inference" = "lightblue"
  )
  colors <- colors[used_levels]
  return(list(line_df = line_df, rect_df = rect_df, colors = colors, used_levels = used_levels))
}