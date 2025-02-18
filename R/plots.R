#' Get COSMIC reference signatures matrix
#'
#' @return matrix
#' @export
#'
#' @examples
#' cosmic <- get_cosmic()
get_cosmic <- function() {
    P <- read.csv(
        "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3.1_SBS_GRCh37.txt", # nolint: line_length_linter.
        sep = "\t"
    )
    rownames(P) <- P$Type
    P <- P[, 2:ncol(P)]
    P <- as.matrix(P)
    return(P)
}

#' Get a named list of colors used for COSMIC signature plots
#'
#' @return list
#' @export
#'
#' @examples
#' COSMIC_colors <- get_cosmic_colors()
get_cosmic_colors <- function() {
    COSMIC_colors <- c(
        'C>A' = rgb(8, 181, 236, maxColorValue = 255),
        'C>G' = rgb(0, 0, 0, maxColorValue = 255),
        'C>T' = rgb(225, 37, 33, maxColorValue = 255),
        'T>A' = rgb(198, 193, 195, maxColorValue = 255),
        'T>C' = rgb(153, 200, 87, maxColorValue = 255),
        'T>G' = rgb(233, 190, 189, maxColorValue = 255)
    )
    return(COSMIC_colors)
}

#' Get theme used for signatures plots
#'
#' @return ggplot2 theme object
#' @noRd
get_sig_theme <- function() {
    ggthemes::theme_few() +
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

#' Plot a single mutational signature.
#'
#' This function can be used to plot a reference mutational signature, an
#' estimated mutational signature with posterior uncertainty as points with
#' error bars, or both. The default reference is COSMIC, and unless otherwise
#' specified, the reference signature chosen to plot is that with highest
#' cosine similarity to the estimated signature.
#'
#' One of the following must be true
#' - res and sig are specified
#' - ref_matrix is specified and ref is integer, a column name, or "best"
#' - ref is a vector
#'
#' @param res bayesNMF object, results list
#' @param sig integer, index of the estimated signature to plot
#' @param ref integer, column name string, numeric vector, or "best"
#' @param ref_matrix matrix or "cosmic", reference signatures for comparison
#' @param title string
#' @param cosine boolean, whether to report cosine similarity between estimated
#' and reference. Ignored if either of sig, ref are NULL.
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' res <- readRDS("examples/plot_example.rds")
#' data <- readRDS("examples/3_64_1_cosmic.rds")
#' plot_sig(res = res, sig = 1, title = "Plotting an estimated signature with the best matched COSMIC signature")
#' plot_sig(res = res, sig = 1, ref_matrix = data$P, title = "Plotting an estimated signature with the best matched signature from a custom reference")
#' plot_sig(ref = "SBS3", title = "Plotting a reference signature alone")
#' plot_sig(ref = runif(96), title = "Plotting a custom reference signature alone")
#' plot_sig(res = res, sig = 1, ref = "SBS3", title = "Plot an estimated signature with a specific COSMIC signature")
# plot_sig(res = res, sig = 1, ref = "SBS88", ref_matrix = data$P, title = "Plot an estimated signature with a specific signature from a custom reference")
plot_sig <- function(
    res = NULL, sig = NULL,
    ref = "best", ref_matrix = "cosmic",
    title = "", cosine = TRUE
) {
    if (is.null(res) & is.null(ref)) {
        stop("Must provide either `res` (results of bayesNMF) or `ref` (vector for reference signature or column name of `ref_matrix`).")
    } else if (is.null(ref) & is.null(res) & is.null(sig)) {
        stop("To plot a signature from `res`, must specify signature index with `sig`.")
    }

    # rescale factors to sum to 1
    if (!is.null(res)) {
        res$credible_intervals$P[[1]] <- sweep(res$credible_intervals$P[[1]], 2, colSums(res$MAP$P), '/')
        res$credible_intervals$P[[2]] <- sweep(res$credible_intervals$P[[2]], 2, colSums(res$MAP$P), '/')
        res$MAP$P <- sweep(res$MAP$P, 2, colSums(res$MAP$P), '/')
    }

    if ('character' %in% class(ref)) {
        if ('character' %in% class(ref_matrix)) {
            if (ref_matrix == 'cosmic') {
                ref_matrix = get_cosmic()
            } else {
                stop("Parameter `ref_matrix` must be a matrix or the string 'cosmic'.")
            }
        }
        if (ref == 'best') {
            sim <- pairwise_sim(res$MAP$P[,sig], ref_matrix)
            assignment_res <- assign_signatures(sim)
            ref = colnames(assignment_res)
            title = glue::glue("{title}\nBest match in reference is {ref}")
        }
        if (!(ref %in% colnames(ref_matrix)) & !is.numeric(ref)) {
            stop(glue::glue("The name {ref} is not found as a column of `ref_matrix`."))
        }
        ref = ref_matrix[,ref]
    }
    cosmic = get_cosmic()

    plot_ref = !is.null(ref)
    plot_est = !is.null(res)

    if (plot_ref & plot_est & cosine) {
        cosine_sim = pairwise_sim(res$MAP$P[,sig], ref)[1,1]
        cosine_sim = round(cosine_sim, 3)
        title = glue::glue("{title}\nCosine similarity = {cosine_sim}")
    }

    if (plot_ref) {
        if (is.null(names(ref))) {
            names(ref) = rownames(cosmic)
        }
    }

    COSMIC_colors <- get_cosmic_colors()
    plot_dat <- data.frame(
        'mutation' = rownames(cosmic)
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

    if (plot_ref) {
        ref_plot_dat <- data.frame(
            'ref' = ref,
            'mutation' = names(ref)
        )
        plot_dat <- merge(plot_dat, ref_plot_dat)
    }
    if (plot_est) {
        est_plot_dat <- data.frame(
            'mutation' = rownames(cosmic),
            'mean' = res$MAP$P[,sig],
            'lower' = res$credible_intervals$P[[1]][,sig],
            'upper' = res$credible_intervals$P[[2]][,sig]
        )
        plot_dat <- merge(plot_dat, est_plot_dat)
    }

    # build bar plot
    plot <- plot_dat %>%
        ggplot2::ggplot(ggplot2::aes(x = x.label, fill = center)) +
        ggplot2::facet_grid(cols = ggplot2::vars(facet.label), scales = 'free') +
        ggplot2::scale_fill_manual(values = COSMIC_colors) +
        get_sig_theme() +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
        ggplot2::ggtitle(title)

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

#' Plot signature contributions
#'
#' @param res_list Named list, containing one or more bayesNMF objects. Names will become identifiers along the top of the plot.
#' @param ref_matrix matrix, "cosmic", or NULL, reference signatures to align to
#' @param title string, title of the produced plot
#' @param return_df boolean, whether to return summary data frame
#'
#' @return ggplot object, or list holding ggplot object and data frame if return_df = TRUE
#' @export
#'
#' @examples
#' res <- readRDS("examples/plot_example.rds")
#' data <- readRDS("examples/3_64_1_cosmic.rds")
#' plot_results(list("Example" = res), title = "Results of a single run")
#' plot_results(list("Example" = res), reference = data$P, title = "Results of a single run with custom reference")
plot_results <- function(res_list, ref_matrix = 'cosmic', title = "", return_df = TRUE) {
    if ('character' %in% class(ref_matrix)) {
        if (ref_matrix == 'cosmic') {
            ref_matrix = get_cosmic()
        } else {
            stop("Parameter `ref_matrix` must be a matrix or 'cosmic'")
        }
    }

    results_df <- get_results_df(res_list, ref_matrix) %>%
        dplyr::filter(Med_Contribution > 0)

    color_range = 2**seq(
        floor(min(log(results_df$Med_Contribution, base = 2))),
        ceiling(max(log(results_df$Med_Contribution, base = 2))),
        by = 1
    )

    results_df <- results_df %>%
        dplyr::mutate(
            Name = paste0(Name, " (", G, ")")
        )

    plot <- results_df %>%
        ggplot2::ggplot(ggplot2::aes(
            x = Name, y = Signature, color = Med_Contribution
        ))

    if (!is.null(ref_matrix)) {
        plot <- plot +
            ggplot2::geom_point(ggplot2::aes(size = Cosine_Similarity)) +
            ggplot2::scale_y_discrete(limits = rev(colnames(ref_matrix))) +
            ggplot2::labs(size = "Posterior Average\nCosine Similarity to\nReference Signature")
    } else {
        plot <- plot +
            ggplot2::geom_point(size = 6)
    }

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
            color = "Median Contribution\n(# Mutations) per Mb",
            title = title
        ) +
        ggplot2::theme(
            text = ggplot2::element_text(size = 11),
            axis.title = ggplot2::element_text(size = 11),
            legend.title = ggplot2::element_text(size = 11),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0.5),
            axis.title.x = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5)
        )

    if (return_df) {
        return(list(
            'plot' = plot,
            'df' = results_df
        ))
    } else {
        return(plot)
    }
}

#' Plot the distribution of mutational counts across signatures for each mutation
#' type across all, a subset, or one subject.
#'
#' @param res bayesNMF object, results list
#' @param subjects integer or vector, indices of subjects to include
#' @param reference matrix, "cosmic", or NULL, reference signatures to align to
#' @param title string, title of the produced plot
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' plot_signature_dist(res, title = "Plot signature distribution across all subjects")
#' plot_signature_dist(res, subject = 1, title = "Plot signature distribution of subject 1")
#' plot_signature_dist(res, subject = 1, reference = NULL, title = "Plot signature distribution of subject 1 w/o reference")
#' plot_signature_dist(res, subject = c(1,4,10), title = "Plot signature distribution of subjects 1, 4, and 10")
plot_signature_dist <- function(
    res, subjects = 1:ncol(res$model$M),
    reference = "cosmic",
    title = "Distribution of Signature Allocation"
) {
    if ('character' %in% class(reference)) {
        if (reference == 'cosmic') {
            reference = get_cosmic()
        } else {
            stop("Parameter `reference` must be a matrix, 'cosmic', or NULL")
        }
    }

    M = matrix(res$model$M[,subjects], ncol = length(subjects))

    all_counts <- do.call(cbind, args = lapply(1:sum(res$MAP$A), function(n) {
        if (length(subjects) == 1) {
            Mhat_n <- matrix(res$MAP$P[,n], ncol = 1) * res$MAP$E[n,subjects]
        } else {
            Mhat_n <- matrix(res$MAP$P[,n], ncol = 1) %*% matrix(res$MAP$E[n,subjects], nrow = 1)
        }
        count_n <- rowSums(Mhat_n)
        return(count_n)
    }))
    rownames(all_counts) <- rownames(get_cosmic())
    all_counts <- cbind(all_counts, rowSums(M - res$MAP$P %*% res$MAP$E[,subjects]))
    colnames(all_counts) <- c(paste0("Signature", 1:sum(res$MAP$A)), 'resid')

    if (!is.null(reference)) {
        sim <- pairwise_sim(res$MAP$P, reference)
        assignment_res <- assign_signatures(sim)
        ref = colnames(assignment_res)

        colnames(all_counts)[1:sum(res$MAP$A)] = ref
    }

    plot_dat <- all_counts %>%
        data.frame() %>%
        tibble::rownames_to_column("mutation") %>%
        dplyr::mutate(
            center = substr(mutation, 3, 5),
            left = substr(mutation, 1, 1),
            right = substr(mutation, 7, 7),
            reference = substr(mutation, 3, 3)
        ) %>%
        tidyr::pivot_longer(2:(sum(res$MAP$A) + 1), names_to = 'signature', values_to = 'count') %>%
        dplyr::mutate(
            signature = dplyr::case_when(
                signature == 'resid' & count < 0 ~ 'neg_resid',
                signature == 'resid' & count >=0 ~ 'pos_resid',
                TRUE ~ signature
            )
        )

    plot <- plot_dat %>%
        ggplot2::ggplot(ggplot2::aes(x = right, y = count, fill = signature)) +
        ggplot2::facet_grid(
            rows = ggplot2::vars(left),
            cols = ggplot2::vars(center),
            switch = 'y'
        ) +
        ggplot2::geom_bar(stat = 'identity') +
        ggplot2::theme(
            strip.placement = 'outside',
            plot.subtitle = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::labs(
            y = 'left',
            subtitle = 'center',
            fill = "Signature",
            title = title
        )

    return(plot)
}
