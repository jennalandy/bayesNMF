#' A single sample of the P matrix votes for signature assignment based on cosine similarity
#'
#' @param P_post matrix, posterior sample of P
#' @param res bayesNMF object
#' @param reference matrix, "cosmic", or NULL, reference signatures to align to
#'
#' @return data frame
#' @noRd
vote <- function(P_post, res, reference) {
    P_post_discovery = P_post[,(res$MAP$A[1,] == 1 & res$final_Theta$recovery == FALSE)]

    if (sum(res$MAP$A[1,] == 1 & res$final_Theta$recovery == FALSE) == 0) {
        return(NULL)
    } else if (sum(res$MAP$A[1,] == 1 & res$final_Theta$recovery == FALSE) == 1) {
        P_post_discovery <- matrix(P_post_discovery, ncol = 1)
    }
    sim <- pairwise_sim(P_post_discovery, reference)
    assigned <- hungarian_algorithm(sim)

    votes <- data.frame(
        cosine_sim = diag(assigned), # voting weight
        sig = colnames(assigned),    # voting for
        # record actual indices in case of A switching
        n = which(res$final_Theta$recovery[res$MAP$A[1,] == 1] == FALSE)
    )
    return(votes)
}

#' Hungarian algorithm for signature alignment
#'
#' @param matrix matrix, difference or similarity matrix
#' @param which boolean, whether we want to maximize or minimize diagonals in the matrix
#'
#' @return aligned matrix
#' @noRd
hungarian_algorithm <- function(matrix, which = "max") {
    if (which == "max") {
        matrix <- -1 * matrix
    }
    hungarian_res <- RcppHungarian::HungarianSolver(matrix)
    hungarian_alignment <- data.frame(hungarian_res$pairs) %>%
        dplyr::filter(X1 != 0 & X2 != 0)

    aligned_matrix <- matrix[
        hungarian_alignment$X1,
        hungarian_alignment$X2
    ]

    if (which == "max") {
        aligned_matrix <- -1 * aligned_matrix
    }

    if (nrow(matrix) == 1) {
        aligned_matrix = matrix(aligned_matrix, dimnames = list(NULL, names(aligned_matrix)))
    }
    return(aligned_matrix)
}

#' Posterior inference on signature assignment
#'
#' This function performs the Hungarian algorithm on each P matrix among the
#' posterior samples. Each sample gets to vote for it's assignment with  cosine
#' similarity as voting weights.
#'
#' @param res bayesNMF object
#' @param reference matrix, "cosmic", or NULL, reference signatures to align to
#'
#' @return list, signature assignments, MAP and credible intervals for cosine similarity
#' @export
#'
#' @examples
#' res <- readRDS("examples/plot_example.rds")
#' assign <- signature_asssignment_inference(res)
#' assign$assignment
#' assign$MAP$cos_sim
#' assign$credible_intervals$cos_sim
signature_asssignment_inference <- function(res, reference = 'cosmic') {
    if ('character' %in% class(reference)) {
        if (reference == 'cosmic') {
            reference = get_cosmic()
        } else {
            stop("Parameter `reference` must be a matrix or 'cosmic'")
        }
    }

    # get votes from all posterior samples (discovery signatures only)
    assignment_votes <- lapply(res$posterior_samples$P, function(P_post) {vote(P_post, res, reference)})
    if (is.null(assignment_votes[[1]])) {
        # happens if there are no discovery signatures to assign
        assignments <- NULL
    } else {
        # sum the votes for each assigned signature to each index
        summarized_votes <- do.call(rbind, assignment_votes) %>%
            dplyr::group_by(sig, n) %>%
            dplyr::summarize(votes = sum(cosine_sim), .groups = 'keep') %>%
            dplyr::ungroup()

        # final assignment based on maximum votes per index
        assignments <- summarized_votes %>%
            dplyr::group_by(n) %>%
            dplyr::summarize(max_votes = max(votes), total = sum(votes), .groups = 'keep') %>%
            merge(summarized_votes) %>%
            dplyr::mutate(score = votes/total) %>%
            dplyr::filter(max_votes == votes) %>%
            dplyr::arrange(n) %>%
            dplyr::select(sig, score, n)
    }

    # append recovery signatures
    if (sum(res$final_Theta$recovery) > 0) {                  # there are recovery
        if (sum(res$MAP$A[1,1:ncol(reference)] == 1) > 0) { # at least one is included
            if (is.null(assignments)) {                       # there are no discovery included
                assignments <- data.frame(                    # assignments is just recovery signatures
                    sig = colnames(reference)[res$MAP$A[1,1:ncol(reference)] == 1],
                    score = 1,
                    n = 1:sum(res$MAP$A[1,1:ncol(reference)] == 1)
                )
            } else {                                          # there is at least one discovery
                assignments <- rbind(                         # append recovery to discovery assignments
                    data.frame(
                        sig = colnames(reference)[res$MAP$A[1,1:ncol(reference)] == 1],
                        score = 1,
                        n = 1:sum(res$MAP$A[1,1:ncol(reference)] == 1)
                    ),
                    assignments
                )
            }
        }
    }

    # compute cosine similarity between each posterior sample's estimated
    # and it's assigned signatures
    assigned_cos_sims <- lapply(res$posterior_samples$P, function(P_post) {
        P_post = P_post[,res$MAP$A[1,] == 1]
        sim <- pairwise_sim(P_post, reference)
        sapply(1:nrow(assignments), function(i) {
            sim[assignments$n[i], assignments$sig[i]]
        })
    }) %>%
        do.call(rbind, .)

    # return assignments, MAP and credible intervals for cosine similarity
    assignment_res <- list(
        'assignment' = assignments,
        'MAP'  = list(
            'cos_sim' = assigned_cos_sims %>% colMeans()
        ),
        'credible_intervals' = list(
            'cos_sim' = list(
                apply(assigned_cos_sims, 2, function(col) {quantile(col, 0.025)}),
                apply(assigned_cos_sims, 2, function(col) {quantile(col, 0.975)})
            )
        )
    )
    return(assignment_res)
}

#' For a named list of results, gets signature assignments, cosine similarities,
#' and distribution of mutations attributed to signatures in each.
#'
#' @param res_list Named list, containing one or more bayesNMF objects. Names will become identifiers along the top of the plot.
#' @param reference matrix, "cosmic", or NULL, reference signatures to align to
#'
#' @return data frame
#' @noRd
get_results_df <- function(res_list, reference) {
    total_counts <- NULL
    total_props <- NULL
    first <- TRUE
    Gs <- c()
    for (name in names(res_list)) {
        res <- res_list[[name]]
        if (is.null(res$converged_at)) {
            print(paste('not done:', name))
            next()
        }

        # record sample size
        Gs <- c(Gs, res$model$dims$G)

        # assign estimated factors to reference, update names
        if (!is.null(reference)) {
            assignment_res <- signature_asssignment_inference(res, reference)
            if (sum(res$MAP$A) == 1) {
                res$MAP$P = matrix(res$MAP$P, ncol = 1)
                res$MAP$E = matrix(res$MAP$E, nrow = 1)
            }
            colnames(res$MAP$P) <- assignment_res$assignment$sig
            rownames(res$MAP$E) <- assignment_res$assignment$sig
        } else {
            assignment_res <- NULL
        }

        # re-scale so factors sum to 1
        res$MAP$E <- sweep(res$MAP$E, 1, colSums(res$MAP$P), '*')
        res$MAP$P <- sweep(res$MAP$P, 2, colSums(res$MAP$P), '/')

        # record number of mutations attributed to each signature
        # median among subjects with the signature
        # vector is length of # reference, we only fill in the ones seen in this sample
        this_counts <- apply(res$MAP$E, 1, function(col) {
            median(col[col > 0])
        })
        if (!is.null(reference)) {
            sig_counts <- rep(0, ncol(reference))
            names(sig_counts) <- colnames(reference)
            sig_counts[names(this_counts)] <- this_counts
        } else {
            names(this_counts) <- paste0(name, ' n', 1:length(this_counts))
            sig_counts <- this_counts
        }


        # record average cosine similarity between estimated and assigned reference
        # vector is length of # reference, we only fill in the ones seen in this sample
        if (!is.null(reference)) {
            this_cos <- assignment_res$MAP$cos_sim
            cosine_sim <- rep(0, ncol(reference))
            names(cosine_sim) <- colnames(reference)
            cosine_sim[assignment_res$assignment$sig] <- this_cos
        } else {
            cosine_sim = NULL
        }

        # add data to results matrix
        if (first) {
            total_counts <- matrix(
                c(name, sig_counts), nrow = 1,
                dimnames = list(c(1), c('Name', names(sig_counts)))
            )
            if (!is.null(reference)) {
                total_cos <- matrix(
                    c(name, cosine_sim), nrow = 1,
                    dimnames = list(c(1), c('Name', names(cosine_sim)))
                )
            } else {
                total_cos = NULL
            }
            first = FALSE
        } else {
            if (is.null(reference)) {
                # no common signature naming scheme
                # add new columns to existing matrix with all 0s
                total_counts <- cbind(
                    total_counts,
                    matrix(
                        0,
                        nrow = nrow(total_counts),
                        ncol = length(sig_counts),
                        dimnames = list(
                            rownames(total_counts),
                            names(sig_counts)
                        )
                    )
                )

                # add old columns to new row with all 0s
                sig_counts <- c(
                    c('Name' = name, rep(0, ncol(total_counts) - length(sig_counts)-1)),
                    sig_counts
                )
                names(sig_counts) = colnames(total_counts)
            } else {
                sig_counts = c('Name' = name, sig_counts)
                cosine_sim = c('Name' = name, cosine_sim)
            }

            # add new data to matrix
            total_counts <- rbind(
                total_counts, sig_counts
            )
            total_cos <- rbind(
                total_cos, cosine_sim
            )
        }
    }

    # reformat total_counts into a long data frame
    # will have Name, Signature, Med_Contribution
    total_counts <- data.frame(total_counts) %>%
        dplyr::mutate(G = Gs) %>%
        tidyr::pivot_longer(
            2:ncol(total_counts),
            names_to = "Signature",
            values_to = "Med_Contribution"
        )

    # reformat total_cos into a long data frame
    # will have Name, Signature, Cosine_Similarity
    if (!is.null(reference)) {
        total_cos <- data.frame(total_cos) %>%
            tidyr::pivot_longer(
                2:ncol(total_cos),
                names_to = "Signature",
                values_to = "Cosine_Similarity"
            )
    }

    # combine results, ensure correct column types
    if (!is.null(reference)) {
        results <- merge(total_counts, total_cos) %>%
            dplyr::mutate(
                Cosine_Similarity = as.numeric(Cosine_Similarity),
                Med_Contribution = as.numeric(Med_Contribution)
            )
        results$Signature = factor(results$Signature, levels = rev(colnames(reference)))
    } else {
        results <- total_counts %>%
            dplyr::mutate(
                Med_Contribution = as.numeric(Med_Contribution)
            )
    }

    return(results)
}
