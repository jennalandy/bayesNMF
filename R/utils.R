#' get RMSE
#'
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#'
#' @return scalar
#' @export
get_RMSE <- function(M, Theta) {
    sqrt(mean((M - Theta$P %*% Theta$E)**2))
}

#' Get KL Divergence
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#'
#' @return scalar
#' @export
get_KLDiv <- function(M, Theta) {
    Mhat <- Theta$P %*% Theta$E
    if (sum(Mhat < 0) > 1) {
        print(Theta$P)
        print(Theta$E)
    }
    Mhat[Mhat == 0] <- 1
    M[M == 0] <- 1
    sum(M * log(M / Mhat) - M + Mhat)
}

#' Reassign signatures in a similarity matrix
#'
#' @param sim_mat similarity matrix between estimated and true signatures
#'
#' @return matrix
reassign_signatures <- function(sim_mat) {
    reassignment <- RcppHungarian::HungarianSolver(-1 * sim_mat)
    reassigned_sim_mat <- sim_mat[, reassignment$pairs[,2]]
    reassigned_sim_mat
}

#' Get similarity matrix
#'
#' @param est_P estimated P (signatures matrix)
#' @param true_P true P (signatures matrix)
#'
#' @return matrix
#' @export
get_sim_mat <- function(est_P, true_P) {
    text2vec::sim2(t(est_P), t(true_P), method = 'cosine')
}

#' Get heatmap
#'
#' @param est_P estimated P (signatures matrix)
#' @param true_P true P (signatures matrix)
#'
#' @return ggplot object
#' @export
get_heatmap <- function(est_P, true_P) {
    sim_mat <- get_sim_mat(est_P, true_P)

    sim_mat_melted <- reshape2::melt(sim_mat) %>%
        dplyr::arrange(-value)

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
        ggplot2::ggplot(ggplot2::aes(x = Var1, y = Var2, fill = value, label = round(value, 2))) +
        ggplot2::geom_tile() +
        ggplot2::geom_text() +
        ggplot2::labs(x = 'Estimated Signatures', y = 'True Signatures', fill = 'Cosine\nSimilarity')

    return(heatmap)
}
