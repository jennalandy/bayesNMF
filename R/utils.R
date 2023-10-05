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

#' Pairwise cosine similarity between rows or columns of matrices
#'
#' @param mat1 matrix, first matrix for comparison
#' @param mat2 matrix, second matrix for comparison
#' @param name1 string, to name rows of similarity matrix
#' @param name2 string, to name columns of similarity matrix
#' @param which string, one of c("rows","cols")
#'
#' @return matrix
#' @export
pairwise_sim <- function(
        mat1, mat2,
        name1 = '',
        name2 = '',
        which = 'rows'
) {
    if (which == 'cols') {
        mat1 = t(mat1)
        mat2 = t(mat2)
    }

    if (ncol(mat1) != ncol(mat2)) {
        overlap_dim = ifelse(which == "rows","cols","rows")
        stop(paste0(
            "Different number of ", overlap_dim, ": ",
            ncol(mat1), " != ", ncol(mat2)
        ))
    }

    rows = nrow(mat1)
    sim_mat = do.call(rbind, lapply(1:rows, function(row_mat1) {
        sapply(1:rows, function(row_mat2) {
            lsa::cosine(mat1[row_mat1,], mat2[row_mat2,])
        })
    }))

    if (name1 != "") {
        rownames(sim_mat) = paste0(name1, 1:rows)
    }

    if (name2 != "") {
        colnames(sim_mat) = paste0(name2, 1:rows)
    }

    return(sim_mat)
}

#' Get heatmap
#'
#' @param est_P estimated P (signatures matrix)
#' @param true_P true P (signatures matrix)
#'
#' @return ggplot object
#' @export
get_heatmap <- function(est_P, true_P) {
    sim_mat <- pairwise_sim(est_P, true_P, which = 'cols')

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
