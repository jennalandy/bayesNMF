#' Log of Sum from Logs
#'
#' @param vec logged values
#'
#' @return log of summed values
sumLog <- function(vec) {
    ord <- sort(vec, decreasing=T)
    s <- ord[1]
    for (i in 2:length(vec)) {
        s <- s + log(1+exp(ord[i]-s))
    }
    return(s)
}

#' Estimate M from current values of Theta
#'
#' @param Theta list
#'
#' @return matrix
get_Mhat <- function(Theta) {
    Theta$P %*% Theta$E
}

#' Estimate multistudy M from current values of Theta
#'
#' @param Theta list of parameters
#' @param dims list of dimensions
#'
#' @return list of matrices
get_Mhat_multistudy <- function(Theta, dims) {
    lapply(1:dims$S, function(s) {
        Theta$P %*% diag(Theta$A[s,]) %*% Theta$E[[s]]
    })
}

#' get RMSE
#'
#' @param M mutational catalog matrix, K x G
#' @param M reconstructed mutational catalog matrix, K x G
#'
#' @return scalar
#' @export
get_RMSE <- function(M, Mhat) {
    sqrt(mean((M - Mhat)**2))
}

#' get RMSE in multistudy setting
#'
#' @param M list of mutational catalog matrices, length S
#' @param Mhat list of reconstructed mutational catalog matrices, length S
#' @param dims list of dimensions
#'
#' @return scalar
#' @export
get_RMSE_multistudy <- function(M, Mhat, dims) {
    M_wide = do.call(cbind, M)
    Mhat_wide = do.call(cbind, Mhat)

    sqrt(mean((M_wide - Mhat_wide)**2))
}

#' Get KL Divergence
#' @param M mutational catalog matrix, K x G
#' @param M reconstructed mutational catalog matrix, K x G
#'
#' @return scalar
#' @export
get_KLDiv <- function(M, Mhat) {
    Mhat[Mhat == 0] <- 1
    M[M == 0] <- 1
    sum(M * log(M / Mhat) - M + Mhat)
}

#' Get KL Divergence in the multistudy setting
#' @param M list of mutational catalog matrices, length S
#' @param Mhat list of reconstructed mutational catalog matrices, length S
#' @param dims list of dimensions
#'
#' @return scalar
#' @export
get_KLDiv_multistudy <- function(M, Mhat, dims) {
    M_wide = do.call(cbind, M)
    Mhat_wide = do.call(cbind, Mhat)

    Mhat_wide[Mhat_wide == 0] <- 1
    M_wide[M_wide == 0] <- 1
    sum(M_wide * log(M_wide / Mhat_wide) - M_wide + Mhat_wide)
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

    sim_mat = do.call(rbind, lapply(1:nrow(mat1), function(row_mat1) {
        sapply(1:nrow(mat2), function(row_mat2) {
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
