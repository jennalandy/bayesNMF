#' Log of Sum from Logs
#'
#' @param vec logged values
#'
#' @return log of summed values
#' @noRd
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
#' @noRd
get_Mhat <- function(Theta) {
    if (sum(Theta$A[1,]) == 0) {
        Mhat <- matrix(0, nrow = nrow(Theta$P), ncol = ncol(Theta$E))
    } else {
        Mhat <- Theta$P %*% diag(Theta$A[1,]) %*% Theta$E
    }
    return(Mhat)
}

#' Estimate M from current values of Theta excluding signature N
#'
#' @param Theta list
#' @param dims list of dimensions
#' @param n integer, signature to exclude
#'
#' @return matrix
#' @noRd
get_Mhat_no_n <- function(Theta, dims, n) {
    if (dims$N > 2) {
        Mhat_no_n = Theta$P[, -n] %*% diag(Theta$A[1, -n]) %*% Theta$E[-n, ]
    } else if (dims$N == 2) {
        Mhat_no_n = Theta$A[1, -n] * matrix(Theta$P[, -n], ncol = 1) %*% matrix(Theta$E[-n, ], nrow = 1)
    } else if (dims$N == 1) {
        Mhat_no_n = matrix(0, nrow = dims$K, ncol = dims$G)
    }
    return(Mhat_no_n)
}

#' Compute log prior
#'
#' @param Theta list of parameters
#' @param likelihood string, one of c('poisson','normal')
#' @param prior string, one of c('exponential','truncnormal','gamma')
#'
#' @return scalar
#' @noRd
get_logprior <- function(
    Theta, likelihood, prior
) {
    logprior = 0
    if (prior == 'truncnormal') {
        logprior <- logprior +
            sum(log(
                truncnorm::dtruncnorm(
                    Theta$P[,Theta$A[1,]==1], a = 0, b = Inf,
                    mean = Theta$Mu_p[,Theta$A[1,]==1],
                    sd = sqrt(Theta$Sigmasq_p[,Theta$A[1,]==1])
                )
            )) +
            sum(log(
                truncnorm::dtruncnorm(
                    Theta$E[Theta$A[1,]==1,], a = 0, b = Inf,
                    mean = Theta$Mu_e[Theta$A[1,]==1,],
                    sd = sqrt(Theta$Sigmasq_e[Theta$A[1,]==1,])
                )
            ))
    } else if (prior == 'exponential') {
        logprior <- logprior +
            sum(log(
                dexp(Theta$P[,Theta$A[1,]==1], Theta$Lambda_p[,Theta$A[1,]==1])
            )) +
            sum(log(
                dexp(Theta$E[Theta$A[1,]==1,], Theta$Lambda_e[Theta$A[1,]==1,])
            ))
    } else if (prior == 'gamma') {
        logprior <- logprior +
            sum(log(
                dgamma(Theta$P[,Theta$A[1,]==1], Theta$Alpha_p[,Theta$A[1,]==1], Theta$Beta_p[,Theta$A[1,]==1])
            )) +
            sum(log(
                dgamma(Theta$E[Theta$A[1,]==1,], Theta$Alpha_e[Theta$A[1,]==1,], Theta$Beta_p[Theta$A[1,]==1,])
            ))
    }

    if (likelihood == 'normal') {
        logprior <- logprior + sum(log(invgamma::dinvgamma(Theta$sigmasq, Theta$Alpha, Theta$Beta)))
    }
    return(logprior)
}

#' get Normal log likelihood
#'
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#'
#' @return scalar
#' @noRd
get_loglik_normal <- function(M, Theta, dims) {
    Mhat = get_Mhat(Theta)
    loglik = (- dims$K * sum(log(2 * pi * Theta$sigmasq)) / 2 -
        sum(colSums((M - Mhat)**2) / (2*Theta$sigmasq)))
    return(loglik)
}


#' get Poisson log likelihood
#'
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param logfac vector, logfac[i] = log(i!)
#'
#' @return scalar
#' @noRd
get_loglik_poisson <- function(M, Theta, dims, logfac) {
    Mhat = get_Mhat(Theta)
    Mhat[Mhat <= 0] <- 1
    M[M <= 0] <- 1
    sum(sapply(1:dims$K, function(k) {
        sum(sapply(1:dims$G, function(g) {
            -Mhat[k,g] +
            M[k,g] * log(Mhat[k,g]) -
            logfac[M[k,g]]
        }))
    }))
}

#' Compute RMSE
#'
#' @param M matrix, K x G
#' @param Mhat reconstructed matrix, K x G
#'
#' @return scalar
#' @noRd
get_RMSE <- function(M, Mhat) {
    sqrt(mean((M - Mhat)**2))
}

#' Compute KL Divergence
#'
#' @param M matrix, K x G
#' @param Mhat reconstructed matrix, K x G
#'
#' @return scalar
#' @noRd
get_KLDiv <- function(M, Mhat) {
    Mhat[Mhat <= 0] <- 1
    M[M <= 0] <- 1
    sum(M * log(M / Mhat) - M + Mhat)
}

#' compute BIC where number of parameters depends on likelihood-prior combination
#'
#' @param loglik scalar, log likelihood at Theta
#' @param Theta list, current values of all unknowns
#' @param dims list, named list of dimensions
#' @param likelihood string, one of c("normal", "poisson")
#' @param prior string, one of c("truncnormal","exponential","gamma")
#'
#' @return scalar, BIC
#' @noRd
get_BIC <- function(loglik, Theta, dims, likelihood, prior) {
    N = sum(Theta$A[1,])
    n_params = N * (dims$G + dims$K)

    return(n_params * log(dims$G) - 2 * loglik)
}

#' Get tempering schedule
#'
#' @param len integer, number of iterations
#'
#' @return vector of gamma values
#' @noRd
get_gamma_sched <- function(len = 5000, n_temp = 1000) {
    nX = round(n_temp / 374)
    gamma_sched <- c(
        rep(0, nX),
        c(sapply(9:5, function(x) {rep(10^(-x), nX)})),
        rep(10**(-4), round(8 * nX)),
        c(sapply(4:1, function(y) {
            c(sapply(seq(0, 8.9, by = 0.1), function(x) {
                rep((1+x)*10**(-y), nX)
            }))
        }))
    )
    gamma_sched <- c(gamma_sched, rep(1, len - length(gamma_sched)))
}


#' Pairwise cosine similarity between rows or columns of matrices
#'
#' @param mat1 matrix, first matrix for comparison
#' @param mat2 matrix, second matrix for comparison
#' @param name1 string, to name rows or cols of similarity matrix
#' @param name2 string, to name rows or cols of similarity matrix
#' @param which string, one of c("rows","cols")
#'
#' @return matrix
#' @export
pairwise_sim <- function(
        mat1, mat2,
        name1 = NULL,
        name2 = NULL,
        which = 'cols'
) {
    row_names = colnames(mat1)
    col_names = colnames(mat2)

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

    if (!is.null(name1)) {
        rownames(sim_mat) = paste0(name1, 1:nrow(sim_mat))
    } else {
        rownames(sim_mat) = row_names
    }

    if (!is.null(name2)) {
        colnames(sim_mat) = paste0(name2, 1:ncol(sim_mat))
    } else {
        colnames(sim_mat) = col_names
    }

    return(sim_mat)
}

#' Plot heatmap of cosine similarities
#' @description Plot a heatmap of cosine similarities between two matrices
#' with `ggplot2`. Can be similarity between rows or columns with the `which`
#' parameter.
#'
#' @param est_P estimated P (signatures matrix)
#' @param true_P true P (signatures matrix)
#' @param est_names names of estimated signatures
#' @param true_names names of true signatures
#' @param which string, one of c("rows","cols")
#'
#' @return ggplot object
#' @export
get_heatmap <- function(
    est_P, true_P,
    est_names = NULL,
    true_names = NULL,
    which = 'cols'
) {
    sim_mat <- pairwise_sim(
        est_P, true_P,
        name1 = est_names,
        name2 = true_names,
        which = which
    )
    sim_mat <- assign_signatures(sim_mat)

    sim_mat_melted <- reshape2::melt(sim_mat)

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

#' Assign signatures based on cosine similarities
#'
#' @param sim_mat similarity matrix
#'
#' @return matrix
#' @export
assign_signatures <- function(sim_mat) {
    reassignment <- RcppHungarian::HungarianSolver(-1 * sim_mat)
    reassigned_sim_mat <- sim_mat[, reassignment$pairs[,2]]
    if (nrow(sim_mat) == 1 | ncol(sim_mat) == 1) {
        reassigned_sim_mat = matrix(reassigned_sim_mat)
        colnames(reassigned_sim_mat) = colnames(sim_mat)[reassignment$pairs[,2]]
        rownames(reassigned_sim_mat) = rownames(sim_mat)
    }
    reassigned_sim_mat
}


#' Get mode of a list of matrices
#'
#' @param matrix_list list of matrices
#'
#' @return named list with mode ('matrix') and indices ('idx')
#' @noRd
get_mode <- function(matrix_list) {
    top_counts <- sapply(matrix_list, function(mat) {
        paste(c(mat), collapse = '')
    })
    str_counts <- table(top_counts)
    str_counts <- sort(str_counts, decreasing = TRUE)

    str_mode = names(str_counts)[1]
    idx_mode = which(top_counts == str_mode)
    matrix_mode = matrix_list[[idx_mode[1]]]

    return(list(
        'matrix' = matrix_mode,
        'top_counts' = str_counts[1:5],
        'idx' = idx_mode
    ))
}

#' Get element-wise mean of a list of matrices
#'
#' @param matrix_list list of matrices
#'
#' @return matrix
#' @noRd
get_mean <- function(matrix_list) {
    return(Reduce(`+`, matrix_list)/length(matrix_list))
}

#' Get element-wise quantiles of a list of matrices
#'
#' @param matrix_list list of matrices
#'
#' @return matrix
#' @noRd
get_quantile <- function(matrix_list, quantiles = c(0.025, 0.975)) {
    if (length(matrix_list) == 0) {
        return()
    }
    quantiles = sort(quantiles)
    quantile_matrices <- list()
    if (!("matrix" %in% class(matrix_list[[1]]))) {
        rows = length(matrix_list[[1]])
        cols = 1
        vector = TRUE
    } else {
        rows = nrow(matrix_list[[1]])
        cols = ncol(matrix_list[[1]])
        vector = FALSE
    }
    for (quantile in quantiles) {
        quantile_matrices[[as.character(quantile)]] <- matrix(
            nrow = rows,
            ncol = cols
        )
    }
    for (row in 1:rows) {
        for (col in 1:cols) {
            quants <- quantile(
                sapply(matrix_list, function(mat) {
                    if(vector) {
                        return(mat[row])
                    } else {
                        return(mat[row, col])
                    }
                }),
                quantiles
            )
            for (i in 1:length(quantiles)) {
                quantile = quantiles[i]
                quantile_matrices[[as.character(quantile)]][row, col] <- quants[i]
            }
        }
    }
    return(quantile_matrices)
}

#' Any named item in `fill_with` that is not specified in `list` gets
#' transferred into `list`. Final `list` is returned.
#'
#' @param list list of user specified values
#' @param fill_with list of default values
#'
#' @return updated list
#' @noRd
fill_list <- function(list, fill_with) {
    for (name in names(fill_with)) {
        if (!(name %in% names(list))) {
            list[[name]] = fill_with[[name]]
        }
    }
    return(list)
}

#' Create a matrix full of element values if element is provided and matrix is not
#'
#' @param Theta list of parameters
#' @param element name of element
#' @param matrix name of matrix
#' @param nrow row dimension of matrix
#' @param ncol column dimension of matrix
#'
#' @return Theta, list of p
#' @noRd
fill_matrix <- function(Theta, element, matrix, nrow, ncol) {
    if (element %in% names(Theta) & !(matrix %in% names(Theta))) {
        Theta[[matrix]] = matrix(Theta[[element]], nrow = nrow, ncol = ncol)
    }
    return(Theta)
}


#' Plot one metric
#'
#' @param x x-axis variable
#' @param y y-axis variable
#' @param vblue x-intercept for blue vertical line (convergence)
#' @param vgreen x-intercept for green vertical line (end of tempering)
#' @param xlab x-axis label
#' @param ylab y-axis label
#'
#' @return ggplot object
#' @noRd
plot_one <- function(x, y, vblue = NULL, vgreen = NULL, xlab = "", ylab = "") {
    if (sum(!is.na(y)) > 0){
        if (sum(y != -Inf, na.rm = TRUE) > 0) {
            plot(x, y, ylab = ylab, xlab = xlab)
            if (!is.null(vblue)) {
                abline(v = vblue, col = 'blue')
            }
            if (!is.null(vgreen)) {
                abline(v = vgreen, col = 'green')
            }
        }
    }
}

#' Plot all metrics and save in a pdf
#'
#' @param metrics list, list of metrics
#' @param plotfile string, file name (pdf) for plot
#' @param stop integer, iteration at which convergence was reached
#' @param learn_A boolean, whether A is being learned
#' @param gamma_sched numeric vector, tempering schedule
#' @param iter integer, current iteration
#' @param true_P matrix, true P matrix
#'
#' @return NULL
#' @noRd
plot_metrics <- function(metrics, plotfile, stop, learn_A, gamma_sched, iter, true_P) {
    if (!is.null(stop)) {
        vblue = stop
    } else { vblue = NULL }
    if (learn_A & gamma_sched[iter] == 1) {
        vgreen = which(gamma_sched == 1)[1]
    } else { vgreen = NULL }


    grDevices::pdf(plotfile)
    graphics::par(mfrow = c(3,1))

    x = unlist(metrics$sample_idx)
    plot_one(x, unlist(metrics$loglik), vblue, vgreen, xlab = "Iteration", ylab = "Log Likelihood")
    plot_one(x, unlist(metrics$logpost), vblue, vgreen, xlab = "Iteration", ylab = "Log Posterior")
    plot_one(x, unlist(metrics$BIC), vblue, vgreen, xlab = "Iteration", ylab = "BIC")
    plot_one(x, unlist(metrics$RMSE), vblue, vgreen, xlab = "Iteration", ylab = "RMSE")
    plot_one(x, unlist(metrics$KL), vblue, vgreen, xlab = "Iteration", ylab = "KL Divergence")

    plot_one(x, unlist(metrics$N), vblue, vgreen, xlab = "Iteration", ylab = "Latent Rank")
    plot_one(x, unlist(metrics$MAP_A_counts), vblue, vgreen, xlab = "Iteration", ylab = "MAP A Counts")
    if (!is.null(true_P)) {
        abline(h = ncol(true_P))
    }

    grDevices::dev.off()
}

#' Compute maximum a-posteriori (MAP) estimate of A, P, E, q
#'
#' @param logs list, list of Gibbs sampler logs
#' @param keep numeric vector, indices to consider for MAP
#' @param final boolean, whether this is the final iteration
#'
#' @return list, MAP estimate of A, P, E, q
#' @noRd
get_MAP <- function(logs, keep, final = FALSE) {

    # get MAP of A matrix (fine to do even if learn_A = FALSE)
    A_MAP = get_mode(logs$A[keep])
    map.idx = keep[A_MAP$idx]
    if (final) {
        keep_sigs <- which(A_MAP$matrix[1,] == 1)
    } else {
        keep_sigs <- 1:ncol(A_MAP$matrix)
    }

    # get MAP of P, E conditional on MAP of A
    MAP <- list(
        A = A_MAP$matrix,
        P = get_mean(logs$P[map.idx]),
        P_acceptance = get_mean(logs$P_acceptance[map.idx]),
        E = get_mean(logs$E[map.idx]),
        E_acceptance = get_mean(logs$E_acceptance[map.idx]),
        q = get_mean(logs$q[map.idx]),
        prob_inclusion = get_mean(logs$prob_inclusion[map.idx]),
        idx = map.idx,
        top_counts = A_MAP$top_counts
    )

    if (ncol(MAP$P) == 1) {
        MAP$P <- matrix(MAP$P, ncol = 1)
        MAP$E <- matrix(MAP$E, nrow = 1)
    }

    MAP$P <- MAP$P[,keep_sigs]
    MAP$E <- MAP$E[keep_sigs,]
    if ("sigmasq" %in% names(logs)) {
        MAP$sigmasq <- get_mean(logs$sigmasq[map.idx])
    }

    return(MAP)
}

#' Update list of metrics
#'
#' @param metrics list, list of metrics
#' @param MAP list, maximum a-posteriori estimates
#' @param iter integer, current iteration
#' @param Theta list, current state of Theta
#' @param M matrix, data matrix
#' @param likelihood string, one of c("normal", "poisson")
#' @param prior string, one of c("truncnormal","exponential","gamma")
#' @param dims list, named list of dimensions
#' @param logfac vector, logfac[i] = log(i!), only needed for `likelihood = 'poisson'`
#'
#' @return list, updated metrics
#' @noRd
update_metrics <- function(
        metrics, MAP, iter, Theta, M,
        likelihood, prior, dims, logfac
) {
    metrics$sample_idx[[iter]] <- iter

    Theta_MAP <- Theta
    Theta_MAP$P = MAP$P
    Theta_MAP$E = MAP$E
    Theta_MAP$A = MAP$A
    Theta_MAP$q = MAP$q
    if (likelihood == 'normal') {
        Theta_MAP$sigmasq = MAP$sigmasq
    }
    Mhat_MAP <- get_Mhat(Theta_MAP)

    # reconstruction errors: RMSE and KL
    metrics$RMSE[[iter]] <- get_RMSE(M, Mhat_MAP)
    metrics$KL[[iter]] <- get_KLDiv(M, Mhat_MAP)

    # likelihood-based metrics
    if (likelihood == 'normal') {
        metrics$loglik[[iter]] <- get_loglik_normal(M, Theta_MAP, dims)
    } else if (likelihood == 'poisson') {
        metrics$loglik[[iter]] <- get_loglik_poisson(M, Theta_MAP, dims, logfac)
    }
    metrics$N[[iter]] <- sum(Theta_MAP$A[1,])
    metrics$n_params[[iter]] <- metrics$N[[iter]] * (dims$G + dims$K + 2)
    metrics$BIC[[iter]] <- get_BIC(
        loglik = metrics$loglik[[iter]],
        Theta = Theta_MAP,
        dims = dims,
        likelihood = likelihood,
        prior = prior
    )

    metrics$logpost[[iter]] <- metrics$loglik[[iter]] + get_logprior(
        Theta_MAP, likelihood, prior
    )

    # top counts for MAP A
    metrics$MAP_A_counts[[iter]] <- MAP$top_counts[1]

    return(list(
        metrics = metrics,
        Theta_MAP = Theta_MAP,
        Mhat_MAP = Mhat_MAP
    ))
}

#' Log progress of Gibbs sampler
#'
#' @param iter integer, current iteration
#' @param done boolean, whether sampler is done
#' @param diff double, time since last log
#' @param convergence_control list, control parameters
#' @param convergence_status list, current status of convergence
#' @param gamma_sched numeric vector, tempering schedule
#' @param MAP list, maximum a-posteriori estimates
#' @param learn_A boolean, whether A is being learned
#'
#' @return NULL
#' @noRd
log_MAP <- function(iter, done, diff, convergence_control, convergence_status, gamma_sched, MAP, learn_A) {
    cat(paste(
        iter, "/", ifelse(done, "", "(up to)"), convergence_control$maxiters,
        "-", round(diff, 4), "seconds",
        ifelse(
            gamma_sched[iter] == 1,
            paste("-", paste0(round(convergence_status$prev_percent_change * 100, 4), "% change"),
                  "-", convergence_status$inarow_no_best, "no best",
                  "-", convergence_status$inarow_no_change, "no change"),
            ""
        ),
        "\n"
    ))

    if (learn_A) {
        print(MAP$top_counts)
        cat("\n")
    }
}

#' Log when sampler converges
#'
#' @param convergence_control list, control parameters
#' @param convergence_status list, current status of convergence
#'
#' @return NULL
#' @noRd
log_converged <- function(convergence_control, convergence_status) {
    cat(paste("\n\nCONVERGED at", convergence_status$best_iter))
    if (convergence_status$why %in% c("no best", "max iters")) {
        cat(paste(
            "\nNo best MAP since sample",
            convergence_status$best_iter, "\n\n"
        ))
    } else {
        cat(paste(
            "\nNo change in MAP over past",
            convergence_control$Ninarow_nochange, "samples\n\n"
        ))
    }
}

#' Compute 95% credible intervals for MAP estimates
#'
#' @param logs list, list of Gibbs sampler logs
#' @param map.idx integer vector, indices used for MAP estiamtes
#'
#' @return list, credible intervals for P, E, q, and sigmasq
#' @noRd
get_credible_intervals <- function(logs, map.idx) {
    credible_intervals <- list()
    credible_intervals[["P"]] <- get_quantile(logs$P[map.idx])
    credible_intervals[["E"]] <- get_quantile(logs$E[map.idx])
    credible_intervals[["q"]] <- get_quantile(logs$q[map.idx])
    if ("sigmasq" %in% names(logs)) {
        credible_intervals[["sigmasq"]] <- get_quantile(logs$sigmasq[map.idx])
    }
    return(credible_intervals)
}

#' Get Posterior pmf of latent rank
#'
#' @param A_list list, list of A matrices
#'
#' @return table, posterior pmf of latent rank
#' @noRd
get_posterior_counts_N <- function(A_list) {
    rank_logs = sapply(A_list, function(A) {sum(A)})
    rank_logs = factor(rank_logs, levels = 1:ncol(A_list[[1]]))
    posterior_counts <- table(rank_logs)
    return(posterior_counts)
}

#' Validate that provided N and max_N are valid
#'
#' @param N integer, fixed rank
#' @param max_N integer, maximum rank
#' @param recovery_priors list of prior parameters
#'
#' @return integer, N or max_N
#' @noRd
validate_N <- function(N, max_N, recovery_priors) {
    if (is.null(N) & is.null(max_N)) {
        stop("Either `N` or `max_N` must be provided.")
    } else if (!is.null(N) & !is.null(max_N)) {
        message("Both `N` and `max_N` provided, using `N`.")
        max_N = NULL
    } else if (is.null(N)) {
        N = max_N
    }
    if (length(recovery_priors) > 0) {
        N = N + recovery_priors$N_r
    }
    return(N)
}

#' Validate likelihood-prior combination
#'
#' @param likelihood string, one of c('normal','poisson')
#' @param prior string, one of c('truncnormal','exponential','gamma')
#' @param fast boolean, if `likelihood == 'poisson'` and `fast = TRUE`, updates
#' from the corresponding `likelihood == 'normal'` model are used as proposals
#' in an efficient Gibb's sampler
#'
#' @return NULL
#' @noRd
validate_model <- function(likelihood, prior, fast) {
    if (!(likelihood %in% c('normal', 'poisson'))) {
        stop("likelihood must be one of c('poisson','normal')")
    } else if (likelihood == 'normal') {
        if (!(prior %in% c('truncnormal','exponential'))) {
            stop("prior must be one of c('truncnormal','exponential') with `likelihood = 'normal'`")
        }
    } else if (likelihood == 'poisson') {
        if (!(prior %in% c('gamma','exponential','truncnormal'))) {
            stop("prior must be one of c('gamma','exponential','truncnormal') with `likelihood = 'poisson'`")
        }
        if (prior == 'gamma' & fast) {
            stop('gamma prior cannot be used with fast sampler')
        }
        if (prior == 'truncnormal' & !fast) {
            stop('truncnormal prior can only be used with fast sampler')
        }
    }
}
