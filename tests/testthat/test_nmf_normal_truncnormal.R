source("setup_normal.R")
library(RcppHungarian)
reassign_signatures <- function(sim_mat) {
    reassignment <- RcppHungarian::HungarianSolver(-1 * sim_mat)
    reassigned_sim_mat <- sim_mat[, reassignment$pairs[,2]]
    reassigned_sim_mat
}

test_that("nmf_normal_truncnormal works with 1 signature", {
    res <- nmf_normal_truncnormal(
        M, N = 5,
        file = "nmf_normal_truncnormal_onesig",
        overwrite = TRUE,
        true_P = P,
        niters = 1000,
        burn_in = 500
    )

    expect_equal(sum(is.na(res$final_values$P)), 0)
    expect_equal(sum(is.na(res$final_values$E)), 0)
})

test_that("nmf_normal_truncnormal works with Normal data generating function", {
    res <- nmf_normal_truncnormal(
        M, N = 5,
        file = "nmf_normal_truncnormal_normal_setup",
        overwrite = TRUE,
        true_P = P,
        niters = 1000,
        burn_in = 500
    )

    expect_equal(sum(is.na(res$final_values$P)), 0)
    expect_equal(sum(is.na(res$final_values$E)), 0)
    reassigned_sim_mat <-reassign_signatures(res$sim_mat)
    expect_gt(min(diag(reassigned_sim_mat)), 0.75)
})

source("setup_poisson.R")

test_that("nmf_normal_truncnormal works with Poisson data generating function", {
    res <- nmf_normal_truncnormal(
        M, N = 5,
        file = "nmf_normal_truncnormal_poisson_setup",
        overwrite = TRUE,
        true_P = P,
        niters = 1000,
        burn_in = 500
    )

    expect_equal(sum(is.na(res$final_values$P)), 0)
    expect_equal(sum(is.na(res$final_values$E)), 0)
    reassigned_sim_mat <- reassign_signatures(res$sim_mat)
    expect_gt(min(diag(reassigned_sim_mat)), 0.75)
})
