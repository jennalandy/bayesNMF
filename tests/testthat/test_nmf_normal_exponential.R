source("setup_poisson.R")
library(RcppHungarian)
reassign_signatures <- function(sim_mat) {
    reassignment <- RcppHungarian::HungarianSolver(-1 * sim_mat)
    reassigned_sim_mat <- sim_mat[, reassignment$pairs[,2]]
    reassigned_sim_mat
}

test_that("nmf_normal_exponential works with one signature", {
    res <- nmf_normal_exponential(
        M, N = 1,
        file = "nmf_normal_exponential_onesig",
        overwrite = TRUE,
        true_P = P,
        niters = 1000,
        burn_in = 500
    )

    expect_equal(sum(is.na(res$final_values$E)), 0)
    expect_equal(sum(is.na(res$final_values$P)), 0)
})

test_that("nmf_normal_exponential works with Poisson data generating function", {
    res <- nmf_normal_exponential(
        M, N = 5,
        file = "nmf_normal_exponential_poisson_setup",
        overwrite = TRUE,
        true_P = P,
        niters = 1000,
        burn_in = 500
    )

    expect_equal(sum(is.na(res$final_values$E)), 0)
    expect_equal(sum(is.na(res$final_values$P)), 0)
    reassigned_sim_mat <-reassign_signatures(res$sim_mat)
    expect_gt(min(diag(reassigned_sim_mat)), 0.75)
})

source("setup_normal.R")

test_that("nmf_normal_exponential works with Normal data generating function", {
    res <- nmf_normal_exponential(
        M, N = 5,
        file = "nmf_normal_exponential_normal_setup",
        overwrite = TRUE,
        true_P = P,
        niters = 5000,
        burn_in = 4000
    )

    expect_equal(sum(is.na(res$final_values$E)), 0)
    expect_equal(sum(is.na(res$final_values$P)), 0)
    reassigned_sim_mat <-reassign_signatures(res$sim_mat)
    expect_gt(min(diag(reassigned_sim_mat)), 0.75)
})

