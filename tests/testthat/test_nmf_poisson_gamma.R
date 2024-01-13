library(RcppHungarian)
reassign_signatures <- function(sim_mat) {
    reassignment <- RcppHungarian::HungarianSolver(-1 * sim_mat)
    reassigned_sim_mat <- sim_mat[, reassignment$pairs[,2]]
    reassigned_sim_mat
}

niters = 1500
burn_in = 1000

source("setup_poisson.R")

test_that("nmf_poisson_gamma works with 1 signature", {
    res <- nmf_poisson_gamma(
        M, N = 1,
        file = "log_files/modelPG_dataP_N1",
        overwrite = TRUE,
        true_P = P,
        niters = niters,
        burn_in = burn_in
    )

    expect_equal(sum(is.na(res$final_values$P)), 0)
    expect_equal(sum(is.na(res$final_values$E)), 0)
})

test_that("nmf_poisson_gamma works with Poisson data generating function", {
    res <- nmf_poisson_gamma(
        M, N = 5,
        file = "log_files/modelPG_dataP_N5",
        overwrite = TRUE,
        true_P = P,
        niters = niters,
        burn_in = burn_in
    )

    expect_equal(sum(is.na(res$final_values$P)), 0)
    expect_equal(sum(is.na(res$final_values$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

source("setup_poisson_sparse.R")

test_that("nmf_poisson_gamma works with sparse Poisson data generating function", {
    res <- nmf_poisson_gamma(
        M, N = 5,
        file = "log_files/modelPG_dataPS_N5",
        overwrite = TRUE,
        true_P = P,
        niters = niters,
        burn_in = burn_in
    )

    expect_equal(sum(is.na(res$final_values$P)), 0)
    expect_equal(sum(is.na(res$final_values$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

source("setup_normal.R")

test_that("nmf_poisson_gamma works with Normal data generating function", {
    res <- nmf_poisson_gamma(
        M, N = 5,
        file = "log_files/modelPG_dataN_N5",
        overwrite = TRUE,
        true_P = P,
        niters = niters,
        burn_in = burn_in
    )

    expect_equal(sum(is.na(res$final_values$P)), 0)
    expect_equal(sum(is.na(res$final_values$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

source("setup_normal_sparse.R")

test_that("nmf_poisson_gamma works with sparse Normal data generating function", {
    res <- nmf_poisson_gamma(
        M, N = 5,
        file = "log_files/modelPG_dataNS_N5",
        overwrite = TRUE,
        true_P = P,
        niters = niters,
        burn_in = burn_in
    )

    expect_equal(sum(is.na(res$final_values$P)), 0)
    expect_equal(sum(is.na(res$final_values$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)libra(✿◕‿◕)
})
