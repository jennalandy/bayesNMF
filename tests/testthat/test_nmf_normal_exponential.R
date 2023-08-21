source("setup_normal.R")

test_that("nmf_normal_exponential works with Normal data generating function", {
    res <- nmf_normal_exponential(
        M, N = 5,
        file = "nmf_normal_exponential_normal_setup",
        overwrite = TRUE,
        true_P = P
    )
    reassigned_sim_mat <-reassign_signatures(res$sim_mat)
    expect_gt(min(diag(reassigned_sim_mat)), 0.75)
})

source("setup_poisson.R")

test_that("nmf_normal_exponential works with Poisson data generating function", {
    res <- nmf_normal_exponential(
        M, N = 5,
        file = "nmf_normal_exponential_poisson_setup",
        overwrite = TRUE,
        true_P = P
    )

    reassigned_sim_mat <-reassign_signatures(res$sim_mat)
    expect_gt(min(diag(reassigned_sim_mat)), 0.75)
})
