source("setup_poisson.R")
source("test_funcs.R")

small_test_convergence_control <- new_convergence_control(maxiters = 1000, MAP_over = 500)
large_test_convergence_control <- new_convergence_control(maxiters = 2000, MAP_over = 500)

test_that("nmf_poisson_gamma works with 1 signature", {
    res <- bayesNMF(
        M, rank = 1,
        likelihood = 'poisson',
        prior = 'gamma',
        file = "log_files/modelPG_dataP_N1",
        overwrite = TRUE,
        true_P = true_P,
        convergence_control = small_test_convergence_control
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)
})

test_that("nmf_poisson_gamma works with 2 signatures", {
    res <- bayesNMF(
        M, rank = 2,
        likelihood = 'poisson',
        prior = 'gamma',
        file = "log_files/modelPG_dataP_N1",
        overwrite = TRUE,
        true_P = true_P,
        convergence_control = small_test_convergence_control
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)
})

test_that("nmf_poisson_gamma learns signatures", {
    res <- bayesNMF(
        M, rank = 5,
        likelihood = 'poisson',
        prior = 'gamma',
        file = "log_files/modelPG_dataP_N5",
        overwrite = TRUE,
        true_P = true_P,
        convergence_control = small_test_convergence_control
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_poisson_gamma learns signatures and rank", {
    res <- bayesNMF(
        M, rank = 1:7,
        likelihood = 'poisson',
        prior = 'gamma',
        file = "log_files/modelPG_dataP_N5",
        overwrite = TRUE,
        true_P = true_P,
        convergence_control = large_test_convergence_control
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})
