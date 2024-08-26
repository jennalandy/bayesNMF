source("setup_poisson.R")
source("test_funcs.R")

small_test_convergence_control <- new_convergence_control(maxiters = 1000, MAP_over = 500)
large_test_convergence_control <- new_convergence_control(maxiters = 2000, MAP_over = 500)

test_that("nmf_normal_exponential works with 1 signature given N", {
    res <- bayesNMF(
        M, rank = 1,
        likelihood = 'normal',
        prior = 'exponential',
        file = "log_files/modelNE_dataP_N1",
        true_P = true_P,
        overwrite = TRUE,
        convergence_control = small_test_convergence_control
    )

    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(sum(is.na(res$MAP$P)), 0)
})

test_that("nmf_normal_exponential works with 2 signature given N", {
    res <- bayesNMF(
        M, rank = 2,
        likelihood = 'normal',
        prior = 'exponential',
        file = "log_files/modelNE_dataP_N2",
        true_P = true_P,
        overwrite = TRUE,
        convergence_control = small_test_convergence_control
    )

    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(sum(is.na(res$MAP$P)), 0)
})

test_that("nmf_normal_exponential learns signatures given N", {
    res <- bayesNMF(
        M, rank = 5,
        likelihood = 'normal',
        prior = 'exponential',
        file = "log_files/modelNE_dataP_N5",
        true_P = true_P,
        overwrite = TRUE,
        convergence_control = small_test_convergence_control
    )

    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(sum(is.na(res$MAP$P)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_exponential learns rank and signatures", {
    res <- bayesNMF(
        M, rank = 1:7,
        likelihood = 'normal',
        prior = 'exponential',
        file = "log_files/modelNE_dataP_maxN7",
        true_P = true_P,
        overwrite = TRUE,
        convergence_control = large_test_convergence_control
    )

    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(sum(is.na(res$MAP$P)), 0)

    expect_lt(abs(sum(res$MAP$A) - 5), 1)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})
