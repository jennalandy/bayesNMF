source("setup_poisson.R")
source("test_funcs.R")

test_that("nmf_poisson_gamma works with 1 signature", {
    res <- bayesNMF(
        M, N = 1,
        likelihood = 'poisson',
        prior = 'gamma',
        file = "log_files/modelPG_dataP_N1",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)
})

test_that("nmf_poisson_gamma works with fixed N", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'poisson',
        prior = 'gamma',
        file = "log_files/modelPG_dataP_N5",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_poisson_gamma works with max N", {
    res <- bayesNMF(
        M, max_N = 7,
        likelihood = 'poisson',
        prior = 'gamma',
        file = "log_files/modelPG_dataP_N5",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_poisson_gamma works with Poisson data generating function given max N and recovery", {
    res <- bayesNMF(
        M, max_N = 7,
        likelihood = 'poisson',
        prior = 'gamma',
        file = "log_files/modelPG_dataP_maxN7_recovery",
        overwrite = TRUE,
        true_P = true_P,
        recovery = TRUE
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_lt(abs(sum(res$MAP$A) - 5), 1)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})
