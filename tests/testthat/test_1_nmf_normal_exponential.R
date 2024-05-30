source("setup_poisson.R")
source("test_funcs.R")

test_that("nmf_normal_exponential works with 1 signature given N", {
    res <- bayesNMF(
        M, N = 1,
        likelihood = 'normal',
        prior = 'exponential',
        file = "log_files/modelNE_dataP_N1",
        true_P = true_P,
        overwrite = TRUE
    )

    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(sum(is.na(res$MAP$P)), 0)
})

test_that("nmf_normal_exponential works with 2 signature given N", {
    res <- bayesNMF(
        M, N = 2,
        likelihood = 'normal',
        prior = 'exponential',
        file = "log_files/modelNE_dataP_N2",
        true_P = true_P,
        overwrite = TRUE
    )

    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(sum(is.na(res$MAP$P)), 0)
})

test_that("nmf_normal_exponential works with Poisson data generating function given N", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'normal',
        prior = 'exponential',
        file = "log_files/modelNE_dataP_N5",
        true_P = true_P,
        overwrite = TRUE
    )

    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(sum(is.na(res$MAP$P)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_exponential works with Poisson data generating function given max_N", {
    res <- bayesNMF(
        M, max_N = 7,
        likelihood = 'normal',
        prior = 'exponential',
        file = "log_files/modelNE_dataP_maxN7",
        true_P = true_P,
        overwrite = TRUE
    )

    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(sum(is.na(res$MAP$P)), 0)

    expect_lt(abs(sum(res$MAP$A) - 5), 1)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_exponential works with Poisson data generating function given max_N and custom a, b", {
    res <- bayesNMF(
        M, max_N = 7,
        likelihood = 'normal',
        prior = 'exponential',
        file = "log_files/modelNE_dataP_maxN7_customab",
        true_P = true_P,
        overwrite = TRUE,
        prior_parameters = list('a' = 0.8, 'b' = 0.4)
    )

    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(sum(is.na(res$MAP$P)), 0)

    expect_lt(abs(sum(res$MAP$A) - 5), 1)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_exponential works with Poisson data generating function given max N and recovery", {
    res <- bayesNMF(
        M, max_N = 7,
        likelihood = 'normal',
        prior = 'exponential',
        file = "log_files/modelNE_dataP_maxN7_recovery",
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
