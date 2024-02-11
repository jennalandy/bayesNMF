source("setup_poisson.R")
source("test_funcs.R")

test_that("nmf_normal_exponential works with one signature", {
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

test_that("nmf_normal_exponential works with Poisson data generating function", {
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

    log_post <- get_proportional_log_posterior(
        Theta = res$final_Theta,
        M = M,
        P = res$MAP$P,
        E = res$MAP$E,
        sigmasq = res$MAP$sigmasq,
        likelihood = 'normal',
        prior = 'exponential'
    )
    expect_true(!is.na(log_post))
})

source("setup_poisson_sparse.R")

test_that("nmf_normal_exponential works with sparse Poisson data generating function", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'normal',
        prior = 'exponential',
        file = "log_files/modelNE_dataPS_N5",
        true_P = true_P,
        overwrite = TRUE
    )

    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(sum(is.na(res$MAP$P)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

source("setup_normal.R")

test_that("nmf_normal_exponential works with Normal data generating function", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'normal',
        prior = 'exponential',
        file = "log_files/modelNE_dataN_N5",
        true_P = true_P,
        overwrite = TRUE
    )

    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(sum(is.na(res$MAP$P)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

source("setup_normal_sparse.R")

test_that("nmf_normal_exponential works with sparse Normal data generating function", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'normal',
        prior = 'exponential',
        file = "log_files/modelNE_dataNS_N5",
        true_P = true_P,
        overwrite = TRUE
    )

    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(sum(is.na(res$MAP$P)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

