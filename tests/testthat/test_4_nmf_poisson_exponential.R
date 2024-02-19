source("setup_poisson.R")
source("test_funcs.R")

test_that("nmf_poisson_exponential works with 1 signature", {
    res <- bayesNMF(
        M, N = 1,
        likelihood = 'poisson',
        prior = 'exponential',
        file = "log_files/modelPE_dataP_N1",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)
})

test_that("nmf_poisson_exponential works with 2 signatures", {
    res <- bayesNMF(
        M, N = 2,
        likelihood = 'poisson',
        prior = 'exponential',
        file = "log_files/modelPE_dataP_N",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)
})

test_that("nmf_poisson_exponential works with Poisson data generating function given N", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'poisson',
        prior = 'exponential',
        file = "log_files/modelPE_dataP_N5",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)

    log_post <- get_proportional_log_posterior(
        Theta = res$final_Theta,
        M = M,
        P = res$MAP$P,
        E = res$MAP$E,
        sigmasq = res$MAP$sigmasq,
        likelihood = 'poisson',
        prior = 'exponential'
    )
    expect_true(!is.na(log_post))
})

test_that("nmf_poisson_exponential works with Poisson data generating function given max_N", {
    res <- bayesNMF(
        M, max_N = 7,
        likelihood = 'poisson',
        prior = 'exponential',
        file = "log_files/modelPE_dataP_maxN7",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)

    log_post <- get_proportional_log_posterior(
        Theta = res$final_Theta,
        M = M,
        P = res$MAP$P,
        E = res$MAP$E,
        sigmasq = res$MAP$sigmasq,
        likelihood = 'poisson',
        prior = 'exponential'
    )
    expect_true(!is.na(log_post))
})

source("setup_poisson_sparse.R")

test_that("nmf_poisson_exponential works with sparse Poisson data generating function", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'poisson',
        prior = 'exponential',
        file = "log_files/modelPE_dataPS_N5",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

source("setup_normal.R")

test_that("nmf_poisson_exponential works with Normal data generating function", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'poisson',
        prior = 'exponential',
        file = "log_files/modelPE_dataN_N5",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

source("setup_normal_sparse.R")

test_that("nmf_poisson_exponential works with sparse Normal data generating function", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'poisson',
        prior = 'exponential',
        file = "log_files/modelPE_dataNS_N5",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})
