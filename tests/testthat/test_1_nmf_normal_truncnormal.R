source("setup_poisson.R")
source("test_funcs.R")

test_that("nmf_normal_truncnormal works with 1 signature given N", {
    res <- bayesNMF(
        M, N = 1,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_N1",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)
})

test_that("nmf_normal_truncnormal works with Poisson data generating function given N", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_N5",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_truncnormal works with Poisson data generating function given max N", {
    res <- bayesNMF(
        M, max_N = 7,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_maxN7",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_lt(abs(sum(res$MAP$A) - 5), 2)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

source("setup_poisson_sparse.R")

test_that("nmf_normal_truncnormal works with sparse Poisson data generating function given N", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataPS_N5",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_truncnormal works with sparse Poisson data generating function given max N", {
    res <- bayesNMF(
        M, max_N = 7,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataPS_maxN7",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_lt(abs(sum(res$MAP$A) - 5), 2)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

source("setup_normal.R")

test_that("nmf_normal_truncnormal works with Normal data generating function given N", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataN_N5",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_truncnormal works with Normal data generating function given max N", {
    res <- bayesNMF(
        M, max_N = 7,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataN_maxN7",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_lt(abs(sum(res$MAP$A) - 5), 2)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

source("setup_normal_sparse.R")

test_that("nmf_normal_truncnormal works with sparse Normal data generating function given N", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataNS_N5",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_truncnormal works with sparse Normal data generating function given max N", {
    res <- bayesNMF(
        M, max_N = 7,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataNS_maxN7",
        overwrite = TRUE,
        true_P = true_P
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_lt(abs(sum(res$MAP$A) - 5), 2)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})
