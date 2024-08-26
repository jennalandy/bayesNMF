source("setup_poisson.R")
source("test_funcs.R")

small_test_convergence_control <- new_convergence_control(maxiters = 1000, MAP_over = 500)
large_test_convergence_control <- new_convergence_control(maxiters = 2000, MAP_over = 500)

test_that("nmf_normal_truncnormal works with 1 signature given N", {
    res <- bayesNMF(
        M, rank = 1,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_N1",
        overwrite = TRUE,
        true_P = true_P,
        convergence_control = small_test_convergence_control
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(class(res$heatmap), c('gg','ggplot'))
})

test_that("nmf_normal_truncnormal works with 2 signatures given N", {
    res <- bayesNMF(
        M, rank = 2,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_N2",
        overwrite = TRUE,
        true_P = true_P,
        convergence_control = small_test_convergence_control
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_equal(ncol(res$MAP$P), 2)
    expect_equal(nrow(res$MAP$P), 96)
    expect_equal(nrow(res$MAP$E), 2)
})

test_that("nmf_normal_truncnormal learns signatures", {
    res <- bayesNMF(
        M, rank = 5,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_N5",
        overwrite = TRUE,
        true_P = true_P,
        store_logs = TRUE,
        convergence_control = small_test_convergence_control
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_equal(ncol(res$MAP$P), 5)
    expect_equal(nrow(res$MAP$P), 96)
    expect_equal(nrow(res$MAP$E), 5)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_truncnormal learns signatures and rank", {
    res <- bayesNMF(
        M, rank = 1:7,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_maxN7",
        overwrite = TRUE,
        true_P = true_P,
        store_logs = TRUE,
        convergence_control = large_test_convergence_control
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_lt(abs(sum(res$MAP$A) - 5), 1)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})
