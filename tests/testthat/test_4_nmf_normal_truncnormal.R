source("setup_poisson.R")
source("test_funcs.R")

test_that("nmf_normal_truncnormal works with 5 signatures given N and fixed P", {
    res <- bayesNMF(
        M, N = 5,
        fixed = list(P = true_P),
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_N5_fixedP",
        overwrite = TRUE,
        true_P = true_P,
        convergence_control = new_convergence_control(maxiters = 2000)
    )

    expect_equal(res$MAP$P, true_P, tolerance = 1e-10)

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_equal(ncol(res$MAP$P), 5)
    expect_equal(nrow(res$MAP$P), 96)
    expect_equal(nrow(res$MAP$E), 5)
})


test_that("nmf_normal_truncnormal works with 5 signatures given N and fixed sigmasq", {
    res <- bayesNMF(
        M, N = 5,
        fixed = list(sigmasq = rep(5, 96)),
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_N5_fixedP",
        overwrite = TRUE,
        true_P = true_P,
        convergence_control = new_convergence_control(maxiters = 2000)
    )

    expect_equal(res$MAP$sigmasq, rep(5, 96), tolerance = 1e-10)

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_equal(ncol(res$MAP$P), 5)
    expect_equal(nrow(res$MAP$P), 96)
    expect_equal(nrow(res$MAP$E), 5)
})


test_that("nmf_normal_truncnormal works with 1 signature given N", {
    res <- bayesNMF(
        M, N = 1,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_N1",
        overwrite = TRUE,
        true_P = true_P,
        convergence_control = new_convergence_control(maxiters = 2000)
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)
    expect_equal(class(res$heatmap), c('gg','ggplot'))
})

test_that("nmf_normal_truncnormal works with 2 signatures given N", {
    res <- bayesNMF(
        M, N = 2,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_N2",
        overwrite = TRUE,
        true_P = true_P,
        convergence_control = new_convergence_control(maxiters = 2000)
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_equal(ncol(res$MAP$P), 2)
    expect_equal(nrow(res$MAP$P), 96)
    expect_equal(nrow(res$MAP$E), 2)
})

test_that("nmf_normal_truncnormal works with Poisson data generating function given N", {
    res <- bayesNMF(
        M, N = 5,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_N5",
        overwrite = TRUE,
        true_P = true_P,
        store_logs = TRUE,
        convergence_control = new_convergence_control(
            metric = 'logposterior'
        )
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

test_that("nmf_normal_truncnormal works with Poisson data generating function given max N", {
    res <- bayesNMF(
        M, max_N = 7,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_maxN7",
        overwrite = TRUE,
        true_P = true_P,
        store_logs = TRUE
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_lt(abs(sum(res$MAP$A) - 5), 1)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_truncnormal works with Poisson data generating function given max N and custom a, b", {
    res <- bayesNMF(
        M, max_N = 7,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_maxN7_customab",
        overwrite = TRUE,
        true_P = true_P,
        prior_parameters = list('a' = 0.8, 'b' = 0.4)
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_lt(abs(sum(res$MAP$A) - 5), 1)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_truncnormal works with Poisson data generating function given max N and recovery", {
    res <- bayesNMF(
        M, max_N = 7,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/modelNT_dataP_maxN7_recovery",
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
