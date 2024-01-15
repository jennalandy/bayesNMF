library(RcppHungarian)
reassign_signatures <- function(sim_mat) {
    if (nrow(sim_mat) == ncol(sim_mat)) {
        return(reassign_one(sim_mat))
    } else {
        transpose = FALSE
        if (ncol(sim_mat) > nrow(sim_mat)) {
            transpose = TRUE
            sim_mat = t(sim_mat)
        }

        # keep subset of rows that gives best final assignment
        to_try = combn(1:nrow(sim_mat), ncol(sim_mat))
        max_avg = 0
        for (i in 1:ncol(to_try)) {
            sim_mat_i = reassign_one(sim_mat[to_try[,i],])
            avg_i = mean(diag(sim_mat_i))
            if (avg_i > max_avg) {
                max_avg = avg_i
                max_sim_mat = sim_mat_i
            }
        }
        sim_mat = max_sim_mat

        if (transpose) {
            return(t(sim_mat))
        } else {
            return(sim_mat)
        }
    }
}
reassign_one <- function(sim_mat) {
    reassignment <- RcppHungarian::HungarianSolver(-1 * sim_mat)
    reassigned_sim_mat <- sim_mat[, reassignment$pairs[,2]]
    reassigned_sim_mat
}

niters = 1500
burn_in = 1000

niters_maxN = 5000
burn_in_maxN = 3500

source("setup_poisson.R")

test_that("nmf_normal_truncnormal works with 1 signature given N", {
    res <- nmf_normal_truncnormal(
        M, N = 1,
        file = "log_files/modelNT_dataP_N1",
        overwrite = TRUE,
        true_P = P,
        niters = niters,
        burn_in = burn_in
    )

    expect_equal(sum(is.na(res$final_Theta$P)), 0)
    expect_equal(sum(is.na(res$final_Theta$E)), 0)
})

test_that("nmf_normal_truncnormal works with Poisson data generating function given N", {
    res <- nmf_normal_truncnormal(
        M, N = 5,
        file = "log_files/modelNT_dataP_N5",
        overwrite = TRUE,
        true_P = P,
        niters = niters,
        burn_in = burn_in
    )

    expect_equal(sum(is.na(res$final_Theta$P)), 0)
    expect_equal(sum(is.na(res$final_Theta$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_truncnormal works with Poisson data generating function given max N", {
    res <- nmf_normal_truncnormal(
        M, max_N = 7,
        file = "log_files/modelNT_dataP_maxN7",
        overwrite = TRUE,
        true_P = P,
        niters = niters_maxN,
        burn_in = burn_in_maxN
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
    res <- nmf_normal_truncnormal(
        M, N = 5,
        file = "log_files/modelNT_dataPS_N5",
        overwrite = TRUE,
        true_P = P,
        niters = niters,
        burn_in = burn_in
    )

    expect_equal(sum(is.na(res$final_Theta$P)), 0)
    expect_equal(sum(is.na(res$final_Theta$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_truncnormal works with sparse Poisson data generating function given max N", {
    res <- nmf_normal_truncnormal(
        M, max_N = 7,
        file = "log_files/modelNT_dataPS_maxN7",
        overwrite = TRUE,
        true_P = P,
        niters = niters_maxN,
        burn_in = burn_in_maxN
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
    res <- nmf_normal_truncnormal(
        M, N = 5,
        file = "log_files/modelNT_dataN_N5",
        overwrite = TRUE,
        true_P = P,
        niters = niters,
        burn_in = burn_in
    )

    expect_equal(sum(is.na(res$final_Theta$P)), 0)
    expect_equal(sum(is.na(res$final_Theta$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_truncnormal works with Normal data generating function given max N", {
    res <- nmf_normal_truncnormal(
        M, max_N = 7,
        file = "log_files/modelNT_dataN_maxN7",
        overwrite = TRUE,
        true_P = P,
        niters = niters_maxN,
        burn_in = burn_in_maxN
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
    res <- nmf_normal_truncnormal(
        M, N = 5,
        file = "log_files/modelNT_dataNS_N5",
        overwrite = TRUE,
        true_P = P,
        niters = niters,
        burn_in = burn_in
    )

    expect_equal(sum(is.na(res$final_Theta$P)), 0)
    expect_equal(sum(is.na(res$final_Theta$E)), 0)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})

test_that("nmf_normal_truncnormal works with sparse Normal data generating function given max N", {
    res <- nmf_normal_truncnormal(
        M, max_N = 7,
        file = "log_files/modelNT_dataNS_maxN7",
        overwrite = TRUE,
        true_P = P,
        niters = niters_maxN,
        burn_in = burn_in_maxN
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)

    expect_lt(abs(sum(res$MAP$A) - 5), 2)

    sig_sims <- diag(reassign_signatures(res$sim_mat))
    sig_sims <- sig_sims[sig_sims != min(sig_sims)]
    expect_gt(min(sig_sims), 0.8)
})
