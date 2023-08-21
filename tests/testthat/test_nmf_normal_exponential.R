source("setup.R")

test_that("nmf_normal_exponential works", {
    res <- nmf_normal_exponential(M, N = 5, overwrite = TRUE, true_Theta = true_Theta)
    reassigned_sim_mat <-reassign_signatures(res$sim_mat)
    expect_gt(min(diag(reassigned_sim_mat)), 0.75)
})
