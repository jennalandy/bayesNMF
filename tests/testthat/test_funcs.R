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
