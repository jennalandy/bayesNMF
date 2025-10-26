# These functions are not a part of the sampler class
# They may be useful for other purposes

#' Broadcast a vector across a matrix
#' 
#' @description
#' read as "broadcast [vec] across [across] of [of] with [with]"
#' e.g., broadcast(v, "rows", M, "*") =
#'       broadcast vector v across rows of matrix M with multiplication
#'
#' @param vec vector to broadcast
#' @param across string, "rows" or "cols"
#' @param of matrix to broadcast across
#' @param with function to apply
#'
#' @return matrix
#' @export
broadcast <- function(vec, across, of, with) {
  margin <- if (across == "rows") 2L else 1L
  sweep(
    x = of,
    MARGIN = margin,
    STATS = vec,
    FUN = with
  )
}

#' Renormalize NMF results for factors to sum to 1 while maintaining product
#' 
#' @param P matrix, dimensions K x N, factors
#' @param E matrix, dimensions N x G, signatures
#' 
#' @return list, P and E renormalized for factors to sum to 1
#' @export
renormalize <- function(P, E) {
  E <- broadcast(
    colSums(P),
    across = "columns",
    of = E,
    with = "*"
  )
  P <- broadcast(
    colSums(P),
    across = "rows",
    of = P,
    with = "/"
  )
  return(list(P = P, E = E))
}

#' Get the mode of a list of matrices
#' 
#' @description
#' Find the mode of a list of matrices
#' by converting to strings and finding the most frequent string.
#' This is a simple way to find the most common matrix in a list,
#' for example the MAP estimate of a binary matrix.
#'
#' @param matrix_list list of matrices
#'
#' @return matrix
#' @export
get_mode <- function(matrix_list) {
  top_counts <- sapply(matrix_list, function(mat) {
    paste(c(mat), collapse = '')
  })
  str_counts <- table(top_counts)
  str_counts <- sort(str_counts, decreasing = TRUE)

  str_mode = names(str_counts)[1]
  idx_mode = which(top_counts == str_mode)
  matrix_mode = matrix_list[[idx_mode[1]]]

  return(list(
    'matrix' = matrix_mode,
    'top_counts' = str_counts[1:5],
    'idx' = idx_mode
  ))
}

#' Convert a table to a string for logging
#'
#' @param tab table to convert
#'
#' @return string
#' @export
log_table <- function(tab) {
  table_str_lines <- as.character(
    knitr::kable(tab, format = "simple", align = 'c')
  )

  # remove any lines with with only spaces and "-"
  table_str_lines <- sapply(table_str_lines, function(line) {
    trimws(gsub("-", "", line))
  })
  table_str_lines <- table_str_lines[table_str_lines != ""]

  table_str <- paste(table_str_lines, collapse = "\n")
  return(table_str)
}

#' Update a list with a new value at a given index, and remove the first element if the list is too long
#' 
#' @param list list to update
#' @param new_value value to update
#' @param index index to update
#' @param max_length maximum length of list
#' 
#' @return list
#' @export
update_list <- function(list, new_value, index, max_length) {
  if (index <= max_length) {
    list[[index]] <- new_value
  } else {
    list[[1]] <- NULL
    list[[max_length]] <- new_value
  }
  return(list)
}

#' Check if two matrices are equal
#' @description
#' Returns TRUE if two matrices are both non-null and have the same dimensions and values.
#' 
#' @param mat1 matrix
#' @param mat2 matrix
#' 
#' @return boolean
#' @export
check_matrix_equal <- function(mat1, mat2) {
  if (is.null(mat1) | is.null(mat2)) {
    return(FALSE)
  }
  if (nrow(mat1) != nrow(mat2) | ncol(mat1) != ncol(mat2)) {
    return(FALSE)
  }
  return(all(mat1 == mat2))
}

#' Check if two vectors are equal
#' @description
#' Returns TRUE if two vectors are both non-null and have the same length and values.
#' 
#' @param vec1 vector
#' @param vec2 vector
#' 
#' @return boolean
#' @export
check_vector_equal <- function(vec1, vec2) {
  if (is.null(vec1) | is.null(vec2)) {
    return(FALSE)
  }
  if (length(vec1) != length(vec2)) {
    return(FALSE)
  }
  return(all(vec1 == vec2))
}

#' Get COSMIC reference SBS signatures matrix v3.3.1
#'
#' @return matrix
#' @export
#'
#' @examples
#' cosmic <- get_cosmic()
get_cosmic <- function() {
  P <- read.csv(system.file("extdata", "COSMIC_v3.3.1_SBS_GRCh37.csv", package = "bayesNMF"), row.names = 1)
  return(P)
}

#' Download COSMIC reference SBS signatures matrix v3.3.1 from the web
#'
#' @return matrix
#' @export
#'
#' @examples
#' cosmic <- download_cosmic()
download_cosmic <- function() {
  P <- read.csv(
    "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3.1_SBS_GRCh37.txt", # nolint: line_length_linter.
    sep = "\t"
  )
  rownames(P) <- P$Type
  P <- P[, 2:ncol(P)]
  P <- as.matrix(P)
  return(P)
}

#' Get a named list of colors used for COSMIC signature plots
#'
#' @return list
#' @export
#'
#' @examples
#' COSMIC_colors <- get_cosmic_colors()
get_cosmic_colors <- function() {
  COSMIC_colors <- c(
    'C>A' = rgb(8, 181, 236, maxColorValue = 255),
    'C>G' = rgb(0, 0, 0, maxColorValue = 255),
    'C>T' = rgb(225, 37, 33, maxColorValue = 255),
    'T>A' = rgb(198, 193, 195, maxColorValue = 255),
    'T>C' = rgb(153, 200, 87, maxColorValue = 255),
    'T>G' = rgb(233, 190, 189, maxColorValue = 255)
  )
  return(COSMIC_colors)
}

#' Pairwise cosine similarity between rows or columns of matrices
#'
#' @param mat1 matrix, first matrix for comparison
#' @param mat2 matrix, second matrix for comparison or "cosmic" to use COSMIC reference SBS signatures v3.3.1
#' @param name1 string, to name rows or cols of similarity matrix
#' @param name2 string, to name rows or cols of similarity matrix
#' @param which string, which dimension to compute similarity over, one of c("rows","cols")
#'
#' @return matrix
#' @export
pairwise_sim <- function(
  mat1, mat2 = "cosmic",
  name1 = NULL,
  name2 = NULL,
  which = "cols"
) {
  if ("character" %in% class(mat2)) {
    if (mat2 == "cosmic") {
      mat2 = get_cosmic()
    } else {
      stop("Parameter `reference_P` must be a matrix or 'cosmic'")
    }
  }

  row_names <- colnames(mat1)
  col_names <- colnames(mat2)

  if (which == "cols") {
    mat1 <- t(mat1)
    mat2 <- t(mat2)
  }

  if (ncol(mat1) != ncol(mat2)) {
    overlap_dim <- ifelse(which == "rows", "cols", "rows")
    stop(paste0(
      "Different number of ", overlap_dim, ": ",
      ncol(mat1), " != ", ncol(mat2)
    ))
  }

  sim_mat = do.call(rbind, lapply(seq_len(nrow(mat1)), function(row_mat1) {
    sapply(seq_len(nrow(mat2)), function(row_mat2) {
      lsa::cosine(mat1[row_mat1, ], mat2[row_mat2, ])
    })
  }))

  if (!is.null(name1)) {
    rownames(sim_mat) <- paste0(name1, seq_len(nrow(sim_mat)))
  } else {
    rownames(sim_mat) <- row_names
  }

  if (!is.null(name2)) {
    colnames(sim_mat) <- paste0(name2, seq_len(ncol(sim_mat)))
  } else {
    colnames(sim_mat) <- col_names
  }

  return(sim_mat)
}

#' Assign two sets of factors based on cosine similarities
#' @description 
#' Relies on the Hungarian algorithm (RcppHungarian::HungarianSolver) to maximize total cosine similarity between estimated and assigned reference signatures.
#'
#' @param estimated_P matrix, estimated factors
#' @param reference_P matrix, reference factors or "cosmic" to use COSMIC reference SBS signatures v3.3.1
#' @param which string, which dimension to compute similarity over, one of c("rows","cols") (default "cols")
#' @param keep_all_est boolean, whether to keep all estimated factors (default TRUE)
#' @param keep_all_ref boolean, whether to keep all reference factors (default FALSE)
#' @param return_mat boolean, whether to return the similarity matrix (default FALSE)
#'
#' @return aligned similarity matrix with rows as estimated signatures and columns as reference signatures if return_mat is TRUE, otherwise data frame with columns:
#' \itemize{
#'   \item sig_est: estimated signature index or name
#'   \item sig_ref: assigned reference signature index or name
#'   \item cos_sim: cosine similarity between the estimated and assigned reference signatures
#' }
#' @export
hungarian_assignment <- function(
  estimated_P,
  reference_P = "cosmic",
  which = "cols",
  keep_all_est = TRUE,
  keep_all_ref = FALSE,
  return_mat = FALSE,
  check_reference_order = TRUE
) {
  if ("character" %in% class(reference_P)) {
    if (reference_P == "cosmic") {
      reference_P = get_cosmic()
    } else {
      stop("Parameter `reference_P` must be a matrix or 'cosmic'")
    }
  }

  # reorder reference_P to match the order of the data if available
  if (check_reference_order) {
    if (!is.null(rownames(estimated_P))) {
      if (!is.null(rownames(reference_P))) {
        if (!all(rownames(estimated_P) == rownames(reference_P))) {
          if (!setequal(rownames(estimated_P), rownames(reference_P))) {
            warning("Row names of estimated_P and reference_P do not overlap. Reference matrix will not be reordered.")
          } else {
            reference_P <- reference_P[rownames(estimated_P), ]
          }
        } # else they're already in the same order
      } else {
        warning("Row names of reference_P are not available. Reference matrix will not be reordered.")
      }
    } else {
      warning("Row names of estimated_P are not available. Reference matrix will not be reordered.")
    }
  }
  

  sim_mat <- pairwise_sim(estimated_P, reference_P, which = which)

  # give column and row names to maintain labeling
  if (is.null(colnames(reference_P))) {
    rename_ref <- TRUE
    colnames(sim_mat) <- paste0("Ref", seq_len(ncol(reference_P)))
  } else {
    rename_ref <- FALSE
    colnames(sim_mat) <- colnames(reference_P)
  }
  if (is.null(colnames(estimated_P))) {
    rename_est <- TRUE
    rownames(sim_mat) <- paste0("Est", seq_len(ncol(estimated_P)))
  } else {
    rename_est <- FALSE
    rownames(sim_mat) <- colnames(estimated_P)
  }

  # use Hungarian algorithm to assign estimated signatures to reference signatures
  reassignment <- RcppHungarian::HungarianSolver(-1 * sim_mat)
  colnames(reassignment$pairs) <- c("row","col") 
  reassignment$pairs <- data.frame(reassignment$pairs) %>% dplyr::filter( row != 0, col != 0 )
  rows <- reassignment$pairs[, 1]
  cols <- reassignment$pairs[, 2]
  if (keep_all_est) {
    for (row in setdiff(seq_len(nrow(sim_mat)), rows)) {
      rows <- c(rows, row)
    }
  } 
  if (keep_all_ref) {
    for (col in setdiff(seq_len(ncol(sim_mat)), cols)) {
      cols <- c(cols, col)
    }
  }
  reassigned_sim_mat <- sim_mat[rows, cols]

  # if one dimension is 1, convert to matrix
  if (nrow(sim_mat) == 1 || ncol(sim_mat) == 1) {
    reassigned_sim_mat <- matrix(reassigned_sim_mat)
    colnames(reassigned_sim_mat) <- colnames(sim_mat)[cols]
    rownames(reassigned_sim_mat) <- rownames(sim_mat)[rows]
  }

  # if one dimension is greater, fill with cos sim of 0 to make square matrix
  if (nrow(reassigned_sim_mat) > ncol(reassigned_sim_mat)) {
    fill <- matrix(0, nrow = nrow(reassigned_sim_mat), ncol = nrow(reassigned_sim_mat) - ncol(reassigned_sim_mat))
    colnames(fill) <- rep("None", nrow(reassigned_sim_mat) - ncol(reassigned_sim_mat))
    reassigned_sim_mat <- cbind(reassigned_sim_mat, fill)
  }
  if (ncol(reassigned_sim_mat) > nrow(reassigned_sim_mat)) {
    fill <- matrix(0, nrow = ncol(reassigned_sim_mat) - nrow(reassigned_sim_mat), ncol = ncol(reassigned_sim_mat))
    rownames(fill) <- rep("None", ncol(reassigned_sim_mat) - nrow(reassigned_sim_mat))
    reassigned_sim_mat <- rbind(reassigned_sim_mat, fill)
  }

  if (return_mat) {
    return(reassigned_sim_mat)
  }

  assignment_df <- data.frame(
    sig_est = rownames(reassigned_sim_mat),
    sig_ref = colnames(reassigned_sim_mat),
    cos_sim = diag(reassigned_sim_mat)
  )

  # rename to original names if needed
  if (rename_est) {
    assignment_df$sig_est <- as.numeric(stringr::str_remove(assignment_df$sig_est, "Est"))
  }
  if (rename_ref) {
    assignment_df$sig_ref <- as.numeric(stringr::str_remove(assignment_df$sig_ref, "Ref"))
  }

  return(assignment_df)
}