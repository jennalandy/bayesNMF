# These functions are copied over to be methods of the bayesNMF_sampler class

####################################
###### PRIVATE METHODS ##############
####################################

#' Fill hyperprior parameters
#' @description
#' Based on user-provided hyperprior parameters and default hyperprior parameters, create full hyperprior parameters matrices
#' 
#' @param self bayesNMF_sampler object
#' 
#' @return None, updates self$hyperprior_params
#' @noRd
fill_hyperprior_params_ <- function(self) {
  if (self$specs$prior == 'truncnormal') {
    # start with defaults
    hyperprior_params <- get_default_truncnorm_hyperprior_params_(self)
    
    # replace default with inits if applicable
    for (name in names(self$hyperprior_params)) {
      hyperprior_params[[name]] <- self$hyperprior_params[[name]]
    }

    # create matrices from values
    for (name in c("M","S","A","B")) {
      for (end in c("_p", "_e")) {
        mat = paste0(name, end)
        element = tolower(mat)
        nrow = ifelse(end == "_p", self$dims$K, self$dims$N)
        ncol = ifelse(end == "_p", self$dims$N, self$dims$G)
        hyperprior_params <- fill_matrix_(
          hyperprior_params,
          element = element, mat = mat,
          nrow = nrow, ncol = ncol
        )
      }
    }
  } else if (self$specs$prior == 'exponential') {
    # start with defaults
    hyperprior_params <- get_default_exponential_hyperprior_params_(self)

    # replace default with inits if applicable
    for (name in names(self$hyperprior_params)) {
      hyperprior_params[[name]] <- self$hyperprior_params[[name]]
    }

    # create matrices from values
    for (name in c("A","B")) {
      for (end in c("_p", "_e")) {
        mat = paste0(name, end)
        element = tolower(mat)
        nrow = ifelse(end == "_p", self$dims$K, self$dims$N)
        ncol = ifelse(end == "_p", self$dims$N, self$dims$G)
        hyperprior_params <- fill_matrix_(
          hyperprior_params,
          element = element, mat = mat,
          nrow = nrow, ncol = ncol
        )
      }
    }
  } else if (self$specs$prior == 'gamma') {
    # start with defaults
    hyperprior_params <- get_default_gamma_hyperprior_params_(self)

    # replace default with inits if applicable
    for (name in names(self$hyperprior_params)) {
      hyperprior_params[[name]] <- self$hyperprior_params[[name]]
    }

    # create matrices from values
    for (name in c("A","B","C","D")) {
      for (end in c("_p", "_e")) {
        mat = paste0(name, end)
        element = tolower(mat)
        nrow = ifelse(end == "_p", self$dims$K, self$dims$N)
        ncol = ifelse(end == "_p", self$dims$N, self$dims$G)
        hyperprior_params <- fill_matrix_(
          hyperprior_params,
          element = element, mat = mat,
          nrow = nrow, ncol = ncol
        )
      }
    }
  }
  # update hyperprior_params
  self$hyperprior_params <- hyperprior_params
}

#' Fill a matrix in hyperprior_params
#' @description
#' Helper function to fill a matrix in fill_hyperprior_params_
#' 
#' @param hyperprior_params list, hyperprior parameters
#' @param element string, element to fill
#' @param mat string, matrix to fill
#' @param nrow integer, number of rows
#' @param ncol integer, number of columns
#' 
#' @return list, updated hyperprior parameters
#' @noRd
fill_matrix_ <- function(hyperprior_params, element, mat, nrow, ncol) {
  if (
    element %in% names(hyperprior_params) &
    !(mat %in% names(hyperprior_params))
  ) {
    hyperprior_params[[mat]] <- matrix(
      hyperprior_params[[element]],
      nrow = nrow, ncol = ncol
    )
  }
  return(hyperprior_params)
}

#' Get default hyperprior parameters for truncated normal prior
#' @description
#' Helper function to fill_hyperprior_params_
#' 
#' @param self bayesNMF_sampler object
#' 
#' @return list, default hyperprior parameters
#' @noRd
get_default_truncnorm_hyperprior_params_ <- function(self) {
  m_p <- 0
  s_p <- sqrt(mean(self$data) / self$dims$N)
  a_p <- self$dims$N + 1
  b_p <- sqrt(self$dims$N)

  m_e <- 0
  s_e <- sqrt(mean(self$data) / self$dims$N)
  a_e <- self$dims$N + 1
  b_e <- sqrt(self$dims$N)

  return(list(
    m_p = m_p, s_p = s_p, a_p = a_p, b_p = b_p,
    m_e = m_e, s_e = s_e, a_e = a_e, b_e = b_e
  ))
}

#' Get default hyperprior parameters for exponential prior
#' @description
#' Helper function to fill_hyperprior_params_
#' 
#' @param self bayesNMF_sampler object
#' 
#' @return list, default hyperprior parameters
#' @noRd
get_default_exponential_hyperprior_params_ <- function(self) {
  a_p <- 10 * sqrt(self$dims$N)
  b_p <- 10 * sqrt(mean(self$data))
  a_e <- 10 * sqrt(self$dims$N)
  b_e <- 10 * sqrt(mean(self$data))

  return(list(
    a_p = a_p, b_p = b_p, a_e = a_e, b_e = b_e
  ))
}

#' Get default hyperprior parameters for gamma prior
#' @description
#' Helper function to fill_hyperprior_params_
#' 
#' @param self bayesNMF_sampler object
#' 
#' @return list, default hyperprior parameters
#' @noRd
get_default_gamma_hyperprior_params_ <- function(self) {
  a_p <- 10 * sqrt(self$dims$N)
  b_p <- 10
  c_p <- 10 * sqrt(mean(self$data))
  d_p <- 10
  a_e <- 10 * sqrt(self$dims$N)
  b_e <- 10
  c_e <- 10 * sqrt(mean(self$data))
  d_e <- 10

  return(list(
    a_p = a_p, b_p = b_p, c_p = c_p, d_p = d_p,
    a_e = a_e, b_e = b_e, c_e = c_e, d_e = d_e
  ))
}

