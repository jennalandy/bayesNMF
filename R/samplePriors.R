#### Truncated Normal Prior

#' Sample Prior Parameter: Mu for Pn
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length K
#' @noRd
sample_Mu_Pn <- function(n, Theta, dims, gamma) {
    num <- Theta$M_p[,n]/Theta$S_p[,n] + gamma*Theta$P[,n]/Theta$Sigmasq_p[,n]
    denom <- 1/Theta$S_p[,n] + gamma*1/Theta$Sigmasq_p[,n]
    rnorm(dims$K, num/denom, 1/denom)
}

#' Sample Prior Parameter: Mu for En
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_Mu_En <- function(n, Theta, dims, gamma) {
    num <- Theta$M_e[n,]/Theta$S_e[n,] + gamma*Theta$E[n,]/Theta$Sigmasq_e[n,]
    denom <- 1/Theta$S_e[n,] + gamma*1/Theta$Sigmasq_e[n,]
    rnorm(dims$G, num/denom, 1/denom)
}

#' Sample Prior Parameter: Sigmasq for Pn
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length K
#' @noRd
sample_Sigmasq_Pn <- function(n, Theta, dims, gamma) {
    invgamma::rinvgamma(
        n = dims$K,
        shape = Theta$A_p[,n] + gamma*1/2,
        rate = Theta$B_p[,n] + gamma*(Theta$P[,n] - Theta$Mu_p[,n])**2/2
    )
}

#' Sample Prior Parameter: Sigmasq for En
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_Sigmasq_En <- function(n, Theta, dims, gamma) {
    invgamma::rinvgamma(
        n = dims$G,
        shape = Theta$A_e[n,] + gamma*1/2,
        rate = Theta$A_e[n,] + gamma*(Theta$E[n,] - Theta$Mu_e[n,])**2/2
    )
}

#### Exponential Prior

#' Sample Prior Parameter: Lambda for Pn
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length K
#' @noRd
sample_Lambda_Pn <- function(n, Theta, dims, gamma) {
    rgamma(dims$K, Theta$A_p[,n] + gamma * 1, Theta$B_p[,n] + gamma * Theta$P[,n])
}

#' Sample Prior Parameter: Lambda for En
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_Lambda_En <- function(n, Theta, dims, gamma) {
    rgamma(dims$G, Theta$A_e[n,] + gamma * 1, Theta$B_e[n,] + gamma * Theta$E[n,])
}


#### Gamma Prior

#' Sample Prior Parameter: Beta for Pn
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length K
#' @noRd
sample_Beta_Pn <- function(n, Theta, dims, gamma) {
    rgamma(
        dims$K,
        Theta$A_p[,n] + gamma*Theta$Alpha_p[,n],
        Theta$B_p[,n] + gamma*Theta$P[,n]
    )
}

#' Sample Prior Parameter: Beta for En
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_Beta_En <- function(n, Theta, dims, gamma) {
    rgamma(
        dims$G,
        Theta$A_e[n,] + gamma*Theta$Alpha_e[n,],
        Theta$B_e[n,] + gamma*Theta$E[n,]
    )
}

#' Sample Prior Parameter: Alpha for Pkn
#'
#' @param k mutation type index
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Alpha_Pkn <- function(k, n, Theta, dims, gamma) {
    logpdf_prop <- function(x) {
        (Theta$C_p[k,n] - 1) * log(x) -
        Theta$D_p[k,n] * x +
        gamma * x * log(Theta$Beta_p[k,n]) +
        gamma * (x - 1) * log(Theta$P[k,n]) -
        gamma * lgamma(x)
    }

    armspp::arms(
        n_samples = 1,
        log_pdf = logpdf_prop,
        lower = 1e-3,
        upper = 10000
    )
}

#' Sample Prior Parameter: Alpha for Eng
#'
#' @param n signature index
#' @param g tumor genome index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Alpha_Eng <- function(n, g, Theta, dims, gamma) {
    logpdf_prop <- function(x) {
        (Theta$C_e[n,g] - 1) * log(x) -
        Theta$D_e[n,g] * x +
        gamma * x * log(Theta$Beta_e[n,g]) +
        gamma * (x - 1) * log(Theta$E[n,g]) -
        gamma * lgamma(x)
    }

    armspp::arms(
        n_samples = 1,
        log_pdf = logpdf_prop,
        lower = 1e-3,
        upper = 10000
    )
}

