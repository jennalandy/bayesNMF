set.seed(123)
dims = list(
    K = 96,
    G = 60,
    N = 5
)

a = 10
b = 5
Beta = matrix(0.1, nrow = dims$N, ncol = dims$G)
r = rep(3, dims$K)
theta = rep(10, dims$K)


P <- read.csv(
    "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3.1_SBS_GRCh37.txt",
    sep = '\t'
)
P <- as.matrix(P[,2:ncol(P)])
sim <- pairwise_sim(P, P)
row_maxes <- sapply(1:nrow(sim), function(i) {
    sum(sim[i,-i] > 0.7)
})

P <- P[,-which(row_maxes > 0)]

sigs = 2:(2+dims$N - 1) # c(2,3,4,5,10)
P <- P[,sigs]


a = 10
b = 5
Beta = matrix(0.1, nrow = dims$N, ncol = dims$G)
r = rep(3, dims$K)
theta = rep(10, dims$K)

E <- matrix(nrow = dims$N, ncol = dims$G)
for (n in 1:dims$N) {
    for (g in 1:dims$G) {
        E[n,g] <- rexp(1, Beta[n,g])
    }
}

sigmasq <- vector(length = dims$K)
for (k in 1:dims$K) {
    sigmasq[k] <- invgamma::rinvgamma(1, r[k], theta[k])
}

M <- matrix(nrow = dims$K, ncol = dims$G)
for (g in 1:dims$G) {
    M[,g] <- rnorm(dims$K, mean = P %*% E[,g], sd = sqrt(sigmasq))
}
M <- round(M)

M[M<=0] <- 1
true_P <- P

rm(list = c('P','E','dims','theta','r','a','b','Beta'))
