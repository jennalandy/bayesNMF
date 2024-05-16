set.seed(123)
N = 5
G = 60

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

sigs = 2:(2+N - 1) # c(2,3,4,5,10)
P <- P[,sigs]
K = nrow(P)

E <- matrix(rexp(N*G, 0.01), nrow = N, ncol = G)

M <- matrix(nrow = K, ncol = G)
for (k in 1:K) {
    M[k,] <- rpois(G, P[k,]%*%E)
}

true_P <- P
rm(list = c('P','E','N','G','sigs'))
