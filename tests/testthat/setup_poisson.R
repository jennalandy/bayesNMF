set.seed(123)
N = 5
G = 20

sigs = c(2,3,4,5,6)
N = length(sigs)
P <- read.csv(
    "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3.1_SBS_GRCh37.txt",
    sep = '\t'
)
P <- as.matrix(P[,sigs])
K = nrow(P)

E <- matrix(rexp(N*G, 0.001), nrow = N, ncol = G)

M <- matrix(nrow = K, ncol = G)
for (k in 1:K) {
    M[k,] <- rpois(G, P[k,]%*%E)
}

true_P <- P
rm(list = c('P','E','N','G','sigs'))
