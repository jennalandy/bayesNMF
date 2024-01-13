---
editor_options: 
  markdown: 
    wrap: 72
---

# bayesNMF: an R package for Bayesian NMF

## Notes

Hello! This package is a **work in progress**. I will be adding more models, but so far, I have 

- Poisson - Gamma: $M \sim Poisson(PE)$, $P$ and $E$ follow gamma priors
- Normal - Exponential: $M_k \sim N((PE)_k, \sigma^2_k I)$, $P$ and $E$ follow Exponential priors, $\sigma^2_k$ follows an Inverse-Gamma prior 
- Normal - Truncated Normal: $M_k \sim N((PE)_k, \sigma^2_k I)$, $P$ and $E$ follow Truncated-Normal priors, $\sigma^2_k$ follows an Inverse-Gamma prior

I also plan to add markdown files with full model specifications and derivations of the Gibbs updates used in these implementations.

While language and simulation examples are in the context of mutational signatures analysis, this package can be used for any application of NMF.

## Quick Start

### Installation

```{r}
library(devtools)
devtools::install_github("jennalandy/bayesNMF")
library(bayesNMF)
```

### Use

The functions for each model are:

-   `nmf_poisson_gamma`
-   `nmf_normal_exponential`
-   `nmf_normal_truncnormal`

All functions have the same structure, so I will just use
`nmf_normal_truncnormal` for the tutorial.

Running the Normal-Truncated Normal model requires a mutational catalog matrix $M$ and the number of signatures or latent factors $N$. We also recommend setting the file name to something related to your analysis. The default file name will be "nmf_normal_truncnormal".

```{r}
res <- nmf_normal_truncnormal(M, N = 5, file = "my_run")
```

Three files will be created and updated every 100 iterations by default (can be controlled with the `logevery` parameter) 

- `my_run.log` will note the current iteration the method is on, which is useful to estimate the total run time if using a large dataset or a lot of iterations. 
- `my_run.save` holds the current results in a R data file, which can be useful if your run is cut short or you need to access final results at a later time. 
- `my_run.pdf` has three plots: RMSE, KL Divergence, and log likelihood over iterations. Note that log likelihood is specific to the likelihood model, so values from the Poisson-Gamma model are not comparable to those from either of the Normal models.

The maximum a-posteriori (MAP) estimates for $P$ and $E$ are stored in `res$MAP$P` and `res$MAP$E`. The full Gibbs sampler chains are stored in `res$logs`. The reconstruction errors and log likelihood for each iteration are stored in `res$metrics`.

### Compare to True Signatures

We also include commands to compare estimated signatures to the true
signatures matrix to evaluate simulation studies. This could also be a set of signatures from literature that we wish to use as a baseline.

The following functions provide a similarity matrix between the true and estimated $P$ matrices, as well as a heatmap to visualize this.

```{r}
sim_mat <- pairwise_sim(res$MAP$P, true_P)
heatmap <- get_heatmap(est_P = res$MAP$P, true_P = true_P)
```

You can also do this all in one call. If `true_P` is provided to `nmf_normal_truncnormal`, then the similarity matrix and heatmap are stored in `res$sim_mat` and `res$heatmap`, respectively.

```{r}
res <- nmf_normal_truncnormal(M, N = 5, file = "my_run", true_P = P)
```

## Simulated Example

### 1. Simulate Data

Here we simulate a mutational catalog based on a subset of COSMIC signatures, available at the link below. Here, we use 5 of the COSMIC signatures.

```{r}
set.seed(123)

sigs = c(2,3,4,5,6)
N = length(sigs)
P <- read.csv(
  "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3.1_SBS_GRCh37.txt",
  sep = '\t'
)
P <- as.matrix(P[,sigs])
K = nrow(P)
```

We then simulate exposure values, for example here from an exponential distribution.

```{r}
G = 20
E <- matrix(rexp(N*G, 0.001), nrow = N, ncol = G)
```

Assuming a Poisson data generating function, we can finally simulate a mutational catalog $M$.

```{r}
M <- matrix(nrow = K, ncol = G)
for (k in 1:K) {
    M[k,] <- rpois(G, P[k,]%*%E)
}
```

### 2. Running Methods

Now we can run the Normal-Exponential model of Bayesain NMF.

```{r}
res <- nmf_normal_truncnormal(
  M, N = 5,
  file = "my_run",
  true_P = P
)
```

### 3. Evaluating Methods

Metrics are always recorded on each iteration. Here we look at the metrics on the last iteration.

```{r}
data.frame(res$metrics)[res$niters,]
```

```        
         loglik     RMSE      KL
2000 -272791295 5.701652 794.889
```

The similarity matrix and corresponding heatmap are only provided if
`true_P` is used.

```{r}
res$sim_mat
```

```         
                true1     true2     true3      true4     true5
estimated1 0.02639892 0.9993882 0.1246383 0.06637452 0.2701752
estimated2 0.95782867 0.2929271 0.1467035 0.06537255 0.3256951
estimated3 0.12814808 0.7505440 0.5363341 0.24709673 0.8122917
estimated4 0.04909267 0.0153710 0.6617648 0.99525126 0.4061263
estimated5 0.12320290 0.2215279 0.9806191 0.65717947 0.8239119
```

```{r}
res$heatmap
```

![](images/example_heatmap.png){width="491"}

Notice that signatures are reordered on the heatmap, assigned using the
Hungarian algorithm. To see this assignment on the similarity matrix,
you can use the `assign_signatures` function.

```{r}
assigned_sim_mat = assign_signatures(res$sim_mat)
assigned_sim_mat
```

```
               true2      true1     true5      true4     true3
estimated1 0.9993882 0.02639892 0.2701752 0.06637452 0.1246383
estimated2 0.2929271 0.95782867 0.3256951 0.06537255 0.1467035
estimated3 0.7505440 0.12814808 0.8122917 0.24709673 0.5363341
estimated4 0.0153710 0.04909267 0.4061263 0.99525126 0.6617648
estimated5 0.2215279 0.12320290 0.8239119 0.65717947 0.9806191
```
