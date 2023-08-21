# bayesNMF: an R package with various models for Bayesian NMF

While data and simulation examples are in the context of mutational signatures analysis, this package can be used for any application of NMF.

## Quick Start

### Installation

```{r}
library(devtools)
devtools::install_github("jennalandy/bayesNMF")
library(bayesNMF)
```

### Use

Running the Normal-Exponential model requires a mutational catalog matrix $M$ and the number of signatures or latent factors $N$. We also recommend setting the file name to something related to your analysis. The default file name will be "nmf_normal_exponential".

```{r}
res <- nmf_normal_exponential(M, N = 5, file = "my_run")
```

Three files will be created and updated every 100 iterations by default (can be controlled with the `logevery` parameter)
- `my_run.log` will note the current iteration the method is on, which is useful to estimate the total run time if using a large dataset or a lot of iterations. 
- `my_run.save` holds the current results in a R data file, which can be useful if your run is cut short or you need to access final results at a later time.
- `my_run.pdf` has two plots: RMSE over iterations and loglikelihood over iterations.

The final estimates for $P$ and $E$ are stored in `res$P.mean` and `res$E.mean`. The full Gibbs sampler chains are stored in `res$P.chain` and `res$E.chain`. The reconstruction error and loglikelihood under the method's respective model are stored in `res$final_metrics`. The chains of these values over the iterations are stored in `res$RMSE.chain` and `res$loglik.chain`. 

### Compare to True Signatures

We also include commands to compare estimated signatures to the true signatures matrix to evaluate simulation studies. This could also be a set of signatures from literature that we wish to use as a baseline.

The following functions provide a similarity matrix between the true and estimated $P$ matrices, as well as a heatmap to visualize this.

```{r}
sim_mat <- get_sim_mat(est_P = res$P.mean, true_P = true_P)
heatmap <- get_heatmap(est_P = res$P.mean, true_P = true_P)
```

You can also do this all in one call. If `true_P` is provided to `nmf_normal_exponential`, then the similarity matrix and heatmap are stored in `res$sim_mat` and `res$heatmap`, respectively.

```{r}
res <- nmf_normal_exponential(M, N = 5, file = "my_run", true_P = P)
```

## Simulated Example

### 1. Simulate Data

Here we simulate a mutational catalog based on a subset of COSMIC signatures, available at the link below. Here, we sample 5 signatures.

```{r}
set.seed(123)

sigs = c(2,3,4,5,6)
N = length(sigs)
P <- read.csv(
  "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3.1_SBS_GRCh37.txt",
  sep = '\t'
)
P <- as.matrix(P[,sigs])
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
res <- nmf_normal_exponential(
  M, N = 5,
  file = "my_run",
  true_P = P
)
```

### 3. Evaluating Methods

Final metrics are always provided.

```{r}
res$final_metrics
```
```
$loglik
[1] -212892210

$RMSE
[1] 5.68804
```

The similarity matrix and corresponding heatmap are only provided if `true_Theta` is used.

```{r}
res$sim_mat
```
```
           SBS1       SBS2      SBS3       SBS4      SBS5
[1,] 0.02825740 0.99935370 0.1253304 0.06603725 0.2698738
[2,] 0.96033001 0.28413701 0.1467032 0.06751348 0.3236676
[3,] 0.12572760 0.76450866 0.5361641 0.27063974 0.8011454
[4,] 0.03723758 0.01831284 0.6612858 0.99558598 0.4039419
[5,] 0.12688657 0.21944638 0.9809169 0.64927609 0.8246318
```

```{r}
res$heatmap
```
![image](https://github.com/jennalandy/bayesNMF/assets/35237833/793b8562-6185-4b0b-ab76-dcf12b1ced0a)
