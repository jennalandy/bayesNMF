# bayesNMF: Gibbs sampling for Bayesian Nonnegative Matrix Factorization with automatically learned rank and optional Metropolis-Hastings updates for computational efficiency.

Introduced in the paper ["bayesNMF: Fast Bayesian Poisson NMF with Automatically Learned Rank Applied to Mutational Signatures"](https://arxiv.org/abs/2502.18674).

For code to reproduce all results in the paper, see the [bayesNMF_PAPER respository](https://github.com/jennalandy/bayesNMF_PAPER).

## Quick Start

### Installation

```{r}
library(devtools)
devtools::install_github("jennalandy/bayesNMF")
library(bayesNMF)
```

Please see our [tutorial](vignettes/bayesNMF_tutorial.pdf)

Model details and advanced options can be found [here](vignettes/advanced.pdf).
