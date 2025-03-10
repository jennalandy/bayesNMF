% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_results}
\alias{plot_results}
\title{Plot signature contributions}
\usage{
plot_results(res_list, ref_matrix = "cosmic", title = "", return_df = TRUE)
}
\arguments{
\item{res_list}{Named list, containing one or more bayesNMF objects. Names will become identifiers along the top of the plot.}

\item{ref_matrix}{matrix, "cosmic", or NULL, reference signatures to align to}

\item{title}{string, title of the produced plot}

\item{return_df}{boolean, whether to return summary data frame}
}
\value{
ggplot object, or list holding ggplot object and data frame if return_df = TRUE
}
\description{
Plot signature contributions
}
\examples{
res <- readRDS("examples/plot_example.rds")
data <- readRDS("examples/3_64_1_cosmic.rds")
plot_results(list("Example" = res), title = "Results of a single run")
plot_results(list("Example" = res), reference = data$P, title = "Results of a single run with custom reference")
}
