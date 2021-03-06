% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_bw_cluster_summary.R
\name{plot_hclust_heatmap.assemble}
\alias{plot_hclust_heatmap.assemble}
\title{Title}
\usage{
plot_hclust_heatmap.assemble(plot_list)
}
\arguments{
\item{plot_list}{}
}
\value{

}
\description{
Title
}
\examples{
library(ssvRecipes)
mat = matrix(runif(100), ncol = 10)
plot_list = plot_hclust_heatmap(mat)
plot_list$heatmap = plot_list$heatmap + scale_fill_viridis_c(limits = c(-1, 1))
plot_hclust_heatmap.assemble(plot_list)
}
