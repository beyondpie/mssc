
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- ```{r, include = FALSE} -->
<!-- knitr::opts_chunk$set( -->
<!--   collapse = TRUE, -->
<!--   comment = "#>", -->
<!--   fig.path = "man/figures/README-", -->
<!--   out.width = "100%" -->
<!-- ) -->
<!-- ``` -->
<!-- # mssc -->
<!-- badges: start -->
<!-- badges: end -->

# MSSC: Multiple-subject differential analysis for single-cell RNA sequencing data

The goal of mssc is to provide a tool for differential analysis for
single-cell RNA sequencing data by considering the batch effects induced
when multiple subjects under different conditions, which is ignored by
pseudo-bulk analysis and cannot be handled well by simple generalized
linear models.

## Installation

``` r
if(!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("beyondpie/mssc")
```
