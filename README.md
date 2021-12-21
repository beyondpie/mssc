
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
## install cmdstanr firslty.
## directly use install.package("cmdstanr") may face error.
install.packages("cmdstanr",
                 repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
## then install mssc
if(!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("beyondpie/mssc")
## mssc => cmdstanr => cmdstan  (=>: depends on)
## so we also need to install cmdstan, which is not an R packge, but cmdline Stan tool.
## Prefer to use cmdstanr to install it since cmdstanr will then remember where the cmdstan is,
## otherwise you have to install cmdstan independently and tell cmdstanr where it is by call
## `cmdstanr::set_cmdstan_path()`.
cmdstanr::install_cmdstan()
```

## Test

The test script can be download here:
<https://github.com/beyondpie/mssc/blob/main/inst/rscript/test.R>. You
can also use the R command below to find the script.

``` r
## test.R script can be got from the package.
test_script = system.file("rscript", "test.R", package = "mssc", mustWork = TRUE)
```

Then you can run this script to test if you can run mssc without errors.

## Rank differentially expressed genes

Currently, we use the function `evalDeltaMean` to get the differential
mean under two conditions. The absolute value could reflect the effect
size.

-   When we use variational inference, the function returns the average
    differential mean from the samples of conditional gene expression
    levels from the approximated joint posterior distribution. For each
    gene, we have one such value.

-   When we use optimization, the function returns the directly
    differential value from the estimated conditional gene expression
    levels.

Currently, we have no the corresponding p-value like statistics to
measure the differential expression.

## Generalized linear model

A simple way to model the batch effect for different cells is to treat
the batch as one covariant in a generalized linear model, which we have
implemented as a control for mssc.
