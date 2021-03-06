% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mssc2.R
\name{new_MSSC2}
\alias{new_MSSC2}
\title{Create MSSC2 model object}
\usage{
new_MSSC2(modelpath = NULL, glmodelpath = NULL, seed = 1L)
}
\arguments{
\item{modelpath}{string, mssc2 stan script path
use system.file("stan", "mssc2.stan", package = "mssc", mustWork=TRUE) to get the script.}

\item{glmodelpath}{string glmodelpath, glm stan script path
use system.file("stan", "glm.stan", package = "mssc", mustWork=TRUE) to get the script.}

\item{seed}{int, default is 1.}
}
\value{
MSSC2 object
}
\description{
Create MSSC2 model object
}
