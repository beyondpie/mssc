#' Create GenewiseNBFit S3 object
#' It's used for MSSC model initialization.
#' 
#' @return GenewiseNBFit object
#' @example
#' new_GenewiseNBFit()
new_GenewiseNBFit <- function() {
  structure(list(
    NBModel = NULL,
    r_alpha = 0.05,
    r_beta = 0.05,
    mu = 0.0,
    r = 20,
    big_r = 500,
    min_varofmu = 4.0,
    min_varofcond = 0.25,
    min_varofind = 0.25,
    min_tau2 = 0.25
  ), class = "GenewiseNBFit")
}

str.GenewiseNBFit <- function(x) {
  print(x)
}

print.GenewiseNBFit <- function(x) {
  message(paste("Negative Binomial Model:", x$NBModel))
  message(paste("NB disperson's prior of Gamma distribution:",
                paste("alpha:", x$r_alpha), paste("beta:", x$r_beta)))
  message(paste("Negative Binomial mean init value:", x$mu))
  message(paste("Negative Binomial disperson init value:", x$r))
}


