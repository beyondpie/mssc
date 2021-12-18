#' Estimate the normalized mean of a numeric vector on log-level,
#' which is similar with log(mean(y/s)) (or better),
#'
#' NOTE previsous name: init_snb_log_mean
#' @param y numeric vector, should be not all zeros.
#' @param s numeric vector, normalizing constant, should no zeros.
#' @return numeric scalar
#' @export
normLogMean <- function(y, s) {
  if (any(y) < 0) {
    stop("numeric vector y should be no less than zero.")
  }
  if(sum(y) == 0){
    stop("numeric vector y should have at least one element larger than zero.")
  }
  if(any(s <= 0) ) {
    stop("normalizing constant s should be all larger than zero.")
  }
  mu <- log(mean(y)) - log(median(s))
  invisible(mu)
}

#' Get the mu (log-level mean) and r (dispersion) in Negative Binomial with scaling.
#' dispersion is defined as the same as in either in R or Stan.
#'
#' NOTE:
#' From observation, mean can be estimate well, but if it's below 1,
#' then r tends to be estimated too small.
#' @param y numeric vector, should be not all zeros.
#' @param s numeric vector, normalizing constant, should no zeros.
#' @param min_r numeric scalar, dispersion lower threshold.
#' @return list of two-named scalar, mu and r
#' @export
initNBParam <- function(y, s, min_r) {
  mu <- normLogMean(y = y, s = s)
  ## v equals to m + m^2/r
  v <- var(y)
  m <- mean(y)
  ## dispersion
  r <- ifelse(v > m, m^2 / (v - m), min_r)
  invisible(list(mu = mu, r = r))
}

#' Fit mu, mucond, muind for Negative Binomial with scaling.
#' No more Stan script based.
#' @param cond vector of integer, start from 1
#' @param ind vector of integer, start from 1
#' @return list of four-named element, mu, r, mucond, muind
#' mu and r are scalars.
#' mucond is a vector length of number of conditions, usually 2.
#' muind is a vector length of number of individuals.
#' @export
initNBParamWithCondBatch <- function(y, s, cond, ind,
                                     default_mu = 0.0,
                                     default_r = 20) {
  nind <- max(ind)
  ncond <- max(cond)
  result <- list(
    mu = default_mu,
    r = default_r,
    mucond = rep(0.0, ncond),
    muind = rep(0.0, nind)
  )
  if (any(y) < 0) {
    stop("numeric vector y should be no less than zero.")
  }
  if(sum(y) == 0){
    ## if y are all zeros, then we use the default ones.
    ## TODO: maybe we can set muind as the mean of muind from our model,
    ## which can be directly estimated by a pure non-mucond effect model.
    ## for example, we can use random effect model?
    return(invisible(result))
  }

  init_mur <- initNBParam(y, s)
  result$mu <- init_mur$mu
  result$r <- init_mur$r
  
  result$mucond <- vapply(1:ncond, function(i) {
      ss <- s[cond == i]
      yy <- y[cond == i]
      if(sum(yy) == 0) {
        return(0.0)
      }
      t <- normLogMean(y = yy, s = ss)
      init_mucond <- t - result$mu
      ## Since fitting is hard when fix r, and limited data
      ## we directly use the init_mucond
      invisible(init_mucond)
    },
    FUN.VALUE = 0.0
  )
  
  result$muind <- vapply(1:nind, function(i) {
      yy <- y[ind == i]
      if (sum(yy) == 0) {
        return(0.0)
      }
      t <- normLogMean(y = yy, s = s[ind == i])
      return(invisible(t - result$mu - result$mucond[cond[ind == i][1]]))
    },
    FUN.VALUE = 0.0
  )
  invisible(result)
}

#' @param mu numeric vector, mean expression for genes
#' @param min_var scalar, lower bound of estimated variance of mu
#' @return numeric vector with four elements:
#' mean, variance, alpha, beta (param of inv-gamma distribution for var of mu)
#' @export
fitGeneGlobalMeanPriorParams <- function(mu, min_var = 4.0) {
  m <- median(mu)
  ngene <- length(mu)
  v <- sum((mu - m)^2) / ngene
  v <- max(v, min_var)
  ## varofmu prior follows a inv-gamma dist
  ## est hy based on posterior
  alpha <- 1.0 + ngene / 2
  beta <- 1.0 + sum((mu - m)^2) / 2
  invisible(c(m, v, alpha, beta))
}

#' Get parameters for the prior of dispersion.
#' We assume the dispersion follows a log-Normal distribution.
#' @param r numeric vector, dispersions for genes
#' @param min_var scalar, lower bound of estimated var of r.
#' @return numeric vector with four elements:
#' mean, variance, alpha, beta (param of inv-gamma distribution for var of logr)
#' @export
fitDispersionParams <- function(r, min_var = 4.0){
  invisible(fitGeneGlobalMeanPriorParams(mu = log(r), min_var = min_var))
}






