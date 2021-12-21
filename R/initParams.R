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
  if (sum(y) == 0) {
    stop("numeric vector y should have at least one element larger than zero.")
  }
  if (any(s <= 0)) {
    stop("normalizing constant s should be all larger than zero.")
  }
  mu <- log(stats::mean(y)) - log(stats::median(s))
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
  v <- stats::var(y)
  m <- mean(y)
  ## dispersion
  r <- ifelse(v > m, m^2 / (v - m), min_r)
  invisible(list(mu = mu, r = r))
}

#' Fit mu, mucond, muind for Negative Binomial with scaling.
#' No more Stan script based.
#' NOTE previous name: fit_mgsnb
#' @param y numeric vector, should be not all zeros.
#' @param s numeric vector, normalizing constant, should no zeros.
#' @param cond vector of integer, start from 1
#' @param ind vector of integer, start from 1
#' @param default_mu double, default is 0.0.
#' @param default_r double, default is 20.
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
  if (sum(y) == 0) {
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
    if (sum(yy) == 0) {
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

#' NOTE previous name: est_varofmu
#' @param mu numeric vector, mean expression for genes
#' @param min_var scalar, lower bound of estimated variance of umu
#' @return numeric vector with four elements:
#' mean, variance, alpha, beta (param of inv-gamma distribution for var of mu)
#' @export
fitGeneGlobalMeanPriorParams <- function(mu, min_var = 4.0) {
  m <- stats::median(mu)
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
#' NOTE previous name: est_varofr
#' @param r numeric vector, dispersions for genes
#' @param min_var scalar, lower bound of estimated var of r.
#' @return numeric vector with four elements:
#' mean, variance, alpha, beta (param of inv-gamma distribution for var of logr)
#' @export
fitDispersionParams <- function(r, min_var = 4.0) {
  invisible(fitGeneGlobalMeanPriorParams(mu = log(r), min_var = min_var))
}

#' Estimate parameters of inv-gamma dist as the prior
#' for genewise conditional expression.
#' NOTE previous name: est_varofcond
#' @param mucond matrix, ngene by ncond, gene expression level under conditions.
#' @param min_varofcond scalar, lower bound for the estimated variance.
#' @return matrix, ngene by 3 params,
#' including variance of the gene expression under that condition,
#' alpha, beta params for the prior of the variance (inv-gamma distribution).
#' We treat initial means of the gene expressions as zeros.
#' @export
fitGenewiseCondPriorParams <- function(mucond,
                                       min_varofcond = 0.25) {
  ncond <- ncol(mucond)
  ngene <- nrow(mucond)
  ## each condiiton has its own variance.
  ## which follows a inv-gamma prior
  t_d <- vapply(
    1:ncond, function(i) {
      t <- max(abs(mucond[, i]))
      ## set a variance not that small
      v <- max(t^2, min_varofcond)
      ## assume the mean of mucond is around 0.0
      ## then use posterior of inv-gamma to set the hyper priors
      alpha <- 1.0 + ngene / 2
      beta <- 1.0 + sum(mucond[, i]^2) / 2
      invisible(c(v, alpha, beta))
    },
    FUN.VALUE = rep(1.0, 3)
  )
  return(invisible(t(t_d)))
}

#' NOTE previous name: est_varofind
#' @param muind matrix, ngene by nind
#' @param min_varofind scalar, lower bound of estimated variance of individuals.
#' @param min_tau2 scalar, lower bound of estimated tau2.
#' @return list of two elements
#' - est_varofind: matrix, nind by 4
#'   - mean, variance, alpha and beta (inv-gamma dist as the prior of the variance)
#' - est_tau2: vector of 3 elements
#'   - tau2, tau2_alpha, tau2_beta (assume muinds follow a N(0.0, tau), taus is std.)
#' @export
fitGenewiseBatchPriorParams <- function(muind,
                                        min_varofind = 0.25,
                                        min_tau2 = 0.25) {
  nind <- ncol(muind)
  ngene <- nrow(muind)
  t_d <- vapply(
    1:nind, function(i) {
      m <- stats::median(muind[, i])
      v <- max(sum((muind[, i] - m)^2) / ngene, min_varofind)
      alpha <- 1.0 + ngene / 2
      beta <- 1.0 + sum((muind[, i] - m)^2) / 2
      invisible(c(m, v, alpha, beta))
    },
    FUN.VALUE = rep(1.0, 4)
  )
  ## shape: nind by 4
  r <- t(t_d)
  ## assume muinds follow a N(0.0, tau) (tau is sd)
  tau2 <- max(max(abs(r[, 1]))^2, min_tau2)
  ## assume tau2 has a inv-gamma prior
  ## use posterior to set up the hp.
  tau2_alpha <- 1.0 + nind / 2
  tau2_beta <- 1.0 + sum(muind[, 1]^2) / 2
  return(invisible(
    list(
      est_varofind = r,
      est_tau2 = c(tau2, tau2_alpha, tau2_beta)
    )
  ))
}

#' A summary fit function for Genewise Negative Binomial Model
#' NOTE previous name: fit_mgsnb
#' @param cnt matrix, count matrix size of ngene by ncell
#' @param s numeric vector, length of ncell, normalzing/scaling factor for cells
#' @param cond integer vector, length of ncell, which cond the cell belongs to,
#' start from 1.
#' @param ind integer vector, length of ncell, which individuals cells belong to,
#' start from 1.
#' @param default_mu double, default is 0.0
#' @param default_r double, default is 20.
#' @param min_var double, default is 4.0
#' @param min_varofcond double, default is 0.25
#' @param min_varofind double, default is 0.25
#' @param min_tau2 double, default is 0.25
#' @return list of four elements
#' - mgsnb: matrix, ngene by 2 + ncond + nind
#' - mu: output of fitGeneGlobalMeanPriorParams
#' - logr: output of fitGeneGlobalMeanPriorParams
#' - cond: output of fitGenewiseCondPriorParams
#' - ind: output of fitGenewiseBatchPriorParams
#' @export
fitGenewiseNBModel <- function(cnt, s, cond, ind,
                               default_mu = 0.0, default_r = 20,
                               min_var = 4.0, min_varofcond = 0.25,
                               min_varofind = 0.25, min_tau2 = 0.25) {
  if (any(s) == 0) {
    stop("Normalizing factor s should have no zeros.")
  }
  ncond <- max(cond)
  nind <- max(ind)
  ngene <- nrow(cnt)
  t_init_mgsnb <- vapply(1:ngene, function(i) {
    r <- initNBParamWithCondBatch(
      y = cnt[i, ], s = s,
      cond = cond, ind = ind,
      default_mu = default_mu, default_r = default_r
    )
    invisible(unlist(r))
  },
  FUN.VALUE = rep(0.0, 2 + ncond + nind)
  )
  ## shape: ngene by 2 + ncond + nind
  init_mgsnb <- t(t_init_mgsnb)
  init_varofmu <- fitGeneGlobalMeanPriorParams(
    mu = init_mgsnb[, 1],
    min_var = min_var
  )
  init_varofr <- fitGeneGlobalMeanPriorParams(
    mu = log(init_mgsnb[, 2]),
    min_var = min_var
  )
  init_varofcond <- fitGenewiseCondPriorParams(
    mucond = init_mgsnb[, 3:(2 + ncond)],
    min_varofcond = min_varofcond
  )
  init_varofind <- fitGenewiseBatchPriorParams(
    muind = init_mgsnb[, (2 + ncond + 1):ncol(init_mgsnb)],
    min_varofind = min_varofind, min_tau2 = min_tau2
  )
  return(invisible(list(
    mgsnb = init_mgsnb,
    mu = init_varofmu,
    logr = init_varofr,
    cond = init_varofcond,
    ind = init_varofind
  )))
}
