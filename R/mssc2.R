#' Create MSSC2 model object
#'
#' @param modelpath string, mssc2 stan script path
#' use system.file("stan", "mssc2.stan", package = "mssc", mustWork=TRUE) to get the script.
#' @param glmodelpath string glmodelpath, glm stan script path
#' use system.file("stan", "glm.stan", package = "mssc", mustWork=TRUE) to get the script.
#' @param seed int, default is 1.
#' @return MSSC2 object
#' @export
new_MSSC2 <- function(modelpath = NULL, glmodelpath = NULL, seed = 1L) {
  r <- structure(list(
    seed = seed,
    modelpath = modelpath,
    model = NULL,
    init_mgsnb = NULL,
    ## predefined params, not model-related hyper params.
    methodparams = list(
      default_mu = 0.0, default_r = 20,
      min_var = 4.0, min_varofcond = 0.25,
      min_varofind = 0.25, min_tau2 = 0.25
    ),
    ncond = NULL,
    nind = NULL,
    ngene = NULL,
    ncell = NULL,
    ## model hyper parameters
    modelhp = NULL,
    ## model parameters' initial values
    modelip = NULL,
    vi = NULL,
    opt = NULL,
    stanparams = list(
      num_iter = 20000,
      vi_refresh = 2000,
      algorithm = "meanfield",
      eval_elbo = 100,
      output_samples = 1000,
      tol_rel_obj = 1e-4,
      adapt_iter = 200,
      adapt_engaged = FALSE,
      eta = 0.1,
      opt_method = "lbfgs",
      init_alpha = 0.001,
      opt_refresh = 10,
      opt_max_iter = 5000,
      tol_obj = NULL,
      tol_grad = NULL,
      tol_rel_grad = NULL,
      tol_param = NULL,
      history_size = 10
    ),
    allparams = c(
      "centerofmu", "varofmu", "mu", "centerofr",
      "varofr", "r", "varofcond", "mucond",
      "tau2", "centerofind", "varofind", "muind"
    ),
    glmodelpath = glmodelpath,
    glmodel = NULL,
    glmip = NULL,
    glmvi = NULL,
    glmopt = NULL,
    allparamsglm = c(
      "mu", "nb_r", "mucond", "muind"
    )
  ), class = "MSSC2")

  ## only compile model. glm need to be explicitly compiled when in use.
  r$model <- compileStanModel(
    model_path = r$modelpath,
    use_thread = NULL,
    use_mpi = NULL, use_opencl = NULL
  )
  invisible(r)
}

#' Initialize MSSC2 Genewise parameters, which is further used for
#' initialize stan model parameters.
#'
#' Call fitGenewiseNBModel to initialize the parameters.
#'
#' @param mssc2 MSSC2 object
#' @param cnt matrix, count matrix size of ngene by ncell
#' @param s numeric vector, length of ncell, normalzing/scaling factor for cells
#' @param cond integer vector, length of ncell, which cond the cell belongs to,
#' start from 1.
#' @param ind integer vector, length of ncell, which individuals cells belong to,
#' start from 1.
#' @return MSSC2 object
#' @export
initMSSC2GenewiseParams <- function(mssc2, cnt, s, cond, ind) {
  if (class(mssc2) != "MSSC2") {
    stop("Input mssc2 is not an objet of MSSC2 S3 class.")
  }
  p <- mssc2$methodparams
  init_mgsnb <- fitGenewiseNBModel(
    cnt = cnt, s = s, cond = cond, ind = ind,
    default_mu = p$defualt_mu, default_r = p$default_r,
    min_var = p$min_var, min_varofcond = p$min_varofcond,
    min_varofind = p$min_varofind, min_tau2 = p$min_tau2
  )
  m <- mssc2
  m$init_mgsnb <- init_mgsnb
  m$ngene <- nrow(cnt)
  m$ncell <- ncol(cnt)
  m$ncond <- max(cond)
  m$nind <- max(ind)
  invisible(m)
}

#' Set hyper parameters for the stan model of mssc2.
#'
#' @param mssc2 MSSC2 object
#' @return MSSC2 object
#' @export
setStanHyperparams <- function(mssc2) {
  if (class(mssc2) != "MSSC2") {
    stop("Input mssc2 is not an objet of MSSC2 S3 class.")
  }
  if (is.null(mssc2$init_mgsnb)) {
    stop("mssc2 init_mgsnb is NULL. Use fitGenewiseNBModel to set it.")
  }
  init_mgsnb <- mssc2$init_mgsnb
  hp <- list(
    hp_varofmu = init_mgsnb$mu[3:4],
    hp_varofr = init_mgsnb$logr[3:4],
    hp_varofcond = init_mgsnb$cond[, 2:3],
    hp_varofind = init_mgsnb$ind$est_varofind[, 3:4],
    hp_tau2 = init_mgsnb$ind$est_tau2[2:3]
  )
  ## create a new mssc2 object
  mssc2$modelhp <- hp
  invisible(mssc2)
}

#' Set initial values for the parameters in the stan model of mssc2.
#'
#' TODO: muind estimated by the genes from DESEQ2.
#' @param mssc2 MSSC2 object
#' @return MSSC2 object
#' @export
setStanInitialParams <- function(mssc2) {
  if (class(mssc2) != "MSSC2") {
    stop("Input mssc2 is not an objet of MSSC2 S3 class.")
  }
  if (is.null(mssc2$init_mgsnb)) {
    stop("mssc2 init_mgsnb is NULL. Use fitGenewiseNBModel to set it.")
  }
  init_mgsnb <- mssc2$init_mgsnb
  ## * set params initiolization
  centerofmu <- init_mgsnb$mu[1]
  varofmu <- init_mgsnb$mu[2]
  mu <- init_mgsnb$mgsnb[, 1]
  raw_mu <- (mu - centerofmu) / sqrt(varofmu)

  ### r in negative binomial
  r <- init_mgsnb$mgsnb[, 2]
  ### NOTE: log of r level
  #### Should change the name of r to logr
  centerofr <- init_mgsnb$logr[1]
  varofr <- init_mgsnb$logr[2]
  raw_r <- (log(r) - centerofr) / sqrt(varofr)

  ### ngene by ncond
  ncond <- mssc2$ncond
  nind <- mssc2$nind
  mucond <- init_mgsnb$mgsnb[, 3:(2 + ncond)]
  ### ncond by 1
  varofcond <- init_mgsnb$cond[, 1]
  raw_mucond <- mucond %*% diag(1 / sqrt(varofcond))

  ### scalar
  tau2 <- init_mgsnb$ind$est_tau2[1]
  ### nind by 1
  centerofind <- init_mgsnb$ind$est_varofind[, 1]
  ### nind by 1
  raw_centerofind <- centerofind / sqrt(tau2)
  ### nind by 1
  varofind <- init_mgsnb$ind$est_varofind[, 2]
  ### ngene by nind
  muind <- init_mgsnb$mgsnb[, (2 + ncond + 1):(2 + ncond + nind)]
  ngene <- mssc2$ngene
  raw_muind <- (muind - repVecRowise(centerofind, n = ngene)) %*%
    diag(1 / sqrt(varofind))
  mssc2$modelip <- list(
    centerofmu = centerofmu,
    varofmu = varofmu,
    raw_mu = raw_mu,
    centerofr = centerofr,
    varofr = varofr,
    raw_r = raw_r,
    varofcond = varofcond,
    raw_mucond = raw_mucond,
    tau2 = tau2,
    raw_centerofind = raw_centerofind,
    varofind = varofind,
    raw_muind = raw_muind
  )
  invisible(mssc2)
}

#' Convert input data to the format of Stan Input.
#'
#' @param cnt matrix, count matrix size of ngene by ncell
#' @param ind integer vector, length of ncell, which individuals cells belong to,
#' start from 1
#' @param cond integer vector, length of ncell, which conditions cells belong to,
#' start from 1.
#' @param s normalizing/scaling factor, length of ncell.
#' @param hp list of 4 elements.
#' @export
toStanInput <- function(cnt, ind, cond, s, hp) {
  invisible(c(list(
    ncell = ncol(cnt),
    nind = max(ind),
    ncond = max(cond),
    ngene = nrow(cnt),
    s = s, cond = cond,
    ind = ind, y = t(cnt)
  ), hp))
}

#' Run variational inference for mssc2.
#' @param mssc2 MSSC2 object
#' @param data list, generated by toStanInput
#' @return MSSC2 object
#' @export
runVI <- function(mssc2, data) {
  if (class(mssc2) != "MSSC2") {
    stop("Input mssc2 is not an objet of MSSC2 S3 class.")
  }

  if (is.null(mssc2$model)) {
    stop("mssc2 model has not been compiled.")
  }

  if (is.null(mssc2$modelip)) {
    stop("mssc2 hasn't initialize the parameters.")
  }

  sp <- mssc2$stanparams
  mssc2$vi <- mssc2$model$variational(
    data = data,
    init = list(mssc2$modelip),
    seed = mssc2$seed,
    refresh = sp$vi_refresh,
    iter = sp$num_iter,
    eval_elbo = sp$eval_elbo,
    adapt_engaged = sp$adapt_engaged,
    adapt_iter = sp$adapt_iter,
    algorithm = sp$algorithm,
    output_samples = sp$output_samples,
    tol_rel_obj = sp$tol_rel_obj,
    eta = sp$eta
  )
  return(mssc2)
}

#' Run MAP optimization for mssc2.
#' @param mssc2 MSSC2 object
#' @param data list, generated by toStanInput
#' @return MSSC2 object
#' @export
runMAP <- function(mssc2, data) {
  if (class(mssc2) != "MSSC2") {
    stop("Input mssc2 is not an objet of MSSC2 S3 class.")
  }

  if (is.null(mssc2$model)) {
    stop("mssc2 model has not been compiled.")
  }
  if (is.null(mssc2$modelip)) {
    stop("mssc2 hasn't initialize the parameters.")
  }
  sp <- mssc2$stanparams
  mssc2$opt <- mssc2$model$optimize(
    data = data,
    init = list(mssc2$modelip),
    seed = mssc2$seed,
    refresh = sp$opt_refresh,
    iter = sp$opt_max_iter,
    algorithm = sp$opt_method,
    init_alpha = sp$init_alpha,
    tol_obj = sp$tol_obj,
    tol_rel_obj = sp$tol_rel_obj,
    tol_grad = sp$tol_grad,
    tol_rel_grad = sp$tol_rel_grad,
    tol_param = sp$tol_param,
    history_size = sp$history_size
  )
  invisible(mssc2)
}

#' Get Pareto Smoothed Importance Sampling weights.
#' Weights are normalized and NOT in scale level (what do I mean?) by default.
#'
#' This is specific for variational inference.
#'
#' @param mssc2 MSSC2 object
#' @param takelog bool, default FALSE
#' @param donormalize bool, default TRUE
#' @return list of two elements:
#' - psis: output of loo:psis
#' - weight: output of loo::weights.importance_sampling
#' @export
PSIS <- function(mssc2, takelog = FALSE, donormalize = TRUE) {
  if (class(mssc2) != "MSSC2") {
    stop("Input mssc2 is not an objet of MSSC2 S3 class.")
  }
  if (is.null(mssc2$vi)) {
    stop("MSSC2 hasn't run variotional inference yet or may have troubles on that.")
  }
  ## TODO: check the meaninig below
  log_ratios <- mssc2$vi$lp() - mssc2$vi$lp_approx()
  utils::capture.output(suppressWarnings(r <- loo::psis(
    log_ratios = log_ratios,
    r_eff = NA
  )))
  ## call loo::weights.importance_sampling
  w <- stats::weights(r, log = takelog, normalize = donormalize)
  invisible(list(psis = r, weight = w))
}

#' Extract draws/samples from the MSSC2 model
#' @param mssc2 MSSC2 object
#' @param param string, name of the parameter we want to extract
#' @param genenms vector of string, length of ngene
#' @param method string, "opt" or "vi", default is "vi"
#' @return NA (when have error) or vector/matrix (depend on param)
#' - samples of the parameter for all the genes
#' @export
extractDraws <- function(mssc2, param, genenms = NULL, method = "vi") {
  if (class(mssc2) != "MSSC2") {
    stop("Input mssc2 is not an objet of MSSC2 S3 class.")
  }
  if (method == "vi") {
    fit <- mssc2$vi
  } else {
    fit <- mssc2$opt
  }
  if (is.null(fit)) {
    stop(paste(method, " is null. Might becuase it hasn't been run or have troubles."))
  }
  if (!(param %in% mssc2$allparams)) {
    stop(paste(
      param, "is not in the list of",
      paste(mssc2$allparams, collapse = ","), "."
    ))
  }
  if ((!is.null(genenms)) & (length(genenms) != mssc2$ngene)) {
    stop(paste(
      length(genenms), "is not equal to ngene",
      mssc2$ngene, "in the model."
    ))
  }
  tryCatch(
    {
      ## when the corresponding parameter is scalar.
      if (param %in% c("centerofmu", "varofmu", "centerofr", "varofr", "tau2")) {
        return(invisible(fit$draws(param)))
      }
      ## when param is vector length of ngene
      if (param %in% c("mu", "r", "nb_r")) {
        t <- fit$draws(getNameOfVectorFromCmdstanr(nm = param, n = mssc2$ngene))
        names(t) <- genenms
        return(invisible(t))
      }
      ## when param is vector length of ncond
      if (param %in% c("varofcond")) {
        return(invisible(fit$draws(getNameOfVectorFromCmdstanr(nm = param, n = mssc2$ncond))))
      }
      ## when param is vector length of nind
      if (param %in% c("varofind", "centerofind")) {
        return(invisible(fit$draws(getNameOfVectorFromCmdstanr(nm = param, n = mssc2$nind))))
      }
      ## when param is matrix size of ngene-by-ncond
      if (param %in% c("mucond")) {
        t <- fit$draws(getNameOfMatrixFromCmdstanr(nm = param, nr = mssc2$ngene, nc = mssc2$ncond))
        return(invisible(split_matrix_col(mat = t, secondim = mssc2$ngene, secondimNms = genenms)))
      }
      ## when param is matrix size of ngene-by-nind
      if (param %in% c("muind")) {
        t <- fit$draws(getNameOfMatrixFromCmdstanr(nm = param, nr = mssc2$ngene, nc = mssc2$nind))
        return(invisible(split_matrix_col(mat = t, secondim = mssc2$ngene, secondimNms = genenms)))
      }
    }, ## end of try block
    error = function(e) {
      warning(e)
      invisible(NA)
    }
  ) ## end of trycatch block
}

#' Extract all parameters' samples from MSSC2 Model.
#' @param mssc2 MSSC2 object
#' @param genenms NULL or vector of strings.
#' @param method string, "vi" or "opt", default is "vi"
#' @return list, samples for different parameters.
#' @export
extractDrawsAll <- function(mssc2, genenms = NULL, method = "vi") {
  if (class(mssc2) != "MSSC2") {
    stop("Input mssc2 is not an objet of MSSC2 S3 class.")
  }
  if (method == "vi") {
    fit <- mssc2$vi
  } else {
    fit <- mssc2$opt
  }
  if (is.null(fit)) {
    stop(paste(method, " is null. Might becuase it hasn't been run or have troubles."))
  }
  est_params <- lapply(mssc2$allparams, function(nm) {
    extractDraws(mssc2 = mssc2, genenms = nm, method = method)
  })
  names(est_params) <- mssc2$allparams
  invisible(est_params)
}

#' Evaluate delta mean between two conditions.
#' @param mucond 3-dim arrays, nsample by ngene by ncond
#' @param twoHotVec integer vector, like (1,-1) or (0, 0, -1, 1, 0),
#' the two conditions we want to compare, one is 1 and the other is -1.
#' @return numeric vector, length of ngene, delta mean for different genes.
#' @export
evalDeltaMean <- function(mucond, twoHotVec) {
  if (sum(twoHotVec) != 0) {
    stop("Set 1 and -1 for the two conditions we are interested.")
  }
  if (dim(mucond)[3] != length(twoHotVec)) {
    stop(paste(
      "Mucond has", dim(mucond)[3], "conditions;",
      "twoHotVec has", length(twoHotVec), "conditions."
    ))
  }
  n <- dim(mucond)[1]
  tmp <- t(vapply(1:n, function(i) {
    mucond[i, , ] %*% twoHotVec
  }, FUN.VALUE = rep(0.0, dim(mucond)[2])))
  r <- colMeans(tmp)
  if (!is.null(dimnames(mucond)[[2]])) {
    names(r) <- dimnames(mucond)[[2]]
  }
  invisible(r)
}

#' Evaluate t statistics and p-value between two conditions.
#' @param mucond 3-dim arrays, nsample by ngene by ncond
#' @param twoHotVec integer vector, like (1,-1) or (0, 0, -1, 1, 0),
#' the two conditions we want to compare, one is 1 and the other is -1.
#' @param alternative string, default "two.sided"
#' @param paired bool, default FALSE
#' @param var.equal bool, default FALSE
#' @return matrix, ngene by 2 (t statistics and pvalue in order)
#' NA will be used for t statistics when some troubles happend in t.test,
#' in this case, p-value is 1.0
#' @export
evalTstat <- function(mucond, twoHotVec,
                      alternative = "two.sided",
                      paired = TRUE,
                      var.equal = FALSE) {
  if (sum(twoHotVec) != 0) {
    stop("Set 1 and -1 for the two conditions we are interested.")
  }
  if (dim(mucond)[3] != length(twoHotVec)) {
    stop(paste(
      "Mucond has", dim(mucond)[3], "conditions;",
      "twoHotVec has", length(twoHotVec), "conditions."
    ))
  }
  group1 <- mucond[, , twoHotVec == 1]
  group2 <- mucond[, , twoHotVec == -1]
  ngene <- dim(mucond)[2]
  tstat <- vapply(1:ngene, function(i) {
    tryCatch(
      {
        s <- stats::t.test(
          x = group1[, i], y = group2[, i],
          alternative = alternative,
          paired = paired,
          var.equal = var.equal
        )
        invisible(c(s$statistic, s$p.value))
      },
      error = function(e) {
        warning(e)
        invisible(c(NA, 1.0))
      }
    )
  }, FUN.VALUE = c(0.0, 1.0))
  r <- t(tstat)
  colnames(r) <- c("tstat", "pvalue")
  if (!is.null(dimnames(mucond)[[2]])) {
    rownames(r) <- dimnames(mucond)[[2]]
  }
  invisible(r)
}
#' Calculate AUCs for stats in colume-wise.
#' @param stats matrix, ngene by n_measurements
#' @param c1 integer vector, index of genes for condition 1
#' @param c2 integer vector, index of genes for condition -1
#' @return numeric vector, within 0-1 region, AUCs for different measurements.
#' @export
calColAUC <- function(stats, c1, c2) {
  t <- stats
  if (!is.matrix(stats)) {
    t <- as.matrix(stats, ncol = 1)
  }
  groundtruth <- c(rep(TRUE, length(c1)), rep(FALSE, length(c2)))
  r <- caTools::colAUC(t[c(c1, c2), ], groundtruth)
  invisible(r)
}

#' Set initial values for the generalized linear model.
#' TODO: muind estimated by the genes from DESEQ2.
#' @param mssc2 MSSC2 object
#' @return MSSC2 object
#' @export
setGLMInitialParams <- function(mssc2) {
  if (class(mssc2) != "MSSC2") {
    stop("Input mssc2 is not an objet of MSSC2 S3 class.")
  }
  if (is.null(mssc2$init_mgsnb)) {
    stop("mssc2 init_mgsnb is NULL. Use fitGenewiseNBModel to set it.")
  }
  init_mgsnb <- mssc2$init_mgsnb
  mu <- init_mgsnb$mgsnb[, 1]
  ### nb_r in negative binomial
  nb_r <- init_mgsnb$mgsnb[, 2]
  log_of_r <- log(nb_r)
  ### ngene by ncond
  mucond <- init_mgsnb$mgsnb[, 3:(2 + mssc2$ncond)]
  ### ngene by nind
  muind <- init_mgsnb$mgsnb[, (2 + mssc2$ncond + 1):(2 + mssc2$ncond + mssc2$nind)]
  mssc2$glmip <- list(
    mu = mu,
    r = log_of_r,
    mucond = mucond,
    muind = muind
  )
  invisible(mssc2)
}

#' Run MAP optimization for generalized lineaer model.
#' @param mssc2 MSSC2 object
#' @param data list, generated by toStanInput
#' @return MSSC2 object
#' @export
runGLMAP <- function(mssc2, data) {
  if (class(mssc2) != "MSSC2") {
    stop("Input mssc2 is not an objet of MSSC2 S3 class.")
  }

  if (is.null(mssc2$glmodel)) {
    stop("mssc2 model has not been compiled.")
  }
  if (is.null(mssc2$gmlip)) {
    stop("mssc2 hasn't initialize the parameters.")
  }
  sp <- mssc2$stanparams
  mssc2$glmopt <- mssc2$glmodel$optimize(
    data = data,
    init = list(mssc2$glmip),
    seed = mssc2$seed,
    refresh = sp$opt_refresh,
    max_iter = sp$opt_max_iter,
    opt_method = sp$opt_metod,
    init_alpha = sp$init_alpha,
    tol_obj = sp$tol_obj,
    tol_rel_obj = sp$tol_rel_obj,
    tol_grad = sp$tol_grad,
    tol_rel_grad = sp$tol_rel_grad,
    tol_param = sp$tol_param,
    history_size = sp$history_size
  )
  return(mssc2)
}

#' Extract draws/samples from the generalized linear model.
#' @param mssc2 MSSC2 object
#' @param param string, name of the parameter we want to extract
#' @param genenms vector of string, length of ngene
#' @param method string, "opt" or "vi", default is "vi"
#' @return NA (when have error) or vector/matrix (depend on param)
#' - samples of the parameter for all the genes
#' @export
extractGLMDraws <- function(mssc2, param, genenms = NULL, method = "opt") {
  if (class(mssc2) != "MSSC2") {
    stop("Input mssc2 is not an objet of MSSC2 S3 class.")
  }
  if (method == "vi") {
    fit <- mssc2$glmvi
  } else {
    fit <- mssc2$glmopt
  }
  if (is.null(fit)) {
    stop(paste(method, " is null. Might becuase it hasn't been run or have troubles."))
  }
  if (!(param %in% mssc2$allparamsglm)) {
    stop(paste(
      param, "is not in the list of",
      paste(mssc2$allparamsglm, collapse = ","), "."
    ))
  }
  if ((!is.null(genenms)) & (length(genenms) != mssc2$ngene)) {
    stop(paste(
      length(genenms), "is not equal to ngene",
      mssc2$ngene, "in the model."
    ))
  }
  tryCatch(
    {
      if (param %in% c("mu", "r", "nb_r")) {
        t <- fit$draws(getNameOfVectorFromCmdstanr(nm = param, n = mssc2$ngene))
        names(t) <- genenms
        return(invisible(t))
      }
      if (param %in% c("mucond")) {
        t <- fit$draws(getNameOfMatrixFromCmdstanr(nm = param, nr = mssc2$ngene, nc = mssc2$ncond))
        return(invisible(split_matrix_col(mat = t, secondim = mssc2$ngene, secondimNms = genenms)))
      }
      if (param %in% c("muind")) {
        t <- fit$draws(getNameOfMatrixFromCmdstanr(nm = param, nr = mssc2$ngene, nc = mssc2$nind))
        return(invisible(split_matrix_col(mat = t, secondim = mssc2$ngene, secondimNms = genenms)))
      }
    },
    error = function(e) {
      warning(e)
      return(invisible(NA))
    }
  ) ## of trycatch block
}

#' Extract all parameters' samples from generalised linear model.
#' @param mssc2 MSSC2 object
#' @param genenms NULL or vector of strings.
#' @param method string, "vi" or "opt", default is "opt"
#' @return list, samples for different parameters.
#' @export
extractGLMDrawsAll <- function(mssc2, genenms = NULL, method = "opt") {
  if (class(mssc2) != "MSSC2") {
    stop("Input mssc2 is not an objet of MSSC2 S3 class.")
  }
  if (method == "vi") {
    fit <- mssc2$glmvi
  } else {
    fit <- mssc2$glmopt
  }
  if (is.null(fit)) {
    stop(paste(method, " is null. Might becuase it hasn't been run or have troubles."))
  }
  est_params <- lapply(mssc2$allparamsglm, function(nm) {
    extractGLMDraws(mssc2 = mssc2, genenms = genenms, method = method)
  })
  names(est_params) <- mssc2$allparamsglm
  invisible(est_params)
}
