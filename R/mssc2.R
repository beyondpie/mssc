#' Create MSSC2 model object
#'
#' @return MSSC2 object
#' @export
new_MSSC2 <- function(modelpath = NULL, glmodelpath = NULL) {
  r <- structure(list(
    modelpath = modelpath,
    model = NULL,
    init_mgsnb = NULL,
    ## predefined params, not model-related hyper params.
    methodparams = list(
      default_mu = 0.0, default_r = 20,
      min_var = 4.0, min_varofcond = 0.25,
      min_varofind = 0.25, min_tau2 = 0.25
    ),
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
      adapt_engaged = FASLE,
      eta = 0.1
    ),
    allparams = c(
      "centerofmu", "varofmu", "mu", "centerofr",
      "varofr", "r", "varofcond", "mucond",
      "tau2", "centerofind", "varofind", "muind"
    ),
    glmodelpath = glmodelpath,
    glmodel = NULL,
    allparamsglm = c(
      "mu", "nb_r", "mucond", "muind"
    )
  ), class = "MSSC2")
  
  ## only compile model. glm need to be explicitly compiled when in use.
  r$model <- compileStanModel(
    model_path = r$model_path,
    use_thread = NULL, use_mpi = NULL,
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
  ## create a new mssc2 since R "copy when modified".
  mssc2$init_mgsnb <- init_mgsnb
  invisible(mssc2)
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
  if(is.null(mssc2$init_mgsnb)) {
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
  mssc2$hp <- hp
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
  if(is.null(mssc2$init_mgsnb)) {
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
  ncond <- max(cond)
  nind <- max(ind)
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
  ngene <- nrow(cnt)
  raw_muind <- (muind - rep_row(centerofind, n = ngene)) %*%
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
      raw_muind = raw_muind)
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
  invisible(list(
    ncell = ncol(cnt),
    nind = max(ind),
    ncond = max(cond),
    ngene = nrow(cnt),
    s = s, cond = cond,
    ind = ind, y = t(cnt)), hp)
}

