#' Create MSSC2 model object
#'
#' @return MSSC2 object
#' @export
new_MSSC2 <- function(modelpath = NULL, glmodelpath = NULL) {
  structure(list(
    nind = NULL,
    ncond = 2,
    modelpath = modelpath,
    model = NULL,
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
}

#' @export
initMSSC2Params <- function(mssc2)
