#' Test MSSC2
#'
#' A modified data set from a public PBMC single-cell RNA sequencing data,
#' which has 200 genes, 6799 cells from 10 individuals and 2 conditions, is used.
#' @format A list with 6 elements
#' \describe{
#'   \item{sumcnt}{total number of counts per cell}
#'   \item{s}{normaling constant per cell}
#'   \item{ind}{which individual a cell belongs to, start from 1.}
#'   \item{cond}{which condition a cell belongs to, start from 1.}
#'   \item{y2c}{count matrix, ngene by ncell}
#'   \item{est_snb_mat}{matrix, estimated parameters for different genes}
#' }
#' @source A public PBMC single-cell RNA sequencing data

library(mssc)

pbmc <- readRDS(system.file("extdata", "refPBMC.rds", package = "mssc", mustWork = TRUE))
nind <- max(pbmc$ind)

## it takes some time (around several minitutes) to let cmdstanr compile mssc2 stan script.
mssc <- new_MSSC2(modelpath = system.file("stan", "mssc2.stan", package = "mssc", mustWork = TRUE),
                   glmodelpath = system.file("stan", "glm.stan", package = "mssc", mustWork = TRUE))

mssc <- initMSSC2GenewiseParams(mssc2 = mssc, cnt = pbmc$y2c[1:10,], s = pbmc$s,
                                cond = pbmc$cond, ind = pbmc$ind)
mssc <- setStanHyperparams(mssc2 = mssc)
mssc <- setStanInitialParams(mssc2 = mssc)
data <- toStanInput(cnt = pbmc$y2c[1:10, ], s = pbmc$s, cond = pbmc$cond,
                    ind = pbmc$ind, hp = mssc$modelhp)

## still need to compile.
mssc$model$compile()
## run variational inference
mssc <- runVI(mssc2 = mssc, data = data)
## run MAP
## A bug: when mssc is updated, we need to re-compile.
mssc$model$compile()
mssc <- runMAP(mssc2 = mssc, data = data)

## analysis
vi_est_params <- extractDrawsAll(mssc2 = mssc,
                                 genenms = rownames(pbmc$y2c[1:10,]),
                                 method = "vi")
str(vi_est_params)
opt_est_params <- extractDrawsAll(mssc2 = mssc,
                                  genenms = rownames(pbmc$y2c[1:10,]),
                                  method = "opt")
vi_mucond <- extractDraws(mssc2 = mssc, param = "mucond",
                       genenms = rownames(pbmc$y2c[1:10, ]),
                       method = "vi")
opt_mucond <- extractDraws(mssc2 = mssc, param = "mucond",
                       genenms = rownames(pbmc$y2c[1:10, ]),
                       method = "opt")



rankings <- model$get_ranking_statistics(
  mucond = mucond,
  two_hot_vec = c(1, -1)
)
str(rankings)

psis <- model$psis()
print(psis$psis)
rankings <- model$get_ranking_statistics(
  mucond = mucond,
  two_hot_vec = c(1, -1)
)
str(rankings)

## test glm
init_params_of_glm <- model$init_glm_params(cnt = pbmc$y2c[1:10,],
                                            s = pbmc$s,
                                            cond = pbmc$cond,
                                            ind = pbmc$ind)
model$run_glm_opt(data = data, list_wrap_ip = list(init_params_of_glm))
est_params_of_glm <- model$extract_draws_all_from_glm
mucond_of_glm <- model$extract_draws_from_glm(
  param = "mucond", ngene = 10,
  genenms = rownames(pbmc$y2c[1:10, ])
)
ranking_from_glm <- model$get_opt_ranking_statistic(
  mucond = mucond_of_glm,
  two_hot_vec = c(1, -1)
)
str(ranking_from_glm)
