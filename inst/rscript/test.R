test <- function() {
  pbmc <- readRDS("snb_pool_ref_pbmc.rds")
  nind <- max(pbmc$ind)
  model <- High2$new(
    stan_snb_path = here::here("src", "mssc", "stan", "snb.stan"),
    stan_high2_path = here::here("src", "mssc", "stan", "mssc_2-0.stan"),
    stan_glm_path = here::here("src", "mssc", "stan", "glm.stan"),
    nind = nind,
    tol_rel_obj = 0.0001,
    adapt_engaged = FALSE
  )

  init_params <- model$init_params(
    cnt = pbmc$y2c[1:10, ],
    s = pbmc$s,
    cond = pbmc$cond,
    ind = pbmc$ind
  )

  data <- model$to_model_data(
    cnt = pbmc$y2c[1:10, ],
    s = pbmc$s,
    cond = pbmc$cond,
    ind = pbmc$ind,
    hp = init_params$hp
  )

  ## variational inference
  model$run(data = data, list_wrap_ip = list(init_params$ip))

  est_params <- model$extract_draws_all(
    ngene = 10,
    genenms = rownames(pbmc$y2c[1:10, ])
  )
  str(est_params)

  mucond <- model$extract_draws(
    param = "mucond", ngene = 10,
    genenms = rownames(pbmc$y2c[1:10, ])
  )
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

  ## optimization
  model$run_opt(data = data, list_wrap_ip = list(init_params$ip))
  est_params <- model$extract_draws_all(
    ngene = 10,
    genenms = rownames(pbmc$y2c[1:10, ]),
    method = "opt"
  )
  str(est_params)

  mucond <- model$extract_draws(
    param = "mucond", ngene = 10,
    genenms = rownames(pbmc$y2c[1:10, ]),
    method = "opt"
  )
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
}

test()
