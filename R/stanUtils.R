#' Compile stan model with cmdstanr
#'
#' NOTE previous name: init_stan_model
#' @param model_path string, stan script
#' @param use_thread bool/NULL, when reduce_sum or map_rect in script, set it as TRUE.
#' When running, remember to export STAN_NUM_THREADS=4. Default is NULL.
#' @param use_mpi bool/NULL, when map_rect in script, set it as TRUE. Deault is NULL.
#' @param use_opencl bool/NUYLL, enable the OpenCL backend. Deault is NULL.
#' When set is as TRUE, remember to set STAN_OPENCL=true when compiling.
#' https://mc-stan.org/docs/2_26/cmdstan-guide/parallelization.html
#' - The Stan model compiled with STAN_OPENCL can also be supplied the
#'  OpenCL platform and device IDs of the target device. These IDs
#'  determine the device on which to run the OpenCL-supported functions
#'  on. You can list the devices on your system using the clinfo
#'  program. If the system has one GPU and no OpenCL CPU runtime, the
#'  platform and device IDs of the GPU are typically 0. In that case you
#'  can also omit the OpenCL IDs as the default 0 IDs are used in that
#'  case.
#'  We supply these IDs when starting the executable as shown below:
#'  path/to/model sample data file=data.json opencl platform=0 device=1
#' 
#' NOTE[BUG]:  use_xxx = NULL when we don't want to use them.
#' When I set them as FALSE, they still treat them as TRUE.
#' @return compiled stan model
#' @export
compileStanModel <- function(model_path,
                             use_thread = NULL,
                             use_mpi = NULL,
                             use_opencl = NULL) {
  if(!file.exists(model_path)){
    stop(paste(model_path, "does not exist."))
  }
  r <- suppressMessages(cmdstanr::cmdstan_model(
    stan_file = model_path,
    compile = TRUE,
    quiet = TRUE,
    pedantic = TRUE,
    cpp_options = list(stan_threads = use_thread,
                       stan_map = use_mpi,
                       stan_opencl = use_opencl),
    stanc_options = list(),
    force_recompile = FALSE
  ))
  invisible(r)
}

#' Get names of vector in cmdstanr.
#' @param nm string, name of variable in stan script
#' @param n integer, length of the variable.
#' @param l string, left bracket, default "["
#' @param r string, right brackt, default "]"
#' @return vector of string
#' @export
getNameOfVectorFromCmdstanr <- function(nm = "MuInd", n = 10L, l = "[", r = "]") {
  invisible(vapply(seq_len(n), function(i) {
    paste0(nm, ls, i, rs)
  }, FUN.VALUE = "MuInd[1]"))
}

#' Get names of matrix in cmdstanr.
#' @param nm string, name of variable in stan script
#' @param nr integer, nrow of the variable.
#' @param nc integer, ncolof the variable.
#' @param l string, left bracket, default "["
#' @param r string, right brackt, default "]"
#' @return vector of string
#' @export
getNameOfMatrixFromCmdstanr <- function(nm = "MuInd", nr, nc, l = "[",r = "]") {
  r <- rep("", nr * nc)
  for (i in 1:nr) {
    for (j in 1:nc) {
      r[nc * (i - 1) + j] <- paste0(nm, ls, i, ",", j, rs)
    }
  }
  invisible(r)
}

#' Replicate a vector in row-wise for n times.
#'
#' NOTE previous name: rep_row
#' @param x vector of numeric
#' @param n integer, times of repeat in row-wise
#' @return matrix
#' @export  
repVecRowise <- function(x, n) {
  matrix(rep(x, each = n), nrow = n, byrow = FALSE)
}

#' Replicate a vector in column-wise for n times.
#'
#' NOTE previous name: rep_col
#' @param x vector of numeric
#' @param n integer, times of repeat in column-wise
#' @return matrix
#' @export
repVecColwise <- function(x, n) {
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

#' Split matrix columns into a matrix, which is row-major order
#' @param mat matrix, size of nsample by (second_dim * third_dim)
#' @param secondim integer, secondim we want
#' @param secondimNms vector of string, names for secondim.
#' @return three-dim array
#' @export
split_matrix_col <- function(mat, secondim, secondimNms = NULL) {
  t <- vapply(seq_len(nrow(mat)), function(i) {
    ## byrow = TRUE, means we order the elements row by row.
    ## so each time we use ncol element to fill a row.
    matrix(mat[i, ], nrow = secondim, byrow = TRUE)
  }, FUN.VALUE = matrix(mat[1, ], nrow = secondim, byrow = TRUE))
  r <- aperm(t, c(3, 1, 2))
  if (!is.null(secondimNms)) {
    dimnames(r)[[2]] <- secondimNms
  }
  return(invisible(r))
}
