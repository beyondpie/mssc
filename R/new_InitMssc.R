##' Create InitMssc object
##'
##' Initialize a InitMssc object in order to initialize the parameters
##' of MSSC.
##' 
##' @return  InitMssc
##' @export
##' @example
##' new_InitMssc()
new_InitMssc <- function() {
  structure(list(
    nb_model = NULL,
    gamma_alpha = 0.05,
    gamma_beta = 0.05,
    mu = 0.0,
    r = 20,
    big_r = 500,
    min_varofmu = 4.0,
    min_varofcond = 0.25,
    min_varofind = 0.25,
    min_tau2 = 0.25
  ), class = "InitMssc")
}


