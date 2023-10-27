#' 5 factor model formula
#'
#' @description
#' Returns the yields associated with a set of maturities, using decay parameters and factors with the 5 factor formula.
#'
#'
#' @param lambda Decay parameters, associated with the slope, curvature 1 and curvature 2 loadings.Numeric values.
#' @param factors Estimated factors (see factors_5f). Numeric Vector.
#' @param tau Time to maturity. Numeric. Must match the length of Y and be in the same base, i.e, annual, monthly etc.
#'
#' @return Vector of yields corresponding to the maturities.
#' @export
#' 
#'
#' @examples
#' 
#' #library('Quandl')
#' 
#' #dados <- Quandl("USTREASURY/YIELD")
#' dados <- ettj_package_data
#' yields <- dados[1,2:ncol(dados)]
#' maturidades <- c(1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30);
#' factors <- factors_5f(c(.9, .035, .035), yields, maturidades)
#' yc <- yc5f(factors[6:8], factors[1:5], maturidades)

yc5f <- function(lambda, factors, tau){
  lam1 = lambda[[1]]
  lam2 = lambda[[2]]
  lam3 = lambda[[3]]
  H = matrix(1,length(tau),5)
  H[,2] = (1-exp(-t(tau)*lam1))/(t(tau)*lam1)
  H[,3] = (1-exp(-t(tau)*lam1))/(t(tau)*lam1) - exp(-t(tau)*lam1)
  H[,4] = (1-exp(-t(tau)*lam2))/(t(tau)*lam2) - exp(-t(tau)*lam2)
  H[,5] = (1-exp(-t(tau)*lam3))/(t(tau)*lam3) - exp(-t(tau)*lam3)
  y_hat = H%*%factors # returns an 1xn vector
  return(t(y_hat))
}
