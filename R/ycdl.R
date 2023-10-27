#' Diebold and Li (2006) 3 factor model formula
#'
#' @description
#' Returns the yields associated with a set of maturities, using decay parameters and factors with the Diebold and Li (2006) formula.
#'
#' @param lambda Decay parameter.Numeric.
#' @param factors Observed yields. Numeric Vector.
#' @param tau Time to maturity. Numeric. Must match the length of Y and be in the same base, i.e, annual, monthly etc.
#'
#' @return Vector of yields corresponding to the maturities.
#' @export
#'
#' @examples
#' 
#' #library('Quandl')
#' 
#' #dados <- Quandl("USTREASURY/YIELD")
#' dados <- ettj_package_data
#' yields <- dados[1,2:ncol(dados)]
#' maturidades <- c(1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30);
#' factors <- factors_dl(c(.9,.035), yields, maturidades)
#' yc <- ycdl(factors[4], factors[1:3], maturidades)
#' 
#' 
ycdl <- function(lambda, factors, tau){
  lam1 = lambda[[1]]
  H = matrix(1,length(tau),3)
  H[,2] = (1-exp(-t(tau)*lam1))/(t(tau)*lam1)
  H[,3] = (1-exp(-t(tau)*lam1))/(t(tau)*lam1) - exp(-t(tau)*lam1)
  y_hat = H%*%factors # returns an 1xn vector
  return(t(y_hat))
}