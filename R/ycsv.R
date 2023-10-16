#' svensson (1994) 4 factor model formula
#'
#' @description
#' Returns the yields associated with a set of maturities, using decay parameters and factors with the Svensson (1994) formula
#'
#'
#' @param lambda Decay parameters
#' @param factors The estimated factors Beta 1 through 4
#' @param tau Time to maturity
#'
#' @return vector of yields corresponding to the maturities
#' @export
#'
#' @examples
ycsv <- function(lambda, factors, tau){
  lam1 = lambda[[1]]
  lam2 = lambda[[2]]
  H = matrix(1,length(tau),4)
  H[,2] = (1-exp(-t(tau)*lam1))/(t(tau)*lam1)
  H[,3] = (1-exp(-t(tau)*lam1))/(t(tau)*lam1) - exp(-t(tau)*lam1)
  H[,4] = (1-exp(-t(tau)*lam2))/(t(tau)*lam2) - exp(-t(tau)*lam2)
  y_hat = H%*%factors # returns an 1xn vector
  return(t(y_hat))
}
