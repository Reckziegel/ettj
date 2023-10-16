#' Estimates the squared root of the mean squared error (rmse) for an OLS regression of the Svensson (1994) factor model
#'
#' @description
#' Performs an OLS regression using the 4 factor yield curve model:
#'
#' \eqn{y = \beta_{1}((1-e^{(-\tau^{'}  \lambda_{1})})/(\tau^{'}\lambda_{1})) +
#'                \beta_{2}((1-e^{(-\tau^{'}  \lambda_{1})})/(\tau^{'}\lambda_{1}) - e^{(-\tau^{'}  \lambda_{1})}) +
#'                \beta_{3}((1-e^{(-\tau^{'}  \lambda_{2})})/(\tau^{'}\lambda_{2}) - e^{(-\tau^{'}  \lambda_{2})})}
#'
#'
#' @param lambda Decay parameters
#' @param Y Observed yields
#' @param tau Time to maturity
#'
#' @return squared root of mean squared error
#' @export
#'
#' @examples
#'
rmse_sv <- function(lambda, Y, tau){
  lam1 = lambda[[1]]
  lam2 = lambda[[2]]
  H = matrix(1,ncol(Y),4)
  H[,2] = (1-exp(-t(tau)*lam1))/(t(tau)*lam1)
  H[,3] = (1-exp(-t(tau)*lam1))/(t(tau)*lam1) - exp(-t(tau)*lam1)
  H[,4] = (1-exp(-t(tau)*lam2))/(t(tau)*lam2) - exp(-t(tau)*lam2)
  B = solve(t(H)%*%H, tol = 1*10^-500)%*%(t(H)%*%t(Y))
  Y_hat = H%*%B
  Y_hat = t(Y_hat)
  erro = Y_hat - Y
  rmse = sqrt(mean(as.matrix(erro)^2))
  return(sum(rmse))
}
