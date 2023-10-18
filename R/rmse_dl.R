#' Estimates the squared root of the mean squared error (rmse) for an OLS regression of the Diebold Li (2006) factor model.
#'
#' @description
#' Performs an OLS regression using the 3 factor yield curve model:
#' \deqn{y_{\tau} = \beta_{1} + 
#'                  \beta_{2}((1-e^{(-\tau^{'}  \lambda_{1})})/(\tau^{'}\lambda_{1})) +
#'                  \beta_{3}((1-e^{(-\tau^{'}  \lambda_{1})})/(\tau^{'}\lambda_{1}) - e^{(-\tau^{'}  \lambda_{1})})} 
#' The independent variable is the yield to maturity, \eqn{y_{\tau}}, the dependent variables are the factor loadings,
#' the bracketed terms that multiply the \eqn{\beta} factors. The results are compared with the empirical yields and summarized 
#' in the rmse results. This function is called by the factors_sv function in order to determine the decay parameters, 
#' \eqn{\lambda}, which minimize the error of the model in relation to the observed values.
#'
#' @param lambda Decay parameters. A pair of numeric values.
#' @param Y Observed yields. Numeric Vector.
#' @param tau Time to maturity. Numeric. Must match the length of Y and be in the same base, i.e, annual, monthly etc.
#'
#' @return Squared root of mean squared error.
#' @export
#'
#' 
rmse_dl <- function(lambda, Y, tau){
  lam1 = lambda[[1]]
  H = matrix(1, ncol(Y), 3)
  H[,2] = (1-exp(-t(tau)*lam1))/(t(tau)*lam1)
  H[,3] = (1-exp(-t(tau)*lam1))/(t(tau)*lam1) - exp(-t(tau)*lam1)
  B = solve(t(H)%*%H, tol = 1*10^-500)%*%(t(H)%*%t(Y))
  Y_hat = H%*%B
  Y_hat = t(Y_hat)
  erro = Y_hat - Y
  rmse = sqrt(mean(as.matrix(erro)^2))
  return(sum(rmse))
}