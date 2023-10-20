#' Estimates the squared root of the mean squared error (rmse) for an OLS regression of the 5 factor model factor model.
#'
#' @description
#' Performs an OLS regression using the 5 factor yield curve model:
#'
#' \deqn{y_{\tau} = \beta_{1} +
#'                  \beta_{2}((1-e^{(-\tau^{'}  \lambda_{1})})/(\tau^{'}\lambda_{1})) +
#'                  \beta_{3}((1-e^{(-\tau^{'}  \lambda_{1})})/(\tau^{'}\lambda_{1}) - e^{(-\tau^{'}  \lambda_{1})}) +
#'                  \beta_{4}((1-e^{(-\tau^{'}  \lambda_{2})})/(\tau^{'}\lambda_{2}) - e^{(-\tau^{'}  \lambda_{2})}) +
#'                  \beta_{5}((1-e^{(-\tau^{'}  \lambda_{3})})/(\tau^{'}\lambda_{3}) - e^{(-\tau^{'}  \lambda_{3})})}
#' The independent variable is the yield to maturity, \eqn{y_{\tau}}, the dependent variables are the factor loadings,
#' the bracketed terms that multiply the \eqn{\beta} factors. The results are compared with the empirical yields and summarized 
#' in the rmse results. This function is called by the factors_sv function in order to determine the decay parameters, 
#' \eqn{\lambda}, which minimize the error of the model in relation to the observed values.
#'
#' @param lambda Decay parameters, associated with the slope, curvature 1 and curvature 2 loadings.Numeric values.
#' @param Y Observed yields. Numeric Vector.
#' @param tau Time to maturity. Numeric. Must match the length of Y and be in the same base, i.e, annual, monthly etc.
#'
#' @return Squared root of mean squared error.
#' @export
#'
#' 
#'
rmse_sv <- function(lambda, Y, tau){
  lam1 = lambda[[1]]
  lam2 = lambda[[2]]
  lam3 = lambda[[3]]
  H = matrix(1,ncol(Y),5)
  H[,2] = (1-exp(-t(tau)*lam1))/(t(tau)*lam1)
  H[,3] = (1-exp(-t(tau)*lam1))/(t(tau)*lam1) - exp(-t(tau)*lam1)
  H[,4] = (1-exp(-t(tau)*lam2))/(t(tau)*lam2) - exp(-t(tau)*lam2)
  H[,5] = (1-exp(-t(tau)*lam3))/(t(tau)*lam3) - exp(-t(tau)*lam3)
  B = solve(t(H)%*%H, tol = 1*10^-500)%*%(t(H)%*%t(Y))
  Y_hat = H%*%B
  Y_hat = t(Y_hat)
  erro = Y_hat - Y
  rmse = sqrt(mean(as.matrix(erro)^2))
  return(sum(rmse))
}
