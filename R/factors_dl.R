#' Estimates the 3 factors and decay parameter in the Diebold and Li (2006) factor model. The decay parameters are optimized using the two step approach.
#'
#' @description
#' Determines the optimal decay parameters calling the rmse function.
#' Then, performs an OLS regression using the 3 factor yield curve model in order to determine the optimal \eqn{\lambda} decay parameters:
#' \deqn{y_{\tau} = \beta_{1}((1-e^{(-\tau^{'}  \lambda_{1})})/(\tau^{'}\lambda_{1})) +
#'                  \beta_{2}((1-e^{(-\tau^{'}  \lambda_{1})})/(\tau^{'}\lambda_{1}) - e^{(-\tau^{'}  \lambda_{1})})}
#'
#' @param lambda Initial guess for the decay parameter. A numeric values.
#' @param yields Observed yields. Numeric Vector.
#' @param maturidades Time to maturity. Numeric. Must match the length of Y and be in the same base, i.e, annual, monthly etc.
#'
#' @return A matrix containing the 3 factors and decay parameter, respectively.
#' @export
#'
#' @examples
#' 
#' library(Quandl)
#' 
#' dados <- Quandl("USTREASURY/YIELD")
#' yields <- dados[1,2:ncol(dados)]
#' maturidades <- c(1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30);
#' factors <- factors_dl(c(.9,.035), yields, maturidades)
#' yc <- ycdl(factors[1], factors[1:3], maturidades)
#'
factors_dl <- function(lambda, yields, maturidades){
  # Treating data -----------------------------------------------------------
  index <- full_join(as_tibble(which(is.na(yields))), as_tibble(which(is.na(maturidades))))
  
  index <- as.double(as.matrix(index))
  
  if(is_empty(index)){
    yields <- yields
    maturidades <- maturidades
  }else{
    yields <- yields[-c(index),]
    maturidades <- maturidades[-c(index),]
  }
  
  
  # Calculating optimal decay parameters --------------------------------
  
  
  optlamfun <- optim(lambda, rmse_dl, Y = yields, tau = maturidades , control = list(maxit = 10000, reltol = 0.00000001))
  
  lam1 <- optlamfun$par[1]
  
  # Estimating Betas --------------------------------------------------------
  
  
  H <-  matrix(NA,length(yields),3)
  B <-  matrix(NA,3,1)
  
  H[1:length(yields),1] <- 1
  H[1:length(yields),2] <- (1-exp(-t(maturidades)*lam1))/(t(maturidades)*lam1)
  H[1:length(yields),3] <- (1-exp(-t(maturidades)*lam1))/(t(maturidades)*lam1) - exp(-t(maturidades)*lam1)
  B[,1] <- solve(t(H)%*%H, tol = 10^-50)%*%(t(H)%*%t(yields))
  
  
  final <- rbind(B, lam1); rownames(final) <- c('Beta 1', 'Beta 2', 'Beta 3', 'lambda 1')
  
  return(final)
}