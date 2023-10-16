#' Estimates the 4 factors and 2 decay parameters in the Svensson (1994) factor model. The decay parameters are optimized using the two step approach.
#'
#' @param lambda initial guess for the decay parameters
#' @param yields A vector of yields
#' @param maturidades A vector of maturities. Be sure to put them in the same base as yield, i.e., per annum, per month etc.
#'
#' @return A matrix containing the 4 factors ad 2 decay parameters in this order
#' @export
#'
#' @examples
#'
factors_sv <- function(lambda, yields, maturidades){
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


  optlamfun <- optim(lambda, rmse_sv, Y = yields, tau = maturidades , control = list(maxit = 10000, reltol = 0.00000001))

  lam1 <- optlamfun$par[1]
  lam2 <- optlamfun$par[2]

  # Estimating Betas --------------------------------------------------------


  H <-  matrix(NA,length(yields),4)
  B <-  matrix(NA,4,1)

  H[1:length(yields),1] <- 1
  H[1:length(yields),2] <- (1-exp(-t(maturidades)*lam1))/(t(maturidades)*lam1)
  H[1:length(yields),3] <- (1-exp(-t(maturidades)*lam1))/(t(maturidades)*lam1) - exp(-t(maturidades)*lam1)
  H[1:length(yields),4] <- (1-exp(-t(maturidades)*lam2))/(t(maturidades)*lam2) - exp(-t(maturidades)*lam2)
  B[,1] <- solve(t(H)%*%H, tol = 10^-50)%*%(t(H)%*%t(yields))


  final <- rbind(B, lam1, lam2); rownames(final) <- c('Beta 1', 'Beta 2', 'Beta 3', 'Beta 4', 'lambda 1', 'lambda 2')

  return(final)
}
