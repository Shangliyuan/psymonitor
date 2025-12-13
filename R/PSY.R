#' @title Estimate PSY's BSADF sequence of test statistics
#'
#' @description \code{PSY} implements the real time bubble detection procedure
#'   of Phillips, Shi and Yu (2015a,b)
#'
#' @param y   A vector. The data.
#' @param swindow0 A positive integer. Minimum window size (default = \eqn{T
#'   (0.01 + 1.8/\sqrt{T})}, where \eqn{T} denotes the sample size)
#' @param IC  An integer. 0 for fixed lag order (default), 1 for AIC and 2 for
#'   BIC (default = 0).
#' @param adflag  An integer, lag order when IC=0; maximum number of
#'   lags when IC>0 (default = 0).
#'
#' @return Vector, BSADF test statistic.
#'
#' @references Phillips, P. C. B., Shi, S., & Yu, J. (2015a). Testing for
#'   multiple bubbles: Historical episodes of exuberance and collapse in the S&P
#'   500. \emph{International Economic Review}, 56(4), 1034--1078.
#' @references Phillips, P. C. B., Shi, S., & Yu, J. (2015b). Testing for
#'   multiple bubbles: Limit Theory for Real-Time Detectors. \emph{International
#'   Economic Review}, 56(4), 1079--1134.
#'
#' @export
#'
#'
#' @examples
#'
#' y     <- rnorm(80)
#' bsadf <- PSY(y, IC = 0, adflag = 1)
#'


PSY <- function(y, swindow0, IC = 0, adflag = 0,
                useParallel = TRUE, nCores) {

  t <- length(y)

  if (missing(swindow0)) {
    swindow0 <- floor(t * (0.01 + 1.8 / sqrt(t)))
  }

  r2_seq <- swindow0:t

  if (useParallel) {
    # Set up parallel backend only when needed
    if (missing(nCores)) {
      nCores <- parallel::detectCores() - 1
    }
    nCores <- max(1L, nCores)  # ensure at least 1 core

    cl <- parallel::makeCluster(nCores)
    doParallel::registerDoParallel(cl)

    bsadf <- foreach::foreach(r2 = r2_seq, .combine = c) %dopar% {
      rwadft <- numeric(r2 - swindow0 + 1)
      for (r1 in 1:(r2 - swindow0 + 1)) {
        rwadft[r1] <- as.numeric(ADF(y[r1:r2], IC, adflag))
      }
      max(rwadft)
    }

    parallel::stopCluster(cl)  # Properly shut down cluster

  } else {
    # Sequential execution: no cluster, no extra ports
    bsadf <- foreach::foreach(r2 = r2_seq, .combine = c) %do% {
      rwadft <- numeric(r2 - swindow0 + 1)
      for (r1 in 1:(r2 - swindow0 + 1)) {
        rwadft[r1] <- as.numeric(ADF(y[r1:r2], IC, adflag))
      }
      max(rwadft)
    }
  }

  return(bsadf)
}

