#' @importFrom stats coef
#' @importFrom glmnet glmnet
#' @export
motl <- function(X, Y, w, lambda, eta, k = 5, ...)
{
  if (length(lambda) == 1 & length(eta) == 1)
    return(fit_motl(X = X, Y = Y, w = w, lambda = lambda, eta = eta))

}
