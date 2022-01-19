#' @importFrom stats coef
#' @importFrom glmnet glmnet
#' @importFrom purrr map map_dbl
#' @export
motl <- function(X, Y, w, lambda, eta, cv = TRUE, k = 5, ...)
{
  n <- nrow(X)
  p <- ncol(X)
  if (length(lambda) == 1 & length(eta) == 1)
    return(fit_motl(X = X, Y = Y, w = w, lambda = lambda, eta = eta, ...))

  ## here is our output from the double map
  raw_estimates <- map(lambda, function(lam)
    map(eta, function(et)
      fit_motl(X = X, Y = Y, w = w, eta = et, lambda = lam, ...)))

  names(raw_estimates) <- as.character(lambda)
  raw_estimates <- lapply(raw_estimates, function(x) {
    names(x) <- as.character(eta)
    return(x)
  })
  if (cv == FALSE)
    return(raw_estimates)

  ## run k fold cross validate to estimate the prediction error, 
  ## for a given lambda,eta pair
  get_error <- function(lam, et) {
    random_order <- sample(n, n)
    subsamples <- split(random_order, 1:k)
    mean(sapply(subsamples, function(idx){
      X.train <- X[-idx,]
      Y.train <- Y[-idx]
      X.test <- X[idx,]
      Y.test <- Y[idx]
      crossprod(Y.test - X.test %*%
                fit_motl(X = X.train, Y = Y.train, w = w,
                         lambda = lam, eta = et, ...)) / length(Y.test)
    }))
  }
  errors <- map(lambda, function(lam) 
   map_dbl(eta, function(et) get_error(lam, et)))
  errors_mat <- matrix(unlist(errors), nrow = length(eta), ncol = length(lambda)) ### rows are etas, columns are lambdas
  errors_mat_t <- t(errors_mat) ## transpose so that lambdas are rows, etas are columns
  rownames(errors_mat_t) <- as.character(lambda)
  colnames(errors_mat_t) <- as.character(eta)
  choice <- which(errors_mat_t == min(errors_mat_t), arr.ind = TRUE) ## coordinate of optimal lambda, eta

  list(estimates = raw_estimates,
       error_mat = errors_mat_t,
       optimal_lambda = choice[1],
       optimal_eta = choice[2])
}
