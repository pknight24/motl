#' Mixed-outcome transfer learning
#' @param X An n-by-p design matrix.
#' @param Y An outcome vector of length n.
#' @param w An estimated coefficient vector of length p.
#' @param lambda A vector of positive tuning parameters corresponding to the penalty on the Euclidean norm of beta.
#' @param eta A vector of positive tuning parameters corresponding to the penalty on the inner product of w and beta.
#' @param cv Logical, indicates whether to evaluate the cross validation error of the estimator at each pair of tuning parameters (lambda, eta).
#' @param k Integer, for k-fold cross validation when cv = TRUE.
#' @param family String indicating which model to fit: currently supports
#' 'gaussian' for continuous outcomes and 'binomial' for binary outcomes.
#' @param ... Additional arguments for the fitting function.
#' @return A list of estimated coefficients corresponding to each pair of tuning parameters. If cv = TRUE, also returns the cross-validation error for each pair of tuning parameters in the form of a matrix, as well as the optimal lambda and eta.
#' @importFrom stats coef optim
#' @importFrom glmnet glmnet
#' @importFrom purrr map map_dbl
#' @export
motl <- function(X, Y, w, lambda, eta, cv = TRUE, k = 5, family = "gaussian",...)
{
  ## cross validation for non-continuous data isn't supported yet
  if (family != "gaussian") cv <- FALSE
  n <- nrow(X)
  p <- ncol(X)
  if (length(lambda) == 1 & length(eta) == 1)
    return(fit_motl(X = X, Y = Y, w = w, lambda = lambda, eta = eta, ...))

  ## here is our output from the double map
  raw_estimates <- map(lambda, function(lam)
    map(eta, function(et)
      fit_motl(X = X, Y = Y, w = w, eta = et, lambda = lam, family = family, ...)))

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
