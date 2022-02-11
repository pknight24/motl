#' Simulate coefficient vectors w and beta according to a multivariate normal distribution
#' @param p Integer, the length of the vectors.
#' @param var_vector A vector of length 2, contains the variance of beta and w (in order).
#' @param rho The correlation between w and beta.
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm rbinom
#' @export
simulate_coef <- function(p, var_vector = c(0.1, 0.1), rho = 0.75)
{
  covariance <- diag(2) %*% diag(var_vector)
  covariance[1,2] <- prod(sqrt(var_vector)) * rho
  covariance[2,1] <- prod(sqrt(var_vector)) * rho
  output <- MASS::mvrnorm(n = p, mu = c(0,0), Sigma = covariance)
  return(list(w = output[,1], beta = output[,2]))
}

#' Simulate data under a linear model
#' @param n Integer, number of samples.
#' @param p Integer, number of features.
#' @param beta A true coefficient vector of length p.
#' @param sigma2 Variance of the error term for Gaussian data.
#' @param family Distribution of the outcome data, defaults to 'gaussian' for
#' continuous data. For binary outcomes, specify 'binomial'.
#' @export
simulate_data <- function(n, p, beta, sigma2 = 0.1, family = "gaussian"){
  X <- matrix(rnorm(n = n*p), nrow = n, ncol = p)
  if (family == "gaussian")
  {
    Y <- X %*% beta + rnorm(n = n, mean = 0, sd = sqrt(sigma2))
  }
  else if (family == "binomial")
  {
    expit <- function(x) exp(x) / (1 + exp(x))
    probs <- apply(X, 1, function(x) expit(crossprod(x, beta)[1,1]))
    Y <- rbinom(n = n, size = 1, prob = probs)
  }
  list(X = X, Y = Y)
}
