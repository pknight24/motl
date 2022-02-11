fit_motl <- function(X, Y, w,lambda = 1, eta = 1, family = "gaussian", rho = 1, n_iter = 5)
{
  n <- nrow(X)
  p <- ncol(X)

  if (family == "gaussian")
  {
    ### initialize the variables; what is the best way to do this?
    ### current approach is to compute the initial b via a standard ridge, and generate u randomly
    b <- coef(glmnet::glmnet(x = X, y = Y, lambda = lambda, intercept = FALSE,
                             standardize = FALSE, alpha = 0))[-1]
    z <- b
    u <- rnorm(p)
    for (i in 1:n_iter)
    {
      ###### stage 1: compute b_k+1
      c <- z - u
      Y_tilde <- Y - X %*% c
      b <- coef(glmnet::glmnet(x = X, y = Y_tilde, lambda = rho/2, standardize = FALSE, intercept = FALSE,
                       alpha = 0))[-1] + c
      ### stage 2 + 3: update other variables
      z <- (2 * eta / (2 * lambda + rho)) * w + (rho / (2 * lambda + rho)) * (b + u)
      u <- u + b - z
    }
    return(b)
  }
  else if (family == "binomial")
  {
    expit <- function(x)exp(x) / (1 + exp(x))
    loss <- function(beta) 
    {
      -1 * sum(sapply(1:n, function(i){
        Y[i] * log(expit(crossprod(X[i,], beta)[1,1])) + 
          (1 - Y[i]) * log(1 - expit(crossprod(X[i,], beta)[1,1]))
      })) + lambda * crossprod(beta) - 2 * eta * t(w) %*% beta
    }
    gradient <- function(beta)
    {
      sum(sapply(1:n, function(i){
        (Y[i] - expit(crossprod(X[i,], beta)[1,1])) * X[i,]
      })) + lambda * beta - 2 * eta * w

    }

    b <- coef(glmnet(x = X, y = Y, family = "binomial", intercept = FALSE,
      lambda = lambda, alpha = 0))[-1]
    optim(par = b, fn = loss, gr = gradient, method = "BFGS")$par

  }

}
