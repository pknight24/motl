fit_motl <- function(X, Y, w,lambda = 1, eta = 1, rho = 1, n_iter = 5)
{
  n <- nrow(X)
  p <- ncol(X)

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
