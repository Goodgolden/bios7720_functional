## Function for doing smoothing parameter selection by
## looping over a grid of candidate (log) smoothing parameters.
## Note that for simplicity of presentation I have skipped normal
## input checking.
## Inputs:
## y = vector of response
## Phi = Full rank design matrix (spline basis) N x K
## S = penalty matrix. Must be K x K
## S_sqrt = S^0.5, not required, but if not specified S must be positive definite
##          or else the cholesky decomposition will fail
## loglamba = vector of log smoothing parameters to search over
## var = logical indicating whether to calculate variance of the coefficients
## tol = tolerance to pass through to R's solve() function
##       we specify the default here as lower than R's default as some of the
##       matrices we will invert are nearly singular
## Outputs:
## xi = vector of estimated coefficients conditional on the optimal smoothing parameter
## var_xi = variance/covariance matrix of estimated coefficients
## sig2e = estimated residual variance
## sp = optimal smoothing parameter
## gcv = vector of GCV values associated with the loglambda input vector
fn_PENSSE_grid <- function(y, Phi, S, S_sqrt = NULL, loglambda,
                           var = TRUE, tol = 1e-30) {
  N <- length(y)
  ## pre-calculate some quantities of interest
  if (is.null(S_sqrt)) {
    S_sqrt <- chol(S) # S^(1/2)
  }
  Phi_sq <- crossprod(Phi) # Phi^t Phi
  # augmented data / design vector/matrix quantities
  K <- ncol(Phi)
  a <- rep(0, K)
  ya <- c(y, a)

  ## if provided multiple lambdas, use gcv to find the best
  nlambda <- length(loglambda)
  lambda <- exp(loglambda)
  if (nlambda > 1) {
    gcv_vec <- rep(NA, nlambda)
    for (l in 1:nlambda) {
      ## get trace(H)
      trH <- sum(diag(Phi_sq %*% solve(Phi_sq + lambda[l] * S, tol = tol)))
      ## get \hat{y}
      B <- sqrt(lambda[l]) * S_sqrt
      PhiB <- rbind(Phi, B)
      PhiB_ya <- crossprod(PhiB, ya)
      PhiB_sq <- crossprod(PhiB)
      xi_hat <- qr.solve(a = PhiB_sq, b = PhiB_ya, tol = tol)
      y_hat <- Phi %*% xi_hat
      ## get GCV criteria
      gcv_vec[l] <- N * sum((y_hat - y)^2 / (N - trH)^2)
    }
    lambda <- lambda[which.min(gcv_vec)]
  }
  ## get the final estimates using either the one lambda supplied
  ## or the "best" selected via GCV
  B <- sqrt(lambda) * S_sqrt
  PhiB <- rbind(Phi, B)
  PhiB_ya <- crossprod(PhiB, ya)
  PhiB_sq_p <- crossprod(PhiB)
  xi_hat <- qr.solve(a = PhiB_sq_p, b = PhiB_ya, tol = tol)
  y_hat <- Phi %*% xi_hat
  ## get variance estimates
  if (var) {
    Phi_sq_p_inv <- solve(Phi_sq + lambda * S, tol = tol)
    trH <- sum(diag(Phi_sq %*% Phi_sq_p_inv))
    sig2e <- sum((y_hat - y)^2) / (N - trH)
    var_xi_hat <- sig2e * Phi_sq_p_inv %*% Phi_sq %*% Phi_sq_p_inv
  } else {
    var_xi_hat <- NULL
  }

  ret <- list(
    xi = xi_hat,
    var_xi = var_xi_hat,
    sig2e = sig2e,
    sp = lambda,
    gcv = gcv_vec
  )
}


## function for selecting smoothing parameter via GCV
## chosen using UNSTABLE numerical optimization
## Inputs:
## loglambda = input parameter vector (in this case, a scalar)
## y = response vector
## ya = augmented response vector
## Phi = spline basis
## Phi_sq = Phi^t Phi
## S = penalty matrix
## S_sqrt = S^0.5
## N = number of data points, should be equal to length(y)
opt_gcv <- function(loglambda, y, Phi, ya, Phi_sq, S, S_sqrt, N, tol = 1e-30) {
  lambda <- exp(loglambda)
  ## get trace(H)
  trH <- sum(diag(Phi_sq %*% solve(Phi_sq + lambda * S, tol = tol)))
  ## get \hat{y} using augmented design matrix
  B <- sqrt(lambda) * S_sqrt
  PhiB <- rbind(Phi, B)
  PhiB_ya <- crossprod(PhiB, ya)
  PhiB_sq <- crossprod(PhiB)
  xi_hat <- qr.solve(a = PhiB_sq, b = PhiB_ya, tol = tol)
  y_hat <- Phi %*% xi_hat
  ## get GCV criteria
  N * sum((y_hat - y)^2 / (N - trH)^2)
}

## Function for doing smoothing parameter selection by (constrained) numerical optimization.
## Modifies fn_PENSSE_grid to choose smoothing paramter via the opt_gcv function.
## Note that for simplicity of presentation I have skipped normal
## input checking.
## Inputs:
## y = vector of response
## Phi = Full rank design matrix (spline basis) N x K
## S = penalty matrix. Must be K x K
## S_sqrt = S^0.5, not required, but if not specified S must be positive definite
##          or else the cholesky decomposition will fail
## loglamba = vector of log smoothing parameters to search over
## var = logical indicating whether to calculate variance of the coefficients
## tol = tolerance to pass through to R's solve() function
##       we specify the default here as lower than R's default as some of the
##       matrices we will invert are nearly singular
## lower = lower bound for log smoothing parameters to search
## upper = upper bound for log smoothing parameters to search
## loglambda_st = starting value for the optimization routine for log(lambda)
## Outputs:
## xi = vector of estimated coefficients conditional on the optimal smoothing parameter
## var_xi = variance/covariance matrix of estimated coefficients
## sig2e = estimated residual variance
## sp = optimal smoothing parameter
## gcv = scalar value of GCV associated with the optimal smoothing parameter
fn_PENSSE_opt <- function(y, Phi, S, S_sqrt = NULL, var = TRUE,
                          tol = 1e-30, lower = -20, upper = 20, loglambda_st = 5) {
  N <- length(y)
  ## pre-calculate some quantities of interest
  if (is.null(S_sqrt)) {
    S_sqrt <- chol(S) # S^(1/2)
  }
  Phi_sq <- crossprod(Phi) # Phi^t Phi
  Phi_y <- t(Phi) %*% y # Phi^t y
  # augmented data / design vector/matrix quantities
  K <- ncol(Phi)
  a <- rep(0, K)
  ya <- c(y, a)
  # do the optimization
  opt <- optim(loglambda_st,
    opt_gcv,
    y = y,
    Phi = Phi,
    ya = ya,
    Phi_sq = Phi_sq,
    S = S,
    S_sqrt = S_sqrt,
    N = N,
    method = "L-BFGS-B",
    control = list(maxit = 10000),
    lower = lower,
    upper = upper
  )
  lambda <- exp(opt$par)
  ## get the final estimates using
  ## either the one lambda supplied
  ## or the "best" selected via GCV
  B <- sqrt(lambda) * S_sqrt
  PhiB <- rbind(Phi, B)
  PhiB_ya <- crossprod(PhiB, ya)
  PhiB_sq_p <- crossprod(PhiB)
  xi_hat <- qr.solve(
    a = PhiB_sq_p,
    b = PhiB_ya,
    tol = tol
  )
  y_hat <- Phi %*% xi_hat
  ## get variance estimates
  if (var) {
    Phi_sq_p_inv <- solve(Phi_sq + lambda * S,
      tol = tol
    )
    trH <- sum(diag(Phi_sq %*% Phi_sq_p_inv))
    sig2e <- sum((y_hat - y)^2) / (N - trH)
    var_xi_hat <- sig2e * Phi_sq_p_inv %*% Phi_sq %*% Phi_sq_p_inv
  } else {
    var_xi_hat <- NULL
  }
  ## return values
  ret <- list(
    xi = xi_hat,
    var_xi = var_xi_hat,
    sig2e = sig2e,
    sp = lambda,
    gcv = opt$value
  )
}
