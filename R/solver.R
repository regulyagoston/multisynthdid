#' contract a 4 dimensional object
contract4 = function(X, v, beta_valid ) {
  #stopifnot(length(dim(X)) == 3, dim(X)[3] == length(v))
  out = array(0, dim = dim(X)[1:3])
  #if (length(v) == 0) { return(out) }
  J = dim(v)[1]
  K = dim(v)[2]
  for (j in 1:J) {
    for (k in 1:K) {
      if ( beta_valid[j,k] == 1 ){ # if NA then just do not update
        out[,,j] = out[,,j] + v[j,k] * X[, ,j,k]
      }
    }
  }
  return(out)
}

#' 3D matrix multiplication
mm3d <- function( x, A, tp = FALSE ){
  J <- dim(A)[3]
  J2 <- dim(x)[3]
  if ( tp ){
    N <- dim(A)[2]
    y <- rep(NA,N*J)
    dim(y) <- c(1,N,J)
    if ( is.na( J2 ) ){
      for( j in 1:J){
        y[,,j] <- x %*% A[,,j]
      }
    } else{
        for( j in 1:J){
          y[,,j] <- x[,,j] %*% A[,,j]
        }
    }
  } else{
    N <- dim(A)[1]
    y <- rep(NA,N*J)
    dim(y) <- c(N,1,J)
    for( j in 1:J){
      y[,,j] <- A[,,j] %*% x
    }
  }
  return(y)
}


#' Updated to 3D arrays with weights given by theta
fw.step = function(A, x, b, theta, eta, alpha = NULL) {
  Ax = mm3d(x,A)
  err <- Ax - b
  eA <- mm3d(err,A,tp=TRUE)
  learning_step <- tcrossprod( x, eta )
  dim( learning_step ) <- c(1,dim(A)[2],dim(A)[3])
  half.grad <- ( eA + learning_step )
  half_grad_v <- half.grad
  dim(half_grad_v) <- c( dim(A)[2],dim(A)[3])
  i = which.min( half_grad_v %*% theta )
  if (!is.null(alpha)) {
    x = x * (1 - alpha)
    x[i] = x[i] + alpha
    return(x)
  } else {
    d.x = -x
    d.x[i] = 1 - x[i]
    if (all(d.x == 0)) { return(x) }
    d.err = ( A[, i, ] - Ax[,,] )
    step_1 <- -mm3d(d.x,half.grad)
    if (dim(A)[3]==1){
      step_2 <- sum(d.err^2)  + eta * sum(d.x^2)
    } else{
      step_2 <- colSums(d.err^2)  + eta * sum(d.x^2)
    }
    step = c(step_1 / step_2) %*% theta
    constrained.step = min(1, max(0, step))
    return( x + constrained.step * d.x )
  }
}

#' A Frank-Wolfe solver for synthetic control weights using exact line search
sc.weight.fw = function(Y, zeta, theta,
                        intercept = TRUE, weight = NULL,
                        min.decrease = 1e-3, max.iter = 1000,
                        scale_vcov = T, smooth = F, time_FE = F ) {
  dimY <- dim(Y)
  T0 = dimY[2] - 1
  N0 = dimY[1]
  J <- dimY[3]
  if ( smooth ){
    stop('Does not work!')
  }
  if ( is.null(weight) ) {
    if ( !smooth ){
      weight = rep(1 / T0, T0)
    } else
      # Infinity means weight all equally
      weight  = 100
  }

  # De-meaning
  if (intercept) {
    if ( !is.na( J ) ){
      for (j in 1:J){Y[,,j] = apply(Y[,,j], 2, function(col) { col - mean(col) }) }
    } else{
      Y <- apply(Y, 2, function(col) { col - mean(col) })
    }
  }

  # In case of smoothing one may want to take time fixed effects as well
  if (smooth && time_FE) {
    if ( !is.na( J ) ){
      for (j in 1:J){Y[,,j] = apply(Y[,,j], 1, function(col) { col - mean(col) }) }
    } else{
      Y <- apply(Y, 1, function(col) { col - mean(col) })
    }
  }

  # Scaling with variance-covariance matrices
  if ( scale_vcov && J > 1){
    mS = array(NA,dim=c(J,J,N0))
    for ( i in 1:N0 ){
      # Covariance
      vcov_time = cov(Y[i,,])
      # take to the power of -1/2
      V <- eigen(vcov_time)$vectors
      L <- eigen(vcov_time)$values
      i2_vcov <- V %*% diag(1/sqrt(L)) %*% t(V)
      # Save vcov-scaler
      mS[,,i] <- i2_vcov
    }
  } else{
    scale_vcov = F
  }

  ###
  # Initialize the iteration process
  t = 0
  vals = rep(NA, max.iter)
  # controls before
  A = Y[, 1:T0, ]
  dim(A) <- c(dim(Y)[1],T0,dim(Y)[3])
  # event-point
  b = Y[, T0 + 1, ]
  dim(b) <- c(dim(Y)[1],1,dim(Y)[3])
  # Learning rate
  eta = N0 * Re(zeta^2)
  # Controls
  YN0 <- Y[1:N0, ,]
  dim(YN0) <- c(N0,dim(Y)[2],dim(Y)[3])

  # Iteration
  while (t < max.iter && (t < 2 || vals[t - 1] - vals[t] > min.decrease^2) ) {
    weight_last = weight
    t = t + 1
    weight.p = fw.step(A, weight, b, theta, eta, alpha = NULL)
    if ( any(is.na(weight.p)) ){
       vals[t] = vals[t-1]
       t = max.iter
       break
    }
    weight = weight.p
    err = mm3d( c(weight, -1), YN0 )
    dim(err) <- c( N0,J )
    # Penalization works for both
    if ( scale_vcov ){
      pen_err = 0
      for ( i in 1:N0 ){
        pen_err = pen_err + t(err[i,]) %*% mS[,,i] %*% err[i,]
      }
      vals[t] = ( Re(zeta^2) * sum(weight^2) )%*% theta + pen_err / (N0*J)

    } else{
      vals[t] = ( Re(zeta^2) * sum(weight^2) + colSums(err^2) / N0 ) %*% theta
    }
  }
  # If newer optimization step is worse use the previous values
  if ( vals[t - 1] - vals[t] <= 0 ){
    weight_out = weight_last
  } else{
    weight_out = weight
  }
  list(weight = weight_out,
       vals = vals,
       convergence = t < max.iter )
}

#' A Frank-Wolfe + Gradient solver for lambda, omega, and beta when there are covariates
#' Uses the exact line search Frank-Wolfe steps for lambda, omega and (1/t)*gradient steps for beta
#' pass update.lambda=FALSE/update.omega=FALSE to fix those weights at initial values, defaulting to uniform 1/T0 and 1/N0
sc.weight.fw.covariates = function(Y, theta, X, beta_valid, zeta.lambda = NULL, zeta.omega = NULL,
                                   lambda.intercept = TRUE, omega.intercept = TRUE,
                                   min.decrease = 1e-3, max.iter = 1000,
                                   lambda = NULL, omega = NULL, beta = NULL, update.lambda = TRUE, update.omega = TRUE ) {


  T0 = ncol(Y) - 1
  N0 = nrow(Y) - 1
  J = dim(Y)[3]
  K = dim(X)[4]
  #if (length(dim(X)) == 2) { dim(X) = c(dim(X), 1) }
  if (is.null(lambda)) {  lambda = rep(1 / T0, T0)   }
  if (is.null(omega)) {  omega = rep(1 / N0, N0)    }
  if (is.null(beta)) {  beta = rep(0, J*K )
                        dim(beta) = c(J,K) }
  if (is.null(zeta.lambda)) {  zeta.lambda = rep(0, J)   }
  if (is.null(zeta.omega)) {  zeta.omega = rep(0, j)   }

  # Optimization function with Xs
  update.weights = function(Y, lambda, omega ) {
    ######
    ## Update Lambda
    # De-mean Y -- observation wise
    Yc_n <- Y[1:N0,,]
    dim(Yc_n) <- c(N0,dim(Y)[2],J)
    Y.lambda = Yc_n
    if (lambda.intercept) {
      if ( !is.na( J ) ){
        for (j in 1:J){Y.lambda[,,j] = apply(Y.lambda[,,j], 2, function(row) { row - mean(row) }) }
      } else{
        Y.lambda <- apply(Y[1:N0,], 2, function(col) { col - mean(col) })
      }
    }
    # forward-step for lambda
    if (update.lambda) {
      Yc_nt0 <- Y.lambda[, 1:T0,]
      dim(Yc_nt0) <- c(N0,T0,J)
      Yc_nt1 <- Y.lambda[, T0+1,]
      dim(Yc_nt1) <- c(N0,1,J)
      lambda = fw.step(Yc_nt0, lambda, Yc_nt1, theta, N0 * Re(zeta.lambda^2))
    }
    # Calculate lambda error
    err.lambda = mm3d(c(lambda, -1),Y.lambda)

    ######
    ## Update omega
    # De-mean Y -- time wise
    Yc_t = rep( NA, T0*dim(Y)[1]*J )
    dim(Yc_t) <- c(T0,dim(Y)[1],J)
    for (j in 1:J){Yc_t[,,j]=t(Y[,1:T0,j])}
    Y.omega = Yc_t
    if (omega.intercept) {
      if ( !is.na( J ) ){
        for (j in 1:J){Y.omega[,,j] = apply(Y.omega[,,j], 2, function(row) { row - mean(row) }) }
      } else{
        Y.omega = apply(t(Y[, 1:T0]), 2, function(row) { row - mean(row) })
      }
    }

    # forward-step for omega
    if (update.omega) {
      Yc_tn0 <- Y.omega[, 1:N0,]
      dim(Yc_tn0) <- c(T0,N0,J)
      Yc_tn1 <- Y.omega[, N0 + 1,]
      dim(Yc_tn1) <- c(T0,1,J)
      omega = fw.step(Yc_tn0, omega, Yc_tn1 , theta, T0 * Re(zeta.omega^2))
    }
    # Calculate error
    err.omega = mm3d(c(omega, -1),Y.omega)

    val = ( Re(zeta.omega^2) * sum(omega^2) +
            Re(zeta.lambda^2) * sum(lambda^2) +
            colSums(err.omega^2) / T0 +
            colSums(err.lambda^2) / N0 ) %*% theta
    list(val = val, lambda = lambda, omega = omega, theta = theta, err.lambda = err.lambda, err.omega = err.omega)
  }

  vals = rep(NA, max.iter)
  t = 0
  Y.beta = Y - contract4(X, beta, beta_valid)
  weights = update.weights(Y.beta, lambda, omega )
  grad.beta = rep(0, J*K )
  dim(grad.beta) = c(J,K)
  # state is kept in weights$lambda, weights$omega, beta
  while (t < max.iter && (t < 2 || vals[t - 1] - vals[t] > min.decrease^2)) {
    t = t + 1
      for ( j in 1 : J ){
        Xjk = X[,,j,]
        dim(Xjk) <- c(N0+1,T0+1,K)
        grad.beta[j,] = -apply(Xjk, 3, function(Xi) {
          t(weights$err.lambda[,,j]) %*% Xi[1:N0, ] %*% c(weights$lambda, -1) / N0 +
            t(weights$err.omega[,,j]) %*% t(Xi[, 1:T0]) %*% c(weights$omega, -1) / T0
        })
      }
    grad.beta = grad.beta * beta_valid


    alpha = 1 / t
    beta = beta - alpha * grad.beta
    Y.beta = Y - contract4(X, beta, beta_valid)
    weights = update.weights(Y.beta, weights$lambda, weights$omega )
    vals[t] = weights$val
  }
  list(lambda = weights$lambda, omega = weights$omega, beta = beta, theta = weights$theta, vals = vals)
}


#' Standardize explanatory covariates X
scale_Xs <- function( X ){
  J = dim(X)[3]
  K = dim(X)[4]
  mX = matrix( NA, nrow = J, ncol = K)
  mS_X = matrix( NA, nrow = J, ncol = K)
  for ( j in 1:J){
    for ( k in 1:K ){
      if (!all(is.na(c(X[,,j,k])))){
        mX[j,k] = mean( X[,,j,k] )
        mS_X[j,k] = sd( X[,,j,k] )
        # Standardize
        X[,,j,k] <- ( X[,,j,k] - mX[j,k] ) / mS_X[j,k]
      }
    }
  }
  return(X=X,mX = mX, mS_X = mS_X)
}




