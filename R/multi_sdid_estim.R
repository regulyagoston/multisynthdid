#' A function mapping a numeric vector to a (presumably sparser) numeric vector of the same shape to
#' be passed onto multi_synthdid_estim
#' @param v a vector
sparsify_function = function(v) { v[v <= max(v)/4] = 0; v/sum(v) }



#' Computes the synthetic diff-in-diff estimate for an average treatment effect on a treated block with multiple outcomes
#'
#' Method is an extension of 'Synthetic Difference in Differences' by Arkhangelsky et al, Algorithm 1.
#'
#' @param Y the observation tensor 3D (unit x time x outcome)
#' @param N0 the number of control units (N_co in the paper). Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps (T_pre in the paper). Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param X an optional 4-D array of time-varying covariates. Shape should be N X T X J x K standing for each j-outcome K covariates.
#' @param noise.level, an estimate of the noise standard deviation sigma. Defaults to the scaled standard deviation of first differences of Y. Scaling is based on Reguly et al (2024)
#' @param eta.omega  determines the tuning parameter zeta.omega = eta.omega * noise.level. Defaults to the value (N_tr T_post)^(1/4).
#' @param eta.lambda analogous for lambda.  Defaults to an 'infinitesimal' value 1e-6.
#' @param zeta.omega if passed, overrides the default zeta.omega = eta.omega * noise.level. Deprecated.
#' @param zeta.lambda analogous for lambda.
#' @param omega.intercept Binary. Use an intercept when estimating omega.
#' @param lambda.intercept Binary. Use an intercept when estimating lambda.
#' @param weights a list with fields lambda, omega and theta. If non-null weights$lambda is passed,
#'        we use them instead of estimating lambda weights. Same for weights$omega.
#'        Default for theta that weights different outcomes during the optimization is equal weights.
#' @param update.omega If true, solve for omega using the passed value of weights$omega only as an initialization.
#'        If false, use it exactly as passed. Defaults to false if a non-null value of weights$omega is passed.
#' @param update.lambda  Analogous.
#' @param min.decrease Tunes a stopping criterion for our weight estimator. Stop after an iteration results in a decrease
#' 		        in penalized MSE smaller than min.decrease^2.
#' @param max.iter A fallback stopping criterion for our weight estimator. Stop after this number of iterations.
#' @param sparsify A function mapping a numeric vector to a (presumably sparser) numeric vector of the same shape, which must sum to one.
#'                  If not null, we try to estimate sparse weights via a second round of Frank-Wolfe optimization
#'                  initialized at sparsify( the solution to the first round ).
#' @param max.iter.pre.sparsify Analogous to max.iter, but for the pre-sparsification first-round of optimization.
#'     		                Not used if sparsify=NULL.
#' @param standardize TRUE or FALSE. Default is TRUE if J > 1, otherwise FALSE.
#'                 If true outcomes are standardized (mean zero unit variance) during the estimation
#' @param scale_vcov TRUE or FALSE. If TRUE: scale during the optimization with specific variance-covariance matrix.
#' @param out Default = 'all'. Options are:
#'              - 'all': provides all possible estimates
#'              - 'tau_curve': only the treatment effects for each outcome after the event (gaps -- JxT1)
#'              - 'tau': only the point estimates after the event (point estimates Jx1)
#'              - 'tau_all': all the curves: synthetic, treated and treatment curves before and after the events (Jx(T0+T1))
#'
#' @return Specified in 'out'. If all, tau, the output is supplemented with 'weights' and 'setup' attached as attributes.
#'         'weights' contains the estimated weights lambda and omega and corresponding intercepts,
#'         as well as regression coefficients beta if X is passed.
#'         'setup' is a list describing the problem passed in: Y, N0, T0, X.
#' @export multi_synthdid_estim
multi_synthdid_estim <- function(Y, N0, T0,
                              X = NULL,
                              noise.level = NULL,
                              eta.omega = NULL,
                              eta.lambda = NULL,
                              zeta.omega  = NULL,
                              zeta.lambda = NULL,
                              omega.intercept = TRUE,
                              lambda.intercept = TRUE,
                              weights = list(omega = NULL, # unit weights
                                             lambda = NULL, # time weights
                                             theta = NULL ), #outcome weights
                              update.omega = is.null(weights$omega),
                              update.lambda = is.null(weights$lambda),
                              min.decrease = NULL,
                              max.iter = 1e4,
                              sparsify = sparsify_function,
                              max.iter.pre.sparsify = 100,
                              standardize = (dim(Y)[3]>1),
                              scale_vcov = T,
                              out = 'all' ) {

  stopifnot(nrow(Y) > N0, ncol(Y) > T0, is.list(weights),
            is.null(weights$lambda) || length(weights$lambda) == T0, is.null(weights$omega) || length(weights$omega) == N0,
            !is.null(weights$lambda) || update.lambda, !is.null(weights$omega) || update.omega ||
            is.logical(standardize) || is.logical(scale_vcov) )

  if ( is.na( dim(Y)[3] ) ){
    dim(Y) <- c(dim(Y)[1],dim(Y)[2],1)
  }
  J <- dim(Y)[3]

  if ( J > 1 && !is.null(X) && scale_vcov ){
    warning('Using explanatory variables X does not support var-covar scaling of Y.')
    scale_vcov = F
  }

  # Standardize Y and X
  if ( J != 1 && standardize ){
    mY = rep( NA, J )
    mS = rep( NA, J )
    for ( j in 1:J){
      mY[j] = mean( Y[,,j] )
      mS[j] = sd( Y[,,j] )
      # Standardize to avoid different gradient scales
      Y[,,j] <- ( Y[,,j] - mY[j] ) / mS[j]
    }
    # Standardize Xs
    if ( !is.null( X ) ){
      scaledX <- scale_Xs( X )
      X <- scaledX$X
    }
  }

  # Set theta weights
  if (is.null(weights$theta) ){
    weights$theta = rep( 1/J, J )
  }

  # Treated units
  N1 = nrow(Y) - N0
  # Post-treatment periods
  T1 = ncol(Y) - T0

  # Set the noise level
  if ( is.null( noise.level ) ){
    eta.omega <- (N1*T1)^(1/4)
    eta.lambda <- 1e-8
    if ( J == 1 ){ # Same as Arkhangelsky (2021)
      noise.level <- sd(apply(Y[1:N0,1:T0,1], 1, diff))
      t_s2_omega = noise.level
      t_s2_lambda = noise.level

    } else{ # Following Chava & Reguly (2024)
      t_s2_omega <- rep(NA,J)
      t_s2_lambda <- t_s2_omega
      # differencing across Ys and calculate the covariance
      for ( j in 1 : J ){
        # Regularizer for omega uses differences across time
        dY_j_om <- apply(Y[1:N0,1:T0,j], 1, diff)
        mdY_j_om <- mean( dY_j_om, na.rm = T )
        t_s2_omega[j] <- mean( ( dY_j_om - mdY_j_om )^2, na.rm = T )
        # Regularizer for lambda uses differences across units
        dY_j_lbd <- apply(Y[1:N0,1:T0,j], 2, diff)
        mdY_j_lbd <- mean( dY_j_lbd, na.rm = T )
        t_s2_lambda[j] <- mean( ( dY_j_lbd - mdY_j_lbd )^2, na.rm = T )
      }
      # Adjust with -1
      scale_N0T0 = N0*T0/(N0*T0-1)
      # Get the trace of the sqrt of sigma^2-s
      t_s2_omega = sqrt( scale_N0T0*t_s2_omega )
      t_s2_lambda = sqrt( t_s2_lambda*t_s2_omega )
      # Scale with number of outcomes
      eta.omega <- eta.omega
    }
    zeta.omega <- eta.omega  * t_s2_omega
    zeta.lambda <- eta.lambda * t_s2_lambda
    # take the minimum for stopping criterion
    min.decrease <- 1e-5 * min( c(t_s2_omega,t_s2_lambda), na.rm = T )
  }

  # Specify sparsify iteration
  if (is.null(sparsify)) {
    max.iter.pre.sparsify = max.iter
  }

  # If there are no covariates in the equation
  if ( is.null(X) ) {
    weights$vals = NULL
    weights$lambda.vals = NULL
    weights$omega.vals = NULL
    Yc   = collapsed.form(Y, N0, T0)
    Yc_n = Yc[1:N0,,]
    dim(Yc_n) <- c(N0,dim(Yc)[2:3])
    if (update.lambda) { # Time weights
      lambda.opt = sc.weight.fw(Yc_n,
                                zeta = zeta.lambda,
                                theta = weights$theta,
                                intercept = lambda.intercept,
                                weight=weights$lambda,
                                min.decrease = min.decrease,
                                max.iter = max.iter.pre.sparsify,
                                scale_vcov = scale_vcov,
                                smooth = F,
                                time_FE = F )
      if ( !is.null( sparsify ) ) { # Sparsify weights
        lambda.opt = sc.weight.fw(Yc_n,
                                  zeta = zeta.lambda,
                                  theta = weights$theta,
                                  intercept = lambda.intercept,
                                  weight=sparsify(lambda.opt$weight),
                                  min.decrease = min.decrease,
                                  max.iter = max.iter,
                                  scale_vcov = scale_vcov,
                                  smooth = F,
                                  time_FE = F )
      }
      weights$lambda      = lambda.opt$weight
      weights$lambda.vals = lambda.opt$vals
      weights$vals        = lambda.opt$vals
    }
    if (update.omega) { # Individual weights
      Yc_n_t = rep( NA, T0*dim(Yc)[1]*J )
      dim(Yc_n_t) <- c(T0,dim(Yc)[1],J)
      for (j in 1:J){Yc_n_t[,,j]=t(Yc[,1:T0,j])}
      omega.opt = sc.weight.fw(Yc_n_t,
                               zeta = zeta.omega,
                               theta = weights$theta,
                               intercept = omega.intercept,
                               weight=weights$omega,
                               min.decrease = min.decrease,
                               max.iter = max.iter.pre.sparsify,
                               scale_vcov = scale_vcov,
                               smooth = F,
                               time_FE = F )
      if(!is.null(sparsify)) {
        omega.opt = sc.weight.fw(Yc_n_t,
                                 zeta = zeta.omega,
                                 theta = weights$theta,
                                 intercept = omega.intercept,
                                 weight=sparsify(omega.opt$weight),
                                 min.decrease = min.decrease,
                                 max.iter = max.iter,
                                 scale_vcov = scale_vcov,
                                 smooth = F,
                                 time_FE = F )
      }
      weights$omega = omega.opt$weight
      weights$omega.vals = omega.opt$vals
      if (is.null(weights$vals)){ weights$vals = omega.opt$vals }
      else { weights$vals = pairwise.sum.decreasing(weights$vals, omega.opt$vals) }
    }
  } else {
    K = dim(X)[4]
    Yc = collapsed.form(Y, N0, T0)
    Xc = rep(NA,(N0+1)*(T0+1)*J*K)
    dim(Xc) <- c(dim(Yc),K)
    for ( k in 1:K ){
      X_k = X[,,,k]
      if ( J == 1){
        dim(X_k) <- c(dim(X_k),1)
      }
      Xc[,,,k] = collapsed.form(X_k,N0,T0)
    }
    # Check if beta is valid
    beta_valid = matrix(NA,nrow=J,ncol=K)
    for ( j in 1:J){
      for ( k in 1:K){
        beta_valid[j,k] = as.integer(!is.na(X[1,1,j,k]))
      }
    }

    weights = sc.weight.fw.covariates(Yc, weights$theta, Xc, beta_valid, zeta.lambda = zeta.lambda, zeta.omega = zeta.omega,
                                      lambda.intercept = lambda.intercept, omega.intercept = omega.intercept,
                                      min.decrease = min.decrease, max.iter = max.iter,
                                      lambda = weights$lambda, omega = weights$omega, update.lambda = update.lambda, update.omega = update.omega)
  }


  # Calculate Y (standardization w partialling out Xs)
  if ( J != 1 && standardize ){
    Yres = Y
    # Convert back Y
    for ( j in 1:J){
      Y[,,j] <- Y[,,j] * mS[j] + mY[j]
    }
    if ( is.null(X) ){
      for ( j in 1:J){
        # Convert back Y
        Yres[,,j] <- Yres[,,j] * mS[j] + mY[j]
      }
    } else{
      for ( j in 1:J){
        X_s = rep(NA,prod(dim(X)[c(1,2,4)]))
        dim(X_s) = c(dim(X)[c(1,2)],1,dim(X)[4])
        for ( k in 1:K ){
          X_s[,,1,k] =  X[,,j,k]*scaledX$mS_X[j,k]+scaledX$mX[j,k]
        }
        X.beta_j = contract4(X_s, matrix( weights$beta[j,], nrow = 1, ncol = K ), matrix( beta_valid[j,], nrow = 1, ncol = K ) )
        # Convert back Y
        Yres[,,j] = ( Yres[,,j] * mS[j] + mY[j] ) - X.beta_j[,,]
      }
    }
  } else if( !is.null( X ) ){
    X.beta = contract4(X, weights$beta, beta_valid)
    Yres = Y - X.beta
  } else{
    Yres = Y
  }

  fp1 = mm3d(t(c(-weights$omega, rep(1 / N1, N1))), Yres, tp = TRUE)
  estimate = rep( NA, J )
  for ( j in 1 : J ){ estimate[j] = fp1[,,j] %*% c(-weights$lambda, rep(1 / T1, T1)) }

  #estimate = t(c(-weights$omega, rep(1 / N1, N1))) %*% (Y - X.beta) %*% c(-weights$lambda, rep(1 / T1, T1))

  ## Add the curve to the estimation
  tau.curve = rep( NA, J * T1 )
  dim(tau.curve) = c(J,T1)
  for ( j in 1 : J ){
    tau.curve[j,] = fp1[,T0 + (1:T1),j] - c(fp1[,1:T0,j] %*% weights$lambda)
  }
  tau.curve

  class(estimate) = 'multi_synthdid_estimate'
  attr(estimate, 'estimator') = "multi_synthdid_estimate"
  attr(estimate, 'weights') = weights
  attr(estimate, 'setup') = list(Y = Y, N0 = N0, T0 = T0, X = X, Y_res = Yres )
  attr(estimate, 'opts') = list(zeta.omega = zeta.omega, zeta.lambda = zeta.lambda,
                                omega.intercept = omega.intercept, lambda.intercept = lambda.intercept,
                                update.omega = update.omega, update.lambda = update.lambda,
                                min.decrease = min.decrease, max.iter=max.iter,
                                noise.level = NULL )

  if ( out == 'all' ){
    return(list( tau = estimate, tau.curve = tau.curve ) )
  } else if ( out == 'tau_curve' ){
    return( tau.curve )
  } else if ( out == 'tau' ){
    return( estimate )
  } else if ( out == 'tau_all' ){
      return( multi_synthdid_curves( estimate, complete = T )$tau_curve )
  } else{
    stop('Wrong input for out use, all, tau_curve or tau.')
  }

}

#' multi_synthdid_estim for synthetic control estimates.
#' Takes all the same parameters, but by default, passes options to use the synthetic control estimator
#' By default, this uses only 'infinitesimal' ridge regularization when estimating the weights.
#' @param Y the observation tensor (N x T x J)
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param eta.omega determines the level of ridge regularization, zeta.omega = eta.omega * noise.level, as in multi_synthdid_estim
#' @param ... additional options for multi_synthdid_estim
#' @return an object like that returned by multi_synthdid_estim
#' @export multi_sc
multi_sc = function(Y, N0, T0, eta.omega = 1e-6, ...) {
  estimate = multi_synthdid_estim(Y, N0, T0, eta.omega = eta.omega,
                               weights = list(lambda = rep(0, T0)), omega.intercept = FALSE, ...)
  attr(estimate, 'estimator') = "multi_sc_estimate"
  estimate
}

#' multi_synthdid_estim for diff-in-diff estimates.
#' Takes all the same parameters, but by default, passes options to use the diff-in-diff estimator
#' @param Y the observation tensor (N x T x J)
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param ... additional  options for multi_synthdid_estim
#' @return an object like that returned by multi_synthdid_estim
#' @export did_estimate
did_estimate = function(Y, N0, T0, ...) {
  estimate = multi_synthdid_estim(Y, N0, T0, weights = list(lambda = rep(1 / T0, T0), omega = rep(1 / N0, N0)), ...)
  attr(estimate, 'estimator') = "multi_did_estimate"
  estimate
}

#' Computes a placebo variant of our estimator using pre-treatment data only
#' @param estimate, as output by multi_synthdid_estim
#' @param treated.fraction, the fraction of pre-treatment data to use as a placebo treatment period
#'        Defaults to NULL, which indicates that it should be the fraction of post-treatment to pre-treatment data
#' @export multi_synthdid_placebo
multi_synthdid_placebo = function(estimate, treated.fraction = NULL) {
  setup = attr(estimate, 'setup')
  opts = attr(estimate, 'opts')
  weights = attr(estimate, 'weights')
  X.beta = contract4(setup$X, weights$beta)
  estimator = attr(estimate, 'estimator')

  if (is.null(treated.fraction)) { treated.fraction = 1 - setup$T0 / ncol(setup$Y) }
  placebo.T0 = floor(setup$T0 * (1 - treated.fraction))

  if( is.null(setup$X) ){
    Xs = NULL
  } else {
    Xs = setup$X
  }

  do.call(estimator, c(list(Y=setup$Y[, 1:setup$T0,], N0=setup$N0, T0=placebo.T0, X=Xs), opts) )
}


#' Outputs the synthetic, treated units and the
#' effect curve that was averaged to produce our point estimate
#'  (at the moment it may give slightly different results
#' if standardization has been used during the optimization )
#' @param estimate, as output by multi_synthdid_estim
#' @param complete, TRUE/FALSE.
#'        - TRUE: it contains pre- and post event periods as well.
#'        - FALSE: it only contains after event values (Default)
#' @return list with:
#'        - synt_curve: synthetic outcome
#'        - tr_curve: treated (averaged) unit outcome
#'        - tau_curve: effect/gap curve
#' @export synthdid_effect_curve
multi_synthdid_curves = function(estimate, complete = F) {
  setup = attr(estimate, 'setup')
  weights = attr(estimate, 'weights')
  N1 = nrow(setup$Y) - setup$N0
  T1 = ncol(setup$Y) - setup$T0
  J = dim(setup$Y)[3]

  #if ( complete ){
  #  Y = setup$Y
  #} else{
    Y = setup$Y_res
  #}


  Y_m = Y[1:setup$N0,,]
  if ( J == 1 ){ dim(Y_m) = c(setup$N0,dim(Y)[2:3]) }
  synt_curve_0 = mm3d(t(weights$omega),Y_m, tp = TRUE)

  Y_t = Y[(setup$N0+1):(setup$N0+N1),,]
  if ( J == 1 || N1 == 1 ){ dim(Y_t) = c(N1,dim(Y)[2:3]) }
  if ( N1 == 1 ){
    tr_curve_0 = Y_t
  } else{
    tr_curve_0 = mm3d(t(rep(1 / N1, N1)),Y_t, tp = TRUE)
  }


  tau_cruve_0 = tr_curve_0 - synt_curve_0

  if ( complete ){
    synt_curve = rep( NA, J * ( T1 + setup$T0 ) )
    dim(synt_curve) = c(J,( T1 + setup$T0 ))
    tr_curve = rep( NA, J * ( T1 + setup$T0 ) )
    dim(tr_curve) = c(J,( T1 + setup$T0 ))
    tau_curve = rep( NA, J * ( T1 + setup$T0 ) )
    dim(tau_curve) = c(J,( T1 + setup$T0 ))
    for ( j in 1 : J ){
      tr_curve[j,] = tr_curve_0[,,j]
      tau_curve[j,] = tau_cruve_0[,,j] - c(tau_cruve_0[,1:setup$T0,j] %*% weights$lambda)
      synt_curve[j,] = tr_curve[j,] + tau_curve[j,]
    }
  } else{
    synt_curve = rep( NA, J * T1 )
    dim(synt_curve) = c(J,T1)
    tr_curve = rep( NA, J * T1 )
    dim(tr_curve) = c(J,T1)
    tau_curve = rep( NA, J * T1 )
    dim(tau_curve) = c(J,T1)
    for ( j in 1 : J ){
      synt_curve[j,] = synt_curve_0[,setup$T0 + (1:T1),j] - c(synt_curve_0[,1:setup$T0,j] %*% weights$lambda)
      tr_curve[j,] = tr_curve_0[,setup$T0 + (1:T1),j] - c(tr_curve_0[,1:setup$T0,j] %*% weights$lambda)
      tau_curve[j,] = tau_cruve_0[,setup$T0 + (1:T1),j] - c(tau_cruve_0[,1:setup$T0,j] %*% weights$lambda)
    }
  }

  return( list( synt_curve = synt_curve,
                tr_curve = tr_curve,
                tau_curve = tau_curve ) )
}
