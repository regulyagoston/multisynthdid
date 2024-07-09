#' Convert classical R regression input to multi synthetic diff-in-diff estimator
#'
#'
#' @param formula_yx classical formula used in R
#'    - multiple outcomes: y1 + y2 + y3 ~ 1
#'    - classical synthetit diff-in-diff: y1 ~ 1
#'    - multiple outcomes with SAME covariates: y1 + y2 + y3 ~ a + b + c
#'    - multiple outcomes with DIFFERENT covariates: cbind(y1,y2,y3)~cbind(a+b,a+c,a+b+c)
#' @param treatment character value for the treatment (logical) variable
#' @param unit character value for the unit identifier variable
#' @param time character value for the time identifier variable
#' @param data data.frame or tibble object
#' @param SE TRUE/FALSE if compute standard errors or not
#' @param se_type type of standard errors
#'    - 'placebo' standard errors (Default)
#'    - 'bootstrap' standard errors
#' @param se_replication positive integer for
#'    the number of replication for SE calculation (default = 200).
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

#'
#'
multi_sdid <- function( formula_yx, treatment, unit, time, data,
                        SE = T, se_type = 'placebo', se_replication = 200,
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
                        standardize = NULL,
                        scale_vcov = T ){

  stopifnot( is.logical(SE) | ( se_type %in% c('placebo','bootstrap') ) |
              ( is.integer( se_replication ) & se_replication > 0 ) )

  # Outcomes
  formula_yx <- as.formula(formula_yx)
  outcome <- all.vars(formula_yx[[2]])
  J = length(outcome)
  # Covariates
  covariates = all.vars(formula_yx[[3]])

  # Check for multivariate equation definition with `cbind`
  if ( length( covariates) > 0 && str_sub(as.character(formula_yx[2]),0,5) == 'cbind'){
    if (str_sub(as.character(formula_yx[3]),0,5) == 'cbind'){
      eq_s = matrix('',nrow=J,ncol=length(covariates))
      eq_0 <- as.character( formula_yx[[3]] )
      eq_1 <- eq_0[2:length(eq_0)]
      for ( j in 1:J ){
        eq_j <- strsplit(unlist(eq_1[j]),' \\+ ')[[1]]
        eq_s[j, 1:length(eq_j)] <- eq_j
      }
    } else{ # Unified covariates -- uses everywhere the same covariates
      eq_s = NULL
    }
  } else{
    covariates = NULL
    eq_s = NULL
  }

  # Create the necessary matrices (tensors) for estimation
  setup0 <- setup_multi_sdid(data, unit = unit, time = time,
                             outcome = outcome, treatment = treatment,
                             covariates = covariates, eq_s = eq_s )


  # Estimate multi_synthetic diff-in-diff
  m_sdid <- multi_synthdid_estim(setup0$Y, setup0$N0, setup0$T0,
                                 X = setup0$X,
                                 noise.level = noise.level,
                                 eta.omega = eta.omega,
                                 eta.lambda = eta.lambda,
                                 zeta.omega  = zeta.omega,
                                 zeta.lambda = zeta.lambda,
                                 omega.intercept = omega.intercept,
                                 lambda.intercept = lambda.intercept,
                                 weights = weights,
                                 update.omega = update.omega,
                                 update.lambda = update.lambda,
                                 min.decrease = min.decrease,
                                 max.iter = max.iter,
                                 sparsify = sparsify_function,
                                 max.iter.pre.sparsify = max.iter.pre.sparsify,
                                 standardize = ifelse( is.null( standardize ), J > 1, standardize ),
                                 scale_vcov = scale_vcov ,
                                 out = 'all' )

  # Restructure the data for user-friendly format
  tau_curves = t(m_sdid$tau.curve)
  colnames( tau_curves ) <- outcome
  tau_df <- as.data.frame( tau_curves )
  time_var <- colnames(setup0$W)
  if ( inherits( data[[time]], 'double' ) ){
    time_var = as.double(time_var)
  } else if ( inherits( data[[time]], 'Date' ) ){
    time_var = as.Date(time_var)
  }
  tau_df$year <- time_var[(setup0$T0 + 1): length(time_var)]

  # Final estimates
  tau <- colMeans(tau_curves)
  # Assign SE-s
  if ( SE ){
    se_tau = se.multi_sdid(m_sdid$tau, method = se_type,
       replications = se_replication )
    # Save for later potential use
    m_sdid$se_all <- se_tau
    # Point estimates
    se_tau_point <- t(as.matrix(se_tau$tau_se))
    colnames( se_tau_point ) <- outcome
    tau <- rbind(tau, se_tau_point)
    rownames( tau ) <- c('Avg. Treatment','SE')
    # Curve estimates
    se_tau_curve <- t( se_tau$tau_curve_se )
    colnames( se_tau_curve ) <- paste0( outcome, '_SE' )
    tau_df <- cbind( tau_df, se_tau_curve[(setup0$T0+1):nrow(se_tau_curve), ] )
  }

  m_sdid$tau_all <- tau
  m_sdid$tau_df <- tau_df
  m_sdid$time <- time_var
  m_sdid$unit <- setup0$unit

  # Weight tables
  weights <- attr(m_sdid$tau,'weights')
  unit_weights <- data.frame( weights$omega, setup0$unit )
  colnames(unit_weights) <- c('weight',unit)
  time_weights <- data.frame( weights$lambda, time_var[1:setup0$T0] )
  colnames(time_weights) <- c('weight',time)
  m_sdid$weights <- list( unit = unit_weights, time = time_weights )

  class(m_sdid) = 'multi_synthdid_obj'
  return( m_sdid )

}
