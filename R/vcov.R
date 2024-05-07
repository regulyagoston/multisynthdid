#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#'
#' Provides variance estimates based on the following two options following Arkangelsky et al.
#' \itemize{
#'   \item The bootstrap, Algorithm 2 in Arkhangelsky et al.
#'   \item Placebo, Algorithm 4 in Arkhangelsky et al.
#'   \item Jacknife is not applicable
#' }
#'
#' "placebo" is the only option that works for only one treated unit.
#'
#' @param object A synthdid model
#' @param method, the CI method. The default is placebo (warning: this may be slow on large
#'  data sets).
#' @param replications, the number of bootstrap replications
#'
#' @references Agoston Reguly, Chava Sudheer.
#'  "The Long-Run Stock Market Performance of Mergers and Acquisitions". SSRN prepint
#'  https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4742401
#'
#' @references Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager.
#'  "Synthetic Difference in Differences". arXiv preprint arXiv:1812.09970, 2019.
#'
#' @method se multi_sdid
#' @export
se.multi_sdid = function(object,
                                  method = c("bootstrap", "placebo"),
                                  replications = 200 ) {
  method = match.arg(method)
  if(method == 'bootstrap') {
    stop('Not implemented')
    se = bootstrap_se(object, replications)
  } else if(method == 'jackknife') {
    stop('Not implemented')
  } else if(method == 'placebo') {
    se = placebo_se(object, replications)
  }
}

#' The bootstrap se, based on Algorithm 2 of Arkhangelsky et al.
#' @param estimate multi_sdid object
#' @param replications number of replicates
bootstrap_se = function(estimate, replications) {
  setup = attr(estimate, 'setup')
  opts = attr(estimate, 'opts')
  opts$out = 'tau_all'
  weights = attr(estimate, 'weights')
  if (setup$N0 == nrow(setup$Y) - 1) { return(NA) }
  theta = function(ind) {
    if(all(ind <= setup$N0) || all(ind > setup$N0)) { NA }
    else {
      weights.boot = weights
      weights.boot$omega = sum_normalize(weights$omega[sort(ind[ind <= setup$N0])])
      if( is.null(setup$X) ){
        Xs = NULL
      } else {
        Xs = setup$X[ind, , ,]
        dim(Xs) <- c(length(ind),dim(setup$X)[2:4])
      }
      do.call(multi_synthdid_estim, c(list(Y=setup$Y[sort(ind),,], N0=sum(ind <= setup$N0), T0=setup$T0, X=Xs, weights=weights.boot), opts))
    }
  }
  J = dim(setup$Y)[3]
  T_all = dim(setup$Y)[2]
  bootstrap.estimates = array(NA, c(J,T_all,replications) )
  count = 0
  while(count < replications) {
    bootstrap.estimates[,,count+1] = theta(sample(1:nrow(setup$Y), replace=TRUE))
  }

  # Overall
  scaler = sqrt((replications-1)/replications)
  tau_se = rep( NA, J )
  for ( j in 1 : J ){
    rep_output_j = bootstrap.estimates[j, ( setup$T0+1 ) : T_all, ]
    tau_se[j] = scaler*sd( colMeans( rep_output_j, na.rm = T ), na.rm = T )
  }

  # For each point in the curve
  tau_curve_se = matrix( NA, nrow = J, ncol = T_all )
  # Get the SE for each point in the curve
  for ( j in 1:J){
    for ( t in 1:T_all){
      tau_curve_se[j,t] = sd( bootstrap.estimates[j,t,], na.rm = T)
    }
  }
  tau_curve_se = scaler * tau_curve_se

  return( list( tau_hat = bootstrap.estimates, tau_se = tau_se, tau_curve_se = tau_curve_se ) )
}


#' The placebo se, based on Algorithm 4 of Arkhangelsky et al.
#' @param estimate multi_sdid object
#' @param replications number of replicates
placebo_se = function(estimate, replications) {
  setup = attr(estimate, 'setup')
  opts = attr(estimate, 'opts')
  opts$out = 'tau_all'
  weights = attr(estimate, 'weights')
  N1 = nrow(setup$Y) - setup$N0
  if (setup$N0 <= N1) { stop('must have more controls than treated units to use the placebo se') }
  theta = function(ind) {
    N0 = length(ind)-N1
    weights.boot = weights
    weights.boot$omega = sum_normalize(weights$omega[ind[1:N0]])
    Y_rnd = setup$Y[ind,,]
    if( is.null(setup$X) ){
      Xs = NULL
    } else {
      Xs = setup$X[ind, , ,]
      dim(Xs) <- c(length(ind),dim(setup$X)[2:4])
    }
    do.call(multi_synthdid_estim, c(list(Y=Y_rnd, N0=N0,  T0=setup$T0, X=Xs, weights=weights.boot), opts))
  }
  rep_output = replicate(replications, theta(sample(1:setup$N0)))
  scaler = sqrt((replications-1)/replications)
  J = dim( rep_output )[1]
  # Overall effect
  T1 = dim(rep_output)[2]
  tau_se = rep( NA, J )
  for ( j in 1 : J ){
    rep_output_j = rep_output[j,(setup$T0+1):T1,]
    tau_se[j] = scaler * sd( colMeans( rep_output_j ) )
  }

  # For each point in the curve
  tau_curve_se = matrix( NA, nrow = J, ncol = T1 )
  # Get the SE for each point in the curve
  for ( j in 1:J){
    for ( t in 1:T1){
      tau_curve_se[j,t] = sd(rep_output[j,t,])
    }
  }
  tau_curve_se = scaler * tau_curve_se

  return( list( tau_hat = rep_output, tau_se = tau_se, tau_curve_se = tau_curve_se ) )
}

#' Auxilary function for bootstrapping
sum_normalize = function(x) {
  if(sum(x) != 0) { x / sum(x) }
  else { rep(1/length(x), length(x)) }
  # if given a vector of zeros, return uniform weights
  # this fine when used in bootstrap and placebo standard errors.
}
