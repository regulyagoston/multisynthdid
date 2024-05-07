#' Collapse Y to an N0+1 x T0+1 x J vector by averaging the last N1=nrow(Y)-N0 rows and T1=ncol(Y)-T0 columns
#' it is the modified version of Arkhangelsky et al for
collapsed.form = function(Y, N0, T0) {
  dimY <- dim(Y)
  N = dimY[1]
  Ta = dimY[2]
  J = dimY[3]
  Y_cf <- rep(NA,(N0+1)*(T0+1)*J)
  dim(Y_cf) <- c(N0+1,T0+1,J)
  for ( j in 1 : J){
    Yc_j <- Y[1:N0, 1:T0, j, drop = FALSE]
    dim( Yc_j ) <- c(N0,T0)
    Y_cf[,,j] <- rbind( cbind( Yc_j, rowMeans(Y[1:N0, (T0 + 1):Ta, j, drop = FALSE]) ),
                        cbind(t(colMeans(Y[(N0 + 1):N, 1:T0, j, drop = FALSE])), mean(Y[(N0 + 1):N, (T0 + 1):Ta, j, drop = FALSE])) )
  }
  return(Y_cf)
}

#' Return the component-wise sum of decreasing vectors in which NA is taken to mean that the vector has stopped decreasing
#' and we can use the last non-na element. Where both are NA, leave as NA.
pairwise.sum.decreasing = function(x, y) {
  na.x = is.na(x)
  na.y = is.na(y)
  x[is.na(x)] = min(x[!na.x])
  y[is.na(y)] = min(y[!na.y])
  pairwise.sum = x + y
  pairwise.sum[na.x & na.y] = NA
  pairwise.sum
}

#' Convert a long (balanced) panel to a wide multi dimensional matrices (tensors)
#'
#' Converts a data set in panel form to matrix format required by multi_synthdid estimators.
#' A typical long panel date set looks like \[unit, time, outcomes, treatment, covariates\].
#' Multi dimensional version of synthdid requires a balanced panel with simultaneous
#' adoption of treatment: each unit must be observed at all times, and all treated units
#' must begin treatment simultaneosly. This function
#' creates num.units x num.time.periods x num.outcome.vars matrices Y and
#' num.units x num.time.periods for treatment matrix W. It also creates a covariance matrix
#' X with 4 dimensions: num.units x num.time.periods x outcome_var x selected.covariates
#' The 3rd dimension specify for which outcome the covariate is used and the 4th specify the
#' selected covariate.
#' In these matrices, columns are sorted by time, and rows for control units
#' appear before those of treated units.
#'
#' @param df A data.frame with columns consisting of units, time, outcome, and treatment indicator.
#' @param unit The column name corresponding to the unit identifier.
#' @param time The column name corresponding to the time identifier.
#' @param outcome The column names corresponding to the outcome identifiers (can be multiple -- 1xJ).
#' @param treatment The column number corresponding to the treatment status.
#' @param covariates name of the used covariates from the data frame. Default is no covariates.
#' @param eq_s character matrix with LxK elements. J is in the same order as `outcome` is given.
#'  The K elements may be blank or containing the covariates referring to the j-th outcome equation.
#' @return A list with entries `Y`: the data tensor, `N0`: the number of control units, `T0`:
#'  the number of time periods before treatment, `W`: the matrix of treatment indicators,
#'  `X` the data tensor for covariates and `unit`: standing for the ordering of control + treated units.
#'
#' @examples
#' \donttest{
#' # Load tobacco sales in long panel format.
#' data("german_reunification")
#' # Transform to N*T matrix format required for multi_synthdid_estim,
#' # where N is the number of units and T the time periods and
#' # J is the used number of outcome variables.
#' setup <- setup_multi_sdid(german_reunification, unit = 'country',
#'                           time = 'year', outcome = c('gdp','infrate','trade'),
#'                           treatment = 'W')
#'
#' # Compute synthdid estimate
#' multi_synthdid_estim(setup$Y, setup$N0, setup$T0)
#' }
#'
#' @export
setup_multi_sdid = function(df, unit, time, outcome, treatment, covariates = NULL, eq_s = NULL ) {

  if (!is.null(eq_s)){
    if (!all( eq_s %in% covariates | c( eq_s == '' ) ) && all( covariates %in% eq_s )){
      stop("Variables defined in `eq_s` must be the same as the variables defined in covariates or `''`.")
    }
    if( length(outcome) != dim(eq_s)[1] ){
      stop('First dimension of `eq_s` must be the same as number of outcomes defined!')
    }
  }

  keep = c(unit, time, outcome, treatment, covariates )
  if (!all( keep %in% colnames(df))) {
    stop("Column identifiers should be in the dataframe.")
  }

  # Remove not-needed variables
  df = df[keep]
  if (!( is.data.frame(df) | is_tibble(df) )){
    stop("Unsupported input type `df` use data.frame or tibble")
  }
  if (is_tibble(df) ){
    df = as.data.frame(df)
  }
  if (anyNA(df)) {
    stop("Missing values in `df`.")
  }
  if (length(unique(df[, treatment])) == 1) {
    stop("There is no variation in treatment status.")
  }
  if (!all(df[, treatment] %in% c(0, 1))) {
    stop("The treatment status should be in 0 or 1.")
  }

  # Convert potential factor/date columns to character
  df = data.frame(
    lapply(df, function(col) {if (is.factor(col) || inherits(col, "Date")) as.character(col) else col}), stringsAsFactors = FALSE
  )
  val <- as.vector(table(df[, unit], df[, time]))
  if (!all(val == 1)) {
    stop("Input `df` must be a balanced panel: it must have an observation for every unit at every time.")
  }

  # set the order
  panel = df[order(df[, unit], df[, time]), ]
  num.years = length(unique(panel[, time]))
  num.units = length(unique(panel[, unit]))

  W = matrix(panel[,treatment], num.units, num.years, byrow = TRUE,
             dimnames = list(unique(panel[,unit]), unique(panel[,time])))
  w = apply(W, 1, any)                         # indicator for units that are treated at any time
  T0 = unname(which(apply(W, 2, any))[1]-1)    # last period nobody is treated
  N0 = sum(!w)

  if(! (all(W[!w,] == 0) && all(W[,1:T0] == 0) && all(W[w, (T0+1):num.years]==1))) {
    stop("The package cannot use this data. Treatment adoption is not simultaneous.")
  }

  unit.order = order(W[,T0+1], unique(panel[,unit]) )

  # Create multi-dimensional objects
  num.outcome = length( outcome )
  Y = rep(NA,num.units*num.years*num.outcome)
  dim(Y) = c(num.units,num.years,num.outcome)
  for ( j in 1 : num.outcome ){
    Yj = matrix(panel[,outcome[j]], num.units, num.years, byrow = TRUE,
                dimnames = list(unique(panel[,unit]), unique(panel[,time])))
    Y[,,j] = Yj[unit.order, ]
  }

  # Covariates
  if ( is.null( covariates ) ){
    X = NULL
  } else{
    num.covars = length( covariates )
    X = rep(NA,num.units*num.years*num.outcome*num.covars)
    dim(X) = c(num.units,num.years,num.outcome,num.covars)
    if ( is.null(eq_s) ){
      if ( num.outcome > 1 ){
        eq_s = rep(covariates,num.outcome)
        dim(eq_s) = c(num.outcome,num.covars)
        if ( num.covars > 1 ){
          eq_s = t(eq_s)
        }
      } else{
        eq_s = matrix(covariates,nrow=1,ncol=num.covars)
      }
    }
    for ( j in 1 : num.outcome ){
      for ( k in 1 : num.covars)
        if ( !(eq_s[j,k] == '' ) ){
          Xjk = matrix(panel[,eq_s[j,k]], num.units, num.years, byrow = TRUE,
                       dimnames = list(unique(panel[,unit]), unique(panel[,time])))
          X[,,j,k] = Xjk[unit.order,]
        }
    }
  }

  list(Y = Y, N0 = N0, T0 = T0, X = X, W = W[unit.order, ], unit = rownames(W[1:N0,]))
}


