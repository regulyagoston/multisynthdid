#' Summarize a multisynthdid object
#' @param object The object to summarize
#' @param weights Report weights or not
#' @param weights_print `non-null` - only report non-null values, `all` - report all values.
#' @method summary.multisynthdid
#' @export summary.multisynthdid
summary.multisynthdid = function(object, weights = F, weight.digits=3, weights_print = 'non-null') {

  stopifnot( class(object) == 'multi_synthdid_obj' )

  N0 = attr(object$tau, 'setup')$N0
  T0 = attr(object$tau, 'setup')$T0
  dims <- dim(attr(object$tau, 'setup')$Y)
  N  = dims[1]
  Tall = dims[2]
  J = dims[3]

  unit_weights <- object$weights$unit
  time_weights <- object$weights$time
  if ( weights_print == 'non-null' ){
    unit_weights <- unit_weights[unit_weights$weight>0,]
    time_weights <- time_weights[time_weights$weight>0,]
  } else if ( weights_print != 'all' ){
    stop('Invalid value for weights_print, use `all` or `non-null`!')
  }
  id_sort <- rev(order(unit_weights[,1]))
  unit_weights <- unit_weights[id_sort,]
  unit_weights$weight <- round(unit_weights$weight, digits=weight.digits)
  time_weights$weight <- round(time_weights$weight, digits=weight.digits)

  dimensions <- c( N-N0,
                   N0,
                   round(1 / sum(unit_weights$weight^2),  weight.digits),
                   Tall-T0,
                   T0,
                   round(1 / sum(time_weights$weight^2), weight.digits),
                   J )
  dimensions <- as.matrix(dimensions)
  rownames( dimensions ) <- c('No. Treated units',
                              'No. Control pool',
                              'Effective no. controls',
                              'No. Post-treatment periods',
                              'No. Pre-treatment periods',
                              'Effective no. pre-treatment periods',
                              'No. Outcomes')

  out <- list(estimate = class(object),
              tau = object$tau_all,
              dimensions = dimensions )

  if ( weights ){
    out$controls = unit_weights
    out$periods  = time_weights
  }


  return( out )
}
