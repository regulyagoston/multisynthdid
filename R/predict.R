#' Predict a multisynthdid object
#' @param object The object to predict
#' @param newdata the new data to use for prediction
#' @param formula_yx formula containing the variables/structure to predict
#' @method predict.multisynthdid
#' @export predict.multisynthdid
#'
#'
#'

predict.multisynthdid <- function( object, formula_yx, treatment, unit, time, newdata,
                     SE = T, se_type = 'placebo', se_replication = 200 ){


opts_0 <- attr(object$tau,"opts")
spec_0 <- attr(object$tau,'specs')
weights_0 <- attr(object$tau,"weights")

pred <- multi_sdid( formula_yx, treatment, unit, time, newdata,
            SE = T, se_type = 'placebo', se_replication = 200,
            noise.level = opts_0$noise.level,
            eta.omega = spec_0$eta.omega,
            eta.lambda = spec_0$eta.lambda,
            zeta.omega  = opts_0$zeta.omega,
            zeta.lambda = opts_0$zeta.lambda,
            omega.intercept = opts_0$omega.intercept,
            lambda.intercept = opts_0$lambda.intercept,
            weights = list(omega = weights_0$omega, # unit weights
                           lambda = weights_0$lambda, # time weights
                           theta = NULL ), #outcome weights
            update.omega = F,
            update.lambda = F,
            min.decrease = NULL,
            max.iter = opts_0$max.iter,
            sparsify = spec_0$sparsify_function,
            max.iter.pre.sparsify = spec_0$max.iter.pre.sparsify,
            standardize = spec_0$standardize,
            scale_vcov = spec_0$scale_vcov )

return(pred)

}
