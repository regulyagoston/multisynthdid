#' Plotting functions for multi_sdid objects
#'
#' @param estimate milti_sdid object
#' @param outcome char vector of outcomes to plot
#' @param placebo TRUE/FALSE to plot the placebo outcomes along with the gaps
#' @param se_type `placebo` or `bootrstrap`
#' @param se_replication how many time SE should be calculated
#' @param max_rep_plot how many replication to plot
#' @import dplyr
#' @import ggplot2
#' @import magrittr
#' @export plot_gaps
plot_gaps <- function( estimate, outcome = NULL, placebo = F,
                       se_type = 'placebo', se_replication = 150,
                       max_rep_plot = 150 ){

  stopifnot( class(estimate) == 'multi_synthdid_obj' )

  if ( is.null( outcome ) ){
    outcome <- colnames( estimate$tau_all )
  }


  synt_gap <- multi_synthdid_curves(estimate$tau, complete = T)
  orig_names <- colnames(estimate$tau_all)
  gap_m <- t( synt_gap$tau_curve )
  colnames( gap_m ) <- orig_names

  df_gap <- as.data.frame(gap_m)
  df_gap$x <- estimate$time
  T0 <- attr( estimate$tau, 'setup' )$T0
  tr_time <- df_gap$x[T0+1]

  # Placebo
  if ( placebo ){
    if (is.null(estimate$se_all)){
      se_tau = se.multi_sdid(estimate$tau, method = se_type,
                             replications = se_replication )
    } else{
      se_tau = estimate$se_all
    }
    plac_vals <- se_tau$tau_hat
    J <- dim(plac_vals)[1]
    Tall <- dim(plac_vals)[2]
    n_simul <- dim(plac_vals)[3]
    dim( plac_vals ) <- c(J*Tall*n_simul)
    df_placebo <- as.data.frame(plac_vals)
    colnames(df_placebo) = 'value'
    df_placebo$variable <- rep( orig_names, Tall * n_simul )
    df_placebo$x <- rep( rep( estimate$time, each = J ), n_simul)
    df_placebo$iter <- rep( 1:n_simul, each = J * Tall )
    df_placebo <- df_placebo %>%
      dplyr::filter( iter <= max_rep_plot ) %>%
      dplyr::mutate( placebo = "T" )
  }

  df_gap_long <- df_gap %>%
    tidyr::pivot_longer(!x,names_to = 'variable', values_to = 'value' ) %>%
    dplyr::filter( variable %in% outcome )

  if ( placebo ){
    df_gap_long <- df_gap_long %>%
      mutate( placebo = "F",
              iter = 0 )
    df_gap_long <- rbind(df_gap_long,df_placebo)
    f <- ggplot( df_gap_long, aes(x=x,y=value, group = iter, colour = placebo, alpha = placebo, linewidth = placebo )) +
      geom_line( ) +
      scale_color_manual( values = c('F'='red','T'='grey') )+
      scale_alpha_manual( values = c('F'=1,'T'=0.2) )+
      scale_linewidth_manual( values = c('F'=1.5,'T'=0.5) )+
      geom_hline(yintercept = 0 )+
      geom_vline(xintercept = tr_time, linetype = 'dashed' )+
      facet_wrap( ~variable, scales = 'free_y' ) +
      labs(y='',x='')+
      theme_bw()+
      theme(legend.position = "none")
  } else{
    f <- ggplot( df_gap_long, aes(x=x,y=value)) +
      geom_line( color = 'red', linewidth = 1.5 )+
      geom_hline(yintercept = 0 )+
      geom_vline(xintercept = tr_time, linetype = 'dashed' )+
      facet_wrap( ~variable, scales = 'free_y' ) +
      labs(y='',x='')+
      theme_bw()
  }

  return( f )


}

#' Plotting functions for multi_sdid objects
#'
#' @param estimate milti_sdid object
#' @param outcome char vector of outcomes to plot
#' @import dplyr
#' @import ggplot2
#' @import magrittr
#' @export plot_outcomes
plot_outcomes <- function( estimate, outcome = NULL ){

  stopifnot( class(estimate) == 'multi_synthdid_obj' )

  if ( is.null( outcome ) ){
    outcome <- colnames( estimate$tau_all )
  }

  synt_gap <- multi_synthdid_curves(estimate$tau, complete = T)
  orig_names <- colnames(estimate$tau_all)
  out_m <- t( synt_gap$tr_curve )
  colnames( out_m ) <- orig_names
  syn_m <- t( synt_gap$tau_curve + synt_gap$tr_curve )
  colnames( syn_m ) <- orig_names

  df_out <- as.data.frame(out_m)
  df_out$x <- estimate$time
  df_out$type = 'Observed'

  df_syn <- as.data.frame( syn_m )
  df_syn$x = estimate$time
  df_syn$type = 'Synthetic'

  df_out <- rbind( df_out, df_syn )

  T0 <- attr( estimate$tau, 'setup' )$T0
  tr_time <- df_syn$x[T0+1]

  df_long <- df_out %>%
    dplyr::select( all_of( outcome ), x, type ) %>%
    tidyr::pivot_longer( any_of( outcome ), names_to = 'variable', values_to = 'value' ) %>%
    dplyr::mutate( type = factor(type) )

  f <- ggplot( df_long, aes(x=x, y=value,colour=type,linetype=type)) +
    geom_line( linewidth = 1.5 )+
    geom_vline(xintercept = tr_time, linetype = 'dashed' )+
    facet_wrap( ~variable, scales = 'free_y' ) +
    labs(y='',x='')+
    theme_bw() +
    theme(legend.position = "none")

  return( f )

}
