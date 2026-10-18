#' Fit diagnostics stored with an FLicc fit (internal)
#'
#' Computed at the end of [fiticc_core()] and stored in `fit$diagnostics`:
#' Newton steps, Hessian check, Newton decrement, maximum gradient, a
#' convergence tier and, for an estimated Dirichlet-multinomial precision,
#' whether it sits at or above its prior centre (near-multinomial).
#' @param fit A fit object from [fiticc_core()].
#' @return A list.
#' @noRd
diagnostics_flicc <- function(fit) {
  nw <- fit$opt$newton
  pd <- if (!is.null(fit$rep$pdHess)) isTRUE(fit$rep$pdHess) else
    if (!is.null(nw$hess_pd)) nw$hess_pd else NA
  dec <- if (!is.null(nw$decrement)) nw$decrement else NA_real_
  mg <- if (!is.null(fit$opt$max_gradient)) fit$opt$max_gradient else NA_real_
  tier <- conv_tier_flicc(fit$opt$objective, pd, dec, mg)
  list(newton_steps = if (!is.null(nw$steps)) nw$steps else 0L,
       hess_pd = pd, decrement = dec, max_gradient = mg,
       optimizer_code = fit$opt$convergence, conv_tier = tier,
       dm = dm_status_flicc(fit))
}

#' Convergence tier from objective, Hessian, Newton decrement and gradient
#' @noRd
conv_tier_flicc <- function(objective, pd, dec, max_grad,
                            dec_ok = 1e-3, dec_fail = 0.1, grad_tol = 0.01) {
  if (!is.finite(objective) || isFALSE(pd) || (is.finite(dec) && dec > dec_fail)) return("fail")
  if ((is.finite(dec) && dec > dec_ok) || (is.finite(max_grad) && max_grad > grad_tol)) return("warning")
  "ok"
}

#' Status of an estimated Dirichlet-multinomial precision (internal)
#'
#' The precision is "at bound" when the estimate is at or above the centre of
#' its prior (beta >= 10 max n for "saturating", theta >= 10 for "linear"):
#' the compositions show no overdispersion beyond the input sample size, the
#' likelihood is close to multinomial and the value is set by the prior.
#' @noRd
dm_status_flicc <- function(fit) {
  td <- fit$tmb_data
  if (!identical(fit$obs_model, "dm") || is.null(td$dm_type) || td$dm_type == 0L) return(NULL)
  pl <- tryCatch(fit$obj$env$parList(par = fit$obj$env$last.par.best), error = function(e) NULL)
  if (is.null(pl) || is.null(pl$log_dm)) return(NULL)
  est <- as.numeric(pl$log_dm)
  mu <- as.numeric(td$dm_prior_mu)
  if (length(mu) != length(est)) mu <- rep(mu[1], length(est))
  data.frame(gear = td$gear_names,
             type = c("fixed", "linear", "saturating")[td$dm_type + 1L],
             estimate = exp(est), prior_centre = exp(mu),
             at_bound = est >= mu - 1e-6, row.names = NULL)
}

#' Print fit messages (internal; switched off with verbose = FALSE)
#' @noRd
messages_flicc <- function(fit) {
  d <- fit$diagnostics
  if (is.null(d)) return(invisible(NULL))
  if (identical(d$conv_tier, "fail"))
    message("FLicc: convergence check failed (Hessian not positive definite or large Newton decrement). ",
            "See flicc_convergence(fit).")
  else if (identical(d$conv_tier, "warning"))
    message("FLicc: convergence warning (max gradient ", signif(d$max_gradient, 2),
            ", Newton decrement ", signif(d$decrement, 2), "). See flicc_convergence(fit).")
  if (!is.null(d$dm) && any(d$dm$at_bound))
    message("FLicc: Dirichlet-multinomial precision at its upper bound for ",
            paste(d$dm$gear[d$dm$at_bound], collapse = ", "),
            ": the compositions show no overdispersion beyond the input sample size, ",
            "so the likelihood is close to multinomial (n_eff about n) and the precision is set by its prior. ",
            "Use raw counts or a larger ESS to estimate the overdispersion.")
  invisible(NULL)
}
