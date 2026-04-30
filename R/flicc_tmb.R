#' Build multi-year TMB data list for FLicc
#'
#' @param lfd named list of FLQuant objects, one per gear
#' @param stklen FLStockLen-like object with matching years
#' @param sel_fun selectivity function per gear
#' @param catch_by_gear numeric vector of relative catch by gear
#' @param settings list with CVL, GL, catch.sd, optional linf.sd, Mk.sd, CVL.sd
#' @param years optional vector of years to include
#'
#' @return list for TMB::MakeADFun
data_tmb_flicc <- function(lfd, stklen, sel_fun, catch_by_gear,
                             settings = list(CVL = 0.1, GL = 50, catch.sd = 0.05),
                             years = NULL) {

  stopifnot(is.list(lfd))
  if (is.null(names(lfd))) {
    names(lfd) <- paste0("gear", seq_along(lfd))
  }
  gear_names <- names(lfd)
  ngear <- length(lfd)

  lens <- as.numeric(dimnames(lfd[[1]])$len)
  yrs_all <- as.numeric(dimnames(lfd[[1]])$year)

  if (is.null(years)) {
    years <- yrs_all
  }
  years <- as.numeric(years)
  nyear <- length(years)
  nlen  <- length(lens)

  # observed length frequencies: obs[len, year, gear]
  obs <- array(NA_real_,
               dim = c(nlen, nyear, ngear),
               dimnames = list(len = lens, year = years, gear = gear_names))

  for (g in seq_len(ngear)) {
    x <- lfd[[g]]
    xyrs <- as.numeric(dimnames(x)$year)
    if (!all(years %in% xyrs)) {
      stop("Not all requested years found in lfd gear: ", gear_names[g])
    }

    for (iy in seq_along(years)) {
      y <- years[iy]
      obs[, iy, g] <- as.numeric(as.data.frame(yearMeans(x[, ac(y)]))$data)
    }
  }

  # stock years
  stk_years <- as.numeric(dimnames(stklen)$year)
  if (!all(years %in% stk_years)) {
    stop("Not all requested years found in stklen")
  }

  stkw <- window(stklen, start = min(years), end = max(years))

  # extract year-varying biology as len x year matrices
  wt_df  <- as.data.frame(catch.wt(stkw))
  mat_df <- as.data.frame(mat(stkw))
  m_df   <- as.data.frame(m(stkw))

  # assume ordering len-fast, then year
  wt_mat  <- matrix(wt_df$data,  nrow = nlen, ncol = nyear)
  mat_mat <- matrix(mat_df$data, nrow = nlen, ncol = nyear)
  m_mat   <- matrix(m_df$data,   nrow = nlen, ncol = nyear)

  # scale M around annual mean so Mk remains the average level
  Mscaler <- sweep(m_mat, 2, colMeans(m_mat, na.rm = TRUE), "/")

  # lengths
  LLB <- lens
  step <- c(diff(LLB), tail(diff(LLB), 1))
  Lmid <- LLB + step / 2

  # selectivity
  if (length(sel_fun) != ngear) {
    stop("length(sel_fun) must equal number of gears")
  }

  sel_type <- ifelse(sel_fun == "logistic", 1L,
                     ifelse(sel_fun == "dsnormal", 2L,
                            ifelse(sel_fun == "normal", 3L, NA_integer_)))
  if (any(is.na(sel_type))) {
    stop("Supported selectivity options are logistic, dsnormal, normal")
  }

  sm_n <- ifelse(sel_fun == "logistic", 2L,
                 ifelse(sel_fun == "dsnormal", 3L,
                        ifelse(sel_fun == "normal", 2L, NA_integer_)))
  sm_start <- cumsum(c(0L, head(sm_n, -1L)))  # 0-based for C++





  # catch weights by gear
  catch_wt <- do.call(
    cbind,
    lapply(catch_by_gear, function(x) as.data.frame(x)$data)
  )
  rownames(catch_wt) <- dimnames(catch_by_gear[[1]])$year

  catch_wt <- catch_wt / rowSums(catch_wt)


  if (ncol(catch_wt) != ngear) {
    stop("length(catch_by_gear) must equal number of gears")
  }
  if (nrow(catch_wt) != nyear) {
    stop("catch_by_gear must have one value per year")
  }

  catch_sd <- if (!is.null(settings$catch.sd)) as.numeric(settings$catch.sd) else 0.05

  # observation model
  obs_model <- if (!is.null(settings$obs_model)) settings$obs_model else "nb"

  obs_model_code <- switch(obs_model,
                           nb = 1L,
                           mn = 2L,
                           dm = 3L,
                           stop("Supported obs_model options are 'nb', 'mn', 'dm'")
  )

  pop_model <- if (!is.null(settings$pop_model)) settings$pop_model else "gamma"

  pop_model_code <- switch(pop_model,
                           gamma = 1L,
                           gtg   = 2L,
                           stop("Supported pop_model options are 'gamma' and 'gtg'")
  )



  # quadrature
  gl <- statmod::gauss.quad(settings$GL, kind = "laguerre")
  node <- as.numeric(gl$nodes)
  quad_wt <- as.numeric(gl$weights)

  # life history initial values
  lhpar <- stklen@lhpar
  Linf_init <- as.numeric(lhpar["linf"])
  Mk_init   <- as.numeric(lhpar["Mk"])

  CVL_init <- if (!is.null(settings$CVL)) as.numeric(settings$CVL) else 0.1
  Galpha_init <- 1 / CVL_init^2

  # GTG controls
  ngtg  <- if (!is.null(settings$ngtg))  as.integer(settings$ngtg)  else 13L
  maxsd <- if (!is.null(settings$maxsd)) as.numeric(settings$maxsd) else 2
  Mpow  <- if (!is.null(settings$Mpow))  as.numeric(settings$Mpow)  else 0

  # GTG distribution over Linf #><> check fixed sd
  SDLinf <- CVL_init * Linf_init

  gtgLinfs <- seq(
    from = Linf_init - maxsd * SDLinf,
    to   = Linf_init + maxsd * SDLinf,
    length.out = ngtg
  )

  recP <- dnorm(gtgLinfs, mean = Linf_init, sd = SDLinf)
  recP <- recP / sum(recP)

  # Bin boundaries for GTG mortality scaling
  LBins <- c(LLB, tail(LLB, 1) + tail(step, 1))
  #LBins <- c(Lmid, tail(Lmid, 1) + tail(step, 1))
  #LBins <- LLB#c(LLB, tail(LLB, 1) + tail(step, 1))


  # LBSPR-style optional length-based M/K scaling
  MKL <- Mk_init * (Linf_init / LBins)^Mpow

  # Matrix: one column per GTG
  MKMat <- matrix(MKL, nrow = length(MKL), ncol = ngtg)


  prior_mu <- numeric()
  prior_sd <- numeric()
  prior_code <- integer()

  if (!is.null(settings$linf.sd)) {
    prior_mu   <- c(prior_mu, log(Linf_init))
    prior_sd   <- c(prior_sd, as.numeric(settings$linf.sd))
    prior_code <- c(prior_code, 1L)
  }

  if (!is.null(settings$Mk.sd)) {
    prior_mu   <- c(prior_mu, log(Mk_init))
    prior_sd   <- c(prior_sd, as.numeric(settings$Mk.sd))
    prior_code <- c(prior_code, 3L)
  }

  if (!is.null(settings$CVL.sd)) {
    prior_mu   <- c(prior_mu, -2 * log(CVL_init))
    prior_sd   <- c(prior_sd, 2 * as.numeric(settings$CVL.sd))
    prior_code <- c(prior_code, 2L)
  }

  list(
    nlen       = nlen,
    nyear      = nyear,
    ngear      = ngear,
    Lmid       = as.numeric(Lmid),
    LLB        = as.numeric(LLB),
    obs        = unname(obs),
    Mscaler    = unname(Mscaler),
    wt         = unname(wt_mat),
    mat        = unname(mat_mat),
    catch_wt   = unname(catch_wt),
    catch_sd   = catch_sd,
    pop_model  = as.integer(pop_model_code),
    obs_model  = as.integer(obs_model_code),
    ngtg       = as.integer(ngtg),
    gtgLinfs   = as.numeric(gtgLinfs),
    recP       = as.numeric(recP),
    MKMat      = unname(MKMat),
    sel_type   = as.integer(sel_type),
    sm_start   = as.integer(sm_start),
    sm_n       = as.integer(sm_n),
    node       = node,
    quad_wt    = quad_wt,
    prior_mu   = as.numeric(prior_mu),
    prior_sd   = as.numeric(prior_sd),
    prior_code = as.integer(prior_code),
    Linf_init   = Linf_init,
    Mk_init     = Mk_init,
    CVL_init    = CVL_init,
    Galpha_init = Galpha_init,
    gear_names  = gear_names,
    year_names  = years
  )
}

#' Generate initial parameter values for FLicc TMB fits
#'
#' Builds a named parameter list suitable for \code{TMB::MakeADFun()} from a
#' TMB data list produced by \code{data_tmb_flicc()}. Initial values are derived
#' from the data where possible, using observed modal length by gear as the
#' starting point for selectivity location parameters.
#'
#' @param tmb_data A list produced by \code{data_tmb_flicc()}.
#'
#' @details
#' Initial values are currently set as follows:
#' \itemize{
#'   \item \code{log_Linf}: initialized from \code{tmb_data$Linf_init}
#'   \item \code{log_Mk}: initialized from \code{tmb_data$Mk_init}
#'   \item \code{log_Galpha}: initialized to \code{log(20)}
#'   \item \code{log_Fk}: initialized to \code{log(0.5)} for all gears
#'   \item \code{log_phi}: initialized to \code{log(20)}
#'   \item \code{Sm}: selectivity parameters initialized by gear using the
#'     observed modal length and default spread values
#' }
#'
#' Selectivity parameters are stored in a packed vector \code{Sm}, with the
#' meaning depending on selectivity type.
#'
#' For logistic selectivity, the packed parameter block is
#' \code{c(log_SL50_rel, log_dSL_rel)}, where \code{SL50 = exp(log_SL50_rel) * Linf}
#' and \code{SL95 = SL50 + exp(log_dSL_rel) * Linf}.
#'
#' For normal selectivity, the packed parameter block is
#' \code{c(log_mode_rel, log_sd_rel)}, where \code{mode = exp(log_mode_rel) * Linf}
#' and \code{sd = exp(log_sd_rel) * Linf}.
#'
#' For double-sided normal selectivity, the packed parameter block is
#' \code{c(log_mode_rel, log_lsd_rel, log_rsd_rel)}, where
#' \code{mode = exp(log_mode_rel) * Linf}, \code{lsd = exp(log_lsd_rel) * Linf},
#' and \code{rsd = exp(log_rsd_rel) * Linf}.
#'
#' @return A named list of parameter vectors and scalars for use in
#'   \code{TMB::MakeADFun()}.
#'
#' @examples
#' \dontrun{
#' data(alfonsino)
#' tmb_data <- data_tmb_flicc(lfd,stklen,sel_fun=c("dsnormal","logistic"),catch_by_gear = c(0.7,0.3))
#' par0 <- init_tmb_flicc(tmb_data)
#' str(par0)
#' }
#'
#' @export
#' Initial values for multi-year FLicc TMB
init_tmb_flicc <- function(tmb_data) {

  nlen     <- tmb_data$nlen
  nyear    <- tmb_data$nyear
  ngear    <- tmb_data$ngear
  sel_type <- tmb_data$sel_type
  Lmid     <- tmb_data$Lmid
  obs      <- tmb_data$obs
  Linf     <- tmb_data$Linf_init

  Sm <- numeric(sum(tmb_data$sm_n))
  pos <- 1L

  for (g in seq_len(ngear)) {
    obs_g <- apply(obs[, , g, drop = FALSE], 1, sum, na.rm = TRUE)
    mode_g <- Lmid[which.max(obs_g)]
    mode_rel <- max(mode_g / Linf, 1e-4)

    if (sel_type[g] == 1L) {
      SL50_rel <- mode_rel
      dSL_rel  <- 0.10
      Sm[pos:(pos + 1L)] <- c(log(SL50_rel), log(dSL_rel))
      pos <- pos + 2L

    } else if (sel_type[g] == 2L) {
      lsd_rel <- 0.10
      rsd_rel <- 0.10
      Sm[pos:(pos + 2L)] <- c(log(mode_rel), log(lsd_rel), log(rsd_rel))
      pos <- pos + 3L

    } else if (sel_type[g] == 3L) {
      sd_rel <- 0.10
      Sm[pos:(pos + 1L)] <- c(log(mode_rel), log(sd_rel))
      pos <- pos + 2L
    }
  }

  list(
    log_Linf   = log(tmb_data$Linf_init),
    log_Galpha = log(tmb_data$Galpha_init),
    log_Mk     = log(tmb_data$Mk_init),
    log_Fk     = matrix(log(0.5), nrow = nyear, ncol = ngear,
                        dimnames = list(tmb_data$year_names, tmb_data$gear_names)),
    Sm         = Sm,
    log_phi    = log(10),
    log_sigmaF = tmb_data$prior_sigmaF_mean
  )
}

#' Fit the FLicc model in a single optimizer pass
#'
#' Internal fitting engine for FLicc. This function performs one model fit only
#' and does not apply any automatic retry or stability refit logic.
#'
#' @param lfd A named list of \code{FLQuant} objects containing observed
#'   length-frequency distributions by gear.
#' @param stklen An \code{FLStockLen} object, or compatible FLR object,
#'   containing stock and life-history information.
#' @param sel_fun Character vector giving the selectivity function for each
#'   gear. Supported options include \code{"logistic"}, \code{"normal"}, and
#'   \code{"dsnormal"}.
#' @param catch_by_gear Relative catch by gear. Can be supplied as a numeric
#'   vector or as an \code{FLQuants} object.
#' @param settings A named list of model settings. Recognized elements are
#'   \code{CVL}, \code{GL}, \code{catch.sd}, and \code{prior_sigmaF}. Missing
#'   elements are filled from internal defaults.
#' @param years Optional numeric vector of years to include in the fit. If
#'   \code{NULL}, all years in \code{lfd} are used.
#' @param iter Integer giving the iteration to fit. Iterations are subset in R
#'   before being passed to TMB.
#' @param compile Logical; if \code{TRUE}, compile the TMB model before fitting.
#' @param silent Logical; passed to \code{TMB::MakeADFun()}.
#' @param dll Character string giving the compiled TMB DLL name.
#' @param sel_fixed Optional object used to fix selectivity parameters.
#' @param start Optional numeric vector of starting values for the estimated
#'   TMB parameters. This is mainly used to restart a fit from
#'   \code{fit$opt$par}. If supplied as a vector, its length must match the
#'   free parameter vector after applying any TMB parameter mapping.
#' @param FLRreport Logical; if \code{TRUE}, convert reported outputs using
#'   \code{as_FLQuants()}.
#'
#' @details
#' \code{fiticc_core()} prepares the TMB data using \code{data_tmb_flicc()},
#' initializes model parameters, applies optional selectivity constraints, builds
#' the TMB objective function with \code{TMB::MakeADFun()}, optimizes it using
#' \code{nlminb()}, and returns both the optimizer output and reported
#' quantities.
#'
#' This function is primarily intended for internal use. Most users should call
#' \code{fiticc()}, which wraps \code{fiticc_core()} and can automatically
#' refit marginal models with a higher \code{GL}.
#'
#' @return
#' An object of class \code{"flicc_tmb_fit"} with optimizer results, sdreport
#' output, parameter summaries, reported quantities, and the processed TMB data.
#'
#' @seealso \code{\link{fiticc}}, \code{\link{need_refit_flicc}}
#'
#' @examples
#' data(alfonsino)
#'
#' fit <- fiticc_core(
#'   lfd_alfonsino,
#'   stklen_alfonsino,
#'   sel_fun = c("dsnormal", "logistic"),
#'   catch_by_gear = c(0.7, 0.3)
#' )
#' @keywords internal
fiticc_core <- function(lfd, stklen,
                     sel_fun = c("logistic", "normal", "dsnormal"),
                     catch_by_gear,
                     settings = list(
                       CVL = 0.1,
                       GL = 50,
                       catch.sd = 0.05,
                       prior_sigmaF = c(log(0.5), 0.3, 1)
                     ),
                     years = NULL,
                     iter = 1,
                     compile = FALSE,
                     silent = TRUE,
                     dll = "FLicc",
                     sel_fixed = NULL,
                     start = NULL,
                     n_restart = 3,
                     grad_tol = 1e-3,
                     control = list(),
                     FLRreport=FALSE ) {

  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Package 'TMB' is required")
  }


  default_settings <- list(
    CVL = 0.1,
    GL = 50,
    catch.sd = 0.05,
    pop_model = "gtg",
    obs_model = "mn",
    ngtg = 13,
    maxsd = 2,
    Mpow = 0,
    prior_sigmaF = c(log(0.5), 0.3, 1),
    linf.sd = NULL,
    Mk.sd = NULL,
    CVL.sd = NULL
  )
  # ---- merge user settings with defaults ----
  if (is.null(settings)) settings <- list()

  if(any(settings$prior_sigmaF == FALSE)) settings$prior_sigmaF <- c(NA_real_, NA_real_, 0)

  # fill missing elements
  for (nm in names(default_settings)) {
    if (is.null(settings[[nm]])) {
      settings[[nm]] <- default_settings[[nm]]
    }
  }

  if (!settings$pop_model %in% c("gamma", "gtg")) {
    stop("settings$pop_model must be 'gamma' or 'gtg'")
  }

  if (!settings$obs_model %in% c("nb", "mn", "dm")) {
    stop("settings$obs_model must be 'nb', 'mn', or 'dm'")
  }

  unknown <- setdiff(names(settings), names(default_settings))
  if (length(unknown) > 0) {
    stop("Unknown settings: ", paste(unknown, collapse = ", "))
  }

  # subset iter in R; do not loop over iter in TMB
  lfd_i <- lapply(lfd, function(x) x[, , iter])
  stk_i <- stklen[, , , iter]

  if (is.null(years)) {
    years <- as.numeric(dimnames(lfd_i[[1]])$year)
  }
  years <- as.numeric(years)

  if(!inherits(catch_by_gear, "FLQuants")){

 catch_by_gear  <-  FLCore::FLQuants(Map(function(x,y){
    flq = FLCore::quantSums(x)
    flq[] <- y
    flq
  },x=lfd_i,y=catch_by_gear))

 }

  tmb_data <- data_tmb_flicc(
    lfd = lfd_i,
    stklen = stk_i,
    sel_fun = sel_fun,
    catch_by_gear = catch_by_gear,
    settings = settings,
    years = years
  )

  # Set priors for random walk on F (if nyear > 0)
  tmb_data$prior_sigmaF_mean <- as.numeric(settings$prior_sigmaF[1])
  tmb_data$prior_sigmaF_sd   <- as.numeric(settings$prior_sigmaF[2])
  tmb_data$prior_sigmaF_use  <- as.integer(settings$prior_sigmaF[3])


  parameters <- init_tmb_flicc(tmb_data)

  fix_res <- apply_sel_fixed_flicc(
    parameters = parameters,
    tmb_data   = tmb_data,
    sel_fixed  = sel_fixed
  )
  parameters <- fix_res$parameters

  map <- list()

  if (is.null(settings$linf.sd)) {
    map$log_Linf <- factor(NA)
  }
  if (is.null(settings$Mk.sd)) {
    map$log_Mk <- factor(NA)
  }
  if (is.null(settings$CVL.sd)) {
    map$log_Galpha <- factor(NA)
  }
  # NEW: phi only relevant for NB
  if (settings$obs_model != "nb") {
    map$log_phi <- factor(NA)
  }

  if (length(unique(years)) <= 1) {
    map$log_sigmaF <- factor(NA)
  }

  if (settings$pop_model == "gtg") {
    map$log_Galpha <- factor(NA)
  } else if (is.null(settings$CVL.sd)) {
    map$log_Galpha <- factor(NA)
  }
  if (!is.null(fix_res$map_Sm)) {
    map$Sm <- fix_res$map_Sm
  }
  if (length(map) == 0) {
    map <- NULL
  }

  if (compile) {
    TMB::compile("src/FLicc.cpp")
    dyn.load(dynlib("src/FLicc"))
  }

  obj <- TMB::MakeADFun(
    data = tmb_data[setdiff(names(tmb_data),
                            c("gear_names", "year_names",
                              "Linf_init", "Mk_init", "CVL_init", "Galpha_init"))],
    parameters = parameters,
    map = map,
    DLL = dll,
    silent = silent
  )

  # Optional restart from user-supplied starting values.
  # This is useful for reducing gradients by restarting from a previous optimum,
  # e.g. start = fit$opt$par. The supplied vector must match obj$par after
  # mapping, because mapped/fixed parameters are not optimized by nlminb().
  opt_start <- obj$par
  if (!is.null(start)) {
    if (is.list(start)) start <- unlist(start)
    start <- as.numeric(start)

    if (length(start) != length(opt_start)) {
      stop(
        "Length of 'start' (", length(start), ") does not match the number ",
        "of estimated TMB parameters (", length(opt_start), ")."
      )
    }

    names(start) <- names(opt_start)
    opt_start <- start
  }

  ctrl <- modifyList(
    list(
      iter.max = 1000,
      eval.max = 2000,
      rel.tol =  1e-9,
      x.tol =  1e-7
    ),
    control
  )

  opt <- nlminb(
    start = opt_start,
    objective = obj$fn,
    gradient = obj$gr,
    control = ctrl
  )

  restart_log <- data.frame(
    restart = 0L,
    objective = opt$objective,
    max_gradient = max(abs(obj$gr(opt$par)), na.rm = TRUE),
    accepted = TRUE
  )

  n_restart_used <- 0L

  for (i in seq_len(n_restart)) {

    old_grad <- max(abs(obj$gr(opt$par)), na.rm = TRUE)

    if (old_grad < grad_tol) break

    opt_new <- nlminb(
      start = opt$par,
      objective = obj$fn,
      gradient = obj$gr,
      control = ctrl
    )

    new_grad <- max(abs(obj$gr(opt_new$par)), na.rm = TRUE)

    accept <- is.finite(opt_new$objective) &&
      is.finite(new_grad) &&
      new_grad < old_grad &&
      opt_new$objective <= opt$objective + 1e-6

    restart_log <- rbind(
      restart_log,
      data.frame(
        restart = i,
        objective = opt_new$objective,
        max_gradient = new_grad,
        accepted = accept
      )
    )

    if (accept) {
      opt <- opt_new
    } else {
      break
    }
  }

  opt$max_gradient <- max(abs(obj$gr(opt$par)), na.rm = TRUE)
  opt$n_restarts <- n_restart_used
  opt$restart_log <- restart_log

  obj$env$last.par.best <- opt$par

  rep <- TMB::sdreport(obj)
  adrep <- summary(rep, "fixed")
  report <- obj$report()
  par_tab <- data.frame(
    par = rownames(adrep),
    mpd = adrep[, 1],
    se  = adrep[, 2],
    row.names = NULL,
    stringsAsFactors = FALSE
  )


  fit <- structure(
    list(
      opt = opt,
      obj = obj,
      rep = rep,
      par = par_tab,
      report = report,
      stklen= stklen,
      tmb_data = tmb_data
    ),
    class = "flicc_tmb_fit"
  )

  fit$settings <- settings
  fit$pop_model <- settings$pop_model
  fit$obs_model <- settings$obs_model
  fit$start_used <- !is.null(start)
  fit$max_gradient <- opt$max_gradient
  fit$n_restarts <- opt$n_restarts
  fit$restart_log <- opt$restart_log

  fit$logLik <- LLflicc(fit)


  if(FLRreport){fit$report <- as_FLQuants(fit,stklen)} else {fit$report$FLReport=FALSE}
  return(fit)

}

#' Fit FLicc length-based model using TMB
#'
#' Fits the FLicc length-based equilibrium model to length-frequency data
#' using Template Model Builder (TMB). The model separates the population
#' process from the observation process, allowing flexible combinations of
#' population structure and likelihood assumptions.
#'
#' @param lfd A named list of \code{FLQuant} objects containing observed
#'   length-frequency data by gear. Values can be counts or proportions.
#'
#' @param stklen An \code{FLStockLen}-like object containing life-history
#'   parameters (e.g. Linf, k, M) and matching dimensions to \code{lfd}.
#'
#' @param sel_fun Character vector specifying selectivity function per gear.
#'   Supported options include \code{"logistic"}, \code{"normal"},
#'   and \code{"dsnormal"}.
#'
#' @param catch_by_gear Numeric vector of relative catch proportions by gear.
#'   Must sum to 1.
#'
#' @param settings Named list of model settings. Key options include:
#'   \describe{
#'     \item{\code{pop_model}}{Population model. Options are:
#'       \code{"gtg"} (default), a growth-type-group approximation inspired
#'       by LBSPR, and \code{"gamma"}, the original FLicc/BLICC continuous
#'       growth formulation.}
#'
#'     \item{\code{obs_model}}{Observation likelihood. Options are:
#'       \code{"mn"} (default, multinomial), \code{"nb"} (negative binomial),
#'       and \code{"dm"} (Dirichlet-multinomial).}
#'
#'     \item{\code{CVL}}{Coefficient of variation in length-at-age.}
#'
#'     \item{\code{GL}}{Number of Gauss–Laguerre quadrature nodes (used only
#'       when \code{pop_model = "gamma"}).}
#'
#'     \item{\code{catch.sd}}{Standard deviation for catch penalty.}
#'
#'     \item{\code{ngtg}}{Number of growth-type groups for
#'       \code{pop_model = "gtg"}. Default is 13.}
#'
#'     \item{\code{maxsd}}{Range (in standard deviations) for GTG Linf
#'       distribution. Default is 2.}
#'
#'     \item{\code{Mpow}}{Optional length-dependent mortality scaling
#'       exponent (GTG only). Default is 0.}
#'
#'     \item{\code{prior_sigmaF}}{Prior for annual fishing mortality variation
#'       (mean, sd, use flag).}
#'
#'     \item{\code{linf.sd, Mk.sd, CVL.sd}}{Optional penalties on life-history
#'       parameters.}
#'   }
#' @param start Optional numeric vector of starting values for the estimated
#'   TMB parameters. A common use is \code{start = fit$opt$par} to restart
#'   from a previous fit and reduce residual gradients. Currently this is used
#'   for ordinary multi-year fits; it is ignored for \code{by_year = TRUE}.
#' @param n_restart Integer. Number of automatic optimizer restarts from the
#'   current optimum. Restarts can reduce residual gradients when \code{nlminb()}
#'   stops by relative convergence before reaching a stationary point. Default
#'   is 3. Use 0 to disable automatic restarts.
#' @param grad_tol Numeric. Maximum absolute gradient tolerance used to stop
#'   automatic restarts and by the optional safe-fit check. Default is
#'   \code{1e-3}.
#' @param control Optional list of control arguments passed to \code{nlminb()}.
#'   Values supplied here override the FLicc defaults.
#' @param by_year Logical. If \code{TRUE}, the model is fitted independently
#'   for each year (LBSPR-style). Only the final fit object is retained, while
#'   annual report quantities are updated and combined across years.
#'   Default is \code{FALSE}.
#'
#' @details
#' The FLicc model separates:
#'
#' \itemize{
#'   \item \strong{Population model} (\code{pop_model}):
#'     \itemize{
#'       \item \code{"gtg"}: discrete growth-type-group approximation where
#'         variability in asymptotic length is represented by a finite set
#'         of groups with associated recruitment proportions.
#'       \item \code{"gamma"}: continuous formulation based on a gamma
#'         distribution of growth, solved using Gauss–Laguerre quadrature.
#'     }
#'
#'   \item \strong{Observation model} (\code{obs_model}):
#'     \itemize{
#'       \item \code{"mn"}: multinomial likelihood for length compositions.
#'         Effective sample size is taken directly from the observed LFD
#'         (sum of counts), consistent with LBSPR and fishblicc.
#'       \item \code{"nb"}: negative binomial likelihood for counts with
#'         overdispersion.
#'       \item \code{"dm"}: Dirichlet-multinomial likelihood for
#'         overdispersed compositions.
#'     }
#' }
#'
#' The default model (\code{pop_model = "gtg"}, \code{obs_model = "mn"})
#' provides a fast and stable formulation closely aligned with LBSPR-style
#' length-based assessment models.
#'
#' @return An object of class \code{"flicc_tmb_fit"} containing:
#' \describe{
#'   \item{\code{opt}}{Optimizer output.}
#'   \item{\code{report}}{Reported model quantities from TMB.}
#'   \item{\code{par}}{Estimated parameters.}
#'   \item{\code{tmb_data}}{Input data passed to TMB.}
#'   \item{\code{pop_model}}{Population model used.}
#'   \item{\code{obs_model}}{Observation model used.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- fiticc(
#'   lfd_alfonsino,
#'   stklen_alfonsino,
#'   sel_fun = c("logistic", "logistic"),
#'   catch_by_gear = c(0.7, 0.3)
#' )
#' }
#'
#' @export
fiticc <- function(lfd, stklen,
                   sel_fun = c("logistic", "normal", "dsnormal"),
                   catch_by_gear,
                   settings = list(
                     CVL = 0.1,
                     GL = 30,
                     catch.sd = 0.05,
                     pop_model = "gtg",
                     obs_model = "mn",
                     ngtg = 13,
                     maxsd = 2,
                     Mpow = 0,
                     prior_sigmaF = c(log(0.5), 0.3, 1),
                     linf.sd = NULL,
                     Mk.sd = NULL,
                     CVL.sd = NULL
                   ),
                   years = NULL,
                   iter = 1,
                   compile = FALSE,
                   silent = TRUE,
                   dll = "FLicc",
                   start = NULL,
                   by_year = FALSE,
                   sel_fixed = NULL,
                   FLRreport = TRUE,
                   safe_fit = TRUE,
                   refit_GL = 100,
                   n_restart = 3,
                   grad_tol = 1e-3,
                   control = list()) {


     if(by_year){
      if (!is.null(start)) {
        warning("Argument 'start' is currently ignored when by_year = TRUE.")
      }
      yrs <- dimnames(stklen)$year
      if(!is.null(years)){
       yrs <-  yrs[an(yrs)%in%years]
      }

      nyrs <- length(yrs)
      # Loops through years
      fity <- lapply(yrs,function(y){

           if (inherits(catch_by_gear, "FLQuants")){
             cbg <- window(catch_by_gear,start= yrs[nyrs],end=yrs[nyrs])
           } else {
             cbg <- catch_by_gear
           }

         fiticc_core(
         lfd = window(lfd,start= y,end=y),
         stklen = window(stklen,start= y,y),
         sel_fun = sel_fun,
         catch_by_gear =cbg,
         settings = settings,
         years = NULL,
         iter = iter,
         compile = compile,
         silent = silent,
         dll = dll,
         sel_fixed = sel_fixed,
         n_restart = n_restart,
         grad_tol = grad_tol,
         control = control,
         FLRreport = TRUE
         )}
      )
      repy <-  fity[[1]]$report
      fit_final = fity[[nyrs]]
      dn <-  names(fity[[1]]$report)[-c(3)]
      qs <- dn[c(1:14)]

      for(y in 2:nyrs){
       for(i in seq(qs)){
        flq <- repy[[qs[i]]]
        if(!inherits(repy[[qs[i]]], "FLQuants")){
        repy[[qs[i]]] <- append(repy[[qs[i]]],fity[[y]]$report[[qs[i]]])
        } else {
        for(f in 1:length(repy[[qs[i]]])){
          repy[[qs[i]]][[f]] <- append(repy[[qs[i]]][[f]],fity[[y]]$report[[qs[i]]][[f]])
        }}
       }
      }
      fit_final <- fity[[nyrs]]
      fit_final$stklen <- stklen
      fit_final$report[qs] <- repy[qs]
  }


  if(!by_year){
  fit1 <- fiticc_core(
    lfd = lfd,
    stklen = stklen,
    sel_fun = sel_fun,
    catch_by_gear = catch_by_gear,
    settings = settings,
    years = years,
    iter = iter,
    compile = compile,
    silent = silent,
    dll = dll,
    sel_fixed = sel_fixed,
    start = start,
    n_restart = n_restart,
    grad_tol = grad_tol,
    control = control,
    FLRreport = FLRreport
  )

  if(pop_model_flicc(fit1)%in%"gtg")  safe_fit = FALSE

  refit_used <- FALSE
  fit_final <- fit1

  if (isTRUE(safe_fit) && need_refit_flicc(fit1, grad_tol = grad_tol)) {

    settings2 <- settings
    settings2$GL <- refit_GL

    fit2 <- tryCatch(
      fiticc_core(
        lfd = lfd,
        stklen = stklen,
        sel_fun = sel_fun,
        catch_by_gear = catch_by_gear,
        settings = settings2,
        years = years,
        iter = iter,
        compile = compile,
        silent = silent,
        dll = dll,
        sel_fixed = sel_fixed,
        start = fit1$opt$par,
        n_restart = n_restart,
        grad_tol = grad_tol,
        control = control,
        FLRreport = FLRreport
      ),
      error = function(e) NULL
    )

    if (!is.null(fit2)) {
      fit_final <- fit2
      refit_used <- TRUE
    }
  }

  fit_final$safe_fit <- list(
    enabled = safe_fit,
    refit_used = refit_used,
    initial_GL = settings$GL,
    refit_GL = if (refit_used) refit_GL else NA_real_,
    final_needs_refit = need_refit_flicc(fit_final, grad_tol = grad_tol)
  )
  }

  fit_final
}


#' Check whether a FLicc fit should be refitted
#'
#' Applies a small set of basic convergence checks to a fitted FLicc model and
#' returns whether a conditional refit is recommended.
#'
#' @param fit A fitted object of class \code{"flicc_tmb_fit"}.
#'
#' @details
#' The current checks are intentionally simple and are designed for stable
#' routine use in \code{fiticc()}. A refit is recommended when:
#' \itemize{
#'   \item the fit object is \code{NULL},
#'   \item the optimizer convergence code is not zero,
#'   \item the objective function is not finite, or
#'   \item the Hessian reported by \code{TMB::sdreport()} is not positive
#'   definite.
#' }
#'
#' These checks are used by \code{fiticc()} to decide whether to rerun the model
#' once with a higher \code{GL} value.
#'
#' @return Logical; \code{TRUE} if a refit is recommended, otherwise
#'   \code{FALSE}.
#'
#' @seealso \code{\link{fiticc}}, \code{\link{fiticc_core}}
#'
#' @examples
#' data(alfonsino)
#'
#' fit <- fiticc(
#'   lfd_alfonsino,
#'   stklen_alfonsino,
#'   sel_fun = c("dsnormal", "logistic"),
#'   catch_by_gear = c(0.7, 0.3)
#' )
#'
#' need_refit_flicc(fit)

need_refit_flicc <- function(fit, grad_tol = 0.02) {

  if (is.null(fit))
    return(TRUE)

  if (is.null(fit$opt$convergence) || fit$opt$convergence != 0)
    return(TRUE)

  if (is.null(fit$opt$objective) || !is.finite(fit$opt$objective))
    return(TRUE)

  if (!is.null(fit$rep$pdHess) && !isTRUE(fit$rep$pdHess))
    return(TRUE)



  FALSE
}

#' Apply fixed selectivity settings to FLicc parameters
#'
#' Updates the packed selectivity parameter vector \code{Sm} and builds a TMB
#' mapping object for fixing selected selectivity parameters by gear. This
#' provides a user-facing interface for constraining selectivity parameters on
#' the natural length scale while retaining the internal relative-to-\code{Linf}
#' parameterization used by the TMB model.
#'
#' @param parameters A named parameter list, typically produced by
#'   \code{init_tmb_flicc()}, containing at least the packed selectivity vector
#'   \code{Sm}.
#' @param tmb_data A TMB data list produced by \code{data_tmb_flicc()}.
#' @param sel_fixed An optional named list specifying fixed selectivity values by
#'   gear. Names may be either gear names or character representations of gear
#'   indices. Each element should itself be a list of fixed parameter values on
#'   the natural length scale.
#'
#' @details
#' The function allows users to fix selectivity parameters by gear without
#' needing to work directly with the internal packed \code{Sm} vector.
#'
#' Supported user-facing selectivity inputs are:
#'
#' \strong{Logistic selectivity} (\code{sel_type == 1})
#' \itemize{
#'   \item \code{SL50}: length at 50 percent selectivity
#'   \item \code{SL95}: length at 95 percent selectivity
#'   \item \code{dSL}: difference \code{SL95 - SL50}
#' }
#'
#' \strong{Normal selectivity} (\code{sel_type == 3})
#' \itemize{
#'   \item \code{mode}: modal length of selectivity
#'   \item \code{sd}: standard deviation of the symmetric dome
#' }
#'
#' \strong{Double-sided normal selectivity} (\code{sel_type == 2})
#' \itemize{
#'   \item \code{mode}: modal length of selectivity
#'   \item \code{lsd}: left-side standard deviation
#'   \item \code{rsd}: right-side standard deviation
#' }
#'
#' Internally, all selectivity parameters are stored relative to \code{Linf} on
#' the log scale. The function converts natural-scale user inputs to this
#' internal parameterization and updates the TMB map so that selected values are
#' held fixed during optimization.
#'
#' This makes it possible, for example, to:
#' \itemize{
#'   \item fix a logistic selectivity curve using \code{SL50} and \code{SL95},
#'   \item fix a symmetric dome using \code{mode} and \code{sd},
#'   \item or approximate a one-sided shoulder by fixing a very large
#'     \code{rsd} for a \code{dsnormal} gear.
#' }
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{parameters}}{The updated parameter list with modified
#'       selectivity values in \code{Sm}.}
#'     \item{\code{map_Sm}}{A factor vector suitable for inclusion in the TMB
#'       \code{map} argument. Elements set to \code{NA} are fixed during
#'       optimization. Returns \code{NULL} if no selectivity parameters are
#'       fixed.}
#'   }
#'
#' @examples
#' \dontrun{
#' tmb_data <- data_tmb_flicc(dl)
#' parameters <- init_tmb_flicc(tmb_data)
#'
#' # Fix logistic selectivity using SL50 and SL95
#' sel_fixed <- list(
#'   "Gill net" = list(SL50 = 22, SL95 = 27)
#' )
#'
#' fix_res <- apply_sel_fixed_flicc(parameters, tmb_data, sel_fixed)
#'
#' # Fix right tail of a double-sided normal to be very broad
#' sel_fixed2 <- list(
#'   "Gill net" = list(rsd = 1000)
#' )
#'
#' fix_res2 <- apply_sel_fixed_flicc(parameters, tmb_data, sel_fixed2)
#'
#' # Fix a symmetric normal dome
#' sel_fixed3 <- list(
#'   "Marine set bagnet" = list(mode = 24, sd = 4)
#' )
#'
#' fix_res3 <- apply_sel_fixed_flicc(parameters, tmb_data, sel_fixed3)
#' }
#'
#' @export

apply_sel_fixed_flicc <- function(parameters, tmb_data, sel_fixed = NULL) {
  if (is.null(sel_fixed) || length(sel_fixed) == 0) {
    return(list(parameters = parameters, map_Sm = NULL))
  }

  Sm_map <- factor(seq_along(parameters$Sm))

  gear_names <- tmb_data$gear_names
  sel_type   <- tmb_data$sel_type
  sm_start   <- tmb_data$sm_start
  sm_n       <- tmb_data$sm_n
  Linf       <- tmb_data$Linf_init

  for (nm in names(sel_fixed)) {
    spec <- sel_fixed[[nm]]
    if (is.null(spec)) next

    if (nm %in% gear_names) {
      g <- match(nm, gear_names)
    } else if (grepl("^[0-9]+$", nm)) {
      g <- as.integer(nm)
      if (g < 1 || g > length(gear_names)) {
        stop("Gear index out of range in sel_fixed: ", nm)
      }
    } else {
      stop("Unknown gear in sel_fixed: ", nm)
    }

    start <- sm_start[g] + 1L

    if (sel_type[g] == 1L) {
      # logistic: c(log_SL50_rel, log_dSL_rel)

      current_SL50 <- exp(parameters$Sm[start]) * Linf
      current_dSL  <- exp(parameters$Sm[start + 1L]) * Linf

      if (!is.null(spec$SL50)) {
        current_SL50 <- as.numeric(spec$SL50)
        if (current_SL50 <= 0) stop("SL50 must be > 0")
        parameters$Sm[start] <- log(current_SL50 / Linf)
        Sm_map[start] <- NA
      }

      if (!is.null(spec$SL95)) {
        dSL <- as.numeric(spec$SL95) - current_SL50
        if (dSL <= 0) stop("For logistic selectivity, SL95 must be > SL50")
        parameters$Sm[start + 1L] <- log(dSL / Linf)
        Sm_map[start + 1L] <- NA
      }

      if (!is.null(spec$dSL)) {
        dSL <- as.numeric(spec$dSL)
        if (dSL <= 0) stop("For logistic selectivity, dSL must be > 0")
        parameters$Sm[start + 1L] <- log(dSL / Linf)
        Sm_map[start + 1L] <- NA
      }

    } else if (sel_type[g] == 2L) {
      # dsnormal: c(log_mode_rel, log_lsd_rel, log_rsd_rel)

      if (!is.null(spec$mode)) {
        mode <- as.numeric(spec$mode)
        if (mode <= 0) stop("mode must be > 0")
        parameters$Sm[start] <- log(mode / Linf)
        Sm_map[start] <- NA
      }

      if (!is.null(spec$lsd)) {
        lsd <- as.numeric(spec$lsd)
        if (lsd <= 0) stop("lsd must be > 0")
        parameters$Sm[start + 1L] <- log(lsd / Linf)
        Sm_map[start + 1L] <- NA
      }

      if (!is.null(spec$rsd)) {
        rsd <- as.numeric(spec$rsd)
        if (rsd <= 0) stop("rsd must be > 0")
        parameters$Sm[start + 2L] <- log(rsd / Linf)
        Sm_map[start + 2L] <- NA
      }

    } else if (sel_type[g] == 3L) {
      # normal: c(log_mode_rel, log_sd_rel)

      if (!is.null(spec$mode)) {
        mode <- as.numeric(spec$mode)
        if (mode <= 0) stop("mode must be > 0")
        parameters$Sm[start] <- log(mode / Linf)
        Sm_map[start] <- NA
      }

      if (!is.null(spec$sd)) {
        sd <- as.numeric(spec$sd)
        if (sd <= 0) stop("sd must be > 0")
        parameters$Sm[start + 1L] <- log(sd / Linf)
        Sm_map[start + 1L] <- NA
      }

    } else {
      stop("Unsupported selectivity type for gear ", g)
    }

    if ((start + sm_n[g] - 1L) > length(parameters$Sm)) {
      stop("Selectivity indexing exceeded Sm length for gear ", g)
    }
  }

  list(parameters = parameters, map_Sm = Sm_map)
}


#' Summarise convergence diagnostics for a FLicc fit
#'
#' Provides a compact, lay-readable convergence summary for a fitted FLicc/TMB
#' model. The function checks whether the optimiser terminated normally, whether
#' the maximum absolute gradient is small, whether the Hessian is positive
#' definite, whether standard errors are available, and whether the numerical
#' Hessian has positive eigenvalues.
#'
#' @param fit A fitted FLicc model object, typically returned by `fiticc()`.
#'   The object is expected to contain `fit$opt`, `fit$obj`, `fit$par`, and
#'   either `fit$rep$pdHess` or `fit$sdrep$pdHess`.
#' @param grad_tol Numeric. Tolerance used to judge whether the maximum absolute
#'   gradient is sufficiently small. Default is `1e-3`.
#' @param eig_tol Numeric. Tolerance used to judge whether Hessian eigenvalues
#'   are positive. Default is `1e-8`.
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{check}{Name of the convergence diagnostic.}
#'   \item{value}{Observed diagnostic value.}
#'   \item{passed}{Logical flag indicating whether the check passed.}
#'   \item{interpretation}{Plain-language explanation of the check.}
#' }
#'
#' @details
#' Optimiser convergence alone does not guarantee a reliable fit. A model may
#' stop with optimiser code `0` but still have a large gradient, indicating that
#' the solution is not close to a stationary point. Conversely, standard errors
#' and a positive-definite Hessian indicate that local curvature is estimable,
#' but should still be interpreted together with the gradient check.
#'
#' The Hessian positive-definiteness check first looks for `fit$rep$pdHess`,
#' then for `fit$sdrep$pdHess`. The eigenvalue check recomputes a numerical
#' Hessian using `stats::optimHess()`.
#'
#' @examples
#' \dontrun{
#' fit <- fiticc(lfds, stock, lhpar)
#' flicc_convergence(fit)
#' }
#'
#' @export
flicc_convergence <- function(fit, grad_tol = 1e-2, eig_tol = 1e-8,optimizer.code=FALSE) {

  opt_code <- fit$opt$convergence
  opt_msg  <- fit$opt$message

  max_grad <- if (!is.null(fit$opt$gradient)) {
    max(abs(fit$opt$gradient), na.rm = TRUE)
  } else if (!is.null(fit$obj$gr)) {
    max(abs(fit$obj$gr(fit$opt$par)), na.rm = TRUE)
  } else {
    NA_real_
  }

  pdHess <- if (!is.null(fit$rep$pdHess)) {
    fit$rep$pdHess
  } else if (!is.null(fit$sdrep$pdHess)) {
    fit$sdrep$pdHess
  } else {
    NA
  }

  n_par <- if (!is.null(fit$par)) nrow(fit$par) else length(fit$opt$par)

  n_nan_se <- if (!is.null(fit$par) && "se" %in% names(fit$par)) {
    sum(is.nan(fit$par$se))
  } else {
    NA_integer_
  }

  nan_pars <- if (!is.null(fit$par) && "se" %in% names(fit$par)) {
    unique(fit$par$par[is.nan(fit$par$se)])
  } else {
    character(0)
  }

  eig <- try({
    H <- stats::optimHess(fit$opt$par, fit$obj$fn, fit$obj$gr)
    eigen(H, symmetric = TRUE, only.values = TRUE)$values
  }, silent = TRUE)

  min_eig <- if (!inherits(eig, "try-error")) min(eig, na.rm = TRUE) else NA_real_

  out <- data.frame(
    check = c(
      "Optimizer finished normally",
      "Gradient is acceptable",
      "Hessian is positive definite",
      "Standard errors available",
      "Hessian eigenvalues are positive"
    ),
    value = c(
      paste0(opt_code, " - ", opt_msg),
      signif(max_grad, 4),
      as.character(pdHess),
      paste0(n_par - n_nan_se, " / ", n_par, " parameters"),
      signif(min_eig, 4)
    ),
    passed = c(
      identical(opt_code, 0L),
      isTRUE(max_grad < grad_tol),
      isTRUE(pdHess),
      isTRUE(n_nan_se == 0),
      isTRUE(min_eig > eig_tol)
    ),
    interpretation = c(
      "The optimizer found a numerical optimum.",
      "The fit is close to a stationary point.",
      "The model curvature is well behaved.",
      "Uncertainty could be estimated for all parameters.",
      "No flat or negative curvature directions were detected."
    ),
    stringsAsFactors = FALSE
  )

  if(!optimizer.code) out <- out[-1,]
  out
}

