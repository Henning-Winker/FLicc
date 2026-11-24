#include <TMB.hpp>

// ---------------------------------------------------------
// Selectivity (unchanged)
// ---------------------------------------------------------
template<class Type>
vector<Type> calc_sel(const vector<Type>& L,
                      int sel_type,
                      const vector<Type>& pars,
                      Type Linf)
{
  int n = L.size();
  vector<Type> sel(n);

  if(sel_type == 1) {
    Type SL50 = exp(pars(0)) * Linf;
    Type dSL  = exp(pars(1)) * Linf;
    Type Delta = dSL;

    for(int i = 0; i < n; i++) {
      sel(i) = Type(1) / (Type(1) + exp(-log(Type(19.0)) * (L(i) - SL50) / Delta));
    }

  } else if(sel_type == 2) {

    Type mode = exp(pars(0)) * Linf;
    Type lsd  = exp(pars(1)) * Linf;
    Type rsd  = exp(pars(2)) * Linf;
    Type two  = Type(2);

    // AD-safe branch: the mode is estimated, so a plain if() would be frozen
    // at the side of the mode each bin had when the tape was recorded.
    for(int i = 0; i < n; i++) {
      Type d = L(i) - mode;
      Type sd_i = CppAD::CondExpLt(L(i), mode, lsd, rsd);
      sel(i) = exp(-(d * d) / (two * sd_i * sd_i));
    }

  } else if(sel_type == 3) {

    Type mode = exp(pars(0)) * Linf;
    Type sd   = exp(pars(1)) * Linf;
    Type two  = Type(2);

    for(int i = 0; i < n; i++) {
      Type d = L(i) - mode;
      sel(i) = exp(-(d * d) / (two * sd * sd));
    }

  } else if(sel_type == 4) {

    // binormal: sum of two normal curves, scaled to a maximum of 1
    //   pars: log(mode1/Linf), log(sd1/Linf), log((mode2 - mode1)/Linf),
    //         log(sd2/Linf), log(h) with h the height of the second curve
    //         relative to the first before scaling
    Type m1  = exp(pars(0)) * Linf;
    Type s1  = exp(pars(1)) * Linf;
    Type m2  = m1 + exp(pars(2)) * Linf;
    Type s2  = exp(pars(3)) * Linf;
    Type h   = exp(pars(4));
    Type two = Type(2);
    Type smax = Type(0);
    for(int i = 0; i < n; i++) {
      Type d1 = L(i) - m1, d2 = L(i) - m2;
      sel(i) = exp(-(d1 * d1) / (two * s1 * s1)) + h * exp(-(d2 * d2) / (two * s2 * s2));
      // AD-safe running maximum
      smax = CppAD::CondExpGt(sel(i), smax, sel(i), smax);
    }
    for(int i = 0; i < n; i++) sel(i) = sel(i) / smax;

  } else if(sel_type == 5) {

    // bidsnormal: normal curve plus a double-sided normal (right limb rsd2,
    // e.g. fixed at 1000 for an asymptotic second mode), scaled to max 1
    //   pars: log(mode1/Linf), log(sd1/Linf), log((mode2 - mode1)/Linf),
    //         log(lsd2/Linf), log(rsd2/Linf), log(h)
    Type m1   = exp(pars(0)) * Linf;
    Type s1   = exp(pars(1)) * Linf;
    Type m2   = m1 + exp(pars(2)) * Linf;
    Type lsd2 = exp(pars(3)) * Linf;
    Type rsd2 = exp(pars(4)) * Linf;
    Type h    = exp(pars(5));
    Type two  = Type(2);
    Type smax = Type(0);
    for(int i = 0; i < n; i++) {
      Type d1 = L(i) - m1, d2 = L(i) - m2;
      Type sd2_i = CppAD::CondExpLt(L(i), m2, lsd2, rsd2);
      sel(i) = exp(-(d1 * d1) / (two * s1 * s1)) + h * exp(-(d2 * d2) / (two * sd2_i * sd2_i));
      smax = CppAD::CondExpGt(sel(i), smax, sel(i), smax);
    }
    for(int i = 0; i < n; i++) sel(i) = sel(i) / smax;

  } else {
    error("Unknown sel_type");
  }

  return sel;
}

// ---------------------------------------------------------
// Negative binomial
// ---------------------------------------------------------
template<class Type>
Type dnbinom_phi(Type x, Type mu, Type phi, int give_log = 0)
{
  Type logres =
    lgamma(x + phi) - lgamma(phi) - lgamma(x + Type(1)) +
    phi * log(phi / (phi + mu)) +
    x   * log(mu  / (phi + mu));

  if(give_log) return logres;
  return exp(logres);
}

// ---------------------------------------------------------
// Population recursion (gamma-blicc)
// ---------------------------------------------------------
template<class Type>
vector<Type> pop_len_glq(const vector<Type>& node,
                         const vector<Type>& quad_wt,
                         const vector<Type>& Len,
                         const vector<Type>& Zki,
                         Type Galpha,
                         Type Gbeta)
{
  int nv  = node.size();
  int LN  = Len.size();
  int LN1 = LN - 1;

  Type lgamma_Galpha = lgamma(Galpha);
  Type Galpha_1      = Galpha - Type(1);

  vector<Type> surv(LN);
  surv.setZero();

  vector<Type> x_beta(nv);
  vector<Type> log_x_beta(nv);

  vector<Type> Zii(LN);
  Zii.setZero();

  vector<Type> NI(LN);
  NI.setZero();

  Zii(0) = -Zki(0);
  for(int i = 1; i < (LN - 1); ++i) {
    Zii(i) = Zki(i - 1) - Zki(i);
  }

  for(int i = 0; i < nv; ++i) {
    x_beta(i)     = node(i) / Gbeta;
    log_x_beta(i) = log(x_beta(i));

    surv(0) += exp(
      log(node(i) + Gbeta * Len(0)) * Galpha_1
    - Gbeta * Len(0)
      - lgamma_Galpha
    ) * quad_wt(i);
  }

  for(int Li = 1; Li < LN; ++Li) {
    Type Ln = Len(Li);

    for(int i = 0; i < nv; ++i) {
      Type lim = Type(0);

      for(int Lii = 0; Lii < Li; ++Lii) {
        Type Lrange = Ln - Len(Lii);
        lim += log(x_beta(i) + Lrange) * Zii(Lii);
      }

      surv(Li) += exp(
        lim
        + log_x_beta(i) * Zki(Li - 1)
        + log(node(i) + Gbeta * Ln) * Galpha_1
        - (Gbeta * Ln + lgamma_Galpha)
      ) * quad_wt(i);
    }
  }

  for(int i = 0; i < LN1; ++i) {
    NI(i) = (surv(i) - surv(i + 1)) / Zki(i);
  }

  NI(LN1) = surv(LN1) / Zki(LN1);

  return NI;
}
// ---------------------------------------------------------
// Population recursion (gtg-lbspr)
// ---------------------------------------------------------
// ---------------------------------------------------------
// Negative log-likelihood of one composition o (counts, length x sex cells)
// given predicted proportions p, for the sexed part of a gear-year. Same
// observation models as the unsexed compositions in the main template.
// ---------------------------------------------------------
template<class Type>
Type comp_nll_flicc(const vector<Type>& o, const vector<Type>& p, int obs_model,
                    int dm_type, Type log_dm_g, Type phi)
{
  Type nll = Type(0);
  int nb = o.size();
  Type o_tot = o.sum();
  if(o_tot <= Type(0)) return nll;
  if(obs_model == 1) {
    for(int k = 0; k < nb; k++) nll -= dnbinom_phi(o(k), o_tot * p(k) + Type(1e-12), phi, 1);
  } else if(obs_model == 2) {
    for(int k = 0; k < nb; k++) {
      Type pp = CppAD::CondExpLt(p(k), Type(1e-12), Type(1e-12), p(k));
      nll -= o(k) * log(pp);
    }
  } else {
    Type alpha0 = o_tot;
    if(dm_type == 1) alpha0 = exp(log_dm_g) * o_tot;
    if(dm_type == 2) alpha0 = exp(log_dm_g);
    nll -= lgamma(alpha0);
    nll += lgamma(o_tot + alpha0);
    for(int k = 0; k < nb; k++) {
      Type pp = CppAD::CondExpLt(p(k), Type(1e-12), Type(1e-12), p(k));
      Type a = alpha0 * pp;
      nll -= lgamma(o(k) + a);
      nll += lgamma(a);
    }
  }
  return nll;
}

template<class Type>
vector<Type> pop_len_gtg(const vector<Type>& LLB,
                         const vector<Type>& gtgLinfs,
                         const matrix<Type>& ZKLMat,
                         const vector<Type>& recP)
{
  int nlen = LLB.size();
  int ngtg = gtgLinfs.size();

  vector<Type> N(nlen);
  N.setZero();

  // Build upper bin boundary vector of length nlen + 1
  vector<Type> LBins(nlen + 1);
  for(int l = 0; l < nlen; l++) {
    LBins(l) = LLB(l);
  }

  if(nlen > 1) {
    Type step_last = LLB(nlen - 1) - LLB(nlen - 2);
    LBins(nlen) = LLB(nlen - 1) + step_last;
  } else {
    LBins(nlen) = LLB(0) + Type(1.0);
  }

  matrix<Type> NPRFished(nlen + 1, ngtg);
  NPRFished.setZero();

  // recruitment entering first boundary
  for(int g = 0; g < ngtg; g++) {
    NPRFished(0, g) = recP(g);
  }

  // recursion across boundaries
  for(int g = 0; g < ngtg; g++) {
    for(int l = 1; l < (nlen + 1); l++) {

      Type num = gtgLinfs(g) - LBins(l);
      Type den = gtgLinfs(g) - LBins(l - 1);

      // AD-safe floor: gtgLinfs depend on estimated Linf and CVL, so a plain
      // if() would be frozen at the values used when the tape was recorded.
      num = CppAD::CondExpLe(num, Type(1e-12), Type(1e-12), num);
      den = CppAD::CondExpLe(den, Type(1e-12), Type(1e-12), den);

      NPRFished(l, g) =
        NPRFished(l - 1, g) * pow(num / den, ZKLMat(l - 1, g));
    }
  }

  // convert boundary survivorship to numbers in bins
  for(int l = 0; l < nlen; l++) {
    for(int g = 0; g < ngtg; g++) {
      N(l) += (NPRFished(l, g) - NPRFished(l + 1, g)) / ZKLMat(l, g);
    }
  }

  return N;
}

// ---------------------------------------------------------
// GTG survivorship at the lower bound of each length class, i.e. at the
// length label, as LBB without the pile-up correction (Pile = 0): sum over
// groups of the numbers reaching LBins(l). Same recursion as pop_len_gtg().
// ---------------------------------------------------------
template<class Type>
vector<Type> pop_bound_gtg(const vector<Type>& LLB,
                           const vector<Type>& gtgLinfs,
                           const matrix<Type>& ZKLMat,
                           const vector<Type>& recP)
{
  int nlen = LLB.size();
  int ngtg = gtgLinfs.size();
  vector<Type> LBins(nlen + 1);
  for(int l = 0; l < nlen; l++) LBins(l) = LLB(l);
  LBins(nlen) = (nlen > 1) ? Type(LLB(nlen - 1) + (LLB(nlen - 1) - LLB(nlen - 2))) : Type(LLB(0) + Type(1.0));
  vector<Type> Nb(nlen);
  Nb.setZero();
  for(int g = 0; g < ngtg; g++) {
    Type surv = recP(g);
    Nb(0) += surv;
    for(int l = 1; l < nlen; l++) {
      Type num = gtgLinfs(g) - LBins(l);
      Type den = gtgLinfs(g) - LBins(l - 1);
      num = CppAD::CondExpLe(num, Type(1e-12), Type(1e-12), num);
      den = CppAD::CondExpLe(den, Type(1e-12), Type(1e-12), den);
      surv = surv * pow(num / den, ZKLMat(l - 1, g));
      Nb(l) += surv;
    }
  }
  return Nb;
}

// ---------------------------------------------------------
// Objective (MULTI-YEAR VERSION)
// ---------------------------------------------------------
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(nlen);
  DATA_INTEGER(nyear);
  DATA_INTEGER(ngear);

  DATA_VECTOR(Lmid);
  DATA_VECTOR(LLB);

  DATA_ARRAY(obs);        // nlen x nyear x ngear
  DATA_MATRIX(Mscaler);   // nlen x nyear
  DATA_MATRIX(wt);        // nlen x nyear
  DATA_MATRIX(mat);       // nlen x nyear

  DATA_MATRIX(catch_wt);
  DATA_SCALAR(catch_sd);
  DATA_INTEGER(obs_model);   // 1 = nb, 2 = multinom, 3 = dirichlet

  DATA_INTEGER(pop_model);   // 1 = gamma, 2 = gtg
  DATA_INTEGER(ngtg);
  DATA_VECTOR(gtg_z);        // standardised GTG Linf nodes, seq(-maxsd, maxsd)
  DATA_VECTOR(recP);         // recruitment weights by GTG (normal on gtg_z)
  DATA_VECTOR(LBins_mk);     // nlen + 1 bin boundaries for length-based M/K
  DATA_SCALAR(Mpow);         // M/K scales as (Linf / L)^Mpow

  DATA_IVECTOR(sel_type);
  DATA_IVECTOR(sm_start);
  DATA_IVECTOR(sm_n);

  DATA_VECTOR(node);
  DATA_VECTOR(quad_wt);

  DATA_VECTOR(prior_mu);
  DATA_VECTOR(prior_sd);
  DATA_IVECTOR(prior_code);
  // Random walk for fk(g) - effectly F
  DATA_SCALAR(prior_sigmaF_mean);
  DATA_SCALAR(prior_sigmaF_sd);
  DATA_INTEGER(prior_sigmaF_use);
  // FM penalty
  DATA_SCALAR(FM_min);
  DATA_SCALAR(FM_max);
  DATA_SCALAR(FMpen_sd);
  // Upper-tail robustness for the composition likelihood
  DATA_SCALAR(rob_eps);      // 0 = off
  DATA_INTEGER(plus_bin);    // 0-based first pooled bin; -1 = off
  DATA_INTEGER(cut_bin);     // 0-based first excluded (tail-cut) bin; -1 = off
  // Per-recruit reference points at the estimates (evaluated only in double
  // mode, i.e. obj$report(), and only when requested; never during fitting)
  DATA_VECTOR(brp_F);        // apical F (or F/M if brp_FMscale) values
  DATA_VECTOR(brp_spr);      // SPR targets (proportions) for Fspr
  DATA_INTEGER(brp_nyears);  // terminal years averaged
  DATA_SCALAR(brp_k);        // von Bertalanffy k (annual scale)
  DATA_INTEGER(brp_FMscale); // 1: brp_F are F/M multipliers of M-at-length
  DATA_SCALAR(brp_spawn);    // spawning time fraction
  DATA_SCALAR(brp_Fmax);     // upper bound for the Fspr search
  DATA_INTEGER(catch_units); // 1: catch_by_gear in weight (default), 0: in numbers
  DATA_INTEGER(dm_type);     // 0: alpha0 = n (fixed), 1: alpha0 = theta n, 2: alpha0 = beta
  DATA_VECTOR(dm_prior_mu);  // weak prior on log_dm by gear (keeps it finite without overdispersion)
  DATA_SCALAR(dm_prior_sd);  // <= 0: no prior
  DATA_SCALAR(pile);         // GTG predicted composition: 1 numbers integrated over the class (default),
                             // 0 survivorship at the class bound (LBB Pile = 0), in between: LBB's mix
  DATA_INTEGER(brp_only);    // 1: per-recruit call, skip years before the
                             //    terminal brp_nyears (objective not used)

  // Sex structure (nsex = 2; nsex = 1 is the one-sex model, unchanged).
  // Sex 1 (female) uses Linf, Mk, wt, mat, Mscaler and is the reference for
  // F/k and SPR; sex 2 (male) has its own Linf, M/k, weight, maturity and M
  // at length. Recruits are split by sex_ratio; F and selectivity are shared.
  DATA_IVECTOR(gear_type);   // ngear: 0 fishery, 1 survey (selectivity only: no F, no catch share)
  DATA_INTEGER(nsex);
  DATA_ARRAY(obs_sex);       // nlen x nyear x ngear x 2: sexed compositions (F, M); zeros = none
  DATA_IMATRIX(sex_type);    // nyear x ngear, as SS3 gender codes for the sexed part:
                             // 0 none, 1 females only, 2 males only, 3 both (joint length x sex)
  DATA_IVECTOR(sex_comb);    // by gear: number of leading (retained) length bins in which the
                             // sexes are combined in sexed samples (0 = none)
  DATA_SCALAR(sex_ratio);    // proportion female at recruitment (1 when nsex = 1)
  DATA_SCALAR(k_ratio);      // k_female / k_male: F/k and time scale of male survival
  DATA_SCALAR(mk_ratio);     // (M/k)_male / (M/k)_female of the inputs (start value of log_Mk_m; not used in the model)
  DATA_MATRIX(Mscaler_m);    // nlen x nyear, male
  DATA_MATRIX(wt_m);         // nlen x nyear, male
  DATA_MATRIX(mat_m);        // nlen x nyear, male

  // FLicc-SPM (spm_use = 1): Pella-Tomlinson state-space surplus production
  // model in spawning biomass (JABBA-Select / JABBA3 Selex state), fitted to
  // catch and CPUE in the same objective as the length compositions. r and m
  // come from the per-recruit equilibrium of the gtg model (with the
  // estimated selectivity) at Fmsy or Fspr; the EB/SB conversion from the
  // same equilibrium on an F grid; SPR by year is linked to SB/SB0 under
  // Beverton-Holt steepness. The SPM part is a translation of the subset of
  // JABBA3::jb3_nll() with constant process error, catch without error and
  // the production lag. spm_use = 0 leaves FLicc unchanged.
  DATA_INTEGER(spm_use);
  DATA_INTEGER(spm_nyear);
  DATA_VECTOR(spm_catch);      // spm_nyear: catch / mean catch (internal units)
  DATA_IVECTOR(spm_lfd_t);     // nyear: SPM year (0-based) of each LFD year
  DATA_INTEGER(spm_nidx);
  DATA_VECTOR(spm_logI);       // log index observations
  DATA_IVECTOR(spm_obs_t);     // SPM year (0-based) of each observation
  DATA_IVECTOR(spm_obs_i);     // index (0-based) of each observation
  DATA_VECTOR(spm_se2);        // fixed observation variance (JABBA: se^2 + min se^2)
  DATA_IVECTOR(spm_idx_gear);  // nidx: gear (0-based) whose selectivity defines the index; -1 SB, -2 joint fishery
  DATA_IVECTOR(spm_idx_units); // nidx: 1 biomass, 0 numbers
  DATA_SCALAR(spm_h);          // Beverton-Holt steepness
  DATA_INTEGER(spm_shape);     // source of r and m: 0 Fmsy, 1 Fspr (spm_spr_tgt)
  DATA_SCALAR(spm_bmsyk);      // > 0: fixed Bmsy/K of the production function (shape m from it); <= 0: per recruit
  DATA_INTEGER(spm_sdb_fl);    // 1: process SD centred on 0.8 sigR (2 Hmsy)^0.45 (FishLife, Schaefer-equivalent r = 2 Fmsy)
  DATA_SCALAR(spm_sigR);       // recruitment SD for spm_sdb_fl
  DATA_SCALAR(spm_spr_tgt);
  DATA_SCALAR(spm_r_sd);       // > 0: log r ~ N(log r(theta), spm_r_sd); <= 0: r = r(theta)
  DATA_SCALAR(spm_m_sd);       // > 0: log m ~ N(log m(theta), spm_m_sd); <= 0: m = m(theta)
  DATA_INTEGER(spm_link_spr);  // 1: log D_spr,y ~ N(log SB_t/SB0, spm_sigma_link)
  DATA_SCALAR(spm_sigma_link);
  DATA_IVECTOR(spm_lagidx);    // spm_nyear: state (0-based) whose depletion drives production in t (lag)
  DATA_VECTOR(spm_bk_pr);      // initial depletion SB_1/SB0: median, log sd, on
  DATA_VECTOR(spm_sdb_pr);     // process error SD: median, log sd, on
  DATA_VECTOR(spm_K_pr);       // K (internal units): median, log sd, on
  DATA_VECTOR(spm_xgrid);      // F grid for the EB/SB curves: multiples of M at Linf, first 0
  DATA_INTEGER(spm_newton_it); // taped Newton iterations for Fmsy / Fspr and m
  DATA_VECTOR(spm_Pbound);     // soft bounds on SB/SB0 (JABBA: 0.02, 1.3)
  DATA_INTEGER(spm_timing);
  DATA_INTEGER(spm_catch_err); // 1: catch observed with error, log C_obs ~ N(log C_t, spm_sdc) (catch deviations)
  DATA_SCALAR(spm_sdc);        // SD of the catch error (log scale)
  DATA_IVECTOR(spm_frw_t);     // years t (0-based, >= 1) with a random-walk increment log H_t - log H_(t-1) (SPiCT-like F)
  DATA_VECTOR(spm_sdf_pr);     // prior on the SD of the F random walk: c(median, log sd, on)    // biomass in the SPM per recruit: 1 start of year (as the discrete SPM), 0 mean of year
  DATA_VECTOR(spm_selw);       // ngear: fixed weights of the gear selectivities for the per-recruit
                               // block; all <= 0: F-weighted joint selectivity of the terminal years
  // conditional age at length (ALK): growth (Linf, k) from age-length keys
  DATA_INTEGER(alk_on);      // 0 = no ALK
  DATA_INTEGER(alk_mkk);     // 1: M fixed at its input, M/k = M / k follows the estimated k
  DATA_SCALAR(alk_logM);     // log M (alk_mkk = 1)
  DATA_SCALAR(alk_t0);       // von Bertalanffy t0: start value (estimated) or fixed value
  DATA_SCALAR(alk_maxsd);    // Linf range in SD (truncated normal)
  DATA_VECTOR(alk_gl_x);     // Gauss-Legendre nodes on [-1, 1]
  DATA_VECTOR(alk_gl_w);     // Gauss-Legendre weights
  DATA_VECTOR(alk_age);      // ages of the keys (rows of alk_obs)
  DATA_IVECTOR(alk_year);    // year index of each key (0-based)
  DATA_VECTOR(alk_tau);      // sampling time within the year of each key
  DATA_ARRAY(alk_obs);       // nage x nlen x nkey: aged fish by age and length class

  PARAMETER(log_Linf);
  PARAMETER(log_Galpha);
  PARAMETER(log_Mk);
  PARAMETER_MATRIX(log_Fk);   // nyear x ngear
  PARAMETER_VECTOR(Sm);
  PARAMETER(log_phi);
  PARAMETER_VECTOR(log_dm);  // DM precision by gear: log theta (linear) or log beta (saturating)
  PARAMETER(log_sigmaF); // random walk F
  PARAMETER(log_Linf_m);   // male Linf   (nsex = 2; fixed unless linf_m.sd)
  PARAMETER(log_Mk_m);     // male M/k    (nsex = 2; fixed unless Mk_m.sd)
  PARAMETER(log_Galpha_m); // male 1/CVL^2 (nsex = 2; fixed unless CVL_m.sd)
  PARAMETER(spm_logK);     // FLicc-SPM (mapped off when spm_use = 0)
  PARAMETER(spm_logr);
  PARAMETER(spm_logm);
  PARAMETER(spm_lbk);
  PARAMETER(spm_logsdb);
  PARAMETER_VECTOR(spm_logq);
  PARAMETER_VECTOR(spm_logsdi);
  PARAMETER_VECTOR(spm_logP);  // log SB_t/SB0 (random effects)
  PARAMETER(spm_trho);         // AR1 correlation of the SPR-link errors (fixed), rho = 2/(1+exp(-x))-1
  PARAMETER_VECTOR(spm_cdev);  // spm_nyear: log catch deviations (catch error; random effects with re = TRUE)
  PARAMETER(spm_logsdf);       // log SD of the F random walk
  PARAMETER(log_k);        // von Bertalanffy k (estimated with an ALK; otherwise fixed)
  PARAMETER(t0g);          // von Bertalanffy t0 (estimated with an ALK unless fixed; enters the ALK only)

  Type nll = 0.0;

  Type Linf   = exp(log_Linf);
  Type Galpha = exp(log_Galpha);
  Type Mk     = (alk_mkk == 1) ? Type(exp(alk_logM - log_k)) : Type(exp(log_Mk));
  Type phi    = exp(log_phi);
  Type Gbeta  = Galpha / Linf;
  Type sigmaF = exp(log_sigmaF);

  // --- GTG growth groups and M/K from current parameters ---
  // Linf groups spread by CVL = 1 / sqrt(Galpha); M/K uses the estimated Mk
  // and Linf, so priors on Linf, Mk and CVL now reach GTG survival.
  Type CVL = Type(1) / sqrt(Galpha);
  vector<Type> gtgLinfs(ngtg);
  for(int g = 0; g < ngtg; g++) {
    gtgLinfs(g) = Linf * (Type(1) + CVL * gtg_z(g));
  }
  matrix<Type> MKMat(nlen + 1, ngtg);
  for(int l = 0; l < (nlen + 1); l++) {
    Type mkl = Mk * pow(Linf / LBins_mk(l), Mpow);
    for(int g = 0; g < ngtg; g++) {
      MKMat(l, g) = mkl;
    }
  }

  // male growth groups and M/K (nsex = 2): own Linf, M/k and CVL, same
  // growth-group grid. Male parameters are separate from the female ones, so
  // either sex can be fixed while the other is estimated.
  Type Linf_m = exp(log_Linf_m);
  Type Mk_m   = exp(log_Mk_m);
  Type CVL_m  = Type(1) / sqrt(exp(log_Galpha_m));
  vector<Type> gtgLinfs_m(ngtg);
  matrix<Type> MKMat_m(nlen + 1, ngtg);
  if(nsex == 2) {
    if(pop_model != 2) error("Sex structure needs pop_model = gtg");
    for(int g = 0; g < ngtg; g++) gtgLinfs_m(g) = Linf_m * (Type(1) + CVL_m * gtg_z(g));
    for(int l = 0; l < (nlen + 1); l++) {
      Type mkl = Mk_m * pow(Linf_m / LBins_mk(l), Mpow);
      for(int g = 0; g < ngtg; g++) MKMat_m(l, g) = mkl;
    }
  } else {
    gtgLinfs_m = gtgLinfs; MKMat_m = MKMat;
  }
  // male unfished numbers per recruit (male time scale), for EPR
  vector<Type> N0_m(nlen);
  N0_m.setZero();
  if(nsex == 2) {
    matrix<Type> ZKL0m(nlen + 1, ngtg);
    for(int g = 0; g < ngtg; g++) {
      for(int l = 0; l < nlen; l++) ZKL0m(l, g) = MKMat_m(l, g) * Mscaler_m(l, 0);
      ZKL0m(nlen, g) = MKMat_m(nlen, g) * Mscaler_m(nlen - 1, 0);
    }
    N0_m = pop_len_gtg(LLB, gtgLinfs_m, ZKL0m, recP);
  }

  // --- COMMON SELECTIVITY ---
  matrix<Type> Sel(nlen, ngear);
  for(int g = 0; g < ngear; g++) {
    vector<Type> pars(sm_n(g));
    for(int j = 0; j < sm_n(g); j++) {
      pars(j) = Sm(sm_start(g) + j);
    }
    Sel.col(g) = calc_sel(Lmid, sel_type(g), pars, Linf);
  }

  // --- COMMON M, Z0, N0 ---
  vector<Type> M(nlen), Z0(nlen), N0(nlen);

  for(int l = 0; l < nlen; l++) {
    M(l)  = Mk * Mscaler(l, 0);   // common across years by design
    Z0(l) = M(l);
  }

  if(pop_model == 1) {
    N0 = pop_len_glq(node, quad_wt, LLB, Z0, Galpha, Gbeta);

  } else if(pop_model == 2) {
    matrix<Type> ZKL0(nlen + 1, ngtg);
    ZKL0.setZero();

    for(int g = 0; g < ngtg; g++) {
      for(int l = 0; l < nlen; l++) {
        ZKL0(l, g) = MKMat(l, g) * Mscaler(l, 0);
      }
      ZKL0(nlen, g) = MKMat(nlen, g) * Mscaler(nlen - 1, 0);
    }

    N0 = pop_len_gtg(LLB, gtgLinfs, ZKL0, recP);

  } else {
    error("Unknown pop_model");
  }





  // --- OUTPUT OBJECTS ---
  array<Type> plen(nlen, nyear, ngear);
  matrix<Type> Fk(nyear, ngear);
  vector<Type> spr_y(nyear);
  vector<Type> epr_y(nyear);  // exploitable potential ratio: vulnerable biomass per recruit, fished / unfished
  matrix<Type> N_y(nlen, nyear);
  matrix<Type> plen_all_y(nlen, nyear);
  matrix<Type> catch_share(nyear, ngear);   // predicted catch shares (units of catch_units)
  matrix<Type> catch_share_n(nyear, ngear); // predicted catch shares in numbers
  matrix<Type> dm_alpha0(nyear, ngear);     // Dirichlet precision used
  matrix<Type> dm_neff(nyear, ngear);       // effective sample size n (1 + a0) / (n + a0)
  dm_alpha0.setZero(); dm_neff.setZero();
  matrix<Type> sel_joint(nlen, nyear);
  matrix<Type> Fk_l(nlen, nyear);
  array<Type> plen_sex(nlen, nyear, ngear, 2); // predicted joint length x sex proportions within gear
  matrix<Type> N_y_m(nlen, nyear);             // male numbers per recruit (male k time scale)
  plen_sex.setZero(); N_y_m.setZero();
  vector<Type> fapic_y(nyear);                 // apical F/k by year (joint selectivity), for the SPM F link
  fapic_y.setZero();

  matrix<Type> ZK_y(nlen, nyear);              // Z/K at length by year (for the ALK)
  ZK_y.setZero();
  sel_joint.setZero();
  epr_y.setZero();
  Fk_l.setZero();
  plen.setZero();
  Fk.setZero();
  spr_y.setZero();
  N_y.setZero();
  plen_all_y.setZero();
  plen.setZero();
  Fk.setZero();
  spr_y.setZero();
  N_y.setZero();
  plen_all_y.setZero();


  // --- LOOP OVER YEARS ---
  int y_first = 0;
  if(brp_only == 1) {
    int nyb = brp_nyears < 1 ? 1 : (brp_nyears > nyear ? nyear : brp_nyears);
    y_first = nyear - nyb;
  }
  for(int y = y_first; y < nyear; y++) {

    vector<Type> Z(nlen), N(nlen), Nobs(nlen), NobsF(nlen), NobsM(nlen);
    NobsF.setZero(); NobsM.setZero();
    matrix<Type> F_len(nlen, ngear), Cpred(nlen, ngear), mu(nlen, ngear);

    vector<Type> obs_gear(ngear);
    vector<Type> pred_gear(ngear);

    obs_gear.setZero();
    pred_gear.setZero();

    for(int g = 0; g < ngear; g++) {
      Fk(y,g) = exp(log_Fk(y,g));
    }




    // Z at length
    for(int l = 0; l < nlen; l++) {

      Z(l) = M(l);

      for(int g = 0; g < ngear; g++) {
        F_len(l,g) = (gear_type(g) == 1) ? Type(0) : Type(Fk(y,g) * Sel(l,g));
        Z(l) += F_len(l,g);
      }
    }

    for(int l = 0; l < nlen; l++) ZK_y(l, y) = Z(l);

    // what each gear samples at length: catch rate F(l) for fisheries, the
    // population through the selectivity curve for surveys
    matrix<Type> Wsel(nlen, ngear);
    for(int l = 0; l < nlen; l++)
      for(int g = 0; g < ngear; g++)
        Wsel(l,g) = (gear_type(g) == 1) ? Type(Sel(l,g)) : Type(F_len(l,g));

    if(pop_model == 1) {
      N = pop_len_glq(node, quad_wt, LLB, Z, Galpha, Gbeta);
      Nobs = N;

    } else if(pop_model == 2) {
      matrix<Type> ZKL(nlen + 1, ngtg);
      ZKL.setZero();

      for(int g = 0; g < ngtg; g++) {
        for(int l = 0; l < nlen; l++) {
          ZKL(l, g) = MKMat(l, g) * Mscaler(l, y) + (Z(l) - M(l));
        }
        ZKL(nlen, g) = MKMat(nlen, g) * Mscaler(nlen - 1, y) + (Z(nlen - 1) - M(nlen - 1));
      }

      N = pop_len_gtg(LLB, gtgLinfs, ZKL, recP);
      if (pile < Type(1)) {
        vector<Type> Nb = pop_bound_gtg(LLB, gtgLinfs, ZKL, recP);
        for(int l = 0; l < nlen; l++) Nobs(l) = Nb(l) + pile * (N(l) - Nb(l));
      } else {
        Nobs = N;
      }

      if(nsex == 2) {
        // male survival in male time units: M/K of males, F/k scaled by k_f / k_m
        matrix<Type> ZKLm(nlen + 1, ngtg);
        for(int g = 0; g < ngtg; g++) {
          for(int l = 0; l < nlen; l++)
            ZKLm(l, g) = MKMat_m(l, g) * Mscaler_m(l, y) + (Z(l) - M(l)) * k_ratio;
          ZKLm(nlen, g) = MKMat_m(nlen, g) * Mscaler_m(nlen - 1, y) + (Z(nlen - 1) - M(nlen - 1)) * k_ratio;
        }
        vector<Type> Nm = pop_len_gtg(LLB, gtgLinfs_m, ZKLm, recP);
        vector<Type> Nmobs = Nm;
        if (pile < Type(1)) {
          vector<Type> Nbm = pop_bound_gtg(LLB, gtgLinfs_m, ZKLm, recP);
          for(int l = 0; l < nlen; l++) Nmobs(l) = Nbm(l) + pile * (Nm(l) - Nbm(l));
        }
        for(int l = 0; l < nlen; l++) {
          N_y_m(l, y) = Nm(l);
          // numbers in real time, relative to the female time scale (x k_f / k_m),
          // and recruitment split by sex
          NobsF(l) = sex_ratio * Nobs(l);
          NobsM(l) = (Type(1) - sex_ratio) * k_ratio * Nmobs(l);
        }
      }

    } else {
      error("Unknown pop_model");
    }


    for(int l = 0; l < nlen; l++) {
      N_y(l,y) = N(l);
    }

    Type spr0 = 0.0;
    Type sprf = 0.0;

    for(int l = 0; l < nlen; l++) {
      spr0 += N0(l) * mat(l,y) * wt(l,y);
      sprf += N(l)  * mat(l,y) * wt(l,y);
    }

    if(spr0 > Type(0)) {
      spr_y(y) = sprf / spr0;
    } else {
      spr_y(y) = Type(0);
    }

    for(int g = 0; g < ngear; g++) {
      for(int l = 0; l < nlen; l++) {
        if(nsex == 2) {
          Cpred(l,g) = (NobsF(l) + NobsM(l)) * Wsel(l,g);
        } else {
          Cpred(l,g) = Nobs(l) * Wsel(l,g);
        }
        pred_gear(g) += Cpred(l,g);
        obs_gear(g)  += obs(l,y,g);
      }
    }

    // predicted gear shares of the catch, compared with catch_by_gear:
    // in weight (catch_units = 1, default; e.g. landings in tonnes) or in
    // numbers (catch_units = 0, as in fishblicc)
    vector<Type> pred_share(ngear);
    for(int g = 0; g < ngear; g++) {
      pred_share(g) = Type(0);
      if(gear_type(g) == 1) continue;          // surveys take no catch
      for(int l = 0; l < nlen; l++)
        if(nsex == 2 && catch_units == 1)
          pred_share(g) += (NobsF(l) * wt(l,y) + NobsM(l) * wt_m(l,y)) * F_len(l,g);
        else
          pred_share(g) += (catch_units == 1) ? Type(Cpred(l,g) * wt(l,y)) : Type(Cpred(l,g));
    }
    vector<Type> pred_pgear(ngear);
    Type pred_total = Type(0);

    pred_pgear.setZero();

    for(int g = 0; g < ngear; g++) {
      pred_total += pred_share(g);
    }

    for(int g = 0; g < ngear; g++) {
      if(pred_total > Type(0)) {
        pred_pgear(g) = pred_share(g) / pred_total;
      } else {
        pred_pgear(g) = Type(0);
      }
    }
    for(int g = 0; g < ngear; g++) catch_share(y, g) = pred_pgear(g);
    Type pred_total_n = Type(0);
    for(int g = 0; g < ngear; g++) if(gear_type(g) == 0) pred_total_n += pred_gear(g);
    for(int g = 0; g < ngear; g++)
      catch_share_n(y, g) = (pred_total_n > Type(0) && gear_type(g) == 0) ? Type(pred_gear(g) / pred_total_n) : Type(0);


    Type fmax = Type(0);
    Type M_at_fmax = M(0);

    for(int l = 0; l < nlen; l++) {

      Type fsum = Type(0);

      for(int g = 0; g < ngear; g++) {
        fsum += F_len(l,g);
      }

      Fk_l(l,y) = fsum;

      // AD-safe running maximum (which bin is apical depends on parameters)
      M_at_fmax = CppAD::CondExpGt(fsum, fmax, M(l), M_at_fmax);
      fmax      = CppAD::CondExpGt(fsum, fmax, fsum, fmax);
    }
    fapic_y(y) = fmax;

    if (FMpen_sd > Type(0)) {

      Type FM_y = fmax / (M_at_fmax + Type(1e-12)) + Type(1e-12);

      // AD-safe one-sided penalties (zero inside [FM_min, FM_max])
      nll += CppAD::CondExpLt(FM_y, FM_min, Type(1e6) * pow(log(FM_y) - log(FM_min), 2), Type(0));
      nll += CppAD::CondExpGt(FM_y, FM_max, Type(1e6) * pow(log(FM_y) - log(FM_max), 2), Type(0));
    }


    for(int l = 0; l < nlen; l++) {
      if(fmax > Type(0)) {
        sel_joint(l,y) = Fk_l(l,y) / fmax;
      } else {
        sel_joint(l,y) = Type(0);
      }
    }

    // Exploitable potential ratio (EPR): biomass vulnerable to the fishery
    // (weighted by the combined selectivity of all gears) per recruit at the
    // year's F relative to the unfished state; equivalent to LBB's B/B0
    // (Froese et al. 2018) for one gear.
    {
      Type epr0 = Type(0);
      Type eprf = Type(0);
      for(int l = 0; l < nlen; l++) {
        if(nsex == 2) {
          // both sexes per total recruit, on the female time scale
          epr0 += (sex_ratio * N0(l) * wt(l,y) + (Type(1) - sex_ratio) * k_ratio * N0_m(l) * wt_m(l,y)) * sel_joint(l,y);
          eprf += (sex_ratio * N(l)  * wt(l,y) + (Type(1) - sex_ratio) * k_ratio * N_y_m(l,y) * wt_m(l,y)) * sel_joint(l,y);
        } else {
          epr0 += N0(l) * wt(l,y) * sel_joint(l,y);
          eprf += N(l)  * wt(l,y) * sel_joint(l,y);
        }
      }
      epr_y(y) = CppAD::CondExpGt(epr0, Type(0), eprf / (epr0 + Type(1e-30)), Type(0));
    }

    // combined expected length composition of the whole catch (numbers):
    // the gears' catch-at-length summed, i.e. each gear's composition weighted
    // by its predicted catch in numbers
    Type ctot = Type(0);
    for(int g = 0; g < ngear; g++) {
      if(gear_type(g) == 0) ctot += pred_gear(g);
    }

    for(int l = 0; l < nlen; l++) {
      Type ctot_l = Type(0);

      for(int g = 0; g < ngear; g++) {
        if(gear_type(g) == 0) ctot_l += Cpred(l,g);
      }

      if(ctot > Type(0)) {
        plen_all_y(l,y) = ctot_l / ctot;
      } else {
        plen_all_y(l,y) = Type(0);
      }
    }



    // soft constraint on predicted gear shares using relative catch
    // (reference: the last fishery gear; surveys have no catch and are skipped)
    int gref = ngear - 1;
    while(gref > 0 && gear_type(gref) == 1) gref--;
    for(int g = 0; g < gref; g++) {
      if(gear_type(g) == 1) continue;
      if(pred_pgear(g) > Type(0) &&
         pred_pgear(gref) > Type(0) &&
         catch_wt(y,g) > Type(0) &&
         catch_wt(y,gref) > Type(0)) {

        Type logratio_pred  = log(pred_pgear(g) / pred_pgear(gref));
        Type logratio_catch = log(catch_wt(y,g) / catch_wt(y,gref));

        nll -= dnorm(logratio_pred, logratio_catch, catch_sd, true);
      }
    }



    // Likelihood function

    // Within-gear predicted proportions
    for(int g = 0; g < ngear; g++) {
      for(int l = 0; l < nlen; l++) {
        if(pred_gear(g) > Type(0)) {
          plen(l,y,g) = Cpred(l,g) / pred_gear(g);
        } else {
          plen(l,y,g) = Type(0);
        }
      }
    }

    // Observation likelihood
    // Optional treatment of the upper tail (all off by default):
    //  - cut_bin >= 0 excludes bins l >= cut_bin from the likelihood
    //    (tail cut); observed and predicted compositions are both
    //    truncated and the predicted proportions renormalised over the
    //    retained bins, so the fit conditions on fish below the cut;
    //  - plus_bin >= 0 pools observed and predicted counts in bins
    //    l >= plus_bin into bin plus_bin (tail compression), within the
    //    retained bins;
    //  - rob_eps > 0 mixes predicted proportions with a uniform:
    //    p = (1 - rob_eps) * p + rob_eps / nb.
    int nkeep = (cut_bin > 0 && cut_bin < nlen) ? cut_bin : nlen;
    int nb = (plus_bin >= 0 && plus_bin < nkeep) ? (plus_bin + 1) : nkeep;

    for(int g = 0; g < ngear; g++) {
      vector<Type> o(nb), p(nb);
      o.setZero(); p.setZero();
      for(int l = 0; l < nkeep; l++) {
        int k = (l < nb) ? l : (nb - 1);
        o(k) += obs(l,y,g);
        if(pred_gear(g) > Type(0)) p(k) += Cpred(l,g) / pred_gear(g);
      }
      // renormalise over retained bins (no-op without a tail cut)
      Type psum = p.sum();
      if(psum > Type(0)) p = p / psum;
      // observed total over retained bins (= obs_gear(g) without a cut)
      Type o_tot = o.sum();
      if(rob_eps > Type(0)) {
        for(int k = 0; k < nb; k++) {
          p(k) = (Type(1) - rob_eps) * p(k) + rob_eps / Type(nb);
        }
      }

      if(obs_model == 1) {
        // Negative binomial on ESS-scaled counts
        for(int k = 0; k < nb; k++) {
          Type muk = o_tot * p(k) + Type(1e-12);
          if(k < nlen) mu(k,g) = muk;
          nll -= dnbinom_phi(o(k), muk, phi, 1);
        }

      } else if(obs_model == 2) {
        // Multinomial on within-gear compositions
        // obs is already ESS-scaled by lfdess()
        if(o_tot > Type(0) && pred_gear(g) > Type(0)) {
          for(int k = 0; k < nb; k++) {
            // AD-safe floor (a plain if() would be frozen at tape time)
            Type p_pred = CppAD::CondExpLt(p(k), Type(1e-12), Type(1e-12), p(k));
            nll -= o(k) * log(p_pred);
          }
        }

      } else if(obs_model == 3) {
        // Dirichlet-multinomial on within-gear compositions. Precision:
        //  dm_type 0: alpha0 = n (n = observed total over retained bins)
        //  dm_type 1: alpha0 = theta n (linear, Thorson et al. 2017)
        //  dm_type 2: alpha0 = beta (saturating: n_eff -> beta + 1 as n grows)
        if(o_tot > Type(0) && pred_gear(g) > Type(0)) {
          Type alpha0 = o_tot;
          if(dm_type == 1) alpha0 = exp(log_dm(g)) * o_tot;
          if(dm_type == 2) alpha0 = exp(log_dm(g));
          dm_alpha0(y, g) = alpha0;
          dm_neff(y, g) = o_tot * (Type(1) + alpha0) / (o_tot + alpha0);
          nll -= lgamma(alpha0);
          nll += lgamma(o_tot + alpha0);
          for(int k = 0; k < nb; k++) {
            // AD-safe floor (a plain if() would be frozen at tape time)
            Type p_pred = CppAD::CondExpLt(p(k), Type(1e-12), Type(1e-12), p(k));
            Type a = alpha0 * p_pred;
            nll -= lgamma(o(k) + a);
            nll += lgamma(a);
          }
        }

      } else {
        error("Unknown obs_model");
      }
    }

    // Sexed compositions (nsex = 2), as SS3 gender codes 1-3:
    //  3 (both sexes): one joint composition over length x sex, so the sex
    //    ratio at length is information;
    //  1 / 2 (one sex sampled): the composition within that sex.
    // Tail cut and plus group as for the unsexed part, within each sex; same
    // observation model and DM precision by gear. An unsexed (combined)
    // composition of the same gear-year enters separately above.
    if(nsex == 2) {
      for(int g = 0; g < ngear; g++) {
        for(int sx = 0; sx < 2; sx++)
          for(int l = 0; l < nlen; l++) {
            Type cs = (sx == 0) ? NobsF(l) * Wsel(l, g) : NobsM(l) * Wsel(l, g);
            if(pred_gear(g) > Type(0)) plen_sex(l, y, g, sx) = cs / pred_gear(g);
          }
        int st = sex_type(y, g);
        if(st <= 0 || st > 3) continue;     // 0 none; 1 F, 2 M, 3 joint (SS3 gender codes)
        // sexes combined below a length (SS3 'combine males into females at or
        // below a bin'): the first cb retained bins are one cell of both sexes
        // (observed: F + M columns; predicted: females + males), the bins above
        // are female and male cells, all in one composition
        int cb = sex_comb(g);
        if(cb > nb) cb = nb;
        if(cb > 0) {
          int nc = cb + 2 * (nb - cb);
          vector<Type> oc(nc), pc(nc);
          oc.setZero(); pc.setZero();
          for(int sx = 0; sx < 2; sx++) {
            for(int l = 0; l < nkeep; l++) {
              int k = (l < nb) ? l : (nb - 1);
              int idx = (k < cb) ? k : (cb + sx * (nb - cb) + (k - cb));
              oc(idx) += obs_sex(l, y, g, sx);
              pc(idx) += plen_sex(l, y, g, sx);
            }
          }
          Type pcs = pc.sum();
          if(pcs > Type(0)) pc = pc / pcs;
          if(rob_eps > Type(0)) {
            for(int k = 0; k < nc; k++) pc(k) = (Type(1) - rob_eps) * pc(k) + rob_eps / Type(nc);
          }
          if(pred_gear(g) > Type(0))
            nll += comp_nll_flicc(oc, pc, obs_model, dm_type, log_dm(g), phi);
          continue;
        }
        int s0 = (st == 2) ? 1 : 0;          // first sex in the composition
        int s1 = (st == 1) ? 0 : 1;          // last sex in the composition
        int ns = s1 - s0 + 1;
        vector<Type> o(ns * nb), p(ns * nb);
        o.setZero(); p.setZero();
        for(int sx = s0; sx <= s1; sx++) {
          for(int l = 0; l < nkeep; l++) {
            int k = (sx - s0) * nb + ((l < nb) ? l : (nb - 1));
            o(k) += obs_sex(l, y, g, sx);
            p(k) += plen_sex(l, y, g, sx);
          }
        }
        Type psum = p.sum();
        if(psum > Type(0)) p = p / psum;
        if(rob_eps > Type(0)) {
          for(int k = 0; k < ns * nb; k++) p(k) = (Type(1) - rob_eps) * p(k) + rob_eps / Type(ns * nb);
        }
        if(pred_gear(g) > Type(0))
          nll += comp_nll_flicc(o, p, obs_model, dm_type, log_dm(g), phi);
      }
    }

  }


  // --- FLicc-SPM: surplus production model in the same objective ---
  if(spm_use == 1 && brp_only == 0) {
    if(pop_model != 2) error("FLicc-SPM needs pop_model = gtg");
    Type kspm = exp(log_k);     // annual k: estimated with an ALK, else fixed at brp_k (log_k mapped)
    // smooth bounds as JABBA3 (.smax, .smin)
    auto smax_ = [](Type x, Type a, Type eps) { return Type(a + Type(0.5) * ((x - a) + sqrt((x - a) * (x - a) + eps))); };
    auto smin_ = [](Type x, Type a, Type eps) { return Type(a - Type(0.5) * ((a - x) + sqrt((a - x) * (a - x) + eps))); };

    // per-recruit inputs as the brp block: terminal brp_nyears averages of the
    // joint selectivity (apical 1), M, maturity and weight; AD-safe
    int nyb = brp_nyears < 1 ? 1 : (brp_nyears > nyear ? nyear : brp_nyears);
    int yb0 = nyear - nyb;
    vector<Type> sA(nlen), Ms(nlen), matA(nlen), wtA(nlen), Ms_m(nlen), wtA_m(nlen), MA(nlen), MA_m(nlen);
    sA.setZero(); Ms.setZero(); matA.setZero(); wtA.setZero(); Ms_m.setZero(); wtA_m.setZero();
    bool selw_on = (spm_selw.size() == ngear) && (spm_selw.sum() > Type(0));
    for(int l = 0; l < nlen; l++) {
      if(selw_on) for(int g = 0; g < ngear; g++) sA(l) += spm_selw(g) * Sel(l, g);
      for(int y = yb0; y < nyear; y++) {
        if(!selw_on) sA(l) += sel_joint(l, y) / Type(nyb);
        Ms(l)    += Mscaler(l, y)   / Type(nyb);
        matA(l)  += mat(l, y)       / Type(nyb);
        wtA(l)   += wt(l, y)        / Type(nyb);
        Ms_m(l)  += Mscaler_m(l, y) / Type(nyb);
        wtA_m(l) += wt_m(l, y)      / Type(nyb);
      }
    }
    // smooth maximum (p-norm, p = 20) instead of a hard max: r, m, BmsyK and
    // Hmsy do not depend on the scale of sA, and a hard max puts kinks into
    // the fixed-effect surface when the apical bin switches
    Type samax = Type(0);
    for(int l = 0; l < nlen; l++) samax += pow(sA(l) + Type(1e-12), Type(20));
    samax = pow(samax, Type(1) / Type(20));
    for(int l = 0; l < nlen; l++) {
      sA(l) = sA(l) / (samax + Type(1e-30));
      MA(l) = Mk * kspm * Ms(l);                    // annual M at length, female
      MA_m(l) = Mk_m * (kspm / k_ratio) * Ms_m(l);  // male
    }
    Type male_w = (nsex == 2) ? Type(Type(1) - sex_ratio) : Type(0);
    Type Mref = Mk * kspm;                          // annual M at Linf: scale of the F grid

    // numbers at length per recruit at annual apical F x (female; male on the
    // female time scale), as popN() / popNm() of the brp block
    auto popF = [&](Type x, vector<Type>& Zl) {
      matrix<Type> ZK(nlen + 1, ngtg);
      for(int l = 0; l < nlen; l++) Zl(l) = MA(l) + x * sA(l);
      for(int g = 0; g < ngtg; g++) {
        for(int l = 0; l < nlen; l++) ZK(l, g) = MKMat(l, g) * Ms(l) + x * sA(l) / kspm;
        ZK(nlen, g) = MKMat(nlen, g) * Ms(nlen - 1) + x * sA(nlen - 1) / kspm;
      }
      return pop_len_gtg(LLB, gtgLinfs, ZK, recP);
    };
    auto popM = [&](Type x, vector<Type>& Zl) {
      matrix<Type> ZK(nlen + 1, ngtg);
      for(int l = 0; l < nlen; l++) Zl(l) = MA_m(l) + x * sA(l);
      for(int g = 0; g < ngtg; g++) {
        for(int l = 0; l < nlen; l++) ZK(l, g) = MKMat_m(l, g) * Ms_m(l) + x * sA(l) * k_ratio / kspm;
        ZK(nlen, g) = MKMat_m(nlen, g) * Ms_m(nlen - 1) + x * sA(nlen - 1) * k_ratio / kspm;
      }
      vector<Type> NI = pop_len_gtg(LLB, gtgLinfs_m, ZK, recP);
      return vector<Type>(NI * k_ratio);
    };
    vector<Type> Zl0(nlen);
    vector<Type> NI00 = popF(Type(0), Zl0);
    Type scale0 = NI00(0);

    // per recruit at x: SBPR (female), YPR (both sexes), EBPR of the joint
    // fishery (both sexes, weight) and of each index (biomass or numbers)
    // out: 0 SBPR (FLicc convention, for SPR, Fref, BmsyK), 1 YPR, 2 EBPR of the
    // fishery, 3.. EBPR of each index, last SBPR for the EB/SB ratios. With
    // spm_timing = 0 (default) the biomasses are means over the year, so that
    // C / EB = F (Hmsy = Fmsy, the common currency with F). spm_timing = 1
    // (legacy) converts them to the start of the year (x Z / (1 - exp(-Z))).
    int nq_out = 4 + spm_nidx;
    auto prq = [&](Type x) {
      vector<Type> out(nq_out);
      out.setZero();
      vector<Type> Zl(nlen), Zlm(nlen);
      vector<Type> NI = popF(x, Zl);
      vector<Type> NIm(nlen);
      NIm.setZero();
      if(nsex == 2) NIm = popM(x, Zlm);
      for(int l = 0; l < nlen; l++) {
        Type n = NI(l) / scale0;
        Type c = n * x * sA(l);                          // catch = N F (time-integrated numbers)
        out(0) += sex_ratio * n * matA(l) * exp(-Zl(l) * brp_spawn) * wtA(l);
        Type cf = (spm_timing == 1) ? Type(Zl(l) / (Type(1) - exp(-Zl(l)))) : Type(1);
        Type cfm = cf;
        out(nq_out - 1) += sex_ratio * n * cf * matA(l) * wtA(l) *
                           ((spm_timing == 1) ? Type(1) : Type(exp(-Zl(l) * brp_spawn)));
        out(1) += sex_ratio * c * wtA(l);
        out(2) += sex_ratio * n * cf * wtA(l) * sA(l);
        Type nm = Type(0), cm = Type(0);
        if(nsex == 2) {
          nm = NIm(l) / scale0;
          cm = nm * x * sA(l);
          out(1) += male_w * cm * wtA_m(l);
          if(spm_timing == 1) cfm = Zlm(l) / (Type(1) - exp(-Zlm(l)));
          out(2) += male_w * nm * cfm * wtA_m(l) * sA(l);
        }
        for(int i = 0; i < spm_nidx; i++) {
          int gi = spm_idx_gear(i);
          if(gi == -1) continue;                         // SB index: filled below
          Type si = (gi >= 0) ? Type(Sel(l, gi)) : Type(sA(l));
          if(spm_idx_units(i) == 1) out(3 + i) += (sex_ratio * n * cf * wtA(l) + male_w * nm * cfm * wtA_m(l)) * si;
          else                      out(3 + i) += (sex_ratio * n * cf + male_w * nm * cfm) * si;
        }
      }
      for(int i = 0; i < spm_nidx; i++) if(spm_idx_gear(i) == -1) out(3 + i) = out(nq_out - 1);
      return out;
    };
    // Beverton-Holt equilibrium recruitment (R/R0) at SPR
    auto recbh = [&](Type spr) {
      Type R = (Type(4) * spm_h * spr - (Type(1) - spm_h)) / ((Type(5) * spm_h - Type(1)) * spr);
      return smax_(R, Type(1e-6), Type(1e-12));
    };

    // equilibrium on the F grid: depletion P = SB/SB0 and the EB/SB ratios
    int nx = spm_xgrid.size();
    vector<Type> q0 = prq(Type(0));
    Type SBPR0 = q0(0);
    vector<Type> Pg(nx), Yg(nx), lxg(nx), sprg(nx);
    matrix<Type> rhog(nx, 1 + spm_nidx);             // column 0: fishery (removals), 1..: indices
    for(int j = 0; j < nx; j++) {
      Type x = spm_xgrid(j) * Mref;
      vector<Type> q = (j == 0) ? q0 : prq(x);
      Type spr = q(0) / SBPR0;
      Type R = recbh(spr);
      sprg(j) = spr;
      Pg(j) = R * spr;
      Yg(j) = R * q(1);
      lxg(j) = log(x + Type(1e-12));
      rhog(j, 0) = q(2) / q(nq_out - 1);
      for(int i = 0; i < spm_nidx; i++) rhog(j, 1 + i) = q(3 + i) / q(nq_out - 1);
    }
    // EB/SB curve at depletion P: Nadaraya-Watson kernel smoother in log P over
    // the grid nodes (Gaussian kernel, bandwidth 0.25 in log P). Smooth in P,
    // which matters because P is a random effect: with linear interpolation
    // the Laplace log-determinant jumps when a state crosses a node. Flat
    // beyond the grid.
    vector<Type> lPg(nx);
    for(int j = 0; j < nx; j++) lPg(j) = log(Pg(j) + Type(1e-4));
    Type bw2 = Type(2) * Type(0.25) * Type(0.25);
    auto rho_at = [&](Type P, int col) {
      Type lp = log(P + Type(1e-4));
      Type sw = Type(0), sr = Type(0);
      for(int j = 0; j < nx; j++) {
        Type w = exp(-(lp - lPg(j)) * (lp - lPg(j)) / bw2);
        sw += w; sr += w * rhog(j, col);
      }
      return Type(sr / (sw + Type(1e-300)));
    };

    // reference F: Fmsy (max equilibrium yield) or Fspr, by taped damped
    // Newton iterations on log F from the best grid node. The iterations are
    // on the tape (TMB tapes once), so the reference point and its
    // derivatives follow the parameters; converged, the derivative is the
    // implicit-function derivative.
    Type uref;
    Type dd = Type(1e-3);
    if(spm_shape == 0) {
      Type ybest = Yg(1);
      uref = lxg(1);
      for(int j = 2; j < nx; j++) {
        uref  = CppAD::CondExpGt(Yg(j), ybest, lxg(j), uref);
        ybest = CppAD::CondExpGt(Yg(j), ybest, Yg(j), ybest);
      }
      auto Yu = [&](Type u) { vector<Type> q = prq(exp(u)); return Type(recbh(q(0) / SBPR0) * q(1)); };
      for(int it = 0; it < spm_newton_it; it++) {
        Type y0 = Yu(uref), yp = Yu(uref + dd), ym = Yu(uref - dd);
        Type g = (yp - ym) / (Type(2) * dd);
        Type H = (yp - Type(2) * y0 + ym) / (dd * dd);
        Type step = CppAD::CondExpLt(H, Type(0), -g / (H - Type(1e-30)),
                                     CppAD::CondExpGt(g, Type(0), Type(0.5), Type(-0.5)));
        step = CppAD::CondExpGt(step, Type(0.5), Type(0.5), step);
        step = CppAD::CondExpLt(step, Type(-0.5), Type(-0.5), step);
        uref += step;
      }
    } else {
      // largest grid F with SPR above the target, then Newton on log SPR
      uref = lxg(1);
      for(int j = 1; j < nx; j++) uref = CppAD::CondExpGt(sprg(j), spm_spr_tgt, lxg(j), uref);
      Type ltgt = log(spm_spr_tgt);
      auto gu = [&](Type u) { vector<Type> q = prq(exp(u)); return Type(log(q(0) / SBPR0) - ltgt); };
      for(int it = 0; it < spm_newton_it; it++) {
        Type g0 = gu(uref);
        Type gd = (gu(uref + dd) - gu(uref - dd)) / (Type(2) * dd);
        Type step = -g0 / CppAD::CondExpLt(gd, Type(-1e-8), gd, Type(-1e-8));
        step = CppAD::CondExpGt(step, Type(1), Type(1), step);
        step = CppAD::CondExpLt(step, Type(-1), Type(-1), step);
        uref += step;
      }
    }
    Type Fref = exp(uref);
    vector<Type> qr = prq(Fref);
    Type spr_ref = qr(0) / SBPR0;
    Type BmsyK_pr = recbh(spr_ref) * spr_ref;          // SB_ref / SB0
    Type Hmsy_pr = qr(1) / qr(2);                      // yield / exploitable biomass at Fref
    Type Yref = recbh(spr_ref) * qr(1);                // equilibrium yield per R0 at Fref
    // Pella-Tomlinson shape from BmsyK = m^(-1/(m-1)): Newton on v = log m
    Type lbk_ref = (spm_bmsyk > Type(0)) ? Type(log(spm_bmsyk)) : Type(log(BmsyK_pr));
    auto fv = [&](Type v) {
      Type ev = exp(v) - Type(1);
      Type av = CppAD::CondExpLt(fabs(v), Type(1e-6), Type(1) - v / Type(2), v / (ev + Type(1e-30)));
      return Type(lbk_ref + av);
    };
    Type vm = log(Type(2));
    for(int it = 0; it < spm_newton_it; it++) {
      Type f0 = fv(vm);
      Type fd = (fv(vm + dd) - fv(vm - dd)) / (Type(2) * dd);
      Type step = -f0 / CppAD::CondExpLt(fd, Type(-1e-8), fd, Type(-1e-8));
      step = CppAD::CondExpGt(step, Type(1), Type(1), step);
      step = CppAD::CondExpLt(step, Type(-1), Type(-1), step);
      vm += step;
    }
    Type m_pr = exp(vm);
    Type r_pr = m_pr * Hmsy_pr;                        // Hmsy = r / m (JABBA3); Hmsy = Fmsy with mean-year EB
    Type sdb_fl = Type(0.8) * spm_sigR * pow(Type(2) * Hmsy_pr, Type(0.45));   // FishLife SD at r_Schaefer = 2 Hmsy

    // ---- state-space surplus production model (JABBA3 Selex, SB state) ----
    Type K = exp(spm_logK);
    Type r, m;
    if(spm_r_sd > Type(0)) { r = exp(spm_logr); nll -= dnorm(spm_logr, log(r_pr), spm_r_sd, true); } else r = r_pr;
    if(spm_m_sd > Type(0)) { m = exp(spm_logm); nll -= dnorm(spm_logm, log(m_pr), spm_m_sd, true); } else m = m_pr;
    // K prior: on = 1 lognormal; on = 2 one-sided (JABBA3 "half8": flat below the
    // median, half-normal on log K above)
    if(spm_K_pr(2) == Type(1)) nll -= dnorm(spm_logK, log(spm_K_pr(0)), spm_K_pr(1), true);
    if(spm_K_pr(2) == Type(2)) {
      Type dk = smax_(spm_logK - log(spm_K_pr(0)), Type(0), Type(1e-4)) / spm_K_pr(1);
      nll += Type(0.5) * dk * dk;
    }
    if(spm_bk_pr(2) > Type(0)) nll -= dnorm(spm_lbk, log(spm_bk_pr(0)), spm_bk_pr(1), true);
    Type sdb = exp(spm_logsdb);
    if(spm_sdb_fl == 1) {
      if(spm_sdb_pr(2) > Type(0)) nll -= dnorm(spm_logsdb, log(sdb_fl), spm_sdb_pr(1), true);
      else sdb = sdb_fl;                               // penalised fit: SD fixed at the FishLife value
    } else if(spm_sdb_pr(2) > Type(0)) nll -= dnorm(spm_logsdb, log(spm_sdb_pr(0)), spm_sdb_pr(1), true);
    vector<Type> sdi = exp(spm_logsdi);
    for(int i = 0; i < spm_nidx; i++) {
      // JABBA: gamma(0.001, 0.001) on 1/sdi^2 (with Jacobian), soft upper bound 1
      nll -= Type(0.001) * (Type(-2) * spm_logsdi(i)) - Type(0.001) * exp(Type(-2) * spm_logsdi(i));
      Type pb = smax_(spm_logsdi(i), Type(0), Type(1e-4)) - Type(0);
      nll += Type(500) * pb * pb;
    }
    Type BmsyK_spm = pow(m, Type(-1) / (m - Type(1)));
    Type Hmsy_spm = r / (m - Type(1)) * (Type(1) - Type(1) / m);

    int NT = spm_nyear;
    vector<Type> P = exp(spm_logP);
    vector<Type> rem(NT), Pmean(NT), H(NT), estC(NT);
    for(int t = 0; t < NT; t++) rem(t) = rho_at(P(t), 0);
    // catch: observed, or estimated with error (log C_t = log C_obs - cdev_t, cdev ~ N(0, sdc))
    for(int t = 0; t < NT; t++) estC(t) = spm_catch(t) * ((spm_catch_err == 1) ? exp(spm_cdev(t)) : Type(1));
    if(spm_catch_err == 1) nll -= dnorm(spm_cdev, Type(0), spm_sdc, true).sum();
    Pmean(0) = spm_lbk;
    for(int t = 1; t < NT; t++) {
      Type Pl = P(spm_lagidx(t - 1));
      Type pm = P(t - 1) + r / (m - Type(1)) * Pl * (Type(1) - pow(Pl, m - Type(1))) - estC(t - 1) / (K * rem(t - 1));
      Pmean(t) = log(smax_(pm, Type(0.001), Type(1e-7)));
    }
    for(int t = 0; t < NT; t++) {
      nll -= dnorm(spm_logP(t), Pmean(t), sdb, true);
      Type pb = smin_(spm_logP(t) - log(spm_Pbound(0)), Type(0), Type(1e-4)) +
                smax_(spm_logP(t) - log(spm_Pbound(1)), Type(0), Type(1e-4));
      nll += Type(500) * pb * pb;
      H(t) = estC(t) / (K * P(t) * rem(t));
    }
    // optional random walk on the log harvest rate (SPiCT-like F process, as JABBA3 frw)
    Type sdf = exp(spm_logsdf);
    if(spm_frw_t.size() > 0) {
      for(int k = 0; k < spm_frw_t.size(); k++) {
        int t = spm_frw_t(k);
        nll -= dnorm(log(H(t)), log(H(t - 1)), sdf, true);
      }
      if(spm_sdf_pr(2) > Type(0)) nll -= dnorm(spm_logsdf, log(spm_sdf_pr(0)), spm_sdf_pr(1), true);
    }
    vector<Type> procdev = spm_logP - Pmean;

    // indices: log I ~ N(log q + log K + log P + log rho_i(P), sqrt(se^2 + sdi^2))
    int nobs = spm_logI.size();
    vector<Type> spm_Imean(nobs);
    for(int k = 0; k < nobs; k++) {
      int t = spm_obs_t(k), i = spm_obs_i(k);
      Type rho_i = (spm_idx_gear(i) == -1) ? Type(1) : Type(rho_at(P(t), 1 + i));
      spm_Imean(k) = spm_logq(i) + spm_logK + spm_logP(t) + log(rho_i);
      nll -= dnorm(spm_logI(k), spm_Imean(k), sqrt(spm_se2(k) + sdi(i) * sdi(i)), true);
    }

    // links between the length model and the SPM. The equilibrium lag of the
    // length model gives persistent link errors, so they are AR1 over the
    // length-data years (rho^gap for gaps): e_0 ~ N(0, sig),
    // e_y ~ N(rho_g e_(y-1), sig sqrt(1 - rho_g^2)), rho_g = rho^gap
    auto ar1_nll = [&](vector<Type> e, Type rho, Type sig) {
      Type out = -dnorm(e(0), Type(0), sig, true);
      for(int y = 1; y < e.size(); y++) {
        Type rg = rho;                       // rho^gap by products (pow() has no derivative at rho = 0)
        for(int k = 1; k < spm_lfd_t(y) - spm_lfd_t(y - 1); k++) rg *= rho;
        out -= dnorm(e(y), rg * e(y - 1), sig * sqrt(Type(1) - rg * rg + Type(1e-10)), true);
      }
      return out;
    };
    Type rho_link = Type(2) / (Type(1) + exp(-spm_trho)) - Type(1);
    // SPR link: depletion implied by the length model's SPR under steepness
    vector<Type> D_spr(nyear), e_spr(nyear);
    D_spr.setZero(); e_spr.setZero();
    for(int y = 0; y < nyear; y++) {
      Type D = (Type(4) * spm_h * spr_y(y) - (Type(1) - spm_h)) / (Type(5) * spm_h - Type(1));
      D_spr(y) = smax_(D, Type(1e-3), Type(1e-8));
      e_spr(y) = log(D_spr(y)) - spm_logP(spm_lfd_t(y));
    }
    if(spm_link_spr == 1) nll += ar1_nll(e_spr, rho_link, spm_sigma_link);

    Type MSY_spm = BmsyK_spm * K * Hmsy_spm * rho_at(BmsyK_spm, 0);
    vector<Type> BBmsy = P / BmsyK_spm;
    vector<Type> HHmsy = H / Hmsy_spm;
    vector<Type> spm_logBBmsy = log(BBmsy);
    vector<Type> spm_logHHmsy = log(HHmsy);
    vector<Type> spm_logB = spm_logP + spm_logK;
    Type spm_logMSY = log(MSY_spm);
    Type spm_r = r, spm_m = m, spm_K = K, spm_sdb = sdb;
    vector<Type> spm_q = exp(spm_logq);
    REPORT(SBPR0); REPORT(Fref); REPORT(Yref); REPORT(spr_ref); REPORT(BmsyK_pr); REPORT(Hmsy_pr); REPORT(m_pr); REPORT(r_pr); REPORT(sdb_fl);
    vector<Type> spm_sdi = sdi;
    REPORT(spm_r); REPORT(spm_m); REPORT(spm_K); REPORT(spm_sdb); REPORT(spm_sdi); REPORT(spm_q);
    REPORT(BmsyK_spm); REPORT(Hmsy_spm); REPORT(MSY_spm);
    REPORT(estC); REPORT(sdf); REPORT(P); REPORT(BBmsy); REPORT(H); REPORT(HHmsy); REPORT(rem); REPORT(Pmean); REPORT(procdev);
    REPORT(spm_Imean); REPORT(D_spr); REPORT(e_spr); REPORT(rho_link);
    vector<Type> spm_sA = sA;
    REPORT(Pg); REPORT(Yg); REPORT(rhog); REPORT(sprg); REPORT(spm_sA);
    ADREPORT(spm_logB); ADREPORT(spm_logBBmsy); ADREPORT(spm_logHHmsy); ADREPORT(spm_logMSY);
    ADREPORT(spm_logP);
  }

  // --- PER-RECRUIT REFERENCE POINTS (double mode only, on request) ---
  // Mirrors the R functions pr_flicc() / spr_flicc(): terminal-year averages
  // of joint selectivity (scaled to apical 1), M, maturity and weight;
  // numbers per recruit scaled to the first bin of the unfished state;
  // catch = N F: the numbers at length are integrated over the time spent in
  // each class (deaths / Z), so no Baranov factor; SBPR = sum(N * mat * exp(-Z * spawn) * w).
  if(isDouble<Type>::value && (brp_F.size() > 0 || brp_spr.size() > 0)) {
    int ny = brp_nyears < 1 ? 1 : (brp_nyears > nyear ? nyear : brp_nyears);
    int y0 = nyear - ny;
    vector<Type> sA(nlen), Ms(nlen), MA(nlen), matA(nlen), wtA(nlen);
    vector<Type> Ms_m(nlen), wtA_m(nlen), MA_m(nlen);
    sA.setZero(); Ms.setZero(); matA.setZero(); wtA.setZero(); Ms_m.setZero(); wtA_m.setZero();
    for(int l = 0; l < nlen; l++) {
      for(int y = y0; y < nyear; y++) {
        sA(l)   += sel_joint(l, y) / Type(ny);
        Ms(l)   += Mscaler(l, y)   / Type(ny);
        matA(l) += mat(l, y)       / Type(ny);
        wtA(l)  += wt(l, y)        / Type(ny);
        Ms_m(l) += Mscaler_m(l, y) / Type(ny);
        wtA_m(l) += wt_m(l, y)     / Type(ny);
      }
    }
    Type smax = Type(0);
    for(int l = 0; l < nlen; l++) if(sA(l) > smax) smax = sA(l);
    for(int l = 0; l < nlen; l++) {
      if(smax > Type(0)) sA(l) = sA(l) / smax;
      MA(l) = Mk * brp_k * Ms(l);
      MA_m(l) = Mk_m * (brp_k / k_ratio) * Ms_m(l);   // male annual M
    }

    // male numbers-at-length per recruit (nsex = 2), on the female time scale
    // (x k_ratio); same fishing at length as females
    auto popNm = [&](Type x, vector<Type>& Fl, vector<Type>& Zl) {
      for(int l = 0; l < nlen; l++) {
        Fl(l) = (brp_FMscale == 1) ? x * MA(l) * sA(l) : x * sA(l);
        Zl(l) = MA_m(l) + Fl(l);
      }
      matrix<Type> ZK(nlen + 1, ngtg);
      for(int g = 0; g < ngtg; g++) {
        for(int l = 0; l < nlen; l++) ZK(l, g) = MKMat_m(l, g) * Ms_m(l) + Fl(l) * k_ratio / brp_k;
        ZK(nlen, g) = MKMat_m(nlen, g) * Ms_m(nlen - 1) + Fl(nlen - 1) * k_ratio / brp_k;
      }
      vector<Type> NI = pop_len_gtg(LLB, gtgLinfs_m, ZK, recP);
      return vector<Type>(NI * k_ratio);
    };
    // per total recruit: females sex_ratio, males 1 - sex_ratio (nsex = 1: sex_ratio = 1)
    Type male_w = (nsex == 2) ? Type(Type(1) - sex_ratio) : Type(0);

    // numbers-at-length per recruit (unscaled) for fishing level x
    auto popN = [&](Type x, vector<Type>& Fl, vector<Type>& Zl) {
      for(int l = 0; l < nlen; l++) {
        Fl(l) = (brp_FMscale == 1) ? x * MA(l) * sA(l) : x * sA(l);
        Zl(l) = MA(l) + Fl(l);
      }
      vector<Type> NI(nlen);
      if(pop_model == 1) {
        vector<Type> Zk(nlen);
        for(int l = 0; l < nlen; l++) Zk(l) = Zl(l) / brp_k;
        NI = pop_len_glq(node, quad_wt, LLB, Zk, Galpha, Gbeta);
      } else {
        matrix<Type> ZK(nlen + 1, ngtg);
        for(int g = 0; g < ngtg; g++) {
          for(int l = 0; l < nlen; l++) {
            ZK(l, g) = MKMat(l, g) * Ms(l) + Fl(l) / brp_k;
          }
          ZK(nlen, g) = MKMat(nlen, g) * Ms(nlen - 1) + Fl(nlen - 1) / brp_k;
        }
        NI = pop_len_gtg(LLB, gtgLinfs, ZK, recP);
      }
      return NI;
    };

    vector<Type> Fl(nlen), Zl(nlen);
    vector<Type> NI0 = popN(Type(0), Fl, Zl);
    Type scale0 = NI0(0);
    Type brp_SBPR0 = Type(0);
    for(int l = 0; l < nlen; l++) {
      brp_SBPR0 += sex_ratio * NI0(l) / scale0 * matA(l) * exp(-Zl(l) * brp_spawn) * wtA(l);
    }
    auto sbpr = [&](Type x) {
      vector<Type> NI = popN(x, Fl, Zl);
      Type sb = Type(0);
      for(int l = 0; l < nlen; l++) {
        sb += sex_ratio * NI(l) / scale0 * matA(l) * exp(-Zl(l) * brp_spawn) * wtA(l);
      }
      return sb;
    };

    int nF = brp_F.size();
    matrix<Type> brp_N(nlen, nF), brp_C(nlen, nF), brp_Fl(nlen, nF);
    matrix<Type> brp_N_m(nlen, nF), brp_C_m(nlen, nF);   // males per male recruit (nsex = 2)
    brp_N_m.setZero(); brp_C_m.setZero();
    vector<Type> brp_SBPR(nF), brp_YPR(nF);
    for(int i = 0; i < nF; i++) {
      vector<Type> NI = popN(brp_F(i), Fl, Zl);
      Type sb = Type(0), yp = Type(0);
      for(int l = 0; l < nlen; l++) {
        Type n = NI(l) / scale0;
        Type c = n * Fl(l);
        brp_N(l, i) = n; brp_C(l, i) = c; brp_Fl(l, i) = Fl(l);
        sb += sex_ratio * n * matA(l) * exp(-Zl(l) * brp_spawn) * wtA(l);
        yp += sex_ratio * c * wtA(l);
      }
      if(nsex == 2) {
        vector<Type> Flm(nlen), Zlm(nlen);
        vector<Type> NIm = popNm(brp_F(i), Flm, Zlm);
        for(int l = 0; l < nlen; l++) {
          Type n = NIm(l) / scale0;
          Type c = n * Flm(l);
          brp_N_m(l, i) = n; brp_C_m(l, i) = c;
          yp += male_w * c * wtA_m(l);
        }
      }
      brp_SBPR(i) = sb; brp_YPR(i) = yp;   // SBPR: females; YPR: both sexes; per total recruit
    }

    // Fspr on log F by the Illinois method (bracketed regula falsi; SPR
    // decreases monotonically with F); -1 when the target is not reached
    // within (1e-8, brp_Fmax)
    int nS = brp_spr.size();
    vector<Type> brp_Fspr(nS);
    Type xlo = log(Type(1e-8)), xhi = log(brp_Fmax);
    Type glo0 = sbpr(exp(xlo)), ghi0 = sbpr(exp(xhi));
    for(int j = 0; j < nS; j++) {
      Type tgt = brp_spr(j) * brp_SBPR0;
      Type a = xlo, b = xhi, fa = glo0 - tgt, fb = ghi0 - tgt;
      if(fa < Type(0) || fb > Type(0)) {
        brp_Fspr(j) = Type(-1);
        continue;
      }
      int side = 0;
      Type c = a;
      for(int it = 0; it < 200; it++) {
        c = (a * fb - b * fa) / (fb - fa);
        Type fc = sbpr(exp(c)) - tgt;
        if(fabs(fc) <= Type(1e-14) * brp_SBPR0 || fabs(b - a) < Type(1e-12)) break;
        if(fc > Type(0)) {          // root lies above c
          a = c; fa = fc;
          if(side == -1) fb /= Type(2);
          side = -1;
        } else {
          b = c; fb = fc;
          if(side == 1) fa /= Type(2);
          side = 1;
        }
      }
      brp_Fspr(j) = exp(c);
    }

    // numbers-at-length per recruit at each Fspr (reference LFD, LBIspr)
    matrix<Type> brp_Nspr(nlen, nS), brp_Nspr_m(nlen, nS);
    brp_Nspr.setZero(); brp_Nspr_m.setZero();
    for(int j = 0; j < nS; j++) {
      if(brp_Fspr(j) > Type(0)) {
        vector<Type> NI = popN(brp_Fspr(j), Fl, Zl);
        for(int l = 0; l < nlen; l++) brp_Nspr(l, j) = NI(l) / scale0;
        if(nsex == 2) {
          vector<Type> Flm(nlen), Zlm(nlen);
          vector<Type> NIm = popNm(brp_Fspr(j), Flm, Zlm);
          for(int l = 0; l < nlen; l++) brp_Nspr_m(l, j) = NIm(l) / scale0;
        }
      }
    }
    REPORT(brp_Nspr_m);
    REPORT(brp_N_m);
    REPORT(brp_C_m);
    REPORT(MA_m);
    REPORT(wtA_m);
    REPORT(wtA);
    REPORT(matA);

    REPORT(brp_Nspr);
    REPORT(brp_SBPR0);
    REPORT(brp_SBPR);
    REPORT(brp_YPR);
    REPORT(brp_N);
    REPORT(brp_C);
    REPORT(brp_Fl);
    REPORT(brp_Fspr);
    REPORT(sA);
    REPORT(MA);
  }

  // Random walk on log F by gear: a penalty when log_Fk are fixed effects,
  // the random-effect density when settings$F_re = TRUE (Laplace)
  if(prior_sigmaF_use == 1 && nyear > 1) {
    for(int g = 0; g < ngear; g++) {
      if(gear_type(g) == 1) continue;      // surveys: no F (fixed), no random walk
      for(int y = 1; y < nyear; y++) {
        nll -= dnorm(log_Fk(y,g), log_Fk(y-1,g), sigmaF, true);
      }
    }

    // optional prior on sigmaF (prior_sigmaF_sd <= 0: none, e.g. when F is a
    // random effect and sigmaF is estimated from the marginal likelihood)
    if(obs_model == 3 && dm_type > 0 && dm_prior_sd > Type(0)) {
      for(int g = 0; g < ngear; g++) nll -= dnorm(log_dm(g), dm_prior_mu(g), dm_prior_sd, true);
    }
    if(prior_sigmaF_sd > Type(0)) {
      nll -= dnorm(log_sigmaF, prior_sigmaF_mean, prior_sigmaF_sd, true);
    }
  }

  // --- CONDITIONAL AGE AT LENGTH (ALK) ---
  // P(a | l) for a key of year y sampled at time tau: fish of age a are in length class l when their Linf lies in
  // [LB_l / g_a, LB_l+1 / g_a], g_a = 1 - exp(-k (a + tau - t0)); integrated over that interval (Gauss-Legendre)
  // of the Linf density (normal, CVL, truncated at +- maxsd) times the survival to length L = Linf g_a along the
  // length classes with the year's Z/K (LBSPR recursion). Selectivity at length cancels in P(a | l). Smooth in
  // k, Linf and CVL (interval limits by conditional expressions).
  Type nll_alk = Type(0);
  array<Type> alk_pexp(alk_obs.dim(0), alk_obs.dim(1), alk_obs.dim(2));   // expected age proportions by length
  alk_pexp.setZero();
  Type kgr = exp(log_k);
  if(alk_on == 1) {
    int na = alk_obs.dim(0), nkey = alk_obs.dim(2), nq = alk_gl_x.size();
    vector<Type> LB(nlen + 1);
    for(int l = 0; l < nlen; l++) LB(l) = LLB(l);
    LB(nlen) = (nlen > 1) ? Type(LLB(nlen - 1) + (LLB(nlen - 1) - LLB(nlen - 2))) : Type(LLB(0) + Type(1));
    Type sdL = CVL * Linf;
    Type Llo = Linf * (Type(1) - alk_maxsd * CVL);
    Type Lhi = Linf * (Type(1) + alk_maxsd * CVL);
    for(int kk = 0; kk < nkey; kk++) {
      int y = alk_year(kk);
      for(int l = 0; l < nlen; l++) {
        Type nobs = Type(0);
        for(int a = 0; a < na; a++) nobs += alk_obs(a, l, kk);
        if(nobs <= Type(0)) continue;
        vector<Type> C(na);
        for(int a = 0; a < na; a++) {
          Type ga = Type(1) - exp(-kgr * (alk_age(a) + alk_tau(kk) - t0g));
          ga = CppAD::CondExpGt(ga, Type(1e-4), ga, Type(1e-4));
          Type lo = LB(l) / ga, hi = LB(l + 1) / ga;
          lo = CppAD::CondExpGt(lo, Llo, lo, Llo);
          hi = CppAD::CondExpLt(hi, Lhi, hi, Lhi);
          Type w = hi - lo;
          w = CppAD::CondExpGt(w, Type(0), w, Type(0));
          Type acc = Type(0);
          for(int q = 0; q < nq; q++) {
            Type Li = lo + w * (alk_gl_x(q) + Type(1)) / Type(2);
            Type L  = Li * ga;
            Type ls = MKMat(0, 0) * Mscaler(0, y) * log((Li - LB(0)) / Li);
            for(int i = 0; i < l; i++) ls += ZK_y(i, y) * log((Li - LB(i + 1)) / (Li - LB(i)));
            ls += ZK_y(l, y) * log((Li - L) / (Li - LB(l)));
            acc += alk_gl_w(q) * dnorm(Li, Linf, sdL, false) * exp(ls);
          }
          C(a) = acc * w / Type(2);
        }
        // small floor on every age before normalising: finite gradients when a length class lies beyond the
        // Linf range (all C = 0), e.g. while Linf is low during the optimisation
        Type eps = Type(1e-20);
        Type Csum = C.sum() + Type(na) * eps;
        for(int a = 0; a < na; a++) {
          alk_pexp(a, l, kk) = (C(a) + eps) / Csum;
          if(alk_obs(a, l, kk) > Type(0)) nll_alk -= alk_obs(a, l, kk) * log((C(a) + eps) / Csum);
        }
      }
    }
    nll += nll_alk;
  }
  REPORT(nll_alk);
  REPORT(alk_pexp);
  REPORT(kgr);
  REPORT(t0g);

  // --- NORMAL PRIOR PENALTIES ---
  for(int i = 0; i < prior_mu.size(); i++) {
    Type x;
    int code = prior_code(i);

    if(code == 1) {
      x = log_Linf;
    } else if(code == 2) {
      x = log_Galpha;
    } else if(code == 3) {
      x = log_Mk;
    } else if(code == 4) {
      x = log_phi;
    } else if(code == 5) {
      x = log_Linf_m;
    } else if(code == 6) {
      x = log_Mk_m;
    } else if(code == 7) {
      x = log_Galpha_m;
    } else if(code == 9) {
      x = t0g;                                    // von Bertalanffy t0 (ALK growth): normal prior on the natural scale
    } else if(code >= 100 && code < 200) {
      int g = code - 100;
      x = log_Fk(0, g);
    } else if(code >= 200) {
      int j = code - 200;
      x = Sm(j);
    } else {
      error("Unknown prior_code");
    }

    nll -= dnorm(x, prior_mu(i), prior_sd(i), true);
  }

  vector<Type> log_spr_y = log(spr_y);
  vector<Type> log_epr_y = log(epr_y + Type(1e-30));

  REPORT(Sel);        // nlen x ngear
  REPORT(plen);       // nlen x nyear x ngear
  REPORT(epr_y);
  REPORT(spr_y);      // nyear x 1
  REPORT(catch_share); // nyear x ngear, predicted catch shares in weight
  REPORT(catch_share_n); // nyear x ngear, predicted catch shares in numbers
  REPORT(dm_alpha0);
  REPORT(dm_neff);
  REPORT(Fk);         // nyear x ngear on natural scale
  REPORT(sigmaF);
  REPORT(N_y);
  REPORT(plen_all_y);
  REPORT(sel_joint);
  REPORT(Fk_l);
  REPORT(Linf);
  REPORT(Galpha);
  REPORT(Gbeta);
  REPORT(Mk);
  REPORT(phi);
  REPORT(CVL);
  REPORT(gtgLinfs);   // GTG Linf groups at the estimates (pop_model = "gtg")
  REPORT(recP);
  REPORT(MKMat);      // (nlen + 1) x ngtg M/K used in GTG survival
  REPORT(nsex);
  REPORT(Linf_m);
  REPORT(Mk_m);
  REPORT(CVL_m);
  REPORT(gtgLinfs_m);
  REPORT(MKMat_m);
  REPORT(plen_sex);   // nlen x nyear x ngear x 2 (nsex = 2)
  REPORT(N_y_m);



  ADREPORT(log_Linf);
  ADREPORT(log_Galpha);
  ADREPORT(log_Mk);
  ADREPORT(log_sigmaF);
  if(nsex == 2) { ADREPORT(log_Linf_m); ADREPORT(log_Mk_m); ADREPORT(log_Galpha_m); }
  ADREPORT(log_spr_y);
  ADREPORT(log_epr_y);

  return nll;
}
