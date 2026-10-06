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

    for(int i = 0; i < n; i++) {
      Type d = L(i) - mode;
      if(L(i) < mode) {
        sel(i) = exp(-(d * d) / (two * lsd * lsd));
      } else {
        sel(i) = exp(-(d * d) / (two * rsd * rsd));
      }
    }

  } else if(sel_type == 3) {

    Type mode = exp(pars(0)) * Linf;
    Type sd   = exp(pars(1)) * Linf;
    Type two  = Type(2);

    for(int i = 0; i < n; i++) {
      Type d = L(i) - mode;
      sel(i) = exp(-(d * d) / (two * sd * sd));
    }

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
  // Per-recruit reference points at the estimates (evaluated only in double
  // mode, i.e. obj$report(), and only when requested; never during fitting)
  DATA_VECTOR(brp_F);        // apical F (or F/M if brp_FMscale) values
  DATA_VECTOR(brp_spr);      // SPR targets (proportions) for Fspr
  DATA_INTEGER(brp_nyears);  // terminal years averaged
  DATA_SCALAR(brp_k);        // von Bertalanffy k (annual scale)
  DATA_INTEGER(brp_FMscale); // 1: brp_F are F/M multipliers of M-at-length
  DATA_SCALAR(brp_spawn);    // spawning time fraction
  DATA_SCALAR(brp_Fmax);     // upper bound for the Fspr search
  DATA_INTEGER(brp_only);    // 1: per-recruit call, skip years before the
                             //    terminal brp_nyears (objective not used)

  PARAMETER(log_Linf);
  PARAMETER(log_Galpha);
  PARAMETER(log_Mk);
  PARAMETER_MATRIX(log_Fk);   // nyear x ngear
  PARAMETER_VECTOR(Sm);
  PARAMETER(log_phi);
  PARAMETER(log_sigmaF); // random walk F

  Type nll = 0.0;

  Type Linf   = exp(log_Linf);
  Type Galpha = exp(log_Galpha);
  Type Mk     = exp(log_Mk);
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
  matrix<Type> N_y(nlen, nyear);
  matrix<Type> plen_all_y(nlen, nyear);
  matrix<Type> sel_joint(nlen, nyear);
  matrix<Type> Fk_l(nlen, nyear);

  sel_joint.setZero();
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

    vector<Type> Z(nlen), N(nlen);
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
        F_len(l,g) = Fk(y,g) * Sel(l,g);
        Z(l) += F_len(l,g);
      }
    }

    if(pop_model == 1) {
      N = pop_len_glq(node, quad_wt, LLB, Z, Galpha, Gbeta);

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
        Cpred(l,g) = N(l) * F_len(l,g);
        pred_gear(g) += Cpred(l,g);
        obs_gear(g)  += obs(l,y,g);
      }
    }

    // predicted gear shares
    vector<Type> pred_pgear(ngear);
    Type pred_total = Type(0);

    pred_pgear.setZero();

    for(int g = 0; g < ngear; g++) {
      pred_total += pred_gear(g);
    }

    for(int g = 0; g < ngear; g++) {
      if(pred_total > Type(0)) {
        pred_pgear(g) = pred_gear(g) / pred_total;
      } else {
        pred_pgear(g) = Type(0);
      }
    }


    Type fmax = Type(0);
    Type M_at_fmax = M(0);

    for(int l = 0; l < nlen; l++) {

      Type fsum = Type(0);

      for(int g = 0; g < ngear; g++) {
        fsum += F_len(l,g);
      }

      Fk_l(l,y) = fsum;

      if(fsum > fmax) {
        fmax = fsum;
        M_at_fmax = M(l);
      }
    }

    if (FMpen_sd > Type(0)) {

      Type FM_y = fmax / (M_at_fmax + Type(1e-12)) + Type(1e-12);

      if (FM_y < FM_min) {
        nll += Type(1e6) * pow(log(FM_y) - log(FM_min), 2);
      }

      if (FM_y > FM_max) {
        nll += Type(1e6) * pow(log(FM_y) - log(FM_max), 2);
      }
    }


    for(int l = 0; l < nlen; l++) {
      if(fmax > Type(0)) {
        sel_joint(l,y) = Fk_l(l,y) / fmax;
      } else {
        sel_joint(l,y) = Type(0);
      }
    }

    // combined expected length composition across gears
    Type ctot = Type(0);
    for(int g = 0; g < ngear; g++) {
      ctot += catch_wt(y,g) * pred_gear(g);
    }

    for(int l = 0; l < nlen; l++) {
      Type ctot_l = Type(0);

      for(int g = 0; g < ngear; g++) {
        ctot_l += catch_wt(y,g) * Cpred(l,g);
      }

      if(ctot > Type(0)) {
        plen_all_y(l,y) = ctot_l / ctot;
      } else {
        plen_all_y(l,y) = Type(0);
      }
    }



    // soft constraint on predicted gear shares using relative catch
    for(int g = 0; g < (ngear - 1); g++) {
      if(pred_pgear(g) > Type(0) &&
         pred_pgear(ngear - 1) > Type(0) &&
         catch_wt(y,g) > Type(0) &&
         catch_wt(y,ngear - 1) > Type(0)) {

        Type logratio_pred  = log(pred_pgear(g) / pred_pgear(ngear - 1));
        Type logratio_catch = log(catch_wt(y,g) / catch_wt(y,ngear - 1));

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
    // Optional robustness for the upper tail (both off by default):
    //  - plus_bin >= 0 pools observed and predicted counts in bins
    //    l >= plus_bin into bin plus_bin (tail compression);
    //  - rob_eps > 0 mixes predicted proportions with a uniform:
    //    p = (1 - rob_eps) * p + rob_eps / nb.
    int nb = (plus_bin >= 0 && plus_bin < nlen) ? (plus_bin + 1) : nlen;

    for(int g = 0; g < ngear; g++) {
      vector<Type> o(nb), p(nb);
      o.setZero(); p.setZero();
      for(int l = 0; l < nlen; l++) {
        int k = (l < nb) ? l : (nb - 1);
        o(k) += obs(l,y,g);
        if(pred_gear(g) > Type(0)) p(k) += Cpred(l,g) / pred_gear(g);
      }
      if(rob_eps > Type(0)) {
        for(int k = 0; k < nb; k++) {
          p(k) = (Type(1) - rob_eps) * p(k) + rob_eps / Type(nb);
        }
      }

      if(obs_model == 1) {
        // Negative binomial on ESS-scaled counts
        for(int k = 0; k < nb; k++) {
          Type muk = obs_gear(g) * p(k) + Type(1e-12);
          if(k < nlen) mu(k,g) = muk;
          nll -= dnbinom_phi(o(k), muk, phi, 1);
        }

      } else if(obs_model == 2) {
        // Multinomial on within-gear compositions
        // obs is already ESS-scaled by lfdess()
        if(obs_gear(g) > Type(0) && pred_gear(g) > Type(0)) {
          for(int k = 0; k < nb; k++) {
            Type p_pred = p(k);
            if(p_pred < Type(1e-12)) p_pred = Type(1e-12);
            nll -= o(k) * log(p_pred);
          }
        }

      } else if(obs_model == 3) {
        // Dirichlet-multinomial on within-gear compositions
        // Here obs_gear(g) acts as fixed ESS / precision
        if(obs_gear(g) > Type(0) && pred_gear(g) > Type(0)) {
          Type alpha0 = obs_gear(g);
          nll -= lgamma(alpha0);
          nll += lgamma(obs_gear(g) + alpha0);
          for(int k = 0; k < nb; k++) {
            Type p_pred = p(k);
            if(p_pred < Type(1e-12)) p_pred = Type(1e-12);
            Type a = alpha0 * p_pred;
            nll -= lgamma(o(k) + a);
            nll += lgamma(a);
          }
        }

      } else {
        error("Unknown obs_model");
      }
    }


  }


  // --- PER-RECRUIT REFERENCE POINTS (double mode only, on request) ---
  // Mirrors the R functions pr_flicc() / spr_flicc(): terminal-year averages
  // of joint selectivity (scaled to apical 1), M, maturity and weight;
  // numbers per recruit scaled to the first bin of the unfished state;
  // Baranov catch; SBPR = sum(N * mat * exp(-Z * spawn) * w).
  if(isDouble<Type>::value && (brp_F.size() > 0 || brp_spr.size() > 0)) {
    int ny = brp_nyears < 1 ? 1 : (brp_nyears > nyear ? nyear : brp_nyears);
    int y0 = nyear - ny;
    vector<Type> sA(nlen), Ms(nlen), MA(nlen), matA(nlen), wtA(nlen);
    sA.setZero(); Ms.setZero(); matA.setZero(); wtA.setZero();
    for(int l = 0; l < nlen; l++) {
      for(int y = y0; y < nyear; y++) {
        sA(l)   += sel_joint(l, y) / Type(ny);
        Ms(l)   += Mscaler(l, y)   / Type(ny);
        matA(l) += mat(l, y)       / Type(ny);
        wtA(l)  += wt(l, y)        / Type(ny);
      }
    }
    Type smax = Type(0);
    for(int l = 0; l < nlen; l++) if(sA(l) > smax) smax = sA(l);
    for(int l = 0; l < nlen; l++) {
      if(smax > Type(0)) sA(l) = sA(l) / smax;
      MA(l) = Mk * brp_k * Ms(l);
    }

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
      brp_SBPR0 += NI0(l) / scale0 * matA(l) * exp(-Zl(l) * brp_spawn) * wtA(l);
    }
    auto sbpr = [&](Type x) {
      vector<Type> NI = popN(x, Fl, Zl);
      Type sb = Type(0);
      for(int l = 0; l < nlen; l++) {
        sb += NI(l) / scale0 * matA(l) * exp(-Zl(l) * brp_spawn) * wtA(l);
      }
      return sb;
    };

    int nF = brp_F.size();
    matrix<Type> brp_N(nlen, nF), brp_C(nlen, nF), brp_Fl(nlen, nF);
    vector<Type> brp_SBPR(nF), brp_YPR(nF);
    for(int i = 0; i < nF; i++) {
      vector<Type> NI = popN(brp_F(i), Fl, Zl);
      Type sb = Type(0), yp = Type(0);
      for(int l = 0; l < nlen; l++) {
        Type n = NI(l) / scale0;
        Type c = (Zl(l) > Type(0)) ? n * Fl(l) / Zl(l) * (Type(1) - exp(-Zl(l))) : Type(0);
        brp_N(l, i) = n; brp_C(l, i) = c; brp_Fl(l, i) = Fl(l);
        sb += n * matA(l) * exp(-Zl(l) * brp_spawn) * wtA(l);
        yp += c * wtA(l);
      }
      brp_SBPR(i) = sb; brp_YPR(i) = yp;
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
    matrix<Type> brp_Nspr(nlen, nS);
    brp_Nspr.setZero();
    for(int j = 0; j < nS; j++) {
      if(brp_Fspr(j) > Type(0)) {
        vector<Type> NI = popN(brp_Fspr(j), Fl, Zl);
        for(int l = 0; l < nlen; l++) brp_Nspr(l, j) = NI(l) / scale0;
      }
    }

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

  // Penalty on random walk F
  if(prior_sigmaF_use == 1 && nyear > 1) {
    for(int g = 0; g < ngear; g++) {
      for(int y = 1; y < nyear; y++) {
        nll -= dnorm(log_Fk(y,g), log_Fk(y-1,g), sigmaF, true);
      }
    }

    nll -= dnorm(log_sigmaF, prior_sigmaF_mean, prior_sigmaF_sd, true);
  }

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

  REPORT(Sel);        // nlen x ngear
  REPORT(plen);       // nlen x nyear x ngear
  REPORT(spr_y);      // nyear x 1
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



  ADREPORT(log_Linf);
  ADREPORT(log_Galpha);
  ADREPORT(log_Mk);
  ADREPORT(log_sigmaF);
  ADREPORT(log_spr_y);

  return nll;
}
