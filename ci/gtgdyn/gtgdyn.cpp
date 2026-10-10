// gtg.dyn trial: dynamic, length-indexed growth-type-group model (standalone TMB).
//
// No age dimension. For growth group g, age at length boundary j is
// A_gj = -log(1 - L_j / Linf_g) / K (data). Numbers in bin j at time t:
//
//   N_gj(t) = int_0^{D_gj} R(b) S_gj(t - tau) exp(-M tau - sel_j [CF(t) - CF(t - tau)]) dtau,
//   b       = t - tau - A_gj                                  (recruitment time)
//   S_gj(s) = exp(-M A_gj - sum_{i<j} sel_i [CF(b + A_g,i+1) - CF(b + A_gi)])
//
// CF = cumulative F (piecewise linear in time; F constant before the first and
// after the last model year). D_gj is the bin residence time (Tcap for the bin
// containing Linf_g). The integral uses fixed quadrature nodes (data), so no
// time point depends on a parameter and the AD tape is valid everywhere.
// With constant F and R this reduces to the equilibrium GTG recursion.
//
// pop = 0: equilibrium GTG (current FLicc 'gtg', continuous form); 1: gtg.dyn.

#include <TMB.hpp>

template<class Type>
Type cfun(double t, const vector<Type> &cumF, const vector<Type> &F, int nF) {
  if (t < 0.0) return F(0) * t;
  if (t >= nF) return cumF(nF) + F(nF - 1) * (t - nF);
  int k = (int) std::floor(t);
  if (k < 0 || k >= nF) Rf_error("cfun index %d (t = %f, nF = %d)", k, t, nF);
  return cumF(k) + F(k) * (t - k);
}

template<class Type>
Type logrfun(double b, const vector<Type> &logR, int nF) {
  double u = b - 0.5;                      // knots at year mid-points
  if (u <= 0.0) return logR(0);
  if (u >= nF - 1) return logR(nF - 1);
  int k = (int) std::floor(u);
  if (k < 0 || k + 1 >= nF) Rf_error("logrfun index %d (b = %f)", k, b);
  double fr = u - k;
  return (1.0 - fr) * logR(k) + fr * logR(k + 1);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(pop);              // 0 eq, 1 dyn
  DATA_MATRIX(obs);               // ny x nb counts
  DATA_IVECTOR(kd);               // model-year index of each data year
  DATA_VECTOR(tq);                // observation times within year
  DATA_VECTOR(lmid);
  DATA_VECTOR(mw);                // maturity x weight at bin mid
  DATA_VECTOR(wg);                // GTG weights
  DATA_MATRIX(A);                 // ng x (nb+1) age at boundary (large if unreachable)
  DATA_MATRIX(dT);                // ng x nb residence time (0 unreachable)
  DATA_IMATRIX(reach);            // 0 unreachable, 1 finite bin, 2 bin containing Linf_g
  DATA_IVECTOR(node_g);           // quadrature nodes (dyn)
  DATA_IVECTOR(node_j);
  DATA_VECTOR(node_tau);
  DATA_VECTOR(node_w);
  DATA_SCALAR(db);                // birth-time grid step (<= 0: exact, no grid)
  DATA_SCALAR(b0);                // first grid point (H constant below)
  DATA_INTEGER(nbg);              // number of grid points
  DATA_SCALAR(M);
  DATA_INTEGER(sel_type);         // 0 logistic, 1 dome (double normal)
  DATA_INTEGER(rmode);            // 0 constant R, 1 iid, 2 random walk
  DATA_VECTOR(prior_sigF);        // mu (log), sd, on
  DATA_VECTOR(prior_sigR);

  PARAMETER_VECTOR(logF);         // nF
  PARAMETER_VECTOR(theta);        // selectivity
  PARAMETER_VECTOR(logR);         // nF (mapped off if rmode 0)
  PARAMETER(log_sigF);
  PARAMETER(log_sigR);

  int nF = logF.size(), ny = obs.rows(), nb = obs.cols(), ng = wg.size();
  int nn = node_g.size(), nq = tq.size();
  if (isDouble<Type>::value) {
    if (kd.size() != ny) Rf_error("kd size %d != ny %d", (int) kd.size(), ny);
    if (lmid.size() != nb || mw.size() != nb) Rf_error("lmid/mw size != nb %d", nb);
    if (A.rows() != ng || A.cols() != nb + 1) Rf_error("A dims %d x %d, expected %d x %d", (int) A.rows(), (int) A.cols(), ng, nb + 1);
    if (dT.rows() != ng || dT.cols() != nb || reach.rows() != ng || reach.cols() != nb) Rf_error("dT/reach dims");
    if (logR.size() != nF) Rf_error("logR size %d != nF %d", (int) logR.size(), nF);
    for (int y = 0; y < ny; y++) if (kd(y) < 0 || kd(y) >= nF) Rf_error("kd(%d) = %d out of range", y, kd(y));
    for (int n = 0; n < nn; n++) if (node_g(n) < 0 || node_g(n) >= ng || node_j(n) < 0 || node_j(n) >= nb) Rf_error("node %d out of range", n);
    if (sel_type == 1 && theta.size() < 3) Rf_error("theta size %d", (int) theta.size());
  }
  Type nll = 0;

  // selectivity
  vector<Type> sel(nb);
  for (int j = 0; j < nb; j++) {
    if (sel_type == 0) {
      sel(j) = Type(1) / (Type(1) + exp(-log(Type(19)) * (lmid(j) - theta(0)) / exp(theta(1))));
    } else {
      Type sd = CppAD::CondExpLt(Type(lmid(j)), theta(0), exp(theta(1)), exp(theta(2)));
      Type z = (lmid(j) - theta(0)) / sd;
      sel(j) = exp(Type(-0.5) * z * z);
    }
  }

  vector<Type> F = exp(logF);
  vector<Type> cumF(nF + 1);
  cumF(0) = 0;
  for (int k = 0; k < nF; k++) cumF(k + 1) = cumF(k) + F(k);

  // equilibrium numbers per unit recruitment (summed over GTG) at F = f
  auto eqN = [&](Type f) {
    vector<Type> N(nb); N.setZero();
    for (int g = 0; g < ng; g++) {
      Type E = 0;
      for (int j = 0; j < nb; j++) {
        if (reach(g, j) == 0) break;
        Type Z = M + f * sel(j);
        Type frac = (reach(g, j) == 2) ? Type(1) : Type(1) - exp(-Z * dT(g, j));
        N(j) += wg(g) * exp(-E) * frac / Z;
        if (reach(g, j) == 2) break;
        E += Z * dT(g, j);
      }
    }
    return N;
  };

  // cumulative fishing exposure H(g, m, j) on reaching boundary j for a cohort
  // born at grid time b0 + m db (running sum over bins; computed once)
  double dbd = asDouble(db), b0d = asDouble(b0);
  bool use_grid = (pop == 1) && (dbd > 0);
  array<Type> H(use_grid ? ng : 1, use_grid ? nbg : 1, nb + 1);
  if (use_grid) {
    for (int g = 0; g < ng; g++) for (int m = 0; m < nbg; m++) {
      double b = b0d + dbd * m;
      H(g, m, 0) = 0;
      Type cprev = cfun(b, cumF, F, nF);
      for (int i = 0; i < nb; i++) {
        double Ai1 = asDouble(A(g, i + 1));
        if (Ai1 > 1e9) { for (int ii = i + 1; ii <= nb; ii++) H(g, m, ii) = H(g, m, i); break; }
        Type cnext = cfun(b + Ai1, cumF, F, nF);
        H(g, m, i + 1) = H(g, m, i) + sel(i) * (cnext - cprev);
        cprev = cnext;
      }
    }
  }

  // dynamic numbers (summed over GTG) at time t
  auto dynN = [&](double t) {
    vector<Type> N(nb); N.setZero();
    Type cft = cfun(t, cumF, F, nF);
    for (int n = 0; n < nn; n++) {
      int g = node_g(n), j = node_j(n);
      double tau = asDouble(node_tau(n)), Aj = asDouble(A(g, j));
      double s = t - tau, b = s - Aj;
      Type E = M * Aj;
      if (use_grid) {
        double u = (b - b0d) / dbd;
        int m = (int) std::floor(u);
        if (m < 0) m = 0;
        if (m > nbg - 2) m = nbg - 2;
        double fr = u - m; if (fr < 0) fr = 0; if (fr > 1) fr = 1;
        E += (1.0 - fr) * H(g, m, j) + fr * H(g, m + 1, j);
      } else {
        Type cprev = cfun(b, cumF, F, nF);
        for (int i = 0; i < j; i++) {
          Type cnext = cfun(b + asDouble(A(g, i + 1)), cumF, F, nF);
          E += sel(i) * (cnext - cprev);
          cprev = cnext;
        }
      }
      Type within = M * tau + sel(j) * (cft - cfun(s, cumF, F, nF));
      Type R = (rmode == 0) ? Type(1) : exp(logrfun(b, logR, nF));
      N(j) += wg(g) * node_w(n) * R * exp(-E - within);
    }
    return N;
  };

  // predicted catch proportions
  matrix<Type> pred(ny, nb);
  for (int y = 0; y < ny; y++) {
    vector<Type> N(nb);
    if (pop == 0) {
      N = eqN(F(kd(y)));
    } else {
      N.setZero();
      for (int q = 0; q < nq; q++) N += dynN(kd(y) + asDouble(tq(q)));
      N /= Type(nq);
    }
    vector<Type> C = N * sel * F(kd(y));
    Type tot = C.sum();
    for (int j = 0; j < nb; j++) {
      pred(y, j) = C(j) / tot;
      nll -= obs(y, j) * log(pred(y, j) + Type(1e-12));
    }
  }

  // random walk on log F
  Type sigF = exp(log_sigF);
  for (int k = 1; k < nF; k++) nll -= dnorm(logF(k), logF(k - 1), sigF, true);
  if (prior_sigF(2) > 0) nll -= dnorm(log_sigF, prior_sigF(0), prior_sigF(1), true);

  // recruitment
  Type sigR = exp(log_sigR);
  if (rmode == 1) {
    for (int k = 0; k < nF; k++) nll -= dnorm(logR(k), Type(0), sigR, true);
  } else if (rmode == 2) {
    nll -= dnorm(logR(0), Type(0), Type(1), true);
    for (int k = 1; k < nF; k++) nll -= dnorm(logR(k), logR(k - 1), sigR, true);
  }
  if (rmode > 0 && prior_sigR(2) > 0) nll -= dnorm(log_sigR, prior_sigR(0), prior_sigR(1), true);

  REPORT(pred);
  REPORT(sel);
  REPORT(F);

  // derived quantities (double evaluation only, keeps the tape small)
  if (isDouble<Type>::value) {
    vector<Type> mwv(nb);
    for (int j = 0; j < nb; j++) mwv(j) = mw(j);
    Type sb0 = (eqN(Type(0)) * mwv).sum();
    vector<Type> spr(ny), ssb_rel(ny);
    for (int y = 0; y < ny; y++) {
      spr(y) = (eqN(F(kd(y))) * mwv).sum() / sb0;
      ssb_rel(y) = (pop == 1) ? (dynN(kd(y) + 0.5) * mwv).sum() / sb0 : spr(y);
    }
    REPORT(spr);
    REPORT(ssb_rel);
  }
  return nll;
}
