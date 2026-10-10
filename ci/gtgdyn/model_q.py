"""gtg.dyn, quadrature form: all evaluation times are data (TMB-friendly).

N_gj(t) = int_0^{D_gj} R(t - tau - A_gj) * S_gj(t - tau)
                      * exp(-M tau - sel_j [CF(t) - CF(t - tau)]) dtau
S_gj(s) = exp(-M A_gj - sum_{i<j} sel_i [CF(s - A_gj + A_g,i+1) - CF(s - A_gj + A_gi)])
D_gj = A_g,j+1 - A_gj (bin residence time), or Tcap for the bin containing Linf_g.
The integral uses composite Gauss-Legendre (panels of width <= h, ngl nodes).
"""
import numpy as np
from common import *

def gl(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return (x + 1) / 2, w / 2

class DynQ:
    def __init__(self, Lg, wg, nq=2, ngl=2, h=2.0, Tcap=25.0):
        self.Lg, self.wg, self.nq = Lg, wg, nq
        ng = len(Lg)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = 1 - LB[None, :] / Lg[:, None]
            A = np.where(r > 0, -np.log(np.where(r > 0, r, 1)) / K, np.inf)
        self.A = A
        self.reach = np.isfinite(A[:, :-1])
        D = A[:, 1:] - A[:, :-1]
        D = np.where(self.reach & ~np.isfinite(D), Tcap, D)
        self.D = np.where(self.reach, D, 0.0)
        x, w = gl(ngl)
        # quadrature nodes per (g, j): lists of (tau, weight)
        self.nodes = {}
        for g in range(ng):
            for j in range(NB):
                if not self.reach[g, j]: continue
                Dj = self.D[g, j]; npan = int(np.ceil(Dj / h))
                hp = Dj / npan
                tau = (np.arange(npan)[:, None] + x[None, :]).ravel() * hp
                wt_ = np.tile(w, npan) * hp
                self.nodes[(g, j)] = (tau, wt_)
        self.tq = (np.arange(nq) + 0.5) / nq
        self.MW = mat(LMID) * wt(LMID)

    @staticmethod
    def cf(t, edges, F):
        cf_e = np.r_[0.0, np.cumsum(F)]
        x = np.r_[edges[0] - 1000.0, edges, edges[-1] + 1.0, edges[-1] + 1001.0]
        y = np.r_[-1000.0 * F[0], cf_e, cf_e[-1] + 1000.0 * F[-1]]
        return np.interp(t, x, y)

    def N(self, t, years, F, sel, logR=None):
        edges = years.astype(float)
        cf = lambda u: self.cf(u, edges, F)
        ng = len(self.Lg)
        out = np.zeros((ng, NB))
        cft = cf(t)
        for (g, j), (tau, w) in self.nodes.items():
            Aj = self.A[g, j]
            s = t - tau; b = s - Aj
            Ai = self.A[g, :j]; Ai1 = self.A[g, 1:j + 1]
            E = M * Aj + (sel[:j][None, :] * (cf(b[:, None] + Ai1[None, :]) - cf(b[:, None] + Ai[None, :]))).sum(1)
            within = M * tau + sel[j] * (cft - cf(s))
            R = 1.0 if logR is None else np.exp(np.interp(b, edges + 0.5, logR))
            out[g, j] = np.sum(w * R * np.exp(-E - within))
        return out

    def catch(self, y, years, F, sel, logR=None):
        k = np.searchsorted(years, y)
        C = sum(self.wg @ self.N(y + tq, years, F, sel, logR) for tq in self.tq) / self.nq
        return C * F[k] * sel


class DynQG(DynQ):
    """As DynQ, but the cumulative fishing exposure H_g(b, j) of a cohort born at
    time b is computed once on a birth-time grid (step db) for all j by a running
    sum, and interpolated linearly in b for each quadrature node.
    Cost per evaluation: ng x nbgrid x nb (grid) + nodes x times (lookup)."""
    def __init__(self, Lg, wg, db=0.25, **kw):
        super().__init__(Lg, wg, **kw)
        self.db = db
        Afin = np.where(np.isfinite(self.A), self.A, np.nan)
        self.b0 = -np.nanmax(np.where(self.reach, self.A[:, :-1], np.nan)) - 1.0

    def H(self, years, F, sel, bmax):
        edges = years.astype(float)
        bg = self.b0 + self.db * np.arange(int(np.ceil((bmax - self.b0) / self.db)) + 2)
        ng = len(self.Lg)
        Hm = np.zeros((ng, len(bg), NB + 1))
        for g in range(ng):
            Ag = np.where(np.isfinite(self.A[g]), self.A[g], np.nan)
            cfv = self.cf(bg[:, None] + np.nan_to_num(Ag, nan=0.0)[None, :], edges, F)  # nbg x NB+1
            d = np.diff(cfv, axis=1) * sel[None, :]
            d = np.where(np.isfinite(Ag[1:])[None, :], d, 0.0)
            Hm[g, :, 1:] = np.cumsum(d, axis=1)
        return bg, Hm

    def N(self, t, years, F, sel, logR=None, _cache={}):
        edges = years.astype(float)
        cf = lambda u: self.cf(u, edges, F)
        # model time origin: years[0]
        bg, Hm = self.H(years, F, sel, t + 1.0)
        ng = len(self.Lg)
        out = np.zeros((ng, NB))
        cft = cf(t)
        for (g, j), (tau, w) in self.nodes.items():
            Aj = self.A[g, j]
            s = t - tau; b = s - Aj
            u = (b - bg[0]) / self.db
            m = np.clip(np.floor(u).astype(int), 0, len(bg) - 2)
            fr = np.clip(u - m, 0.0, 1.0)
            Hj = (1 - fr) * Hm[g, m, j] + fr * Hm[g, m + 1, j]
            E = M * Aj + Hj
            within = M * tau + sel[j] * (cft - cf(s))
            R = 1.0 if logR is None else np.exp(np.interp(b, edges + 0.5, logR))
            out[g, j] = np.sum(w * R * np.exp(-E - within))
        return out
