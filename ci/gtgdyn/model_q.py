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
