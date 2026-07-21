"""
NOTE 12, stage A-v2 (fixes) + stage B (EXACT ENVELOPE MASSES - the decisive test).

Stage-A fixes (honesty ledger):
 (a) MC coefficient slip [-2,0,1,0] -> [-2,1,0,0] (tau = t^2 - 2t^3).
 (b) P3 interval [-0.99,1.5] was the WRONG regime: param convergence holds on
     |t| <= rho < 1 only (cusp transients live in the boundary layer |t| -> 1).
 (c) three_pt used la.solve on a 3x2 matrix -> lstsq.
 (d) contact pairs printed min/max (d=9 tuple order had flipped; values were right).

REVISED THEORY (after I = 16.2888% vs census 8.726% - naive shadow FAIL):
 missed_n = cubic_shadow INTERSECT un-rescued, where the rescue envelope is the
 second real branch of p_d(t)=s (and far-field ones): r_env(s) = min over real
 roots t of tau_d(t). As d->inf (odd d), t2-cusp migrates to -1^+, rescue recedes,
 cone mass crawls UP toward the full shadow mass I_inf = 16.2888%.

PREDICTIONS, LOCKED:
 (A2) rho-fit s(d)=s*+B R^d via lstsq: s* in (-1.02,-0.94); corner(11) s in (-0.950,-0.935).
 (A3) sup on [-0.95,0.95] of ||(p_d,tau_d)-(p_inf,tau_inf)||: DECREASES with d,
      ratio sup(5)/sup(9) in (1.2, 3.0).
 (B1) exact envelope masses m(d) (no map-MC): m(5,7,9) within 0.15% ABS of censuses
      (8.69, 8.74, 8.726); m(9) in [8.60%, 8.80%].
 (B2) m(d) odd d to 29: tail APPROACHES I_inf = 16.2888% from below, eventually
      increasing; m(29) in (10%, 14%); m(d) monotone for d >= 13.
 (B3) fixed MC of I_inf agrees with quad to 3 sigma: I in (16.25%, 16.35%).
"""
import json, math
import numpy as np
import mpmath as mp
import sympy as sp
from sympy import symbols, expand, integrate, diff, Rational as R

w = symbols("w")
mp.mp.dps = 60
out = {}

def seed(d):
    return expand(2*w - 3*w**2 + w*(1-w)*(w**(sp.Integer(d)-2) - R(6, d*(d+1))))

# ---------- A2: corner fit (lstsq) ----------
print("="*84); print("(A2) crunode corner migration fit"); print("-"*84)
s5, s7, s9 = -0.8818834118879952, -0.9094571136757867, -0.928716289105
best = None
for Rg in np.linspace(1e-4, 0.99995, 60000):
    A = np.array([[1, Rg**5], [1, Rg**7], [1, Rg**9]])
    sol, res, *_ = np.linalg.lstsq(A, np.array([s5, s7, s9]), rcond=None)
    sstar, B = sol
    err = abs(sstar + B*Rg**5 - s5) + abs(sstar + B*Rg**7 - s7) + abs(sstar + B*Rg**9 - s9)
    if best is None or err < best[0]: best = (err, sstar, B, Rg)
err, sstar, B, Rg = best
s11 = sstar + B*Rg**11
print(f"  fit s(d)=s*+B R^d:  s*={sstar:.5f}, B={B:.5f}, R={Rg:.5f}, design-residual {err:.2e}")
print(f"  s* in (-1.02,-0.94): {-1.02 < sstar < -0.94}    corner(11): s11 = {s11:.5f} in (-0.950,-0.935): {-0.950 < s11 < -0.935}")
out["A2"] = {"sstar": float(sstar), "B": float(B), "R": float(Rg), "s11": float(s11), "resid": float(err)}

# ---------- A3: param convergence on |t|<=0.95 ----------
print("="*84); print("(A3) param convergence, |t| <= 0.95 only (honest regime)"); print("-"*84)
p_inf = 2*w - 3*w**2; Phi_inf = expand(integrate(p_inf, w))
p_inff = sp.lambdify(w, p_inf, "numpy"); phi_inff = sp.lambdify(w, Phi_inf, "numpy")
TT = np.linspace(-0.95, 0.95, 30001)
sups = {}
for d in (5, 7, 9):
    p = seed(d); Phi = expand(integrate(p, w))
    pf = sp.lambdify(w, p, "numpy"); phif = sp.lambdify(w, Phi, "numpy")
    dmax = float(np.max(np.hypot(pf(TT)-p_inff(TT), (TT*pf(TT)-phif(TT)) - (TT*p_inff(TT)-phi_inff(TT)))))
    sups[d] = dmax
    print(f"  d={d}: sup = {dmax:.6f}")
print(f"  decreasing: {sups[5]>sups[7]>sups[9]}   ratio 5/9 = {sups[5]/sups[9]:.3f}  [predict (1.2,3.0)]")
out["A3"] = {"sups": sups, "ratio": sups[5]/sups[9]}

# ---------- B3: fixed MC of the shadow mass ----------
print("="*84); print("(B3) I_inf fixed Monte Carlo vs quad"); print("-"*84)
sig = 1.5
def shadow_mask(S, Rv):
    t = (2 + np.sqrt(np.clip(4 - 12*S, 0, None)))/6
    env = -2*t**3 + t**2                      # FIXED coefficients
    return (S <= 1/3) & (Rv < env)
rng = np.random.default_rng(7)
NMC = 4_000_000
S = rng.normal(0, sig, NMC); Rv = rng.normal(0, sig, NMC)
mc = float(np.mean(shadow_mask(S, Rv))); mc_err = 3*math.sqrt(mc*(1-mc)/NMC)
Iq = 0.162887856959152
print(f"  I_inf MC = {100*mc:.4f}% +- {100*mc_err:.4f}%   quad = {100*Iq:.4f}%   agree(3sig): {abs(mc-Iq) < mc_err}")
out["B3"] = {"I_MC": mc, "I_MC_err": mc_err, "I_quad": Iq}

# ---------- B1/B2: EXACT envelope masses per chamber ----------
print("="*84); print("(B1/B2) exact envelope masses m(d)  -- NO map Monte Carlo --"); print("-"*84)
def envelope_mass(d, ngrid=4001):
    p = seed(d); Phi = expand(integrate(p, w))
    pc = np.array([float(c) for c in sp.Poly(p, w).all_coeffs()])            # p (desc), len d+1
    # tau(w) = w p(w) - Phi(w), coefficients once, padded to common length:
    tau = sp.Poly(sp.expand(w*p - Phi), w).all_coeffs()
    tc = np.array([float(c) for c in tau])                                   # len d+2?  w p has deg d+1
    S = np.linspace(-12, 6, ngrid)
    Ren = np.full(ngrid, np.inf)
    for k, sv in enumerate(S):
        coef = pc.copy(); coef[-1] -= sv                                     # p(w) - s = 0
        roots = np.roots(coef)
        treal = roots[np.abs(roots.imag) < 1e-9].real
        if len(treal):
            Ren[k] = np.polyval(tc, treal).min()                             # min tau over real t with p(t)=s
        else:
            Ren[k] = np.inf   # no real critical point -> h monotone -> never missed
    from math import erf
    cdf = 0.5*(1+np.vectorize(erf)(Ren/(sig*math.sqrt(2))))
    pdf = np.exp(-S**2/(2*sig**2))/(sig*math.sqrt(2*math.pi))
    return float(np.trapezoid(pdf*cdf, S)), S, Ren

ms = {}
for d in (5, 7, 9, 11, 13, 15, 17, 19, 23, 29):
    m, Sg, Ren = envelope_mass(d)
    ms[d] = m
    print(f"  d={d:2d} (n={d+1:2d}): exact envelope mass m = {100*m:.4f}%", flush=True)
cen = {5: 8.6900, 7: 8.7400, 9: 8.7260}
for d in (5,7,9):
    dev = abs(100*ms[d]-cen[d])
    print(f"     vs census d={d}: |m - census| = {dev:.4f}%  [predict <= 0.15%]")
print(f"  approach: m(29) = {100*ms[29]:.4f}%  in (10%,14%): {10 < 100*ms[29] < 14}")
tail = [ms[d] for d in (13,15,17,19,23,29)]
print(f"  monotone for d>=13: {all(tail[i] < tail[i+1] for i in range(len(tail)-1))};  all < I_inf: {all(v < Iq for v in ms.values())}")
out["B"] = {"m": {str(d): ms[d] for d in ms}, "I_inf": Iq,
            "dev": {str(d): abs(100*ms[d]-cen[d]) for d in cen}}
json.dump(out, open("limits_stageB.json", "w"), indent=1)
print("\nsaved limits_stageB.json")
