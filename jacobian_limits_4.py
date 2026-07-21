"""
NOTE 12, stage C-v2: the corrected clipped-shadow limit and the exact m*.

THEORY (after stage-C-v1's LINE-only integral (14.57%) overdrew Region s<-1):
 limit real wall in the Gaussian window =
   [main cubic branch t in (-1,1]: envelope tau_+(s) of t^2 - 2t^3],
   [LEFT pin layer t = -1 - x/d, odd d: (S_-(x), R_-(x)) = (-5+2e^x, 3-2e^x) = LINE r=-s-2],
   [RIGHT pin layer t = 1 + y/d: p ~ -1 - delta(3 + e^y): tau -> s (line r = s)],
 with algebraic crossing certificates:
   corner (LINE vs main): -2(t+1)^2 (t-1) >= 0 for t<=1, equality t=1 -> corner -> (-1,-1);
   partA integral exact: int phi(s) Phi(s) ds = Phi(a)^2/2.

 m* = int_{-1}^{1/3} phi(s) Phi(-s-2) ds + Phi(-1/sigma)^2 / 2,   sigma = 1.5.

PREDICTIONS, LOCKED:
 (D1) identities: (t^2-2t^3) - (-(2t-3t^2)-2) = -2 (t+1)^2 (t-1) EXACT (sympy);
      int_-inf^a phi Phi = Phi(a)^2/2: numerical check to 12 digits.
 (D2) m* in (7.4%, 8.9%); and m* <= m(d) for all computed d (finite-d cones
      sit above the limit: layer corrections are one-sided).
 (D3) polished-root masses m(d), d = 9..121 odd: monotone decreasing for
      d >= 11, tail m(d) - m* ~ C d^-alpha with alpha in (0.35, 1.2);
      m(121) in (m*, m* + 0.30pp).
 (D4) envelope audit grid: at s = -2, d=9 env in (-2.25,-2.0); at s=0 env in
      (-2.05,-1.90); both converge toward LINE/right-layer limits from above.
"""
import json, math
import numpy as np
import mpmath as mp
import sympy as sp
from sympy import symbols, expand, integrate, Rational as R

mp.mp.dps = 80
w = symbols("w")
out = {}
sig = 1.5

def seed(d):
    return expand(2*w - 3*w**2 + w*(1-w)*(w**(sp.Integer(d)-2) - R(6, d*(d+1))))

# ---------- (D1) algebraic certificates ----------
print("="*84); print("(D1) algebraic certificates"); print("-"*84)
t = symbols("t")
cert = sp.expand((t**2 - 2*t**3) - (-(2*t - 3*t**2) - 2))
print(f"  main-vs-LINE difference:  {sp.factor(cert)}   [predict -2*(t+1)^2*(t-1)]")
d1ok = sp.factor(cert) == sp.factor(-2*(t+1)**2*(t-1))
phi = lambda x: mp.exp(-x**2/(2*sig**2))/(sig*mp.sqrt(2*mp.pi))
PHI = lambda x: mp.mpf("0.5")*(1+mp.erf(x/(sig*mp.sqrt(2))))
lhs = mp.quad(lambda s: phi(s)*PHI(s), [-mp.inf, -1])
rhs = PHI(-1)**2/2
print(f"  partA identity: quad = {mp.nstr(lhs, 15)}  Phi(-1)^2/2 = {mp.nstr(rhs, 15)}  agree: {abs(lhs-rhs) < mp.mpf('1e-14')}")
out["D1"] = {"diff_factor": str(sp.factor(cert)), "partA": mp.nstr(rhs, 15), "partA_ok": bool(abs(lhs-rhs) < mp.mpf('1e-14'))}

# ---------- (D2) m* ----------
print("="*84); print("(D2) exact limit mass m*"); print("-"*84)
partB = mp.quad(lambda s: phi(s)*PHI(-s-2), [-1, mp.mpf(1)/3])
mstar = rhs + partB
print(f"  partB = {mp.nstr(partB, 15)}")
print(f"  m*    = {mp.nstr(mstar, 15)} = {mp.nstr(100*mstar, 10)} %   [predict (7.4, 8.9)]")
out["D2"] = {"partB": mp.nstr(partB, 15), "mstar": mp.nstr(mstar, 15), "mstar_pct": mp.nstr(100*mstar, 10)}

# ---------- (D3/D4) polished-root envelope masses ----------
print("="*84); print("(D3) polished-root masses m(d), d = 9..121"); print("-"*84)
def poly_data(d):
    p = seed(d); Phi = expand(integrate(p, w))
    pcn = np.array([float(c) for c in sp.Poly(p, w).all_coeffs()])
    tcn = np.array([float(c) for c in sp.Poly(sp.expand(w*p - Phi), w).all_coeffs()])
    coef = [mp.mpf(str(c)) for c in sp.Poly(p, w).all_coeffs()]
    return p, Phi, pcn, tcn, coef

def env_at(sv, pcn, tcn, coef):
    pc = pcn.copy(); pc[-1] -= sv
    cand = np.roots(pc)
    coef_s = list(coef); coef_s[-1] = coef_s[-1] - mp.mpf(float(sv))   # shift: p(z) - sv
    f = lambda z: mp.polyval(coef_s, z)                                # LEDGER: was 'coef' -> polished to wrong roots, invisible at s=0
    df = lambda z: mp.polyval([c*(len(coef_s)-1-i0) for i0, c in enumerate(coef_s[:-1])], z)
    taus = []
    for r0v in cand:
        if abs(r0v.imag) > 1e-2:      # genuinely complex: skip (Newton stays off-axis)
            continue
        z = mp.mpc(mp.mpf(float(r0v.real)), mp.mpf(float(r0v.imag)))
        for _ in range(25):
            try: dz = f(z)/df(z)
            except ZeroDivisionError: break
            z -= dz
            if abs(dz) < mp.mpf(10)**-50: break
        if abs(mp.im(z)) < mp.mpf(10)**-10 and abs(mp.re(z)) < 6:
            tt = mp.re(z)
            taus.append(float(np.polyval(tcn, float(tt))))
    return min(taus) if taus else np.inf

def masses_polished(dmax=121):
    S = np.linspace(-10, 5, 2601)
    pdf = np.exp(-S**2/(2*sig**2))/(sig*math.sqrt(2*math.pi))
    from math import erf as merf
    ms = {}
    for d in list(range(9, 20, 2)) + [23, 29, 35, 45, 61, 91, 121]:
        if d > dmax: break
        p, Phi, pcn, tcn, coef = poly_data(d)
        Ren = np.array([env_at(sv, pcn, tcn, coef) for sv in S])
        cdf = np.array([0.5*(1+merf(rv/(sig*math.sqrt(2)))) if np.isfinite(rv) else 1.0 for rv in Ren])
        ms[d] = float(np.trapezoid(pdf*cdf, S))
        print(f"  m({d:3d}) = {100*ms[d]:.4f}%   m-m* = {100*(ms[d]-float(mstar)):+.4f}pp", flush=True)
    return ms, S

ms, Sgrid = masses_polished()
mono = all(ms[d1] >= ms[d2] for d1, d2 in zip(sorted(ms), sorted(ms)[1:]) if d1 >= 11)
print(f"  monotone decreasing for d >= 11: {mono}")
ds = np.array([d for d in sorted(ms) if d >= 13], float)
dv = np.array([ms[d] for d in ds]) - float(mstar)
ok = dv > 0
lg = np.polyfit(np.log(ds[ok]), np.log(dv[ok]), 1)
print(f"  tail exponent alpha = {-lg[0]:.3f}   [predict (0.35,1.2)];  m(121) - m* = {100*(ms[121]-float(mstar)):+.4f}pp  [predict < 0.30pp]")
out["D3"] = {"m": {str(d): ms[d] for d in ms}, "alpha": float(-lg[0]), "monotone": bool(mono)}

# ---------- (D4) envelope audit ----------
print("="*84); print("(D4) envelope audit grid"); print("-"*84)
for d in (9, 29, 61):
    p, Phi, pcn, tcn, coef = poly_data(d)
    e2 = env_at(-2.0, pcn, tcn, coef); e0 = env_at(0.0, pcn, tcn, coef)
    print(f"  d={d:3d}: env(-2) = {e2:+.5f}   env(0) = {e0:+.5f}   [d=9 predict (-2.25,-2.0)/(-2.05,-1.90); -> -2 from above]")
out["D4"] = "see prints"
json.dump(out, open("limits_stageD.json", "w"), indent=1)
print("\nsaved limits_stageD.json")
