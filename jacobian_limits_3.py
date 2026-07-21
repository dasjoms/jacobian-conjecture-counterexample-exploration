"""
NOTE 12, stage C: THE BOUNDARY-LAYER THEOREM (lim锥 cone, exact).

After B1 CONFIRMED (envelope masses reproduce censuses 8.5908/8.7687/8.7736 vs
8.69/8.74/8.726 without any map MC) and B2 REFUTED the naive approach to the
16.29% shadow (tail m(d) DECREASES from d=11: 8.7755 -> 8.6634 @ d=29),
the layer analysis says:
  limit real wall in window = cubic param curve t in [-1,1]  (pins: p_d(1)=-1,
  Phi_d(1)=0 give the fixed point (-1,-1) EXACTLY for all d; p_d(-1)->-3)
  PLUS the left-layer segment  (S_-(x),R_-(x)) = (-5+2 e^x, 3-2 e^x):  r = -s-2.
  ==> LIMIT CONE: r < min(...) with envelope = LINE r=-s-2 for s in (-3,1/3],
      right-layer line r=s for s<-3 (both beat cubic there), negligible tails.
  ==> m* = int_{-inf}^{-3} phi Phi(s) + int_{-3}^{1/3} phi Phi(-s-2).
  Corner = main-env vs LINE crossing: (t+1)^2 (t-1) = 0 -> corner -> (-1,-1)
  exactly (rho-fit of finite corners: s*=-0.9733, consistent crawl).

PREDICTIONS, LOCKED:
 (C1) m* in (8.35%, 8.75%); m(29)=8.6634 within [m*, m*+0.15%].
 (C2) tail m(d)-m* for odd d=13..45 decays ~1/d^alpha with alpha in (0.55, 1.15)
      (least-squares log-log fit, using m* from C1).
 (C3) LEFT-LINE CERTIFICATE: layer functions: p_d(-1-x/d) - (-5+2 e^x) -> 0 and
      tau_d(-1-x/d) - (3-2e^x) -> 0 as d->inf: numerically max |.| over x in
      [0, 3] shrinks (d=9,15,25,45); slope dr = -ds EXACT on the segment:
      R_-(x) = -S_-(x) - 2 (symbolic identity check).
 (C4) cusp2 migration: s_c2(d) -> -3, r_c2(d) INCREASES (collect d=5..15 odd);
      corner contacts t2(d): model t2 = -1 - (ln c)/d with c -> 2.0..2.2:
      print c_eff(d) = exp(d*(-1-t2)): trend toward ~2.
 (C5) prediction lock for round n=12 (d=11): corner s in (-0.950,-0.937);
      census in [m(11)-0.13%, m(11)+0.13%] = [8.63%, 8.89%]; cusp2-r(11) > 1.55.
"""
import json, math
import numpy as np
import mpmath as mp
import sympy as sp
from sympy import symbols, expand, integrate, Rational as R

mp.mp.dps = 60
w = symbols("w")
out = {}
sig = 1.5

def seed(d):
    return expand(2*w - 3*w**2 + w*(1-w)*(w**(sp.Integer(d)-2) - R(6, d*(d+1))))

# ---------- (C1) the exact limit integral ----------
print("="*84); print("(C1) m* = integral of the clipped-shadow limit"); print("-"*84)
phi = lambda x: mp.exp(-x**2/(2*sig**2))/(sig*mp.sqrt(2*mp.pi))
PHI = lambda x: mp.mpf("0.5")*(1+mp.erf(x/(sig*mp.sqrt(2))))
part1 = mp.quad(lambda s: phi(s)*PHI(s), [-mp.inf, -3])
part2 = mp.quad(lambda s: phi(s)*PHI(-s-2), [-3, mp.mpf(1)/3])
mstar = part1 + part2
print(f"  part1 (s<-3, env=s)        = {mp.nstr(part1, 12)}")
print(f"  part2 (-3..1/3, env=-s-2)  = {mp.nstr(part2, 12)}")
print(f"  m* = {mp.nstr(mstar, 15)} = {mp.nstr(100*mstar, 10)} %   [predict (8.35,8.75)]")
# finite-d masses from stage B (envelope data, exact recompute here for extension)
def envelope_mass(d, ngrid=4001):
    p = seed(d); Phi = expand(integrate(p, w))
    pc = np.array([float(c) for c in sp.Poly(p, w).all_coeffs()])
    tc = np.array([float(c) for c in sp.Poly(sp.expand(w*p - Phi), w).all_coeffs()])
    S = np.linspace(-12, 6, ngrid)
    Ren = np.full(ngrid, np.inf)
    for k, sv in enumerate(S):
        coef = pc.copy(); coef[-1] -= sv
        roots = np.roots(coef)
        treal = roots[np.abs(roots.imag) < 1e-9].real
        if len(treal): Ren[k] = np.polyval(tc, treal).min()
    cdf = 0.5*(1+np.vectorize(math.erf)(Ren/(sig*math.sqrt(2))))
    pdf = np.exp(-S**2/(2*sig**2))/(sig*math.sqrt(2*math.pi))
    return float(np.trapezoid(pdf*cdf, S))
ms = {}
for d in (5, 7, 9, 11, 13, 15, 17, 19, 23, 29, 35, 45):
    ms[d] = envelope_mass(d)
    print(f"  m({d:2d}) = {100*ms[d]:.4f}%   m-m* = {100*(ms[d]-mstar):+.4f}pp   [within m*..m*+0.15: {float(mstar) <= ms[d] <= float(mstar)+0.0015}]", flush=True)
out["C1"] = {"mstar": mp.nstr(mstar, 15), "mstar_pct": mp.nstr(100*mstar, 10),
             "part1": mp.nstr(part1, 12), "part2": mp.nstr(part2, 12),
             "m": {str(d): ms[d] for d in ms}}

# ---------- (C2) tail decay rate ----------
print("="*84); print("(C2) tail decay  m(d) - m*  ~  d^(-alpha)"); print("-"*84)
ds = np.array([13, 15, 17, 19, 23, 29, 35, 45], float)
dv = np.array([ms[d] for d in ds]) - float(mstar)
ok = dv > 0
lg = np.polyfit(np.log(ds[ok]), np.log(dv[ok]), 1)
print(f"  log-log slope alpha = {-lg[0]:.3f}   [predict (0.55,1.15)]   (resid rms {np.std(np.polyval(lg, np.log(ds[ok])) - np.log(dv[ok])):.3f})")
out["C2"] = {"alpha": float(-lg[0]), "tail6": list(dv)}

# ---------- (C3) layer certificates ----------
print("="*84); print("(C3) left-layer certificates:  p_d(-1-x/d) -> -5+2 e^x,  tau -> 3-2 e^x"); print("-"*84)
XX = np.linspace(0, 3, 401)
for d in (9, 15, 25, 45):
    p = seed(d); Phi = expand(integrate(p, w))
    pf = sp.lambdify(w, p, "numpy"); phif = sp.lambdify(w, Phi, "numpy")
    T = -1 - XX/d
    dp = np.max(np.abs(pf(T) - (-5 + 2*np.exp(XX))))
    dtau = np.max(np.abs((T*pf(T) - phif(T)) - (3 - 2*np.exp(XX))))
    print(f"  d={d:2d}: max|Dp| = {dp:.6f},  max|Dtau| = {dtau:.6f}   over x in [0,3]  (should shrink)")
x = sp.symbols("x")
print(f"  line identity: 3-2e^x = -(-5+2e^x) - 2  ->  {sp.simplify((3-2*sp.exp(x)) + (-5+2*sp.exp(x)) + 2) == 0}  (r = -s-2 EXACT)")
out["C3"] = "see prints"

# ---------- (C4) cusp2 migration + contact layer-model ----------
print("="*84); print("(C4) cusp2 (t~ -0.89) migration + corner contact layer model"); print("-"*84)
mp.mp.dps = 120
p_inf = 2*w - 3*w**2
for d in range(5, 16, 2):
    p = seed(d); pp = sp.diff(p, w); Phi = expand(integrate(p, w))
    coef = [mp.mpf(str(c)) for c in sp.Poly(pp, w).all_coeffs()]
    f = lambda z: mp.polyval(coef, z)
    roots = []
    for s0 in [mp.mpf(-2)+mp.mpf(i)/40 for i in range(121)]:
        try:
            z = mp.findroot(f, s0, tol=mp.mpf(10)**-90, maxsteps=80)
            if abs(mp.im(z)) < mp.mpf(10)**-40 and all(abs(z-q) > mp.mpf(10)**-15 for q in roots):
                roots.append(mp.re(z))
        except Exception: pass
    roots.sort()
    t2 = roots[0]
    pd_ = sp.lambdify(w, p, "numpy"); phid_ = sp.lambdify(w, Phi, "numpy")
    sc2 = float(pd_(float(t2))); rc2 = float(float(t2)*pd_(float(t2)) - phid_(float(t2)))
    print(f"  d={d:2d}: t2 = {mp.nstr(t2, 18)}   cusp2 (s, r) = ({sc2:+.6f}, {rc2:+.6f})   [s->-3? r increasing: {d==5 or rc2>1.5526:}]")
    out.setdefault("C4", {})[str(d)] = {"t2": mp.nstr(t2, 40), "s_c2": sc2, "r_c2": rc2}
# contact t2 layer model c_eff = exp(d*(-1-t2))
con = {5: -1.271214, 7: -1.155970, 9: -1.108918}
print("  corner contact layer model t2 = -1 - (ln c)/d,  c_eff = e^{d(-1-t2)}:")
for d, t2v in con.items():
    print(f"    d={d}: c_eff = {math.exp(d*(-1-t2v)):.4f}   [trend -> 2.0..2.2]")
out["C4"]["c_eff"] = {str(d): math.exp(d*(-1-t2v)) for d, t2v in con.items()}

# ---------- (C5) round-12 locks ----------
print("="*84); print("(C5) LOCKS for the n=12 round (d=11)"); print("-"*84)
m11 = ms[11]
print(f"  corner(11): s in (-0.950, -0.937)   [rho-fit -0.94217]")
print(f"  exact envelope m(11) = {100*m11:.4f}%   census(11) lock: [8.63%, 8.89%]")
print(f"  cusp2(11): r > 1.55 (printed above)")
out["C5"] = {"m11": m11, "census_lock": [8.63, 8.89], "corner_lock": [-0.950, -0.937]}
json.dump(out, open("limits_stageC.json", "w"), indent=1)
print("\nsaved limits_stageC.json")
