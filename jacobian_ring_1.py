"""
NOTE 13 ("THE RINGING CONSTANT"), stage A: the exact 1/d layer corrections.

The tower's missed-cone mass m(d) rings around m* = 8.64446465377901%:
  overshoot (m(9)  = 8.7736)  ->  zero-cross (~d ~ 45-60)  ->  undershoot
  (m(121) = 8.6329 < m*)  ->  locked rebound m(161) > m(121)  ->  home at m*.
That is TWO decay scales fighting: a +A/d term (left pin layer) and a
negative, slower term from the right pin layer (root t_d(s) = 1 + c/d with
c solving  c(4 + e^c) = d(|s|-1)  =>  c ~ ln d - ln ln d).

TODAY'S JOB: derive symbolically, audit numerically at 100 digits.

  (A1) LEFT LAYER THEOREM (machine-algebra certificate, symbolic d):
       with eps = 1/d, t = -1 - x*eps, odd d, x = x0 + x1*eps + x2*eps^2:
       tau_L(s; eps) = -s - 2 + A_L(s)*eps + B_L(s)*eps^2 + O(eps^3)
       A_L(s) = (s+5) - 5*ln((s+5)/2)            [both sides symbolic]
       Prediction: sp.simplify(A_L - candidate) == 0 EXACTLY, and B_L > 0
       on (-5, 1/3] (plotted + sampled, not claimed as proof).
  (A2) Numerical audits: exact root solves p_d(t) = s (100-digit Newton),
       d = 21, 61, 121, 201;  s in {-2, -1.2, -0.5, 0, 0.3}:
       |d*(tau_exact - tau_limit) - A_L| -> 0 like 1/d (left), and the
       right-layer model residual (implied c vs exact) small & shrinking.
  (A3) THE MASS MODEL (quadrature, 50 digits, sigma = 1.5):
       dm_L(d) = (1/d) * INT_{-1}^{oo} phi(s) phi(-s-2) A_L(s) ds      > 0
       dm_R(d) = (1/d) * INT_{-oo}^{-1} phi(s)^2 * [c(s;d)*(s+1) + (c-1)e^c/d] ds < 0
       c(s; d) solves c*(4+e^c) = d*(|s|-1).
       Model m_mod(d) = m* + dm_L + dm_R  vs  measured series (stage D of
       note 12): predict max |m_mod - m| < 1.5e-5 for d >= 29  (0.0015pp),
       better than 8e-6 for d >= 61.
  (A4) PREDICTIONS LOCKED HERE for stage B (exact envelope masses):
       m(131), m(141), m(151), m(161), m(171), m(181), m(201) intervals,
       min location d* in (121, 161), rebound m(161) > m(121) = 8.6329.
"""
import json, math
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 100
out = {}

# ---------------------------------------------------------------- (A1) LEFT
print("=" * 88); print("(A1) LEFT LAYER, symbolic eps-expansion (odd d)"); print("-" * 88)
x, Y, s = sp.symbols("x Y s")
S5 = (s + 5) / 2                       # e^{x0}
L0 = sp.log(S5)                        # x0

def pow_series(numer_shift, order):
    """(1+xY)^(1/Y + numer_shift) expanded in Y to O(Y^order), symbolic x."""
    lser = sp.series(sp.log(1 + x*Y), Y, 0, order).removeO().expand()/Y
    e    = sp.series(sp.exp(lser), Y, 0, order).removeO().expand()
    fac  = sp.series((1 + x*Y)**numer_shift, Y, 0, order).removeO().expand()
    return sp.expand(e * fac)

ORD = 3
E1   = pow_series(0, ORD)              # (1+xY)^(1/Y)
Em2  = sp.expand(pow_series(0, ORD) * sp.series((1+x*Y)**(-2), Y, 0, ORD).removeO())  # ^(1/Y-2)
E1p1 = sp.expand(pow_series(1, ORD))   # (1+xY)^(1/Y+1)

t  = -1 - x*Y
Yo1  = sp.series(Y/(1+Y), Y, 0, ORD+1).removeO().expand()       # Y/(1+Y) as series
cd   = sp.expand(6*Y*Yo1)                                        # c_d = 6Y^2/(1+Y)
p  = 2*t - 3*t**2 + t*(1 - t)*((-Em2) - cd)          # odd d: t^(d-2) = -(1+xY)^(1/Y-2)
Phi  = t**2 - t**3 + Y*(-1)*E1 - E1p1*Yo1 - cd*(t**2/2 - t**3/3)

x0 = sp.symbols("x0"); x1 = sp.symbols("x1"); x2 = sp.symbols("x2")
xsub = x0 + x1*Y + x2*Y**2

E0 = sp.symbols("E0")
def Yseries(e, order=ORD):
    """substitute x -> xsub inside e, series-expand all exp() atoms in Y to
    O(Y^order), mark exp(x0) as E0, expand fully."""
    e = sp.expand(e.subs(x, xsub))
    for a in list(e.atoms(sp.exp)):
        arg  = sp.expand(a.args[0])
        arg0 = arg.subs(Y, 0)
        rest = sp.expand(arg - arg0)
        ser  = sp.series(sp.exp(rest), Y, 0, order).removeO().expand()
        mark = E0 if arg0 == x0 else sp.exp(arg0)
        e = sp.expand(e.subs(a, mark*ser))
    return sp.expand(e)

def coeff_Y(e, k):
    return sp.expand(e).coeff(Y, k)

pexp = Yseries(p - s)
eq0  = coeff_Y(pexp, 0)
eq1  = coeff_Y(pexp, 1)
eq2  = coeff_Y(pexp, 2)
sol1 = sp.solve(sp.Eq(eq1, 0), x1)[0]
sol2 = sp.solve(sp.Eq(eq2.subs(x1, sol1), 0), x2)[0]
print("x0 =", L0)
print("x1 =", sp.simplify(sol1.subs(E0, S5).subs(x0, L0)))
print("x2 =", sp.simplify(sol2.subs(E0, S5).subs(x0, L0)))

tauexp = Yseries(t*s - Phi)
tau0 = sp.simplify(coeff_Y(tauexp, 0).subs(E0, S5).subs(x0, L0))
tau1 = sp.simplify(coeff_Y(tauexp, 1).subs(E0, S5).subs(x0, L0))
tau2 = sp.simplify(coeff_Y(tauexp, 2).subs(E0, S5).subs(x0, L0))
print("\ntau_L series:")
print("  O(1)    =", tau0, "   [predict -s-2]")
print("  O(eps)  =", tau1)
print("  O(eps^2)=", tau2)

cand = (s + 5) - 5*sp.log((s+5)/2)
AL_cert = sp.simplify(tau1 - cand) == 0
print(f"\n(A1) CERTIFICATE:  A_L(s) = (s+5) - 5*ln((s+5)/2)  EXACTLY:  {AL_cert}")
lim_cert = sp.simplify(tau0 + s + 2) == 0
print(f"(A1) CERTIFICATE:  limit line  tau0 = -s - 2        EXACTLY:  {lim_cert}")
out["A1"] = {"AL_cert": bool(AL_cert), "limit_cert": bool(lim_cert),
             "AL": str(tau1), "BL": str(tau2)}

# positivity sweep of A_L on the effective window (numeric, not a proof)
xs_ = np.linspace(-4.999, 8, 400)
AL_ = (xs_+5) - 5*np.log((xs_+5)/2)
print(f"  A_L(s) range on window: [{AL_.min():.6f}, {AL_.max():.6f}]  (min at s={xs_[np.argmin(AL_)]:.3f})")

# ---------------------------------------------------------------- (A2) AUDITS
print("=" * 88); print("(A2) 100-digit audits of both layers"); print("-" * 88)
def seed_co(d):
    """exact rational coefficient lists (mpf) of p_d and tau=w p-Phi."""
    w = sp.symbols("w")
    p_ = sp.expand(2*w - 3*w**2 + w*(1-w)*(w**(sp.Integer(d)-2) - sp.Rational(6, d*(d+1))))
    Ph = sp.expand(sp.integrate(p_, w))
    ta = sp.expand(w*p_ - Ph)
    pc = [mp.mpf(str(c)) for c in sp.Poly(p_, w).all_coeffs()]
    tc = [mp.mpf(str(c)) for c in sp.Poly(ta, w).all_coeffs()]
    return pc, tc

def exact_branch(d, sv, t0):
    pc, tc = seed_co(d)
    pcs = list(pc); pcs[-1] -= mp.mpf(sv)
    dps = [c*(len(pcs)-1-i) for i, c in enumerate(pcs[:-1])]
    f  = lambda z: mp.polyval(pcs, z)
    df = lambda z: mp.polyval(dps, z)
    z = mp.mpf(t0)
    for _ in range(60):
        dz = f(z)/df(z); z -= dz
        if abs(dz) < mp.mpf(10)**(-90): break
    return mp.re(z), mp.polyval(tc, z)

audL = []
for d in (21, 61, 121, 201):
    for sv in (-0.5, 0, 0.3):
        x0v = float(np.log((sv+5)/2)); t0 = -1 - x0v/d
        tt, tau = exact_branch(d, sv, t0)
        ALv  = (sv+5) - 5*np.log((sv+5)/2)
        err  = float(mp.fabs(d*(tau - mp.mpf(-sv-2)) - mp.mpf(ALv)))
        audL.append((d, sv, err))
        print(f"  LEFT  d={d:3d} s={sv:+.1f}:  d*dtau = {mp.nstr(d*(tau+mp.mpf(sv)+2), 12)}  vs A_L={ALv:.10f}   |diff|={err:.3e}")

print()
audR = []
for d in (21, 61, 121, 201):
    for sv in (-1.2, -2, -3):
        # model c
        g = lambda c: c*(4+mp.e**c) - d*(abs(sv)-1)
        c  = mp.findroot(g, 2.5)
        t0 = 1 + float(c)/d
        tt, tau = exact_branch(d, sv, t0)
        resid_t = float(mp.fabs(tt - (1 + c/d)))            # how good t = 1+c/d
        # model tau:  t s - Phi ; expand: tau_mod = s + [c(s+1)]/d + (c-1)e^c/d^2
        tau_mod = sv + float(c*(sv+1))/d + float((c-1)*mp.e**c)/d**2
        err  = float(mp.fabs(tau - mp.mpf(tau_mod)))
        audR.append((d, sv, float(tau-sv), err, float(c)))
        print(f"  RIGHT d={d:3d} s={sv:+.1f}:  tau-s = {float(tau-sv):+.8f}   model {tau_mod-sv:+.8f}   |mod err|={err:.3e}   c={float(c):.4f}   |t-(1+c/d)|={resid_t:.2e}")
out["A2"] = {"audLerr": audL, "audRerr": audR}

# ---------------------------------------------------------------- (A3) MODEL
print("=" * 88); print("(A3) mass model vs measured series"); print("-" * 88)
mp.mp.dps = 50
sig = mp.mpf("1.5")
phi = lambda v: mp.e**(-(v**2)/(2*sig**2))/(sig*mp.sqrt(2*mp.pi))
PHI = lambda v: mp.mpf("0.5")*(1 + mp.erf(v/(sig*mp.sqrt(2))))
mstar = PHI(-1)**2/2 + mp.quad(lambda u: phi(u)*PHI(-u-2), [-1, mp.inf])
print(f"  m* rebuilt = {mp.nstr(mstar, 18)}  (stored 0.08644464653779007893290411735)")

ALf = lambda u: (u+5) - 5*mp.log((u+5)/2)
A0  = mp.quad(lambda u: phi(u)*phi(-u-2)*ALf(u), [-1, mp.inf])
print(f"  A0 = INT phi(s)phi(-s-2)A_L ds = {mp.nstr(A0, 15)}")

B0  = mp.quad(lambda u: phi(u)**2*(mp.fabs(u)-1), [-mp.inf, -1])
print(f"  B0 = INT phi^2 (|s|-1) ds        = {mp.nstr(B0, 15)}  (c=1 limit coefficient)")

def c_of(u, d):
    lam = d*(mp.fabs(u)-1)
    if lam < mp.mpf("1e-12"): return mp.mpf(0)
    return mp.findroot(lambda c: c*(4+mp.e**c) - lam, 1.5)

def right_term(d):
    # includes the (c-1)e^c/d^2 piece: integrand [c(s+1) + (c-1)e^c/d]
    def integ(u):
        if mp.fabs(u)-1 < mp.mpf("1e-9"): return mp.mpf(0)
        c = c_of(u, d)
        return phi(u)**2*(c*(u+1) + (c-1)*mp.e**c/d)
    return mp.quad(integ, [-7, -1]) / d     # gaussian tail <-7 is ~1e-11

stage = json.load(open("limits_stageD.json"))
meas = {int(k): v for k, v in stage["D3"]["m"].items()}
meas[9]  = 0.08773584175648119      # stage-B value (same machinery)
meas[11] = 0.08755547157788661
meas[13] = 0.08736426802531863
meas[23] = 0.086790154156153
meas[29] = 0.08663354848051617

rows = {}
print("   d |     measured      model (A0/d + R(d))     diff(pp)")
for d in sorted(meas):
    Rt  = right_term(d)
    mdl = float(mstar) + float(A0)/d + float(Rt)
    df  = (float(meas[d]) - mdl)
    rows[d] = {"meas": meas[d], "model": mdl, "diff": df, "Rt": float(Rt), "Lft": float(A0)/d}
    print(f"  {d:3d} | {100*meas[d]:12.8f}  {100*mdl:12.8f}      {100*df:+9.5f}")
out["A3"] = {"mstar": str(mstar), "A0": str(A0), "B0": str(B0), "rows": rows}

# ---------------------------------------------------------------- (A4) LOCKS
print("=" * 88); print("(A4) PREDICTIONS LOCKED for stage B (exact masses)"); print("-" * 88)
locks = {}
for d in (131, 141, 151, 161, 171, 181, 201):
    Rt  = right_term(d)
    mdl = float(mstar) + float(A0)/d + float(Rt)
    # tolerance: last-model-error trend + 1/d^2 floor
    locks[d] = (mdl, mdl - 4.0e-5, mdl + 4.0e-5)
    print(f"  m({d:3d}) model = {100*mdl:.6f}%   LOCK: [{100*(mdl-4e-5):.6f}, {100*(mdl+4e-5):.6f}]   (m-m* = {100*(mdl-float(mstar)):+.6f}pp)")
locks["min_d"]    = "d* in (121, 161)"
locks["rebound"]  = "m(161) > m(121) = 8.63288597 (note-12 lock, re-locked)"
locks["min_val"]  = "m(d*) > 8.624pp?? see stage C"   # refined in stage C
out["A4"] = {str(k): v for k, v in locks.items() if k != "min_val"}
out["A4"]["rebound"] = locks["rebound"]; out["A4"]["min_d"] = locks["min_d"]
json.dump(out, open("ring_stageA.json", "w"), indent=1, default=str)
print("\nsaved ring_stageA.json")
