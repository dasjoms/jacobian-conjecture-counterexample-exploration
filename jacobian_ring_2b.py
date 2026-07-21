"""
NOTE 13, stage A2b: min-refine + LOCKS for the rebound (stage B checks them).

Model = GOLD model of note-13 stage A2 (validated: diffs +0.0011pp @ d=45
shrinking smoothly to +0.0004pp @ d=121).  ODD d only (cone chambers).
Tasks:
  - re-audit dm(d) for d in {45, 61, 91, 121} (recomputed, compared to run-1);
  - refine the minimum over odd d = 113..131 (predict d* in (121, 131));
  - quadratic fit on the neighborhood -> d*, m(d*);
  - LOCK exact masses: d in {131,141,151,161,171,181,201}:
      window = m_mod(d) + [bias(121)-drift, bias(45) as cap] style band; PLUS
      correlated-difference locks:
        m(161) > m(121) - 3e-6 (prob units)     [strong rebound lock]
        m(171) > m(161) - 1.5e-6                [slope lock]
saves ring_stageA2.json
"""
import json, math
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 50
sig = mp.mpf("1.5")
phi = lambda v: mp.e**(-(v**2)/(2*sig**2))/(sig*mp.sqrt(2*mp.pi))
PHI = lambda v: mp.mpf("0.5")*(1 + mp.erf(v/(sig*mp.sqrt(2))))
mstar = PHI(-1)**2/2 + mp.quad(lambda u: phi(u)*PHI(-u-2), [-1, mp.inf])

AL = lambda u: (u+5)*(1 - mp.log((u+5)/2))
BLs = sp.sympify(json.load(open("ring_probe_sym.json"))["BL"])
_s = sp.symbols("s")
BLsy = sp.lambdify(_s, BLs, "mpmath")
BL = lambda u: mp.mpf(BLsy(mp.mpf(u)))

_pcache = {}
def seed_co(d):
    if d in _pcache: return _pcache[d]
    w = sp.symbols("w")
    p_ = sp.expand(2*w - 3*w**2 + w*(1-w)*(w**(sp.Integer(d)-2) - sp.Rational(6, d*(d+1))))
    Ph = sp.expand(sp.integrate(p_, w)); ta = sp.expand(w*p_ - Ph)
    pc = [mp.mpf(str(c)) for c in sp.Poly(p_, w).all_coeffs()]
    tc = [mp.mpf(str(c)) for c in sp.Poly(ta, w).all_coeffs()]
    _pcache[d] = (pc, tc); return pc, tc

def c_of(u, d):
    lam = d*(mp.mpf(u) + 1)
    if abs(lam) < mp.mpf("1e-6"): return lam/5
    return mp.findroot(lambda c: c*(4+mp.e**c) - lam, 0.0)

def rightE(sv, d):
    pc, tc = seed_co(d)
    pcs = list(pc); pcs[-1] -= mp.mpf(sv)
    dps = [c*(len(pcs)-1-i) for i, c in enumerate(pcs[:-1])]
    f, df = (lambda z: mp.polyval(pcs, z)), (lambda z: mp.polyval(dps, z))
    nrm = lambda z: sum(abs(c)*abs(z)**k for k, c in enumerate(reversed(pcs)))
    disc = 4 - 12*float(sv)
    t_cub = (2 + math.sqrt(max(disc, 0.0)))/6
    for z0 in (1 + float(c_of(sv, d))/d, t_cub, 1.6):
        z = mp.mpf(z0)
        try:
            for _ in range(80):
                fz, dfz = f(z), df(z)
                dz = fz/dfz; z -= dz
                if abs(dz) < mp.mpf(10)**(-46)*(1+abs(z)): break
            if abs(f(z))/nrm(z) < mp.mpf(10)**(-42) and mp.re(z) > 0.5:
                return mp.polyval(tc, z)
        except (ZeroDivisionError, ValueError, OverflowError):
            continue
    # guaranteed bracket: f(1) = -1-s > 0, f(3) << 0  =>  manual bisect
    a, b, fa = mp.mpf(1), mp.mpf(3), f(mp.mpf(1))
    for _ in range(300):
        m = (a+b)/2; fm = f(m)
        if fa*fm > 0: a, fa = m, fm
        else: b = m
        if b - a < mp.mpf(10)**(-42): break
    z = (a+b)/2
    if abs(f(z))/nrm(z) < mp.mpf(10)**(-40):
        return mp.polyval(tc, z)
    raise RuntimeError(f"no right-branch root d={d} s={sv} resid {f(z)}")

_out = {}
def save_out():
    json.dump(_out, open("ring_stageA2.json", "w"), indent=1, default=str)

_WMIN = mp.mpf("1e-11")
def weight_ok(sv):
    return phi(sv)*phi(env_star(sv)) > _WMIN

def leftF(sv, d):  return -sv - 2 + AL(sv)/d + BL(sv)/d**2
def env_star(sv):  return min(-sv-2, sv)

def dm(d):
    sL = np.linspace(-5.4, -1.0, 1401)
    sR = np.linspace(-1.0, 2.8, 1401)
    tot = mp.mpf(0)
    for sv in sL:
        if not weight_ok(sv): continue
        tot += phi(sv)*phi(env_star(sv))*(rightE(mp.mpf(float(sv)), d) - env_star(sv))*(sL[1]-sL[0])
    for sv in sR:
        if not weight_ok(sv): continue
        tot += phi(sv)*phi(env_star(sv))*(leftF(mp.mpf(float(sv)), d) - env_star(sv))*(sR[1]-sR[0])
    return tot

anchors_all = {9: 0.08773584175648119, 11: 0.08755547157788661, 13: 0.08736426802531863,
               15: 0.08719990720393687, 17: 0.08706461572198995, 19: 0.08695451490949664,
               23: 0.086790154156153, 29: 0.08663354848051617, 35: 0.08653775409398155,
               45: 0.08644568837824819, 61: 0.08637768869767828, 91: 0.08633686882830088,
               121: 0.08632885968162476}
run1_diffs = {9: -5.593e-4, 11: -2.729e-4, 13: -1.443e-4, 15: -7.900e-5, 17: -4.31e-5,
              19: -2.21e-5, 23: -1.5e-6, 29: 8.9e-6, 35: 1.15e-5, 45: 1.09e-5,
              61: 8.7e-6, 91: 5.5e-6, 121: 3.8e-6}   # prob units, run-1 stdout
print("="*92); print("(re-audit of the model bias at d = 45, 61, 91, 121)")
diffs = dict(run1_diffs)
for d in (45, 61, 91, 121):
    mm = float(mstar) + float(dm(d))
    df = anchors_all[d] - mm
    print(f"  d={d:3d}: m_mod = {100*mm:.7f}%  meas-model = {100*df:+.6f}pp  (run1: {100*run1_diffs[d]:+.6f}pp)")
    diffs[d] = df

print("="*92); print("(min refine, odd d = 113..131; locks)")
sweep = {}
for d in [113, 117, 119, 121, 123, 125, 127, 129, 131, 141, 151, 161, 171, 181, 201]:
    sweep[d] = float(mstar) + float(dm(d))
    print(f"  m_mod({d:3d}) = {100*sweep[d]:.7f}%", flush=True)
    _out["sweep"] = {str(k): v for k, v in sweep.items()}; save_out()
xs = np.array([113,117,119,121,123,125,127,129,131], float)
ys = np.array([sweep[int(v)] for v in xs])
cf = np.polyfit(xs, ys*100, 2)
dmin = -cf[1]/(2*cf[0]); mmin = np.polyval(cf, dmin)
print(f"  quadratic fit:  d* = {dmin:.2f},  m(d*) = {mmin:.7f}%   (predict d* in (121,131))")

bias45, bias121 = diffs[45], diffs[121]
locks, diflocks = {}, {}
for d in (131, 141, 151, 161, 171, 181, 201):
    lo = sweep[d] + min(bias121 - 8e-6, 0)      # bias shrinking; allow small drift
    hi = sweep[d] + bias45                      # cap by largest smooth bias seen
    locks[d] = (100*lo, 100*hi)
    print(f"  LOCK m({d:3d}) in [{100*lo:.7f}, {100*hi:.7f}] pp-scaled")
diflocks["rebound_strong"] = "m(161) > m(121) - 3e-6  (prob units; m(121) = 0.08632885968162476)"
diflocks["slope_171_161"]  = "m(171) > m(161) - 1.5e-6"
diflocks["slope_201_161"]  = "m(201) > m(161) - 1.5e-6"
diflocks["min_lock"]       = f"d* in (119, 131); m(121) within 6e-7 of sweep min + bias"
print("  DIFF-LOCKS:", *diflocks.values(), sep="\n    ")
_out.update({"mstar": str(mstar), "validation_diffs": diffs,
             "sweep": {str(d): sweep[d] for d in sweep},
             "quadfit_min": {"d": float(dmin), "m_pct": float(mmin)},
             "locks": {str(d): locks[d] for d in locks}, "diflocks": diflocks})
save_out()
print("\nsaved ring_stageA2.json")
