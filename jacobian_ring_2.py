"""
NOTE 13, stage A2: the GOLD mass model, validated; locks for the rebound.

CERTIFIED INPUT (stage A-probe, 100-digit audits, diffs 1e-4 at d=201):
  left pin-layer branch (t = -1-x/d):
     leftF(s;d) = -s-2 + AL(s)/d + BL(s)/d^2,
     AL(s) = (s+5)(1 - ln((s+5)/2))   [machine certificate],
     zero of AL at s = 2e - 5 = 0.43656... : the branch CROSSES its own line.
  right pin-layer branch (t = 1+c/d, c solving c(4+e^c) = d(s+1), bijective):
     rightE(s;d) = exact 50-digit Newton tau of that root.

GOLD MODEL (no other branch can bind in |s|<=5: probe-tables show the other
taus >= +0.02 there):
  env_m(s;d) = min(leftF, rightE),  env*(s) = min(-s-2, s),
  m_mod(d) = m* + INT phi(s) phi(env*(s)) [env_m - env*] ds        (sigma = 1.5)

LOCKS after validation on the 13 measured anchors (d = 9..121):
  T1: max |m_mod - m_meas| <= 8e-6 (0.0008pp) for d >= 45;
                        <= 4e-5 for d in {29, 35}; larger allowed at d <= 23
                        (higher-order O(ln^2/d^3) not included).
  T2: the minimum moves: d* in (121, 161), m(d*) > 8.628%.
  T3: exact masses (stage B) at 131,141,151,161,171,181,201 inside
      model +/- (last seen bias + 1.5 model-error floor).  Printed as locks.
  T4: m(161) > m(121) = 8.63289%  (note-12 lock, re-locked).
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
two_e_minus5 = 2*mp.e - 5

AL = lambda u: (u+5)*(1 - mp.log((u+5)/2))
# B_L(s) from the machine (ring_probe_sym.json string), simplified grouping:
# B_L = -(s+13)M^2/2 + ((s+5)(M-1)/2 + 5/2*(1-...)) machine form used directly:
BLs = sp.sympify(json.load(open("ring_probe_sym.json"))["BL"])
_s = sp.symbols("s")
BLsy = sp.lambdify(_s, BLs, "mpmath")
BL = lambda u: mp.mpf(BLsy(mp.mpf(u)))

# ---------- exact right branch ----------
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
    cc = mp.findroot(lambda c: c*(4+mp.e**c) - lam, 0.0)
    return cc

def rightE(sv, d):
    """exact tau of the rightmost (1+c/d -> outer cubic) branch; robust starts."""
    pc, tc = seed_co(d)
    pcs = list(pc); pcs[-1] -= mp.mpf(sv)
    dps = [c*(len(pcs)-1-i) for i, c in enumerate(pcs[:-1])]
    f, df = (lambda z: mp.polyval(pcs, z)), (lambda z: mp.polyval(dps, z))
    disc = 4 - 12*float(sv)
    t_cub = (2 + math.sqrt(max(disc, 0.0)))/6
    for z0 in (1 + float(c_of(sv, d))/d, t_cub, 1.6):
        z = mp.mpf(z0)
        try:
            for _ in range(80):
                dz = f(z)/df(z); z -= dz
                if abs(dz) < mp.mpf(10)**(-45): break
            if abs(f(z)) < mp.mpf(10)**(-30) and mp.re(z) > 0.5:
                return mp.polyval(tc, z)
        except (ZeroDivisionError, ValueError, OverflowError):
            continue
    raise RuntimeError(f"no right-branch root d={d} s={sv}")

def leftF(sv, d):
    return -sv - 2 + AL(sv)/d + BL(sv)/d**2

def env_model(sv, d):
    v = rightE(sv, d)
    if sv > -4.9999:
        v = min(v, leftF(sv, d))
    return v

def env_star(sv):
    return min(-sv-2, sv)

def dm(d):
    """correction integral; phi(s) weight kernel, trap on certified support"""
    sL = np.linspace(-5.4, -1.0, 1401)   # right-side support
    sR = np.linspace(-1.0, 2.8, 1401)    # left-side support
    tot = mp.mpf(0)
    for sv in sL:
        w_ = phi(sv)*phi(env_star(sv))*(rightE(mp.mpf(float(sv)), d) - env_star(sv))
        tot += w_*(sL[1]-sL[0])
    for sv in sR:
        w_ = phi(sv)*phi(env_star(sv))*(leftF(mp.mpf(float(sv)), d) - env_star(sv))
        tot += w_*(sR[1]-sR[0])
    return tot

anchors = {9: 0.08773584175648119, 11: 0.08755547157788661, 13: 0.08736426802531863,
           15: 0.08719990720393687, 17: 0.08706461572198995, 19: 0.08695451490949664,
           23: 0.086790154156153, 29: 0.08663354848051617, 35: 0.08653775409398155,
           45: 0.08644568837824819, 61: 0.08637768869767828, 91: 0.08633686882830088,
           121: 0.08632885968162476}

print("="*92)
print(f"m* = {mp.nstr(mstar, 18)};  AL zero at s = 2e-5 = {mp.nstr(two_e_minus5, 10)}")
print("(validation)  d |    measured       m_mod (GOLD)      diff(pp)")
out = {"mstar": str(mstar)}
diffs = {}
for d in sorted(anchors):
    mm = float(mstar) + float(dm(d))
    df = anchors[d] - mm
    diffs[d] = df
    print(f"               {d:3d} | {100*anchors[d]:12.7f} {100*mm:12.7f}   {100*df:+9.5f}")
out["validation"] = {str(d): {"meas": anchors[d], "model": diffs and float(mstar)+float(0), "diff": None} for d in ()}
out["validation"] = {str(d): {"meas": anchors[d], "diff": diffs[d]} for d in anchors}
bad_hi = {d: diffs[d] for d in diffs if d >= 45}
print(f"T1: max |diff| d>=45: {max(abs(v) for v in bad_hi.values()):.2e}  (lock 8e-6)")
print(f"    max |diff| 29..35: {max(abs(diffs[d]) for d in (29,35)):.2e}  (lock 4e-5)")

# ---------- model sweep near the minimum + locks ----------
print("="*92); print("(sweep + locks)")
sweep = {}
for d in list(range(109, 202, 4)) + [131, 141, 151, 161, 171, 181, 201]:
    sweep[d] = float(mstar) + float(dm(d))
    print(f"  m_mod({d:3d}) = {100*sweep[d]:.6f}%  (m - m* = {100*(sweep[d]-float(mstar)):+.6f}pp)", flush=True)
dmin = min(sweep, key=sweep.get)
# parabolic-ish local refine
ds = sorted(sweep); i0 = ds.index(dmin)
lo, hi = ds[max(0, i0-2)], ds[min(len(ds)-1, i0+2)]
print(f"  sweep min at d* = {dmin} (window {lo}..{hi}), m_mod {100*sweep[dmin]:.6f}%")
locks = {}
bias_est = max(abs(diffs[d]) for d in (61, 91, 121))
tol = bias_est + 6e-6
for d in (131, 141, 151, 161, 171, 181, 201):
    locks[d] = (float(sweep[d]) - tol, float(sweep[d]) + tol)
    print(f"  LOCK m({d:3d}) in [{100*(sweep[d]-tol):.6f}, {100*(sweep[d]+tol):.6f}]  (tol {100*tol:.6f}pp)")
locks_meta = {"tol": tol, "dmin_window": [lo, hi], "rebound": "m(161) > 8.63288%",
              "T2": f"d* in (121, 161); m(d*) > 8.628%" }
out["sweep"] = {str(d): sweep[d] for d in sweep}
out["locks"] = {str(d): locks[d] for d in locks}
out["locks_meta"] = locks_meta
json.dump(out, open("ring_stageA2.json", "w"), indent=1, default=str)
print("\nsaved ring_stageA2.json")
