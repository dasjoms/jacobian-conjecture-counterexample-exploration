"""
NOTE 13, stage B: EXACT envelope masses at d = 131..201 -- the rebound verdict.

Locks (ring_stageA2.json, made BEFORE these runs):
  windows:  m(131) in [8.6321043, 8.6336200]  ... m(201) in [8.6330844, 8.6346001]
  diffs:    m(161) > m(121) - 3e-6   [note-12's m(161) > m(121), sharpened]
            m(171) > m(161) - 1.5e-6 ;  m(201) > m(161) - 1.5e-6
Machinery: env_at(s) = min over real roots t of p_d(t)=s of tau_d(t), each root
np.roots-guess + 80-digit Newton polish, with resid/nrm acceptance audit and a
dropped-root counter (dropped roots can only bias env UP; report count == 0
goal).  Trapezoid on S in [-8, 4], 2201 points, sigma = 1.5.
"""
import json, sys, math
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 80
sig = 1.5
w = sp.symbols("w")

def poly_data(d):
    p = sp.expand(2*w - 3*w**2 + w*(1-w)*(w**(sp.Integer(d)-2) - sp.Rational(6, d*(d+1))))
    Phi = sp.expand(sp.integrate(p, w))
    pcn = np.array([float(c) for c in sp.Poly(p, w).all_coeffs()])
    tcn = np.array([float(c) for c in sp.Poly(sp.expand(w*p - Phi), w).all_coeffs()])
    coef = [mp.mpf(str(c)) for c in sp.Poly(p, w).all_coeffs()]
    return pcn, tcn, coef

DROPS = {"roots": 0, "kept": 0}
def env_at(sv, pcn, tcn, coef):
    pc = pcn.copy(); pc[-1] -= sv
    cand = np.roots(pc)
    coef_s = list(coef); coef_s[-1] = coef_s[-1] - mp.mpf(float(sv))
    f  = lambda z: mp.polyval(coef_s, z)
    df = lambda z: mp.polyval([c*(len(coef_s)-1-i) for i, c in enumerate(coef_s[:-1])], z)
    nrm = lambda z: sum(abs(c)*abs(z)**k for k, c in enumerate(reversed(coef_s)))
    taus = []
    for r0v in cand:
        if abs(r0v.imag) > 1e-2:
            continue
        z = mp.mpc(mp.mpf(float(r0v.real)), mp.mpf(float(r0v.imag)))
        try:
            for _ in range(40):
                dz = f(z)/df(z); z -= dz
                if abs(dz) < mp.mpf(10)**(-60): break
        except (ZeroDivisionError, ValueError, OverflowError):
            DROPS["roots"] += 1; continue
        ok = (abs(mp.im(z)) < mp.mpf(10)**(-9) and abs(mp.re(z)) < 6
              and abs(f(z))/nrm(z) < mp.mpf(10)**(-45))
        if ok:
            DROPS["kept"] += 1
            taus.append(float(np.polyval(tcn, float(mp.re(z)))))
        elif abs(mp.im(z)) < mp.mpf(10)**(-9):
            DROPS["roots"] += 1
    return min(taus) if taus else np.inf

def envelope_mass(d, ngrid=2201):
    pcn, tcn, coef = poly_data(d)
    S = np.linspace(-8, 4, ngrid)
    Ren = np.array([env_at(sv, pcn, tcn, coef) for sv in S])
    from math import erf as merf
    cdf = np.array([0.5*(1+merf(rv/(sig*math.sqrt(2)))) if np.isfinite(rv) else 1.0 for rv in Ren])
    pdf = np.exp(-S**2/(2*sig**2))/(sig*math.sqrt(2*math.pi))
    return float(np.trapezoid(pdf*cdf, S)), Ren, S

if __name__ == "__main__":
    ds = [int(v) for v in sys.argv[1:]]
    j = json.load(open("ring_stageA2.json"))
    locks = {int(k): v for k, v in j["locks"].items()}
    m121 = 0.08632885968162476
    out = json.load(open("ring_stageB.json")) if __import__("os").path.exists("ring_stageB.json") else {}
    for d in ds:
        DROPS["roots"] = 0; DROPS["kept"] = 0
        m, Ren, S = envelope_mass(d)
        lo, hi = locks[d]
        inwin = float(lo) <= 100*m <= float(hi)
        out[str(d)] = {"mass": m, "mass_pct": 100*m, "lock_window": [lo, hi],
                       "in_window": inwin, "drops": DROPS["roots"], "kept": DROPS["kept"]}
        json.dump(out, open("ring_stageB.json", "w"), indent=1)
        print(f"d={d:3d}:  m = {100*m:.7f}%   window [{lo:.6f}, {hi:.6f}]   INSIDE: {inwin}"
              f"   (drops {DROPS['roots']}/{DROPS['kept']})", flush=True)
    # diff-verdicts when the needed partners exist
    def g(k): return out.get(str(k), {}).get("mass")
    if g(161) and g(171):
        print(f"DIFF m(171)-m(161) = {1e6*(g(171)-g(161)):+.3f}e-6   lock > -1.5e-6: {(g(171)-g(161)) > -1.5e-6}")
    if g(161) and g(201):
        print(f"DIFF m(201)-m(161) = {1e6*(g(201)-g(161)):+.3f}e-6   lock > -1.5e-6: {(g(201)-g(161)) > -1.5e-6}")
    if g(161):
        print(f"DIFF m(161)-m(121) = {1e6*(g(161)-m121):+.3f}e-6   "
              f"note-12 lock m(161) > m(121): {g(161) > m121};  strong lock > -3e-6: {(g(161)-m121) > -3e-6}")
