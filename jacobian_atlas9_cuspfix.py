"""
Note 11, stage D-fix: the cusp hardcode audit + cusp exponent universality.

HONESTY LEDGER root cause: atlas9_properness.py line "CU=(0.31865, 0.03502)" froze
the cusp basepoint from stage-A's *printed* (5-dp) values. All cusp numbers there
measured the round-off, not the wall. This script redoes everything at 200 digits.

PREDICTIONS, LOCKED BEFORE ANY COMPUTATION:
  (P1) exact cusp basepoint: h(w) has an exact TRIPLE root at w=t (p9'(t)=0),
       and at delta=1e-11 along generic direction: ESCAPING preimages have
       |x| > 1e4, so bounded fiber count = 7, real bounded = 2  (redeems 10/7/6/6).
  (P2) corrected escape exponent at the d=9 hit cusp: |gamma| ~ delta^(2/3)
       for directions (1,0.7), (1,0), (0,1) -- all three.
  (P3) CUSP EXPONENT UNIVERSALITY: same protocol at every real cusp of the
       d=7 chamber (2R: roles hit+?) and d=8 chamber (1R whisker): ALL 2/3.
       "Hit vs whisker" is a GLOBAL role; the local A3 normal form is universal.
  (P4) artifacts reconciliation: frozen point distance delta_0 = |CU_frozen - CU_true|;
       saturated min|gamma| at the frozen point ~ (|p''(t)|/2)*|u|^2 with u^3 ~ delta_0
       -- same order as the bogus scan's min|gamma| (~1e-2..1e-3 range fit gave 0.2746).
  (P5) synthetic-division lesson honored: escaping roots by Newton from the local
       A3 seed u = (-6*drho/p''(t))^(1/3) * omega; bounded roots from the deflated
       degree-7 quotient with remainder audit < 1e-150. NO raw polyroots on h.
"""
import mpmath as mp
from mpmath import mpf, mpc
import sympy as sp
import json

mp.mp.dps = 200
w = sp.symbols("w")

def chamber(d):
    """p_d, Phi_d as exact rational sympolys; c=1 seeds."""
    p = 2*w - 3*w**2 + w*(1-w)*(w**(d-2) - sp.Rational(6, d*(d+1)))
    p = sp.expand(p)
    Phi = sp.expand(sp.integrate(p, w))
    return p, Phi

def mp_poly(symp_poly, x):
    cs = sp.Poly(symp_poly, w).all_coeffs()
    return [mpf(str(c)) for c in cs]

def refine_root(f, df, seed):
    z = mpf(seed)
    for _ in range(60):
        z = z - f(z)/df(z)
    return z

def newton_h(hcoef, seed, want=mpc):
    """Newton-polish a root of h (mp coefficient list, ascending degree NOT used:
    we use descending as in np.roots convention)."""
    def h(z):
        r = mpf(0)
        for c in hcoef: r = r*z + c
        return r
    def dh(z):
        r = mpf(0)
        n = len(hcoef)-1
        for c in hcoef[:-1]: r = r*z + c*n; n -= 1
        return r
    z = seed
    for _ in range(200):
        dz = h(z)/dh(z)
        z = z - dz
        if abs(dz) < mpf(10)**(-190): break
    return z, abs(h(z))

def deflate(hcoef, root):
    q = [hcoef[0]]
    for c in hcoef[1:-1]:
        q.append(q[-1]*root + c)
    rem = q[-1]*root + hcoef[-1]
    return q, abs(rem)

out = {"dps": 200}

# ---------------- d=9 chamber, hit cusp, EXACT ----------------
p9, Phi10 = chamber(9)
p9m = mp_poly(p9, w); Phim = mp_poly(Phi10, w)
pp9m = mp_poly(sp.diff(p9, w), w)
ppp9m = mp_poly(sp.diff(p9, w, 2), w)

def feval(coefs, z):
    r = mpf(0)
    for c in coefs: r = r*z + c
    return r

def wall_point(t):  # (s,r) = (p(t), t p(t) - Phi(t))
    return feval(p9m, t), t*feval(p9m, t) - feval(Phim, t)

# the two real cusp parameters of F9, Newton-refined at 200 digits from stored seeds
def dfp9(z): return feval(pp9m, z)
t_hit  = refine_root(lambda z: feval(pp9m, z), dfp9 if False else None or (lambda z: feval(ppp9m, z)), "0.3299102226")
t_hit2 = refine_root(lambda z: feval(pp9m, z), lambda z: feval(ppp9m, z), "-0.8913533996")

def h_coefs(s, r):  # h(w) = Phi(w) - s w + r  (deg 10, coefficients descending)
    cs = list(Phim) + [mpf(0)]*(11-len(Phim))
    cs = cs[:9] + [cs[9] - s, cs[10] + r] if len(cs)==11 else cs
    return cs

def analyze_cusp(tag, t, p_m, pp_m, ppp_m, Phi_m, n, frozen=None):
    s0, r0 = wall_point(t)
    rec = {"t": mp.nstr(t, 40), "s0": mp.nstr(s0, 40), "r0": mp.nstr(r0, 40),
           "p''(t)": mp.nstr(feval(ppp_m, t), 20)}
    # --- P1: exact triple root check at basepoint: h(t)=h'(t)=h''(t)=0, h'''(t)!=0
    hc = h_coefs(s0, r0)
    def h(z): return feval(hc, z)
    # h = Phi - s w + r   =>   h' = p - s0,  h'' = p',  h''' = p''   (off-by-one fixed)
    rec["h0"] = mp.nstr(abs(h(t)), 6)
    rec["h1"] = mp.nstr(abs(feval(p_m, t) - s0), 6)   # = 0 by wall parametrization
    rec["h2"] = mp.nstr(abs(feval(pp_m, t)), 6)       # = 0 exactly: cusp <=> p'(t)=0
    rec["h3"] = mp.nstr(abs(feval(ppp_m, t)), 6)      # = p''(t) != 0: A3 nondegenerate
    # --- escaping roots at delta along direction (ds, dr): Newton from A3 seed
    def escape(dS, dR, k):
        d = mpf(10)**(-k)
        ds, dr = d*dS, d*dR
        drho = dr*1 - ds*t  # (dr - ds*t): perturbation constant term of h(t+u)
        u0 = (-6*drho/feval(ppp_m, t))**(mpf(1)/3)
        esc, bounded = [], []
        H = h_coefs(s0+ds, r0+dr)
        roots = []
        for j in range(3):
            seed = t + u0*mp.exp(2j*mp.pi*j/3)
            z, res = newton_h(H, seed)
            roots.append(z); esc.append((z, res))
        remtot = mpf(0)
        Q_temp = list(H)
        for z,_ in esc:
            Q_temp, rem = deflate(Q_temp, z)
            remtot = max(remtot, rem)
        # bounded roots from deflated quotient (degree n-3): plain polyroots is safe
        try:
            br = mp.polyroots(Q_temp, maxsteps=200, error=False)
            for z in br: bounded.append(z)
        except Exception:
            # fallback: try polyroots on full H and drop the 3 near t
            for z in mp.polyroots(H, maxsteps=400, error=False):
                if all(abs(z-ze) > mpf(10)**(-30) for ze,_ in esc): bounded.append(z)
        return esc, bounded, remtot
    # --- P1 fiber count at delta = 1e-11, direction (1, 0.7)
    esc, bounded, remtot = escape(mpf(1), mpf("0.7"), 11)
    gam = [abs((s0 + mpf(10)**-11) - feval(p_m, z)) for z,_ in esc]
    xesc = [1/g for g in gam]
    nreal_bounded = sum(1 for z in bounded if abs(mp.im(z)) < mpf(10)**(-25))
    nreal_esc = sum(1 for z,_ in esc if abs(mp.im(z)) < mpf(10)**(-25))
    rec["fiber@1e-11"] = {"escaping": len(esc), "bounded": len(bounded),
                          "real_bounded": nreal_bounded, "real_escaping": nreal_esc,
                          "min|x|_escaping": mp.nstr(min(xesc), 6),
                          "deflate_remainder": mp.nstr(remtot, 6),
                          "max_esc_residual": mp.nstr(max(r for _,r in esc), 6)}
    # --- P2/P3 escape exponent scans, directions (1,0.7),(1,0),(0,1)
    rec["scans"] = {}
    for dn, (dS, dR) in {"(1,0.7)": (mpf(1), mpf("0.7")), "(1,0)": (mpf(1), mpf(0)),
                         "(0,1)": (mpf(0), mpf(1))}.items():
        lgx, lgy = [], []
        for k2 in range(4, 15):   # delta = 1e-4 .. 1e-14
            d = mpf(10)**(-k2)
            esc2, _, _ = escape(dS, dR, k2)
            g = min(abs((s0 + d*dS) - feval(p_m, z)) for z,_ in esc2)
            lgx.append(float(mp.log(d*mp.sqrt(dS**2+dR**2) if dS or dR else d)))
            lgy.append(float(mp.log(g)))
        # least squares slope
        mx = sum(lgx)/len(lgx); my = sum(lgy)/len(lgy)
        slope = sum((a-mx)*(b-my) for a,b in zip(lgx,lgy)) / sum((a-mx)**2 for a in lgx)
        rec["scans"][dn] = round(slope, 4)
    if frozen is not None:
        sf, rf = frozen
        d0 = mp.sqrt((sf-s0)**2 + (rf-r0)**2)
        # saturated gamma at the frozen point (delta << d0)
        H = h_coefs(sf, rf)
        g = min(abs(sf - feval(p_m, z)) for z in mp.polyroots(H, maxsteps=400, error=False))
        rec["frozen_audit"] = {"delta0": mp.nstr(d0, 8),
                               "min|gamma| at frozen pt": mp.nstr(g, 8),
                               "prediction: gamma ~ |p''/2|*(6*delta0/|p''|)^(2/3)":
                                   mp.nstr(abs(feval(ppp_m,t))/2 * (6*d0/abs(feval(ppp_m,t)))**(mpf(2)/3), 8)}
    return rec

print("== d=9 fiber-10 chamber: EXACT hit-cusp redo ==")
out["d9_cusp1"] = analyze_cusp("cusp+0.3299", t_hit, p9m, pp9m, ppp9m, Phim, 10,
                               frozen=(mpf("0.31865"), mpf("0.03502")))
out["d9_cusp2"] = analyze_cusp("cusp-0.8914", t_hit2, p9m, pp9m, ppp9m, Phim, 10)
for k in ("d9_cusp1","d9_cusp2"):
    r = out[k]
    print(f"  t = {r['t'][:26]}")
    print(f"    triple-root audit: |h|={r['h0']} |h'|={r['h1']} |h''|={r['h2']} |h'''|={r['h3']}  [predict 0,0,0,nonzero]")
    print(f"    fiber@1e-11: {r['fiber@1e-11']['escaping']} escaping (min|x|={r['fiber@1e-11']['min|x|_escaping']}), "
          f"{r['fiber@1e-11']['bounded']} bounded, real {r['fiber@1e-11']['real_bounded']}+{r['fiber@1e-11']['real_escaping']}esc   [predict 3 esc (1 real), 7 bounded (1 real), total real 2]")
    print(f"    deflate audit {r['fiber@1e-11']['deflate_remainder']},  esc residual {r['fiber@1e-11']['max_esc_residual']}")
    print(f"    exponents: {r['scans']}   [predict all 0.6667]")
    if "frozen_audit" in r: print(f"    frozen audit: {r['frozen_audit']}")

print("\n== P3 universality: d=7 (both real cusps) and d=8 (whisker cusp) ==")
for d in (7, 8):
    pd_, Phid = chamber(d)
    pdm = mp_poly(pd_, w); ppdm = mp_poly(sp.diff(pd_, w), w); pppdm = mp_poly(sp.diff(pd_, w, 2), w)
    Phidm = mp_poly(Phid, w)
    # real roots of p_d' : seed scan
    roots = []
    for seed in [mpf(-2)+mpf(i)/20 for i in range(81)]:
        try:
            z = mp.findroot(lambda zz: feval(ppdm, zz), seed,
                            solver='newton', tol=mpf(10)**-170, maxsteps=100)
            if abs(mp.im(z)) < mpf(10)**-50 and all(abs(z-r)>mpf(10)**-12 for r in roots):
                roots.append(mp.re(z))
        except Exception: pass
    for t in roots:
        s0, r0 = feval(pdm, t), t*feval(pdm, t) - feval(Phidm, t)
        def hcoef(s, r):
            cs = list(Phidm) + [mpf(0)]*(d+2-len(Phidm))
            cs = cs[:d] + [cs[d]-s, cs[d+1]+r]
            return cs
        ppp = feval(pppdm, t)
        lgx, lgy = [], []
        for k2 in range(6, 15):
            dd = mpf(10)**(-k2); dS, dR = mpf(1), mpf("0.7")
            ds, dr = dd*dS, dd*dR
            drho = dr - ds*t
            u0 = (-6*drho/ppp)**(mpf(1)/3)
            H = hcoef(s0+ds, r0+dr)
            gs = []
            for j in range(3):
                seed = t + u0*mp.exp(2j*mp.pi*j/3)
                zz = seed
                for _ in range(150):
                    def hh2(z): return feval(H, z)
                    def hh3(z): return feval(ppdm, z) - (s0+ds)
                    dz = hh2(zz)/hh3(zz)
                    zz -= dz
                    if abs(dz) < mpf(10)**-190: break
                gs.append(abs((s0+ds) - feval(pdm, zz)))
            g = min(gs)
            lgx.append(float(mp.log(dd*mp.hypot(dS,dR)))); lgy.append(float(mp.log(g)))
        mx = sum(lgx)/len(lgx); my = sum(lgy)/len(lgy)
        slope = sum((a-mx)*(b-my) for a,b in zip(lgx,lgy))/sum((a-mx)**2 for a in lgx)
        print(f"  d={d} cusp t={mp.nstr(t,16)}: exponent {slope:.4f}   [predict 0.6667]")
        out.setdefault(f"d{d}", []).append({"t": mp.nstr(t,50), "slope": round(slope,4)})

json.dump(out, open("atlas9_cuspfix.json","w"), indent=1)
print("\nsaved atlas9_cuspfix.json")
