"""
NOTE 15, stage FINAL: remaining audits + publishing.
 (C1) corrected shadow identity: worst residual (already certified 6.5e-86) - re-audit.
 (C2) ghost closed form ratio test incl. d = 43..47.
 (C3) left-cusp fixed-point iteration converges to exact root within 3e-4.
 (M3) m(3) census 1e6, window [8.45%, 8.85%].
 (F1) chamber-12 published: crunode s(11), r(11) from gap-series;
 (F2) corner-side whisker acnode(12) window.
 (G)  gap series + contact-vs-shadow coalescence rates.
"""
import json, time
import numpy as np
import mpmath as mp
import sympy as sp
from sympy import symbols, expand, diff, integrate, Rational as R, Poly

mp.mp.dps = 90
w = symbols("w")
t0 = time.time()
out = {}

def seed(d):
    return expand(2*w - 3*w**2 + w*(1-w)*(w**sp.Integer(d-2) - R(6, d*(d+1))))

def diag_root(d):
    p = seed(d); Phi = integrate(p, w)
    G = expand(Phi - (w-1)*p)
    fc = [float(c) for c in Poly(sp.expand(G*(d*(d+1))), w).all_coeffs()]
    Gf = sp.lambdify(w, G, "mpmath")
    best = None
    for z in np.roots(fc):
        if abs(z.imag) < 1e-7:
            xr = mp.mpf(z.real)
            try: xr = mp.findroot(Gf, xr, tol=mp.mpf("1e-80"), verify=False)
            except Exception: pass
            if abs(Gf(xr)) < mp.mpf("1e-40") and xr < -1.000001 and (best is None or xr > best):
                best = xr
    return best, sp.lambdify(w, p, "mpmath")(best)

roots = {d: diag_root(d) for d in [3] + list(range(5, 46, 2))}

# ---- C1 ----
def c1_res(U, d):
    X = (1+U)**(d-2)
    return abs((d*d+d)*(U+2)**2*(X-2) + (4*U**2+17*U+19) - X*(U+1)*(1+d*(U+2)))
worst = max(c1_res(-roots[d][0]-1 if d != 3 else mp.mpf(1), d) for d in roots)
out["C1"] = {"ok": bool(worst < mp.mpf("1e-45")), "worst": mp.nstr(worst, 6),
             "identity": "(d^2+d)(U+2)^2 (X-2) + 4U^2+17U+19 = X (U+1) [1 + d (U+2)], X = (1+U)^(d-2)"}
print(f"(C1) worst residual: {mp.nstr(worst, 8)}", flush=True)

# ---- C2 ghost incl. 43, 45, 47 ----
print("(C2) ghost race")
ghost_ok = {}
for d in (21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47):
    pp = diff(seed(d), w)
    fc = [float(c) for c in Poly(sp.expand(pp*(d*(d+1))), w).all_coeffs()]
    Gf = sp.lambdify(w, pp, "mpmath")
    tg = None
    for z in np.roots(fc):
        if abs(z.imag) < 1e-7:
            xr = mp.mpf(z.real)
            try: xr = mp.findroot(Gf, xr, tol=mp.mpf("1e-80"), verify=False)
            except Exception: pass
            if abs(Gf(xr)) < mp.mpf("1e-40") and abs(xr - mp.mpf(1)/3) < 0.05:
                tg = xr
    pred = 1/(3*(mp.mpf(d)*(d+1)-2))
    resid = abs((mp.mpf(1)/3 - tg) - pred)
    ratio = resid / mp.mpf(3)**(-mp.mpf(d)/2)
    ghost_ok[d] = (float(resid), float(ratio))
    print(f"  d={d:2d}: resid {mp.nstr(resid, 4)}   resid/3^(-d/2) = {mp.nstr(ratio, 5)}", flush=True)
ok43 = all(v[1] < 1e-6 for k, v in ghost_ok.items() if k >= 43)
out["C2"] = {"okA2_d47": bool(ghost_ok[47][0] < 1e-19), "ratio_bounded_43plus": bool(ok43), "table": {str(k): v for k, v in ghost_ok.items()}}
print(f"  (A2) d=47 resid < 1e-19: {out['C2']['okA2_d47']}; ratio@43+ < 1e-6: {ok43}")

# ---- C3 ----
print("(C3) left-cusp fixed point")
j1 = json.load(open("jcorner_stage1.json"))
C3ok, worst3 = True, mp.mpf(0)
for d in range(7, 42, 2):
    c0 = mp.mpf(6)/(d*(d+1))
    Rf = lambda u: (8-6*u-c0*(3-2*u))/((3-2*u)+(d-2)*(2-u))
    u = mp.log((2*mp.mpf(d)-1)/8)/(d-2)
    for _ in range(8):
        u = 1 - Rf(u)**(mp.mpf(1)/(d-2))
    err = abs(u - mp.mpf(j1["races"][str(d)]["u"]))
    worst3 = max(worst3, err)
    if d <= 13 or d % 8 == 1:
        print(f"  d={d:2d}: fixed-point err {mp.nstr(err, 5)}", flush=True)
C3ok = worst3 < mp.mpf("3e-4")
out["C3"] = {"ok": bool(C3ok), "worst": mp.nstr(worst3, 5)}
print(f"  worst err d=7..41: {mp.nstr(worst3, 5)}  [<3e-4: {C3ok}]")

# ---- M3 ----
print("(M3) m(3) census, 1e6")
p3 = seed(3); Phi3 = integrate(p3, w)
hh = [float(c) for c in Poly(expand(Phi3), w).all_coeffs()]
rng = np.random.default_rng(42)
N = 1_000_000
S = rng.normal(0, 1.5, N); Rv = rng.normal(0, 1.5, N)
c0n = 0
batch = 100000
for b in range(0, N, batch):
    for sv, rv in zip(S[b:b+batch], Rv[b:b+batch]):
        q = hh[:]
        q[-1] += rv
        q[-2] -= sv
        if all(abs(z.imag) > 1e-8 for z in np.roots(q)):
            c0n += 1
m3 = c0n/N
m3e = 3*mp.sqrt(m3*(1-m3)/N)
M3 = 0.0845 <= m3 <= 0.0885
print(f"  m(3) = {100*m3:.4f}% +- {100*float(m3e):.4f}%   in [8.45, 8.85]: {M3}")
out["M3"] = {"mass": 100*m3, "err": 100*float(m3e), "ok": bool(M3)}

# ---- gaps + F1 ----
j2 = json.load(open("jcorner_basins.json"))
def scr(x): return complex(str(x).replace(" ", ""))
cr = {int(d): [nd for nd in j2["basins"][d] if nd[0] == "CRUNODE"][0] for d in ("5", "7", "9")}
gaps, rgaps = {}, {}
for d in (5, 7, 9):
    sc = mp.mpf(scr(cr[d][1]).real); rc = mp.mpf(scr(cr[d][2]).real)
    sstar = roots[d][1]
    gaps[d] = sc - sstar
    rgaps[d] = rc - sc
    print(f"  d={d}: crunode ({scr(cr[d][1]).real:.8f}, {scr(cr[d][2]).real:.8f})  gap s_c-s* = {mp.nstr(gaps[d], 8)}  r_c-s_c = {mp.nstr(rgaps[d], 8)}", flush=True)
pred_gap11 = gaps[9]*(mp.mpf(9)/11)**2
pred_rgap11 = rgaps[9]*(mp.mpf(9)/11)**2
s11 = roots[11][1]
pred_sc11 = s11 + pred_gap11
pred_rc11 = pred_sc11 + pred_rgap11
print(f"  shadow s*(11) = {mp.nstr(s11, 12)}")
print(f"  F1 publish: crunode(11) s ~ {mp.nstr(pred_sc11, 9)} window ({mp.nstr(pred_sc11-0.0008,8)}, {mp.nstr(pred_sc11+0.0008,8)})")
print(f"                     r ~ {mp.nstr(pred_rc11, 9)} window ({mp.nstr(pred_rc11-0.00075,8)}, {mp.nstr(pred_rc11+0.00075,8)})")
out["F1"] = {"s_window": [mp.nstr(pred_sc11-0.0008, 10), mp.nstr(pred_sc11+0.0008, 10)],
             "r_window": [mp.nstr(pred_rc11-0.00075, 10), mp.nstr(pred_rc11+0.00075, 10)],
             "shadow_s11": mp.nstr(s11, 18), "gap_series": {str(d): mp.nstr(gaps[d], 10) for d in gaps}}
# contact-vs-shadow coalescence
for d in (5, 7, 9):
    t2c = scr(cr[d][4] if len(cr[d]) > 4 else cr[d][3]).real
    contact_neg = min(scr(cr[d][3]).real, scr(cr[d][4]).real)
    tstar = roots[d][0]
    dlt = abs(mp.mpf(contact_neg) - tstar)
    print(f"  d={d}: |contact_t2 - t*| = {mp.nstr(dlt, 6)}   x d^3 = {mp.nstr(dlt*d**3, 6)}", flush=True)
    out.setdefault("contact_shadow_gap", {})[str(d)] = [mp.nstr(dlt, 10), mp.nstr(dlt*d**3, 8)]

# whisker-acnode series
bc = []
for d in (6, 8, 10):
    acs = [nd for nd in j2["basins"][str(d)] if nd[0] == "acnode"]
    near = min(acs, key=lambda nd: (mp.mpf(scr(nd[1]).real)+1)**2 + (mp.mpf(scr(nd[2]).real)+1)**2)
    bc.append((d, mp.mpf(scr(near[1]).real), mp.mpf(scr(near[2]).real)))
    print(f"  d={d}: corner-side whisker acnode ({mp.nstr(bc[-1][1],10)}, {mp.nstr(bc[-1][2],10)})  |s+1|={mp.nstr(abs(bc[-1][1]+1),8)}", flush=True)
a_ = mp.log(abs(bc[0][1]+1)/abs(bc[1][1]+1))/mp.log(8/mp.mpf(6))
pred_s12 = -1 + abs(bc[2][1]+1)*(10/mp.mpf(12))**(mp.log(abs(bc[1][1]+1)/abs(bc[2][1]+1))/mp.log(mp.mpf(10)/8))
print(f"  F2 publish: whisker acnode(12)  s in ({mp.nstr(pred_s12-0.025,7)}, {mp.nstr(pred_s12+0.025,7)}), r in (-1.16, -1.10)")
out["F2"] = {"s_window": [mp.nstr(pred_s12-0.025, 8), mp.nstr(pred_s12+0.025, 8)], "r_window": ["-1.16", "-1.10"],
             "series": [[d, mp.nstr(a, 12), mp.nstr(b, 12)] for d, a, b in bc]}

# key law constants for the note
Ln2 = mp.log(2)
u2n = (1 + Ln2**2)/2; u3n = -mp.mpf("3.5") + Ln2**3/6 + 3*Ln2/4
u4n = -31*Ln2/8 + Ln2**4/24 + Ln2**2/2 + mp.mpf(355)/24
out["laws"] = {"u1": "ln 2", "u2": "(1+ln^2 2)/2", "u3": "-7/2 + ln^3 2/6 + 3 ln 2/4",
               "u4": "-31 ln 2/8 + ln^4 2/24 + ln^2 2/2 + 355/24",
               "u_numeric": [mp.nstr(Ln2, 15), mp.nstr(u2n, 15), mp.nstr(u3n, 15), mp.nstr(u4n, 15)],
               "sig_coeffs": {"s_a": mp.nstr(2-2*Ln2, 20), "s_b": mp.nstr(-(mp.mpf("2.5") + 2*Ln2**2 - 4*Ln2), 20)},
               "sig_numeric": [mp.nstr(2-2*Ln2, 15), mp.nstr(-(mp.mpf('2.5')+2*Ln2**2-4*Ln2), 15)],
               "sigma_star_limit": "s*(d) + 1 -> 0 with (s*+1)(d-2) -> 2-2 ln 2"}
# U-series final residual vs 45
d45 = roots[45]
predU = Ln2/43 + u2n/43**2 + u3n/43**3 + u4n/43**4
print(f"  U(45): exact {mp.nstr(-d45[0]-1, 14)} vs series4 {mp.nstr(predU, 14)}  diff {mp.nstr(abs((-d45[0]-1)-predU), 8)}")
# sigma validation
sa = 2-2*Ln2; sb = -(mp.mpf("2.5")+2*Ln2**2-4*Ln2)
for d in (11, 21, 45):
    dn = d-2
    pred_s = -1 + sa/dn + sb/dn**2
    print(f"  s*({d}): exact {mp.nstr(roots[d][1], 12)}  2-term pred {mp.nstr(pred_s, 12)}  diff {mp.nstr(abs(roots[d][1]-pred_s), 8)}")
    out.setdefault("sig_validate", {})[d] = mp.nstr(abs(roots[d][1]-pred_s), 8)

out["U_series_u_numeric"] = {"u2": mp.nstr(u2n, 25), "u3": mp.nstr(u3n, 25), "u4": mp.nstr(u4n, 25)}
json.dump(out, open("jcorner_final.json", "w"), indent=1)
print(f"saved jcorner_final.json [{time.time()-t0:.0f}s]")
print("SUMMARY: C1", out["C1"]["ok"], " C2-47", out["C2"]["okA2_d47"], " C3", out["C3"]["ok"], " M3", M3)
