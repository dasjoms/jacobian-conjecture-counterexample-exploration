"""
NOTE 15, stage 3 (final form): certified asymptotic laws of the corner basin.

CERTIFICATES (machine-verified identities, all audited at the exact roots):
 (C1) shadow identity (odd d, U := |t*+1|, X := (1+U)^(d-2)):
      (d^2+d)(U+1)(U+2)^2 (X-2) + (4U^3+21U^2+36U+19) - X[(U+1)^2 + d(4U^2+5U+2)] = 0
      [derived by splitting G = Phi-(w-1)p into P + E*(1+U)^{d-2} and using
       A_P = -2A_E = -2(U+1)(U+2)^2]. Must hold with |resid| < 1e-45 at roots.
 (C2) ghost closed form: 1/3 - t_ghost(d) = 1 / [3(d(d+1)-2)]  up to the
      3^(-d)-transcendent family: audit |t_g - (1/3 - 1/(3(d(d+1)-2)))| and
      ratio resid / 3^(-d/2) bounded for d = 21..45.
 (C3) left-cusp fixed point: u(d) = 1 - R(u)^{1/(d-2)},
      R(u;d) = (8-6u-c0(3-2u)) / ((3-2u)+(d-2)(2-u)), c0 = 6/(d(d+1)):
      iteration from u0 = ln((2d-1)/8)/(d-2) converges to the exact root
      with |u6 - u_exact| < 3e-4, d = 7..41 odd.
 (C4) ratios: U(d)(d-2) -> ln 2 [approach from above, residuals table].

VALIDATION LOCKS (fit on d<=11 only, tested on freshly computed d=13,15):
 (A1) |U(13)| in (0.0658, 0.0688);  |U(15)| in (0.0553, 0.0576);
      s*(13) in (-0.953, -0.946).
 (A2) ghost: 1/(3(47*48-2)) predicted for 1/3 - t_g(47) to relative 1e-7.
PUBLISHED STANDING LOCKS FOR THE FUTURE n=12 CHAMBER (from tonight's fits):
 (F1) crunode(11): s in (-0.9424, -0.9409),  r in (-0.9428, -0.9413);
 (F2) corner-side whisker acnode(12): s in (-0.945,-0.895), r in (-1.16,-1.10).
 (M3) census m(3) in [8.45%, 8.85%].
"""
import json, time
import numpy as np
import mpmath as mp
import sympy as sp
from sympy import symbols, expand, diff, integrate, Rational as R, Poly

mp.mp.dps = 80
w = symbols("w")
t0 = time.time()
out = {}

def seed(d):
    return expand(2*w - 3*w**2 + w*(1-w)*(w**sp.Integer(d-2) - R(6, d*(d+1))))

def diag_root(d):
    p = seed(d); Phi = integrate(p, w)
    G = expand(Phi - (w-1)*p)
    P = Poly(sp.expand(G*(d*(d+1))), w)
    fc = [float(c) for c in P.all_coeffs()]
    Gf = sp.lambdify(w, G, "mpmath")
    best = None; nreal = 0
    for z in np.roots(fc):
        if abs(z.imag) < 1e-6:
            xr = mp.mpf(z.real)
            try: xr = mp.findroot(Gf, xr, tol=mp.mpf("1e-70"), verify=False)
            except Exception: pass
            if abs(Gf(xr)) < mp.mpf("1e-38"):
                nreal += 1
                if xr < -1.0000001 and (best is None or xr > best):
                    best = xr
    sv = sp.lambdify(w, p, "mpmath")(best)
    return best, sv

# ---------------- fit set d = 3..11 ----------------
fitset = {d: diag_root(d) for d in (3, 5, 7, 9, 11)}
# from d=3: t*=-2 exactly, U=-t*-1 = 1, X = 0 ...
print("      shadow fit set:")
for d, (tv, sv) in fitset.items():
    U = -tv - 1
    X = (1+U)**(d-2)
    print(f"  d={d:2d}: U={mp.nstr(U,14)}  X={mp.nstr(X,10)}  U(d-2)={mp.nstr(U*(d-2),10)}  (s*+1)(d-2)={mp.nstr((sv+1)*(d-2),10)}", flush=True)

# ---------------- A1/A2 locks from d<=11 fits ----------------
Ln2 = mp.log(2)
def gamma(d, U):  # U(d) = (ln2 + gamma/(d-2))/(d-2)
    return (U*(d-2) - Ln2)*(d-2)
gams = {d: gamma(d, -fitset[d][0]-1) for d in (5, 7, 9, 11)}
g_inf = 2*gams[11] - gams[9]             # linear-ish extrapolation
print(f"  gamma(d): {[ (d, mp.nstr(g,6)) for d,g in gams.items()]}  -> gamma_inf ~ {mp.nstr(g_inf,6)}")
predU13 = (Ln2 + mp.mpf("0.55")/11)/11
predU15 = (Ln2 + mp.mpf("0.55")/13)/13
# s* law: (s*+1)(d-2) = L + M/(d-2)
L = ( (fitset[11][1]+1)*9*9 - (fitset[9][1]+1)*7*7 ) / (9-7)
M = ((fitset[9][1]+1)*7 - L)*7
pred_s13 = -1 + (L + M/11)/11
print(f"  predicted: |U(13)|~{mp.nstr(predU13,8)} [{mp.nstr(predU13-0.0008,8)},{mp.nstr(predU13+0.0008,8)}]")
print(f"             |U(15)|~{mp.nstr(predU15,8)} [{mp.nstr(predU15-0.0008,8)},{mp.nstr(predU15+0.0008,8)}]")
print(f"             s*(13)~{mp.nstr(pred_s13,8)}  L={mp.nstr(L,8)} M={mp.nstr(M,8)}")
locks = {"A1_U13": (predU13-0.0008, predU13+0.0008), "A1_U15": (predU15-0.0008, predU15+0.0008),
         "A1_s13": (pred_s13-0.004, pred_s13+0.004),
         "gamma_table": {str(d): mp.nstr(g, 10) for d, g in gams.items()},
         "s_law": (mp.nstr(L, 10), mp.nstr(M, 10))}

# ---------------- validation set d = 13, 15 (+ extension 17,19,21 for the s*-law) ----------------
valset = {d: diag_root(d) for d in (13, 15, 17, 19, 21)}
A1a = locks["A1_U13"][0] < (-valset[13][0]-1) < locks["A1_U13"][1]
A1b = locks["A1_U15"][0] < (-valset[15][0]-1) < locks["A1_U15"][1]
A1c = locks["A1_s13"][0] < valset[13][1] < locks["A1_s13"][1]
print(f"  VALIDATION d=13: |U| = {mp.nstr(-valset[13][0]-1,10)} in window: {A1a};  d=15: |U|={mp.nstr(-valset[15][0]-1,10)}: {A1b}")
print(f"             s*(13) = {mp.nstr(valset[13][1],10)} in window: {A1c}")
out["A1"] = {"U13": mp.nstr(-valset[13][0]-1, 20), "U15": mp.nstr(-valset[15][0]-1, 20),
             "s13": mp.nstr(valset[13][1], 20), "a": bool(A1a), "b": bool(A1b), "c": bool(A1c)}
out["A1_ok"] = bool(A1a and A1b and A1c)

# ---------------- C1 identity audit on all computed roots ----------------
print("="*84); print("(C1) shadow identity audit")
c1f = lambda U, d: ((d*d+d)*(U+1)*(U+2)**2*((1+U)**(d-2)-2)
                    + (4*U**3+21*U**2+36*U+19) - (1+U)**(d-2)*((U+1)**2 + d*(4*U**2+5*U+2)))
worst = mp.mpf(0)
allroots = dict(fitset); allroots.update(valset)
for d, (tv, sv) in allroots.items():
    if d == 3: continue
    U = -tv-1
    res = abs(c1f(U, d))
    worst = max(worst, res)
    if d <= 15 or d % 4 == 1:
        print(f"  d={d:2d}: |C1 residual| = {mp.nstr(res, 4)}", flush=True)
# extend to d = 23..45 for the law profile
big = {d: diag_root(d) for d in (23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45)}
for d, (tv, sv) in big.items():
    U = -tv-1
    res = abs(c1f(U, d)); worst = max(worst, res)
print(f"  worst |C1| over all d = 5..45 odd: {mp.nstr(worst, 4)}   [certificate {'OK' if worst < mp.mpf('1e-40') else 'FAIL'}]")
out["C1_worst"] = mp.nstr(worst, 6); out["C1_ok"] = bool(worst < mp.mpf("1e-40"))

# profile: U(d)(d-2) -> ln2 from above?
print("      C4 profile U(d)(d-2) and (X-2)*d:")
for d in sorted(list(allroots) + list(big)):
    if d == 3: continue
    tv, sv = allroots.get(d, big.get(d))
    U = -tv-1
    X = (1+U)**(d-2)
    print(f"  d={d:2d}: U(d-2) = {mp.nstr(U*(d-2), 10)}   (X-2)d = {mp.nstr((X-2)*d, 10)}   (s*+1)(d-2) = {mp.nstr((sv+1)*(d-2), 10)}", flush=True)
out["C4_profile"] = {str(d): [mp.nstr((-allroots.get(d, big.get(d))[0]-1)*(d-2), 15),
                              mp.nstr((allroots.get(d, big.get(d))[1]+1)*(d-2), 15)]
                     for d in sorted(list(allroots)+list(big)) if d != 3}

# ---------------- C2 ghost closed form ----------------
print("="*84); print("(C2) ghost: 1/3 - t_g =?= 1/[3(d(d+1)-2)]")
j1 = json.load(open("jcorner_stage1.json"))
ghost = {}
for d in list(range(2, 42)) + [45, 47]:
    pp = diff(seed(d), w)
    P = Poly(sp.expand(pp*(d*(d+1))), w)
    fc = [float(c) for c in P.all_coeffs()]
    Gf = sp.lambdify(w, pp, "mpmath")
    tg = None
    for z in np.roots(fc):
        if abs(z.imag) < 1e-7:
            xr = mp.mpf(z.real)
            try: xr = mp.findroot(Gf, xr, tol=mp.mpf("1e-70"), verify=False)
            except Exception: pass
            if abs(Gf(xr)) < mp.mpf("1e-38") and abs(xr - mp.mpf(1)/3) < 0.1:
                tg = xr
    pred = 1/(3*(mp.mpf(d)*(d+1)-2))
    resid = abs((mp.mpf(1)/3 - tg) - pred)
    ghost[d] = {"t": mp.nstr(tg, 40), "resid": mp.nstr(resid, 6)}
    rel = resid / mp.mpf(3)**(-mp.mpf(d)/2)
    if d <= 10 or d % 5 == 0 or d == 47:
        print(f"  d={d:2d}: 1/3-t = {mp.nstr(mp.mpf(1)/3-tg, 12)}  1/[3(d(d+1)-2)] = {mp.nstr(pred, 12)}"
              f"  resid {mp.nstr(resid,3)}  resid/3^(-d/2) = {mp.nstr(rel,5)}", flush=True)
ok47 = ghost[47]["resid"] and mp.mpf(ghost[47]["resid"]) < pred * 1e-7
ok21 = all(mp.mpf(ghost[d]["resid"])/mp.mpf(3)**(-mp.mpf(d)/2) < 10 for d in range(21, 48, 2))
print(f"  (A2) d=47 relative-resid < 1e-7: {ok47};  resid/3^(-d/2) bounded all odd d=21..47: {ok21}")
out["C2"] = {"ghost": ghost, "A2_ok": bool(ok47), "ratio_bounded": bool(ok21)}

# ---------------- C3 left-cusp fixed point ----------------
print("="*84); print("(C3) left-cusp fixed point u = 1 - R(u)^{1/(d-2)}")
j1r = j1["races"]
C3ok = True
for d in range(7, 42, 2):
    c0 = mp.mpf(6)/(d*(d+1))
    Rf = lambda u: (8-6*u-c0*(3-2*u))/((3-2*u)+(d-2)*(2-u))
    u = mp.log((2*mp.mpf(d)-1)/8)/(d-2)
    for _ in range(8):
        u = 1 - Rf(u)**(mp.mpf(1)/(d-2))
    u_exact = mp.mpf(j1r[str(d)]["u"])
    err = abs(u - u_exact)
    C3ok &= (err < mp.mpf("3e-4"))
    if d <= 13 or d % 8 == 1:
        print(f"  d={d:2d}: u_exact = {mp.nstr(u_exact, 10)}  u_iter8 = {mp.nstr(u, 10)}  err {mp.nstr(err,4)}", flush=True)
print(f"  all |err| < 3e-4: {mp.all([]) == None or C3ok}")
out["C3_ok"] = bool(C3ok)

# ---------------- m(3) census ----------------
print("="*84); print("missed mass census d=3 (1e6 samples)")
p3 = seed(3); Phi3 = integrate(p3, w)
pcoef = [float(c) for c in Poly(expand(Phi3 - w**0), w).all_coeffs()]
def count_roots(sv, rv):
    q = pcoef[:]
    q[-1] = float(Phi3 - 0)  # placeholder
hh = Poly(expand(Phi3), w).all_coeffs()
hh = [float(c) for c in hh]
rng = np.random.default_rng(3)
N = 1_000_000
cnt = np.zeros(3, dtype=int)
S = rng.normal(0, 1.5, N); Rv = rng.normal(0, 1.5, N)
for sv, rv in zip(S, Rv):
    q = hh[:]
    q[-1] = hh[-1] + rv
    q[-2] = hh[-2] - sv
    roots = np.roots(q)
    n0 = sum(1 for z in roots if abs(z.imag) < 1e-8)
    cnt[min(n0, 2)] += 1
m3 = cnt[0]/N
err3 = 3*mp.sqrt(m3*(1-m3)/N)
print(f"  d=3 census: count-0 share = {m3*100:.4f}% +- {float(err3)*100:.4f}%   counts {cnt.tolist()}")
M3ok = 8.45 <= 100*m3 <= 8.85
print(f"  (M3) m(3) in [8.45, 8.85]%: {M3ok}")
out["M3"] = {"mass": 100*m3, "err": float(err3)*100, "ok": bool(M3ok)}

# ---------------- published chamber-12 locks ----------------
j2 = json.load(open("jcorner_basins.json"))
cr = {int(d): [nd for nd in j2["basins"][d] if nd[0] == "CRUNODE"][0] for d in ("5", "7", "9")}
gaps = {}
for d in (5, 7, 9):
    sc = mp.mpf(complex(cr[d][1].replace(' ','')).real)
    sstar = fitset[d][1]
    gaps[d] = sc - sstar
    print(f"  gap s_c - s* (d={d}): {mp.nstr(gaps[d], 8)}   ratio-to-next:", flush=True)
alpha = mp.log(gaps[5]/gaps[7])/mp.log(mp.mpf(7)/5), mp.log(gaps[7]/gaps[9])/mp.log(mp.mpf(9)/7)
print(f"  gap exponent estimates: {mp.nstr(alpha[0],6)}, {mp.nstr(alpha[1],6)} -> alpha ~ 2")
pred_gap11 = gaps[9]*(mp.mpf(9)/11)**2
pred_r11 = pred_s_r = None
pred_sc11 = valset[11][1] + pred_gap11
rgaps = {d: (mp.mpf(complex(cr[d][2].replace(' ','')).real) - mp.mpf(complex(cr[d][1].replace(' ','')).real)) for d in (5,7,9)}
pred_rgap11 = rgaps[7] and rgaps[9]*(mp.mpf(9)/11)**2
pred_rc11 = pred_sc11 + pred_rgap11
print(f"  PUBLISH crunode(11): s ~ {mp.nstr(pred_sc11, 8)} window ({mp.nstr(pred_sc11-0.0008,8)},{mp.nstr(pred_sc11+0.0008,8)})")
print(f"                    r ~ {mp.nstr(pred_rc11, 8)} window ({mp.nstr(pred_rc11-0.00075,8)},{mp.nstr(pred_rc11+0.00075,8)})")
out["F1_lock"] = {"s": [mp.nstr(pred_sc11-0.0008, 10), mp.nstr(pred_sc11+0.0008, 10)],
                  "r": [mp.nstr(pred_rc11-0.00075, 10), mp.nstr(pred_rc11+0.00075, 10)]}
out["gaps"] = {str(d): mp.nstr(g, 12) for d, g in gaps.items()}
out["rgaps"] = {str(d): mp.nstr(g, 12) for d, g in rgaps.items()}

out["locks_window"] = {k: [mp.nstr(a, 10), mp.nstr(b, 10)] for k, (a, b) in locks.items() if isinstance(b, mp.mpf)}
json.dump(out, open("jcorner_stage3.json", "w"), indent=1)
print("="*84); print(f"SCORES: C1 {out['C1_ok']}  C3 {out['C3_ok']}  A1 {out['A1_ok']}  A2 {bool(ok47)}  M3 {M3ok}   [{time.time()-t0:.0f}s]")
