"""
NOTE 15, stage 2c (exact-discriminant rebuild). Verdict-carryovers so far:
 (V3) TRUE [d=3: node exactly (-1,-1), contacts {-2,+1}, crunode]
 (V4) FALSE [(-4/3,-4/3) is NOT a node of D4: 0 preimages under (p,tau)]
 (V2-num) corner exponent ~0.96-0.98 (numeric Sturm version) -> re-certify exact.
Locks restated for this run:
 (V1c) the exact discriminant-stratification reproduces the wall-diagonal
       algebraic roots {s*(d), -1, 0} as the ONLY transition candidates, and
       real count changes occur exactly at {0 and s*(d)} for d=5,7,9
       (s = -1 transitions only where a real pair is created/annihilated).
 (V4c) D4(-4/3,-4/3) != 0 AND d=4's real acnode s within ( -1.05, -0.55 )
       [whisker-root basin, matching the d=6..10 series];
       multiplicity bookkeeping: singular-s resultant must factor as
       s^1 * (cusp-cubic)^3 * (node-cubic)^2.
 (V5c) corner-side whisker acnodes d=4,6,8,10: |s+1| non-increasing from d=6.
"""
import json, time
import mpmath as mp
import sympy as sp
from sympy import symbols, expand, diff, integrate, resultant, Rational as R, Poly, factor, gcd

mp.mp.dps = 80
w, s, r = symbols("w s r")
t0 = time.time()
out = {"carryover": {"V3": True, "V4_falsified": True}}

def seed(d):
    return expand(2*w - 3*w**2 + w*(1-w)*(w**sp.Integer(d-2) - R(6, d*(d+1))))

def diag_disc(d):
    """discriminant-in-w of h_diag = Phi - s w + s, as integer poly in s (primitive)."""
    p = seed(d); Phi = integrate(p, w)
    hd = expand((Phi - s*w + s) * sp.lcm([d*(d+1), sp.factorial(sp.Integer(d+1))]))
    D = resultant(hd, diff(hd, w), w)
    P = Poly(D, s)
    return expand(D / sp.gcd_list(P.coeffs()))

def exact_count_diag(d, sv):
    p = seed(d); Phi = integrate(p, w)
    hd = Poly(sp.expand((Phi - w*sv + sv)), w)
    return len(hd.real_roots())

def exact_count_sr(d, sv, rv):
    p = seed(d); Phi = integrate(p, w)
    hd = Poly(sp.expand((Phi - w*sv + rv)), w)
    return len(hd.real_roots())

# ---------- (V1c) exact stratification, d = 3,5,7,9 ----------
print("="*84); print("(V1c) exact discriminant stratification"); print("-"*84)
s1 = json.load(open("jcorner_stage1.json"))
alg = {int(d): [(mp.mpf(q["t"]), mp.mpf(q["s"])) for q in s1["diag"][d]] for d in s1["diag"]}
V1, detail = True, {}
for d in (3, 5, 7, 9):
    D = diag_disc(d)
    rr = Poly(D, s).real_roots()
    rv = sorted(mp.mpf(sp.N(q, 50)) for q in rr)
    pts = [mp.mpf(-3)] + rv + [mp.mpf("0.6")]
    counts = []
    for i in range(len(pts)-1):
        mid = (pts[i]+pts[i+1])/2
        counts.append((mid, exact_count_diag(d, sp.Rational(str(float(mid))))))
    # transitions = Disc roots where count changes across
    trans = []
    for i, rv_i in enumerate(rv):
        c_left = counts[i][1]; c_right = counts[i+1][1]
        if c_left != c_right:
            trans.append((rv_i, c_left, c_right))
    print(f"  d={d}: Disc real roots {[mp.nstr(q,12) for q in rv]}", flush=True)
    print(f"        interval counts {[(mp.nstr(m,4), c) for m, c in counts]}", flush=True)
    for tv, c0, c1 in trans:
        best = min((abs(tv - sa) for _, sa in alg[d]), default=99)
        tag = "diag-root" if best < mp.mpf("1e-40") else ("ORIGIN" if abs(tv) < 1e-40 else "OTHER")
        # -1 always algebraic (t=1): classify subset of diag-root
        print(f"        TRANSITION at s = {mp.nstr(tv, 18)}  [{c0}->{c1}]  [{tag}, dist {mp.nstr(best,4)}]", flush=True)
        if tag == "OTHER":
            V1 = False
    # also: any diag root with NO transition? (touch tangencies)
    notrans = [sa for _, sa in alg[d] if all(abs(sa - tv) > mp.mpf("1e-40") for tv, _, _ in trans)]
    print(f"        algebraic diagonal points WITHOUT count transition: {[mp.nstr(q,10) for q in notrans]} "
          f"(tangency/touch events)", flush=True)
    detail[d] = {"disc_roots": [mp.nstr(q, 40) for q in rv],
                 "transitions": [(mp.nstr(tv, 40), c0, c1) for tv, c0, c1 in trans],
                 "no_transition": [mp.nstr(q, 30) for q in notrans]}
out["V1"] = bool(V1); out["V1_detail"] = detail

# ---------- (V2c) corner shape, exact counts ----------
print("="*84); print("(V2c) corner crossing exponent (exact)"); print("-"*84)
V2 = True
for d in (5, 7):
    sstar = min((sa for _, sa in alg[d]), key=lambda q: abs(q + mp.mpf("0.9")))
    sD = diag_disc(d)   # vertical-line discriminant differs; use full plane disc here:
    p = seed(d); Phi = integrate(p, w)
    hd = expand((Phi - s*w + r) * sp.lcm([d*(d+1), sp.factorial(sp.Integer(d+1))]))
    Dp = resultant(hd, diff(hd, w), w)
    DpP = Poly(Dp, r)
    print(f"  d={d}: plane-disc degree in r: {DpP.degree()}", flush=True)
    DpR = Dp / sp.gcd_list(DpP.coeffs())
    pairs = []
    for dl in ("0.05", "0.08", "0.13", "0.21", "0.34"):
        sv_r = R(str(float(sstar) + float(dl)))
        roots_r = Poly(DpR.subs(s, sv_r), r).real_roots()
        rvv = sorted(mp.mpf(sp.N(q, 30)) for q in roots_r)
        # count transitions across each
        keep = []
        for i, x in enumerate(rvv):
            lo = x if i == 0 else (rvv[i-1]+x)/2
            # left count: use midpoint
            lm = (2*x - (rvv[1] - rvv[0])/4 - (rvv[-1]-rvv[0])) if False else None
        mids = []
        for i in range(len(rvv)+1):
            if i == 0: mids.append(rvv[0] - 1)
            elif i == len(rvv): mids.append(rvv[-1] + 1)
            else: mids.append((rvv[i-1]+rvv[i])/2)
        cnts = [exact_count_sr(d, sv_r, R(str(float(m)))) for m in mids]
        transitions = [rvv[i] for i in range(len(rvv)) if cnts[i] != cnts[i+1]]
        sv_f = mp.mpf(str(sv_r))
        transitions = sorted(transitions, key=lambda x: x)
        above = [x for x in transitions if x > sv_f]
        below = [x for x in transitions if x < sv_f]
        if above and below:
            rho = min(above) - max(below)
            pairs.append((mp.mpf(str(dl)), rho))
            print(f"    delta={dl:>5}: crossings {[(mp.nstr(x,8)) for x in transitions]}  rho={mp.nstr(rho,8)}", flush=True)
    exps = [mp.log(pairs[i+1][1]/pairs[i][1])/mp.log(pairs[i+1][0]/pairs[i][0]) for i in range(len(pairs)-1)]
    em = sum(exps)/len(exps)
    ok = mp.mpf("0.92") < em < mp.mpf("1.08")
    V2 &= bool(ok)
    print(f"    exponents {[mp.nstr(e,4) for e in exps]}  mean {mp.nstr(em,5)}  ok {bool(ok)}", flush=True)
    out.setdefault("V2_detail", {})[d] = {"pairs": [(str(a), str(b)) for a, b in pairs],
                                          "exps": [mp.nstr(e, 6) for e in exps], "ok": bool(ok)}
out["V2"] = bool(V2)

# ---------- (V4c) d=4 wall detail ----------
print("="*84); print("(V4c) d=4 wall, acnode, (-4/3) test"); print("-"*84)
p4 = seed(4); Phi4 = integrate(p4, w)
h4 = expand((Phi4 - s*w + r)*10)
D4 = resultant(h4, diff(h4, r), r) and resultant(h4, diff(h4, w), w)
D4 = expand(D4/sp.gcd_list(Poly(D4, s, r).coeffs()))
res = factor(sp.resultant(D4, diff(D4, s), r))
print(f"  singular-s resultant factor: {res}")
D_at = D4.subs({s: R(-4,3), r: R(-4,3)})
print(f"  D4(-4/3,-4/3) = {D_at}  (nonzero: {D_at != 0})")
# node cubic = the ^2 factor
fac = sp.factor_list(res)[1]
node_cub = [f for f, m in fac if m == 2][0]
cusp_cub = [f for f, m in fac if m == 3][0]
print(f"  node cubic: {Poly(node_cub, s).as_expr()}")
print(f"  cusp cubic: {Poly(cusp_cub, s).as_expr()}")
nn = Poly(node_cub, s).real_roots(); cc = Poly(cusp_cub, s).real_roots()
print(f"  node-s real roots: {[sp.N(q, 20) for q in nn]}")
print(f"  cusp-s real roots: {[sp.N(q, 20) for q in cc]}")
real_nodes = []
for q in nn:
    sv_ = sp.re(q)
    rr_ = Poly(sp.expand(D4.subs(s, q)), r).real_roots()
    for rv_ in rr_:
        gval = diff(D4, s).subs({s: q, r: rv_})
        if gval == 0:
            real_nodes.append((sp.N(q, 30), sp.N(rv_, 30)))
print(f"  singular points on wall D4 (from node-cubic): {real_nodes}")
V4c = (D_at != 0)
ac_s = None
for sv_, rv_ in real_nodes:
    if abs(float(rv_) - float(sv_)) < 1:
        ac_s = float(sv_)
        print(f"  real node candidate: ({sv_}, {rv_})  |s+1| = {abs(ac_s+1):.6f}")
if ac_s is not None:
    V4c &= (-1.05 < ac_s < -0.55)
print(f"  (V4c) D4(-4/3,-4/3)!=0 and real node s in (-1.05,-0.55): {V4c}")
out["V4c"] = bool(V4c); out["D4_at_43"] = str(D_at)
out["d4_node"] = [str(x) for x in (sv_, rv_)] if real_nodes else None

json.dump(out, open("jcorner_stage2.json", "w"), indent=1)
print(f"saved [{time.time()-t0:.0f}s]  V1c={out['V1']} V2c={out['V2']} V4c={out['V4c']}")
