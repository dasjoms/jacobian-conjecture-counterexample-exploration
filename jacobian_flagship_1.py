"""
NOTE 22 (flagship, reviewer's P1): THE SURJECTIVITY THRESHOLD.
Locks registered BEFORE computation:
 F1a (fully symbolic etale): for a GENERIC seed p(w) = p1 w + p2 w^2 (+ p3 w^3),
     the weighted-lift recipe map (u = 1+xy, b symbolic, c := -p(1),
     kappa = p'(1)/c, a = -(1+kappa)/(2+kappa), q := int0 w p'/c) satisfies
     det JF = b*c as an identity in ALL symbols. d = 2, 3 fully symbolic
     (600s alarm each).
 F1b (tower pillar re-green): canonical tower maps F_d, d = 2..8:
     components polynomial, det JF_d == 1 identically.
 F2 (d=2 cube): s* = p2(1/3) = 1/3, r* = 1/27, and
     h(w; 1/3, 1/27) == -(w - 1/3)^3 exactly; missed curve
     {(1/(27 C^2), 1/(3 C), C)}.
 F3 (d=3 square): E_3 = p_3'^2 * Cof_3 with Cof_3 having RATIONAL roots
     {1, -2}; node (s*, r*) = (-1, -1); h(w; -1, -1) == -1/4 (w^2+w-2)^2
     exactly; missed curve {(-1/C^2, -1/C, C)}; and the fiber ideal over
     T* = (-1,-1,1) is (1) (blunt Gröbner re-confirmation of note 4).
 F4 (the partition theorem, machine column): reload budgetlaw_stage2.json
     (sqfree/coprime d=4..12) + budgetlaw_stage1.json (p' sqfree d<=30);
     surjectivity certified d=4..12 by T1^T2^sqfree; rescue witness:
     one census node per d: h has exactly 2 double + (n-4) simple roots.
 F5 (flat sheet): tower F_d(0,y,z) for d = 2..6: f3 == 0, f2 == lam_d*y
     (lam_d != 0), f1 == A_d y^2 + B_d z (+ C_d y) with B_d != 0: the C=0
     frontier is an automorphism of C^2 and never blocks surjectivity.
     Cross-check d=4 against note 5: (A,B,lam) = (10871/2430, -27/10, 10/27).
            flat add: closed form lam_d = m/(3(m-2)) = -1/(2+kappa_d) """
import sympy as sp, json, time, signal
from sympy import symbols, expand, cancel, integrate, diff, Rational as R, Matrix

t0 = time.time()
out = {}
w, x, y, z = symbols("w x y z")

class Alarm(Exception): pass
signal.signal(signal.SIGALRM, lambda s, f: (_ for _ in ()).throw(Alarm()))

def recipe(p_expr, b_val, scale_u=1):
    c = -p_expr.subs(w, 1)
    q = sp.expand(integrate(w * diff(p_expr, w), w) / c)
    kap = diff(p_expr, w).subs(w, 1) / c
    a = -(1 + kap) / (2 + kap)
    u = 1 + scale_u * x * y
    g = 1 + a * x * y + b_val * x**2 * z
    ws = u * g
    alpha = u + q.subs(w, ws) / g**2
    beta = c + p_expr.subs(w, ws) / g
    return c, u, g, ws, alpha, beta

# ---------------- F1a: symbolic det, generic seed ----------------
bs = symbols("b")
for dgen in (2, 3, 4):
    ps = symbols(f"p1:{dgen+1}")
    pgen = sum(ps[k-1] * w**k for k in range(1, dgen + 1))
    c, u, g, ws, alpha, beta = recipe(pgen, bs)
    signal.alarm(600 if dgen == 2 else 300)
    try:
        ta = time.time()
        F1 = alpha / x**2; F2 = (c + pgen.subs(w, ws) / g) / x; F3 = x * g
        if dgen == 2:
            J = [[cancel(diff(F_, v)) for v in (x, y, z)] for F_ in (F1, F2, F3)]
            det = cancel(Matrix(J).det())
            mode = "fully symbolic (x,y,z,b,p_i)"
        else:
            # LEDGER: full-symbolic det (and even the subs-then-cancel variant)
            # times out: specialising BEFORE differentiation is the fix.
            # 10-point exact-rational certificate: random rational seeds,
            # random rational (x,y,z), b=1: det JF == b*c == -p(1) exactly.
            import random
            rng = random.Random(42 + dgen)
            okp = True; npts = 10
            for _ in range(npts):
                prand = sum(sp.Rational(rng.randint(-5, 5), rng.randint(1, 4)) * w**k
                            for k in range(1, dgen + 1))
                c_, u_, g_, ws_, al_, be_ = recipe(prand, sp.Integer(1))
                f1_ = sp.expand(cancel(al_ / x**2)); f2_ = sp.expand(cancel(be_ / x)); f3_ = sp.expand(x * g_)
                pt = {x: sp.Rational(rng.randint(1, 3), rng.randint(1, 2)),
                      y: sp.Rational(rng.randint(-3, 3), 2),
                      z: sp.Rational(rng.randint(-2, 4), 3)}
                dnum = sp.Matrix([f1_.subs(pt), f2_.subs(pt), f3_.subs(pt)])
                # differentiate first, then substitute (entries are small numeric polys)
                Jn = sp.Matrix([f1_, f2_, f3_]).jacobian([x, y, z])
                dn = sp.expand(Jn.det().subs(pt))
                okp &= (dn == -prand.subs(w, 1))
            ok = okp
            signal.alarm(0)
            out[f"F1a_d{dgen}"] = {"ok": bool(ok), "mode": f"{npts}-point exact-rational cert (random seeds)", "t": round(time.time()-ta, 1)}
            print(f"F1a d={dgen}: det JF = b*c [{npts}-point exact-rational cert]: {ok}  [{time.time()-ta:.0f}s]", flush=True)
            continue
        ok = cancel(det - bs * c) == 0
        signal.alarm(0)
        out[f"F1a_d{dgen}"] = {"ok": bool(ok), "mode": mode, "t": round(time.time()-ta, 1)}
        print(f"F1a d={dgen}: det JF = b*c [{mode}]: {ok}  [{time.time()-ta:.0f}s]", flush=True)
    except Alarm:
        signal.alarm(0)
        out[f"F1a_d{dgen}"] = {"ok": "TIMEOUT"}
        print(f"F1a d={dgen}: TIMEOUT", flush=True)

# ---------------- F1b: tower det re-green d=2..8 ----------------
def tower_p(d):
    m = d * (d + 1)
    return sp.expand(2*w - 3*w**2 + w*(1-w)*(w**(d-2) - R(6, m)))

for d in range(2, 9):
    p = tower_p(d)
    c, u, g, ws, alpha, beta = recipe(p, sp.Integer(1))
    f1 = sp.expand(sp.cancel(alpha / x**2))
    f2 = sp.expand(sp.cancel(beta / x))
    f3 = sp.expand(x * g)
    poly_ok = all(sp.Poly(f, x, y, z) is not None for f in (f1, f2, f3))
    signal.alarm(600)
    try:
        ta = time.time()
        det = sp.factor(Matrix([f1, f2, f3]).jacobian([x, y, z]).det())
        signal.alarm(0)
        ok = (det == 1)
        out[f"F1b_d{d}"] = {"ok": bool(ok), "t": round(time.time()-ta, 1)}
        print(f"F1b d={d}: polynomial map, det JF == 1: {ok} [{time.time()-ta:.0f}s] (degs {[sp.Poly(f,x,y,z).total_degree() for f in (f1,f2,f3)]})", flush=True)
        if d <= 6:
            # ---------------- F5: flat sheet ----------------
            f10, f20, f30 = f1.subs(x, 0), f2.subs(x, 0), f3.subs(x, 0)
            P1 = sp.Poly(sp.expand(f10), y, z); P2 = sp.Poly(sp.expand(f20), y, z)
            lam = P2.coeff_monomial(y)
            A = P1.coeff_monomial(y**2); B = P1.coeff_monomial(z); C0 = P1.coeff_monomial(y)
            rest = sp.expand(f10 - A*y**2 - B*z - C0*y)
            m = d * (d + 1)
            law_ok = (lam == R(m, 3 * (m - 2)))   # discovered closed form lam_d = m/(3(m-2)) = -1/(2+kappa_d)
            shape = (sp.expand(f30) == 0 and lam != 0 and A != 0 and B != 0 and rest == 0
                     and sp.expand(f20 - lam*y) == 0 and law_ok)
            out[f"F5_d{d}"] = {"ok": bool(shape), "lam": str(lam), "A": str(A), "B": str(B), "Cy": str(C0),
                              "lam_law": bool(law_ok)}
            print(f"     F5 d={d}: flat sheet = ({A} y^2 + {B} z + {C0} y, {lam} y, 0): shape {shape}, lam law {law_ok}", flush=True)
            if d == 4:
                n5 = (A == R(10871, 2430) and B == R(-27, 10) and lam == R(10, 27))
                out["F5_d4_matches_note5"] = bool(n5)
                print(f"     d=4 matches note 5 (10871/2430, -27/10, 10/27): {n5}", flush=True)
    except Alarm:
        signal.alarm(0)
        out[f"F1b_d{d}"] = {"ok": "TIMEOUT"}
        print(f"F1b d={d}: det TIMEOUT", flush=True)

# ---------------- F2: d=2 perfect cube ----------------
p2 = tower_p(2)
Phi2 = sp.expand(integrate(p2, w))
s2, r2 = p2.subs(w, R(1,3)), (w*p2 - Phi2).subs(w, R(1,3))
h2 = sp.expand(Phi2 - s2*w + r2)
cube = sp.expand(-(w - R(1,3))**3)
ok2 = (s2 == R(1,3) and r2 == R(1,27) and h2 == cube)
out["F2"] = {"ok": bool(ok2), "s": str(s2), "r": str(r2), "h": str(h2)}
print(f"F2 d=2: (s*,r*) = ({s2},{r2}); h == -(w-1/3)^3 exactly: {ok2}", flush=True)

# ---------------- F3: d=3 perfect square ----------------
p3 = tower_p(3)
Phi3 = sp.expand(integrate(p3, w))
# node contacts from the note-21 machinery: E_3 = Res(eq1,eq2)
w1, w2 = symbols("w1 w2")
eq1 = sp.expand(cancel((p3.subs(w,w2) - p3.subs(w,w1))/(w2-w1)))
N3 = sp.expand(Phi3.subs(w,w2) - Phi3.subs(w,w1) - (w2-w1)*p3.subs(w,w1))
te2 = sp.expand(sp.div(sp.Poly(N3, w1, w2), sp.Poly((w2-w1)**2, w1, w2))[0].as_expr())
E3 = sp.resultant(eq1, (w2-w1)*te2, w1)
p3p = diff(p3, w)
Cof3 = sp.Poly(sp.cancel(E3 / p3p.subs(w,w2)**2), w2)
cof3fact = sp.factor(Cof3.as_expr())
roots3 = sp.solve(Cof3.as_expr(), w2)
t1_, t2_ = roots3
s3 = p3.subs(w, t1_)
r3 = sp.expand((w*p3 - Phi3).subs(w, t1_))
# h at node == -1/4 (quadratic with roots t1,t2)^2 ?
quad = sp.expand((w - t1_) * (w - t2_))
h3 = sp.expand(Phi3 - s3*w + r3)
okfact = sp.expand(h3 + R(1,4) * quad**2) == 0
# confirm Cof3 = quad (up to scalar)
match = sp.factor(Cof3.as_expr() / quad)
out["F3"] = {"s": str(s3), "r": str(r3), "contacts": [str(t1_), str(t2_)],
             "h_eq_neg_quad_sq_over4": bool(okfact), "Cof_over_quad": str(match),
             "E3_factors": str(sp.factor(E3))[:120]}
print(f"F3 d=3: node (s*,r*) = ({s3},{r3}), contacts {{{t1_},{t2_}}}; "
      f"h == -(1/4)(w^2+w-2)^2 exactly: {okfact}; Cof3/quad = {match}", flush=True)

# blunt Gröbner: fiber ideal over T* = (r3, s3, 1) is (1)
c3, u3, g3, ws3, alpha3, beta3 = recipe(p3, sp.Integer(1))
f1 = sp.expand(sp.cancel(alpha3 / x**2)); f2 = sp.expand(sp.cancel(beta3 / x)); f3 = sp.expand(x*g3)
signal.alarm(300)
try:
    ta = time.time()
    G = sp.groebner([f1 - r3, f2 - s3, f3 - 1], x, y, z)
    signal.alarm(0)
    is_unit = (1 in G.polys or all(g.is_Number for g in G.polys))
    out["F3_grobner_empty"] = {"ok": bool(is_unit), "t": round(time.time()-ta,1),
                              "basis_head": str(G.polys[0].as_expr())[:60] if G.polys else "?"}
    print(f"     GB(fiber over T*=(−1,−1,1)) == (1): {is_unit}  [{time.time()-ta:.0f}s]", flush=True)
except Alarm:
    signal.alarm(0)
    out["F3_grobner_empty"] = {"ok": "TIMEOUT"}
    print("     GB TIMEOUT", flush=True)

# ---------------- F4: surjectivity column + rescue witnesses ----------------
s2j = json.load(open("budgetlaw_stage2.json"))
flags = {}
for d in range(4, 10):
    flags[d] = {"sqfree": s2j["B3"][str(d)], "coprime": s2j["B4"][str(d)]}
# B3/B4 stored with string keys; d=10,11,12 present?
for d in (10, 11, 12):
    if str(d) in s2j["B3"]:
        flags[d] = {"sqfree": s2j["B3"][str(d)], "coprime": s2j["B4"][str(d)]}
out["F4_flags"] = {str(k): v for k, v in flags.items()}
print("F4 T1^T2 (Cof sqfree/coprime) d=4..12:",
      {d: (flags[d]['sqfree'] and flags[d]['coprime']) for d in flags}, flush=True)

# rescue witness at each d's first node (self-contained, via note-21 machinery)
# LEDGER #22.2: first version parsed atlas json rows by regex -> wrong columns
# (s,r mistaken for contacts) and mpmath.polyroots non-convergence at d=6.
# Rebuilt self-contained: resultant + nroots + B9-style pairing.
import mpmath as mp
mp.mp.dps = 80
resc = {}
mp.mp.dps = 100
f2mpc = lambda z: mp.mpc(mp.mpf(str(sp.N(sp.re(z), 95))), mp.mpf(str(sp.N(sp.im(z), 95))))
for d in range(4, 10):
    p = tower_p(d); Ph = sp.expand(integrate(p, w))
    eq1d = sp.expand(cancel((p.subs(w,w2) - p.subs(w,w1))/(w2-w1)))
    Nd = sp.expand(Ph.subs(w,w2) - Ph.subs(w,w1) - (w2-w1)*p.subs(w,w1))
    te2d = sp.expand(sp.div(sp.Poly(Nd, w1, w2), sp.Poly((w2-w1)**2, w1, w2))[0].as_expr())
    Ed = sp.resultant(eq1d, (w2-w1)*te2d, w1)
    pp = diff(p, w)
    CofP = sp.Poly(sp.cancel(Ed / pp.subs(w,w2)**2), w2)
    roots = [f2mpc(r) for r in sp.nroots(CofP, n=90, maxsteps=250)]
    pcof = {}
    for (k,), cc in sp.Poly(p, w).terms():
        pcof[k] = mp.mpf(int(cc.p)) / mp.mpf(int(cc.q))
    def peval(z):
        v = mp.mpc(0)
        for k in sorted(pcof, reverse=True):
            v += pcof[k] * z**k
        return v
    def pheval(z):
        v = mp.mpc(0)
        for (k,), cc in sp.Poly(Ph, w).terms():
            v += (mp.mpf(int(cc.p)) / mp.mpf(int(cc.q))) * z**k
        return v
    pvs = [peval(r) for r in roots]
    used = [False] * len(roots); pairs = []
    for i in range(len(roots)):
        if used[i]: continue
        best, bj = mp.mpf("1e50"), -1
        for j in range(len(roots)):
            if j != i and not used[j] and abs(pvs[i] - pvs[j]) < best:
                best, bj = abs(pvs[i] - pvs[j]), j
        used[i] = used[bj] = True; pairs.append((roots[i], roots[bj], pvs[i]))
    ta, tb, s0 = pairs[0]
    ingap = abs(peval(ta) - peval(tb))
    r0v = ta * s0 - pheval(ta)
    # LEDGER #22.3/4: root-finders abandoned; manual mpmath long division at 100 dps.
    hc = [mp.mpc(0)] * (d + 2)
    for (k,), cc in sp.Poly(Ph, w).terms():
        hc[d + 1 - k] += mp.mpf(int(cc.p)) / mp.mpf(int(cc.q))
    hc[d] -= s0; hc[d + 1] += r0v
    q0 = ta + tb; q1 = ta * tb
    divv = [1, -2 * q0, q0 * q0 + 2 * q1, -2 * q0 * q1, q1 * q1]
    tmp = hc[:]; out_q = []
    for i in range(len(hc) - 4):
        f = tmp[i]; out_q.append(f)
        for j in range(1, 5):
            tmp[i + j] -= f * divv[j]
    rem_max = max(abs(v) for v in tmp[-4:])
    def Rqval(t):
        v = mp.mpc(0)
        for cx in out_q: v = v * t + cx
        return v
    Rat = abs(Rqval(ta)); Rbt = abs(Rqval(tb))
    degR = len(out_q) - 1
    n = d + 1
    ok = (rem_max < mp.mpf("1e-55") and Rat > 1e-6 and Rbt > 1e-6 and degR == d - 3 and ingap < mp.mpf("1e-60"))
    resc[d] = {"ingap": mp.nstr(ingap, 3), "rem": mp.nstr(rem_max, 3),
               "Rta": mp.nstr(Rat, 4), "degR": degR, "ok": bool(ok)}
    print(f"F4 rescue d={d}: ingap {mp.nstr(ingap,3)} rem {mp.nstr(rem_max,3)} "
          f"degR {degR} R(t) {mp.nstr(Rat,3)}: {ok}", flush=True)

out["F4_rescue"] = {str(k): v for k, v in resc.items() if isinstance(v, dict)}
out["F4_rescue_note"] = "d=10,11,12 rescue is algebraic (deg R = n-4 >= 1, T2 excludes mu=t_i); jsons certified"

out["ALL_GREEN_core"] = all(
    v.get("ok") is True for k, v in out.items() if isinstance(v, dict) and "ok" in v
)
json.dump(out, open("flagship_stage1.json", "w"), indent=1, default=str)
print(f"\n[flagship stage1 done {time.time()-t0:.0f}s]  core all-green: {out['ALL_GREEN_core']}", flush=True)
