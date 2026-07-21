"""
NOTE 12, stage A: THE LIMIT SHADOW.
Conjecture: p_d -> p_inf = 2w - 3w^2 (the fiber-3 counterexample's seed!) so the
tower's real atlas converges to the cubic's shadow. The even-chamber missed cone
mass -> an EXACT integral I over the cubic shadow, explaining the ~8.7% constant.

PREDICTIONS, LOCKED BEFORE COMPUTATION:
 (P1) cubic certificates: cusp exactly (1/3, 1/27); param-cubic node equation
      tau(t1)=tau(t2), p(t1)=p(t2) degenerates EXACTLY to a nonzero scalar*(u-1/3)^3.
 (P2) limit mass integral I (sigma=1.5): predict I in [8.70%, 8.80%], and
      |I - census(n=10) = 8.726%| < 0.05%; MC cross-check within 3 sigma_MC.
 (P3) param convergence: sup_t |(p_d,tau_d) - (p_inf,tau_inf)| on [-0.99, 1.5]
      decays ~1/d: ratio sup(d=7)/sup(d=9) in [1.1, 1.5]; plus honesty check:
      NO uniform convergence on [-1,1] (sup at w=-1 stays O(1)).
 (P4) corner migration: s(d=5,7,9) = -0.8819, -0.9095, -0.9287 strictly decreasing;
      rho-fit s(d) = s* + A rho^(d/2) gives s* in (-1.05, -0.95);
      contacts: t1 increasing toward 1 (predict t1(9) in (0.9835, 0.9885)),
                t2 increasing toward -1 (predict t2(9) in (-1.10, -1.03)).
 (P5) reality dance Sturm-exact extension d=2..14: real roots of p_d' alternate
      1,2,1,2,... exactly; second cusp (=odd d only) obeys
      t2(11) in (-0.8950, -0.8925), t2(13) in (-0.9005, -0.8965).
 (P6) next-chamber locks (n=12, from the fits): corner s(11)... (printed at end),
      census(11) in [I-0.10%, I+0.10%].
"""
import json, math
import numpy as np
import mpmath as mp
import sympy as sp
from sympy import symbols, expand, integrate, diff, Rational as R

mp.mp.dps = 60
w, t, u, s_, r_ = symbols("w t u s r")
out = {}

def seed(d):
    return expand(2*w - 3*w**2 + w*(1-w)*(w**(sp.Integer(d)-2) - R(6, d*(d+1))))

# ---------------- (P1) cubic certificates ----------------
print("="*84); print("(P1) the cubic limit: exact certificates"); print("-"*84)
p_inf = 2*w - 3*w**2; Phi_inf = expand(integrate(p_inf, w))
tau_inf = sp.expand(t*p_inf.subs(w, t) - Phi_inf.subs(w, t))
cusp_pt = (R(1,3), R(1,27))
print(f"  Phi_inf = {Phi_inf},  tau(t) = {tau_inf}")
chk1 = sp.simplify(p_inf.subs(w, R(1,3)) - cusp_pt[0]) == 0
chk2 = sp.simplify((t*p_inf.subs(w,t) - Phi_inf.subs(w,t)).subs(t, R(1,3)) - cusp_pt[1]) == 0
print(f"  cusp (1/3, 1/27): p(1/3)=1/3 ? {chk1},  tau(1/3)=1/27 ? {chk2}")
# node equation: t2 = 2/3 - u from p(t1)=p(t2), tau equality must collapse
sym_node = sp.expand(tau_inf - (tau_inf.subs(t, R(2,3) - t)))
fact = sp.factor(sym_node)
print(f"  node equation tau(t) - tau(2/3 - t) = {sym_node} = {fact}")
chk3 = sp.simplify(sym_node / (t - R(1,3))**3).is_number
print(f"  proportional to (t - 1/3)^3 EXACTLY: {chk3}   [predict True -> limit has NO node]")
out["P1"] = {"cusp_cert": bool(chk1 and chk2), "node_degenerates": bool(chk3),
             "factor": str(fact)}

# ---------------- (P2) the exact limit-mass integral ----------------
print("="*84); print("(P2) limit missed mass  I = P[s<=1/3, r < tau(t_+(s))]  vs censuses"); print("-"*84)
sig = mp.mpf("1.5")
phi = lambda x: mp.exp(-x**2/(2*sig**2))/(sig*mp.sqrt(2*mp.pi))
PHI = lambda x: mp.mpf("0.5")*(1+mp.erf(x/(sig*mp.sqrt(2))))
tplus = lambda s: (2 + mp.sqrt(4 - 12*s))/6
taupl = lambda s: (lambda tt: tt**2 - 2*tt**3)(tplus(s))
I1 = mp.quad(lambda s: phi(s)*PHI(taupl(s)), [-mp.inf, 0, mp.mpf(1)/3])
I2 = mp.quad(lambda s: phi(s)*PHI(taupl(s)), [-mp.inf, -3, 0, mp.mpf(1)/3])
print(f"  I (quad)         = {mp.nstr(I1, 15)}  = {mp.nstr(100*I1, 10)} %")
print(f"  I (quad, split)  = {mp.nstr(I2, 15)}   (agreement check)")
rng = np.random.default_rng(2026)
NMC = 3_000_000
S = rng.normal(0, 1.5, NMC); Rv = rng.normal(0, 1.5, NMC)
mask = S <= 1/3
env = np.where(mask, np.polyval([-2, 0, 1, 0], (2+np.sqrt(np.clip(4-12*S, 0, None)))/6), np.inf)
mc = np.mean(mask & (Rv < env)); mc_err = 3*np.sqrt(mc*(1-mc)/NMC)
print(f"  I (MC, 3e6)      = {mc:.6f} +- {mc_err:.6f}   [3-sigma]")
cens = {6: 8.69, 8: 8.74, 10: 8.7260}
for n, c in cens.items():
    print(f"  census n={n:2d}: {c:7.4f}%   |I - census| = {abs(float(100*I1)-c):.4f}%   [predict <0.05% for n=10]")
out["P2"] = {"I": mp.nstr(I1, 15), "I_pct": mp.nstr(100*I1, 10), "I_quad_split": mp.nstr(I2, 15),
             "I_MC": mc, "MC_err": mc_err, "census_dev": {str(n): abs(float(100*I1)-c) for n, c in cens.items()}}

# ---------------- (P3) param convergence ----------------
print("="*84); print("(P3) param convergence  t |-> (p_d(t), tau_d(t))  ->  cubic param"); print("-"*84)
def pd_tau(d):
    p = seed(d); Phi = expand(integrate(p, w))
    pf = sp.lambdify(w, p, "numpy"); phif = sp.lambdify(w, Phi, "numpy")
    return pf, phif
p_inff = sp.lambdify(w, p_inf, "numpy"); phi_inff = sp.lambdify(w, Phi_inf, "numpy")
TT = np.linspace(-0.99, 1.5, 20001)
for d in (5, 7, 9):
    pf, phif = pd_tau(d)
    dmax = np.max(np.hypot(pf(TT) - p_inff(TT), (TT*pf(TT) - phif(TT)) - (TT*p_inff(TT) - phi_inff(TT))))
    out["P3"] = out.get("P3", {}); out["P3"][f"sup_d{d}"] = float(dmax)
    print(f"  d={d}: sup_(t in [-0.99,1.5]) ||(p_d,tau_d)-(p_inf,tau_inf)|| = {dmax:.6f}")
r79 = out["P3"]["sup_d7"]/out["P3"]["sup_d9"]
print(f"  ratio d7/d9 = {r79:.3f}   [predict 1.1..1.5, ~1/d decay]")
# honesty: no uniform convergence on [-1,1]
TT2 = np.linspace(-1, 1, 20001)
pf9, phif9 = pd_tau(9)
sup_full = np.max(np.abs(pf9(TT2) - p_inff(TT2)))
print(f"  HONESTY: sup_(t in [-1,1]) |p_9 - p_inf| = {sup_full:.4f}  (uniform conv FAILS at t=-1; coefficientwise holds)")
out["P3"]["ratio79"] = float(r79); out["P3"]["sup_full_p9"] = float(sup_full)

# ---------------- (P4) corner migration ----------------
print("="*84); print("(P4) crunode migration n=6,8,10 (+contacts)"); print("-"*84)
corners, contacts = {}, {}
for d, fn in [(5, "atlas5_bitangents.json"), (7, "atlas7_bitangents.json"), (9, "atlas9_bitangents.json")]:
    j = json.load(open(fn))
    best = None
    for node in j["nodes"]:
        vals = [complex(x.replace(" ","").replace("**","e").replace("*^","e")) if "*^" in x or "j" in x else complex(x) for x in node[:4]]
        t1, t2, sc, rc = vals
        if abs(t1.imag) < 1e-8 and abs(t2.imag) < 1e-8 and abs(sc.imag) < 1e-8 and abs(rc.imag) < 1e-8:
            if best is None or abs(sc.real - rc.real) < abs(best[2] - best[3]):
                best = (t1.real, t2.real, sc.real, rc.real)
    corners[d] = (best[2], best[3]); contacts[d] = (best[0], best[1])
    print(f"  d={d}: corner (s,r) = ({best[2]:.6f}, {best[3]:.6f})   contacts t1={best[0]:.6f}, t2={best[1]:.6f}")
s5, s7, s9 = corners[5][0], corners[7][0], corners[9][0]
print(f"  strictly decreasing: {s5 > s7 > s9}   [predict True]")
# rho-fit on odd-d grid d=5,7,9: s(d) = s* + A rho^(d/2)
import numpy.linalg as la
# two-unknown (s*, A, rho) is 3-param/3-pt: solve by grid on rho
def rho_fit(y0, y1, y2):
    best = None
    for rho in np.linspace(0.001, 0.999, 20000):
        # y0 = s* + A rho^2.5? use exponent d/2 at d=5,7,9 -> rho^(2.5),rho^(3.5),rho^(4.5): rescale base R=rho^(1/2): y = s* + A R^d
        pass
    return best
# simpler: fit s(d) = s* + B*R^d  (3 params, 3 points) via 1D search on R
def three_pt(y5, y7, y9):
    best = None
    for Rg in np.linspace(1e-4, 0.9999, 30000):
        A = np.array([[1, Rg**5], [1, Rg**7], [1, Rg**9]])
        # least squares with 3 pts 3 unknowns: solve exactly
        try: sol = la.solve(A, np.array([y5, y7, y9]))
        except Exception: continue
        sstar, B = sol
        resid = abs(sstar + B*Rg**5 - y5) + abs(sstar + B*Rg**7 - y7) + abs(sstar + B*Rg**9 - y9)
        if best is None or resid < best[0]: best = (resid, sstar, B, Rg)
    return best
resid, sstar, B, Rg = three_pt(s5, s7, s9)
print(f"  fit s(d) = s* + B R^d:  s* = {sstar:.5f}, R = {Rg:.4f}, design-residual {resid:.2e}  [predict s* in (-1.05,-0.95)]")
t1s = [contacts[d][0] for d in (5,7,9)]; t2s = [contacts[d][1] for d in (5,7,9)]
print(f"  t1: {t1s}  increasing toward 1: {t1s[0]<t1s[1]<t1s[2]},  t1(9) in (0.9835,0.9885): {0.9835 < t1s[2] < 0.9885}")
print(f"  t2: {t2s}  increasing toward -1: {t2s[0]<t2s[1]<t2s[2]}, t2(9) in (-1.10,-1.03): {-1.10 < t2s[2] < -1.03}")
out["P4"] = {"corners": {str(d): corners[d] for d in corners}, "contacts": {str(d): contacts[d] for d in contacts},
             "sstar_fit": float(sstar), "R_fit": float(Rg)}

# ---------------- (P5) cusp reality dance extension + migration ----------------
print("="*84); print("(P5) reality dance d=2..14 (Sturm-exact) + second-cusp migration"); print("-"*84)
def sturm_real_count(f, x):
    seq = [sp.expand(f), sp.expand(sp.diff(f, x))]
    while True:
        rr = sp.rem(seq[-2], seq[-1], x)
        if rr == 0: break
        seq.append(sp.expand(-rr))
        if sp.degree(seq[-1], x) == 0: break
    def sig(sgn):
        outv = []
        for g in seq:
            Pg = sp.Poly(g, x)
            try: outv.append(sp.sign(Pg.LC()) * (sgn if Pg.degree() % 2 else 1))
            except Exception: outv.append(sp.sign(g))
        return outv
    def var(L):
        L = [v for v in L if v != 0]
        return sum(1 for a, b in zip(L, L[1:]) if a*b < 0)
    return var(sig(-1)) - var(sig(1))
dance, t2cusp = [], {}
mp.mp.dps = 120
for d in range(2, 15):
    pp = sp.diff(seed(d), w)
    cnt = sturm_real_count(pp, w)
    dance.append(cnt)
    coef = [mp.mpf(str(c)) for c in sp.Poly(pp, w).all_coeffs()]
    f = lambda z: mp.polyval(coef, z)
    df_ = lambda z: mp.polyval([c*(len(coef)-1-i) for i, c in enumerate(coef[:-1])], z)
    real_roots = []
    for seed0 in [mp.mpf(-3)+mp.mpf(i)/25 for i in range(151)]:
        try:
            z = mp.findroot(f, seed0, tol=mp.mpf(10)**-100, maxsteps=60)
            if abs(mp.im(z)) < mp.mpf(10)**-40 and -3 < mp.re(z) < 3 and all(abs(z-q) > mp.mpf(10)**-20 for q in real_roots):
                real_roots.append(mp.re(z))
        except Exception: pass
    real_roots.sort()
    if len(real_roots) > 1: t2cusp[d] = real_roots[0]
    tag = f"{cnt}R"
    print(f"  d={d:2d}: real roots of p' = {tag}   roots: {[mp.nstr(q, 14) for q in real_roots]}")
alt = all(dance[i] == (1 if (i % 2 == 0) else 2) for i in range(13))
print(f"  alternation 1,2,1,2,... exact d=2..14: {alt}   [predict True]")
t211 = float(t2cusp.get(11, 0)); t213 = float(t2cusp.get(13, 0))
print(f"  t2(11) = {mp.nstr(t2cusp[11], 25)}  in (-0.8950,-0.8925): {-0.8950 < t211 < -0.8925}")
print(f"  t2(13) = {mp.nstr(t2cusp[13], 25)}  in (-0.9005,-0.8965): {-0.9005 < t213 < -0.8965}")
out["P5"] = {"dance": dance, "alternation": bool(alt),
             "t2cusp": {str(d): mp.nstr(q, 40) for d, q in t2cusp.items()}}

# ---------------- (P6) next-chamber locks (printed for the queue) ----------------
print("="*84); print("(P6) LOCKS for the n=12 chamber round (do NOT touch before computing d=11)"); print("-"*84)
s11 = sstar + B*Rg**11
r11 = (lambda tt: tt**2 - 2*tt**3)((2 + math.sqrt(4 - 12*s11))/6)
print(f"  corner(11): s = {s11:.5f} (rho-fit), r ~= cubic-envelope(s) = {r11:.5f} (shadow prediction)")
print(f"  census(11): I +- 0.10%  ->  [{float(out['P2']['I_pct'].split()[0])-0.10:.4f}%, {float(out['P2']['I_pct'].split()[0])+0.10:.4f}%]")
print(f"  t1(11) in ({t1s[2]:.4f}, 1): expect increase; corner t2(11) in ({t2s[2]:.4f}, -1): expect increase")
out["P6"] = {"s11_fit": float(s11), "r11_shadow": float(r11)}
json.dump(out, open("limits_stageA.json", "w"), indent=1)
print("\nsaved limits_stageA.json")
