"""
Note 14, stage E: the certificate bundle for chamber n=11 + the 10-chamber table.
  (a) det JF = bc = 1 for F10 at 5/5 exact rational points (explainer map, a = -217/162)
  (b) Sturm-EXACT real-cusp counts d = 2..10 (reality dance 1,2,1,2,1,2,1,2,1)
  (c) class check for d=10: tangent pencil deg = 11 (leading (1-n)*a_n != 0)
  (d) the refined TERM LAW certificate for n = 10, 11 (from stored walls):
      support = cone \ {const, s, r^{n-2}, r^n [fictitious: deg_r = n-1]} + {s^n}
  (e) assemble the 10-row fingerprint/content/K/Sturm table (jsons of stages A,B).
"""
import sympy as sp, json, random, time
from sympy import symbols, diff, integrate, expand, Rational as R

t0 = time.time()
out = {}
x, y, z, w = symbols("x y z w")

# ---------- (a) det ----------
def det_at_points(d, npts=5, seed=99):
    p = expand(2*w - 3*w**2 + w*(1-w)*(w**(sp.Integer(d)-2) - R(6, d*(d+1))))
    c = sp.Integer(1); b = sp.Integer(1)
    q = expand(integrate(w*diff(p, w)/c, w))
    kap = R(diff(p, w).subs(w, 1))/c
    a = R(-(1+kap)/(2+kap))
    u = 1 + x*y
    g = 1 + a*x*y + b*x**2*z
    ws = u*g
    alpha = u + q.subs(w, ws)/g**2
    beta = c + p.subs(w, ws)/g
    M = sp.Matrix([alpha/x**2, beta/x, x*g]).jacobian([x, y, z])
    random.seed(seed + d)
    ok = 0
    for _ in range(npts):
        while True:
            pt = {x: sp.Rational(random.randint(-4,4), random.randint(1,5)),
                  y: sp.Rational(random.randint(-4,4), random.randint(1,5)),
                  z: sp.Rational(random.randint(-4,4), random.randint(1,5))}
            if pt[x] != 0 and sp.simplify(g.subs(pt)) != 0: break
        vals = [entry.subs(pt) for entry in M]
        det = sp.Matrix(3, 3, lambda ii, jj: vals[3*ii+jj]).det()
        if sp.simplify(det - 1) == 0: ok += 1
    return ok

ok = det_at_points(10)
print(f"(a) d = 10 (fiber 11): det JF = 1 at {ok}/5 exact rational points  [{time.time()-t0:.0f}s]", flush=True)
out["det_10"] = ok

# ---------- (b) Sturm dance (exact identical implementation as note 10, proven) ----------
def sturm_real_count(f, x):
    seq = [expand(f), expand(diff(f, x))]
    while sp.degree(seq[-1], x) > 0 or seq[-1].is_number is False:
        rr = sp.rem(seq[-2], seq[-1], x)
        if rr == 0 or rr.is_zero: break
        seq.append(expand(-rr))
        if sp.degree(seq[-1], x) == 0: break
    def sig(inf_sign):
        out = []
        for g in seq:
            lc = sp.LC(sp.Poly(g, x)); dg = sp.degree(g, x)
            out.append(sp.sign(lc) * (inf_sign if dg % 2 else 1) if dg > 0 else sp.sign(g))
        return out
    def variations(L):
        L = [v for v in L if v != 0]
        return sum(1 for a, b in zip(L, L[1:]) if a*b < 0)
    return variations(sig(-1)) - variations(sig(1))

dance = []
for d in range(2, 11):
    p = expand(2*w - 3*w**2 + w*(1-w)*(w**(sp.Integer(d)-2) - R(6, d*(d+1))))
    cnt = sturm_real_count(diff(p, w), w)
    sqf = sp.degree(sp.gcd(diff(p,w), diff(p,w,2)), w) == 0
    dance.append((d, cnt, sqf))
    print(f"(b) d={d:2d}: Sturm real cusps = {cnt}   p' squarefree: {sqf}", flush=True)
out["sturm_dance"] = dance

# ---------- (c) class check d=10 ----------
n = 11
p10 = -w**10 + w**9 - R(162,55)*w**2 + R(107,55)*w
Phi11 = expand(integrate(p10, w))
u0, v0 = sp.symbols("u0 v0")
gpen = expand(Phi11 + u0*p10 - w*p10 - v0)          # tangent-through-(u0,v0) equation
lc = sp.LC(sp.Poly(gpen, w))
class_ok = (sp.degree(gpen, w) == n) and (lc == (n-1)*abs(R(-1, n)))
print(f"(c) d=10: tangent pencil deg = {sp.degree(gpen, w)} (predict 11), lc = {lc} = (1-n)*lead(Phi) != 0: {class_ok}")
out["class_10"] = bool(class_ok)

# ---------- (d) refined TERM LAW ----------
s, r = sp.symbols("s r")
law_cert = {}
for dd, fn in ((9, "atlas9_wall.txt"), (10, "atlas10_wall.txt")):
    nn = dd + 1
    D = sp.sympify(open(fn).read())
    P = sp.Poly(D, s, r)
    supp = {m[0] for m in P.terms()}
    cone = {(i,j) for i in range(nn+1) for j in range(nn+1) if (nn-1)*j + nn*i <= nn*(nn-1)}
    pred_support = (cone - {(0,0), (1,0), (0,nn-2), (0,nn)}) | {(nn,0)}
    okk = (supp == pred_support)
    law_cert[nn] = okk
    print(f"(d) n={nn}: support == cone\\{{const,s,r^(n-2),r^n}} + s^n  EXACT: {okk}   terms {len(supp)} = {nn*(nn+1)//2 - 2}")
out["term_law_cert"] = law_cert

# ---------- (e) table ----------
A = json.load(open("atlas10_stageA.json")); B = json.load(open("atlas10_bitangents.json"))
out["stage_summary"] = {"terms": A["terms"], "holes": A["holes"], "s11": A["s11"],
    "content11": A["content11"], "real_cusps": A["real_cusps"], "n_cusps": A["n_cusps"],
    "irreducible": A["irreducible"], "param": A["param_on_wall"], "nodes": B["n_nodes"],
    "crunodes": B["crunodes"], "acnodes": B["acnodes"], "delta": B["delta"],
    "K_scale": B["eliminant_scale"], "sqfree": B["cofactor_squarefree"], "coprime": B["cofactor_coprime_p10p"]}
json.dump(out, open("atlas10_theorem.json", "w"), indent=1, default=str)
print(f"saved atlas10_theorem.json  [{time.time()-t0:.0f}s]")
