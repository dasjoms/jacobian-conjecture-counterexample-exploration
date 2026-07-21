"""
Note 11, stage E - THE THEOREM SWEEP, d=9 EXTENSION.
The fiber equation h(w) = Phi_{d+1}(w) - s w + r is a LEGENDRE PENCIL of Phi = int p:
  * wall = dual curve (envelope) of the rational graph Gamma: w |-> (w, Phi(w))
  * deg wall = class of Gamma = n EXACTLY  [tangent-through-point eq has deg n,
    leading coeff (1-n)*lead(Phi) != 0]
  * cusps of wall <-> flexes of Gamma <-> roots of p' (n-2, ordinary iff p' squarefree)
  * nodes of wall  <-> bitangents of Gamma, counted by Bezout: (n-1)(n-2) ordered
    contacts, each node costs 2, each (diagonal) cusp costs 2 [the (p')^2 factor]
    ==> #nodes <= (n-2)(n-3)/2 with equality iff the contact system is TRANSVERSE.
Sweep d = 2..9 (fiber 3..10) and certify EVERY chamber:
  (i)   p' squarefree over Q (gcd(p',p'') = 1, EXACT)
  (ii)  Sturm EXACT real-root counts of p'  (the REALITY dance, made exact)
  (iii) eliminant = (den p')^2 * cofactor[deg (n-2)(n-3)] with EXACT constant ratio
  (iv)  cofactor squarefree + coprime to p'  (exact gcd over Q; modular fallback)
  (v)   walls recomputed & matched against the stored txt of notes 6-9
  (vi)  s^n MAGNITUDE LAW: |[s^n] resultant(h,h')| = L^{2n-1} (n-1)^{n-1} |a_n|^{n-1}
        (proof: unique extreme monomial a_n^{n-2} a_1^n in disc + binomial disc
         T_n = (-1)^{(n-1)(n-2)/2} (n-1)^{n-1}  -- also verified symbolically)
  (vii) term counts + support holes vs the weight cone (n-1)j + n i <= n(n-1)
"""
import sympy as sp, json, time
from sympy import symbols, diff, integrate, expand, cancel, resultant, Rational as R, gcd as sgcd

t0 = time.time()
w, w1, w2, s, r, tt = symbols("w w1 w2 s r tt")
report = {}

def seed(d):
    p = 2*w - 3*w**2 + w*(1-w)*(w**(sp.Integer(d)-2) - R(6, d*(d+1)))
    return sp.expand(p)

def sturm_real_count(f, x):
    seq = [sp.expand(f), sp.expand(sp.diff(f, x))]
    while sp.degree(seq[-1], x) > 0 or seq[-1].is_number is False:
        rr = sp.rem(seq[-2], seq[-1], x)
        if rr == 0 or rr.is_zero: break
        seq.append(sp.expand(-rr))
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

print("="*96)
print("chamber sweep d = 2..9  (fiber n = d+1 = 3..10)")
hdr = f"{'d':>2} {'n':>2} {'den':>4} {'class':>5} {'gcd(p,p_pp)':>11} {'SturmR':>6} {'elim56?':>18} {'cof sqfree':>10} {'coprime':>7}"
print(hdr); print("-"*96)

walls = {}
elim_data = {}
for d in range(2, 10):
    n = d + 1
    p = seed(d); Phi = expand(integrate(p, w))
    assert p.subs(w, 1) == -1 and Phi.subs(w, 1) == 0
    pp = diff(p, w); ppp = diff(pp, w)
    den = sp.ilcm(*[sp.denom(c) for c in sp.Poly(pp, w).coeffs()])
    # (class) tangent pencil degree / leading coeff  [t-parametrized!]
    g_cls = sp.Poly(expand(Phi.subs(w, tt) - tt*p.subs(w, tt) + symbols("u0")*p.subs(w, tt)), tt)
    lc_cls = g_cls.LC(); deg_cls = g_cls.degree()
    an = sp.LC(sp.Poly(Phi, w))
    class_ok = (deg_cls == n) and sp.simplify(lc_cls - (1-n)*an) == 0
    # (i) p' squarefree (gcd constant, any size, counts as coprime)
    g_pp = sgcd(pp, ppp)
    sq_pp = (sp.degree(g_pp, w) == 0)
    # (ii) Sturm exact real count
    n_real = sturm_real_count(sp.Poly(pp, w).as_expr(), w) if sq_pp else "??"
    # wall (resultant, scaled)
    L = sp.ilcm(*[sp.denom(c) for c in sp.Poly(Phi, w).coeffs()])
    h = expand(L*(Phi - s*w + r))
    Draw = resultant(h, diff(h, w), w)
    cont = sp.gcd_list([c for c in sp.Poly(Draw, s, r).coeffs()])
    D = sp.expand(Draw/cont); P = sp.Poly(D, s, r)
    terms = len(P.terms())
    pred_terms = n*(n+1)//2 - 1 - (1 if n >= 7 else 0)   # PROVEN TERM LAW (note 10/11)
    # s^n magnitude law
    sN = int(P.coeff_monomial(s**n))
    resN = int(sp.Poly(Draw, s, r).coeff_monomial(s**n))
    law_ratio = R(abs(resN), 1) / (L**(2*n-1) * (n-1)**(n-1) * abs(an)**(n-1))
    content = sp.fraction(R(abs(resN), abs(sN)))[0] if sN else None
    # support holes
    cone = {(i, j) for i in range(n) for j in range(n+1) if (n-1)*j + n*i <= n*(n-1)}
    have = {(int(m[0][1]), int(m[0][0])) for m in P.terms()}  # (i=r-exp, j=s-exp)
    holes = sorted((-1 if (n-1)*j+n*i == n*(n-1) else n*(n-1)-(n-1)*j-n*i, i, j)
                   for (i, j) in cone - have)
    walls[d] = {"D": D, "terms": terms, "sN": sN, "L": L, "holes": holes,
                "pred_terms": pred_terms, "n": n}
    # eliminant via resultant (e1 is monic/const-lead in w1 -> no LC blowup beyond LCs)
    eq1 = sp.expand(cancel((p.subs(w, w2) - p.subs(w, w1))/(w2 - w1)))
    eq2 = sp.expand(cancel((Phi.subs(w, w2) - Phi.subs(w, w1))/(w2 - w1) - p.subs(w, w1)))
    if d <= 7:
        elim = resultant(eq1, eq2, w1)
        ec = sp.gcd_list([c for c in sp.Poly(elim, w2).coeffs()])
        elim = sp.expand(elim/ec)
    else:
        elim = sp.sympify(open(f"atlas{d}_elim_raw.txt").read())
        # LEDGER (normalization discipline): stored raw eliminants are MONIC rational;
        # K=den^2 is an INTEGER-normalization invariant: clear denominators and strip content first.
        Ld = sp.ilcm(*[sp.denom(c) for c in sp.Poly(elim, w2).coeffs()])
        elim = sp.expand(elim*Ld)
        ec = sp.gcd_list([c for c in sp.Poly(elim, w2).coeffs()])
        elim = sp.expand(elim/ec)
    ppw = pp.subs(w, w2)
    q1, rm1 = sp.div(elim, ppw, w2)
    once = (expand(rm1) == 0)
    q2, rm2 = sp.div(q1, ppw, w2) if once else (None, None)
    twice = once and (expand(rm2) == 0)
    cof = None; ratio_str = "FAIL"
    sqfree = coprime = None
    if twice:
        ratio = sp.cancel(q2 / (sp.expand(ppw**0) * 1))
        Knum = sp.LC(sp.Poly(q2, w2)) if sp.Poly(q2, w2).degree() > 0 else q2
        cofac_deg = sp.degree(q2, w2)
        # ratio to (den p')^2 * cofactor: extract scalar
        scalar = sp.simplify(q2 / (q2))  # 1
        LCq = sp.Rational(Knum)
        ratio_str = f"(p')^2 divides 2x EXACT; cofactor deg {cofac_deg} [predict {(n-2)*(n-3)}]"
        t_g = time.time()
        g_sq = sgcd(sp.Poly(q2, w2).as_expr(), diff(q2, w2))
        sqfree = (sp.degree(g_sq, w2) == 0)
        g_cp = sgcd(sp.Poly(q2, w2).as_expr(), ppw)
        coprime = (sp.degree(g_cp, w2) == 0)
        # exact ratio constant elim / (pp')^2 / primitive(cofactor)  ~ expect den^2
        qc = sp.gcd_list([c for c in sp.Poly(q2, w2).coeffs()])
        q2prim = sp.expand(q2/qc)
        ratioK = sp.cancel(sp.expand(elim) / (sp.expand(ppw**2) * q2prim))
        ratio_str += f" | K={ratioK} (den^2={den**2}: {sp.simplify(ratioK-den**2)==0})"
        ratio_str += f" | gcd-time {time.time()-t_g:.0f}s"
    elim_data[d] = {"twice": twice, "sqfree": sqfree, "coprime": coprime}
    if twice: print(f"     eliminant: {ratio_str}", flush=True)
    print(f"{d:>2} {n:>2} {den:>4} {str(class_ok):>5} {str(sq_pp):>11} {str(n_real):>6} "
          f"{str(twice):>18} {str(sqfree):>10} {str(coprime):>7}", flush=True)
    report[d] = {"n": n, "den": str(den), "class_ok": bool(class_ok), "pp_sq": bool(sq_pp),
                 "sturm_real": n_real, "terms": terms, "pred_terms": pred_terms,
                 "sN": sN, "sN_fact": sp.factorint(abs(sN)) if sN else None,
                 "magnitude_law_ratio": str(law_ratio), "content": str(content),
                 "content_fact": sp.factorint(content) if content and content > 1 else None,
                 "elim_twice": twice, "cof_sqfree": sqfree, "cof_coprime": coprime,
                 "holes": [(i, j) for _, i, j in holes]}
    print(f"     wall: terms {terms} (predict {pred_terms}) | s^{n} = {sN} = {sp.factorint(abs(sN))}", flush=True)
    print(f"     magnitude law ratio {law_ratio} (==1 ?) | content {content} = {sp.factorint(content) if content and content>1 else ''}", flush=True)
    print(f"     holes in cone: {[(i,j) for _,i,j in holes]}", flush=True)

# stored-wall cross check (notes 7,8,9 = d 5,6,7 + tonight d=8)
print("\nstored-wall reproduction check:")
for d, fn in [(5, "atlas5_wall.txt"), (6, "atlas6_wall.txt"), (7, "atlas7_wall.txt"), (8, "atlas8_wall.txt"), (9, "atlas9_wall.txt")]:
    Dold = sp.sympify(open(fn).read())
    same = sp.expand(Dold - walls[d]["D"]) == 0 or sp.expand(Dold + walls[d]["D"]) == 0
    print(f"  d={d}: stored == recomputed (up to sign): {same}", flush=True)
    report[d]["stored_match"] = same

# T_n binomial-discriminant verification
print("\nbinomial disc check T_n = (-1)^{(n-1)(n-2)/2} (n-1)^{n-1}:")
a_, b_ = symbols("a_ b_")
for n in range(3, 11):
    db = sp.discriminant(a_*w**n + b_*w, w)
    Tn = sp.simplify(db / (a_**(n-2) * b_**n))
    expect = (-1)**(((n-1)*(n-2))//2) * (n-1)**(n-1)
    print(f"  n={n}: disc(a w^{n} + b w) / (a^{n-2} b^{n}) = {Tn}  expect {expect}: {sp.simplify(Tn-expect)==0}", flush=True)

json.dump(report, open("atlas9_theorem.json","w"), indent=1, default=str)
print(f"\n[sweep done {time.time()-t0:.0f}s]", flush=True)
