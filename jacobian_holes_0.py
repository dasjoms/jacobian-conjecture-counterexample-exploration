"""
Note 11, stage 0 - THE THIRD HOLE, settled before the decic wall is computed.
Wall D(s,r) = prim(resultant(h,h',w)), h = Phi(w) - s w + r.
Hole census (notes 6-10):  constant=0 and s^1=0 ALWAYS;  r^{n-2} = 0 for n>=7 only.

DEFICIT CLASSIFICATION (hand proof, to certify): in the abstract discriminant a monomial
prod a_k^{e_k} obeys sum e_k = 2n-2, sum (n-k)e_k = n(n-1). The r^{n-2}s^0-coefficient needs
e_0 = n-2, e_1 = 0, so on {a_2..a_n}: sum e = n and sum (n-k)e = n. Restricting to the
TOWER support (only a_2, a_3, a_{n-1}, a_n != 0), write (alpha,beta3,gamma,delta) for
(e_2, e_3, e_{n-1}, e_n):  delta = (n-3)alpha + (n-4)beta3 >= 0 forces exactly:
  (0,0,n,0):   a_{n-1}^n * r^{n-2}
  (1,0,2,n-3): a_2 a_{n-1}^2 a_n^{n-3} * r^{n-2}
  (0,1,3,n-4): a_3 a_{n-1}^3 a_n^{n-4} * r^{n-2}
for ALL n >= 7 (others give negative gamma). So beta := [r^{n-2}] disc = C1(n) a_{n-1}^n
+ C2(n) a_2 a_{n-1}^2 a_n^{n-3} + C3(n) a_3 a_{n-1}^3 a_n^{n-4} with universal C's.
Compute C1,C2,C3(n) from ONE symbolic disc per n (spares f = w^n + w^{n-1} + a3 w^3 + a2 w^2 + R),
evaluate the tower specialization (a_n=-1/n, a_{n-1}=1/(n-1), a_2=1-3/(n(n-1)),
a_3=-(1-2/(n(n-1)))) and PRINT beta(7..14) as EXACT rationals BEFORE the n=10 wall exists.
"""
import sympy as sp, json, time
from sympy import symbols, Rational as R

t0 = time.time()
w, Rr, a2, a3 = symbols("w Rr a2 a3")
out = {"C": {}, "beta_tower": {}}

def beta_poly(n):
    f = w**n + w**(n-1) + a3*w**3 + a2*w**2 + Rr
    disc = sp.discriminant(f, w)
    P = sp.Poly(disc, Rr)
    return sp.expand(P.coeff_monomial(Rr**(n-2)))

print(f"{'n':>3} | {'beta poly (in a2,a3 with a_n=a_{n-1}=1)':<55} | time")
for n in range(7, 13):
    bp = beta_poly(n)
    P = sp.Poly(bp, a2, a3)
    C1 = sp.Rational(P.coeff_monomial(a2**0*a3**0)) if (0,0) in [tuple(m) for m,_ in P.terms()] else sp.Rational(0)
    monos = {tuple(m[0]): sp.Rational(m[1]) for m in P.terms()}
    C1 = monos.get((0,0), sp.Rational(0)); C2 = monos.get((1,0), sp.Rational(0)); C3 = monos.get((0,1), sp.Rational(0))
    other = {k: v for k, v in monos.items() if k not in ((0,0),(1,0),(0,1))}
    out["C"][n] = {"C1": str(C1), "C2": str(C2), "C3": str(C3),
                   "C1f": sp.factorint(abs(C1)) if C1 != 0 else 0,
                   "C2f": sp.factorint(abs(C2)) if C2 != 0 else 0,
                   "C3f": sp.factorint(abs(C3)) if C3 != 0 else 0,
                   "classification_3_monomials": other == {}}
    print(f"{n:>3} | C1={C1} C2={C2} C3={C3} | extra monomials: {other if other else 'NONE (classification OK)'} | {time.time()-t0:.0f}s", flush=True)
    # tower specialization, EXACT
    M = n*(n-1)
    A_n = R(-1, n); A_nm1 = R(1, n-1); A2 = 1 - R(3, M); A3 = -(1 - R(2, M))
    beta_t = C1*A_nm1**n + C2*A2*A_nm1**2*A_n**(n-3) + C3*A3*A_nm1**3*A_n**(n-4)
    out["beta_tower"][n] = str(sp.nsimplify(beta_t))
    print(f"      beta_tower({n}) = {sp.nsimplify(beta_t)}   {'== 0 HOLE OPEN' if beta_t==0 else '!= 0 HOLE CLOSES'}", flush=True)

print("\n== verdict lock ==")
for n in range(7, 13):
    print(f"  n={n}: hole r^" + "{" + f"n-2" + "}" + f" = {out['beta_tower'][n]}")
terms10 = None
b10 = sp.sympify(out["beta_tower"][10])
terms10 = 56 - 2 - (1 if b10 == 0 else 0)
print(f"\nLOCKED: beta(10) = {b10}  -> wall terms(10) predicted = {terms10} "
      f"(cone C(10)=56, minus const,s^1 always{', minus r^8 if hole open' if b10==0 else ''})")
out["terms10_prediction"] = terms10
json.dump(out, open("holes0_beta.json","w"), indent=1)
print(f"[{time.time()-t0:.0f}s]")
