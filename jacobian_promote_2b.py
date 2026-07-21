#!/usr/bin/env python3
# LAB NOTE 16 — STAGE 2b: corrected basis (d = 1/ee + 2), bulletproof zero tests,
# even-d sign theorem, ghost-envelope audit vs note-15 data.
# LOCKS (before compute):
#   M1: odd shadow identity H_odd - (U+1)F_odd is the ZERO poly in (X,U) over Q(d).
#   M2: H_even = (U+1)F_even with
#       F_even = 4U^2+17U+19 - 2(d^2+d)(U+2)^2 + X[(U+1) - d(U+2) - d^2(U+2)^2],
#       AND for every even d>=4: A(U,d)=d^2(U+2)^2+d(U+2)-(U+1) > 0 and
#       B(U,d)=2(d^2+d)(U+2)^2-(4U^2+17U+19) > 0 on U>=-1
#       => F_even = -X*A - B < 0 on (-1, inf) => NO even shadow near the corner.
#   M3: eps-machine in the (d-2)-basis reproduces published u1..u4 and s_a,s_b EXACTLY
#       (sympy simplify == 0 against note-15's closed forms).
#   M4: ghost proven envelope err(d) dominates note-15's actual residuals (d=21..47).
import sympy as sp
from sympy import Rational as Q, symbols, log
import json, time
t0 = time.time()
w,U,X,dd,ee = symbols('w U X dd ee')
q = symbols('q', nonnegative=True)
OUT = {"M1": {}, "M2": {}, "M3": {}, "M4": {}}
c0s = Q(6,1)/(dd*(dd+1))
par = symbols('par')
tU = -(1+U)
p_expr = 2*tU - 3*tU**2 - c0s*tU + c0s*tU**2 + (-par)*(1+U)*X - par*(1+U)**2*X
Phi_expr = tU**2 - tU**3 - c0s*tU**2/2 + c0s*tU**3/3 + (par*(1+U)**2*X)/dd + par*(1+U)**3*X/(dd+1)
H = sp.expand((Phi_expr - (tU-1)*p_expr) * dd * (dd+1))

def poly_zero(expr):
    P = sp.Poly(sp.expand(expr), X, U)
    return all(cf == 0 for cf in P.coeffs())

F_odd = (dd**2+dd)*(U+2)**2*(X-2) + 4*U**2+17*U + 19 - X*(U+1)*(1+dd*(U+2))
F_even = 4*U**2+17*U+19 - 2*(dd**2+dd)*(U+2)**2 + X*((U+1) - dd*(U+2) - dd**2*(U+2)**2)
M1_odd = poly_zero(sp.expand(H.subs(par,-1)) - (U+1)*F_odd)
M2_even_div = poly_zero(sp.expand(H.subs(par,1)) - (U+1)*F_even)
OUT["M1"]["odd_identity_zero_poly"] = M1_odd
OUT["M2"]["even_identity_zero_poly"] = M2_even_div
print("M1 odd:", M1_odd, " M2 even:", M2_even_div, flush=True)

# M2 certificates: A>0, B>0 on U>=0 for dd>=4; and on U=-y in [0,1]
A = dd**2*(U+2)**2 + dd*(U+2) - (U+1)
B = 2*(dd**2+dd)*(U+2)**2 - (4*U**2 + 17*U + 19)
cert = {}
for name, expr, base in (("A(U) dd=4+q", A, 4), ("B(U) dd=2+q", B, 2)):
    PU = sp.Poly(sp.expand(expr.subs(dd, base+q)), U)
    cert[name] = all(sp.LC(sp.Poly(cf, q)) is not None and all(c2 >= 0 for c2 in sp.Poly(cf, q).coeffs()) for cf in PU.coeffs())
# U = -y, y in [0,1]: monotone decreasing w/ min at y=1 for both, values at y=1 positive
y = symbols('y', nonnegative=True)
Ay = sp.expand(A.subs(U, -y)); By = sp.expand(B.subs(U, -y))
vertA = sp.solve(sp.Eq(sp.diff(Ay.subs(dd, 2+q), y), 0), y)[0]
vertB = sp.solve(sp.Eq(sp.diff(By.subs(dd, 2+q), y), 0), y)[0]
def posfrac(z):    # z rational in q: certify num>=0, den>0 for q>=0 via coefficients
    n, dn = sp.fraction(sp.together(z))
    return (all(c >= 0 for c in sp.Poly(n, q).coeffs()) and
            all(c > 0 for c in sp.Poly(sp.expand(dn), q).coeffs()) and bool(dn.subs(q,0) > 0))
cert["A(-y): vertex y* >= 1 for dd>=2"] = posfrac(vertA - 1)
cert["A(-y)_at_y=1"] = sp.expand(Ay.subs({y:1}))           # dd^2+dd > 0
cert["B(-y): vertex form"] = str(sp.factor(vertB))
cert["B(-y): vertex y* >= 1 for dd>=2"] = posfrac(vertB - 1)
cert["B(-y)_at_y=1"] = sp.expand(By.subs({y:1}))            # 2dd^2+2dd-6 >0 for dd>=2: coeffs in q
cert["B(-y)_at_y=1_poscoeff_dd2+q"] = all(c >= 0 for c in sp.Poly(By.subs({y:1, dd:2+q}), q).coeffs())
OUT["M2"]["sign_certificates"] = {k: (str(v) if not isinstance(v,(bool,list)) else v) for k,v in cert.items()}
print("M2 certs:", cert, flush=True)

# ---------------- M3: eps-machine in (d-2)-basis ----------------
d_expr = 1/ee + 2                                     # ee = 1/(d-2)
L2s = log(2)
u1,u2,u3,u4 = symbols('u1 u2 u3 u4')
Us = u1*ee + u2*ee**2 + u3*ee**3 + u4*ee**4
log1p = sp.expand(Us - Us**2/2 + Us**3/3 - Us**4/4)
log1p = sum(log1p.coeff(ee,k)*ee**k for k in range(1,5))
expo = sp.expand(log1p/ee)                            # (d-2) ln(1+U) = ln(1+U)/ee
E0, E1, E2, E3, E4 = [sp.expand(expo).coeff(ee,k) for k in range(5)]
Xs_full = sp.series(sp.exp(E0 + E1*ee + E2*ee**2 + E3*ee**3 + E4*ee**4), ee, 0, 5).removeO().expand()
P2 = sp.series((Us+2)**2, ee, 0, 4).removeO().expand()
d2pd = sp.expand(d_expr**2 + d_expr)                  # = ee^-2 + 5 ee^-1 + 6
Fs = sp.expand((ee**-2 + 5*ee**-1 + 6)*P2*(Xs_full-2) + 4*Us**2 + 17*Us + 19 - Xs_full*(1+Us)*(1 + d_expr*(Us+2)))
Fs = sp.expand(Fs)
eq = {k: Fs.coeff(ee, k) for k in (-2,-1,0,1)}
s1 = sp.solve(sp.Eq(eq[-2],0), u1, dict=True)
u1v = log(2)
u2v = sp.expand(sp.solve(sp.Eq(sp.expand(eq[-1].subs(u1,u1v)),0), u2)[0])
u3v = sp.expand(sp.solve(sp.Eq(sp.expand(eq[0].subs({u1:u1v, u2:u2v})),0), u3)[0])
u4v = sp.expand(sp.solve(sp.Eq(sp.expand(eq[1].subs({u1:u1v, u2:u2v, u3:u3v})),0), u4)[0])
pub2 = (1 + L2s**2)/2
pub3 = -Q(7,2) + L2s**3/6 + 3*L2s/4
pub4 = -31*L2s/8 + L2s**4/24 + L2s**2/2 + Q(355,24)
OUT["M3"]["eq1_branch"] = [str(s) for s in s1]
OUT["M3"]["u2_match"] = bool(sp.simplify(u2v-pub2)==0)
OUT["M3"]["u3_match"] = bool(sp.simplify(u3v-pub3)==0)
OUT["M3"]["u4_match"] = bool(sp.simplify(u4v-pub4)==0)
OUT["M3"]["u_resolved"] = [str(u2v), str(u3v), str(u4v)]
print("M3 u matches:", OUT["M3"]["u2_match"], OUT["M3"]["u3_match"], OUT["M3"]["u4_match"], flush=True)

# s*-law: p(t*(U)) with odd parity + c0 series; s+1 = p+1 at ee^1, ee^2
c0_ser = sp.series(c0s.subs(dd, d_expr), ee, 0, 4).removeO()      # 6 ee^2/(1+ee)... check form
Usol = Us.subs({u1:u1v, u2:u2v, u3:u3v, u4:u4v})
Xsol = sp.expand(Xs_full.subs({u1:u1v, u2:u2v, u3:u3v, u4:u4v}))
pU = 2*tU - 3*tU**2 - c0_ser*tU + c0_ser*tU**2 + (1+U)*X + (1+U)**2*X     # odd parity
s_ser = sp.expand(pU.subs({U: Usol, X: Xsol}) + 1)
s_ser = sp.series(s_ser, ee, 0, 3).removeO().expand()
s_a = sp.expand(s_ser).coeff(ee,1); s_b = sp.expand(s_ser).coeff(ee,2); s_0 = sp.expand(s_ser).coeff(ee,0)
OUT["M3"]["s_lead0_zero"] = bool(sp.simplify(s_0)==0)
OUT["M3"]["s_a_match"] = bool(sp.simplify(s_a-(2-2*L2s))==0)
OUT["M3"]["s_b_match"] = bool(sp.simplify(s_b+(Q(5,2)+2*L2s**2-4*L2s))==0)
OUT["M3"]["s_b_resolved"] = str(sp.expand(s_b))
print("M3 s:", OUT["M3"]["s_lead0_zero"], OUT["M3"]["s_a_match"], OUT["M3"]["s_b_match"], s_b, flush=True)

# ---------------- M4: ghost envelope audit vs note-15 C2 ----------------
J = json.load(open('/home/user/jcorner_final.json'))
tab = J['C2']['table']
def proven_err(d_):
    c = Q(6, d_*(d_+1)); rho = Q(1, d_*d_)
    Ap = (Q(1,3) + rho)**(d_-3)
    Emax = Ap * (Q(2*d_,3) - 1 + Q(1,d_))
    return Emax/(6-2*c)
audit = {}
dom_ok = True
for k, v in tab.items():
    d_ = int(k); resid_actual = float(v[0])
    env = float(proven_err(d_))
    audit[d_] = dict(actual=resid_actual, proven_env=env, dominates=env >= resid_actual)
    if env < resid_actual: dom_ok = False
OUT["M4"]["audit"] = audit
OUT["M4"]["envelope_dominates_actual_all"] = dom_ok
OUT["M4"]["proven_err_scaled_3^d"] = {str(d_): float(proven_err(d_)*3**d_) for d_ in (12,15,21,31,47)}
print("M4 dominates:", dom_ok, flush=True)

json.dump(OUT, open('/home/user/promote_stage2b.json','w'), indent=1, default=str)
print("2b DONE %.1fs" % (time.time()-t0))
