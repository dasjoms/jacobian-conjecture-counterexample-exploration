#!/usr/bin/env python3
# LAB NOTE 16 — STAGE 2: contractions (L3,L4), symbolic shadow identity (L5), series (L6)
# =====================================================================
# PREDICTIONS LOCKED BEFORE COMPUTE (2026-07-21, off the back of stage 1 L1/L2 green):
#   P3: ghost contraction T(v)=[-c0/3+E(v)]/(6-2c0), E(v)=(1/3+v)^(d-2)(2d/3-1-dv):
#       for every d=12..200: L_g(d) < 1/2, invariance margin>0, and the proven
#       envelope err(d) = Emax/(6-2c0) satisfies  err(d)*3^d <= d^1 * 1  (roughly O(d)).
#       Existence double-door: sign(p'(1/3-rho)) != sign(p'(1/3+rho)), rho=1/d^2, d=12..200.
#       d=5..11: report how early contraction starts (expected ~d=6-8).
#   P4: left cusp, odd d=7..99, I=[0.85 u0, 1.15 u0], u0=ln((2d-1)/8)/(d-2):
#       (i) R'<0 throughout I; (ii) 0<R_lo<=R_hi<=1; (iii) L_L(d) < 1 (expect ~<0.07 by d=15);
#       (iv) H sign change in I (existence, no invariance needed);
#       (v) invariance envelope R(u_hi)>=(1-u_hi)^(d-2) and R(u_lo)<=(1-u_lo)^(d-2);
#       (vi) 8-iteration bound L_L^8*diam < 1e-12 for d>=11.
#   P5: H = d(d+1) G(t=-(1+U)) with odd parity is EXACTLY (U+1)*F_odd(d,U,X)
#       identically in symbols (d,U,X); even parity gives (U+1)*F_even: NEW identity.
#   P6: eps re-solve reproduces published u2,u3,u4 EXACTLY (sympy zero diff);
#       and s-law: s_a = 2-2ln2, s_b = -(5/2+2ln^2 2-4ln2) exactly.
#       Even-d shadow: positive root of F_even(.,X=(1+U)^(d-2)) exists for d=4..30.
import sympy as sp
from sympy import Rational as Q, symbols, expand, log
import mpmath as mp
import json, time
mp.mp.dps = 60
t0 = time.time()
w, U, X, dd, ee = symbols('w U X dd ee')
OUT = {"L3": {}, "L4": {}, "L5": {}, "L6": {}}
def c0(d_): return Q(6, d_*(d_+1))

# ---------------- L3: GHOST CONTRACTION (exact rationals) ----------------
def ghost_stats(d_):
    rho = Q(1, d_*d_)
    c = c0(d_)
    den = 6 - 2*c
    Ap = (Q(1,3) + rho)**(d_-3)                        # (1/3+rho)^{d-3} exact
    Emax = Ap * (Q(2*d_,3) - 1 + Q(1,d_))              # |E| bound
    Epmax = Ap * (Q(2*d_*d_,3) - Q(8*d_,3) + 2 + Q(d_-1, d_))  # |E'| bound (d>=4)
    Lg = Epmax / den
    inv_num = (c/3 + Emax) / den                        # |T| worst inside ball
    margin = rho - inv_num
    err = Emax / den                                    # proven |v + c0/(3(6-2c0))| bound
    return dict(rho=rho, Lg=Lg, margin=margin, err=err,
                err_scaled=err * Q(3)**d_)
L3rows = {}
first_contract = None
sign_ok = True
for d_ in range(5, 201):
    st = ghost_stats(d_)
    contracted = bool(st['Lg'] < Q(1,2)) and bool(st['margin'] > 0)
    if d_ >= 12 and not contracted:
        L3rows[d_] = "FAIL"
    if first_contract is None and contracted: first_contract = d_
    if d_ in (5,6,7,8,9,10,11,12,15,20,30,47,100,200):
        L3rows[d_] = dict(Lg=float(st['Lg']), margin=float(st['margin']),
                          err_scaled=float(st['err_scaled']), contracted=contracted)
for d_ in range(12, 201):
    pp = (2 - c0(d_)) + (2*c0(d_)-6)*w + (d_-1)*w**(d_-2) - d_*w**(d_-1)
    lo = pp.subs(w, Q(1,3) - Q(1,d_*d_)); hi = pp.subs(w, Q(1,3) + Q(1,d_*d_))
    if not (bool(lo > 0) and bool(hi < 0)):  # p'(1/3+v)=6?  sign convention: p'(1/3-rho)>0? ghost left of 1/3
        sign_ok = False; L3rows[f'signfail_{d_}'] = (str(lo), str(hi))
OUT["L3"]["rows"] = L3rows
OUT["L3"]["first_contract_d"] = first_contract
OUT["L3"]["sign_door_ok_d12..200"] = sign_ok
OUT["L3"]["all_d12..200_contracted"] = not any(isinstance(v,str) for v in L3rows.values()) and not any(str(k).startswith('signfail') for k in L3rows)
print("L3 done %.1fs first_contract=%s sign_ok=%s" % (time.time()-t0, first_contract, sign_ok), flush=True)

# ---------------- L4: LEFT CUSP CONTRACTION (exact over rational envelope) ----------------
def left_stats(d_):
    c = c0(d_)
    u0 = mp.log(mp.mpf(2*d_-1)/8)/(d_-2)
    half = mp.mpf(1)/(2**80)
    u0_lo = Q(int(mp.floor((u0-half)*(2**80))), 2**80)
    u0_hi = Q(int(mp.floor((u0+half)*(2**80)))+1, 2**80)
    u_lo = u0_lo * Q(17,20); u_hi = u0_hi * Q(23,20)
    N  = lambda u: 8 - 6*u - c*(3 - 2*u)
    Dn = lambda u: (2*d_ - 1) - d_*u
    M  = lambda u: (-6 + 2*c)*Dn(u) + d_*N(u)          # R' numerator (affine)
    m_lo, m_hi = M(u_lo), M(u_hi)
    monotone_neg = bool(m_lo < 0) and bool(m_hi < 0)
    R_lo, R_hi = N(u_hi)/Dn(u_hi), N(u_lo)/Dn(u_lo)
    in01 = bool(R_lo > 0) and bool(R_hi <= 1)
    Rpmax = max(-m_lo, -m_hi, m_lo, m_hi) / Dn(u_hi)**2 # |R'| max via endpoints
    L_L = Rpmax / ( (d_-2) * R_lo )
    H = lambda u: (1 - u)**(d_-2) - N(u)/Dn(u)
    sign_change = bool(H(u_lo) > 0) != bool(H(u_hi) > 0)
    inv_lo = R_hi <= (1 - u_lo)**(d_-2)
    inv_hi = R_lo >= (1 - u_hi)**(d_-2)
    diam = u_hi - u_lo
    err8 = L_L**8 * diam
    return dict(u_lo=u_lo, u_hi=u_hi, monotone_neg=monotone_neg, in01=in01,
                L_L=L_L, sign_change=sign_change, inv=(bool(inv_lo), bool(inv_hi)),
                err8=err8)
L4rows = {}; L4fail = []
for d_ in range(7, 100, 2):
    st = left_stats(d_)
    ok = st['monotone_neg'] and st['in01'] and bool(st['L_L'] < 1) and st['sign_change']
    if not ok: L4fail.append((d_, {k:(str(v)[:40]) for k,v in st.items()}))
    if d_ in (7,9,11,15,21,31,51,99):
        L4rows[d_] = dict(L_L=float(st['L_L']), err8=float(st['err8']),
                          sign_change=st['sign_change'], inv=st['inv'],
                          u_lo=float(st['u_lo']), u_hi=float(st['u_hi']))
OUT["L4"]["rows"] = L4rows
OUT["L4"]["fails"] = L4fail
OUT["L4"]["PASS_odd_7..99"] = (len(L4fail)==0)
print("L4 done %.1fs fails=%d" % (time.time()-t0, len(L4fail)), flush=True)

# ---------------- L5: SYMBOLIC SHADOW IDENTITY, both parities ----------------
c0s = Q(6,1)/(dd*(dd+1))
par = symbols('par')   # par = (-1)^d
tU = -(1+U)
base_p = 2*tU - 3*tU**2 - c0s*tU + c0s*tU**2
base_Phi = tU**2 - tU**3 - c0s*tU**2/2 + c0s*tU**3/3
# w^(d-1) = (-par)*(1+U)*X ; w^d = par*(1+U)^2*X ; w^(d+1) = (-par)*(1+U)^3*X
p_expr = base_p + (-par)*(1+U)*X - par*(1+U)**2*X
Phi_expr = base_Phi + (par*(1+U)**2*X)/dd - ((-par)*(1+U)**3*X)/(dd+1)
G_expr = Phi_expr - (tU-1)*p_expr
H = sp.expand(G_expr * dd * (dd+1))
H_par = sp.expand(H.subs(par, -1))          # odd d
F_odd = (dd**2+dd)*(U+2)**2*(X-2) + 4*U**2 + 17*U + 19 - X*(U+1)*(1+dd*(U+2))
rem_odd = sp.expand(H_par - (U+1)*F_odd)
L5_odd_zero = (rem_odd == 0)
H_E = sp.expand(H.subs(par, 1))             # even d
q_even, r_even = sp.div(sp.Poly(H_E, U), sp.Poly(U+1, U))
L5_even_div = (sp.expand(r_even.as_expr()) == 0)
F_even = sp.expand(q_even.as_expr())
OUT["L5"]["odd_identity_exact_zero"] = L5_odd_zero
OUT["L5"]["even_divisible_by_(U+1)"] = L5_even_div
OUT["L5"]["F_even"] = str(sp.factor(F_even))
print("L5 done %.1fs odd_zero=%s even_div=%s" % (time.time()-t0, L5_odd_zero, L5_even_div), flush=True)

# even-d shadow roots numeric (d=4..30)
even_roots = {}
for d_ in range(4, 31, 2):
    Fe = sp.sympify(F_even.subs(dd, d_))
    f = sp.lambdify((U,X), Fe, 'mpmath')
    try:
        g = lambda u: f(mp.mpf(u), (1+mp.mpf(u))**(d_-2))
        # scan for sign change in (1e-6, 0.5)
        grid = [mp.mpf(10)**(mp.mpf(-6) + i*mp.mpf('0.05')) for i in range(80)]
        vals = [g(x) for x in grid]
        root = None
        for i in range(len(grid)-1):
            if vals[i]*vals[i+1] < 0:
                root = mp.findroot(g, (grid[i], grid[i+1])); break
        even_roots[d_] = str(root) if root is not None else None
    except Exception as e:
        even_roots[d_] = f"ERR {e}"
OUT["L5"]["even_shadow_roots_U"] = even_roots
print("L5 even roots done %.1fs" % (time.time()-t0), flush=True)

# ---------------- L6: eps re-solve, exact match vs published ----------------
L2sym = log(2)
u1, u2, u3, u4 = symbols('u1 u2 u3 u4')
Us = u1*ee + u2*ee**2 + u3*ee**3 + u4*ee**4
log1p = (Us - Us**2/2 + Us**3/3 - Us**4/4)
log1p = sp.expand(log1p)
log1p = sum(c*ee**k for k, c in [(k, log1p.coeff(ee,k)) for k in range(1,5)])
expo = sp.expand((1/ee - 2) * log1p)
expo = sum(sp.expand(expo).coeff(ee,k)*ee**k for k in range(-1,4))
A1 = sp.expand(expo).coeff(ee,-1); B0 = sp.expand(expo).coeff(ee,0)
C1 = sp.expand(expo).coeff(ee,1); D2 = sp.expand(expo).coeff(ee,2); E3 = sp.expand(expo).coeff(ee,3); F4 = sp.expand(expo).coeff(ee,4)
OUT["L6"]["expo_coeff_-1"] = str(A1)   # must be '0': log1p starts at ee^1, so (1/ee-2)*log1p has no ee^-1 term
Xseries = sp.exp(B0 + C1*ee + D2*ee**2 + E3*ee**3)  # drop F4 (needed only for u5)
Xser = sp.series(Xseries, ee, 0, 4).removeO().expand()
P2 = sp.series((Us+2)**2, ee, 0, 4).removeO().expand()
Fs = (ee**-2 + ee**-1)*P2*(Xser-2) + 4*Us**2 + 17*Us + 19 - Xser*(1+Us)*(1 + ee**-1*(Us+2))
Fs = sp.expand(Fs)
eq1 = sp.expand(Fs).coeff(ee, -2)
eq2 = sp.expand(Fs).coeff(ee, -1)
eq3 = sp.expand(Fs).coeff(ee, 0)
eq4 = sp.expand(Fs).coeff(ee, 1)
s1 = sp.solve(sp.Eq(eq1, 0), u1, dict=True)
OUT["L6"]["eq1_solutions"] = [str(s) for s in s1]
# pick u1 = ln 2 branch (e^{u1}=2)
u1v = log(2)
eq2s = sp.simplify(eq2.subs(u1, u1v))
s2 = sp.solve(sp.Eq(eq2s, 0), u2)[0]
u2v = sp.expand(s2)
eq3s = sp.simplify(eq3.subs({u1: u1v, u2: u2v}))
s3 = sp.solve(sp.Eq(eq3s, 0), u3)[0]
u3v = sp.expand(s3)
eq4s = sp.simplify(eq4.subs({u1: u1v, u2: u2v, u3: u3v}))
s4 = sp.solve(sp.Eq(eq4s, 0), u4)[0]
u4v = sp.expand(s4)
pub2 = (1 + L2sym**2)/2
pub3 = -Q(7,2) + L2sym**3/6 + 3*L2sym/4
pub4 = -31*L2sym/8 + L2sym**4/24 + L2sym**2/2 + Q(355,24)
OUT["L6"]["u1"] = str(u1v)
OUT["L6"]["u2_exact_match"] = bool(sp.simplify(u2v - pub2) == 0)
OUT["L6"]["u3_exact_match"] = bool(sp.simplify(u3v - pub3) == 0)
OUT["L6"]["u4_exact_match"] = bool(sp.simplify(u4v - pub4) == 0)
OUT["L6"]["u2_resolved"] = str(u2v); OUT["L6"]["u3_resolved"] = str(u3v); OUT["L6"]["u4_resolved"] = str(u4v)
print("L6 u-solve done %.1fs match=%s%s%s" % (time.time()-t0, OUT["L6"]["u2_exact_match"], OUT["L6"]["u3_exact_match"], OUT["L6"]["u4_exact_match"]), flush=True)

# s*-law via p(-1-U) with odd parity, expand s+1 at orders ee,ee^2 using Xseries+u values
pU = sp.expand(p_expr.subs({par: -1, dd: 1/ee}))
# substitute X series (with solved u) and expand; need pU+1 at ee^1, ee^2
Xfull = sp.series(sp.exp(B0 + C1*ee + D2*ee**2 + E3*ee**3 + F4*ee**4), ee, 0, 5).removeO().expand()
Xsub = Xfull.subs({u1: u1v, u2: u2v, u3: u3v, u4: u4v})
Usol = Us.subs({u1: u1v, u2: u2v, u3: u3v, u4: u4v})
s_expr = sp.expand(pU.subs({U: Usol, X: Xsub}) + 1)
s_a = sp.expand(s_expr).coeff(ee, 1)
s_b = sp.expand(s_expr).coeff(ee, 2)
OUT["L6"]["s_a_exact_match_2-2ln2"] = bool(sp.simplify(s_a - (2 - 2*L2sym)) == 0)
OUT["L6"]["s_b_exact_match"] = bool(sp.simplify(s_b + (Q(5,2) + 2*L2sym**2 - 4*L2sym)) == 0)
OUT["L6"]["s_b_resolved"] = str(sp.expand(s_b))

json.dump(OUT, open('/home/user/promote_stage2.json','w'), indent=1, default=str)
print("ALL DONE %.1fs" % (time.time()-t0))
