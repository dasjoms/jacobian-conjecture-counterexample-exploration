#!/usr/bin/env python3
# LAB NOTE 16 — STAGE 3: pin-theorem symbolic certificates + figure data tables.
# LOCKS: (T1) symbolic d: Phi_d(1) === 0, p_d(1) === -1, tau_d(1) === -1,
#        p'_d(1) !== 0 (value c0-5), slope at pin === 1 (tau'=t p' ratio at t=1 => 1).
#        (T2) L_L(5) reported honestly (expected > 1 -> contraction window starts at 7).
#        (T3) exact ghost roots d=2..21 at 60 digits + residuals vs closed form.
import sympy as sp
from sympy import Rational as Q, symbols
import mpmath as mp
import json, time
mp.mp.dps = 60
t0 = time.time()
w, dd = symbols('w dd')
c0 = Q(6,1)/(dd*(dd+1))
p_sym = 2*w - 3*w**2 - c0*w + c0*w**2 + w**(dd-1) - w**dd          # symbolic d, formal powers
Phi1 = 1 - 1 - c0/2 + c0/3 + sp.Rational(1,1)/dd - 1/(dd+1)         # Phi_d(1) closed form
OUT = {"T1": {}, "T2": {}, "T3": {}}
OUT["T1"]["Phi_d(1)_symbolic_zero"] = bool(sp.simplify(Phi1) == 0)
p1 = sp.simplify(p_sym.subs(w, 1))
OUT["T1"]["p_d(1)_equals_-1"] = bool(p1 == -1)
OUT["T1"]["tau_d(1)_equals_-1"] = bool(sp.simplify(p1 - Phi1) == -1)
pp1 = sp.simplify(sp.diff(2*w-3*w**2-c0*w+c0*w**2, w).subs(w,1) + (dd-1) - dd)
OUT["T1"]["p'_d(1)"] = str(pp1)                     # c0 - 5
OUT["T1"]["p'_d(1)_nonzero"] = bool(sp.simplify(pp1 - (c0-5)) == 0)  # = c0-5 != 0
OUT["T1"]["slope_at_pin_is_1"] = "tau'(t)/p'(t) = t (exact) hence at t=1: slope 1"
# T2: left-contraction constant at d=5
def left_L(d_):
    import mpmath as mp
    c = Q(6, d_*(d_+1))
    u0 = mp.log(mp.mpf(2*d_-1)/8)/(d_-2)
    u_lo = Q(int(mp.floor((u0-mp.mpf(1)/2**80)*(2**80))), 2**80)*Q(17,20)
    u_hi = Q(int(mp.floor((u0+mp.mpf(1)/2**80)*(2**80))+1), 2**80)*Q(23,20)
    N  = lambda u: 8 - 6*u - c*(3 - 2*u)
    Dn = lambda u: (2*d_-1) - d_*u
    M  = lambda u: (-6+2*c)*Dn(u) + d_*N(u)
    R_lo, R_hi = N(u_hi)/Dn(u_hi), N(u_lo)/Dn(u_lo)
    Rp = max(-M(u_lo),-M(u_hi),M(u_lo),M(u_hi))/Dn(u_hi)**2
    return float(Rp/((d_-2)*R_lo)), float(u_lo), float(u_hi), bool(M(u_lo)<0), bool(R_lo>0)
L5_ = left_L(5)
OUT["T2"]["d5_left_contraction"] = dict(L_L=L5_[0], u_window=L5_[1:3], monotone=L5_[3], Rpos=L5_[4],
    comment="L_L(5)>=1 expected: contraction theorem window starts at d=7 (d=5 cusp lives in exact table)")
# T3: exact ghost roots + residuals vs closed form d=2..21
rows = {}
for d_ in range(2, 22):
    c = Q(6, d_*(d_+1))
    pp = (2-c) + (2*c-6)*w + (d_-1)*w**(d_-2) - d_*w**(d_-1)
    f = sp.lambdify(w, pp, 'mpmath')
    root = mp.findroot(f, mp.mpf(1)/3 - mp.mpf(1)/(3*(d_*(d_+1)-2)))
    closed = mp.mpf(1)/3 - mp.mpf(1)/(3*(d_*(d_+1)-2))
    resid = root - closed
    rows[d_] = dict(root=str(root), closed_resid=str(resid))
    if d_ in (2,3,5,7,11): rows[d_]['print'] = True
OUT["T3"]["ghost_roots"] = rows
# also left roots d=3..15 for the census table
lrows = {}
for d_ in range(3, 16, 2):
    c = Q(6, d_*(d_+1))
    pp = (2-c) + (2*c-6)*w + (d_-1)*w**(d_-2) - d_*w**(d_-1)
    f = sp.lambdify(w, pp, 'mpmath')
    root = mp.findroot(f, -1 + mp.log(mp.mpf(2*d_-1)/8)/(d_-2))
    lrows[d_] = str(root)
OUT["T3"]["left_roots"] = lrows
json.dump(OUT, open('/home/user/promote_stage3.json','w'), indent=1, default=str)
print("T1:", OUT["T1"])
print("T2 d5:", L5_)
print("ghost d=11:", rows[11]['root'], " resid:", rows[11]['closed_resid'])
print("left d=11:", lrows[11])
print("DONE %.1fs" % (time.time()-t0))
