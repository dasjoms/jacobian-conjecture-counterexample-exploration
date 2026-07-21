#!/usr/bin/env python3
# LAB NOTE 16 figure: four panels of the promotion.
import json, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fractions import Fraction as Fr

S3 = json.load(open('/home/user/promote_stage3.json'))
S2 = json.load(open('/home/user/promote_stage2.json'))
JF = json.load(open('/home/user/jcorner_final.json'))

ghost = {int(k): float(v['root']) for k, v in S3['T3']['ghost_roots'].items()}
left  = {int(k): float(v) for k, v in S3['T3']['left_roots'].items()}

fig, ax = plt.subplots(2, 2, figsize=(12.6, 9.6))
fig.suptitle("NOTE 16 — THE PROMOTION: five laws become theorems", fontsize=15, fontweight='bold')

# ---------- A: the dance panorama ----------
a = ax[0,0]
for d, r in ghost.items():
    a.plot(r, d, 'o', ms=6, color='#0a8a3a', zorder=5)
for d, r in left.items():
    a.plot(r, d, 's', ms=6, color='#c03030', zorder=5)
for d in range(2, 22):
    if d not in left: a.plot([], [])
a.axvline(1/3, color='#0a8a3a', ls='--', lw=1.2, alpha=0.7)
a.axvline(-1, color='#c03030', ls='--', lw=1.2, alpha=0.7)
a.axvline(0.5, color='gray', ls=':', lw=1)
a.axvspan(0, 0.5, color='#0a8a3a', alpha=0.06)
a.axvspan(-1, 0, color='#c03030', alpha=0.05)
a.text(1/3+0.012, 20.2, "t=1/3", color='#0a8a3a', fontsize=9)
a.text(-1.04, 2.0, "t=-1", color='#c03030', fontsize=9, rotation=90, va='bottom')
a.set_xlim(-1.7, 0.62); a.set_ylim(1.2, 21.8)
a.set_xlabel("t"); a.set_ylabel("d")
a.set_title("(A) the reality dance, pinned: exactly one root in (0,1/2),\none extra in (-1,0) for odd d  [THEOREM, all d]", fontsize=10)
a.plot([], [], 'o', color='#0a8a3a', label='ghost cusp (all d)')
a.plot([], [], 's', color='#c03030', label='left cusp (odd d)')
a.legend(loc='lower left', fontsize=8)

# ---------- B: contraction constants ----------
b = ax[0,1]
def ghost_L(d):
    rho = Fr(1, d*d); c = Fr(6, d*(d+1)); den = 6 - 2*c
    Ap = (Fr(1,3)+rho)**(d-3)
    Ep = Ap * (Fr(2*d*d,3) - Fr(8*d,3) + 2 + Fr(d-1,d))
    return float(Ep/den)
def left_L(d):
    u0 = math.log((2*d-1)/8)/(d-2)
    u_lo = u0*0.85; u_hi = u0*1.15
    c = Fr(6, d*(d+1)); ul = Fr(u_lo); uh = Fr(u_hi)
    N  = lambda u: 8 - 6*u - c*(3 - 2*u)
    Dn = lambda u: (2*d-1) - d*u
    M  = lambda u: (-6+2*c)*Dn(u) + d*N(u)
    R_lo = N(uh)/Dn(uh)
    Rp = max(-M(ul), -M(uh)) / Dn(uh)**2
    return float(Rp/((d-2)*R_lo))
ds = list(range(5, 61))
b.semilogy(ds, [ghost_L(d) for d in ds], '.-', ms=5, lw=1, color='#2040c0', label='ghost:  L_g(d)')
dos = list(range(5, 60, 2))
b.semilogy(dos, [left_L(d) for d in dos], 's-', ms=4, lw=1, color='#a05000', label='left cusp:  L_L(d)')
b.axhline(1.0, color='red', lw=1.2, ls='--'); b.axhline(0.5, color='red', lw=0.8, ls=':')
b.axvspan(6.7, 60, ymin=0, ymax=1, color='#2040c0', alpha=0.04)
b.text(30, 1.25, "Banach zone: L < 1/2, d >= 7", fontsize=9, color='#2040c0')
b.text(6.0, 1.25, "L = 1", color='red', fontsize=8)
b.set_xlabel("d"); b.set_ylabel("proven Lipschitz constant")
b.set_title("(B) both fixed-point maps are contractions:\nexistence + uniqueness + lightning iteration", fontsize=10)
b.legend(fontsize=8, loc='upper right'); b.set_ylim(1e-7, 3)

# ---------- C: ghost residual vs proven envelope ----------
c = ax[1,0]
def proven_err(d):
    cc = Fr(6, d*(d+1)); rho = Fr(1, d*d)
    Ap = (Fr(1,3)+rho)**(d-3)
    return float(Ap*(Fr(2*d,3)-1+Fr(1,d))/(6-2*cc))
dA = list(range(7, 52))
c.semilogy(dA, [proven_err(d)*1.0 for d in dA], '--', color='black', lw=1.4, label='proven envelope  err(d)')
tab = JF['C2']['table']
da = sorted(int(k) for k in tab)
c.semilogy(da, [float(tab[str(d)][0]) for d in da], 'o', color='#0a8a3a', ms=5, label='actual  |t_g - closed form|  (d=21..47)')
# d=12..21 actual from stage 3
d12 = list(range(12, 22))
res = []
for d in d12:
    root = float(S3['T3']['ghost_roots'][str(d)]['root'])
    closed = 1/3 - 1/(3*(d*(d+1)-2))
    res.append(abs(root-closed))
c.semilogy(d12, res, 'o', color='#7ac0a0', ms=5, label='actual (d=12..21)')
c.semilogy(dA, [d*3.0**2/3.0**d for d in dA], ':', color='gray', lw=1, label='9d * 3^-d  reference')
c.set_xlabel("d"); c.set_ylabel("error")
c.set_title("(C) the ghost's closed form:  1/3 - t_g = 1/[3(d(d+1)-2)] + eps\ntheorem envelope always above reality (audit: dominates at every d)", fontsize=10)
c.legend(fontsize=8, loc='lower left')
c.set_ylim(1e-24, 1e-1)

# ---------- D: parity dichotomy ----------
d4 = ax[1,1]
Us = np.linspace(0.0, 0.36, 1201)
def Fn_odd(U, d):
    X = (1+U)**(d-2)
    return ((d*d+d)*(U+2)**2*(X-2) + 4*U*U + 17*U + 19 - X*(U+1)*(1+d*(U+2))) / (d*(d+1)*(U+2)**2)
def Fn_even(U, d):
    X = (1+U)**(d-2)
    return (4*U*U + 17*U + 19 - 2*(d*d+d)*(U+2)**2 + X*((U+1) - d*(U+2) - d*d*(U+2)**2)) / (d*(d+1)*(U+2)**2)
cols = {5:'#003399', 7:'#2255cc', 9:'#5588ee', 11:'#88aaff'}
for d, col in cols.items():
    d4.plot(Us, [Fn_odd(u, d) for u in Us], color=col, lw=1.6, label=f'odd d={d}')
cole = {4:'#8a0050', 6:'#b03070', 8:'#d06090'}
for d, col in cole.items():
    d4.plot(Us, [Fn_even(u, d) for u in Us], '--', color=col, lw=1.6, label=f'even d={d}')
d4.axhline(0, color='black', lw=0.8)
roots = {5: 0.27116049, 7: 0.15595096, 9: 0.10890925, 11: 0.08345862}
for d, r in roots.items():
    d4.plot(r, 0, 'o', color=cols[d], ms=7, zorder=6, mec='black')
d4.text(0.242, -2.6, "even: F = -XA - B < 0\nno corner shadow (THEOREM)", fontsize=9, color='#8a0050')
d4.text(0.13, 2.2, "odd: one crossing each,\n U(d) -> ln2/(d-2) (THEOREM + series)", fontsize=9, color='#003399')
d4.set_xlabel("U = |t*+1|"); d4.set_ylabel("F / [d(d+1)(U+2)^2]")
d4.set_title("(D) the parity split of the corner shadow", fontsize=10)
d4.set_ylim(-8, 6); d4.legend(fontsize=7, loc='lower left', ncol=2)

fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig('/home/user/promote_figure.png', dpi=150)
print("figure saved")
