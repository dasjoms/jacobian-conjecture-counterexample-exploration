"""
NOTE 15, stage 1 (rerun, instrument build): THE SINGULARITY CENSUS.
===================================================================
Fast exact cusp machinery: np.roots prefilter + 100-digit mpmath Newton polish
on the EXACT rational p_d', audits: residual + Sturm cross-count.
No predictions locked here (instrument); the honesty-ledger records that the
first stage-1 attempt (killed at 1500s by CRootOf isolation) falsified two of
my pre-committed locks BEFORE completion:
  * (L3) "cusp1 race: d*(1-t1) -> 5/2"  -- FALSIFIED by d<=15 data: the
    positive real cusp does not migrate to 1, it is the GHOST cusp -> 1/3.
    (I had conflated cusp parameters with bitangent contact parameters.)
  * (L4a) "u*d/ln d in (0.8,1.5) for all odd d=11..41" -- FALSIFIED at d=11
    (value ~0.49); the correct numerator is ln((2d-1)/8), refined in stage 3.
Confirmed before the kill: dance 2,1 pattern d<=20; t2(11) = -0.8940298 in
note-12's window (-0.8950,-0.8925).
===================================================================
"""
import json, time
import numpy as np
import mpmath as mp
import sympy as sp
from sympy import symbols, expand, diff, integrate, Rational as R, Poly

mp.mp.dps = 100
w = symbols("w")
t0 = time.time()
out = {"honesty": {"L3_falsified": "cusp->1 race was a contact-cusp confusion",
                   "L4a_falsified": "u*d/lnd ~ 0.49 at d=11, not in (0.8,1.5)"}}

def seed(d):
    return expand(2*w - 3*w**2 + w*(1-w)*(w**sp.Integer(d-2) - R(6, d*(d+1))))

def sturm_real_count(f, x):
    seq = [expand(f), expand(diff(f, x))]
    while True:
        rr = sp.rem(seq[-2], seq[-1], x)
        if rr == 0 or getattr(rr, "is_zero", False):
            break
        seq.append(expand(-rr))
        if seq[-1].is_number:
            break
    def sig(inf_sign):
        L = []
        for g in seq:
            P = sp.Poly(g, x); lc, dg = sp.LC(P), sp.degree(P, x)
            L.append(sp.sign(lc)*(inf_sign if dg % 2 else 1) if dg > 0 else sp.sign(g))
        return L
    def var(L):
        L = [v for v in L if v != 0]
        return sum(1 for a, b in zip(L, L[1:]) if a*b < 0)
    return var(sig(-1)) - var(sig(1))

# ---------------- dance d = 2..21 + spots ----------------
print("="*84); print("reality dance (Sturm exact) d=2..21, spots {23,29,37,41}"); print("-"*84)
dance = {}
for d in range(2, 22):
    cnt = sturm_real_count(diff(seed(d), w), w)
    dance[d] = int(cnt); print(f"  d={d:2d}: {cnt}", flush=True)
spots = {}
for d in (23, 29, 37, 41):
    cnt = sturm_real_count(diff(seed(d), w), w)
    spots[d] = int(cnt); print(f"  spot d={d}: {cnt}", flush=True)
dance["spots"] = spots
ok_dance = [dance[d] for d in range(11, 21)] == [2,1,2,1,2,1,2,1,2,1] and dance[21] == 2
ok_spots = spots == {23: 2, 29: 2, 37: 2, 41: 2}
print(f"  odd-d count 2 / even-d count 1 through d=21: {ok_dance}; spots all 2: {ok_spots}")
out["dance"] = dance; out["dance_ok"] = bool(ok_dance and ok_spots)

# ---------------- fast exact roots: np prefilter + Newmap polish, audits ----------------
def real_roots_exact(d):
    """real roots of p_d' at 100 digits, audited."""
    p = seed(d); pp = diff(p, w)
    P = Poly(sp.expand(pp * (d*(d+1))), w)     # integer coeffs
    coeffs = [mp.mpf(int(c)) for c in P.all_coeffs()]
    fc = [float(c) for c in P.all_coeffs()]
    cand = np.roots(fc)
    ppf = sp.lambdify(w, pp, "mpmath")
    pppf = sp.lambdify(w, diff(pp, w), "mpmath")
    nrm = sum(abs(c) for c in coeffs)
    roots = []
    for z in cand:
        if abs(z.imag) < 1e-7:
            xr = mp.mpf(z.real)
            try:
                xr = mp.findroot(ppf, xr, tol=mp.mpf("1e-90"), verify=False)
            except Exception:
                pass
            resid = abs(ppf(xr))
            if resid < nrm * mp.mpf("1e-60"):
                if all(abs(xr - q) > 1e-30 for q in roots):
                    roots.append(xr)
    roots.sort()
    aud = max((abs(ppf(r_)) for r_ in roots), default=mp.mpf(0))
    ok = aud < nrm * mp.mpf("1e-60")
    return roots, aud, ok

roots_tab, races, ghosts = {}, {}, {}
print("="*84); print("real cusp parameters d=2..41 (audited, 60-digit strings stored)"); print("-"*84)
bad = []
for d in range(2, 42):
    roots, aud, ok = real_roots_exact(d)
    if not ok: bad.append(d)
    roots_tab[d] = [mp.nstr(r_, 60) for r_ in roots]
    p = seed(d); Phi = integrate(p, w)
    ta = sp.lambdify(w, expand(w*p - Phi), "mpmath")
    if d % 2 == 1 and d >= 3:
        t1 = max(roots, key=abs)          # ghost near +1/3 (largest |.| of the pair? see below)
        # robust species: ghost = root nearest 1/3; left = the other
        t1 = min(roots, key=lambda r_: abs(r_ - mp.mpf(1)/3))
        t2 = [r_ for r_ in roots if r_ != t1]
        t2 = t2[0] if t2 else None
        races[d] = {"t_ghost": mp.nstr(t1, 60), "t_left": mp.nstr(t2, 60),
                    "u": mp.nstr(t2+1, 60), "s_ghost": mp.nstr(sp.lambdify(w, p, "mpmath")(t1), 50),
                    "tau_ghost": mp.nstr(ta(t1), 50), "s_left": mp.nstr(sp.lambdify(w, p, "mpmath")(t2), 50),
                    "tau_left": mp.nstr(ta(t2), 50)}
    else:
        tv = min(roots, key=lambda r_: abs(r_ - mp.mpf(1)/3))
        ghosts[d] = {"t_ghost": mp.nstr(tv, 60),
                     "s_ghost": mp.nstr(sp.lambdify(w, p, "mpmath")(tv), 50),
                     "tau_ghost": mp.nstr(ta(tv), 50)}
    if d <= 11 or d % 4 == 1:
        print(f"  d={d:2d}: {[mp.nstr(r_,10) for r_ in roots]}  audit {mp.nstr(aud,3)}", flush=True)
# Sturm cross-count audit (d<=21 + spots)
cross_ok = True
for d in list(range(2, 22)) + [23, 29, 37, 41]:
    cnt = dance[d] if d <= 21 else dance["spots"][d]
    if len(roots_tab[d]) != cnt:
        cross_ok = False
print(f"  cross-audit root-count == Sturm everywhere checked: {cross_ok};  bad-audit d: {bad}")
out["roots"] = roots_tab; out["races"] = races; out["ghosts"] = ghosts
out["cross_ok"] = bool(cross_ok); out["root_audits_ok"] = (not bad)

print("="*84); print("ghost race (all d): t_ghost -> 1/3, image -> (1/3, 1/27)")
for d in range(2, 42):
    tg = mp.mpf((races if d % 2 else ghosts)[d]["t_ghost"])
    print(f"  d={d:2d}: t_ghost={mp.nstr(tg,12)}  d^2*(1/3 - t)={mp.nstr(d*d*(mp.mpf(1)/3-tg),8)}"
          + (f"   u_left={races[d]['u'][:12]}  u(d)*(d-2)/ln((2d-1)/8)=" +
             mp.nstr((mp.mpf(races[d]["u"])*(d-2))/mp.log((2*d-1)/8), 6) if d % 2 else ""), flush=True)

print("="*84); print("diagonal crossings Phi = (t-1)p, d = 2..12")
diag = {}
for d in range(2, 13):
    p = seed(d); Phi = integrate(p, w)
    G = Poly(sp.expand((Phi - (w-1)*p) * (d*(d+1))), w)
    fc = [float(c) for c in G.all_coeffs()]
    pf = sp.lambdify(w, expand(Phi - (w-1)*p), "mpmath")
    pts = []
    for z in np.roots(fc):
        if abs(z.imag) < 1e-8:
            xr = mp.mpf(z.real)
            try: xr = mp.findroot(pf, xr, tol=mp.mpf("1e-85"), verify=False)
            except Exception: pass
            if abs(pf(xr)) < mp.mpf("1e-50") and all(abs(xr-q[0]) > 1e-25 for q in pts):
                spoint = sp.lambdify(w, p, "mpmath")(xr)
                pts.append((xr, spoint))
    diag[d] = [{"t": mp.nstr(a, 50), "s": mp.nstr(b, 50)} for a, b in sorted(pts)]
    print(f"  d={d:2d}: " + "  ".join(f"t={mp.nstr(a,10)} s={mp.nstr(b,10)}" for a, b in sorted(pts)), flush=True)
out["diag"] = diag

# t=1 certificate: (p_d(1), tau_d(1)) == (-1,-1) for all d = 2..41
allcert = True
for d in range(2, 42):
    p = seed(d); Phi = integrate(p, w)
    allcert &= (p.subs(w, 1) == -1) and (sp.expand(w*p - Phi).subs(w, 1) == -1)
print("="*84); print(f"CERT: (-1,-1) = (p_d(1), tau_d(1)) is a wall point for EVERY d = 2..41: {allcert}")
out["corner_fixed_point_cert"] = bool(allcert)

json.dump(out, open("jcorner_stage1.json", "w"), indent=1)
print(f"saved jcorner_stage1.json   [{time.time()-t0:.0f}s]")
