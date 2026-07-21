"""
NOTE 13, stage C: the law extraction + final audits.

Exact series (anchors d=9..121 note-12 + 131..201 tonight, all measured by the
same envelope machinery) vs the two-quadrature law

    m(d) = m* + [A0 - B0*(ln d - ln ln d)]/d + O((ln d)^2/d^2),
    A0 = INT_{-1}^{inf} phi(s) phi(-s-2) AL(s) ds     (AL = (s+5)(1-ln((s+5)/2)))
    B0 = INT_{-inf}^{-1} phi(s)^2 (|s|-1) ds          (= 0.0185094013627119)

Tasks:
 C1: quadratures A0, B0, and the second-order cross coefficient
     K2 = INT phi(s)^2 * (|s|/(2 sig^2)) * (|s|-1)^2 ds  (predicts the model's
     own residual bias ~ K2*(ln d)^2/d^2; fitted bias from exact-minus-model).
 C2: fit exact series d >= 61:  delta_m * d = A_eff - B_fit L(d); report
     B_fit vs B0 (predict within 5%); A_eff vs A0 (absorbs log-constants).
 C3: self-consistency: d* predicted by the fitted law vs model sweep 123.6.
 C4: d0 (m = m* crossing) from the law vs exact data (between 45 and 61).
 C5: assemble the verdict table for the note.
"""
import json, math
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 60
sig = mp.mpf("1.5")
phi = lambda v: mp.e**(-(v**2)/(2*sig**2))/(sig*mp.sqrt(2*mp.pi))
PHI = lambda v: mp.mpf("0.5")*(1 + mp.erf(v/(sig*mp.sqrt(2))))
mstar = PHI(-1)**2/2 + mp.quad(lambda u: phi(u)*PHI(-u-2), [-1, mp.inf])
mstar_f = float(mstar)

AL = lambda u: (u+5)*(1 - mp.log((u+5)/2))
A0 = mp.quad(lambda u: phi(u)*phi(-u-2)*AL(u), [-1, mp.inf])
B0 = mp.quad(lambda u: phi(u)**2*(mp.fabs(u)-1), [-mp.inf, -1])
K2 = mp.quad(lambda u: phi(u)**2*(mp.fabs(u)/(2*sig**2))*(mp.fabs(u)-1)**2, [-mp.inf, -1])
print(f"(C1) A0 = {mp.nstr(A0, 15)}")
print(f"(C1) B0 = {mp.nstr(B0, 15)}")
print(f"(C1) K2 = {mp.nstr(K2, 15)}")

a2 = json.load(open("ring_stageA2.json"))
sweep = {int(k): v for k, v in a2["sweep"].items()}
sb = json.load(open("ring_stageB.json"))
exact = {int(k): v["mass"] for k, v in sb.items()}
anchors = {9: 0.08773584175648119, 11: 0.08755547157788661, 13: 0.08736426802531863,
           15: 0.08719990720393687, 17: 0.08706461572198995, 19: 0.08695451490949664,
           23: 0.086790154156153, 29: 0.08663354848051617, 35: 0.08653775409398155,
           45: 0.08644568837824819, 61: 0.08637768869767828, 91: 0.08633686882830088,
           121: 0.08632885968162476}
exact_full = dict(anchors); exact_full.update(exact)

# fitted model-bias (exact - model) vs K2 (ln d)^2 / d^2
print("\n(C1b) model-bias  /=~  K2 (ln d)^2/d^2 :")
for d in (45, 61, 91, 121, 131, 141, 151, 161, 171, 181, 201):
    if d in sweep or d in (45, 61, 91, 121):
        if d in sweep: mod = sweep[d]
        else:
            # model values recomputed in stage A2 re-audit
            mod = {45: None, 61: None, 91: None, 121: None}[d]
val_diffs = {45: 1.094e-6, 61: 8.68e-7, 91: 5.49e-7, 121: 3.78e-7,
             131: exact[131]-sweep[131], 141: exact[141]-sweep[141],
             151: exact[151]-sweep[151], 161: exact[161]-sweep[161],
             171: exact[171]-sweep[171], 181: exact[181]-sweep[181],
             201: exact[201]-sweep[201]}
# NOTE: diffs at 45..121 from re-audit stdout (prob units).
predK = float(K2)
print("   d |    bias(prob)     K2(ln d)^2/d^2   ratio")
for d, b in val_diffs.items():
    pk = predK*(math.log(d))**2/d**2
    print(f"  {d:3d} | {b:+.4e}     {pk:.4e}      {b/pk:+.3f}")

# C2: fit exact series d >= 61
print("\n(C2) fit:  delta_m * d = A_eff - B_fit * L(d),  L = ln d - ln ln d")
ds = np.array(sorted(d for d in exact_full if d >= 61), float)
dd = np.array([exact_full[int(v)] for v in ds]) - mstar_f
L  = np.log(ds) - np.log(np.log(ds))
Y  = dd*ds                       # delta_m * d   (prob units)
cf = np.polyfit(L, Y, 1)
Bfit, Aeff = -cf[0], cf[1]
resid = np.polyval(cf, L) - Y
print(f"  A_eff = {Aeff:.6e}   B_fit = {Bfit:.6e}   rms resid = {np.sqrt((resid**2).mean()):.3e}")
print(f"  vs first-principles:  A0 = {float(A0):.6e}   B0 = {float(B0):.6e}")
print(f"  B-fit/B0 = {Bfit/float(B0):.4f}   [predict within 5%]")
for v in ds:
    pred = (Aeff - Bfit*(math.log(v)-math.log(math.log(v))))/v + mstar_f
    print(f"    d={int(v):3d}: exact {100*exact_full[int(v)]:.7f}  law {100*pred:.7f}  d_m-d_law {1e6*(exact_full[int(v)]-pred):+.3f}e-6")

# C3: min from the law
print("\n(C3) minimum self-consistency")
dss = np.linspace(100, 200, 40001)
mm = (Aeff - Bfit*(np.log(dss)-np.log(np.log(dss))))/dss
imin = int(np.argmin(mm)); dlaw = dss[imin]
print(f"  law minimum at d = {dlaw:.2f}   (model sweep: 123.60; predict within 6)")

# C4: zero crossing d0
Lc = Aeff/Bfit
d0 = math.exp(Lc + math.log(Lc))   # ln d0 - ln ln d0 = Lc  =>  d0 = Lc * e^Lc approx, iterate
for _ in range(30): d0 = math.exp(Lc)*math.log(d0)
print(f"(C4) m = m* crossing: law d0 = {d0:.1f}   exact data: m(45) > m*, m(61) < m*  => in (45,61)")

json.dump({"A0": str(A0), "B0": str(B0), "K2": str(K2),
           "fit": {"Aeff": Aeff, "Bfit": Bfit, "rms": float(np.sqrt((resid**2).mean()))},
           "law_min_d": float(dlaw), "crossing_d0": float(d0),
           "exact_full": {str(d): exact_full[d] for d in sorted(exact_full)},
           "model_bias": {str(d): val_diffs[d] for d in val_diffs}},
          open("ring_stageC.json", "w"), indent=1, default=str)
print("\nsaved ring_stageC.json")
