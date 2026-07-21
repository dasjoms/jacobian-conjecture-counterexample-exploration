"""Assemble NOTE 13 (The Ringing Constant) as one markdown file with the figure embedded."""
import json, base64

fig = base64.b64encode(open("ring_figure.png", "rb").read()).decode()
a2 = json.load(open("ring_stageA2.json")); sb = json.load(open("ring_stageB.json"))
sc = json.load(open("ring_stageC.json"))
fit = sc["fit"]

def L(d):
    import math
    return math.log(d) - math.log(math.log(d))
mstar = 0.08644464653779008
Aeff, Bfit = fit["Aeff"], fit["Bfit"]
m251 = 100*(mstar + (Aeff - Bfit*L(251))/251)

rowsB = "\n".join(
    f"| {d} | {sb[str(d)]['mass_pct']:.7f} | [{sb[str(d)]['lock_window'][0]:.6f}, {sb[str(d)]['lock_window'][1]:.6f}] "
    f"| {'**yes**' if sb[str(d)]['in_window'] else 'NO'} | {sb[str(d)]['drops']} / {sb[str(d)]['kept']} |"
    for d in (131, 141, 151, 161, 171, 181, 201))

exact_full = {int(k): v for k, v in sc["exact_full"].items()}
rowsFull = "\n".join(
    f"| {d} | {100*exact_full[d]:.7f} | {100*(exact_full[d]-mstar):+.7f} |"
    for d in sorted(exact_full))

sweep = {int(k): v for k, v in a2["sweep"].items()}
rowsSweep = "\n".join(f"| {d} | {100*sweep[d]:.7f} |" for d in sorted(sweep))

md = f"""# JACOBIAN LAB — Note 13: **The Ringing Constant** 🚂🔔
### The tower overshoots its limit, undershoots it, then comes home — and the home address is $[A - B(\\ln d-\\ln\\ln d)]/d$ with $B$ predicted *a priori*

*Date: 2026-07-20/21 (last shift of the night). Instruments: sympy 1.14 (symbolic layer series with
$\\varepsilon = 1/d$ and solved perturbations substituted), mpmath 50–100 digits (Newton, bisection, quadrature),
numpy. All locks were frozen in `jacobian_ring_2b.py` BEFORE the exact masses were measured
(`jacobian_ring_3.py`); data in `ring_stage{{A2,B,C}}.json`; figure `ring_figure.png`; scripts
`jacobian_ring_{{1,1b,2,2b,3,4,figure}}.py`.*

---

## 0. The question note 12 left hanging

Note 12 proved the missed-cone mass of the tower's even-fiber chambers converges to an exact constant

  $$m^* = \\Phi_\\sigma(-1)^2/2 + \\int_{{-1}}^{{\\infty}}\\!\\phi_\\sigma(s)\\,\\Phi_\\sigma(-s-2)\\,ds
        = 0.0864446465377901\\ldots \\qquad(\\sigma = 1.5)$$

— and ended with a locked falsifiable prediction: *"the tail must rebound: m(161) > m(121) = 8.6329,
both tending to m\\* from below."* The measured series was monotone down through d = 121 and sat
**below** m\\* — 8.63289 < 8.64446. For the theory to be complete, the descent had to turn around.

Tonight: it turns around, on exact data, precisely where the two-scale theory said it would.

| d | 121 | 131 | 141 | 151 | **161** | 171 | 181 | 201 |
|---|---|---|---|---|---|---|---|---|
| exact m(d) [%] | 8.6328860 | 8.6328505 | 8.6328959 | 8.6329839 | **8.6331010** | 8.6332363 | 8.6333759 | 8.6336791 |

$$m(161) - m(121) = +2.150\\times10^{{-6}}\\ \\text{{(prob units)}} > 0\\quad\\checkmark$$

---

## 1. The two scales (the whole story in one picture, panel (c))

The missed mass is $m(d) = \\int\\phi_\\sigma(s)\\,\\Phi_\\sigma(\\mathrm{{env}}_d(s))\\,ds$ where
$\\mathrm{{env}}_d(s) = \\min_t\\tau_d(t)$ over real roots of $p_d(t) = s$; the limit envelope is
$\\min(s,\\,-s-2)$ with the corner at the pin $(-1,-1)$. Finite $d$ perturbs the two binding branches
**with opposite signs and different decay rates**:

**LEFT pin layer** (t = −1 − x/d, odd d): the symbolic ε-series with the solved perturbations
$x = x_0 + x_1\\varepsilon + x_2\\varepsilon^2$, $x_0 = \\ln\\frac{{s+5}}2$, substituted back *by the machine*, gives

$$\\boxed{{\\ \\tau_L(s; d) = -s-2 + \\frac{{A_L(s)}}{{d}} + \\frac{{B_L(s)}}{{d^2}} + O(d^{{-3}}),\\qquad
A_L(s) = (s+5)\\Big(1 - \\ln\\tfrac{{s+5}}2\\Big)\\ }}$$

with $B_L$ an explicit log-quadratic (machine string in `ring_probe_sym.json`). Certificates:
$\\tau_0 = -s-2$ exact; $A_L$ exact; 100-digit root audits at d = 21, 61, 121, 201 track
$A_L + B_L/d$ to $10^{{-4}}$ at d = 201. **$A_L \\ge 0$ on the whole weight window**, with equality
only at s = 2e − 5 — and then it goes *slightly negative*: near $s = 0.437$, the finite-d branch
crosses its own asymptote (red dotted marker in panel (c)); Euler's number fell out of the wall.

**RIGHT pin layer** (t = 1 + c/d): here c = O(1) fails — c *grows with d*. The balance
$p_d(1+c/d) = s$ closes as

$$c\\,(4 + e^c) = d\\,|s+1|\\quad\\Longrightarrow\\quad c = \\ln d - \\ln\\ln d + \\ln|s+1| + o(1),$$

$$\\tau_R(s; d) = s + \\frac{{c\\,(s+1)}}{{d}} + \\frac{{(c-1)e^c}}{{d^2}} + \\cdots\\ <\\ s$$

(the branch pokes *below* the limit line r = s; exact-vs-model audits: 2%–4% at d = 21 shrinking
to 0.2%–3% at d = 201; the Lambert-W iterates converge to the measured branch).

So the ring is a **fight between a fast positive term and a slow negative term**:

$$m(d) \\approx m^* + \\underbrace{{\\frac{{A_0}}{{d}}}}_{{>0,\\ \\text{{left}}}}\\ -\\ \\underbrace{{\\frac{{B\\,(\\ln d - \\ln\\ln d)}}{{d}}}}_{{>0,\\ \\text{{right, slower}}}}$$

which is *exactly* the anatomy of an overshoot → zero-crossing → undershoot → rebound.

---

## 2. The GOLD model and its interrogation

Replace the branches by the certified expansions (left) and an exact 50-digit solve of the right
branch per quadrature node (right), take the min — the result is a measurement-free predictor
`m_mod(d)`:

| d | exact | Gold model | diff [pp] |
|---|---|---|---|
| 9 | 8.7735842 | 8.8295098 | −0.05593 |
| 15 | 8.7199907 | 8.7278876 | −0.00790 |
| 23 | 8.6790154 | 8.6791676 | −0.00015 |
| 29 | 8.6633548 | 8.6624648 | +0.00089 |
| 45 | 8.6445688 | 8.6434749 | +0.00109 |
| 61 | 8.6377689 | 8.6369006 | +0.00087 |
| 91 | 8.6336869 | 8.6331384 | +0.00055 |
| 121 | 8.6328860 | 8.6325077 | +0.00038 |

The small-d drift is the dropped $O((\\ln d)^2/d^3)$ tail; the smooth +bias at d ≥ 29 **is itself
theory-shaped**: (exact − model) / [$(\\ln d)^2/d^2$] = 0.25–0.28, flat over d = 131…201 — the
second-order cross term $\\tfrac12\\Phi''\\,\\delta\\mathrm{{env}}^2$ (predicted coefficient family
$K_2 = 0.00925$, $K_2^{{L}} = 0.01416$; the observed 0.0024-coefficient's exact composition is
tonight's leftover, queued).

The same model, swept densely, placed the minimum at d\\* ≈ **123.6**, m ≈ 8.63251 — and every one
of the seven exact measurements then landed inside the pre-locked windows:

| d | exact m(d) [%] | locked window [%] | inside? | dropped roots |
|---|---|---|---|---|
{rowsB}

---

## 3. The law

Least squares on the ten exact points d ≥ 61, in coordinates $(L,\\ d(m-m^*))$,
$L(d) = \\ln d - \\ln\\ln d$:

$$\\boxed{{\\ m(d) = m^* + \\frac{{0.046546 - 0.018764\\,(\\ln d - \\ln\\ln d)}}{{d}},\\qquad
\\text{{rms residual }}1.2\\times10^{{-5}}\\ }}$$

Against **first principles** (pure quadrature, no fit):

$$B_0 = \\int_{{-\\infty}}^{{-1}}\\phi_\\sigma(s)^2(|s|-1)\\,ds = 0.0185094013627119
\\quad\\Longrightarrow\\quad B_{{fit}}/B_0 = 1.0138\\ \\ (\\textbf{{1.4%}})$$

$$A_0 = \\int_{{-1}}^{{\\infty}}\\phi_\\sigma(s)\\,\\phi_\\sigma(-s-2)\\,A_L(s)\\,ds = 0.0311269250795446,
\\qquad A_0 + B_{{const}} = 0.03531$$

with $B_{{const}} = \\int\\phi^2(s+1)\\ln(|s|-1)\\,ds = 0.0041816$ the right layer's constant-in-d
piece. The effective intercept $A_{{eff}} = 0.0465$ reabsorbs the remaining slow constants along
the fit window (intercepts are not invariant under sub-leading rearrangement; the gap
$A_{{eff}} - A_0 - B_{{const}} = 0.0112$ decomposes into the $(1-1/c)|s|/d$ family and is queued).

Consequences of the law (all checked against exact data):

* **zero-crossing** $d_0 = 45.7$ — exact data bracket: m(45) > m\\*, m(61) < m\\* ✓;
* **minimum** $d^* = 128.3$ (fit), depth $m^* - m(d^*) \\approx 0.0119$pp — model sweep says 123.6;
  exact min-bracket d ∈ [121, 131]: m(121) = 8.63289, m(131) = 8.63285 both consistent ✓;
* **rebound slopes**: m(171)−m(161) = +1.353e-6, m(201)−m(161) = +5.781e-6 (locks: > −1.5e-6) ✓;
* **next far lock**: m(251) = **{m251:.6f}%** per the law — frozen here for a future round.

The full 20-chamber exact series (nine anchors from note 12 + seven tonight + four re-audited):

| d | m(d) [%] | m − m* [pp] |
|---|---|---|
{rowsFull}

---

## 4. The figure

![The Ringing Constant](data:image/png;base64,{fig})

**(a)** the ring itself; **(b)** the log-correction straightened, fit vs first principles (two nearly
parallel lines — the scientific content of this note); **(c)** the two scales made visible at
d = 201: measured $d\\cdot\\delta\\mathrm{{env}}$ tracks the gold right-layer curve for s < −1 and the
green certified $A_L + B_L/d$ for s > −1, switching at the corner; **(d)** scoreboard.

---

## 5. HONESTY LEDGER

1. **The phantom certificate.** First A_L run *printed False* — correctly! — because I had certified
   the series with the perturbation x₁, x₂ still symbolic (never substituted the solution). The
   audits caught it the same hour: s = 0 matched by luck (A_L(0) agrees for both wrong and right
   formulas), s = ±0.5/0.3 did not. New house rule, in ink: **an expansion is only as certified as
   its last substitution.**
2. **My own hand algebra lost a factor.** Hand-derived $A_L = (s+5) - 5\\ln((s+5)/2)$; machine-correct
   $(s+5)(1 - \\ln((s+5)/2))$. The cramped x₁-cancellation in my head was exactly wrong.
3. **λ-sign flip.** In the figure's $c_of$ I fed the *signed* s+1 instead of its modulus — a curve
   blew up to 290 on the draft render, which is how it got caught. (The mass model was immune: it
   only used $c$ as a Newton start and audited residuals.)
4. **Solver pathology tour.** mp.findroot stall at λ→0 (linear branch added); mp.findroot-bisect
   tolerance blowup in a wide bracket (manual 300-iteration bisect); np.roots+polish acceptance
   upgraded with an explicit $|$resid$|/$norm test — 0 drops in 28 585 polished roots.
5. **The amber.** The "≤ 0.0008pp at d ≥ 45" tightness lock missed by 0.0003pp — the second-order
   $(\\ln d)^2/d^2$ model-bias *is* the miss, with the right shape (ratio flat 0.25–0.28). Logged
   amber, not red: a measurement-free model good to a few parts in $10^6$ of its own target.
6. Five model-bias numbers (d = 45..121) were transcribed from run-1 stdout into the final script
   after a crash — re-computed identically on the re-audit run (0.001094 vs 0.001090 etc.); noted.

---

## 6. Scoreboard

| # | Prediction (locked before measurement) | Verdict |
|---|---|---|
| 1 | seven exact masses in locked windows (d = 131…201) | 🟢 7/7 inside |
| 2 | **m(161) > m(121)** (note-12's falsifiable lock) | 🟢 +2.150e-6 |
| 3 | strong: m(161) > m(121) − 3e-6 | 🟢 |
| 4 | slopes: m(171), m(201) > m(161) − 1.5e-6 | 🟢 +1.353e-6, +5.781e-6 |
| 5 | layer law $\\tau_L = -s-2 + A_L/d + B_L/d^2$ certified + 100-digit audits | 🟢 ≤ 1e-4 @201 |
| 6 | $A_L = (s+5)(1-\\ln((s+5)/2))$; zero at s = 2e − 5 = 0.436563657 | 🟢 exact |
| 7 | right-layer Lambert model, shrinking residuals | 🟢 0.2%–3% @201 vs 2%–4% @21 |
| 8 | $B_{{fit}}$ within 5% of first-principles $B_0$ | 🟢 1.38% |
| 9 | law min d\\* ∈ (119, 131) | 🟢 128.3 (fit) / 123.6 (sweep) |
| 10 | crossing d₀ ∈ (45, 61) | 🟢 45.7 |
| 11 | windows tightness d ∈ {{29, 35}} ≤ 0.004pp | 🟢 0.00115pp |
| 12 | windows tightness d ≥ 45 ≤ 0.0008pp | 🟡 0.00109pp = diagnosed 2nd-order bias |

**11 green, 1 amber.**

---

## 7. Queue (the train rolls on)

1. **Chamber n = 12 (d = 11)** — the full atlas we promised: wall deg 12, **terms(12) = 78 − 2 = 76**,
   K = den²(p₁₁′) → compute (queue's K(11) = 3025 is for d = 10/n = 11 — both chambers now fully
   locked: n = 11: terms 64, K = 3025, nodes 36 = 0 crunodes + 4 acnodes, census ∈ [8.63, 8.89],
   corner ∈ (−0.950, −0.937), cusp2 r(11) > 1.5526, S₁₁).
2. **Second-order closure**: compose the 0.0024 coefficient from $K_2$, $K_2^L$, the $c^2$-terms and
   the $\\Phi''$-crosses — the natural next note.
3. **Corner $\\sigma(d)$ law** vs note-12's $s^* = -0.9733$ fit: reconcile the finite-window fit with
   the tie-point analysis (model gives $s_c(d)$ cheaply; the crunode track converges to the pin).
4. General-n proof of $C_2, C_3$; content-chain law (needs the n ≥ 11 chambers); small-d K quirks.
5. Far lock frozen now: **m(251) = {m251:.6f}%**.
6. Swallowtail (5,1)-contact 4/5-slope; generic-seed deg-5 control atlas; Moh 2-D chamber.

*The tower dreamed of the cubic, overshot its dream, dipped below it, and tonight — on exact data,
inside locked windows — it turned around and started walking home, at exactly
$0.0188(\\ln d - \\ln\\ln d)/d$ per step. Good night, Jacobian lab. 🌙🚂*
"""

open("jacobian_ringing.md", "w").write(md)
print(f"saved jacobian_ringing.md  ({len(md)/1024:.0f} KB)")
