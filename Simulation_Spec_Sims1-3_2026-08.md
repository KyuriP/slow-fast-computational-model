# Simulation Specification — Revised Uncentered Model (Sims 1–3)

Code-ready spec for the first implementation pass. Feedback, shocks, and hysteresis (old
Simulation 4 / dynamic extension) are explicitly out of scope until 1–3 are stable, per plan.

## Core model (all three simulations)

    logit Pr(S_i = 1 | S_-i, P) = τ_i + Σ_{j≠i} ω_ij S_j + γ_i P

- `τ_i` — symptom-specific baseline activation tendency (threshold)
- `ω_ij` — symptom–symptom coupling (symmetric, `ω_ij = ω_ji`)
- `P` — slow contextual field (a person-level scalar, not time-indexed in Sims 1–3)
- `γ_i` — symptom-specific sensitivity to context

No β. No `J(m − 1/2)`. No shared threshold. Sampling is via the same exact-enumeration Gibbs/joint
sampler already implemented in `13_network_full_check.R` (`build_dgp()` / `sample_states()`),
adapted to take `τ_i`, `ω_ij`, `γ_i` directly instead of deriving them from the mean-field input —
this reuses existing, validated code rather than starting fresh.

## Design decisions to confirm before coding (flagging, not assuming)

1. **N = 12** (carries over from the existing appendix/scripts, keeps exact enumeration cheap:
   2^12 = 4,096 states). Alternative: N = 14 to match Cramer et al. (2016)'s symptom count exactly
   — stronger literature tie, ~4x the enumeration cost (2^14 = 16,384 states), still tractable.
   Recommend starting at N = 12 for the pilot and only moving to 14 if you want the figure to
   visually mirror Cramer's network one-to-one.
2. **τ_i calibration source.** Cramer et al. (2016) estimated real thresholds and weights for 14
   depression symptoms from the VATSPUD data (n = 8,973) via IsingFit (L1-regularized logistic
   regression, EBIC selection) — confirmed from the paper text, but the actual numeric table isn't
   in machine-readable form on the PLOS page (thresholds are shown as a figure — node fill in Fig.
   2a and a separate panel in Fig. 3, not a text table). Fastest path to the real numbers: ask
   Denny directly (he's a co-author and likely has the fitted values on hand), rather than trying
   to extract them from a rendered image. Until then, the τ_i range below (Uniform(−3.0, −1.0)) is
   a literature-informed placeholder reflecting the qualitative pattern the paper reports
   (low-base-rate symptoms like thoughts of death get much more negative thresholds than
   high-base-rate symptoms like fatigue) — swap in real values once obtained. This is a discrete,
   trackable action item, not a blocker for building/testing the code.
3. **"High burden" cutoff** for the probability-of-high-burden outcome in Sim 1/3: proposed
   `D ≥ N/2` (half or more symptoms active), reported alongside the full sumscore distribution so
   the cutoff choice doesn't hide anything. Flag if you'd rather use a DSM-style count (e.g., ≥5
   symptoms, mirroring MDD's "5 of 9" criterion, rescaled for N=12).

## Estimator and software (addresses Denny's "underdescribed" methods complaint directly)

Both `Simulation 2` and `Simulation 3` fit networks via **nodewise logistic regression**, base R
`stats::glm(family = binomial())`, one model per node regressing `S_i` on all other symptoms (plus
`P` in the context-adjusted version), symmetrized by averaging the two directed coefficient
estimates for each edge: `ω̂_ij = (β̂_{i←j} + β̂_{j←i}) / 2`. This is the same estimator already
implemented and validated in `13_network_full_check.R`'s `fit_edges()`. No new estimator or
package needs to be introduced — this should be stated explicitly in the Methods text, by name,
which is exactly what Denny flagged as missing.

## Outcome definitions (confounding framing, not estimation-bias framing)

- **Global strength:** `GS = Σ_{i<j} |ω̂_ij|`
- **Recovery error:** `MAE = C(N,2)^-1 Σ_{i<j} |ω̂_ij − ω_ij|` — read as *misspecification relative
  to the true causal coupling*, not estimator bias (per Denny's A/B/C distinction: this compares
  the causal/data-generating weights (A) against the estimated weights (C), which is a confounding
  question, not an estimation-error question in the statistical sense).

---

## Simulation 1 — Mean context shift under fixed coupling

**Question.** Can two groups show different symptom burden even when their symptom–symptom
coupling is identical?

**Data-generating model.** Core model above, with `P` a *fixed constant per group* (no
between-person variation within a group in this simulation — that's Simulation 2's job).

**Design.**
- Two groups: low-context (`P_g = P_low`) and high-context (`P_g = P_high`)
- Held fixed across groups: `τ_i`, `ω_ij`, `γ_i` (one data-generating network, shared)
- `n = 1,000–2,000` simulated individuals per group
- `100` replicate networks (redraw `τ_i`, `ω_ij`, `γ_i` each replicate) to average out
  network-specific idiosyncrasy, matching how Figure 7's numbers are already computed

**Outputs.**
- Mean symptom sumscore `D = Σ S_i`, by group
- Variance of `D`, by group
- `Pr(D ≥ N/2)`, by group ("high burden" rate)

**Figure.** Same network (drawn once, shown for reference), two symptom-sum distributions
side-by-side or overlaid (low- vs. high-context group) — replaces part of the old Figures 3–5,
now under the uncentered formulation and without requiring the fast→slow loop.

---

## Simulation 2 — Omitted context variation creates apparent connectivity

**Question.** When context varies across individuals and is omitted from the network model, do
estimated symptom–symptom edges become inflated even though the true coupling is unchanged?

**Data-generating model.** `P_n ~ Normal(μ_P, σ_P²)`, drawn independently per individual `n`; same
core model, with `τ_i`, `ω_ij`, `γ_i` fixed (one data-generating network, shared across the whole
sweep so `σ_P` is the only thing that changes between conditions).

**Design.**
- Sweep `σ_P` across a small grid (see pilot grid below) at fixed `μ_P`
- For each `σ_P`: fit symptom-only (`ω̂^omit`) and context-adjusted (`ω̂^adj`) models
- `n = 1,000–2,000` individuals per condition, `100` replicates per `σ_P` value, `3` independently
  drawn data-generating networks (matches the existing σ_P × W design's replicate/network
  structure in `13_network_full_check.R`, minus the window-length axis, which isn't needed once
  there's no time-indexed slow process in this simulation)

**Outputs.** `GS(ω̂^omit)` vs. `GS(ω̂^adj)`, paired gap with SE; `MAE^omit` vs. `MAE^adj`, paired
gap with SE; example estimated-vs-true network diagrams at the largest `σ_P` (same style as
current Figure 8, averaged across replicates as already fixed in that figure).

**Expected result.** As `σ_P` increases, the symptom-only estimator increasingly absorbs
context variation into apparent edges (`GS^omit` and `MAE^omit` both rise with `σ_P`); the
context-adjusted estimator stays close to `ω_ij` regardless of `σ_P`.

**Figure.** Directly replaces current Figures 7–8, now with a clean causal interpretation (no A vs.
C conflation) and no window-length axis to explain, since there's no time dimension in this
version.

---

## Simulation 3 — Between-group comparisons under mean vs. variance differences in context

**Question.** This is the one that speaks directly to the empirical situation the paper opens with
— what typically happens when researchers compare symptom networks across groups.

**Data-generating model.** Two groups, same true network: `ω_ij^(A) = ω_ij^(B)`, `τ_i^(A) =
τ_i^(B)`, `γ_i^(A) = γ_i^(B)`. Context distributions differ:

    P_n^(A) ~ Normal(μ_A, σ_A²)
    P_n^(B) ~ Normal(μ_B, σ_B²)

**Three cases, crossed with symptom-only vs. context-adjusted estimation:**

| Case | μ_A vs μ_B | σ_A vs σ_B | Expected pattern |
|---|---|---|---|
| 1 | differ | equal | Burden/threshold-level differences (Sim 1's effect); minimal apparent connectivity difference between groups even when omitted |
| 2 | equal | differ | Apparent connectivity differs between groups when context is omitted (higher-σ group looks more connected), despite identical true coupling; gap shrinks under context-adjustment |
| 3 | differ | differ | Both effects superimposed — the case closest to real between-group comparisons in the literature |

**Outputs.** For each case and each group: mean `D`, `Pr(D ≥ N/2)`, `GS^omit`, `GS^adj`, `MAE^omit`,
`MAE^adj`, plus the *between-group difference* in each of these (which is the quantity an applied
researcher would actually report as "groups A and B differ in network structure").

**Why this matters for the paper.** This is the direct answer to Denny's "what people usually find
is no connectivity differences between groups while they do find threshold differences — you might
look at the inverse" comment: Case 1 is exactly that inverse scenario (real threshold-level
separation, no real connectivity difference), and shows whether the symptom-only estimator
correctly reports "no connectivity difference" or falsely detects one.

---

## Pilot grid (run before committing to final parameter values)

Purpose: confirm the model operates in a psychologically plausible range (no saturation to
all-symptoms-off or all-symptoms-on) before spending compute on the full replicate/network sweep.

**Starting parameter values (placeholders pending real Cramer-derived τ_i — see decision #2):**

    N            = 12
    density      = 0.30–0.50 (fraction of possible edges nonzero)
    ω_ij         ~ Uniform(0.20, 0.40) where nonzero
    τ_i          ~ Uniform(−3.0, −1.0), one draw per symptom per replicate network
    γ_i          ~ Uniform(0.6, 1.2)
    μ_P grid     = {−1.0, 0.0, 1.0}         (low / medium / high)
    σ_P grid     = {0, 0.25, 0.5, 1.0}      (0 = degenerate sanity check: omitted-context
                                              distortion should be ≈0 here)
    n per cell   = 1,000 (pilot), 2,000 (final)
    replicates   = 20 (pilot), 100 (final, if compute allows)

**Pilot outputs to check, across the full μ_P × σ_P grid:**
- Mean symptom sumscore `D` and its variance
- Proportion of individuals with all-symptoms-off or all-symptoms-on (should be small — large
  values indicate saturation and mean the τ_i/ω_ij range needs adjusting)
- `GS^omit` and `GS^adj` at each grid point, to confirm the σ_P = 0 row shows ≈0 omitted-context
  gap (this is a correctness check on the code, not just a calibration check — if σ_P = 0 doesn't
  give ≈0 gap, something is wrong in the estimator or the sampler)

## Explicitly out of scope for this pass

Time-varying `P_t`, drift, acute stressors/shocks, fast→slow feedback, and hysteresis — all
deferred to the optional Simulation 4 / dynamic extension, which only gets built once Sims 1–3 are
validated and the team has confirmed the main paper doesn't need it.
