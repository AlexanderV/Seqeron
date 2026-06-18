# Validation Report: ONCO-SIG-004 — Mutational Process Classification

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.ClassifyMutationalProcess(exposures, contributionCutoff)`, `OncologyAnalyzer.GetMutationalProcess(signatureLabel)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (test-coverage gaps closed in this session; no algorithm defect)

## Stage A — Description

### Sources opened this session (independent of the repo artefacts)

| Source | Retrieval | What it confirms |
|--------|-----------|------------------|
| deconstructSigs `whichSignatures.R` (raw on GitHub) | WebFetch of `https://raw.githubusercontent.com/raerose01/deconstructSigs/master/R/whichSignatures.R` | Default `signature.cutoff = 0.06`; cutoff applied verbatim as `weights[weights < signature.cutoff ] <- 0` — comparison operator is strict `<` (so exactly 0.06 is **retained**). |
| Wikipedia "Mutational signatures" (with primary refs) | WebFetch | SBS1 = spontaneous deamination of 5-methylcytosine (age-correlated/clock-like); SBS2 & SBS13 = AID/APOBEC cytidine-deaminase activity; SBS4 = tobacco smoking; SBS5 = age-related (clock-like); SBS6/SBS15/SBS20/SBS26 = DNA mismatch-repair deficiency / MSI; SBS7a–d = ultraviolet radiation. |
| COSMIC SBS catalogue + primary literature (SBS20, SBS5/SBS1) | WebSearch | SBS20 proposed aetiology verbatim "Concurrent POLD1 mutations and defective DNA mismatch repair" (one of the MSI/MMR-deficiency set); SBS5 aetiology "Unknown (clock-like signature)"; SBS1 "spontaneous deamination of 5-methylcytosine", both clock-like (mutation count ∝ age). |

### Formula check
- Normalized relative contribution `wᵢ = eᵢ / Σe` — matches deconstructSigs "weights W normalized between 0 and 1". ✓
- Presence rule `wᵢ ≥ τ`, `τ = 0.06`, strict-`<` exclusion — matches the fetched reference line `weights[weights < signature.cutoff] <- 0` with `signature.cutoff = 0.06`. ✓
- Per-process aggregation by summation of surviving member contributions — additive deconstructSigs weights; reasonable and documented assumption (ASM-01). ✓
- Dominant process = argmax of per-process aggregated contribution. ✓

### SBS → process map check (every mapped label re-grounded externally)
SBS1→Aging, SBS5→Aging (clock-like), SBS2→APOBEC, SBS13→APOBEC, SBS4→Tobacco, SBS7a/7b/7c/7d→UV, SBS6/SBS15/SBS20/SBS26→MMR deficiency — **all confirmed** against COSMIC aetiology strings / Wikipedia primary refs. No mis-mapping found. (SBS20 is correctly grouped under MMR deficiency: its aetiology is concurrent POLD1 proofreading loss *and* defective MMR.)

### Edge-case semantics check
- Σ exposure = 0 (or empty list) ⇒ normalization undefined ⇒ no active processes, dominant = Unknown — sound. ✓
- Sub-cutoff signatures dropped ⇒ surviving contributions may sum < 1 (remainder = "unknown") — matches Rosenthal 2016. ✓
- Unmapped labels (e.g. SBS99) contribute to no recognized process — sound, COSMIC-consistent. ✓
- Negative / NaN exposure invalid; cutoff constrained to `[0,1)` — defensible contract. ✓

### Independent cross-check (hand computation against the sourced rules)
Dataset {SBS2:50, SBS13:30, SBS1:15, SBS4:5}, Σ=100 ⇒ w = {0.50, 0.30, 0.15, 0.05}. SBS4 0.05 < 0.06 → dropped. APOBEC = 0.50+0.30 = **0.80**; Aging = **0.15**; Tobacco = 0. Active = {APOBEC, Aging}; dominant = **APOBEC (0.80)**. Matches the algorithm. Boundary: SBS4 6/100 = 0.06 is **not** `< 0.06` ⇒ retained (matches strict-`<` source rule).

### Findings / divergences
None. Every formula, threshold, comparison operator, and aetiology mapping traces to an external source retrieved this session.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`:
- `MutationalProcess` enum (3138–3166), `SbsToProcess` map (3176–3192).
- `GetMutationalProcess` (3203–3209): case-insensitive dictionary lookup, Unknown fallback, null guard.
- `ClassifyMutationalProcess` (3271–3339): validates cutoff ∈ [0,1) (NaN rejected first), validates each exposure (NaN/negative → throw), computes `total`, short-circuits `total ≤ 0` → empty/Unknown, then per-signature `wᵢ = eᵢ/total`, `continue` if `wᵢ < cutoff` (strict, exactly matching `weights[weights<signature.cutoff]<-0`), `continue` if process Unknown, accumulates per process, orders descending by contribution then by process enum, returns dominant = first.

### Formula realised correctly?
Yes. The code is a faithful, exact realisation of the validated description: `eᵢ/Σe`, strict-`<` 0.06 cutoff, COSMIC map, additive per-process aggregation, argmax dominant. The map entries match the externally-confirmed aetiologies one-for-one.

### Cross-verification table (recomputed vs code, via the green test run)

| Input | External-sourced expectation | Code result | Match |
|-------|------------------------------|-------------|-------|
| {SBS2:50,SBS13:30,SBS1:15,SBS4:5} | APOBEC 0.80, Aging 0.15, Tobacco dropped; dominant APOBEC | same | ✓ |
| {SBS1:94,SBS4:6} | SBS4 w=0.06 retained, Tobacco=0.06 | same | ✓ |
| {SBS1:94.001,SBS4:5.999} | SBS4 w≈0.05999 < 0.06 excluded | same | ✓ |
| {SBS6:25,SBS15:25,SBS20:25,SBS26:25} | one MMR process = 1.0 | same | ✓ |
| {SBS7a:40,SBS7b:60} | UV = 1.0 | same | ✓ |
| {SBS1:50,SBS2:50} | tie 0.50/0.50, dominant = Aging (lower enum) | same | ✓ |
| GetMutationalProcess(SBS1/5/2/13/4/7a/6/15/20/26) | Aging,Aging,APOBEC,APOBEC,Tobacco,UV,MMR×4 | same | ✓ |

### Variant/delegate consistency
`ClassifyMutationalProcess` delegates label→process resolution to `GetMutationalProcess`; the two are consistent (the classification test set and the lookup `[TestCase]`s agree).

### Test quality audit (HARD gate)
- **Sourced expectations, not code echoes:** all expected values (0.80, 0.15, 0.06, 1.0, the 10-label map, the strict-`<` boundary, tie→Aging) trace to the deconstructSigs cutoff rule and COSMIC aetiologies — re-confirmed externally this session, not read off the implementation. A deliberately-wrong impl (e.g. `<=` cutoff, summing instead of dropping sub-cutoff, mis-mapping SBS20) would fail M3/M4/M8/M9.
- **No green-washing:** exact `Is.EqualTo(...).Within(1e-10)` on all numeric contributions; `Is.EquivalentTo` on the exact active set; exact enum equality on dominant/lookup. No `Greater`/`AtLeast`/range/`Contains` where an exact value is known; no widened tolerance; no skipped/ignored test.
- **Coverage gaps found and closed (Stage-B defect, test-quality only):** the original 30-test fixture left two real code branches unexercised — (1) the documented **NaN-exposure** guard (`double.IsNaN(exposure) → ArgumentException`, contract §3.3) and the **NaN-cutoff** guard (`double.IsNaN(contributionCutoff) → ArgumentOutOfRangeException`), and (2) the **tie-break ordering** branch `.ThenBy(p => p.Process)` (two processes with equal aggregated contribution). Added three tests: `ClassifyMutationalProcess_NaNExposure_Throws`, `ClassifyMutationalProcess_NaNCutoff_Throws`, and `ClassifyMutationalProcess_EqualContributions_OrderedByProcessEnum` (Aging/APOBEC tie at 0.50 → dominant Aging by enum). Fixture 30 → 33.
- **Honest green:** full unfiltered `dotnet test` = **6644 passed, 0 failed** (1 pre-existing benchmark skipped, unrelated). `dotnet build` 0 errors; the changed test file is warning-free (the 4 build warnings are pre-existing, in `ApproximateMatcher_EditDistance_Tests.cs`).

### Findings / defects
No algorithm defect. The implementation matches the externally validated description exactly. The only finding is the test-coverage gap above, fully fixed in this session (logged FINDINGS A47).

## Verdict & follow-ups
- **Stage A: PASS** — biology/maths independently confirmed against deconstructSigs source (cutoff 0.06, strict `<`) and COSMIC/primary-literature aetiologies (all 13 mapped labels).
- **Stage B: PASS-WITH-NOTES** — code is correct; two untested branches (NaN guards, tie-break ordering) were covered with exact-sourced/derivation-locked tests.
- **End-state: CLEAN** — no defect outstanding; suite 6644 green, build clean.
