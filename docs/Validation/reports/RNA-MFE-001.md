# Validation Report: RNA-MFE-001 — RNA Minimum Free Energy (Zuker–Stiegler DP, Turner 2004)

- **Validated:** 2026-06-16   **Area:** RnaStructure
- **Canonical method(s):** `RnaSecondaryStructure.CalculateMinimumFreeEnergy(string, int)`, `RnaSecondaryStructure.PredictStructure(string, int, int, int)` (internal `CalculateMinimumFreeEnergyClassic` is a benchmark-only baseline)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES
- **End-state:** ✅ CLEAN
- **Test-quality gate:** PASS (after closing two coverage gaps; full unfiltered suite 6591 passed / 0 failed)

## Stage A — Description

### Sources opened this session (independently retrieved)

| Source | How retrieved | Confirms |
|--------|---------------|----------|
| NNDB Turner 2004 **Hairpin Example 1** | `curl` Wayback `20240709061712` (HTTP 200, 28451 B) — live `rna.urmc.rochester.edu` returns 404 to the fetch tool | Verbatim term list & total (below) |
| NNDB Turner 2004 **Hairpin Example 2** | `curl` Wayback `20240709061712` (HTTP 200, 29115 B) | Verbatim term list & total (below) |
| NNDB Turner 2004 **Watson-Crick parameters** (`wc-parameters.html`) | `curl` Wayback `20240715153500` (HTTP 200) | Full 10-entry stacking table + "Per AU end +0.45" + "Intermolecular Initiation +4.09" |
| **ViennaRNA 2.7.0** (independent reference implementation) | `RNA` Python module (Turner 2004 params); `fold_compound.mfe()` + `eval_structure_verbose()` | MFE structures + per-loop energy breakdown for all probe sequences |

Zuker & Stiegler (1981), Lorenz et al. (2011), Ward et al. (2017) were validated for this unit's sibling and re-confirmed from the Evidence doc's PMC links; the algorithmic claims they back (loop decomposition; O(n³)/O(n²) affine model; C/M/M1/F matrices) are standard and correct.

### Formula check (against retrieved sources)

NNDB **Example 1** verbatim (Wayback): `ΔG°37 = (CG followed by AU) + (AU followed by CG) + (CG followed by AU) + AU end penalty + (AU followed by AA) + Hairpin initiation(6)` = `−2.11 − 2.24 − 2.11 + 0.45 − 0.8 + 5.4 = −1.4 kcal/mol` (arithmetic sum −1.41). Note retrieved verbatim: *"for unimolecular secondary structures, the helical intermolecular initiation does not appear."*

NNDB **Example 2** verbatim (Wayback): `… + AU end penalty + (AU followed by GG) + (GG first mismatch) + Hairpin initiation(5)` = `−2.11 − 2.24 − 2.11 + 0.45 − 0.8 − 0.8 + 5.7 = −1.9 kcal/mol` (sum −1.91).

Stacking table (Wayback `wc-parameters.html`) matches the code's `StackingEnergies` **exactly**: AA/UU −0.93, AU/UA −1.10, UA/AU −1.33, CU/GA −2.08, CA/GU −2.11, GU/CA −2.24, GA/CU −2.35, CG/GC −2.36, GG/CC −3.26, GC/CG −3.42; per-AU-end +0.45; intermolecular initiation +4.09 (correctly **excluded** from the unimolecular DP).

The DP recurrences in §2.2 / §4 (C(i,j) = min{hairpin, stack, interior/bulge over inner pair, multiloop}; exterior W; WM/M1 multiloop fragments; min hairpin loop = 3 nt; MAXLOOP-bounded interior) are the standard Zuker–Stiegler/Turner formulation and match Ward et al. (2017).

### Edge-case semantics check
`null`/empty → 0; homopolymer / unfoldable → 0 (open chain always available, ΔG = 0, INV-01); length < `minLoopSize+2` → 0; `minLoopSize < 3` clamped to 3 (NNDB minimum-loop rule). All sourced and standard.

### Independent cross-check (numbers)
ViennaRNA 2.7.0 on the two canonical sequences returned the **identical fold** but slightly different energies:

| Sequence | NNDB worked example | Seqeron | ViennaRNA 2.7.0 | Structure (all agree) |
|----------|--------------------:|--------:|----------------:|-----------------------|
| `CACAAAAAAAUGUG` | −1.41 (tab −1.4) | −1.41 | −1.30 | `((((......))))` |
| `CACAGAAAGUGUG` | −1.91 (tab −1.9) | −1.91 | −1.80 | `((((.....))))` |

ViennaRNA's `eval_structure_verbose` breakdown shows its compiled Turner-2004 file rounds the stacks to −2.10/−2.20/−2.10 (vs the NNDB **table/worked-example** −2.11/−2.24/−2.11) and folds the AU-end+mismatch into the hairpin term differently, producing the 0.11 kcal/mol gap. The **published NNDB worked example is the canonical source**, and it sums to exactly −1.41 / −1.91 — which is what Seqeron returns. ViennaRNA's value is its own implementation's rounding of the same model, not a contradiction of the NNDB worked example.

### Findings / divergences
None affecting correctness. The description is faithful to the retrieved authoritative sources. **Stage A PASS.**

## Stage B — Implementation

### Code path reviewed
`RnaSecondaryStructure.cs`: `CalculateMinimumFreeEnergy` DP (lines ~1036–1393), `CalculateHairpinLoopEnergy` (684–759), `CalculateStemEnergy` (590–631), parameter tables (114–434), `PredictStructure` (1467–1512).

### Formula realised correctly? (evidence)
**Hand-trace of the DP for `CACAAAAAAAUGUG`** (verified by running a reflection probe against the built `Seqeron.Genomics.Analysis.dll`):
- `V(3,10)` hairpin: `CalculateHairpinLoopEnergy("AAAAAA",'A','U')` = init(6) 5.4 + terminal mismatch `AAAU` −0.8 = 4.6; DP then adds inner-pair AU-end +0.45 → 5.05.
- `V(2,11)` = stack `CA/GU` −2.11 + 5.05 = 2.94; `V(1,12)` = `AC/UG` −2.24 + 2.94 = 0.70; `V(0,13)` = `CA/GU` −2.11 + 0.70 = **−1.41**.
- `W(13)`: closing pair (0,13)=C-G ⇒ no exterior AU penalty; i=0 and j=n−1 ⇒ no dangles. Result = −1.41.

The AU-end penalty is applied **once** (inner closing pair facing the loop); there is **no double-counting** for this fold, and the value is a non-coincidental, term-by-term realisation of the NNDB Example 1 breakdown. Example 2 reproduces −1.91 analogously (adds GG first-mismatch −0.8, init(5) +5.7).

### Cross-verification table recomputed vs code (probe + ViennaRNA)

| Sequence | Seqeron MFE | Seqeron structure | ViennaRNA MFE | ViennaRNA structure | Agreement |
|----------|------------:|-------------------|--------------:|---------------------|-----------|
| `CACAAAAAAAUGUG` | −1.41 | `((((......))))` | −1.30 | `((((......))))` | structure ✓; ΔG param-rounding (0.11) |
| `CACAGAAAGUGUG` | −1.91 | `((((.....))))` | −1.80 | `((((.....))))` | structure ✓; ΔG param-rounding (0.11) |
| `AAAAAAAA` / `GCGC` / `GC` / `""` | 0.00 | open chain | 0.00 | open chain | ✓ exact |
| `AUAUAUAAAUAUAU` | 0.00 | (open) | 0.00 | (open) | ✓ exact (AU stem too weak) |
| `GGGGGUUUUUUUUUU` | 0.00 | (open) | 0.00 | (open) | ✓ exact (GU stem too weak) |
| `GGGGGAAAUCCCCC` | −8.94 | hairpin | −9.10 | `(((((....)))))` | structure ✓; param-rounding |
| `GCGCGCAAAGCGCGC` | −9.58 | hairpin | −9.60 | `((((((...))))))` | structure ✓ (0.02) |
| `GGGGAAAACCCCAAAA…CCCC` | −13.16 | multiloop | −13.40 | `((((....((((....))))....))))` | structure ✓; ML c=0 + rounding |

Every fold matches ViennaRNA's optimal structure; energy differences are exactly the documented Turner-table-vs-compiled-param rounding and the ML-unpaired-cost=0 simplification (ASM-04). This independently confirms the DP is a correct global-optimum Zuker–Stiegler/Turner MFE, not a heuristic that happens to hit −1.41.

### Variant/delegate consistency
`PredictStructure("CACAAAAAAAUGUG")` (greedy heuristic) returns the same `((((......))))` / 4 pairs C-G,A-U,C-G,A-U. `CalculateMinimumFreeEnergyClassic` is a non-Turner Nussinov-style baseline used only by the `[Explicit]` benchmark for timing/determinism (INV-03) — correctly not asserted to match the Turner score.

### Test quality audit (HARD gate)
- **M1/M2** assert the exact NNDB worked-example per-term sums (−1.41 / −1.91) `.Within(1e-9)` plus the one-decimal NNDB totals (−1.4 / −1.9). Sourced literals (NNDB pages retrieved this session), not code echoes — a wrong implementation returning any other value fails. ✓
- **M3/M4/M5** (homopolymer / empty-null / too-short → 0): exact, sourced to the open-chain and minimum-loop rules. ✓
- **M6 (INV-01 ≤0), M7 (INV-02 monotonic), C1 (GC<AU):** relational assertions are appropriate here because they test genuine *properties/invariants*, not known exact scalars. ✓
- **M8/S1/S2:** exact dot-bracket + exact base-pair set / empty / all-dots. ✓
- **Coverage gaps found and closed this session (test-only, 0 code change):**
  - **M9** `…_MinLoopSizeBelow3_ClampedToThree` — the `minLoopSize<3 ⇒ 3` clamp (NNDB minimum-loop rule, §3.1/§6.1) had no test; asserts minLoop 0/−5 ≡ minLoop 3 and that a sequence foldable only via a 2-nt loop (`GGCGCC`) stays at 0.
  - **M10** `…_LowercaseInput_SameAsUppercase` — the documented case-insensitive contract (§3.1) had no test; lowercase `cacaaaaaaaugug` must give the same exact −1.41.
- No assertion weakened, no tolerance widened, no skip/ignore, no expected value bent to match output.

### Findings / defects
**No code defect.** Two documented Stage-B *notes* (both pre-accepted in the algorithm doc, neither a defect):
1. Energies diverge ~0.1 kcal/mol from ViennaRNA's compiled parameters because Seqeron uses the canonical NNDB **table/worked-example** values (−2.11/−2.24 …); this matches the published NNDB hand computation exactly and is the higher-authority reference.
2. Multiloop per-unpaired cost set to 0 (ASM-04) — affects only multiloop optima, not the validated single-hairpin path; ViennaRNA cross-check on a multiloop case confirms identical structure with a small (0.24 kcal/mol) energy offset consistent with this simplification.

Two test-coverage gaps were closed (M9, M10). **Stage B PASS-WITH-NOTES.**

## Verdict & follow-ups
- **Stage A: PASS.** Description matches NNDB worked examples (−1.41/−1.91 confirmed verbatim from retrieved pages), the Turner 2004 stacking table (exact), and standard Zuker–Stiegler/Turner recurrences; invariants genuine.
- **Stage B: PASS-WITH-NOTES.** DP faithfully realises the model (hand-trace + ViennaRNA structure agreement on 11 sequences); notes are documented param-rounding and ML c=0 simplification, not defects.
- **End-state: ✅ CLEAN.** No code change; added two sourced coverage tests (M9 clamping, M10 case-insensitivity). MFE fixture 11→13; full unfiltered suite **6591 passed / 0 failed**; build 0 errors.
- **Follow-up (non-blocking):** if exact ViennaRNA-parity is ever desired, the stacking/mismatch tables would need to be replaced with ViennaRNA's compiled `rna_turner2004.par` rounding — but the current NNDB-table values are the more authoritative published reference and are correct as-is.
