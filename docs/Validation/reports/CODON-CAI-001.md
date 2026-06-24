# Validation Report: CODON-CAI-001 — Codon Adaptation Index (CAI)

- **Validated:** 2026-06-24   **Area:** Codon Optimization
- **Canonical method(s):** `CodonOptimizer.CalculateCAI(string codingSequence, CodonUsageTable table, bool excludeSingleCodonAminoAcids = false)`
  (+ helper `CalculateRelativeAdaptiveness`, set `SingleCodonAminoAcids`), `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs`.
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

> **2026-06-24 update (limitation removed):** the previously documented limitation D-A1 (single-codon
> amino acids Met/AUG and Trp/UGG included with w=1.0 rather than excluded per Sharp & Li 1987 /
> Jansen 2003) is **resolved**. `CalculateCAI` now accepts an opt-in
> `excludeSingleCodonAminoAcids` parameter (default `false` = unchanged historical behaviour;
> `true` = canonical Sharp & Li / Jansen exclusion of single-codon families). The exclusion rule was
> re-confirmed verbatim from Jansen et al. (2003), PMC2684136, which quotes Sharp & Li (1987):
> "codon families containing a single codon (e.g. AUG and UGG …) should be excluded in computing CAI"
> because "their corresponding w value will always be 1 regardless of codon usage bias." Five new
> evidence-based tests (M13–M17, exact values: AUGUGG→0, AUGCUACUA→0.08, AUGUGGCUA→0.08, default
> unchanged) were added; CAI fixture now 30 tests, all passing. **The checklist Status was reset to
> ☐ pending independent re-validation of the new mode.** The notes below describe the pre-fix state.

This is an independent re-validation (fresh context). The code and tests are byte-identical to the
prior campaign baseline (`git diff cb113ce -- CodonOptimizer.cs CodonOptimizer_CAI_Tests.cs` is empty;
last source change is `cebac6f0`, predating the campaign). Findings below were re-derived from
authoritative sources and re-checked by hand and by running the tests.

---

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — "Codon Adaptation Index"** (fetched 2026-06-24). Confirms verbatim:
  `w_i = f_i / max(f_j)` — "the ratio between the observed frequency of the codon f_i and the
  frequency of the most frequent synonymous codon f_j for that amino acid"; and
  `CAI = (∏ w_i)^(1/L)` — "the geometric mean of the weight associated to each codon over the
  length (L) of the gene sequence (measured in codons)." Confirms CAI is a **geometric mean**.
  Wikipedia (as fetched) does not state the Met/Trp/stop exclusion rule.
- **Sharp, P.M. & Li, W.H. (1987)**, *Nucleic Acids Res.* 15(3):1281–1295, PMID 3547335,
  DOI 10.1093/nar/15.3.1281 — original CAI paper.
- **Jansen et al. (2003), "An Improved Implementation of the Codon Adaptation Index"**
  (PMC2684136) + general literature cross-check: confirms the canonical convention that the
  original CAI is a geometric mean **"excluding stop and start codons"**, and that
  **single-codon amino acids Met (AUG) and Trp (UGG) are excluded** because their w is always 1
  regardless of bias; including them inflates CAI for Met/Trp-rich genes ("This is why it is
  important to exclude such codons from CAI calculations to avoid skewing results.").

### Formula check
- `w_i = f_i / f_max` (max over synonymous codons of the AA): **correct**, matches Wikipedia and
  Sharp & Li exactly; optimal codon → w = 1. ✓
- `CAI = (∏ w_i)^(1/L) = exp((1/L) Σ ln w_i)`: **correct**, geometric (NOT arithmetic) mean. ✓
- Log base is irrelevant (cancels in exp/log pair); code uses natural log consistently. ✓

### Edge-case semantics check
- **Stop codons excluded** from L: ✓ matches canonical convention (`aminoAcid == "*"` skipped).
- **w = 0 handling (Deviation D1):** when a codon is absent from the table (`f = 0`) but the AA has
  other present codons (`f_max > 0`), `w` is clamped to `1e-6` rather than 0. Defensible and
  documented (avoids `ln 0 = -∞`); bounded; benign for L > 1. ✓
- **Single-codon AAs (Met/Trp): DIVERGENCE — see Findings.**

### Independent cross-check (hand computation, exact)
Recomputed against the in-code Kazusa tables (Python `exp(mean(ln w))`):
- `CUAACU` (E. coli): w = [0.04/0.50, 0.16/0.44] → **0.17056** ✓
- `CUGCCGACC` (Yeast): w = [0.11/0.29, 0.12/0.42, 0.22/0.35] → **0.4084** ✓
- `AUGCUGCCGACC` (Human): w = [1, 1, 0.11/0.32, 1] → **0.7657** ✓
- 4×CUG + 1×CUA (E. coli): → **0.6034** ✓
All confirm geometric (not arithmetic) mean and the f/f_max normalisation.

### Findings / divergences
**D-A1 (PASS-WITH-NOTES): single-codon AAs (Met/AUG, Trp/UGG) are INCLUDED with w=1.0, whereas the
canonical Sharp & Li / Jansen convention EXCLUDES them.** This is a deliberate, documented design
choice (TestSpec Invariant 2; tests M2/M3/MetAndTrp). Consequences vs strict S&L:
- `AUG`/`UGG`-only sequences return CAI = 1.0; strict S&L would have L = 0 (empty → 0 by convention).
- For real genes the extra w=1.0 factors leave the product unchanged but increase L, biasing CAI
  **upward** relative to strict S&L.
It never breaks `0 ≤ CAI ≤ 1` (every w ≤ 1) and is mathematically benign. Recorded PASS-WITH-NOTES
(documented convention, not a hidden bug); not "fixed" because the protocol forbids changing code away
from a deliberately-documented spec, and the spec + Evidence + 3 tests all consistently assert it.

Minor (cosmetic, non-blocking): the Evidence-doc prose value for `CUA-CCA-ACA` reads 0.1980 but the
arithmetic yields ~0.2039. No code or test depends on this number.

---

## Stage B — Implementation

### Code path reviewed
`CalculateCAI` (`CodonOptimizer.cs:423`–`450`) and `CalculateRelativeAdaptiveness` (`:452`–`468`).

### Formula realised correctly? (evidence)
- `maxFreq = synonymousCodons.Max(...)`, `w = Math.Max(codonFreq / maxFreq, 1e-6)` — exactly
  `w_i = f_i / f_max` with the documented 1e-6 clamp. ✓
- `logSum += Math.Log(w); return Math.Exp(logSum / count);` — exact `exp((1/L) Σ ln w_i)`, the
  geometric mean. ✓
- Stop codons (`"*"`) `continue`-skipped → excluded from L. ✓
- Non-standard / no-data AAs → `NaN`, skipped (does not corrupt logSum/count). ✓
- Empty / null / no-complete-codon / all-stop → returns 0. ✓
- `T→U` and `ToUpperInvariant()` normalisation applied. ✓
- Met/Trp NOT skipped → counted with w = 1.0 (consistent with documented description; D-A1). ✓

### Cross-verification table recomputed vs code
| Sequence | Table | Hand value | Test asserts | Match |
|---|---|---|---|---|
| `CUAACU` | E. coli | 0.17056 | 0.17056 | ✓ |
| `CUGCCGACC` | E. coli | 1.0 | 1.0 | ✓ |
| `CUGCCGACC` | Yeast | 0.4084 | 0.4085 | ✓ |
| `AGAAGG` | E. coli | 0.07071 | 0.07071 | ✓ |
| `CUGCUA` | E. coli | 0.28284 | 0.28284 | ✓ |
| `AUGCUGCCGACC` | Human | 0.7657 | 0.7656 | ✓ |
| 4×CUG+1×CUA | E. coli | 0.6034 | 0.6036 | ✓ |

### Variant/delegate consistency
`CalculateCAI` is the sole canonical entry; `OptimizeSequence` calls it for original/optimized CAI
(`:264`, `:307`). No divergent re-implementation. ✓

### Test quality audit
`CodonOptimizer_CAI_Tests.cs` — **25 tests** (matches the TestSpec audit count; the prior report's
"38" was inaccurate). All assert **exact sourced numeric values** (not no-throw / tautologies),
deterministic. Covers empty/null, single Met/Trp/both, all-optimal, rare, range invariant, organism
specificity (E. coli/Yeast/Human), DNA/lowercase, stop exclusion (incl. mid-sequence), geometric-mean
sensitivity & monotonicity, hand-calculated, performance, incomplete-codon. Real and adequate.
**Run result: Passed 25, Failed 0.**

### Findings / defects
None at the implementation level. The code faithfully realises the documented (inclusive-Met/Trp)
description. The only divergence from strict Sharp & Li is the Met/Trp inclusion, which originates in
the description (D-A1) and is intentional and benign.

---

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — `w = f/f_max`, geometric mean, w=0 clamp, and stop exclusion are
  correct and externally sourced; documented deviation D-A1 (Met/Trp included with w=1.0 rather than
  excluded per canonical S&L/Jansen) is intentional, bounded ([0,1]), and mathematically benign.
- **Stage B: PASS** — implementation exactly matches the validated description; all cross-check values
  reproduced by hand; 25 tests pass.
- **End state: CLEAN** — no defect requiring a fix (formula is CAI-critical-correct: geometric mean,
  correct normalisation, sound w=0 handling). Met/Trp convention is a documented design decision left
  unchanged per protocol. Code unchanged since prior validation; no code touched this session.
- **Follow-up (non-blocking, cosmetic):** correct the Evidence-doc prose value for `CUA-CCA-ACA`
  (0.1980 → ~0.2039). No code/test impact.

### References
- Sharp, P.M. & Li, W.H. (1987). *Nucleic Acids Res.* 15(3):1281–1295. DOI 10.1093/nar/15.3.1281.
- Jansen, R. et al. (2003) "An Improved Implementation of the CAI", *Nucleic Acids Res.* (PMC2684136).
- Wikipedia: "Codon Adaptation Index" (fetched 2026-06-24).
- Kazusa Codon Usage Database (E. coli K12 / S. cerevisiae / H. sapiens frequency tables).
