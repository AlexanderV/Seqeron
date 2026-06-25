# Validation Report: CODON-CAI-001 — Codon Adaptation Index (CAI)

- **Validated:** 2026-06-24   **Area:** Codon Optimization
- **Canonical method(s):** `CodonOptimizer.CalculateCAI(string codingSequence, CodonUsageTable table, bool excludeSingleCodonAminoAcids = false)`
  (+ helper `CalculateRelativeAdaptiveness`, set `SingleCodonAminoAcids`), `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

> **2026-06-25 re-validation (post limitation-elimination campaign).** This unit was reset to ⬜ in the
> 2026-06-25 re-reset because the opt-in `excludeSingleCodonAminoAcids` parameter was added during the
> campaign. It has now been **independently re-validated in a fresh context** against external first
> sources retrieved this session, with all CAI values re-derived by hand (Python) and re-run against
> the code. The `excludeSingleCodonAminoAcids` exclusion rule is confirmed verbatim from Jansen et al.
> (2003), PMC2684136, quoting Sharp & Li (1987): "codon families containing a single codon (e.g. AUG
> and UGG in the standard genetic code) should be excluded in computing CAI" because "their
> corresponding w value will always be 1 regardless of codon usage bias of the gene." Stage A is now
> **PASS** (not PASS-WITH-NOTES): the former divergence D-A1 is no longer a divergence — strict
> Sharp & Li/Jansen exclusion is selectable, and the historical inclusive behaviour is the documented
> default. Four new edge-case tests were added this session to close a test-coverage gap on the
> zero-frequency-codon clamp and the no-data-amino-acid skip (see Stage B); CAI fixture now **34 tests,
> all passing**, full `Seqeron.Genomics.Tests` suite **18787 passed / 0 failed**.

---

## Stage A — Description

### Sources opened & what they confirm (re-fetched this session)
- **Wikipedia — "Codon Adaptation Index"** (re-fetched 2026-06-25). Confirms verbatim:
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
  documented (avoids `ln 0 = -∞`); bounded; benign for L > 1. ✓ — now covered by added tests.
- **No-data amino acid** (`f_max ≤ 0`): `w = NaN`, codon skipped, not counted in L. ✓ — now tested.
- **Single-codon AAs (Met/Trp):** selectable via `excludeSingleCodonAminoAcids`. With `true`,
  Met/AUG and Trp/UGG are excluded exactly as Sharp & Li (1987)/Jansen (2003) prescribe; with the
  default `false`, the historical inclusive (w=1.0) behaviour is preserved. **No longer a divergence.**

### Independent cross-check (hand computation, exact — Python, this session)
Recomputed against the in-code Kazusa E. coli K12 table (`exp(mean(ln w))`):
- `CUAACU`: w = [0.04/0.50=0.08, 0.16/0.44=0.36364] → **0.17056057308448833** ✓
- `AGAAGG`: w = [0.04/0.40=0.10, 0.02/0.40=0.05] → **0.07071067811865474** ✓
- `CUGCUA`: w = [1.0, 0.08] → **0.282842712474619** ✓
- 4×CUG + 1×CUA: → **0.6034176336545163** ✓

**Single-codon-AA exclusion cross-check (inclusive vs exclusive):**
| Sequence | Inclusive (default) | Exclusive (`true`) | Effect |
|---|---|---|---|
| `AUGUGG` (Met+Trp only) | **1.0** | **0.0** (L=0, no scored codons) | exclusion fired |
| `AUGCUACUA` | **0.18566355334451112** | **0.08** (Met dropped, two CUA) | exclusion fired |
| `AUGUGGCUA` | **0.43088693800637673** | **0.08** (Met+Trp dropped, one CUA) | exclusion fired |
| `CUGCUA` (no Met/Trp) | **0.282842712474619** | **0.282842712474619** | unchanged (correct) |
All confirm geometric (not arithmetic) mean, f/f_max normalisation, and that exclusion only removes
single-codon families. These exact values match the fixture assertions to ≤1e-10.

**Zero-frequency clamp cross-check** (custom table from reference `"CUG"`, Leu={CUG:1.0}):
- `CUACUG` → w=[1e-6, 1.0] → exp((ln1e-6+ln1)/2) = **0.001** ✓
- `CUA` → w=[1e-6] → **1e-6** ✓
- `UUUCUG` (UUU=Phe, no data → skipped) → only CUG → **1.0** ✓

### Findings / divergences
None at the description level. The single-codon-AA exclusion is now a first-class, sourced option
(former D-A1 resolved). The `1e-6` clamp and `NaN` skip are documented, bounded fallbacks for partial
custom tables (a real-world case Sharp & Li did not face) and are mathematically benign.

---

## Stage B — Implementation

### Code path reviewed
`CalculateCAI` (`CodonOptimizer.cs:473`–`504`), `CalculateRelativeAdaptiveness` (`:506`–`522`),
and the derived `SingleCodonAminoAcids` set built in the static ctor (`:131`–`144`).

### Formula realised correctly? (evidence)
- `maxFreq = synonymousCodons.Max(...)`, `w = Math.Max(codonFreq / maxFreq, 1e-6)` — exactly
  `w_i = f_i / f_max` with the documented 1e-6 clamp. ✓
- `logSum += Math.Log(w); return Math.Exp(logSum / count);` — exact `exp((1/L) Σ ln w_i)`, the
  geometric mean. ✓
- Stop codons (`"*"`) `continue`-skipped → excluded from L. ✓
- Non-standard / no-data AAs → `NaN`, skipped (does not corrupt logSum/count). ✓ (now tested)
- Empty / null / no-complete-codon / all-stop → returns 0. ✓
- `T→U` and `ToUpperInvariant()` normalisation applied. ✓
- `SingleCodonAminoAcids` is **derived** (`AminoAcidToCodons` groups of size 1, excluding `"*"`),
  not hard-coded → {M, W} for the standard code. ✓
- `excludeSingleCodonAminoAcids` (default `false`): when `true`, Met/Trp `continue`-skipped before
  scoring (`:494`), exactly the Sharp & Li/Jansen exclusion; when `false`, counted with w=1.0. ✓

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
`CodonOptimizer_CAI_Tests.cs` — **34 tests** (30 pre-existing + 4 added this session). All assert
**exact sourced numeric values** (Sharp & Li / Jansen / Kazusa hand-computation), not no-throw /
tautologies / code-echoes; deterministic. Covers empty/null, single Met/Trp/both, all-optimal, rare,
range invariant, organism specificity (E. coli/Yeast/Human), DNA/lowercase, stop exclusion (incl.
mid-sequence), geometric-mean sensitivity & monotonicity, hand-calculated, performance,
incomplete-codon, and the full exclusion mode (default==explicit-false, all-single-codon→0, Met
dropped, Met+Trp dropped, no-Met/Trp-unchanged).

**Coverage gap closed this session** — the previous suite had no test for the documented `1e-6`
zero-frequency clamp nor the `NaN` no-data-AA skip (both Stage-A edge cases). Added 4 tests using a
public partial table (`CreateCodonTableFromSequence("CUG", …)`):
`AbsentCodonWithPresentSynonym_ClampsWeightToEpsilon` (CUACUG→0.001),
`AllCodonsAbsentFromFamily_ClampsToEpsilon` (CUA→1e-6),
`AminoAcidWithNoFrequencyData_IsSkipped` (UUUCUG→1.0),
`AllCodonsHaveNoFrequencyData_ReturnsZero` (UUUUUU→0). All four hand-derived, exact.
**Run result: CAI fixture Passed 34, Failed 0; full Seqeron.Genomics.Tests 18787 passed, 0 failed.**

### Findings / defects
None at the implementation level. The code faithfully realises the validated description: w=f/f_max,
geometric mean, stop & (optionally) single-codon-AA exclusion, clamp/skip fallbacks. The only change
this session was additive (4 edge-case tests); no production code was touched.

---

## Verdict & follow-ups
- **Stage A: PASS** — `w = f/f_max`, geometric mean, stop exclusion, and the opt-in single-codon-AA
  exclusion are all correct and externally sourced (Wikipedia + Jansen 2003/PMC2684136 quoting
  Sharp & Li 1987). The clamp/skip fallbacks for partial tables are documented and benign. The former
  D-A1 divergence is resolved: strict Sharp & Li/Jansen behaviour is selectable, historical inclusive
  behaviour is the documented default.
- **Stage B: PASS** — implementation exactly matches the validated description; every cross-check value
  reproduced by independent hand-computation (≤1e-10); 34 tests pass; full suite 18787/0.
- **End state: CLEAN** — no defect. Added 4 edge-case tests to close a coverage gap on the
  zero-frequency clamp and no-data-AA skip. No production code touched.

### References
- Sharp, P.M. & Li, W.H. (1987). *Nucleic Acids Res.* 15(3):1281–1295. DOI 10.1093/nar/15.3.1281.
- Jansen, R. et al. (2003) "An Improved Implementation of the CAI", *Nucleic Acids Res.* (PMC2684136).
- Wikipedia: "Codon Adaptation Index" (fetched 2026-06-24).
- Kazusa Codon Usage Database (E. coli K12 / S. cerevisiae / H. sapiens frequency tables).
