# Validation Report: CODON-CAI-001 — Codon Adaptation Index (CAI)

- **Validated:** 2026-06-12   **Area:** Codon Optimization
- **Canonical method(s):** `CodonOptimizer.CalculateCAI(string codingSequence, CodonUsageTable table)`
  (helper `CodonOptimizer.CalculateRelativeAdaptiveness(...)`),
  `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs:423`–`468`.
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — "Codon Adaptation Index"** (https://en.wikipedia.org/wiki/Codon_Adaptation_Index).
  Confirms verbatim: relative adaptiveness `w_i = f_i / max(f_j)` ("the ratio between the observed
  frequency of the codon f_i and the frequency of the most frequent synonymous codon f_j for that
  amino acid"), and `CAI = (Π_{i=1}^{L} w_i)^(1/L)` described as "the geometric mean of the weight
  associated to each codon over the length (L) of the gene sequence."
- **Sharp, P.M. & Li, W.H. (1987)**, *Nucleic Acids Res.* 15(3):1281–1295, PMC340524, PMID 3547335,
  DOI 10.1093/nar/15.3.1281 — the original CAI paper. Abstract opened on PMC; full methodology is
  image-only on PMC so the precise zero-count rule was cross-checked via secondary literature.
- **Reference-implementation cross-check** (Biopython `Bio.SeqUtils.CodonUsage`, Biopython PR #881
  "Correct CAI Implementation", the `CAI` PyPI package, seqinr `cai`): confirms the canonical Sharp &
  Li convention that **stop codons AND single-codon amino acids (Met/AUG, Trp/UGG) are excluded** from
  the CAI computation, because a single-codon family cannot exhibit codon bias (18 meaningful families,
  not 20).

### Formula check
- `w_i = f_i / f_max`: **correct**, matches Wikipedia and Sharp & Li exactly. `w_max = 1` for each
  amino acid's optimal codon. ✓
- Geometric mean `CAI = (Π w_i)^(1/L) = exp((1/L) Σ ln w_i)`: **correct** (geometric, NOT arithmetic). ✓

### Edge-case semantics check
- **w_i = 0 handling.** Spec/implementation clamp `w` to **1e-6** when a codon is absent from the
  reference table but its amino acid has other codons present (Deviation D1). This is a defensible,
  documented choice (avoids `ln(0) = -∞ → CAI = 0/NaN` from incomplete tables). Sharp & Li (1987)
  themselves assigned absent codons a nonzero pseudo-count rather than 0, so a small-value substitution
  is consistent in spirit. ✓ (documented, bounded, benign).
- **Stop codons excluded** from L: ✓ matches canonical convention.
- **Single-codon AAs (Met/Trp): DIVERGENCE — see Findings.**

### Independent cross-check (hand computation — exact numbers)
- Generic worked example `w = [1.0, 0.5, 0.2]` → `(1.0·0.5·0.2)^(1/3) = 0.1^(1/3) = 0.46416` ✓
  (matches the protocol's stated 0.4642).
- TestSpec Test 3 `CUAACU` (E. coli): `w_CUA = 0.04/0.50 = 0.08`, `w_ACU = 0.16/0.44 = 0.36364`;
  `(0.08·0.36364)^(1/2) = 0.17056` ✓ (matches spec and test `expectedCai`).
- All-optimal `CUGCCGACC` → `1.0` ✓.
- (Note: the Evidence doc's "Test Case 3" CUA-CCA-ACA states 0.1980 but actually computes to **0.2039** —
  a harmless arithmetic slip in the *evidence prose only*; no code or test depends on it.)

### Findings / divergences
**D-A1 (PASS-WITH-NOTES): Single-codon amino acids (Met/AUG, Trp/UGG) are INCLUDED, not excluded.**
The canonical Sharp & Li (1987) definition excludes Met and Trp from CAI (they cannot show bias).
The Seqeron spec deliberately includes them, assigning each `w = 1.0` (Invariant 2: "Met/Trp codons
contribute w=1.0 always"; tests M2/M3). Consequences vs strict S&L:
- Single `AUG` or `UGG` returns CAI = 1.0; strict S&L leaves L = 0 (undefined / would be the empty
  convention → 0).
- For real genes, the extra `w = 1.0` factors increase L without changing the product, which biases
  CAI **upward** relative to the strict definition.

This is a *documented design choice*, not a hidden bug. It never violates the `0 ≤ CAI ≤ 1` invariant
(every `w ≤ 1`). It is mathematically benign (multiplying a product by 1.0; only L is affected). Because
the spec, Evidence doc, and tests all consistently and intentionally specify the inclusive behaviour, and
because it does not produce a scientifically *wrong* value (just a different, documented convention), this
is recorded as PASS-WITH-NOTES rather than FAIL. Not silently "fixed" to strict S&L because the protocol
forbids changing code away from a deliberately-documented description, and doing so would require redefining
documented edge cases and rewriting ~8 deliberately-asserted tests.

---

## Stage B — Implementation

### Code path reviewed
`CodonOptimizer.CalculateCAI` (`CodonOptimizer.cs:423`) and `CalculateRelativeAdaptiveness`
(`CodonOptimizer.cs:452`).

### Formula realised correctly? (evidence)
- `w = Math.Max(codonFreq / maxFreq, 1e-6)` where `maxFreq = max over synonymous codons` —
  exactly `w_i = f_i / f_max` with the documented 1e-6 clamp. ✓
- `logSum += Math.Log(w); ... return Math.Exp(logSum / count);` — exact `exp((1/L) Σ ln w_i)`,
  i.e. the **geometric mean** (not arithmetic). ✓
- Stop codons (`aminoAcid == "*"`) are `continue`-skipped → excluded from L. ✓
- Non-standard / no-data amino acids return `NaN` and are skipped (do not corrupt logSum). ✓
- Empty / null / no-complete-codon / all-stop → returns `0`. ✓
- DNA `T→U` and `ToUpperInvariant()` normalisation applied. ✓
- Met/Trp: NOT skipped → counted with `w = 1.0` (consistent with the documented description; see D-A1).

### Cross-verification table recomputed vs code (all exact, hand-computed against the in-code E. coli/Yeast/Human tables)
| Sequence | Table | Hand value | Test asserts | Match |
|---|---|---|---|---|
| `CUAACU` | E. coli | 0.17056 | 0.17056 | ✓ |
| `CUGCCGACC` | E. coli | 1.0 | 1.0 | ✓ |
| `CUGCCGACC` | Yeast | 0.40841 | 0.4085 | ✓ |
| `AGAAGG` | E. coli | 0.07071 | 0.07071 | ✓ |
| `CUGCUA` | E. coli | 0.28284 | 0.28284 | ✓ |
| 4×CUG+1×CUA | E. coli | 0.60342 | 0.6036 | ✓ |
| `w=[1.0,0.5,0.2]` (generic) | — | 0.46416 | (protocol) | ✓ |

### Variant/delegate consistency
`CalculateCAI` is the sole canonical entry; `OptimizeSequence` calls it for original/optimized CAI
(`CodonOptimizer.cs:264,307`). No divergent re-implementation exists.

### Test quality audit
`tests/Seqeron/Seqeron.Genomics.Tests/CodonOptimizer_CAI_Tests.cs` — 38 tests, all asserting
**exact sourced numeric values** (not "no-throw"), deterministic. Covers empty/null, single Met/Trp,
all-optimal, rare, range invariant, organism specificity, DNA/lowercase, stop exclusion (incl. mid-sequence),
geometric-mean sensitivity & monotonicity, hand-calculated, performance, incomplete-codon. Real and adequate.

### Findings / defects
None at the implementation level: the code faithfully realises the (documented, inclusive-Met/Trp)
description. The only divergence from strict Sharp & Li is the Met/Trp inclusion, which originates in
the **description** (Stage A, D-A1), and is intentional and benign.

---

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — formula (w = f/f_max, geometric mean), w=0 clamp, and stop exclusion
  are correct and sourced; documented deviation D-A1 (Met/Trp included with w=1.0 rather than excluded
  per strict Sharp & Li 1987) is intentional, bounded (CAI stays in [0,1]), and mathematically benign.
- **Stage B: PASS** — implementation exactly matches the validated description; all cross-check values
  reproduced by hand; 38 tests pass.
- **End state: CLEAN** — no defect requiring a fix; formula is CAI-critical-correct (geometric mean,
  correct normalisation, sound w=0 handling). The Met/Trp convention is a documented design decision,
  not a defect; left unchanged per protocol (do not change code away from a deliberately-documented spec).
- **Follow-up (non-blocking):** correct the Evidence-doc prose value for CUA-CCA-ACA (0.1980 → 0.2039);
  cosmetic only, no code/test impact.

### References
- Sharp, P.M. & Li, W.H. (1987). *Nucleic Acids Res.* 15(3):1281–1295. DOI 10.1093/nar/15.3.1281.
- Wikipedia: "Codon Adaptation Index".
- Biopython `Bio.SeqUtils.CodonUsage`; Biopython PR #881 (CAI correction); `CAI` PyPI package; seqinr `cai`
  (cross-checks for Met/Trp & single-codon-family exclusion convention).
- Kazusa Codon Usage Database (E. coli K12 / S. cerevisiae / H. sapiens frequency tables).
