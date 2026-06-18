# Validation Report: ONCO-CNA-001 — Copy-Number Alteration Classification

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.Log2RatioToCopyNumber`, `CallCopyNumber`,
  `ClassifyCopyNumber`, `ClassifyCopyNumbers` (file `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm

- **CNVkit `cnvlib/call.py`** (fetched this session from
  `https://raw.githubusercontent.com/etal/cnvkit/master/cnvlib/call.py`). Retrieved verbatim:
  - `do_call` default `thresholds = (-1.1, -0.25, 0.2, 0.7)`.
  - `_log2_ratio_to_absolute_pure`: `ncopies = ref_copies * 2**log2_ratio`.
  - `absolute_threshold` binning loop:
    ```python
    cnum = 0
    for cnum, thresh in enumerate(thresholds):
        if row.log2 <= thresh:
            if ref_copies != ploidy:
                cnum = int(cnum * ref_copies / ploidy)
            break
    else:
        cnum = int(np.ceil(_log2_ratio_to_absolute_pure(row.log2, ref_copies)))
    ```
  - `absolute_threshold` docstring (verbatim): `DEL(0) < -1.1`, `LOSS(1) < -0.25`,
    `GAIN(3) >= +0.2`, `AMP(4) >= +0.7`, derived from `R> log2(2:6 / 4) = -1.0, -0.4150375,
    0.0, 0.3219281, 0.5849625` under 50% clonality, ±0.1 noise.
  - NaN handling (verbatim): `if np.isnan(row.log2): ... absolutes[idx] = ref_copies; continue`
    — replaced with **`ref_copies`** (the neutral reference copies), with a logged warning.
- **GISTIC2.0 (Mermel et al. 2011)** — corroborates the neutral ±0.1 noise band and high-amplitude
  amp/del cutoffs (0.848 / −0.737), consistent with the CNVkit extremes. (Background, not the
  formula source.)

### Formula check

- Absolute copy number `n = ploidy · 2^log2` matches `_log2_ratio_to_absolute_pure` exactly
  (with `ref_copies = ploidy = 2` for autosomes).
- Hard-threshold integer call (first cutoff with `log2 ≤ cutoff`, counting from 0; else
  `ceil(ploidy·2^log2)`) matches the binning loop exactly, including the inclusive `<=` boundary.
- State mapping 0→DeepDeletion, 1→Loss, 2→Neutral, 3→Gain, ≥4→Amplification matches the
  DEL/LOSS/neutral/GAIN/AMP labels.

### Edge-case semantics check

- Boundary inclusivity (`log2 == thresh` → lower state): matches `if row.log2 <= thresh`.
- NaN log2 → neutral reference copy number: matches CNVkit (`ref_copies`); for the diploid
  default `ref_copies = ploidy = 2`.
- Above last threshold → `ceil`, not `round`: matches `np.ceil`.

### Independent cross-check (hand computation, this session)

R-example reproduced exactly with Python: `log2(2:6/4) = -1.0, -0.4150375, 0.0, 0.3219281,
0.5849625`. Default-threshold binning reproduced for every TestSpec value:

| log2 | n = 2·2^log2 | CN (hand) | CN (TestSpec) |
|------|--------------|-----------|---------------|
| −2.0 | 0.5 | 0 | 0 ✓ |
| −1.1 | 0.467 | 0 (boundary) | 0 ✓ |
| −1.0 | 1.0 | 1 | 1 ✓ |
| −0.25 | 1.682 | 1 (boundary) | 1 ✓ |
| 0.0 | 2.0 | 2 | 2 ✓ |
| 0.2 | 2.297 | 2 (boundary) | 2 ✓ |
| log2(3/2)=0.5849625 | 3.0 | 3 | 3 ✓ |
| 0.7 | 3.249 | 3 (boundary) | 3 ✓ |
| 0.8 | 3.482 | ceil→4 | 4 ✓ |
| 1.0 | 4.0 | ceil→4 | 4 ✓ |
| 2.0 | 8.0 | ceil→8 | 8 ✓ |

Absolute-CN formula: 0→2.0, 1→4.0, −1→1.0, log2(3/2)→3.0 — all confirmed.

### Findings / divergences

- The description (TestSpec §1.3 / Evidence) states NaN → "neutral reference copy number (CN 2)".
  CNVkit actually assigns **`ref_copies`** (not a literal 2, not a rounded ploidy). For the
  diploid default these coincide at 2. The doc phrasing "rounded ploidy" is a faithful description
  of the diploid case and is correct for the only ploidy the unit exposes by default. Minor wording
  nuance; does not affect any output. **No defect.**
- TestSpec S1 cites log2 = −0.3 (CN 1 under both default and custom), but the implemented test uses
  log2 = −0.2 (CN 1 custom vs CN 2 default) — a *stronger*, discriminating choice. The test is
  better than the spec row; not a defect.

Stage A: **PASS** — formula, thresholds, boundary convention, NaN handling, and the five-state
mapping all match the cited primary reference implementation, verified by independent hand
computation.

## Stage B — Implementation

### Code path reviewed

- `Log2RatioToCopyNumber` (OncologyAnalyzer.cs:3982) — `ploidy * Math.Pow(2, log2)`, ploidy>0 guard.
- `CallCopyNumber` (4010) — validate thresholds, ploidy guard, NaN → `round(ploidy)`, ascending
  scan returning first `log2 <= cutoff`, else `ceil(Log2RatioToCopyNumber)`.
- `ClassifyCopyNumber` (4053) — integer CN + continuous absolute (ploidy for NaN) + state.
- `ClassifyCopyNumbers` (4077) — null guard, per-element map, order/length preserving.
- `StateFromCopyNumber` (4095) — ≥4→Amplification, 0/1/2/3 → DeepDeletion/Loss/Neutral/Gain.
- `ValidateThresholds` (4117) — exactly 4, no NaN, strictly ascending.

### Formula realised correctly?

Yes. Every line matches the CNVkit reference: `n = ploidy·2^log2`; inclusive ascending bin scan;
`ceil` (not round) above the last threshold; state labels. The full suite re-confirms the
hand-computed cross-check table value-for-value.

### Cross-verification table recomputed vs code

The full unfiltered suite passed (6649/0), including the parametrised `CallCopyNumber` and
`ClassifyCopyNumber` cases covering every row of the table above. All match.

### Variant/delegate consistency

`ClassifyCopyNumbers` delegates per element to `ClassifyCopyNumber`, which delegates the integer
to `CallCopyNumber` and the continuous value to `Log2RatioToCopyNumber` — single source of truth,
no divergence. Batch test M13 confirms order/length preservation.

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** all asserts use exact externally-derived values (CN integers,
  absolute-CN doubles, state enums) traceable to the CNVkit source and hand computation.
- **Mutation-killing:** S3 (`0.8 → CN 4`) specifically kills a `round`-instead-of-`ceil` mutant;
  S1 (`-0.2`) kills a "thresholds ignored / defaults always used" mutant; boundary cases kill an
  exclusive-`<` mutant; M12 kills a "NaN → throw / 0" mutant.
- **Coverage:** all four public methods exercised; M1–M13, S1–S3, C1, E1–E5 all present; every
  Stage-A branch (each bin, both boundaries, NaN, above-last-ceil, custom thresholds) and every
  documented error (null/empty/wrong-count/non-ascending/NaN-threshold/non-positive-ploidy) covered.
- **No green-washing:** the only inequality assert is S2 monotonicity (`GreaterThanOrEqualTo`),
  which is the correct shape for an order-invariant property and is backed by exact-value tests
  (M1–M6); permitted by protocol for invariants.
- **Honest green:** full unfiltered `dotnet test` = **Failed: 0, Passed: 6649**; `dotnet build`
  0 errors (4 pre-existing warnings in unrelated files, none in CNA code/tests).

Minor coverage note (not a defect): no `ClassifyCopyNumber` parametrised row reaches the high
amplification CN 8 (log2 2.0); that path is covered via `CallCopyNumber` M6 plus the CN≥4 state
mapping. The continuous absolute-CN for the ceil case (0.8 → 3.482) is not separately asserted,
but the formula is covered by M11.

### Findings / defects

None. No code or test change required.

## Verdict & follow-ups

- **Stage A: PASS. Stage B: PASS.**
- **End state: ✅ CLEAN** — no defect found; algorithm fully functional and faithfully matches the
  validated CNVkit `absolute_threshold` / `_log2_ratio_to_absolute_pure` description.
- **Test-quality gate: PASS** — exact sourced expectations, full branch/error coverage, no
  green-washing, honest full-suite green.
- No defects logged.
