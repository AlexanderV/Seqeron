# Hydrophobicity Analysis (Kyte-Doolittle GRAVY and Hydropathy Profile)

| Field | Value |
|-------|-------|
| Algorithm Group | Statistics |
| Test Unit ID | SEQ-HYDRO-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Quantifies the hydropathic character of a protein from its amino-acid sequence using the Kyte-Doolittle hydropathy scale [1]. Two outputs are produced: the GRAVY index (Grand Average of Hydropathy), a single scalar summarizing whether a protein is overall hydrophobic (positive) or hydrophilic (negative) [4], and a sliding-window hydropathy profile that locates hydrophobic stretches such as transmembrane segments [5]. The computation is exact and specification-driven: every residue maps to a fixed published constant and the outputs are deterministic averages.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Each amino-acid side chain has a characteristic affinity for water. Kyte & Doolittle combined water-vapor transfer free energies with interior/exterior side-chain distributions to assign each of the 20 standard residues a single "hydropathy" index ranging from +4.5 (most hydrophobic, Ile) to −4.5 (most hydrophilic, Arg) [1]. Averaging these indices reveals hydrophobic protein regions, which tend to be buried or membrane-spanning [5].

### 2.2 Core Model

Let `kd(r)` be the Kyte-Doolittle index of residue `r` [1][2]. For a sequence `S = s_1 … s_n`:

- **GRAVY** = (Σ kd(s_i)) / n — the sum of hydropathy values divided by the number of residues [4]. Biopython's `gravy()` implements this as `total_gravy / length` [3].
- **Hydropathy profile** for window size `W`: for each window start `i = 1 … n−W+1`, the value is Σ w_j·kd(s_{i+j}) / Σ w_j. With the default edge weight 1 all `w_j = 1`, i.e. the unweighted mean (1/W)·Σ kd of the original Kyte-Doolittle method [1]. With edge weight `e < 1` (ExPASy ProtScale "linear" weight-variation model, Biopython `protein_scale(kd, W, e)`) the centre residue has weight 1, the two ends weight `e`, and `w_j = e + j·2(1−e)/(W−1)` for the j-th position from either end [3]. This yields exactly `n − W + 1` values; value `k` (0-based) is the score of the window's central residue `k + (W−1)/2` for odd `W`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | GRAVY of a single recognized residue equals that residue's `kd` value | GRAVY = sum/length with length 1 [3][4] |
| INV-02 | Profile produces exactly `n − W + 1` values for `W ≤ n`, else none | Window loop `range(n−W+1)` [3] |
| INV-03 | Each profile value is the unweighted mean of its window | Default equal weighting (`edge=1.0`) [3] |
| INV-04 | GRAVY is case-insensitive | Scale defined on uppercase; input is uppercased before lookup |
| INV-05 | Empty/null sequence → GRAVY 0 and empty profile | No residues to average (contract) |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | Kyte-Doolittle GRAVY | Per-window profile |
|--------|----------------------|--------------------|
| Output | one scalar for the whole protein | a series, one value per window |
| Typical use | bulk hydrophobicity / solubility screening [4] | locating transmembrane / surface regions (W=19 peaks > 1.6; W=9 surface) [5] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| proteinSequence | string | required | one-letter amino-acid sequence | case-insensitive; non-standard residues skipped/0 |
| windowSize | int | 9 | sliding-window length for the profile | ≥ 1 (else `ArgumentOutOfRangeException`); profile empty when windowSize > length |
| edgeWeight | double | 1.0 | relative weight of the window edges (linear model) | in [0, 1] (else `ArgumentOutOfRangeException`); `< 1` requires an odd window (else `ArgumentException`) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| GRAVY | double | sum of `kd` values of recognized residues / recognized-residue count; 0 if none |
| profile | IEnumerable\<double\> | `n − W + 1` unweighted window means (lazy) |

### 3.3 Preconditions and Validation

Input is case-insensitive (uppercased before lookup). Only the 20 standard residues are defined in the scale [1]; any other character (ambiguity codes B/Z/X, gaps, stop) is skipped by GRAVY (not added to the count) and contributes 0 to a profile window. Null or empty input returns GRAVY 0 and an empty profile without throwing. No exceptions are raised for any character input.

## 4. Algorithm

### 4.1 High-Level Steps

1. Uppercase the input.
2. **GRAVY:** sum `kd` over recognized residues, count them, return sum/count (0 if count is 0).
3. **Profile:** validate `windowSize ≥ 1` and `edgeWeight ∈ [0,1]` eagerly; if `windowSize > length`, yield nothing; otherwise build the position weights and slide the window across all `n − W + 1` positions, yielding Σ w·kd / Σ w (= sum / W for edge weight 1).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

Kyte-Doolittle hydropathy index (the materially behavior-defining lookup table) [1][2]:

| A 1.8 | R −4.5 | N −3.5 | D −3.5 | C 2.5 | Q −3.5 | E −3.5 | G −0.4 | H −3.2 | I 4.5 |
|-------|--------|--------|--------|-------|--------|--------|--------|--------|-------|
| **L 3.8** | **K −3.9** | **M 1.9** | **F 2.8** | **P −1.6** | **S −0.8** | **T −0.7** | **W −0.9** | **Y −1.3** | **V 4.2** |

Recommended windows: 9 (surface regions), 19 (transmembrane, peaks > 1.6) [5].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| GRAVY | O(n) | O(1) | single pass, constant scale lookup |
| Profile | O(n·W) | O(1) | recomputes each window sum; W is the window size |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceStatistics.CalculateHydrophobicity(string)`: GRAVY index over recognized residues.
- `SequenceStatistics.CalculateHydrophobicityProfile(string, int windowSize = 9, double edgeWeight = 1.0)`: lazy sliding-window hydropathy profile (argument validation is eager).

### 5.2 Current Behavior

The scale is a static `Dictionary<char,double>` matching Biopython `kd` exactly [2]. GRAVY divides by the count of *recognized* residues (not raw string length); the profile recomputes each window's sum directly (O(n·W)) and yields lazily via `yield return`. This is not a search/matching unit, so the repository suffix tree is not applicable (no occurrence enumeration or pattern lookup is performed).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Kyte-Doolittle `kd` values for all 20 standard residues [1][2].
- GRAVY = Σ kd / residue count [3][4].
- Profile = `n − W + 1` unweighted window means with equal per-position weight [3].

**Intentionally simplified:**

- (none).

**Not implemented:**

- ProtScale's "exponential" weight-variation model and scale normalisation (only the linear model, as in Biopython, is implemented).
- Alternative scales (Hopp-Woods, Eisenberg); **users should rely on:** no current alternative in this class — only Kyte-Doolittle is provided.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Non-standard residues skipped (GRAVY divides by recognized count; profile adds 0 and keeps the position's weight in the divisor) | Deviation | Biopython `gravy()` raises `KeyError`; `protein_scale` adds 0 for an unknown centre residue (same as here) but also drops the symmetric partner of an off-centre unknown (not done here) | accepted | Sources define only the 20 standard residues [1][4]; canonical values unchanged. Tracked in Evidence Assumptions. |
| 2 | Even window with edge weight 1 → plain mean over W residues | Deviation | Biopython counts residue W/2 twice and divides by W+1 (e.g. "AV", W=2: 3.4 vs 3.0 here) | accepted | ProtScale accepts odd windows only; Biopython's value is an artefact of its odd-window loop. Weighted (edge < 1) even windows are rejected. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty sequence | GRAVY 0; empty profile | nothing to average (INV-05) |
| windowSize > length | empty profile | `n − W + 1 ≤ 0` (INV-02) [3] |
| windowSize < 1 | `ArgumentOutOfRangeException` | no window; previously W=0 yielded NaN (0/0) |
| lowercase input | same GRAVY as uppercase | input uppercased (INV-04) |
| non-standard residue (e.g. X) | skipped in GRAVY; 0 in profile window | undefined in scale [1] |

### 6.2 Limitations

Only the Kyte-Doolittle scale is supported. The profile supports the linear edge-weight model only. Biological interpretation thresholds (e.g. transmembrane peaks > 1.6 at W=19) are the caller's responsibility — the method returns raw averages, not classifications.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
double gravy = SequenceStatistics.CalculateHydrophobicity("FLIV");      // 3.825
var profile = SequenceStatistics.CalculateHydrophobicityProfile("FLIV", 3).ToList();
// profile == [ (2.8+3.8+4.5)/3, (3.8+4.5+4.2)/3 ] == [ 3.7, 4.1666666667 ]
```

**Numerical walk-through:** "RKDE" → (−4.5 + −3.9 + −3.5 + −3.5)/4 = −15.4/4 = −3.85 (hydrophilic) [1][4].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_CalculateHydrophobicity_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceStatistics_CalculateHydrophobicity_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [SEQ-HYDRO-001-Evidence.md](../../../docs/Evidence/SEQ-HYDRO-001-Evidence.md)
- Related algorithms: [Molecular_Weight_Calculation](Molecular_Weight_Calculation.md), [Isoelectric_Point](Isoelectric_Point.md)

## 8. References

1. Kyte, J., Doolittle, R.F. 1982. A simple method for displaying the hydropathic character of a protein. J Mol Biol 157(1):105–132. https://doi.org/10.1016/0022-2836(82)90515-0
2. Biopython. Bio/SeqUtils/ProtParamData.py (kd scale), master branch. https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/ProtParamData.py
3. Biopython. Bio/SeqUtils/ProtParam.py (gravy, protein_scale), master branch. https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/ProtParam.py
4. Expasy. ProtParam documentation (GRAVY). https://web.expasy.org/protparam/protparam-doc.html
5. GCAT (Davidson College). Kyte-Doolittle background. https://gcat.davidson.edu/DGPB/kd/kyte-doolittle-background.htm
