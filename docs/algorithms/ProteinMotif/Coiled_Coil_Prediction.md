# Coiled-Coil Prediction (Heptad-Repeat a/d Hydrophobic-Core Detection)

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinMotif |
| Test Unit ID | PROTMOTIF-CC-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Coiled coils are bundles of α-helices wound around one another, encoded by a seven-residue **heptad repeat** in which the hydrophobic core positions drive helix association [1][2]. This algorithm scans a protein sequence and predicts coiled-coil regions by measuring how strongly the hydrophobic-core positions (**a** and **d**) of the heptad are occupied by hydrophobic residues. It is a heuristic, sequence-only predictor: for each window it computes the fraction of a/d positions filled by core residues, maximised over the seven possible heptad registers, and reports contiguous high-scoring spans. It intentionally does **not** use the COILS position-specific scoring matrix (whose weights were not obtainable from authoritative sources); it uses the fully-specified a/d occupancy model instead [2].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The coiled coil is one of the most common protein oligomerization motifs. Its sequence signature is a repeating block of seven residues, conventionally labelled **a b c d e f g**, that recurs as (abcdefg)n [1][2]. When folded into an α-helix (3.6 residues/turn), residues a and d fall on the same face, forming a hydrophobic stripe whose burial against another helix ("knobs into holes") stabilises the assembly [1].

### 2.2 Core Model

For a sequence `S` of length `n` (uppercased), a window of length `W` starting at index `i`, and a register `r ∈ {0,…,6}`:

- The heptad position of residue index `k` in register `r` is `p(k,r) = (k − r) mod 7`.
- The **core positions** are `a = 0` and `d = 3` [1][2].
- A residue is a **hydrophobic-core residue** if it is one of `{I, L, V}` — "a and d … often being occupied by isoleucine, leucine, or valine" [3].
- Window occupancy in register `r`: `occ(i,r) = (#core positions in window with residue ∈ {I,L,V}) / (#core positions in window)`.
- Window score: `score(i) = max over r of occ(i,r)` — the register is unknown a priori, so all seven frames are tried and the best is taken [2].

A window is coiled-coil if `score(i) ≥ threshold`. Contiguous coiled windows `[i₀ … i₁]` map to the residue region `[i₀, i₁ + W − 1]`; the region is reported (with `Score` = peak `score` in the run) only if its length is at least the minimum (3 heptads = 21 residues) [1].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every `Score` ∈ [0, 1]. | Score is a count ratio `hydrophobic / core` [3]. |
| INV-02 | Each region spans ≥ 21 residues (3 heptads). | Enforced by `MinCoiledCoilRegion`; multi-heptad requirement [1]. |
| INV-03 | `0 ≤ Start ≤ End ≤ n − 1`; regions are non-overlapping, increasing in `Start`. | Region construction from a single forward scan. |
| INV-04 | Sequences shorter than `windowSize` produce no regions. | No full window exists [2]. |
| INV-05 | A region exists only if some covering window scores ≥ threshold (max over 7 registers). | Definition of `score(i)` [2]. |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | This algorithm (a/d occupancy) | COILS (Lupas 1991) |
|--------|-------------------------------|---------------------|
| Scoring | Fraction of a/d positions in {I,L,V} | Geometric mean of per-position residue frequencies vs reference set [2] |
| Parameter source | Residue set {I,L,V} + heptad geometry | 21×20 position-specific weight matrix |
| Reproducibility from public spec | Fully specified | Requires the published PSSM (not retrieved here) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `proteinSequence` | `string` | required | Amino-acid sequence | Case-insensitive; null/empty → no result |
| `windowSize` | `int` | 28 | Sliding-window length (4 heptads) [2] | ≥ 1; sequences shorter than this yield nothing |
| `threshold` | `double` | 0.5 | Minimum a/d hydrophobic-core occupancy [1] | Typically in [0,1] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Start` | `int` | 0-based inclusive start residue index of the predicted region |
| `End` | `int` | 0-based inclusive end residue index |
| `Score` | `double` | Peak a/d occupancy fraction in the region, ∈ [0,1] |

### 3.3 Preconditions and Validation

Input is uppercased (case-insensitive). `null`, empty, or any sequence shorter than `windowSize` returns an empty enumeration (no exception). Indexing is 0-based; `End` is inclusive. The accepted alphabet is the 20 amino acids; non-{I,L,V} residues simply do not count toward occupancy.

## 4. Algorithm

### 4.1 High-Level Steps

1. Reject null/empty/sub-window input.
2. For each window start `i`, compute `score(i)` = max over the 7 heptad registers of the a/d hydrophobic-core occupancy fraction.
3. Scan the score profile; group contiguous windows with `score ≥ threshold` into runs.
4. Map each run `[i₀,i₁]` to residue region `[i₀, i₁ + W − 1]`; keep it if length ≥ 21; emit `(Start, End, peakScore)`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

- Hydrophobic-core residue set: `{I, L, V}` [3].
- Core heptad positions: `a = 0`, `d = 3` [1][2].
- Window = 28, registers = 7, minimum region = 21 [1][2].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `PredictCoiledCoils` | O(n · W · 7) | O(n) | n = length, W = window. With constants W=28, 7 registers, this is O(n). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProteinMotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)

- `ProteinMotifFinder.PredictCoiledCoils(sequence, windowSize, threshold)`: public entry; returns coiled-coil regions.
- `BestHeptadOccupancy(...)` (private): max a/d occupancy over the 7 registers for one window.
- `BuildRegion(...)` (private): maps a window run to a residue region with the min-length filter.

### 5.2 Current Behavior

The score profile is computed once into an array, then scanned in a single forward pass; regions are yielded lazily. No substring search is performed (the algorithm is a per-position numeric scan over a fixed-size window), so the repository suffix tree is **not** applicable — there is no exact-match occurrence enumeration here, only windowed scoring.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Heptad repeat (abcdefg)n with hydrophobic core at positions a and d [1][2].
- Hydrophobic-core residue set {I, L, V} [3].
- Seven heptad registers tried per window, best taken; window of 28 residues [2].
- Minimum reported coiled coil of multiple (≥3) heptads [1].

**Intentionally simplified:**

- Scoring is the a/d hydrophobic-core occupancy fraction rather than a position-specific frequency model; **consequence:** scores are occupancy fractions in [0,1], not COILS P-scores, and e/g electrostatic-edge residues do not contribute.

**Not implemented:**

- The COILS position-specific scoring matrix (21×20 residue-frequency weights) [2]; **users should rely on:** the dedicated COILS/Paircoil2 tools for PSSM-based probabilities — the matrix weights were not retrievable from authoritative sources and are not fabricated here.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Core residue set limited to {I,L,V} | Assumption | A,M,F-rich coiled coils may score lower | accepted | Exactly the set named in [3]; avoids untraceable constants |
| 2 | COILS PSSM omitted | Deviation | No probabilistic P-score output | accepted | PSSM weights not retrievable; see 5.3 "Not implemented" |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `null` / empty | empty result | validation |
| length < windowSize | empty result | no full window [2] (INV-04) |
| no {I,L,V} residues | empty result | occupancy 0 < threshold [3] |
| off-frame coiled coil | found via 7-register max | register unknown a priori [2] (INV-05) |
| lowercase input | recognised | sequence is uppercased |

### 6.2 Limitations

Sequence-only heuristic: it detects the a/d hydrophobic periodicity but does not model helix propensity, e/g electrostatics, oligomerization state, or register breaks/stutters. It will not reproduce COILS/Paircoil probabilities and may over-predict generic hydrophobic-periodic sequences.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
// A perfect 5-heptad repeat with L at every a and d position.
var regions = ProteinMotifFinder.PredictCoiledCoils(string.Concat(Enumerable.Repeat("LAALAAA", 5)));
// → one region (Start=0, End=34, Score=1.0)
```

**Numerical walk-through:** For `"LAALAAA"×5` (35 aa), register 0 places `a` at indices 0,7,14,… and `d` at 3,10,17,… — every such index holds `L ∈ {I,L,V}`, so occupancy = 1.0 for all 8 windows (W=28). The single run [0..7] maps to residues [0, 7+27] = [0,34], length 35 ≥ 21 → `(0, 34, 1.0)`.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ProteinMotifFinder_PredictCoiledCoils_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_PredictCoiledCoils_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [PROTMOTIF-CC-001-Evidence.md](../../../docs/Evidence/PROTMOTIF-CC-001-Evidence.md)
- Related algorithms: [Transmembrane Helix Prediction](../ProteinMotif/Transmembrane_Helix_Prediction.md)

## 8. References

1. Mason JM, Arndt KM. 2004. Coiled coil domains: stability, specificity, and biological implications. ChemBioChem 5(2):170–176. https://doi.org/10.1002/cbic.200300781
2. Lupas A, Van Dyke M, Stock J. 1991. Predicting coiled coils from protein sequences. Science 252(5009):1162–1164. https://doi.org/10.1126/science.252.5009.1162
3. Wikipedia. Coiled coil. https://en.wikipedia.org/wiki/Coiled_coil (accessed 2026-06-14)
4. Wikipedia. Heptad repeat. https://en.wikipedia.org/wiki/Heptad_repeat (accessed 2026-06-14)
