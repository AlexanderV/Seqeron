# Transmembrane Helix Prediction (Kyte-Doolittle Hydropathy)

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinMotif |
| Test Unit ID | PROTMOTIF-TM-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Predicts transmembrane α-helices in a protein sequence using the Kyte-Doolittle hydropathy method [1]. The mean hydropathy is computed over a sliding window of residues; contiguous windows whose mean exceeds a threshold are reported as candidate membrane-spanning segments. The method is a deterministic heuristic (a hydropathy filter), not a probabilistic model; it is the classical, fully-specified alternative to hidden-Markov membrane predictors and is appropriate for fast, transparent screening of single-pass and multi-pass membrane proteins.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Membrane-spanning α-helices are built from stretches of hydrophobic residues that integrate into the hydrophobic interior of the lipid bilayer. A single α-helix needs ≈18–21 residues to cross the ≈30 Å bilayer (each residue rises ≈1.5 Å along the helix axis) [4]. Kyte and Doolittle assigned each amino acid a hydropathy value and showed that averaging these values over a sliding window reveals such hydrophobic stretches as peaks [1].

### 2.2 Core Model

Each residue *r* has a hydropathy value *h(r)* on the Kyte-Doolittle scale [1][3]:

```
I 4.5  V 4.2  L 3.8  F 2.8  C 2.5  M 1.9  A 1.8  G -0.4  T -0.7  S -0.8
W -0.9 Y -1.3 P -1.6 H -3.2 E -3.5 Q -3.5 D -3.5 N -3.5  K -3.9  R -4.5
```

For a window of width *w* starting at position *i*, the profile value is the **arithmetic mean** of the per-residue values in the window [2][4]:

```
P(i) = (1/w) · Σ_{j=i}^{i+w-1} h(s[j])
```

A residue position is part of a transmembrane segment when its window mean exceeds the threshold *T*. Kyte and Doolittle found *w = 19* and *T = 1.6* optimal for identifying membrane-spanning helices [2]. A maximal run of consecutive windows with *P(i) ≥ T* is reported as one segment.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported segment's peak score ≥ *T*. | A segment is, by construction, a run of windows each with `P(i) ≥ T`; the peak is the max of those. |
| INV-02 | `0 ≤ Start ≤ End ≤ length−1`. | `Start` is the first above-threshold window's first residue; `End` is the last covered residue (`lastProfileIndex + windowSize − 1`), clamped to `length−1`. |
| INV-03 | A uniform run of one residue (length ≥ *w*) gives `P(i) = h(r)` everywhere. | The mean of *w* identical values is that value [1]. |
| INV-04 | Null / empty / shorter-than-window input yields no segments. | No window can be formed below *w* residues. |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | Kyte-Doolittle hydropathy | TMHMM (HMM) |
|--------|---------------------------|-------------|
| Model type | Deterministic sliding-window filter | Probabilistic hidden-Markov model |
| Parameters | Fully published (scale, window 19, threshold 1.6) [1][2] | Learned emission/transition probabilities |
| Topology (in/out) | Not predicted | Predicted |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `proteinSequence` | `string` | required | Amino-acid sequence, one-letter codes | case-insensitive; non-standard residues (X/B/Z/*) excluded from the mean |
| `windowSize` | `int` | 19 | Sliding-window width [2] | must be > 0 and ≤ sequence length to yield output |
| `threshold` | `double` | 1.6 | Mean-hydropathy cutoff [2] | any real value |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Start` | `int` | 0-based inclusive index of the first residue of the segment |
| `End` | `int` | 0-based inclusive index of the last residue covered by an above-threshold window (`lastProfileIndex + windowSize − 1`, clamped to `length−1`) |
| `Score` | `double` | Peak (maximum) window mean within the segment |

### 3.3 Preconditions and Validation

Input is uppercased internally (case-insensitive). 0-based, inclusive coordinates. Null, empty, non-positive `windowSize`, or sequences shorter than `windowSize` yield an empty sequence (no exception). Non-standard residues contribute nothing to the window mean (the mean is taken over the residues that have scale values).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate input; uppercase the sequence.
2. Compute the sliding-window hydropathy profile `P(i)` (arithmetic mean over each window).
3. Scan the profile; open a region at the first `P(i) ≥ T`, track the peak, close it when `P(i) < T`.
4. Report a closed region as a segment if its residue span ≥ the minimum helix length; map the closing window index to a residue `End`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

- Kyte-Doolittle scale (20 values) — §2.2, from [1][3].
- Default window = 19, threshold = 1.6 [2].
- Minimum reported helix span = the window width (19), reflecting the ≈18–21 residues needed to span the bilayer [4]. For the canonical window this filter is always satisfied by any region containing at least one above-threshold window.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Profile + scan | O(n·w) | O(n) | n = sequence length, w = window; scan is O(n) over the profile |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProteinMotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)

- `ProteinMotifFinder.PredictTransmembraneHelices(string, int, double)`: returns `(int Start, int End, double Score)` tuples for predicted segments.
- `ProteinMotifFinder.CalculateHydropathyProfile(...)` (private): the sliding-window arithmetic mean.

### 5.2 Current Behavior

The method computes the full profile then performs a single linear scan to delimit threshold-crossing runs. No search/matching over a text is involved (it is a numeric window scan, not occurrence enumeration), so the repository suffix tree is not applicable. The window mean divides by the count of residues that have scale values, so a window containing non-standard residues is averaged over its standard residues only.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Kyte-Doolittle hydropathy scale values (all 20 residues) [1][3].
- Sliding-window arithmetic-mean profile [2][4].
- Default window 19 and threshold 1.6 [2].

**Intentionally simplified:**

- Segment-to-residue boundary mapping: the source defines the window run, not exact residue boundaries; **consequence:** the reported `End` is the last residue actually covered by an above-threshold window, `lastProfileIndex + windowSize − 1` (clamped) — i.e. the union of all passing windows' residue spans (§5.4). It does not change which windows pass detection.

**Not implemented:**

- Membrane topology (in/out orientation) and loop prediction; **users should rely on:** a dedicated topology predictor (e.g. an HMM-based tool) for orientation.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `End = lastProfileIndex + windowSize − 1` (clamped) | Assumption | Affects reported segment end coordinate only | accepted | Last residue covered by an above-threshold window (union of passing windows' spans); detection (which windows pass) is unchanged |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty | empty result | no window possible |
| length < window | empty result | no window possible |
| `windowSize ≤ 0` | empty result | guarded invalid input |
| all-hydrophilic | no segments | every window mean < threshold |
| non-standard residue in window | excluded from the mean | scale defined for 20 standard residues [4] |

### 6.2 Limitations

A hydropathy filter detects hydrophobic stretches; it does not distinguish true membrane helices from hydrophobic cores of soluble proteins, does not predict topology, and may merge or split closely-spaced helices depending on window/threshold. Results are sensitive to the chosen window and threshold; defaults follow Kyte & Doolittle [2].

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
var segments = ProteinMotifFinder.PredictTransmembraneHelices(
    new string('D', 10) + new string('L', 20) + new string('D', 10));
// → one segment: (Start: 5, End: 34, Score: 3.8)
```

**Numerical walk-through:** For `D×10 L×20 D×10` (window 19, threshold 1.6): the profile has 22 points; window means rise above 1.6 first at profile index 5 and last at index 16; any all-Leu window has mean 3.8 (peak). The last passing window (start 16) covers residues 16..34, so the single run maps to residues (5, 34).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ProteinMotifFinder_PredictTransmembraneHelices_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/ProteinMotifFinder_PredictTransmembraneHelices_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [PROTMOTIF-TM-001-Evidence.md](../../../docs/Evidence/PROTMOTIF-TM-001-Evidence.md)

## 8. References

1. Kyte J, Doolittle RF. 1982. A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology* 157(1):105-132. https://doi.org/10.1016/0022-2836(82)90515-0
2. Davidson College, Department of Biology — Genomics Project. Kyte-Doolittle background. https://gcat.davidson.edu/DGPB/kd/kyte-doolittle-background.htm
3. QIAGEN CLC Genomics Workbench manual — Hydrophobicity scales. https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html
4. Biopython — Bio.SeqUtils.ProtParam (`protein_scale`, `gravy`). https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParam.py
