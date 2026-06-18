# Low-Complexity Region Detection (SEG)

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinPred |
| Test Unit ID | DISORDER-LC-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Low-complexity regions (LCRs) in proteins are segments of strongly biased residue composition — homopolymers, short-period repeats, and compositionally skewed mosaics — that are frequent in intrinsically disordered proteins and break the statistical assumptions of sequence-comparison tools [1]. This unit detects LCRs with the SEG algorithm of Wootton & Federhen (1993): a fixed-length window slides over the sequence and the *local compositional complexity* of each window is measured as Shannon entropy in bits per residue [1][3][4]. It is a deterministic, specification-driven heuristic: it reports a segment wherever a window's complexity falls to or below a trigger cutoff (K1), extending the segment over neighbouring residues that keep complexity at or below an extension cutoff (K2) [4].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

SEG measures local complexity from "an information measure of the complexity state vector, which reflects residue composition appearing on a sliding window, with no regard of the patterns or periodicity of sequence repetitiveness" [3][4]. The detector ignores the physicochemical identity of the residues; only their counts within the window matter.

### 2.2 Core Model

For a window of length `L`, let `nᵢ` be the count of residue type `i` and `pᵢ = nᵢ / L`. The window complexity is the Shannon entropy of the composition, in bits per residue [1][3]:

> H = − Σᵢ pᵢ · log₂(pᵢ)

This is exactly what the NCBI BLAST reference implementation computes in `s_Entropy`: it sums over the composition state vector, normalizes counts by the window length, and converts to base-2 via `NCBIMATH_LN2` [3]. The maximum value for the 20-letter amino-acid alphabet is `log₂(20) ≈ 4.322` bits/residue; a window of a single residue type has `H = 0` [4].

The detector applies three parameters [3][4]: trigger window length `W`, trigger complexity `K1`, extension complexity `K2`. Stage 1 marks every length-`W` window with `H ≤ K1`; stage 2 extends each triggered segment while window/segment complexity stays `≤ K2`; overlapping or adjacent segments are merged into contigs [4].

### 2.3 Modeling Assumptions (Optional)

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Composition entropy alone captures "low complexity" | periodic but compositionally rich repeats (high entropy) are not flagged — by design [3] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Each segment satisfies 0 ≤ Start ≤ End < n (0-based inclusive) | spans are built from in-range window indices |
| INV-02 | Reported segments are non-overlapping and non-adjacent | merge step combines spans where `start ≤ lastEnd + 1` [4] |
| INV-03 | A window of `L` identical residues has H = 0 ≤ K1 (always triggers) | single-type composition ⇒ `p=1`, `−1·log₂1 = 0` [3] |
| INV-04 | A window of `W` distinct residues has H = log₂(W) > K2 for default W=12 | log₂12 ≈ 3.585 > 2.5 [4] |
| INV-05 | Every reported segment has length ≥ `minLength` | minimum-length post-filter |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | Protein sequence | case-insensitive; non-letter chars excluded from composition counts |
| triggerWindow | int | 12 | Trigger window length W | ≥ 1; from NCBI `kSegWindow` [3] |
| triggerThreshold | double | 2.2 | Trigger complexity K1 (bits/residue) | from NCBI `kSegLocut` [3] |
| extensionThreshold | double | 2.5 | Extension complexity K2 (bits/residue) | from NCBI `kSegHicut` [3] |
| minLength | int | 1 | Minimum reported segment length | ≥ 1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Start | int | 0-based inclusive start index of the segment |
| End | int | 0-based inclusive end index of the segment |
| Type | string | Convenience composition label (`"X-rich"` or `"X/Y-rich"`); see §5.4 |

### 3.3 Preconditions and Validation

Input is upper-cased (case-insensitive). A `null` sequence raises `ArgumentNullException` on enumeration. A sequence shorter than `triggerWindow` yields no segments (no full trigger window exists) [4]. Coordinates are 0-based inclusive. Only ASCII letters A–Z contribute to composition counts; the entropy denominator is the full window length.

## 4. Algorithm

### 4.1 High-Level Steps

1. Upper-case the sequence; if shorter than `W`, return empty.
2. **Trigger (stage 1):** for each length-`W` window, compute `H`; if `H ≤ K1`, mark all positions in the window.
3. Collect contiguous marked positions into spans.
4. **Extension (stage 2):** grow each span left/right while the growing segment's `H ≤ K2`.
5. Merge overlapping/adjacent segments.
6. Drop segments shorter than `minLength`; classify and emit `(Start, End, Type)`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

| Parameter | Value | Source |
|-----------|-------|--------|
| Window W | 12 | NCBI `kSegWindow` [3] |
| Trigger K1 | 2.2 bits | NCBI `kSegLocut` [3]; GCG `-LOWcut` [4] |
| Extension K2 | 2.5 bits | NCBI `kSegHicut` [3]; GCG `-HIGhcut` [4] |
| Max complexity (20 aa) | log₂(20) ≈ 4.322 | GCG SEG help [4] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Detect LCRs | O(n·W) | O(n) | n windows, each entropy over W residues; recomputed per window |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)

- `DisorderPredictor.PredictLowComplexityRegions(sequence, triggerWindow, triggerThreshold, extensionThreshold, minLength)`: SEG trigger + extension detector.
- `CalculateShannonEntropy` (private): H = −Σ pᵢ·log₂(pᵢ) over a window (matches `s_Entropy` [3]).
- `ClassifyLowComplexityType` (private): convenience `Type` label.

### 5.2 Current Behavior

The detector is a single forward sliding-window scan, not a substring-occurrence search, so the repository suffix tree is **not applicable** (no exact-match enumeration is performed; complexity is a per-window numeric statistic). Entropy is recomputed for each window rather than incrementally maintained. Extension grows one residue at a time while the whole growing segment's entropy stays ≤ K2 (a contig-growth realization of the merge-extension-windows rule).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Complexity = Shannon entropy of window composition in bits/residue, H = −Σ pᵢ·log₂(pᵢ) [1][3].
- Default parameters W=12, K1=2.2, K2=2.5 from the NCBI reference implementation [3][4].
- Two-stage trigger (≤ K1) then extension (≤ K2) with segment merging [4].

**Intentionally simplified:**

- Extension: per-residue contig growth checking the whole segment entropy, instead of merging discrete length-W extension windows; **consequence:** identical boundaries for the homopolymer/biased inputs tested, but boundary placement can differ from NCBI SEG on mixed-complexity edges.
- Entropy recomputed per window (O(n·W)) rather than via Wootton & Federhen's exact compositional-complexity probability `P0` significance test; **consequence:** no per-segment probability/E-value is reported.

**Not implemented:**

- Wootton & Federhen significance measure `P0` (equation-based combinatorial probability) and overlap optimization; **users should rely on:** NCBI `seg` / BLAST for publication-grade masking.
- Region `Type` label is a repository extension, not part of SEG; **users should rely on:** the `Start`/`End` coordinates for spec-conformant segment positions.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `Type` label `"X-rich"`/`"X/Y-rich"` | Assumption | Cosmetic; not in SEG spec | accepted | Dominant residue if fraction > 0.5, else top two |
| 2 | Per-residue greedy extension | Deviation | Boundary placement on mixed edges | accepted | See §5.3; ASM-01 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Sequence shorter than W | empty | no full trigger window [4] |
| Empty string | empty | length 0 < W |
| Window of W distinct residues | not flagged | H = log₂W ≈ 3.585 > K2 [4] |
| Homopolymer ≥ W | one segment over whole run | H = 0 ≤ K1 [3] |
| `minLength` exceeds segment length | segment dropped | minimum-length post-filter |

### 6.2 Limitations

Detects compositionally biased regions only; periodic repeats made of many residue types are not flagged. No probability/E-value output. Per-window entropy recomputation is O(n·W). The `Type` label is heuristic and not part of the SEG specification.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
var regions = DisorderPredictor.PredictLowComplexityRegions(new string('Q', 26)).ToList();
// regions == [ (0, 25, "Q-rich") ] — every 12-window has entropy 0 ≤ 2.2
```

**Numerical walk-through:**

For `AAABBBCCCDDD` (L=12, four types ×3): pᵢ = 3/12 = 0.25 each; H = −4·(0.25·log₂0.25) = −4·(0.25·−2) = 2.0 bits. 2.0 ≤ K1 (2.2) ⇒ the window triggers; with a strict K1 = 0.5, 2.0 > 0.5 ⇒ no trigger.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [DisorderPredictor_LowComplexity_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_LowComplexity_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [DISORDER-LC-001-Evidence.md](../../../docs/Evidence/DISORDER-LC-001-Evidence.md)
- Related algorithms: [Disorder_Prediction](Disorder_Prediction.md)

## 8. References

1. Wootton J.C., Federhen S. 1993. Statistics of local complexity in amino acid sequences and sequence databases. Computers & Chemistry 17(2):149–163. https://doi.org/10.1016/0097-8485(93)85006-X
2. Wootton J.C., Federhen S. 1996. Analysis of compositionally biased regions in sequence databases. Methods in Enzymology 266:554–571. https://doi.org/10.1016/S0076-6879(96)66035-2
3. NCBI C++ Toolkit. `blast_seg.c` (SEG reference implementation; `kSegWindow=12`, `kSegLocut=2.2`, `kSegHicut=2.5`; `s_Entropy`). https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/blast__seg_8c.html
4. SEG program help (GCG/Weizmann mirror) and `ncbi-seg` manpage. https://bip.weizmann.ac.il/education/materials/gcg/seg.html ; https://manpages.ubuntu.com/manpages/focal/man1/ncbi-seg.1.html
