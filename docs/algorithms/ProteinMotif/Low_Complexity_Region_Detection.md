# Low-Complexity Region Detection (SEG)

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinMotif |
| Test Unit ID | PROTMOTIF-LC-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Low-complexity regions (LCRs) are protein segments with strongly biased residue composition — homopolymers, short-period repeats, and compositionally skewed mosaics — which break the statistical assumptions of sequence-comparison tools. The SEG algorithm of Wootton & Federhen (1993) detects them by sliding a fixed-length window along the sequence and measuring local compositional complexity in bits per residue [1][2]. This implementation is specification-driven: it reports a region wherever a window's complexity falls to or below a trigger cutoff, extending it over neighbouring windows that stay below an extension cutoff [2][6].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

About half of SWISS-PROT entries contain at least one nonrandomly low-complexity segment [1]. SEG measures local complexity "based on an information measure of the complexity state vector, which reflects residue composition appearing on a sliding window, with no regard of the patterns or periodicity of sequence repetitiveness" [6].

### 2.2 Core Model

For a window of length `L` containing residue counts `n₁…n_k` (Σnᵢ = L), let pᵢ = nᵢ/L. The SEG local complexity, measured in bits per residue, is the Shannon entropy of the composition [2][3][5]:

> K = −Σᵢ pᵢ·log₂(pᵢ)

K ranges from 0 (a homopolymer window) to log₂(N) for an alphabet of size N; for the 20 amino acids the maximum is log₂(20) ≈ 4.322 bits/residue [2]. This is the bits/residue form of "equation (3) of Wootton & Federhen (1993)" referenced by the SEG specification [2]; the NCBI reference implementation computes it in `s_Entropy` (Shannon entropy of the count vector normalized to bits) [3], and the SeqComplex `ce` routine encodes the identical formula `ce -= (n/tot)·log₂(n/tot)` [4].

Detection is a two-pass procedure [3][6]:

1. **Trigger (pass 1):** every window with complexity ≤ K1 (the trigger/locut cutoff) marks a raw low-complexity segment.
2. **Extension (pass 2):** a triggered segment is extended over adjacent windows whose complexity is ≤ K2 (the extension/hicut cutoff). Only K2 > K1 has any extending effect [2].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ K ≤ log₂(20) for amino-acid windows | Shannon entropy of an N-symbol distribution is bounded by log₂N; N=20 [2] |
| INV-02 | A homopolymer window has K = 0 | a single symbol has p=1, −1·log₂1 = 0 [4][5] |
| INV-03 | A reported region contains at least one window with K ≤ K1 and every window with K ≤ K2 | two-pass trigger/extension rule [2][3] |
| INV-04 | Region boundaries are 0-based inclusive and lie within the sequence | windowed scan; span = [firstStart, lastStart + W − 1] |
| INV-05 | Output is deterministic, ordered by Start | single left-to-right scan |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| proteinSequence | string | required | protein sequence | case-insensitive; null/empty → no regions |
| windowSize | int | 12 | sliding-window length W [2][3] | must be > 0 (else `ArgumentOutOfRangeException`) |
| triggerComplexity | double | 2.2 | trigger cutoff K1, bits/residue [2][3] | typically in [0, log₂20] |
| extensionComplexity | double | 2.5 | extension cutoff K2, bits/residue [2][3] | K2 > K1 to extend [2] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Start | int | 0-based inclusive start residue of the region |
| End | int | 0-based inclusive end residue of the region |
| Complexity | double | minimum window complexity (bits/residue) observed inside the region |

### 3.3 Preconditions and Validation

`windowSize ≤ 0` throws `ArgumentOutOfRangeException`. Null, empty, or sequences shorter than `windowSize` yield an empty result (no complete trigger window exists [2]). Input is upper-cased before counting (case-insensitive); residue identity is taken literally, so any non-standard characters are simply additional symbols in the composition.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `windowSize`; return empty for sequences shorter than the window.
2. For every window start position, compute K = −Σ pᵢ·log₂ pᵢ over the window composition.
3. Scan left to right; group maximal runs of windows with K ≤ K2.
4. Emit a run as a region only if at least one of its windows has K ≤ K1 (triggered); the region span is `[runStart, runEnd + W − 1]` and its reported complexity is the run's minimum window complexity.

### 4.2 Decision Rules, Scoring, Reference Tables

| Parameter | Value | Source |
|-----------|-------|--------|
| Window W | 12 | NCBI SEG man page; `blast_seg.c` `kSegWindow = 12` [2][3] |
| Trigger K1 (locut) | 2.2 bits | man page; `kSegLocut = 2.2` [2][3] |
| Extension K2 (hicut) | 2.5 bits | man page; `kSegHicut = 2.5` [2][3] |
| Max complexity (20 aa) | log₂(20) = 4.322 | man page [2] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindLowComplexityRegions | O(n·W) | O(n) | n = sequence length; each of n−W+1 windows costs O(W) for composition counting; per-window complexities stored in an array |

A suffix tree was **not** used: the operation is a single linear scan computing a per-window numeric score, not an exact substring/occurrence search, so the repository suffix tree does not apply (see "Reuse existing infrastructure").

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProteinMotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)

- `ProteinMotifFinder.FindLowComplexityRegions(string, int, double, double)`: public SEG entry point; returns `(Start, End, Complexity)` regions.
- `ProteinMotifFinder.CalculateSegComplexity(ReadOnlySpan<char>)`: internal per-window Shannon-entropy complexity in bits/residue.

### 5.2 Current Behavior

Per-window complexity is computed with a `stackalloc` count array indexed by character, so the routine is allocation-light. Each distinct residue is tallied once (the count slot is cleared after consumption) to avoid double-counting duplicates. Regions are emitted via iterator (`yield return`) ordered by Start.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Sliding-window complexity K = −Σ pᵢ·log₂ pᵢ in bits/residue, range 0…log₂20 [2][3][5].
- Two-pass trigger (K1) / extension (K2) detection [2][3][6].
- Default parameters W=12, K1=2.2, K2=2.5 [2][3].

**Intentionally simplified:**

- Complexity uses the Shannon-entropy bits/residue form of equation (3); **consequence:** the optional exact pass-2 local optimization that minimizes the multinomial probability P0 of a raw segment (using the log-factorial / permutation form, NCBI `s_LnPerm`/`lnfact[]`) is not performed, so region boundaries match the window-run span rather than the P0-minimized sub-segment.

**Not implemented:**

- The amino-acid-frequency-weighted P0 significance ranking of segments; **users should rely on:** the reported per-region complexity for relative ranking, or an external SEG/NCBI tool for P0 statistics.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Prior "dominant single-AA frequency ≥ 0.4" rule | Deviation | invented, non-source-traceable threshold and output fields | fixed | replaced by SEG complexity; return type changed `(…,DominantAa,Frequency)` → `(…,Complexity)` |
| 2 | Bits/residue entropy form of eq.(3) | Assumption | which interconvertible form of eq.(3) is used | accepted | man-page units ("bits", max log₂20) select this form (ASM in Evidence) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty sequence | empty result | nothing to window |
| length < windowSize | empty result | no complete trigger window [2] |
| windowSize ≤ 0 | `ArgumentOutOfRangeException` | invalid window |
| homopolymer tract | single region, complexity 0 | INV-02 |
| fully diverse sequence | empty result | every window K > K2 [2] |
| lowercase input | detected (upper-cased) | case-insensitive |

### 6.2 Limitations

SEG ignores periodicity/patterns — a perfect short-period repeat with balanced composition (e.g. `ABABAB…`) has high entropy and is not flagged, by design [6]. The implementation does not compute P0-based statistical significance, and reports the window-run span rather than a P0-minimized optimal segment (see 5.3).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Detect poly-Q low-complexity tract with SEG defaults (W=12, K1=2.2, K2=2.5).
var regions = ProteinMotifFinder
    .FindLowComplexityRegions("MKLPRDST" + new string('Q', 20) + "MKLPRDST")
    .ToList();
// -> one region spanning the Q tract; Complexity == 0 for the homopolymer core.
```

**Numerical walk-through:** for a 12-residue window of 11 A and 1 B, p_A = 11/12, p_B = 1/12, so
K = −[(11/12)·log₂(11/12) + (1/12)·log₂(1/12)] = 0.413817 bits/residue — well below K1 = 2.2, so the window triggers a low-complexity segment.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ProteinMotifFinder_FindLowComplexityRegions_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindLowComplexityRegions_Tests.cs) — covers `INV-01`…`INV-05`
- Evidence: [PROTMOTIF-LC-001-Evidence.md](../../../docs/Evidence/PROTMOTIF-LC-001-Evidence.md)
- Related algorithms: [Coiled_Coil_Prediction](Coiled_Coil_Prediction.md)

## 8. References

1. Wootton JC, Federhen S. 1993. Statistics of local complexity in amino acid sequences and sequence databases. Computers & Chemistry 17(2):149–163. https://doi.org/10.1016/0097-8485(93)85006-X
2. NCBI. SEG program man page (`ncbi-seg`). https://manpages.ubuntu.com/manpages/focal/man1/ncbi-seg.1.html
3. NCBI C++ Toolkit. `blast_seg.c` (s_Entropy, kSegWindow/kSegLocut/kSegHicut). https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/blast__seg_8c.html
4. Caballero J. SeqComplex (`SeqComplex.pm`, subs `cwf`, `ce`, `log_k`). https://github.com/caballero/SeqComplex
5. Mier P, et al. 2006. Novel/Shannon-entropy complexity (−Σ pᵢ log pᵢ). Bioinformatics 22(24):2980. https://academic.oup.com/bioinformatics/article/22/24/2980/208627
6. Pei J, Grishin NV. 2005. A new algorithm for detecting low-complexity regions in protein sequences. Bioinformatics 21(2):160–166. https://academic.oup.com/bioinformatics/article/21/2/160/187330
