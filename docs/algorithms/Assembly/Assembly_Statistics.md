# Assembly Statistics (N50 / L50 / Nx / Lx / auN)

| Field | Value |
|-------|-------|
| Algorithm Group | Assembly |
| Test Unit ID | ASSEMBLY-STATS-001 |
| Related Projects | Seqeron.Genomics.Chromosome |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Assembly statistics summarize the contiguity and composition of a genome assembly from the lengths of its contigs or scaffolds. The headline metric is **N50** — the length of the shortest contig in the smallest set of largest contigs whose combined length reaches at least 50% of the assembly — together with its companion count **L50**, the generalizations **Nx/Lx** (e.g. N90/L90), the area-under-the-Nx-curve metric **auN**, and gap (N-run) detection. The computation is exact and deterministic: contigs are sorted longest-first and a cumulative length threshold is applied [1][2][3].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A *de novo* assembly is a set of contigs (gap-free runs) or scaffolds (contigs joined by N-gaps). Because a single global "assembly length" hides how fragmented the result is, contiguity metrics report a weighted length distribution. N50 is the de facto standard contiguity statistic [1].

### 2.2 Core Model

Let the contig lengths be `L₁ ≥ L₂ ≥ … ≥ Lₙ` (sorted descending) and `T = Σ Lᵢ` the total assembly length.

- **Nx / Lx.** For a threshold `x ∈ [0,100]`, take the smallest prefix `L₁…Lₖ` whose cumulative length `Σᵢ₌₁ᵏ Lᵢ ≥ (x/100)·T`. Then `Nx = Lₖ` (the length of the last/shortest contig in that prefix) and `Lx = k` (the count) [1][2]. The cumulative threshold is **inclusive** ("at least x%") [1][3].
- **N50 / L50.** The special case `x = 50`: "the length of the smallest contig in the set that contains the fewest (largest) contigs whose combined length represents at least 50% of the assembly" [1]. N50 is a length; L50 is a count [2].
- **auN.** The area under the Nx curve, `auN = Σᵢ Lᵢ·(Lᵢ / Σⱼ Lⱼ) = Σᵢ Lᵢ² / Σⱼ Lⱼ` — a continuous, more stable contiguity measure than the discrete N50 [4][3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Nx is the length of the shortest contig in the largest-first prefix with cumulative length ≥ (x/100)·T | Definition in §2.2 [1] |
| INV-02 | Lx is the number of contigs in that prefix (a count, not a length) | Definition [2] |
| INV-03 | Nx is non-increasing and Lx non-decreasing in x; hence N90 ≤ N50 and L90 ≥ L50 | Larger x extends the prefix [2] |
| INV-04 | auN = Σ Lᵢ² / T | Formula in §2.2 [4][3] |
| INV-05 | The threshold boundary is inclusive: cumulative reaching exactly (x/100)·T selects that contig | "at least 50%" [1]; QUAST `s <= limit` [3] |
| INV-06 | `FindGaps` yields one entry per maximal run of `N`/`n` with 0-based inclusive `[Start,End]` and `Length = End − Start + 1` | Implementation contract |

### 2.5 Comparison with Related Methods

| Aspect | N50 | auN |
|--------|-----|-----|
| Nature | discrete (a single contig length) | continuous (integral of the Nx curve) |
| Stability | sensitive to small jumps near L50 | stable; joining any two contigs always increases it [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequences | `IEnumerable<(string Id, string Sequence)>` | required | Contigs/scaffolds | `CalculateStatistics`, `FindGaps` |
| sortedLengths | `IReadOnlyList<int>` | required | Lengths sorted **descending** | `CalculateNx` core overload |
| lengths | `IEnumerable<int>` | required | Lengths in any order (sorted internally) | `CalculateNx(lengths,threshold)`, `CalculateN50`, `CalculateAuN` |
| totalLength | `long` | required | Σ of lengths (assembly size) | core `CalculateNx` |
| threshold | `int` | required | Percentage x ∈ [0,100] (50 for N50) | Nx/Lx |
| minGapLength | `int` | 1 | Minimum N-run length reported | `FindGaps` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `AssemblyStatistics` | record struct | `N50, L50, N90, L90, LargestContig, SmallestContig, TotalLength, TotalLengthNoGaps, MeanLength, MedianLength, GcContent, TotalGaps, TotalGapLength, GapPercentage`, `TotalSequences` |
| `NxStatistics` | record struct | `Threshold, Nx, Lx, CumulativeLength` |
| `CalculateN50` | `int` | N50 length |
| `CalculateAuN` | `double` | auN = Σl²/Σl |
| `GapInfo` | record struct | `SequenceId, Start, End, Length, GapType` (0-based inclusive) |

### 3.3 Preconditions and Validation

Lengths/positions are 0-based; gap coordinates are inclusive `[Start, End]`. N-runs are matched case-insensitively (`N`/`n`). GC fraction is over non-N bases. Empty input is handled without throwing: `CalculateStatistics` returns an all-zero record, `CalculateNx` returns `Nx=Lx=0`, `CalculateAuN` returns `0` (see §6.1). No exceptions are raised for valid (possibly empty) collections.

## 4. Algorithm

### 4.1 High-Level Steps

1. Sort contig lengths descending; compute `T = Σ Lᵢ`.
2. Walk the sorted lengths, accumulating cumulative length and a count.
3. At the first contig where `cumulative ≥ (x/100)·T`, return its length as `Nx` and the running count as `Lx`.
4. `auN = Σ Lᵢ² / T`. `FindGaps` scans each sequence for maximal `N`/`n` runs.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The cumulative threshold is evaluated with **integer arithmetic** — `cumulative * 100 ≥ totalLength * threshold` — rather than `(long)(total * x / 100.0)`. This reproduces QUAST's exact stop test `total − cumulative ≤ total·(100−x)/100` (`s <= limit`) without floating-point rounding of the cutoff on odd totals [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateNx` (pre-sorted) | O(n) | O(1) | single pass over lengths |
| `CalculateNx(lengths,·)` / `CalculateN50` / `CalculateAuN` / `CalculateStatistics` | O(n log n) | O(n) | dominated by the descending sort |
| `FindGaps` | O(m) | O(g) | m = total bases, g = gaps found |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomeAssemblyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs)

- `GenomeAssemblyAnalyzer.CalculateStatistics(sequences)`: aggregate N50/L50/N90/L90/largest/smallest/totals/GC/gaps.
- `GenomeAssemblyAnalyzer.CalculateNx(sortedLengths, totalLength, threshold)`: core Nx/Lx.
- `GenomeAssemblyAnalyzer.CalculateNx(lengths, threshold)`: sorts + totals, delegates to the core overload.
- `GenomeAssemblyAnalyzer.CalculateN50(lengths)`: `CalculateNx(lengths, 50).Nx`.
- `GenomeAssemblyAnalyzer.CalculateAuN(lengths)`: Σl²/Σl.
- `GenomeAssemblyAnalyzer.FindGaps(sequences, minGapLength)`: maximal N-run detection.

### 5.2 Current Behavior

N50/N90/L50/L90 in `CalculateStatistics` are produced by the same `CalculateNx` core as the standalone methods, guaranteeing identical results. **Search reuse:** N/A — this is not a substring/pattern-matching unit; `FindGaps` is a single linear scan for `N`/`n` runs (a suffix tree offers no benefit for a one-pass character classification), so the repository suffix tree was not used.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- N50/L50 as "the smallest contig … whose combined length represents at least 50% of the assembly" with inclusive cumulative threshold [1].
- General Nx/Lx for arbitrary x [1][2].
- auN = Σl²/Σl [4][3].
- Integer-exact inclusive cutoff matching QUAST `NG50_and_LG50` [3].

**Intentionally simplified:**

- `AssemblyStatistics.MedianLength`: reports the upper median (`lengths[count/2]` of the descending list); the N50 literature does not define an assembly "median contig length", so this auxiliary field is not part of the cited contract. **Consequence:** for even contig counts the median may be the upper of the two central values.

**Not implemented:**

- NG50/LG50 (genome-size-relative variants); **users should rely on:** providing an explicit reference length and the Nx core — no current dedicated NG50 method (out of scope for this unit) [2].

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Empty input returns zeros instead of `None` (QUAST) | Assumption | No valid N50 exists for an empty assembly; zeros avoid an exception | accepted | §6.1; non-correctness-affecting (ASSUMPTION 1) |
| 2 | `MedianLength` upper-median convention | Assumption | Auxiliary field, outside the N50 contract | accepted | §5.3 "Intentionally simplified" (ASSUMPTION 2) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty contig set | all-zero `AssemblyStatistics`; `Nx=Lx=0`; `auN=0` | No defined N50 for an empty assembly [3] |
| Single contig | `Nx = that length`, `Lx = 1` | It is the whole assembly |
| Cumulative exactly = x% of total | that contig is selected (inclusive) | INV-05 [1][3] |
| Leading / trailing N-run | reported as a gap at `[0,k-1]` / `[…,len-1]` | INV-06 |
| All-N contig | `TotalGapLength = length`, `GapPercentage = 100` | Whole sequence is gap |

### 6.2 Limitations

Operates on lengths and base composition only; it does not assess base-level accuracy, misassembly, or completeness (BUSCO-like analysis lives in separate methods). NG50/LG50 against an external genome size is not provided as a dedicated entry point.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var lengths = new[] { 80, 70, 50, 40, 30, 20 }; // total 290
int n50 = GenomeAssemblyAnalyzer.CalculateN50(lengths);              // 70
var nx  = GenomeAssemblyAnalyzer.CalculateNx(lengths, 50);           // Nx=70, Lx=2
double aun = GenomeAssemblyAnalyzer.CalculateAuN(lengths);           // 16700/290 ≈ 57.5862
```

**Numerical walk-through (Assembly A [2]):** lengths {80,70,50,40,30,20}, T=290, 50%·T=145. Cumulative: 80 (<145), 80+70=150 (≥145) → stop. N50 = 70, L50 = 2.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomeAssemblyAnalyzer_AssemblyStatistics_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/GenomeAssemblyAnalyzer_AssemblyStatistics_Tests.cs) — covers `INV-01`…`INV-06`
- Evidence: [ASSEMBLY-STATS-001-Evidence.md](../../../docs/Evidence/ASSEMBLY-STATS-001-Evidence.md)
- TestSpec: [ASSEMBLY-STATS-001.md](../../../tests/TestSpecs/ASSEMBLY-STATS-001.md)

## 8. References

1. Miller JR, Koren S, Sutton G. 2010. Assembly algorithms for next-generation sequencing data. *Genomics* 95(6):315-327. https://pmc.ncbi.nlm.nih.gov/articles/PMC2874646/ (DOI: 10.1016/j.ygeno.2010.03.001)
2. Wikipedia contributors. N50, L50, and related statistics. https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics (accessed 2026-06-13)
3. Gurevich A, Saveliev V, Vyahhi N, Tesler G. 2013. QUAST — `quast_libs/N50.py`. https://raw.githubusercontent.com/ablab/quast/master/quast_libs/N50.py (accessed 2026-06-13; DOI: 10.1093/bioinformatics/btt086)
4. Li H. 2020. auN: a new metric to measure assembly contiguity. https://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity (accessed 2026-06-13)
