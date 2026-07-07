# Coverage (Depth) Calculation

| Field | Value |
|-------|-------|
| Algorithm Group | Assembly |
| Test Unit ID | ASSEMBLY-COVER-001 |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Coverage (sequencing depth) calculation reports, for every position of a reference sequence, how many reads cover that position. It is the per-base depth array that underlies downstream summaries such as average depth and breadth of coverage. NGS coverage "describes the average number of reads that align to, or 'cover,' known reference bases" [1]; the per-position form retains the full distribution so that uneven coverage is visible rather than collapsed into a single mean. The computation is exact arithmetic over read placements: it is not heuristic or probabilistic once the read positions are fixed.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

In shotgun sequencing a genome is read as many short fragments (reads). After reads are aligned to a reference, the depth at a position is the number of reads spanning it [3]. Depth quantifies how much evidence supports each base and is the foundation for variant calling and assembly quality assessment. Two summary statistics are derived from the per-base depth: **average depth** = sum of per-base depths / genome size [2], and **breadth of coverage** = fraction of reference bases covered by at least one read [2][3].

### 2.2 Core Model

For a reference of length G and a set of placed reads, let each placed read r start at position pᵣ and have length Lᵣ. The per-base depth is

> depth[i] = #{ r : pᵣ ≤ i < pᵣ + Lᵣ }, for 0 ≤ i < G

i.e. the number of reads whose span (half-open interval) contains position i [3]. Summary forms [1][2][3]:

- Average depth (coverage) C = Σᵢ depth[i] / G, equivalently the Lander-Waterman C = LN/G with L = read length, N = number of reads, G = genome length [1].
- Breadth = #{ i : depth[i] ≥ 1 } / G [2][3].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Reads are placed at their (single) best ungapped match position before counting depth. | If reads multi-map or align with gaps, the placement — and therefore which positions are incremented — may differ; the counting arithmetic is unchanged. |

Under the Lander-Waterman uniform-random placement model the per-base depth is Poisson with rate C, so P(position uncovered) = e^−C and breadth ≈ 1 − e^−C [4][5]. This is a statistical expectation about the *inputs*, not a rule the depth array obeys — the returned array is exact for whatever placements occur.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Output length = reference length. | One depth value per reference position [3]. |
| INV-02 | depth[i] ≥ 0 for all i. | Depth is a count of reads [2][3]. |
| INV-03 | Σ depth[i] = Σ over placed reads of overlap length with reference. | Sum of per-base depths = total bases mapped [2]. |
| INV-04 | A read placed at p of length L increments exactly [p, min(p+L, G)). | Per-base depth = reads spanning the position, clipped at the reference end [3]. |
| INV-05 | A read that fails to place (best match < minOverlap) adds 0 everywhere. | An unaligned read contributes no depth [2][3]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| reference | string | required | Reference sequence; depth reported per position. | Non-null; 0-based positions. |
| reads | IReadOnlyList&lt;string&gt; | required | Reads to map. | Non-null; individual entries non-null. |
| minOverlap | int | 20 | Minimum matching characters to place a read. | Usability default, not a biological constant. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | int[] | Array of length `reference.Length`; element i = number of placed reads covering reference position i (0-based). |

### 3.3 Preconditions and Validation

Null `reference` or `reads` throws `ArgumentNullException`. An empty reference yields an empty array; an empty reads list yields an all-zero array of length `reference.Length`. Matching is case-insensitive (both sides upper-cased). Placement is ungapped; a read longer than the reference cannot be placed and contributes 0. Reads are not normalized for alphabet (T↔U); characters are compared literally after case folding.

## 4. Algorithm

### 4.1 High-Level Steps

1. Allocate a depth array of length `reference.Length` (all zeros).
2. For each read, find its best ungapped match position via a sliding window requiring ≥ `minOverlap` matching characters.
3. If placed at position p, increment depth[i] for i in [p, min(p + read.Length, reference.Length)).
4. Return the depth array.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Read placement (`FindBestAlignment`) scans positions 0 … (referenceLength − readLength), counting case-insensitive character matches, and keeps the position with the most matches provided it is ≥ `minOverlap`; ties keep the leftmost. Because the scan only considers positions where the read fits entirely, a placed read never extends past the reference end (the `i < reference.Length` guard in the depth loop is defensive).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateCoverage | O(r · n · m) | O(n) | r = #reads, n = reference length, m = read length (placement is the dominant cost; depth increment is O(m) per read). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAssembler.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

- `SequenceAssembler.CalculateCoverage(string, IReadOnlyList<string>, int)`: returns the per-base depth array.
- `SequenceAssembler.FindBestAlignment(string, string, int)` (private): places a read at its best ungapped match.

### 5.2 Current Behavior

Returns the exact per-base depth array; average depth and breadth are left to the caller to derive from it (sum / length, and covered-count / length respectively), preserving the full coverage distribution. The repository suffix tree (`SuffixTree.FindAllOccurrences`) was evaluated for read placement: it finds *exact* occurrences only, whereas `CalculateCoverage` places each read at its single *best* (most-matching, possibly imperfect) position with a `minOverlap` floor — a scoring-based placement the suffix tree does not provide. The existing ungapped best-match scan is therefore used, consistent with the sibling assembly methods in the same class.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Per-base depth[i] = number of placed reads covering position i [3].
- Average depth and breadth are recoverable as Σ depth / G and (covered count) / G [2][3].

**Intentionally simplified:**

- Read placement: ungapped best-match scan; **consequence:** reads requiring gapped/spliced alignment may be placed at a position different from a full aligner, changing which positions are incremented (not the counting rule).

**Not implemented:**

- Quality- or mapping-quality-weighted depth, and multi-mapping read distribution; **users should rely on:** dedicated aligners/`samtools depth` for production BAM-derived depth.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Read placement model (best ungapped match) | Assumption | Determines where reads map; depth arithmetic is unaffected | accepted | ASM-01; tests use exact-match reads to isolate the source-defined counting rule. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null reference or reads | `ArgumentNullException` | Input validation (sibling convention). |
| Empty reads list | All-zero array of length `reference.Length` | No aligned reads → 0 depth everywhere [2][3]. |
| Read below minOverlap | Contributes 0 | Unaligned read adds no depth [2][3] (INV-05). |
| Read longer than reference | Contributes 0 (cannot be placed) | Best-match scan requires the read to fit (INV-04). |
| Lowercase read vs uppercase reference | Maps; depth counted | Case-insensitive comparison. |

### 6.2 Limitations

The depth array reflects the repository's ungapped best-match placement, not a full aligner; it does not model indels, soft-clipping, mapping quality, or multi-mapping. It is suited to exact/near-exact short reads against a reference, and to deriving average depth and breadth, rather than replacing BAM-based depth tools.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var reference = "ACGTTGCAAT";              // length 10
var reads = new[] { "ACGTT", "TTGCA", "GCAAT" }; // unique 5-mers placed at 0, 3, 5
int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);
// depth == [1,1,1,2,2,2,2,2,1,1]
double avg = depth.Average();              // 1.5  (sum 15 / 10)
double breadth = depth.Count(d => d >= 1) / (double)depth.Length; // 1.0
```

**Numerical walk-through:** reads cover [0,5), [3,8), [5,10). Overlaps: positions 3,4 covered by reads 1&2; positions 5,6,7 by reads 2&3; the rest once. Σ depth = 5+5+5 = 15, average = 15/10 = 1.5, breadth = 10/10 = 1.0 [2][3].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceAssembler_CalculateCoverage_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAssembler_CalculateCoverage_Tests.cs) — covers `INV-01`–`INV-05`.
- Evidence: [ASSEMBLY-COVER-001-Evidence.md](../../../docs/Evidence/ASSEMBLY-COVER-001-Evidence.md)
- Related algorithms: [Overlap-Layout-Consensus](../Assembly/Overlap_Layout_Consensus.md)

## 8. References

1. Illumina, Inc. (n.d.). Sequencing Coverage for NGS Experiments. https://sapac.illumina.com/science/technology/next-generation-sequencing/plan-experiments/coverage.html
2. Cook, D.E. (n.d.). Calculate Depth and Breadth of Coverage From a bam File. https://www.danielecook.com/calculate-depth-and-breadth-of-coverage-from-a-bam-file/
3. Metagenomics Wiki. (n.d.). SAMtools: get breadth of coverage. https://www.metagenomics.wiki/tools/samtools/breadth-of-coverage
4. Daley, T. et al. (2020). Predicting the Number of Bases to Attain Sufficient Coverage in High-Throughput Sequencing Experiments. PMC7398442. https://pmc.ncbi.nlm.nih.gov/articles/PMC7398442/
5. Lander, E.S., Waterman, M.S. (1988). Genomic mapping by fingerprinting random clones: a mathematical analysis. Genomics 3(2):231-239. https://doi.org/10.1016/0888-7543(88)90007-9
