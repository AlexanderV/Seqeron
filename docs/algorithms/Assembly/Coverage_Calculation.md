# Coverage (Depth) Calculation

| Field | Value |
|-------|-------|
| Algorithm Group | Assembly |
| Test Unit ID | ASSEMBLY-COVER-001 |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Production |
| Last Reviewed | 2026-11-25 |

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
| ASM-01 | Each read is placed once, end-to-end and ungapped (Bowtie 2 default `--end-to-end`: no soft clipping), at its minimum-Hamming window; optionally on the reverse strand too; optionally skipped when multi-mapped. | Reads that need gaps or clipping may be placed elsewhere or not at all; the counting arithmetic is unchanged. |

Under the Lander-Waterman uniform-random placement model the per-base depth is Poisson with rate C, so P(position uncovered) = e^−C and breadth ≈ 1 − e^−C [4][5]. This is a statistical expectation about the *inputs*, not a rule the depth array obeys — the returned array is exact for whatever placements occur.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Output length = reference length. | One depth value per reference position [3]. |
| INV-02 | depth[i] ≥ 0 for all i. | Depth is a count of reads [2][3]. |
| INV-03 | Σ depth[i] = Σ over placed reads of overlap length with reference. | Sum of per-base depths = total bases mapped [2]. |
| INV-04 | A read placed at p of length L increments exactly [p, p+L) ⊆ [0, G). | Per-base depth = reads whose aligned span covers the position [3][6]; end-to-end placement never overhangs the reference. |
| INV-05 | A read that fails to place (best match < minOverlap) adds 0 everywhere. | An unaligned read contributes no depth [2][3]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| reference | string | required | Reference sequence; depth reported per position. | Non-null; 0-based positions. |
| reads | IReadOnlyList&lt;string&gt; | required | Reads to map. | Non-null; a null element throws `ArgumentException`. |
| minOverlap | int | 20 | Minimum matching characters (length − Hamming distance) to place a read. | ≥ 1, else `ArgumentOutOfRangeException`. Usability default, not a biological constant. |
| includeReverseComplement | bool | true | Also align each read's reverse complement (both strands; Bowtie 2/BWA default) [7]; `false` = forward strand only (Bowtie 2 `--norc`). | 3-argument overload = `true`; 5-argument overload `(reference, reads, minOverlap, includeReverseComplement, excludeMultiMapped = false)` takes it explicitly. |
| excludeMultiMapped | bool | false | Skip reads whose best score is reached at > 1 distinct locus (Bowtie 2 MAPQ ≤ 1 → `samtools depth -Q 2`) [6][8]. | |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | int[] | Array of length `reference.Length`; element i = number of placed reads covering reference position i (0-based). |

### 3.3 Preconditions and Validation

Null `reference` or `reads` throws `ArgumentNullException`; a null read throws `ArgumentException`; `minOverlap < 1` throws `ArgumentOutOfRangeException` (a placement with zero matching bases is not an alignment). An empty reference yields an empty array; an empty reads list yields an all-zero array of length `reference.Length`. Matching is case-insensitive (both sides upper-cased). Placement is ungapped; a read longer than the reference cannot be placed and contributes 0. Reads are not normalized for alphabet (T↔U); characters are compared literally after case folding.

## 4. Algorithm

### 4.1 High-Level Steps

1. Allocate a depth array of length `reference.Length` (all zeros).
2. For each read, find its minimum-Hamming end-to-end window with the canonical `ApproximateMatcher.FindBestMatch` (and, if requested, that of its reverse complement); keep it if it has ≥ `minOverlap` matching characters; if `excludeMultiMapped`, drop the read when `ApproximateMatcher.FindWithMismatches` finds the best distance at more than one distinct start (both strands pooled).
3. If placed at position p, increment depth[i] for i in [p, p + read.Length).
4. Return the depth array.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Placement = minimum Hamming distance over windows [p, p+L), p ∈ 0 … G − L (`ApproximateMatcher.FindBestMatch`: case-insensitive, leftmost on ties). With both strands the reverse-complement placement (Core `TryGetReverseComplement`) replaces the forward one only when its distance is strictly smaller (Bowtie 2 picks pseudo-randomly among equal hits; the forward-first, leftmost choice makes the result deterministic). The read is placed iff L − distance ≥ `minOverlap`. Loci for the multi-mapping test are distinct start positions, so a reverse-palindromic read hitting both strands at the same span is unique. Ungapped placements have no deletions, so `samtools depth` with or without `-J` gives the same depth.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateCoverage | O(r · n · m) | O(n) | r = #reads, n = reference length, m = read length (placement is the dominant cost; depth increment is O(m) per read). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAssembler.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

- `SequenceAssembler.CalculateCoverage(string, IReadOnlyList<string>, int)`: returns the per-base depth array.
- `SequenceAssembler.CalculateCoverage(string, IReadOnlyList<string>, int, bool, bool)`: both strands / multi-mapping options.
- Placement delegates to `ApproximateMatcher.FindBestMatch` / `FindWithMismatches` (PAT-APPROX-001); the former private `FindBestAlignment` scan was removed.

### 5.2 Current Behavior

Returns the exact per-base depth array; average depth and breadth are left to the caller to derive from it (sum / length, and covered-count / length respectively), preserving the full coverage distribution. The repository suffix tree (`SuffixTree.FindAllOccurrences`) was evaluated for read placement: it finds *exact* occurrences only, whereas `CalculateCoverage` places each read at its single *best* (most-matching, possibly imperfect) position with a `minOverlap` floor — a scoring-based placement the suffix tree does not provide. The canonical approximate matcher in the same project is therefore used.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Per-base depth[i] = number of placed reads covering position i [3].
- Average depth and breadth are recoverable as Σ depth / G and (covered count) / G [2][3].

**Intentionally simplified:**

- Read placement: ungapped end-to-end best-Hamming placement (Bowtie 2 `--end-to-end` without gaps, mismatch cost 1, `--score-min` replaced by `minOverlap`); **consequence:** reads requiring gaps, splicing or soft clipping (e.g. overhanging a contig end) are placed differently or not at all (the counting rule is unchanged).

**Not implemented:**

- Graded MAPQ (only the equal-best multi-mapping filter is modelled), base-quality filtering (`samtools depth -q`), duplicate flags; **users should rely on:** dedicated aligners/`samtools depth` for production BAM-derived depth.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Read placement model (end-to-end ungapped best Hamming) | Assumption | Determines where reads map; depth arithmetic is unaffected | accepted | ASM-01; cross-checked against `samtools depth` on 2014 cases (Evidence). |
| 2 | 3-argument overload searches both strands (default flipped from forward-only, audit A6-1) | Resolved | Matches Bowtie 2/BWA defaults; forward-only via `includeReverseComplement: false` | resolved | `cover_oracle.py` (pysam `samtools depth -a`) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null reference or reads | `ArgumentNullException` | Input validation (sibling convention). |
| Empty reads list | All-zero array of length `reference.Length` | No aligned reads → 0 depth everywhere [2][3]. |
| Read below minOverlap | Contributes 0 | Unaligned read adds no depth [2][3] (INV-05). |
| Read longer than reference | Contributes 0 (cannot be placed) | End-to-end alignment requires the whole read to align (Bowtie 2 `--end-to-end`) (INV-04). |
| minOverlap < 1 | `ArgumentOutOfRangeException` | Zero matching bases is not an alignment. |
| Reverse-strand read | Placed by default; unplaced with `includeReverseComplement: false` | Bowtie 2 `--nofw/--norc`: both strands by default [7]. |
| Equal-best hits at several loci | Counted once (forward-first, leftmost); skipped with `excludeMultiMapped` | Bowtie 2 one reported alignment; MAPQ ≤ 1 filtered by `samtools depth -Q 2` [6][8]. |
| Lowercase read vs uppercase reference | Maps; depth counted | Case-insensitive comparison. |

### 6.2 Limitations

The depth array reflects the repository's ungapped best-match placement, not a full aligner; it does not model indels, soft-clipping or graded mapping quality (multi-mapping reads: one placement, or excluded on request). It is suited to exact/near-exact short reads against a reference, and to deriving average depth and breadth, rather than replacing BAM-based depth tools.

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
6. samtools-depth(1) manual. https://raw.githubusercontent.com/samtools/samtools/develop/doc/samtools-depth.1
7. Bowtie 2 manual (`--end-to-end`, `--nofw/--norc`, default reporting, MAPQ). https://raw.githubusercontent.com/BenLangmead/bowtie2/master/MANUAL.markdown
8. Bowtie 2 `unique.h` (`BowtieMapq2`). https://raw.githubusercontent.com/BenLangmead/bowtie2/master/unique.h
