# Breakpoint Detection from Split Reads

| Field | Value |
|-------|-------|
| Algorithm Group | StructuralVar |
| Test Unit ID | SV-BREAKPOINT-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Breakpoint detection localizes the junctions of structural variants from *split reads* — reads whose alignment is interrupted by an SV breakpoint so that one segment maps to the reference and the other (soft-clipped) segment does not [2]. The aligned/clipped junction of each read is a single-base estimate of the breakpoint [3]. The algorithm clusters split reads whose junctions agree within a small positional tolerance and reports a breakpoint for each cluster that meets a minimum read-support threshold [3][4]. It is a heuristic, signature-then-cluster method (not probabilistic); given soft-clip-derived junctions it is deterministic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

When a sequenced fragment crosses an SV breakpoint, an aligner maps the portion matching the reference and soft-clips the rest. In SAM/BAM, a soft clip is the CIGAR `S` operation: it consumes the read but not the reference, and the clipped bases are retained in `SEQ` [1]. A read with one anchored end and one imprecisely-mapping (clipped) end is the *split-read* signature of a breakpoint [2]. Split-read methods achieve single-base-pair resolution because the breakpoint is the exact aligned/clipped junction [2][3].

### 2.2 Core Model

- **Per-read breakpoint coordinate.** The breakpoint is "the marginal point between a clipped sequence and matched sequence" [3] — the reference coordinate at the boundary of the aligned segment. Because POS is "the 1-based leftmost mapping position of the first CIGAR operation that consumes a reference base" and `S` does not consume the reference [1], a left clip places the junction at the read's anchored start and a right clip at the anchored end plus the aligned (M/=/X/D/N-consumed) length [1][3].
- **Clustering.** Junctions are sorted by chromosome then position and grouped while consecutive junctions stay within a tolerance window — ClipCrop clusters breakpoints "within 5-base differences" [3]; mapping imprecision otherwise spreads one true breakpoint across nearby positions [3]. A breakpoint coordinate is per-chromosome (SAM POS is contig-local [1]).
- **Support / calling rule.** "A putative breakpoint is defined when there is at least x soft-clipped reads beginning at position y" [4]; the support of a breakpoint is the number of clipped reads in its cluster [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported breakpoint has `SupportingReads ≥ minSupport`. | Clusters below `minSupport` are dropped before emission [4]. |
| INV-02 | Two reads share a breakpoint only if same chromosome AND junction gap ≤ `clusterTolerance`. | Linear scan splits a cluster on chromosome change or a gap above tolerance [1][3]. |
| INV-03 | A reported position lies within `[min, max]` member junction (≤ tolerance of every member). | Position = rounded mean of member junctions, all inside one tolerance window [3] (ASM-01). |
| INV-04 | `SupportingReads` equals the cluster size. | Support is defined as the count of clipped reads at the position [4]. |
| INV-05 | Number of breakpoints ≤ number of input split reads. | Each read joins exactly one cluster (a partition). |

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The reported cluster coordinate is the rounded mean of member junctions. The sources fix the per-read junction [3] and the tolerance window [3] but not the summary statistic. | A different statistic (mode/min) shifts the reported coordinate by at most the tolerance; cluster membership and support are unchanged. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `splitReads` | `IEnumerable<SplitRead>` | required | Split reads; `SupplementaryPosition` is the junction coordinate, `Chromosome` the contig. | non-null; coordinates are reference positions |
| `clusterTolerance` | `int` | 5 | Max junction-position gap to keep adjacent reads in one breakpoint [3]. | ≥ 0 (bases) |
| `minSupport` | `int` | 2 | Minimum supporting split reads to report a breakpoint [4]. | ≥ 1 |
| `chromosome` / `regionStart` / `regionEnd` (RefineBreakpoint) | `string` / `int` / `int` | required | Candidate region (inclusive bounds) to refine. | `regionStart ≤ regionEnd` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Breakpoint.Position1` / `Position2` | `int` | Consensus breakpoint coordinate (rounded mean of member junctions). |
| `Breakpoint.Chromosome1` / `Chromosome2` | `string` | Contig of the breakpoint (same for both ends; intra-contig clustering). |
| `Breakpoint.SupportingReads` | `int` | Number of split reads in the cluster [4]. |
| `RefineBreakpoint` return | `int?` | Consensus junction (mode; tie → rounded mean) inside the region, or `null` if unsupported. |

### 3.3 Preconditions and Validation

Null `splitReads` throws `ArgumentNullException` (both methods). Empty input yields an empty breakpoint sequence; an unsupported refine region yields `null`. Coordinates are taken as-is (the caller supplies reference positions); region bounds are inclusive. Comparison of chromosome names is ordinal/case-sensitive, consistent with sibling methods in the class.

## 4. Algorithm

### 4.1 High-Level Steps

1. Take each read's junction coordinate (`SupplementaryPosition`) as the single-base breakpoint estimate [3].
2. Sort reads by chromosome, then by junction position [3].
3. Scan linearly; start a new cluster whenever the chromosome changes or the gap to the previous junction exceeds `clusterTolerance` [1][3].
4. Emit a `Breakpoint` for each cluster with `≥ minSupport` reads; its position is the rounded mean of member junctions and its support is the cluster size [4] (ASM-01).
5. `RefineBreakpoint`: among reads whose junction lies in `[regionStart, regionEnd]` on `chromosome`, return the modal junction (tie → rounded mean), else `null` [3][4].

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Parameter | Value | Source |
|-----------|-------|--------|
| Junction cluster tolerance | 5 bases (default) | ClipCrop "clustered within 5-base differences" [3] |
| Minimum supporting reads | 2 (default) | SoftSearch configurable minimum [4]; matches sibling BreakDancer -r = 2 |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindBreakpoints` | O(n log n) | O(n) | Dominated by the sort of n split reads; the cluster scan is O(n). |
| `RefineBreakpoint` | O(n) | O(k) | Single filter + frequency count over n reads; k = junctions in region. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [StructuralVariantAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs)

- `StructuralVariantAnalyzer.FindBreakpoints(splitReads, clusterTolerance, minSupport)`: clusters split-read junctions per chromosome and reports breakpoints meeting minimum support (canonical).
- `StructuralVariantAnalyzer.RefineBreakpoint(chromosome, regionStart, regionEnd, splitReads)`: returns the consensus junction inside a candidate region.

### 5.2 Current Behavior

Clustering is performed on the junction coordinate (`SupplementaryPosition`), which differs from the sibling `ClusterSplitReads`, which groups on the anchored `PrimaryPosition`. The junction is the breakpoint per ClipCrop [3], so it is the correct clustering key for single-base breakpoint localization. The reported coordinate uses the rounded mean of member junctions (ASM-01). Validation is eager (`ArgumentNullException` thrown before iteration) while clustering is deferred via an iterator, matching the sibling `DetectSVs` pattern in the same class.

**Search-reuse note:** the repository suffix tree was evaluated and is not applicable. This unit performs no substring/pattern search; it sorts and groups numeric junction coordinates within a tolerance window, for which a sort + linear scan is the appropriate algorithm.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Per-read breakpoint = aligned/clipped junction; single-base resolution [2][3].
- Positional clustering within a tolerance window [3].
- Support = number of clipped reads at the position; minimum-support gate [4].
- Per-chromosome (contig-local) breakpoint coordinates [1].

**Intentionally simplified:**

- Cluster summary coordinate is the rounded mean of member junctions; **consequence:** the reported position may differ from the modal junction by at most `clusterTolerance` (ASM-01). Membership and support are unaffected.
- Strand/orientation (`Strand1='+'`, `Strand2='-'`) is set to a fixed convention rather than derived from left/right clip side; **consequence:** the `Breakpoint` strands do not encode the L/R clip orientation that ClipCrop [3] uses to distinguish breakpoint sides.

**Not implemented:**

- Junction-sequence assembly / microhomology at the breakpoint; **users should rely on:** `AssembleBreakpointSequence` and `FindMicrohomology` in the same class.
- SV-type classification from paired breakpoints (deletion/duplication/inversion/translocation); **users should rely on:** `DetectSVs` / `ClassifySV` (SV-DETECT-001) for paired-end-signature typing.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Cluster coordinate = rounded mean of junctions | Assumption | Reported position within ≤ tolerance of true junction | accepted | ASM-01; sub-tolerance, no effect on membership/support |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty input | Empty result | No junctions to cluster. |
| Null input | `ArgumentNullException` | Input-validation contract (sibling methods). |
| Cluster below `minSupport` | Not reported | Minimum-support gate [4]. |
| Junctions > tolerance apart | Separate clusters | Cluster only within tolerance [3]. |
| Same junction, different chromosomes | Separate breakpoints | POS is contig-local [1]. |
| Refine region with no member junction | `null` | No support to form a consensus [4]. |

### 6.2 Limitations

Uses only the junction coordinate stored on `SplitRead`; it does not re-derive the junction from a CIGAR string, does not pair two breakpoints into a typed SV, and does not encode left/right clip orientation in the output strands. Support thresholds and tolerance are heuristic defaults [3][4] and should be tuned to library/coverage. The method is intended for short-read soft-clip signatures, not long-read or assembly-based breakpoint refinement.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var reads = new[]
{
    new StructuralVariantAnalyzer.SplitRead("r1", "chr1", 4500, 5000, 40, "ACGT..."),
    new StructuralVariantAnalyzer.SplitRead("r2", "chr1", 4502, 5002, 38, "ACGT..."),
    new StructuralVariantAnalyzer.SplitRead("r3", "chr1", 4504, 5004, 42, "ACGT..."),
};

// Junctions {5000, 5002, 5004} agree within tolerance 5 → one breakpoint at 5002, support 3.
var breakpoints = StructuralVariantAnalyzer.FindBreakpoints(reads, clusterTolerance: 5, minSupport: 2);
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [StructuralVariantAnalyzer_FindBreakpoints_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/StructuralVariantAnalyzer_FindBreakpoints_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [SV-BREAKPOINT-001-Evidence.md](../../../docs/Evidence/SV-BREAKPOINT-001-Evidence.md)
- Related algorithms: [SV_Detection](../StructuralVar/SV_Detection.md)

## 8. References

1. The SAM/BAM Format Specification Working Group. 2024. Sequence Alignment/Map Format Specification (SAMv1) — CIGAR operations and POS. samtools/hts-specs. https://samtools.github.io/hts-specs/SAMv1.pdf
2. Tattini L, D'Aurizio R, Magi A. 2015. Detection of Genomic Structural Variants from Next-Generation Sequencing Data. Front Bioeng Biotechnol 3:92. https://doi.org/10.3389/fbioe.2015.00092
3. Suzuki S, Yasuda T, Shiraishi Y, Miyano S, Nagasaki M. 2011. ClipCrop: a tool for detecting structural variations with single-base resolution using soft-clipping information. BMC Bioinformatics 12(Suppl 14):S7. https://doi.org/10.1186/1471-2105-12-S14-S7
4. Hart SN, Sarangi V, Moore R, et al. 2013. SoftSearch: Integration of Multiple Sequence Features to Identify Breakpoints of Structural Variations. PLoS ONE 8(12):e83356. https://doi.org/10.1371/journal.pone.0083356
