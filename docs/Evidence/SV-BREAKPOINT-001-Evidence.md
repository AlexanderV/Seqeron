# Evidence Artifact: SV-BREAKPOINT-001

**Test Unit ID:** SV-BREAKPOINT-001
**Algorithm:** Breakpoint Detection from Split (soft-clipped) Reads
**Date Collected:** 2026-06-13

---

## Online Sources

### Sequence Alignment/Map Format Specification (SAM/BAM spec, samtools/hts-specs) — CIGAR semantics

**URL:** https://raw.githubusercontent.com/samtools/hts-specs/master/SAMv1.tex (also mirrored at https://davetang.org/wiki/tiki-index.php?page=SAM)
**Accessed:** 2026-06-13
**Authority rank:** 2 (official file-format specification)
**Retrieved how:** WebSearch query `SAM format specification CIGAR consumes query consumes reference S soft clipping present in SEQ M I D N table`; WebFetch of the raw `SAMv1.tex` from the official `samtools/hts-specs` GitHub; cross-confirmed the operation descriptions with WebFetch of davetang.org SAM wiki.

**Key Extracted Points:**

1. **CIGAR consume table (verbatim):** `M alignment match — consumes query yes / reference yes`; `I insertion to the reference — query yes / reference no`; `D deletion from the reference — query no / reference yes`; `N skipped region from the reference — query no / reference yes`; `S soft clipping (clipped sequences present in SEQ) — query yes / reference no`; `H hard clipping (clipped sequences NOT present in SEQ) — query no / reference no`. → soft-clipped (S) bases consume the read but NOT the reference; aligned (M/=/X) and reference-only (D/N) operations advance the reference position.
2. **POS definition:** "1-based leftmost mapping POSition of the first CIGAR operation that 'consumes' a reference base." → the alignment's leftmost reference coordinate is where the first reference-consuming operation lands.
3. **SEQ length rule:** "Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ." → soft-clipped bases ARE present in SEQ (so the clipped sequence can be extracted), confirming S is recoverable from the read.
4. **Clip placement:** "S may only have H operations between them and the ends of the CIGAR string." → soft clips occur only at the read ends (left clip = leading S, right clip = trailing S).

### Tattini L, D'Aurizio R, Magi A (2015) — "Detection of Genomic Structural Variants from Next-Generation Sequencing Data", Front Bioeng Biotechnol 3:92

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4479793/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed review)
**Retrieved how:** WebSearch `breakpoint detection split reads clustering supplementary alignment SV review Tattini Magi structural variant single base resolution`; WebFetch of the open-access PMC article.

**Key Extracted Points:**

1. **Split-read signature:** a split read is one where "one end is anchored to the reference genome and the other end maps imprecisely owing to the presence of an underlying structural variant or indel breakpoint." → an anchored segment plus a clipped/displaced segment marks a breakpoint.
2. **Single-base resolution:** "SR methods allow for the detection of SVs with single base-pair resolution." → the breakpoint is localizable to a single base, at the anchored/clipped junction.
3. **Clustering paradigm:** SR callers (e.g. Splitread) "searches for clusters of split reads" to call a breakpoint. → individual split reads are grouped by position to form a supported breakpoint.

### Suzuki S et al. (2011) — "ClipCrop: a tool for detecting structural variations with single-base resolution using soft-clipping information", BMC Bioinformatics 12(Suppl 14):S7

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3287472/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed)
**Retrieved how:** WebSearch `split read structural variant breakpoint detection soft clipping single base resolution`; WebFetch of the PMC article (followed 301 redirect from ncbi.nlm.nih.gov to pmc.ncbi.nlm.nih.gov).

**Key Extracted Points:**

1. **Clip encoding:** in CIGAR, "31S69M means 31 bases from the left end are clipped, and the rest 69 bases are matched." → leading `S` = left clip; the length before/after the clip locates the junction.
2. **Breakpoint = clip boundary:** "The marginal point between a clipped sequence and matched sequence is denoted as a breakpoint." → the single-base breakpoint is exactly the junction coordinate between aligned and clipped portions.
3. **Left vs right breakpoints:** "L-breakpoints and R-breakpoints are distinguished based on which side is clipped." → orientation of the clip (leading vs trailing S) defines the breakpoint side/strand.
4. **Positional clustering:** "breakpoints are sorted and clustered within 5-base differences." → candidate breakpoints from multiple reads are merged when their positions agree within a small tolerance.

### Hart SN et al. (2013) — "SoftSearch: Integration of Multiple Sequence Features to Identify Breakpoints of Structural Variations", PLoS ONE 8(12):e83356

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3865185/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed)
**Retrieved how:** WebSearch `SoftSearch breakpoints structural variation soft-clipped reads minimum`; WebFetch of the PMC article (followed 301 redirect to pmc.ncbi.nlm.nih.gov).

**Key Extracted Points:**

1. **Breakpoint definition / minimum support:** "A putative breakpoint is defined when there is at least x soft-clipped reads beginning at position y." Default x = 5 soft-clipped reads (the paper notes "In SoftSearch_1 we decreased the minimum number of softclipped reads from 5 to 2", confirming 5 is the default baseline). → a breakpoint is called only when ≥ x clipped reads share a clip position; support = number of clipped reads at that position.
2. **Orientation-consistent combining:** "soft-clipped reads supporting a break point event are combined if the left/right orientation is in the same direction." → clipped reads are clustered into one breakpoint only when they share clip side (left vs right).
3. **Minimum clip length:** "Soft-clipped reads with more than 5 unmapped bases are passed through for further analysis." → reads with ≤ 5 clipped bases are discarded; a minimum clip length filters spurious clips.

---

## Documented Corner Cases and Failure Modes

### From SoftSearch (Hart et al. 2013)

1. **Below-support positions:** a clip position with fewer than x soft-clipped reads is not reported as a breakpoint (default x = 5; configurable down to 2). → low-support clip stacks are filtered out.
2. **Short clips:** reads with ≤ 5 unmapped (clipped) bases are not used — too short to be a reliable breakpoint signal.

### From ClipCrop (Suzuki et al. 2011)

1. **Position jitter:** mapping imprecision spreads true breakpoints across nearby positions; clustering within a small tolerance (5 bases) merges them into one call.

### From SAM spec

1. **Clip side:** a soft clip only occurs at a read end; a leading `S` (left clip) and a trailing `S` (right clip) define opposite-sided breakpoints. A read with no `S` operation carries no breakpoint signal.

---

## Test Datasets

### Dataset: Synthetic split-read breakpoint clusters (derived from the cited clip/junction rules)

**Source:** ClipCrop (junction = clip boundary; cluster within tolerance) + SoftSearch (support = clipped reads beginning at a position; default min support 5; min clip length > 5) + SAM spec (POS / CIGAR consume rules).

The repository `SplitRead` record stores the junction directly: `PrimaryPosition` = the read's anchored reference start (SAM POS), `SupplementaryPosition` = the breakpoint coordinate at the aligned/clipped junction, `ClipLength` = number of clipped bases.

| Parameter | Value |
|-----------|-------|
| Clip-position cluster tolerance | 5 bases (ClipCrop "clustered within 5-base differences") |
| Minimum supporting clipped reads (default) | 2 (SoftSearch configurable minimum; default baseline 5) |
| Minimum clip length to keep a read | > 5 clipped bases (SoftSearch) |
| Breakpoint position | aligned/clipped junction coordinate (ClipCrop "marginal point between clipped and matched sequence") |
| 3 split reads at junction 5000 ± ≤5 b, same chr, support ≥ 2 | one breakpoint at ~5000 |
| 1 isolated split read, min support 2 | no breakpoint |
| 2 split reads 5000 vs 5100 (gap > tolerance) | two separate clip groups (each below support 2) → no breakpoint |

---

## Assumptions

1. **ASSUMPTION: Breakpoint position estimator within a cluster.** The cited sources fix the per-read breakpoint as the clip/aligned junction (ClipCrop) and cluster reads sharing a clip position within a tolerance (ClipCrop 5 b; SoftSearch "beginning at position y"), but they do not prescribe a single summary statistic (mean vs mode vs min) for the cluster's reported coordinate. This unit reports the rounded mean of the member junction coordinates, mirroring the existing sibling `ClusterSplitReads` in the same class. Justification: with reads clustered inside a ≤ tolerance window the mean lies within the same single-base neighbourhood the sources define as one breakpoint; the choice does not change which reads form the cluster or the support count, only the sub-tolerance reported coordinate.

---

## Recommendations for Test Coverage

1. **MUST Test:** Split reads whose junction positions agree within the cluster tolerance and meet minimum support yield exactly one breakpoint at that junction. — Evidence: ClipCrop ("clustered within 5-base differences"; junction = breakpoint); SoftSearch ("at least x soft-clipped reads beginning at position y").
2. **MUST Test:** A clip stack with fewer than the minimum supporting reads yields no breakpoint. — Evidence: SoftSearch (default min support; configurable to 2).
3. **MUST Test:** Two clip groups separated by more than the cluster tolerance are not merged into one breakpoint. — Evidence: ClipCrop (cluster only within tolerance).
4. **MUST Test:** The reported breakpoint carries the correct support count = number of clipped reads in the cluster. — Evidence: SoftSearch (support = clipped reads at the position).
5. **MUST Test:** Reads on different chromosomes are not clustered into the same breakpoint. — Evidence: a breakpoint position is per-chromosome (SAM POS is chromosome-local; ClipCrop sorts/clusters by position within a contig).
6. **SHOULD Test:** Boundary of the cluster tolerance — junctions exactly `tolerance` apart cluster; `tolerance + 1` apart do not. — Rationale: tolerance-window correctness (ClipCrop 5 b).
7. **SHOULD Test:** `RefineBreakpoint` narrows a region to the consensus junction supported by reads inside the region. — Rationale: refinement method contract.
8. **COULD Test:** Empty input yields empty output; null input throws. — Rationale: defined trivial behavior / input-validation contract consistent with sibling methods.

---

## References

1. Li H, et al. (the SAM/BAM Format Specification Working Group). 2024. Sequence Alignment/Map Format Specification (SAMv1). samtools/hts-specs. https://samtools.github.io/hts-specs/SAMv1.pdf (source: https://raw.githubusercontent.com/samtools/hts-specs/master/SAMv1.tex)
2. Tattini L, D'Aurizio R, Magi A. 2015. Detection of Genomic Structural Variants from Next-Generation Sequencing Data. Front Bioeng Biotechnol 3:92. https://doi.org/10.3389/fbioe.2015.00092 (open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC4479793/)
3. Suzuki S, Yasuda T, Shiraishi Y, Miyano S, Nagasaki M. 2011. ClipCrop: a tool for detecting structural variations with single-base resolution using soft-clipping information. BMC Bioinformatics 12(Suppl 14):S7. https://doi.org/10.1186/1471-2105-12-S14-S7 (https://pmc.ncbi.nlm.nih.gov/articles/PMC3287472/)
4. Hart SN, Sarangi V, Moore R, et al. 2013. SoftSearch: Integration of Multiple Sequence Features to Identify Breakpoints of Structural Variations. PLoS ONE 8(12):e83356. https://doi.org/10.1371/journal.pone.0083356 (https://pmc.ncbi.nlm.nih.gov/articles/PMC3865185/)

---

## Change History

- **2026-06-13**: Initial documentation.
