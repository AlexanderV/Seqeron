# Focal Amplification Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-CNA-002 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Simplified |
| Last Reviewed | 2026-09-28 |

## 1. Overview

Focal amplification detection separates highly amplified, *focal* copy-number segments from *broad*
(arm-level) amplifications and maps the focal events to recurrently amplified oncogenes. It is
specification-driven: GISTIC2.0 classifies a copy-number event as focal or arm-level purely by its
length relative to its chromosome arm, and an event is "amplified" when its log2 gain exceeds an
amplitude threshold [1][2]. Focal high-amplitude amplifications are the therapeutically actionable
targets, so the algorithm filters segments to that subset and reports the panel oncogenes resident on
their arms [3][4]. It is a deterministic O(n) filter, not a probabilistic peak caller.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Somatic copy-number alterations (SCNAs) in cancer occur at two scales: broad events that span a whole
chromosome arm (or more), and focal events confined to a sub-arm region. GISTIC2.0 observed that events
occupying exactly one chromosome arm are so frequent that arm-level and focal alterations must be modeled
separately; their lengths form a reproducible bimodal distribution that gives "a natural basis for
classifying events as 'arm-level' and 'focal' based purely on length" [1].

### 2.2 Core Model

For a segment of length L on a chromosome arm of length A with mean log2 copy ratio r:

- **Amplified** ⇔ r > t_amp, with t_amp = 0.1 (GISTIC2 `t_amp`) [2].
- **Focal** ⇔ (L / A) < c, with c = 0.98 (GISTIC2 `broad_len_cutoff`); an event occupying ≥ 98% of an
  arm is arm-level [1][2].
- **Focal amplification** ⇔ Amplified ∧ Focal.

A single-copy gain has copy ratio 3/2, i.e. log2(3/2) = 0.585, which is far above t_amp = 0.1, so the
amplitude test admits any genuine gain while rejecting low-level artifactual segments [3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported amplification satisfies L/A < 0.98 (strict) | focal test is `ArmFraction < BroadLengthCutoff` [1][2] |
| INV-02 | Every reported amplification satisfies r > t_amp | amplitude test is `Log2Ratio > AmplificationLog2Threshold` [2] |
| INV-03 | Output is a subset of the input in input order | the detector is a filter; it constructs no new segments |
| INV-04 | An oncogene is reported only for an arm carrying a focal amplification | the mapper consumes focal amplifications only [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| segments | `IEnumerable<CopyNumberArmSegment>` | required | Arm-anchored copy-number segments | not null; each ArmLength > 0, End > Start |
| thresholds | `FocalAmplificationThresholds?` | GISTIC2 defaults | Amplitude + length cutoffs | null ⇒ t_amp 0.1, broad_len_cutoff 0.98 |
| amplifications | `IEnumerable<CopyNumberArmSegment>` | required | Focal amplifications to map to oncogenes | not null |
| segment.MarkerCount / ArmMarkerCount | `int?` (init) | null | Optional marker counts: arm fraction in marker units (GISTIC2 default) | both or neither; 1 ≤ MarkerCount ≤ ArmMarkerCount |
| amplifiedRegions | `IEnumerable<CopyNumberRegion>` | required | Locus-overlap overload: amplified regions, 1-based closed | not null; End ≥ Start; non-empty chromosome |
| genePanel | `IReadOnlyList<GeneLocus>?` | `DefaultOncogeneLoci` (GRCh38) | Gene loci for the locus-overlap overload | Start ≤ End; non-empty symbol/chromosome |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (DetectFocalAmplifications) | `IReadOnlyList<CopyNumberArmSegment>` | Input segments that are focal amplifications, in input order |
| (IdentifyAmplifiedOncogenes) | `IReadOnlyList<string>` | Distinct panel oncogene symbols on amplified arms, in panel order |
| (IdentifyAmplifiedOncogenes, `CopyNumberRegion` overload) | `IReadOnlyList<string>` | Distinct panel genes whose locus overlaps an amplified region (GISTIC2 `genes_at`), in panel order |

### 3.3 Preconditions and Validation

Null `segments`/`amplifications` ⇒ `ArgumentNullException`. Thresholds outside the GISTIC2 reference ranges
(`t_amp` ∈ [0, ∞), `broad_len_cutoff` ∈ [0, 2], NaN rejected — GISTIC2 `gp_gistic2_from_seg.m` `numeric_arg`
ranges) ⇒ `ArgumentOutOfRangeException`, checked before any segment (also for empty input). A NaN `Log2Ratio` is a
no-call (not above `t_amp`, never reported). Coordinates are half-open. A segment with non-positive `ArmLength` or
with `End ≤ Start` ⇒ `ArgumentException`. Arm labels are matched case-insensitively (Ordinal-ignore-case).
Coordinates are base-pair counts; segment length is `End − Start`. The boundary fraction exactly equal to
the cutoff (0.98) is arm-level (the focal test is strictly less-than).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate each segment (positive arm length, End > Start).
2. For each segment compute amplitude test `r > t_amp` and focal test `L/A < broad_len_cutoff`.
3. Emit the segment iff both tests pass, preserving input order.
4. For oncogene mapping, collect the set of arms carrying a focal amplification, then emit each panel
   oncogene whose arm is in that set, in panel order.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Parameter / Table | Value | Source |
|-------------------|-------|--------|
| t_amp (amplitude) | 0.1 (log2) | GISTIC2 `t_amp` [2] |
| broad_len_cutoff (focal/broad) | 0.98 (fraction of arm) | Mermel 2011; GISTIC2 `broad_len_cutoff` [1][2] |
| ERBB2 | 17q | NCBI Gene 2064 (17q12) [4] |
| MYC | 8q | NCBI Gene 4609 (8q24.21) [4] |
| EGFR | 7p | NCBI Gene 1956 (7p11.2) [4] |
| CCND1 | 11q | NCBI Gene 595 (11q13.3) [4] |
| MDM2 | 12q | NCBI Gene 4193 (12q15) [4] |
| CDK4 | 12q | NCBI Gene 1019 (12q14.1) [4] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| DetectFocalAmplifications | O(n) | O(k) | n segments, k focal amplifications |
| IdentifyAmplifiedOncogenes | O(n + g) | O(n) | g = fixed panel size (6) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.CopyNumberPloidy.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.CopyNumberPloidy.cs)

- `OncologyAnalyzer.DetectFocalAmplifications(segments, thresholds?)`: filters segments to focal amplifications.
- `OncologyAnalyzer.IdentifyAmplifiedOncogenes(amplifications)`: maps focal amplifications to panel oncogenes (arm-level).
- `OncologyAnalyzer.IdentifyAmplifiedOncogenes(IEnumerable<CopyNumberRegion>, IReadOnlyList<GeneLocus>? genePanel = null)`: GISTIC2 locus-overlap mapping (F49); default panel `DefaultOncogeneLoci`.
- `OncologyAnalyzer.IsFocalAmplification(segment, thresholds)`: single-segment predicate (internal helper, public for reuse).

### 5.2 Current Behavior

Detection is a single-pass filter; no segmentation is performed here (segmentation is upstream in
`StructuralVariantAnalyzer.SegmentCopyNumber`, SV-CNV-001 / ONCO-CNA-001). No substring/pattern search is
involved, so the repository suffix tree is **not applicable** to this unit. Arm matching uses
case-insensitive ordinal comparison; arms outside the six-gene panel map to no oncogene.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Length-based focal/arm-level split at 98% of chromosome arm (Mermel 2011; GISTIC2 `broad_len_cutoff` 0.98) [1][2].
- Amplitude threshold t_amp = 0.1 for calling a gain amplified (GISTIC2 `t_amp`) [2].
- Oncogene→arm panel from NCBI Gene cytogenetic locations [4].

**Intentionally simplified:**

- Arm boundaries / arm length are supplied by the caller rather than derived from a bundled cytoband
  file; **consequence:** the caller must provide each arm's length (GISTIC2 reads this from the genome
  assembly). The 0.98 rule and amplitude test are unchanged.
- Oncogene mapping of the `CopyNumberArmSegment` overload is arm-level (any focal amplification on the arm flags the
  gene); **consequence:** a focal amplification elsewhere on the same arm also flags the gene. The GISTIC2 rule —
  locus overlap — is available as the `CopyNumberRegion` overload (below, F49).

**Arm-fraction units (F48).** GISTIC2 measures the arm fraction in **markers** by default: `make_sample_B.m` calls
`normalize_by_arm_length(D,B,cyto,1,2,…)` with `norm_type = 1` ("by number of snps"), i.e.
`fract = (en − st + 1) ./ armlengths_by_snp`, where `armlengths_by_snp` = number of markers whose position lies in the
arm's cytoband span (`find_snps(D, chrn, band.start, band.end)`, `band.start = cyto.start + 1`); the event fraction is
the sum of its segments' fractions (`perform_deconstruction.m` `sum(Bt(:,6))`). A `CopyNumberArmSegment` that carries
`MarkerCount` and `ArmMarkerCount` uses `ArmFraction = MarkerCount / ArmMarkerCount`; without them the fraction stays in
bp (`Length / ArmLength`, GISTIC2 `norm_type = 2`). Octave running the original `normalize_by_arm_length.m`: 3/10 →
0.29999999999999999, 39/40 → 0.97499999999999998 (focal under 0.98), 40/40 → 1. Centromere-spanning segments
(GISTIC2 `ref_length = 2`: p-fraction + q-fraction, e.g. 6/10 + 20/40 = 1.1) are outside the arm-anchored model.

**Locus-overlap gene mapping (F49).** GISTIC2 maps genes to peaks in `genetables.m` via
`genes_at(rg, chr, start, end, 1, partial_hits)` with default `partial_hits = 1`: gene `idx` is reported iff
`rg.chrn == chr ∧ rg.start ≤ pos_end ∧ rg.end ≥ pos_start` (closed intervals; touching at one base counts), in
reference-gene order, deduplicated (`filt_uniq_gene`). `IdentifyAmplifiedOncogenes(IEnumerable<CopyNumberRegion>, …)`
implements this rule (panel order, each gene once; chromosome names match ignoring a leading "chr" and case). The
default panel `DefaultOncogeneLoci` holds GRCh38 loci from GISTIC2's own hg38 reference-gene source
(`refgenes/Gencode.v22.170324/gencode_genes.tsv`, GENCODE v22). Not applied: the "[closest gene]" fallback for
gene-less peaks and widening of peak boundaries to the flanking markers (`genomic_location(…,1)`) — regions are used
as given.

**Ziggurat deconstruction — per-sample building blocks (F50).** Ported from GISTIC2 (internal, Octave-locked):
`make_sample_B.m` (a sample's segments → B rows `[chrn st en amp sample fract]` at value changes and chromosome ends;
`fract` in marker units with `ref_length = 2`, so a **centromere-spanning** row gets p-part ÷ p markers + q-part ÷ q
markers — e.g. 3/8 + 4/12 = 0.70833333333333326, a whole chromosome 2) and `deconstruct_sample.m` /
`deconstruct_chr.m` / `prepare_B.m` / `atomic_zigg_deconstruction.m` / `add_broad_levels_to_zigg.m` (events relative to
given p / q broad levels; the reference quirk "breakpoint on the first of several rows ⇒ all rows use the q level" is
kept). Input model: `ZigguratChromosome(Chromosome, PArmMarkerCount, QArmMarkerCount)` (markers 1..P on p, P+1..P+Q on
q) and `ZigguratSegment(Chromosome, StartMarker, EndMarker, Value)` tiling every chromosome of a sample.

**Ziggurat deconstruction — likelihood table and arm broad level (F51).** `generate_2d_hists.m` builds the cohort
event table: a 51 × 51 histogram of event amplitude (edges −2:0.08:2) × arm fraction (edges 0:0.04:2), plus a
pseudocount of 0.01 % of the event total per bin, normalised and log-transformed (`log_hd`). `find_max_broad_level_by_table.m`
chooses an arm's broad level among its distinct segment values: for each candidate it deconstructs the arm relative to
that level (`ziggurat_on_extremes` + `iterative_ziggurat`) and scores every event and the broad event by table lookup
(`score_ziggs_by_table`: last edge strictly below the value); the highest total wins (first on ties).

**Reference-implementation cross-check (2026-09 review):** the predicate equals the GISTIC2 focal-event filter in
`snputil/reconstruct_genomes.m` (`broad_or_focal='focal'`: `Q(:,8) < broad_len_cutoff` and amplitude vs `t_amp`)
and `score_genome.m` (`Qs.del(:,8) < broad_len_cutoff`). The length test is strict `<`, as in GISTIC2. The amplitude
test is strict `>` per the GISTIC2 docs ("above") and `gene_calls.m` (`A.dat > t_amp`); `reconstruct_genomes.m`
uses `>=` — the two differ only at exact floating-point equality.

**Not implemented:**

- GISTIC2 **ziggurat deconstruction** (`perform_deconstruction.m`, `atomic_zigg_deconstruction.m`,
  cohort-learned broad levels via `find_max_broad_level_by_table`). GISTIC2 applies the focal filter to deconstructed
  SCNA events whose amplitude is relative to the underlying level; this unit treats every input segment as one event
  with amplitude = its log2. **Consequence:** raw segments of an arm-level gain interrupted by a focal peak
  (0.5 | 1.5 | 0.5) report the flanks as focal amplifications, whereas GISTIC2 calls one broad (0.5) + one focal (+1.0)
  event. Callers should pass deconstructed events or arm-merged segments.

- GISTIC2's probabilistic peak/q-value boundary estimation and background-rate modeling; **users should
  rely on:** the full GISTIC2 tool for genome-wide significance peaks. This unit implements only the
  deterministic focal/broad length rule and the oncogene panel mapping.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Amplitude test combined with length rule | Assumption | Defines which gains count as amplifications | accepted | t_amp from GISTIC2 docs (0.1); length rule from paper |
| 2 | Caller-supplied arm length | Assumption | Caller must provide cytoband-derived arm length | accepted | No bundled cytoband table |
| 3 | Arm fraction units | Option | bp unless marker counts supplied | resolved (F48) | GISTIC2 default = markers (`norm_type = 1`) |
| 4 | Gene mapping | Option | arm-level overload kept; locus overlap via `CopyNumberRegion` overload | resolved (F49) | GISTIC2 `genes_at` `partial_hits = 1` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Segment = 0.98 of arm | Not focal (arm-level) | Focal test is strict `< 0.98` [1][2] |
| Segment > 98% of arm | Not focal | Arm-level by GISTIC2 [1] |
| log2 ≤ t_amp | Not amplified | GISTIC2 `t_amp` [2] |
| Null input | ArgumentNullException | Guard |
| t_amp < 0 / NaN; broad_len_cutoff ∉ [0, 2] / NaN | ArgumentOutOfRangeException | GISTIC2 `gp_gistic2_from_seg.m` ranges |
| log2 = NaN | Not reported (no-call) | NaN is not above t_amp |
| Empty input | Empty result | Guard |
| ArmLength ≤ 0 or End ≤ Start | ArgumentException | Validation |
| Only one of MarkerCount/ArmMarkerCount, MarkerCount < 1, or MarkerCount > ArmMarkerCount | ArgumentException | Validation (F48) |
| 39/40 markers covering 99% of arm bp | Focal (0.975 < 0.98) | GISTIC2 marker units (F48) |
| Region touching a gene end at one base | Gene reported | `genes_at` closed overlap (F49) |
| Region End < Start, empty chromosome | ArgumentException | `genes_at` errors on End < Start |

### 6.2 Limitations

No ziggurat deconstruction, no significance testing, no background-rate modeling, and no sub-arm peak localization (these are
GISTIC2's probabilistic stages). The arm-level overload maps the six-gene registry panel by arm; the
`CopyNumberRegion` overload applies GISTIC2 locus overlap (F49) to the default or a caller panel. Arm fraction is in
markers when marker counts are supplied (F48), else bp. Deletions are out of scope (ONCO-CNA-003).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var segments = new[]
{
    new OncologyAnalyzer.CopyNumberArmSegment("17q", 100_000, 600_000, 1_000_000, 1.0), // focal amp
    new OncologyAnalyzer.CopyNumberArmSegment("8q", 0, 990_000, 1_000_000, 1.5),        // arm-level
};
var focal = OncologyAnalyzer.DetectFocalAmplifications(segments); // [17q segment]
var genes = OncologyAnalyzer.IdentifyAmplifiedOncogenes(focal);   // ["ERBB2"]
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_DetectFocalAmplifications_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_DetectFocalAmplifications_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ONCO-CNA-002-Evidence.md](../../../docs/Evidence/ONCO-CNA-002-Evidence.md)
- Related algorithms: [Copy_Number_Alteration_Classification](./Copy_Number_Alteration_Classification.md)

## 8. References

1. Mermel CH, Schumacher SE, Hill B, Meyerson ML, Beroukhim R, Getz G. 2011. GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. Genome Biology 12:R41. https://pmc.ncbi.nlm.nih.gov/articles/PMC3218867/
2. Broad Institute. GISTIC2 documentation (`broad_len_cutoff`, `t_amp`, `t_del`). https://broadinstitute.github.io/gistic2/ ; GISTIC2 MATLAB source https://github.com/broadinstitute/gistic2 (`source/gp_gistic2_from_seg.m`, `source/score_genome.m`, `source/gene_calls.m`, `snputil/reconstruct_genomes.m`).
3. Talevich E, Shain AH, Botton T, Bastian BC. CNVkit — Calling copy number gains and losses. https://cnvkit.readthedocs.io/en/stable/calling.html
4. NCBI Gene: ERBB2 (2064), MYC (4609), EGFR (1956), CCND1 (595), MDM2 (4193), CDK4 (1019). https://www.ncbi.nlm.nih.gov/gene/
5. GISTIC2 source (broadinstitute/gistic2 master 26c590bd): `source/normalize_by_arm_length.m`, `source/make_sample_B.m`, `source/find_snps.m`, `source/genes_at.m`, `source/genetables.m`; `refgenes/Gencode.v22.170324/gencode_genes.tsv` (GRCh38 gene loci).
