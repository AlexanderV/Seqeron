# seqeron-structural-variants — tool map (all 12)

Server: **annotation**. One backing class: `StructuralVariantAnalyzer.*`.
This skill is **not** in `domain-map.json`, so it has **no** generated `_generated/tools.md` —
**this curated map is the index.** Verify schemas in `docs/mcp/tools/annotation/<tool>.md`.

> Coordinates: positions are genomic; SV `Start/End` and probe/breakpoint positions follow the tool
> docs (confirm 0- vs 1-based per doc). Copy number is an integer per segment
> (`round(2·2^meanLogR)`, clamped 0–10). Genotype qualities cap at 99. Always confirm exact I/O in
> the tool doc.

## Paired-end (discordant) evidence

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_discordant_pairs` | `StructuralVariantAnalyzer.FindDiscordantPairs` | Flag read pairs inconsistent with the library: interchromosomal (translocation), insert size outside `mean ± 3·sd` (del/ins), non-FR orientation (inv/dup), or above `maxInsertSize`. Concordant FR pairs excluded. O(n). | `find_discordant_pairs.md` |
| `cluster_discordant_pairs` | `StructuralVariantAnalyzer.ClusterDiscordantPairs` | Group discordant pairs sharing a chromosome pair within `clusterDistance`; emit one SV per cluster with ≥ `minSupport` pairs; SV type inferred from the paired-end signature. | `cluster_discordant_pairs.md` |

## Split-read evidence + breakpoint

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_split_reads` | `StructuralVariantAnalyzer.FindSplitReads` | Parse CIGAR soft-clips (`S`) ≥ `minClipLength`; emit primary + inferred supplementary position, clip length, clipped sequence. Single-base breakpoint evidence. | `find_split_reads.md` |
| `cluster_split_reads` | `StructuralVariantAnalyzer.ClusterSplitReads` | Group same-chr split reads within `clusterDistance` of primary position; emit breakpoint per cluster ≥ `minSupport`; position = mean(primary,supplementary), quality = `min(support·15, 100)`. | `cluster_split_reads.md` |
| `assemble_breakpoint_sequence` | `StructuralVariantAnalyzer.AssembleBreakpointSequence` | Heuristic junction sequence = clipped sequence of the read with the largest `clipLength`; `null` if no reads. `minOverlap` reserved (not used today). | `assemble_breakpoint_sequence.md` |
| `find_microhomology` | `StructuralVariantAnalyzer.FindMicrohomology` | Longest (≤ `maxLength`) sequence that is both a suffix of `leftFlank` and a prefix of `rightFlank`; case-insensitive, upper-cased result. MMBIR/MMEJ signature. | `find_microhomology.md` |

## Read-depth copy number

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `segment_copy_number` | `StructuralVariantAnalyzer.SegmentCopyNumber` | CBS-*like* segmentation of sorted probes: new segment on chr change or `|logR − segMean| > changeThreshold` once segment has ≥ `minProbes`. Reports span, meanLogR, `copyNumber = round(2·2^meanLogR)` (0–10), mean BAF, probeCount. Short segments dropped. | `segment_copy_number.md` |
| `identify_cnvs` | `StructuralVariantAnalyzer.IdentifyCNVs` | Segments with copyNumber ≠ `normalCopyNumber` and length ≥ `minLength` → **Deletion** (below) / **Duplication** (above). id `CNV1…`, quality = `|logRatio|·50`, supportingReads = probeCount. | `identify_cnvs.md` |

## SV post-processing

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `genotype_sv` | `StructuralVariantAnalyzer.GenotypeSV` | Diploid genotype from altFraction = `altReads/totalReads`: `<0.1`→`0/0`, `>0.9`→`1/1`, `0.3–0.7`→`0/1`, else `0/1`; quality capped 99; `totalReads=0`→`./.` q0. Non-negative read counts required. | `genotype_sv.md` |
| `filter_svs` | `StructuralVariantAnalyzer.FilterSVs` | Keep SVs with quality ≥ `minQuality`, support ≥ `minSupport`, and length in `[minLength, maxLength]` — **all four** must hold. | `filter_svs.md` |
| `merge_overlapping_svs` | `StructuralVariantAnalyzer.MergeOverlappingSVs` | Sort by chr+start; merge consecutive same-chr, same-type SVs whose overlap ÷ smaller-length ≥ `overlapFraction`. Merged = union span, first id, summed support, max quality. Only adjacent-in-sort merges. | `merge_overlapping_svs.md` |
| `annotate_svs` | `StructuralVariantAnalyzer.AnnotateSVs` | Overlapping genes (`sv.Start ≤ gene.End && sv.End ≥ gene.Start`) + hit exons (`GENE:exonN`). Impact: exon hit → HIGH (Del/Inv/Transloc) or MODERATE (Dup/other); gene-only → MODIFIER; none → LOW. `isPathogenic` = HIGH∪MODERATE; `populationFrequency` always 0. | `annotate_svs.md` |

## Envelope

- **No guarded unit.** `docs/Validation/LIMITATIONS.md` has **no** SV/CNV row; Test-Unit IDs
  `SV-DETECT-001`, `SV-BREAKPOINT-001`, `SV-CNV-001` do not appear there. Nothing throws a
  `SeqeronLimitationException` — do not fabricate a STOP rule. All three algorithms are
  Implementation-Status **Simplified** (heuristic signature-then-cluster; CBS-*like*, not
  statistical change-point; longest-clip assembly, not true OLC; no purity/GC model). See
  [`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md) and the envelope
  note in [`../SKILL.md`](../SKILL.md).

## Not a tool here (route elsewhere)

- **Tumor allele-specific / major-minor copy number, ASCAT purity-ploidy, LOH/HRD, clonal structure**
  → [`seqeron-oncology`](../../seqeron-oncology/SKILL.md) (C#-API only). `segment_copy_number` here is
  single-copy-number, germline-diploid, no purity model.
- **Chromosome-arm / whole-chromosome amplification-deletion, aneuploidy (trisomy/monosomy), karyotype,
  band-scale rearrangements** → [`bio-chromosome`](../../bio-chromosome/SKILL.md). CNVs here are focal
  segment-level, not arm-scale ploidy.
- **SNP / indel calling (base-level ref-vs-query variants), VEP-like effect, ACMG classification**
  → [`bio-annotation`](../../bio-annotation/SKILL.md). This skill is large rearrangements from read
  signatures, not base-resolution variant calling.
