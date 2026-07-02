# seqeron-transcriptome — tool map (all 15)

Server: **annotation**. One backing class: `TranscriptomeAnalyzer.*`.
This skill is **not** in `domain-map.json`, so it has **no** generated `_generated/tools.md` —
**this curated map is the index.** Verify schemas in `docs/mcp/tools/annotation/<tool>.md`.

> Inputs are **caller-supplied numeric matrices** (counts, expression vectors, PSI tables), not
> sequences — no genomic coordinates unless a splicing tool carries exon `start`/`end`. Units: TPM
> sums to 1,000,000 across the gene set; FPKM = count·1e9/(length·N); log2FC (group2 vs group1);
> BH-adjusted p-values; PSI ∈ [0,1]; Pearson r ∈ [-1,1]. Always confirm exact I/O in the tool doc.

## Quantification & normalization

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `calculate_tpm` | `TranscriptomeAnalyzer.CalculateTPM` | Raw counts + gene lengths → per-gene TPM and FPKM. TPM sums to 1e6; all-zero counts → all 0. | `calculate_tpm.md` |
| `quantile_normalize` | `TranscriptomeAnalyzer.QuantileNormalize` | Bolstad (2003) quantile normalization: rank-mean across equal-length sample vectors. | `quantile_normalize.md` |
| `log2_transform` | `TranscriptomeAnalyzer.Log2Transform` | `log2(value + pseudocount)` per value (pseudocount default 1) — variance stabilization. | `log2_transform.md` |

## Differential expression & enrichment

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `differential_expression` | `TranscriptomeAnalyzer.AnalyzeDifferentialExpression` | Two-group log2FC + t-test p + BH-adjusted p; `isSignificant` iff \|log2FC\|≥`foldChangeThreshold`(1.0) & adjP<`pValueThreshold`(0.05); regulation Up/Down/Unchanged. | `differential_expression.md` |
| `over_representation_analysis` | `TranscriptomeAnalyzer.PerformOverRepresentationAnalysis` | Hypergeometric pathway ORA among DE genes vs `backgroundGeneCount`; per-pathway overlap, enrichmentScore(>1=enriched), p-value. | `over_representation_analysis.md` |
| `enrichment_score` | `TranscriptomeAnalyzer.CalculateEnrichmentScore` | GSEA-like running sum over a ranked gene list vs a gene set → max deviation (score). | `enrichment_score.md` |

## Sample / gene structure

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `perform_pca` | `TranscriptomeAnalyzer.PerformPCA` | Project samples onto PC1/PC2 using the `topGenes`(500) most variable genes. **Lightweight approximation, NOT SVD-PCA**; <2 samples → origin. | `perform_pca.md` |
| `cluster_genes_by_expression` | `TranscriptomeAnalyzer.ClusterGenesByExpression` | k-means-like clustering of gene expression profiles → clusters with members, meanCorrelation, representativeGene (`numClusters`=5, `correlationThreshold`=0.5 informational). | `cluster_genes_by_expression.md` |
| `pearson_correlation` | `TranscriptomeAnalyzer.CalculatePearsonCorrelation` | Pearson r ∈ [-1,1] between two equal-length expression vectors. | `pearson_correlation.md` |
| `build_coexpression_network` | `TranscriptomeAnalyzer.BuildCoExpressionNetwork` | Upper-triangle pairwise Pearson over gene profiles; edge iff \|r\|≥`correlationThreshold`(0.7). | `build_coexpression_network.md` |

## Splicing (PSI / read-count tables)

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_skipped_exon_events` | `TranscriptomeAnalyzer.FindSkippedExonEvents` | Per-exon PSI = inclusion/(inclusion+skipping) from read counts; single-sample (`deltaPsi`=0), eventType `SkippedExon`. | `find_skipped_exon_events.md` |
| `detect_differential_splicing` | `TranscriptomeAnalyzer.DetectDifferentialSplicing` | deltaPSI = psi2−psi1 between two conditions; event when \|deltaPSI\|≥`deltaPsiThreshold`(0.1); eventType IncreasedInclusion/IncreasedSkipping. | `detect_differential_splicing.md` |

## Isoforms

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_dominant_isoforms` | `TranscriptomeAnalyzer.FindDominantIsoforms` | Per gene, highest-expression isoform + dominanceRatio (dominant / total gene expression). | `find_dominant_isoforms.md` |
| `detect_isoform_switching` | `TranscriptomeAnalyzer.DetectIsoformSwitching` | Per-condition usage proportions per gene; switch when \|delta-usage\|≥`switchThreshold`(0.3); switchScore = \|deltaUp\|+\|deltaDown\|. | `detect_isoform_switching.md` |

## QC

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `rnaseq_quality_metrics` | `TranscriptomeAnalyzer.CalculateQualityMetrics` | Library QC rates: mappingRate=mapped/total, exonicRate=exonic/mapped, rRnaRate=rRNA/mapped, detectedGenes (count>0). | `rnaseq_quality_metrics.md` |

## Envelope

- **No guarded unit.** `docs/Validation/LIMITATIONS.md` has **no** transcriptome /
  `TranscriptomeAnalyzer` entry — no `SeqeronLimitationException` is thrown by these tools. The only
  model caveat is documentary: `perform_pca` is a top-variable-gene **approximation, not SVD** — do
  not report SVD-style variance-explained. See
  [`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md).

## Not a tool here (route elsewhere)

- **Splice-*site* prediction from a genomic sequence** (donor/acceptor motifs, gene structure,
  `detect_alternative_splicing` / `SpliceSitePredictor.*`) → [`bio-annotation`](../../bio-annotation/SKILL.md).
  This skill's splicing tools consume **PSI / read-count tables**, not sequence.
- **RNA secondary structure / MFE / stem-loops** → [`seqeron-rna-structure`](../../seqeron-rna-structure/SKILL.md).
- **miRNA seed/target pairing, pre-miRNA hairpins** → [`bio-annotation`](../../bio-annotation/SKILL.md) (guarded).
