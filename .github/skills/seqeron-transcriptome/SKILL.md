---
name: seqeron-transcriptome
description: >-
  Analyze RNA-seq / transcriptome expression matrices with Seqeron (MCP tools OR
  the C# API). Use to quantify expression (TPM / FPKM from counts + gene
  lengths), normalize a multi-sample expression matrix (quantile-normalize,
  log2-transform), run two-group differential expression (log2 fold change +
  t-test + BH-adjusted p), PCA / k-means clustering / Pearson correlation /
  co-expression networks of samples or genes, pathway over-representation (ORA)
  and GSEA-like enrichment, alternative & differential splicing (PSI / deltaPSI,
  skipped-exon events), isoform switching + dominant isoforms, and RNA-seq
  library QC metrics. Triggers: "compute TPM / FPKM", "differential expression
  between conditions", "PCA of these samples", "cluster genes by expression",
  "isoform switching", "differential splicing / PSI", "co-expression network",
  "quantile normalize / log2 transform expression", "RNA-seq QC metrics", "which
  pathways are enriched in my DE genes". Server: annotation
  (TranscriptomeAnalyzer.*).
allowed-tools: Read, Bash, Grep, Glob
---

# seqeron-transcriptome — expression quantification, DE, PCA/clustering, splicing, isoforms

Routing + orchestration skill for the transcriptome / RNA-seq family on the **annotation** server
(15 tools, one backing class: `TranscriptomeAnalyzer.*`). It picks the right tool for a
quantify / normalize / differential-expression / PCA-cluster / splicing / isoform question and gives
a **dual-mode** recipe (MCP tool name **and** the equivalent `Seqeron.Genomics` C# `Method ID`).

- **Rigor is delegated.** Tool-only computation, provenance, envelope, cross-check, units, and the
  alpha / not-for-clinical caveat are owned by **[`bio-rigor`](../bio-rigor/SKILL.md)** — it applies
  by default; do not restate.
- **Don't know a tool name?** Use **[`seqeron-discovery`](../seqeron-discovery/SKILL.md)**
  (`python3 scripts/skills/find-tool.py <kw> --server annotation`) — never guess.
- **Point, don't duplicate.** Full I/O schemas live in `docs/mcp/tools/annotation/*.md`;
  algorithm invariants in `docs/algorithms/Transcriptome/*.md`.

## Scope boundary — this family vs neighbours

- **This skill OWNS** the RNA-seq *expression-matrix* workflow: TPM/FPKM quantification, cross-sample
  normalization, differential expression, PCA/clustering/correlation/co-expression, ORA/GSEA
  enrichment, differential splicing (PSI/deltaPSI + skipped-exon), isoform switching / dominance, and
  RNA-seq QC. **[`bio-annotation`](../bio-annotation/SKILL.md)** name-drops this family shallowly —
  route here.
- **Distinguishing feature: it is multi-sample / multi-gene.** Inputs are matrices
  (`{geneId, group1[], group2[]}`, `{sampleId, expression[]}`, per-gene profiles) — *not* a single
  DNA/protein string. Single-sequence annotation (ORFs, promoters, variants, motifs) stays in
  **[`bio-annotation`](../bio-annotation/SKILL.md)**.
- **Splice-*site* prediction from a genomic sequence** (donor/acceptor motifs, gene structure,
  `detect_alternative_splicing`/`SpliceSitePredictor.*`) is **bio-annotation**. This skill's splicing
  tools operate on **PSI / read-count tables**, not on sequence.
- **RNA secondary structure / MFE / stem-loops** → **[`seqeron-rna-structure`](../seqeron-rna-structure/SKILL.md)**.

## Decision guide — which tool for which question

| Question | Tool ([MCP] / `Method ID`) |
|---|---|
| Counts + gene lengths → TPM & FPKM | `calculate_tpm` / `TranscriptomeAnalyzer.CalculateTPM` |
| Make sample columns comparable (rank-based) | `quantile_normalize` / `TranscriptomeAnalyzer.QuantileNormalize` |
| Variance-stabilize expression values | `log2_transform` / `TranscriptomeAnalyzer.Log2Transform` |
| Two-group differential expression (log2FC + t-test + BH FDR) | `differential_expression` / `TranscriptomeAnalyzer.AnalyzeDifferentialExpression` |
| Project samples onto PC1/PC2 (approx.) | `perform_pca` / `TranscriptomeAnalyzer.PerformPCA` |
| k-means-like clustering of gene profiles | `cluster_genes_by_expression` / `TranscriptomeAnalyzer.ClusterGenesByExpression` |
| Pearson r between two expression vectors | `pearson_correlation` / `TranscriptomeAnalyzer.CalculatePearsonCorrelation` |
| Correlation-edge gene network | `build_coexpression_network` / `TranscriptomeAnalyzer.BuildCoExpressionNetwork` |
| Pathway over-representation among DE genes (hypergeometric) | `over_representation_analysis` / `TranscriptomeAnalyzer.PerformOverRepresentationAnalysis` |
| GSEA-like running-sum enrichment score | `enrichment_score` / `TranscriptomeAnalyzer.CalculateEnrichmentScore` |
| Differential splicing (deltaPSI between two conditions) | `detect_differential_splicing` / `TranscriptomeAnalyzer.DetectDifferentialSplicing` |
| Skipped-exon PSI from inclusion/skipping reads | `find_skipped_exon_events` / `TranscriptomeAnalyzer.FindSkippedExonEvents` |
| Isoform switching between two conditions (usage change) | `detect_isoform_switching` / `TranscriptomeAnalyzer.DetectIsoformSwitching` |
| Per-gene dominant isoform + dominance ratio | `find_dominant_isoforms` / `TranscriptomeAnalyzer.FindDominantIsoforms` |
| RNA-seq library QC (mapping/exonic/rRNA rates, detected genes) | `rnaseq_quality_metrics` / `TranscriptomeAnalyzer.CalculateQualityMetrics` |

## Canonical dual-mode pipelines

> Units: TPM sums to 1e6 across the gene set; FPKM = count·1e9/(length·N). Expression matrices are
> caller-supplied numeric vectors (no coordinates). PCA is a **lightweight top-variable-gene
> approximation, not SVD** (see envelope note). Verify exact I/O in each tool doc.

### (a) Quantify → normalize → differential expression (the main matrix workflow)
1. **[MCP]** `calculate_tpm`(geneCounts=[{geneId,rawCount,length}]) → per-gene `{tpm,fpkm}`.
2. **[MCP]** `quantile_normalize`(samples=number[][]) → matched sample distributions; then
   `log2_transform`(values, pseudocount?=1) → variance-stabilized values.
3. **[MCP]** `differential_expression`(expressionData=[{geneId,group1[],group2[]}], foldChangeThreshold?=1.0, pValueThreshold?=0.05) → per-gene `{log2FoldChange, pValue, adjustedPValue, isSignificant, regulation}`.
- **[C# API]** `TranscriptomeAnalyzer.CalculateTPM` → `.QuantileNormalize` → `.Log2Transform` → `.AnalyzeDifferentialExpression`.
```
Provenance
1) calculate_tpm(geneCounts) → tpm,fpkm (TPM sums to 1e6)
2) quantile_normalize(samples); log2_transform(values,pseudocount=1)
3) differential_expression(expressionData,foldChangeThreshold=1.0,pValueThreshold=0.05) → log2FC,pValue,adjustedPValue(BH),isSignificant,regulation
Cross-check: pearson_correlation of replicates within a group should be high (see (c)).
Caveat: alpha — validate before use.
```

### (b) DE genes → pathway enrichment (ORA + GSEA-like)
1. **[MCP]** take significant genes from pipeline (a).
2. **[MCP]** `over_representation_analysis`(differentiallyExpressedGenes, pathways=[{pathwayId,pathwayName,genes[]}], backgroundGeneCount) → per-pathway `{overlappingGenes, enrichmentScore, pValue}` (hypergeometric).
3. **[MCP]** cross-check the trend with `enrichment_score`(rankedGenes, geneSet) → GSEA-like running-sum score (rank genes by log2FC from (a)).
- **[C# API]** `TranscriptomeAnalyzer.PerformOverRepresentationAnalysis` · `.CalculateEnrichmentScore`.
```
Provenance
1) over_representation_analysis(deGenes,pathways,backgroundGeneCount) → enrichmentScore,pValue per pathway (hypergeometric)
2) enrichment_score(rankedGenes,geneSet) → running-sum score (independent, rank-based)
Cross-check: pathways ORA-enriched should trend positive under enrichment_score.
Caveat: alpha.
```

### (c) Sample structure: PCA + clustering + correlation
1. **[MCP]** `perform_pca`(samples=[{sampleId,expression[]}], topGenes?=500) → per-sample `{pc1,pc2}` (**approximation**, see envelope).
2. **[MCP]** `cluster_genes_by_expression`(geneProfiles=[{geneId,expression[]}], numClusters?=5, correlationThreshold?=0.5) → clusters `{genes[], meanCorrelation, representativeGene}`.
3. **[MCP]** corroborate independently: `pearson_correlation`(expression1, expression2) between samples/genes; `build_coexpression_network`(geneProfiles, correlationThreshold?=0.7) → edges `{gene1,gene2,correlation}`.
- **[C# API]** `TranscriptomeAnalyzer.PerformPCA` · `.ClusterGenesByExpression` · `.CalculatePearsonCorrelation` · `.BuildCoExpressionNetwork`.
```
Provenance
1) perform_pca(samples,topGenes=500) → pc1,pc2 (top-variable-gene approx, NOT SVD)
2) cluster_genes_by_expression(geneProfiles,numClusters=5) → clusters + meanCorrelation
3) pearson_correlation / build_coexpression_network(correlationThreshold=0.7) → independent correlation check
Cross-check: same-condition samples co-cluster and correlate highly.
Caveat: alpha; PCA is an approximation — do not report as SVD-PCA variance-explained.
```

### (d) Differential splicing (PSI / deltaPSI)
1. **[MCP]** `find_skipped_exon_events`(exonData=[{geneId,exonStart,exonEnd,inclusionReads,skippingReads}]) → per-exon `{inclusionLevel(PSI), deltaPsi=0}` (single-sample PSI).
2. **[MCP]** `detect_differential_splicing`(splicingData=[{geneId,start,end,psiCondition1,psiCondition2}], deltaPsiThreshold?=0.1) → events `{eventType(IncreasedInclusion/IncreasedSkipping), inclusionLevel, deltaPsi}`.
- **[C# API]** `TranscriptomeAnalyzer.FindSkippedExonEvents` · `.DetectDifferentialSplicing`.
- Note: these operate on **PSI / read-count tables**. Splice-*site* motif prediction from sequence → `bio-annotation` (`detect_alternative_splicing`/`SpliceSitePredictor.*`).
```
Provenance
1) find_skipped_exon_events(exonData) → PSI = inclusion/(inclusion+skipping)
2) detect_differential_splicing(splicingData,deltaPsiThreshold=0.1) → deltaPSI = psi2 − psi1, eventType
Caveat: alpha.
```

### (e) Isoform switching + dominant isoforms
1. **[MCP]** `find_dominant_isoforms`(isoforms) → per-gene `{dominantIsoform, dominanceRatio}` (dominant expr / total gene expr).
2. **[MCP]** `detect_isoform_switching`(isoformData=[{isoform,expression1,expression2}], switchThreshold?=0.3) → switches `{transcriptId1(down),transcriptId2(up),switchScore=|deltaUp|+|deltaDown|}`.
- **[C# API]** `TranscriptomeAnalyzer.FindDominantIsoforms` · `.DetectIsoformSwitching`.
```
Provenance
1) find_dominant_isoforms(isoforms) → dominantIsoform,dominanceRatio
2) detect_isoform_switching(isoformData,switchThreshold=0.3) → switchScore = |deltaUp|+|deltaDown|
Cross-check: a switch means the dominant isoform differs between conditions.
Caveat: alpha.
```

## End-to-end grounded example — "compare two conditions from raw counts"

**Task.** Given per-gene raw counts + lengths for replicate samples in two conditions,
(1) quantify TPM, (2) QC the libraries, (3) normalize + log2, (4) run differential expression,
(5) confirm sample structure (PCA/correlation), (6) find enriched pathways in the DE genes.

Tool / `Method ID` chain (MCP names; C# path in parentheses):
1. `calculate_tpm`(geneCounts) → per-gene TPM/FPKM. (`TranscriptomeAnalyzer.CalculateTPM`)
2. `rnaseq_quality_metrics`(totalReads, mappedReads, exonicReads, rRnaReads, geneCounts) → mapping/exonic/rRNA rates + detectedGenes. (`.CalculateQualityMetrics`)
3. `quantile_normalize`(samples) → `log2_transform`(values) → comparable, stabilized matrix. (`.QuantileNormalize` → `.Log2Transform`)
4. `differential_expression`(expressionData, foldChangeThreshold=1.0, pValueThreshold=0.05) → significant up/down genes. (`.AnalyzeDifferentialExpression`)
5. `perform_pca`(samples) + `pearson_correlation`(rep vs rep) → conditions separate; replicates correlate. (`.PerformPCA` / `.CalculatePearsonCorrelation`)
6. `over_representation_analysis`(deGenes, pathways, backgroundGeneCount) → enriched pathways. (`.PerformOverRepresentationAnalysis`)
```
Provenance
1) calculate_tpm(geneCounts) → tpm,fpkm (TPM sums to 1e6 per sample)
2) rnaseq_quality_metrics(totalReads,mappedReads,exonicReads,rRnaReads,geneCounts) → mappingRate,exonicRate,rRnaRate,detectedGenes
3) quantile_normalize(samples); log2_transform(values,pseudocount=1)
4) differential_expression(expressionData,foldChangeThreshold=1.0,pValueThreshold=0.05) → log2FC,adjustedPValue(BH),regulation
5) perform_pca(samples,topGenes=500) [approx, not SVD]; pearson_correlation(repA,repB) → replicate concordance
6) over_representation_analysis(deGenes,pathways,backgroundGeneCount) → enrichmentScore,pValue (hypergeometric)
Units: TPM/FPKM (depth+length normalized); log2FC; BH-adjusted p; PSI∈[0,1]; Pearson r∈[-1,1].
Cross-check: DE direction (regulation) consistent with raw TPM means; ORA pathways trend under enrichment_score.
Envelope: no guarded TranscriptomeAnalyzer unit (LIMITATIONS.md has no transcriptome entry).
Caveat: alpha software; not for clinical use. PCA is an approximation — do not report SVD variance-explained.
```

## Reference

- **This family's tool map (all 15 — curated index; NOT in domain-map.json, so there is NO
  generated slice):** [`reference/tool-map.md`](reference/tool-map.md)
- **Fuller recipes + parameter guidance:** [`reference/pipelines.md`](reference/pipelines.md)
- **Algorithm background (invariants/formulas — link, don't copy):**
  [`Expression_Quantification.md`](../../../docs/algorithms/Transcriptome/Expression_Quantification.md) ·
  [`Differential_Expression.md`](../../../docs/algorithms/Transcriptome/Differential_Expression.md) ·
  [`Alternative_Splicing.md`](../../../docs/algorithms/Transcriptome/Alternative_Splicing.md)
- **Operating envelope:** [`LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) — **no
  transcriptome/`TranscriptomeAnalyzer` guarded unit;** nothing here throws `SeqeronLimitationException`.
- **Cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor guardrail) ·
  [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup) ·
  [`bio-annotation`](../bio-annotation/SKILL.md) (overlap: it name-drops transcriptome/RNA-seq and owns splice-*site* prediction from sequence).
