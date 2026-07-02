# seqeron-transcriptome — fuller pipelines & parameter guidance

Dual-mode recipes for the transcriptome / RNA-seq family (annotation server, `TranscriptomeAnalyzer.*`).
Rigor (tool-only, provenance, cross-check, alpha caveat) is delegated to
[`bio-rigor`](../../bio-rigor/SKILL.md). Schemas: `docs/mcp/tools/annotation/<tool>.md`. Tool index:
[`tool-map.md`](tool-map.md).

## Parameter defaults (verify per tool doc)

- `differential_expression`: `foldChangeThreshold` = 1.0 (|log2FC|), `pValueThreshold` = 0.05
  (BH-adjusted). Both must hold for `isSignificant`.
- `perform_pca`: `topGenes` = 500 (most-variable genes selected). **Approximation, not SVD.**
- `cluster_genes_by_expression`: `numClusters` = 5, `correlationThreshold` = 0.5 (informational only).
- `build_coexpression_network`: `correlationThreshold` = 0.7 (|Pearson r| for an edge).
- `log2_transform`: `pseudocount` = 1.0.
- `detect_differential_splicing`: `deltaPsiThreshold` = 0.1.
- `detect_isoform_switching`: `switchThreshold` = 0.3 (|delta-usage|).

## Input shapes (matrices, not sequences)

- `calculate_tpm.geneCounts` = `[{ geneId, rawCount, length }]`.
- `differential_expression.expressionData` = `[{ geneId, group1[], group2[] }]` (replicate arrays).
- `perform_pca.samples` = `[{ sampleId, expression[] }]` (all equal length).
- `cluster_genes_by_expression.geneProfiles` / `build_coexpression_network.geneProfiles` =
  `[{ geneId, expression[] }]`.
- `quantile_normalize.samples` = `number[][]`; `log2_transform.values` = `number[]`.
- `detect_differential_splicing.splicingData` = `[{ geneId, start, end, psiCondition1, psiCondition2 }]`.
- `find_skipped_exon_events.exonData` = `[{ geneId, exonStart, exonEnd, inclusionReads, skippingReads }]`.
- `detect_isoform_switching.isoformData` = `[{ isoform, expression1, expression2 }]`.
- `over_representation_analysis` = `(differentiallyExpressedGenes[], pathways[{pathwayId,pathwayName,genes[]}], backgroundGeneCount)`.
- `enrichment_score` = `(rankedGenes[], geneSet[])`.
- `rnaseq_quality_metrics` = `(totalReads, mappedReads, exonicReads, rRnaReads, geneCounts[])`.

## 1. Quantify → normalize → differential expression

**Goal.** Turn raw counts into normalized expression and a trustworthy DE call.

1. **[MCP]** `calculate_tpm`(geneCounts=[{geneId,rawCount,length}]) → per-gene `{tpm,fpkm}`. TPM sums
   to 1,000,000 across the gene set; all-zero counts emit all-0 (not an error).
2. **[MCP]** `quantile_normalize`(samples=number[][]) → identical per-sample distributions; then
   `log2_transform`(values, pseudocount=1) → variance-stabilized.
3. **[MCP]** `differential_expression`(expressionData=[{geneId,group1[],group2[]}], foldChangeThreshold=1.0,
   pValueThreshold=0.05) → per-gene `{log2FoldChange, pValue, adjustedPValue, isSignificant, regulation}`.
   log2FC is group2 vs group1 (with a pseudocount); adjustedPValue is Benjamini-Hochberg across genes.

- **[C# API]** `TranscriptomeAnalyzer.CalculateTPM` → `.QuantileNormalize` → `.Log2Transform` →
  `.AnalyzeDifferentialExpression`.

```
Provenance
1) calculate_tpm(geneCounts) → tpm,fpkm (TPM sums to 1e6; FPKM = count·1e9/(length·N))
2) quantile_normalize(samples); log2_transform(values,pseudocount=1)
3) differential_expression(expressionData,foldChangeThreshold=1.0,pValueThreshold=0.05)
   → log2FoldChange, pValue(t-test), adjustedPValue(BH), isSignificant, regulation(Up/Down/Unchanged)
Cross-check: regulation direction matches raw TPM group means; replicate pearson_correlation high.
Envelope: no guarded unit. Caveat: alpha — validate before use.
```

## 2. DE genes → pathway enrichment (ORA + GSEA-like)

**Goal.** Interpret a DE gene list at the pathway level, with an independent enrichment view.

1. Take `isSignificant` genes from pipeline 1.
2. **[MCP]** `over_representation_analysis`(differentiallyExpressedGenes, pathways=[{pathwayId,pathwayName,genes[]}],
   backgroundGeneCount) → per-pathway `{genesInPathway, overlappingGenes, enrichmentScore(>1=enriched),
   pValue, genes[]}` (hypergeometric over-representation).
3. **[MCP]** rank genes by log2FC and run `enrichment_score`(rankedGenes, geneSet) → GSEA-like
   running-sum max deviation. Two independent enrichment signals should agree in direction.

- **[C# API]** `TranscriptomeAnalyzer.PerformOverRepresentationAnalysis` · `.CalculateEnrichmentScore`.

```
Provenance
1) over_representation_analysis(deGenes,pathways,backgroundGeneCount) → overlappingGenes,enrichmentScore,pValue (hypergeometric)
2) enrichment_score(rankedGenes,geneSet) → running-sum score (rank-based, independent)
Cross-check: ORA-enriched pathways trend positive under enrichment_score.
Caveat: alpha.
```

## 3. Sample structure: PCA + clustering + correlation

**Goal.** Confirm that samples group by condition and genes co-vary as expected.

1. **[MCP]** `perform_pca`(samples=[{sampleId,expression[]}], topGenes=500) → `{pc1,pc2}` per sample.
   **Approximation** (sum of first/second half of the top-variable-gene values), **not SVD** — do not
   report variance-explained. <2 samples → origin (0,0).
2. **[MCP]** `cluster_genes_by_expression`(geneProfiles=[{geneId,expression[]}], numClusters=5,
   correlationThreshold=0.5) → clusters with member genes, meanCorrelation, representativeGene.
3. **[MCP]** corroborate: `pearson_correlation`(expression1, expression2) between two samples/genes;
   `build_coexpression_network`(geneProfiles, correlationThreshold=0.7) → `{gene1,gene2,correlation}` edges.

- **[C# API]** `TranscriptomeAnalyzer.PerformPCA` · `.ClusterGenesByExpression` ·
  `.CalculatePearsonCorrelation` · `.BuildCoExpressionNetwork`.

```
Provenance
1) perform_pca(samples,topGenes=500) → pc1,pc2  [top-variable-gene approximation, NOT SVD]
2) cluster_genes_by_expression(geneProfiles,numClusters=5,correlationThreshold=0.5) → clusters+meanCorrelation
3) pearson_correlation / build_coexpression_network(correlationThreshold=0.7) → independent correlation check
Cross-check: same-condition samples co-cluster; replicate r high; network edges match cluster co-membership.
Caveat: alpha; PCA approximate — report as PC1/PC2 projection, not SVD variance-explained.
```

## 4. Differential splicing (PSI / deltaPSI)

**Goal.** Quantify exon inclusion and its change between conditions — from read-count / PSI tables.

1. **[MCP]** `find_skipped_exon_events`(exonData=[{geneId,exonStart,exonEnd,inclusionReads,skippingReads}])
   → per-exon PSI = inclusion/(inclusion+skipping); single-sample so `deltaPsi`=0; eventType `SkippedExon`.
2. **[MCP]** `detect_differential_splicing`(splicingData=[{geneId,start,end,psiCondition1,psiCondition2}],
   deltaPsiThreshold=0.1) → `{eventType(IncreasedInclusion/IncreasedSkipping), inclusionLevel(PSI2), deltaPsi=psi2−psi1}`.

- **[C# API]** `TranscriptomeAnalyzer.FindSkippedExonEvents` · `.DetectDifferentialSplicing`.
- These consume **PSI / read counts**, not sequence. Splice-*site* motif prediction from a genomic
  sequence (`detect_alternative_splicing`/`SpliceSitePredictor.*`) → [`bio-annotation`](../../bio-annotation/SKILL.md).

```
Provenance
1) find_skipped_exon_events(exonData) → PSI = inclusion/(inclusion+skipping)
2) detect_differential_splicing(splicingData,deltaPsiThreshold=0.1) → deltaPSI = psi2 − psi1, eventType
Caveat: alpha.
```

## 5. Isoform switching + dominant isoforms

**Goal.** Find genes whose dominant transcript changes between conditions.

1. **[MCP]** `find_dominant_isoforms`(isoforms) → per-gene `{dominantIsoform, dominanceRatio}`
   (dominant expression / total gene expression).
2. **[MCP]** `detect_isoform_switching`(isoformData=[{isoform,expression1,expression2}], switchThreshold=0.3)
   → `{transcriptId1(usage decreased), transcriptId2(usage increased), switchScore=|deltaUp|+|deltaDown|}`.
   Usage = isoform expression ÷ total gene expression per condition; genes need ≥2 isoforms.

- **[C# API]** `TranscriptomeAnalyzer.FindDominantIsoforms` · `.DetectIsoformSwitching`.

```
Provenance
1) find_dominant_isoforms(isoforms) → dominantIsoform,dominanceRatio
2) detect_isoform_switching(isoformData,switchThreshold=0.3) → transcriptId1(down),transcriptId2(up),switchScore
Cross-check: a reported switch implies the per-condition dominant isoform differs.
Caveat: alpha.
```

## 6. RNA-seq library QC

**Goal.** Sanity-check a library before trusting its expression values.

- **[MCP]** `rnaseq_quality_metrics`(totalReads, mappedReads, exonicReads, rRnaReads, geneCounts) →
  `{mappingRate=mapped/total, exonicRate=exonic/mapped, rRnaRate=rRNA/mapped, detectedGenes(count>0)}`.
- **[C# API]** `TranscriptomeAnalyzer.CalculateQualityMetrics`.
- Use as a gate: low mappingRate / high rRnaRate / few detectedGenes → distrust downstream DE.

```
Provenance
rnaseq_quality_metrics(totalReads,mappedReads,exonicReads,rRnaReads,geneCounts) → mappingRate,exonicRate,rRnaRate,detectedGenes
Caveat: alpha.
```

## Envelope

- **No guarded `TranscriptomeAnalyzer` unit** in
  [`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md) — no
  `SeqeronLimitationException` from these tools; no MinimumMode STOP applies. The only caveat is
  documentary: **`perform_pca` is a top-variable-gene approximation, not SVD** — present PC1/PC2 as a
  lightweight projection, never as SVD variance-explained.

## Scope reminders

- **Splice-*site* prediction from sequence** (donor/acceptor, gene structure,
  `detect_alternative_splicing`) → [`bio-annotation`](../../bio-annotation/SKILL.md).
- **RNA secondary structure / MFE / stem-loops** → [`seqeron-rna-structure`](../../seqeron-rna-structure/SKILL.md).
- **miRNA seed/target pairing, pre-miRNA hairpins** → [`bio-annotation`](../../bio-annotation/SKILL.md) (guarded units).
