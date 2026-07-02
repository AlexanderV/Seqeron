---
name: bio-annotation
description: >-
  Annotate and characterize genomic sequences with Seqeron (MCP tools OR the C#
  API). Use to find ORFs / predict genes / promoters / RBS, call + annotate +
  classify variants (SNPs, indels, structural variants, CNVs, VEP-like effect,
  ACMG-like pathogenicity), discover and scan motifs (exact/degenerate/PROSITE/
  PWM), analyze repeats & low-complexity / masking, profile k-mers & composition,
  plus splicing, methylation/epigenetics, miRNA targeting (seed/target pairing,
  pre-miRNA hairpins), and transcriptome/RNA-seq. (RNA secondary-structure folding →
  seqeron-rna-structure; protein-feature prediction → seqeron-protein-features.)
  Triggers: "annotate this
  sequence", "find ORFs / genes / promoters", "call variants", "what motifs are
  in…", "classify this variant", "predict the effect", "mask low-complexity",
  "find repeats", "k-mer profile". Servers: annotation + analysis (~188 tools).
allowed-tools: Read, Bash, Grep, Glob
---

# bio-annotation — structural annotation, variants, motifs, repeats, k-mers

Routing + orchestration skill for the **Annotation** (97 tools) and **Analysis** (91 tools) servers.
It picks the right tool for an annotation/characterization question and gives a **dual-mode** recipe
(MCP tool calls **and** the equivalent `Seqeron.Genomics` C# `Method ID`s).

- **Rigor is delegated.** Parse-with-a-tool, envelope, provenance, cross-check, units/0-based
  coordinates, and the alpha / not-for-clinical-use caveat are all owned by **`bio-rigor`** — it
  applies here by default; do not restate its rules. Variant pathogenicity is decision-relevant →
  always surface the clinical caveat.
- **Don't know the exact tool name?** Use **`seqeron-discovery`**
  (`python3 scripts/skills/find-tool.py <kw> --server annotation|analysis`) — never guess.
- **Point, don't duplicate.** Full I/O schemas live in `docs/mcp/tools/{annotation,analysis}/*.md`;
  algorithm invariants in `docs/algorithms/*`. This skill links; it does not copy.

## Decision guide — question → workflow family → entry tool(s)

With ~188 tools, route by **family** first, then open the family table in
[`reference/tool-map.md`](reference/tool-map.md) for the full per-tool list.

| If the task is about… | Family | Key entry tool(s) ([MCP] / `Method ID`) |
|---|---|---|
| ORFs, genes, promoters, RBS, coding potential, GFF3 | **Structural annotation** | `find_orfs`/`GenomeAnnotator.FindOrfs` · `predict_genes`/`GenomeAnnotator.PredictGenes` · `find_promoter_motifs`/`GenomeAnnotator.FindPromoterMotifs` |
| SNP/indel calling, effect, classification, pathogenicity, Ti/Tv, VCF | **Variant calling + annotation** | `call_variants`/`VariantCaller.CallVariants` · `annotate_variants`/`VariantCaller.AnnotateVariants` · `classify_variant`/`VariantAnnotator.ClassifyVariant` · `predict_pathogenicity`/`VariantAnnotator.PredictPathogenicity` |
| SVs, CNVs, breakpoints, discordant/split reads | **Structural variants** | `find_discordant_pairs`/`StructuralVariantAnalyzer.FindDiscordantPairs` · `identify_cnvs`/`StructuralVariantAnalyzer.IdentifyCNVs` · `annotate_svs`/`StructuralVariantAnalyzer.AnnotateSVs` |
| Motif discovery / exact / degenerate / PROSITE / PWM scan | **Motif discovery & scan** | `discover_motifs`/`MotifFinder.DiscoverMotifs` · `find_exact_motif`/`MotifFinder.FindExactMotif` · `create_pwm`+`scan_with_pwm`/`MotifFinder.CreatePwm`+`ScanWithPwm` |
| Tandem/inverted/direct repeats, microsatellites, palindromes | **Repeat analysis** | `find_tandem_repeats`/`GenomicAnalyzer.FindTandemRepeats` · `find_microsatellites`/`RepeatFinder.FindMicrosatellites` · `find_inverted_repeats`/`RepeatFinder.FindInvertedRepeats` |
| Low-complexity regions, DUST/SEG, masking, entropy | **Complexity / masking** | `find_low_complexity_regions`/`SequenceComplexity.FindLowComplexityRegions` · `mask_low_complexity`/`SequenceComplexity.MaskLowComplexity` · `dust_score`/`SequenceComplexity.CalculateDustScore` |
| k-mer counts/frequencies/spectrum/positions/distance | **k-mer & composition** | `count_kmers`/`KmerAnalyzer.CountKmers` · `most_frequent_kmers`/`KmerAnalyzer.FindMostFrequentKmers` · `kmer_distance`/`KmerAnalyzer.KmerDistance` · `analyze_gc_content`/`GcSkewCalculator.AnalyzeGcContent` |
| Splice sites, introns, gene structure, alt-splicing | **Splicing** | `find_donor_sites`/`SpliceSitePredictor.FindDonorSites` · `predict_gene_structure`/`SpliceSitePredictor.PredictGeneStructure` |
| Methylation, CpG islands, DMRs, chromatin, epigenetic age | **Epigenetics** | `find_cpg_islands`/`EpigeneticsAnalyzer.FindCpGIslands` · `find_dmrs`/`EpigeneticsAnalyzer.FindDMRs` |
| miRNA seeds, targets, hairpins | **miRNA** ⚠ guarded | `find_mirna_target_sites`/`MiRnaAnalyzer.FindTargetSites` (see envelope) |
| RNA secondary structure / MFE / stem-loops | → **seqeron-rna-structure** | see [../seqeron-rna-structure/SKILL.md](../seqeron-rna-structure/SKILL.md) |
| TPM/FPKM, differential expression, PCA, clustering | **Transcriptome / RNA-seq** | `calculate_tpm`/`TranscriptomeAnalyzer.CalculateTPM` · `differential_expression`/`TranscriptomeAnalyzer.AnalyzeDifferentialExpression` |
| Protein motifs/domains, disorder, TM, signal peptide, hydrophobicity profile | → **seqeron-protein-features** | see [../seqeron-protein-features/SKILL.md](../seqeron-protein-features/SKILL.md) |
| ANI, orthologs, synteny, rearrangements, dot-plot | **Comparative genomics** | `calculate_ani`/`ComparativeGenomics.CalculateANI` · `find_orthologs`/`ComparativeGenomics.FindOrthologs` |

**⚠ Envelope note.** miRNA targeting (**MIRNA-TARGET-001**, **MIRNA-CLEAVAGE-001**) is guarded.
(Disorder **DISORDER-REGION-001** and RNA-structure **RNA-STRUCT-001** envelopes now live in
seqeron-protein-features / seqeron-rna-structure.) If a call throws `SeqeronLimitationException` below its MinimumMode, **STOP and
report the limitation** — do not raise the mode to force output (bio-rigor rule 2). Details:
[`reference/pipelines.md`](reference/pipelines.md) → Envelope.

## Canonical dual-mode pipelines

Coordinates: Seqeron variant/annotation outputs are **0-based** on the reference unless a tool doc
says otherwise (`call_variants` `position` is 0-based; `predict_variant_effect` `variantPosition` is
0-based within the CDS). Report the base explicitly.

### (a) Structural annotation: sequence → ORFs → genes → promoters/RBS → GFF3
1. **[MCP]** `find_orfs`(dnaSequence, minLength=100, searchBothStrands=true) → `orfs[{start,end,frame,isReverseComplement,proteinSequence}]` (0-based, `end` incl. stop).
2. **[MCP]** `predict_genes`(dnaSequence, …) → ORF-based gene models.
3. **[MCP]** `find_promoter_motifs`(dnaSequence) → −10/−35 boxes; `find_ribosome_binding_sites`(dnaSequence) → Shine–Dalgarno.
4. **[MCP]** `to_gff3`(genes) → GFF3 lines for downstream tools.
- **[C# API]** `GenomeAnnotator.FindOrfs` → `GenomeAnnotator.PredictGenes` → `GenomeAnnotator.FindPromoterMotifs` / `FindRibosomeBindingSites` → `GenomeAnnotator.ToGff3`.
- **Cross-check:** `longest_orfs_per_frame`/`GenomeAnnotator.FindLongestOrfsPerFrame` should agree with the longest ORF per frame from step 1; `coding_potential`/`GenomeAnnotator.CalculateCodingPotential` corroborates that a called ORF is coding.

### (b) Variant workflow: reference+query → call → annotate → classify → effect
1. **[MCP]** `call_variants`(reference, query) → `variants[{position(0-based),referenceAllele,alternateAllele,type,queryPosition}]`. (Pre-aligned inputs: `call_variants_from_alignment`.)
2. **[MCP]** `annotate_variants`(reference, query, isCodingSequence=true) → per-variant `{variant, effect, mutationType}` in one shot.
3. **[MCP]** `classify_variant`(reference, alternate) → SNV/Insertion/Deletion/MNV/Indel/Complex; `predict_variant_effect`(cds, variantPosition, alternate) → protein consequence.
4. **[MCP]** `variant_statistics`(reference, query) → totals + Ti/Tv + density; `variants_to_vcf`(variants) → VCF v4.2 lines.
- **[C# API]** `VariantCaller.CallVariants` → `VariantCaller.AnnotateVariants` → `VariantAnnotator.ClassifyVariant` / `VariantCaller.PredictEffect` → `VariantCaller.CalculateStatistics` / `VariantCaller.ToVcfLines`.
- **Pathogenicity (decision-relevant):** `predict_pathogenicity`/`VariantAnnotator.PredictPathogenicity` (ACMG-like) — surface the alpha / not-for-clinical caveat and independently validate.
- **Cross-check:** `find_snps`+`find_indels` should reconcile with `call_variants`; `classify_mutation`/`VariantCaller.ClassifyMutation` (Ti/Tv/other) should match the type in `variant_statistics`.

### (c) Motif discovery → scan for occurrences
1. **[MCP]** `discover_motifs`(sequences, motifLength, …) → overrepresented k-mer motifs + a consensus. (Known set: `find_known_motifs`.)
2. **[MCP]** `create_pwm`(alignedSequences) → log-odds PWM → `scan_with_pwm`(sequence, pwm, threshold) → scored hits (0-based positions).
3. **[MCP]** exact/degenerate follow-up: `find_exact_motif`(sequence, motif) or `find_degenerate_motif`(sequence, iupacMotif).
- **[C# API]** `MotifFinder.DiscoverMotifs` → `MotifFinder.CreatePwm` → `MotifFinder.ScanWithPwm`; `MotifFinder.FindExactMotif` / `MotifFinder.FindDegenerateMotif`.
- **Cross-check (bio-rigor):** re-scan the reverse complement, or confirm an exact hit's position via `kmer_positions`/`KmerAnalyzer.FindKmerPositions`.

### (d) Repeat / low-complexity masking
1. **[MCP]** `find_low_complexity_regions`(sequence, …) and/or `find_tandem_repeats` / `find_microsatellites` / `find_inverted_repeats` → interval lists (0-based).
2. **[MCP]** `mask_low_complexity`(sequence, …) → sequence with low-complexity windows masked (feed to motif/alignment steps to avoid spurious hits).
- **[C# API]** `SequenceComplexity.FindLowComplexityRegions` / `RepeatFinder.*` → `SequenceComplexity.MaskLowComplexity`.
- **Cross-check:** `dust_score`/`SequenceComplexity.CalculateDustScore` on a flagged window should exceed the DUST threshold; `tandem_repeat_summary`/`RepeatFinder.GetTandemRepeatSummary` aggregates the microsatellite calls.

### (e) k-mer profile / composition analysis
1. **[MCP]** `count_kmers`(sequence, k) → per-k-mer counts; `kmer_frequencies` → normalized; `kmer_spectrum` → frequency-of-frequencies.
2. **[MCP]** `most_frequent_kmers` / `find_clumps` (windowed enrichment) → candidate signals (e.g. DnaA boxes near an origin).
3. **[MCP]** `kmer_distance`(seq1, seq2, k) → alignment-free composition distance; `analyze_gc_content` → GC summary.
- **[C# API]** `KmerAnalyzer.CountKmers` → `KmerAnalyzer.GetKmerFrequencies` / `GetKmerSpectrum` → `KmerAnalyzer.FindMostFrequentKmers` / `FindClumps` → `KmerAnalyzer.KmerDistance`; `GcSkewCalculator.AnalyzeGcContent`.
- **Cross-check:** GC-skew origin prediction `predict_replication_origin`/`GcSkewCalculator.PredictReplicationOrigin` should co-localize with a `find_clumps` skew-minimum signal.

## End-to-end grounded example — "characterize this locus"

**Task.** Given a bacterial DNA locus (reference) and a re-sequenced variant of it (query),
(1) find the coding ORF, (2) confirm it is coding, (3) locate its promoter/RBS, (4) call and classify
variants between reference and query, (5) predict the protein effect, (6) flag any low-complexity
that could confound the motif scan.

Tool / `Method ID` chain (MCP names; C# path in parentheses):
1. `find_orfs`(reference, minLength=50, searchBothStrands=true) → longest ORF interval + protein. (`GenomeAnnotator.FindOrfs`)
2. `coding_potential`(orf.sequence) → CPAT-style log-likelihood confirms coding. (`GenomeAnnotator.CalculateCodingPotential`)
3. `find_promoter_motifs`(reference) + `find_ribosome_binding_sites`(reference) → −10/−35 + Shine–Dalgarno upstream of the ORF. (`GenomeAnnotator.FindPromoterMotifs` / `FindRibosomeBindingSites`)
4. `annotate_variants`(reference, query, isCodingSequence=true) → `annotated[{variant{position(0-based),referenceAllele,alternateAllele,type,queryPosition}, effect, mutationType}]`. (`VariantCaller.AnnotateVariants`)
5. `variant_statistics`(reference, query) → totals, Ti/Tv, density; corroborate types via `classify_variant`. (`VariantCaller.CalculateStatistics` / `VariantAnnotator.ClassifyVariant`)
6. `find_low_complexity_regions`(reference) → intervals to exclude before trusting any motif hit. (`SequenceComplexity.FindLowComplexityRegions`)

Expected-shape output (values illustrative; **compute with the tools, do not eyeball**):
```
| feature        | coord (0-based) | detail                          |
|----------------|-----------------|---------------------------------|
| ORF            | start..end      | frame, protein MKKK…*           |
| promoter −10   | pos             | TATAAT-like box                 |
| RBS (SD)       | pos             | AGGAGG-like, upstream of ATG    |
| variant #1     | position        | ref>alt, type, effect=Missense  |
| Ti/Tv          | —               | ratio from variant_statistics   |
| low-complexity | start..end      | mask before motif interpretation|

Provenance
1) find_orfs(reference, minLength=50, searchBothStrands=true) → ORF interval + protein
2) coding_potential(orf.sequence) → CPAT log-likelihood (coding vs non-coding)
3) find_promoter_motifs(reference); find_ribosome_binding_sites(reference) → regulatory motifs
4) annotate_variants(reference, query, isCodingSequence=true) → variants (0-based) + effect + mutationType
5) variant_statistics(reference, query) → totals, Ti/Tv, density  [cross-check types: classify_variant]
6) find_low_complexity_regions(reference) → masked intervals
Coordinates: 0-based reference positions. Units: counts / bp / Ti:Tv ratio.
Cross-check: longest_orfs_per_frame agrees with step 1; find_snps+find_indels reconcile with the calls.
Envelope: none of these units guarded. (miRNA/disorder/RNA-structure ARE — see reference/pipelines.md.)
Caveat: alpha software; not for clinical use — independently validate before any decision use.
```

## Reference

- **Full domain tool index (all ~188, generated — do NOT hand-edit):** [`_generated/tools.md`](_generated/tools.md)
  (produced by `scripts/skills/gen-catalog.py`; if absent, use `seqeron-discovery`).
- **Tool map (~188 tools by family, one-liners + Method ID):** [`reference/tool-map.md`](reference/tool-map.md)
- **Fuller recipes + parameter/coordinate guidance + envelope STOP rules:** [`reference/pipelines.md`](reference/pipelines.md)
- **Algorithm background (invariants/formulas — link, don't copy):**
  [`docs/algorithms/Annotation/`](../../../docs/algorithms/Annotation/) ·
  [`Variants/`](../../../docs/algorithms/Variants/) ·
  [`Motif_Discovery/`](../../../docs/algorithms/Motif_Discovery/) ·
  [`Motif_Analysis/`](../../../docs/algorithms/Motif_Analysis/) ·
  [`Repeat_Analysis/`](../../../docs/algorithms/Repeat_Analysis/) ·
  [`Complexity/`](../../../docs/algorithms/Complexity/) ·
  [`K-mer_Analysis/`](../../../docs/algorithms/K-mer_Analysis/) ·
  [`Splicing/`](../../../docs/algorithms/Splicing/) ·
  [`Epigenetics/`](../../../docs/algorithms/Epigenetics/) ·
  [`StructuralVar/`](../../../docs/algorithms/StructuralVar/) ·
  [`MiRNA/`](../../../docs/algorithms/MiRNA/) ·
  [`Transcriptome/`](../../../docs/algorithms/Transcriptome/) ·
  [`Comparative_Genomics/`](../../../docs/algorithms/Comparative_Genomics/)
- **Cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor guardrail) · [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup).
