---
name: bio-metagenomics
description: >-
  Profile and compare microbial communities with Seqeron (MCP tools OR the C#
  API). Use for taxonomic classification of reads/contigs (Kraken k-mer/LCA),
  community/abundance profiling, alpha diversity (Shannon, Simpson, inverse
  Simpson, Chao1, Pielou evenness), beta diversity between samples (Bray-Curtis,
  Jaccard), differential abundance, metagenomic binning into MAGs + completeness/
  contamination, functional prediction, antibiotic-resistance gene detection, and
  pan-genome analysis (core/accessory, Heaps' law, phylogenetic markers). Triggers:
  "classify these reads", "what's the community composition", "profile this
  metagenome", "compute Shannon/Simpson diversity", "compare these two samples",
  "bin these contigs", "build a pan-genome". Server: metagenomics.
allowed-tools: Read, Bash, Grep, Glob
---

# bio-metagenomics — taxonomic classification, community profiling, diversity, binning

Routing + orchestration skill for the **Metagenomics** server (19 tools). It picks the right tool
for a community-ecology / metagenome question and gives a **dual-mode** recipe (MCP tool calls
**and** the equivalent `Seqeron.Genomics` C# `Method ID`s). Two method families back the server:
`MetagenomicsAnalyzer` (classification, profiling, diversity, binning, function) and
`PanGenomeAnalyzer` (multi-genome pan-genome analysis).

- **Rigor is delegated.** Parse-with-a-tool, envelope, provenance, cross-check, units, and the
  alpha / not-for-clinical-use caveat are owned by **`bio-rigor`** — it applies here by default; do
  not restate its rules.
- **Don't know the tool name?** Use **`seqeron-discovery`**
  (`python3 scripts/skills/find-tool.py <kw> --server metagenomics`) — never guess.
- **Point, don't duplicate.** Full I/O schemas live in `docs/mcp/tools/metagenomics/*.md`; algorithm
  invariants in `docs/algorithms/Metagenomics/*.md`. This skill links, it does not copy.

## ⚠ Envelope — one guarded unit in scope

**`bin_contigs` / `MetagenomicsAnalyzer.BinContigs` is guarded by `META-BIN-001`** (domain-level
CheckM completeness/contamination), minimum access mode **`Moderate`** (the default). Under
**`Strict`** it throws `SeqeronLimitationException`. **STOP rule:** if a binning call throws, report
the limitation and do not force output — surface it and (in C#) bootstrap `Moderate`/`Permissive`
only when the caller accepts the domain-level approximation. All other 18 tools are unguarded.

## Decision guide — which tool for which question

| Question | Tool ([MCP] / `Method ID`) |
|---|---|
| Build a Kraken k-mer→taxon **database** from references | `build_kmer_database` / `MetagenomicsAnalyzer.BuildKmerDatabase` |
| **Classify reads** to taxa (Kraken k-mer / LCA) | `classify_reads` / `MetagenomicsAnalyzer.ClassifyReads` |
| Aggregate classifications → **abundance profile** (+ Shannon/Simpson) | `taxonomic_profile` / `MetagenomicsAnalyzer.GenerateTaxonomicProfile` |
| **Alpha diversity** of one sample (Shannon, Simpson, Chao1, evenness) | `alpha_diversity` / `MetagenomicsAnalyzer.CalculateAlphaDiversity` |
| **Beta diversity** between two samples (Bray-Curtis, Jaccard) | `beta_diversity` / `MetagenomicsAnalyzer.CalculateBetaDiversity` |
| Per-taxon **differential abundance** between two condition groups | `differential_abundance` / `MetagenomicsAnalyzer.DifferentialAbundance` |
| **Bin contigs** into MAGs (completeness/contamination) ⚠guarded | `bin_contigs` / `MetagenomicsAnalyzer.BinContigs` |
| Predict **functions** for proteins | `predict_functions` / `MetagenomicsAnalyzer.PredictFunctions` |
| **Functional diversity** (richness + Shannon over functions) | `functional_diversity` / `MetagenomicsAnalyzer.CalculateFunctionalDiversity` |
| Detect **antibiotic-resistance** genes | `find_resistance_genes` / `MetagenomicsAnalyzer.FindResistanceGenes` |
| **Pan-genome** from many genomes / core/accessory / markers | `construct_pangenome`, `cluster_genes`, … / `PanGenomeAnalyzer.*` |

Rule of thumb: **classification** needs a database first (`build_kmer_database` → `classify_reads`);
**profiling** aggregates those classifications; **alpha** is within-sample, **beta** is between-sample;
**binning** needs assembled contigs with coverage; **pan-genome** tools need per-genome gene sets
(see [`reference/tool-map.md`](reference/tool-map.md) for all 19).

## Canonical dual-mode pipelines

### (a) Reads → taxonomic classification → abundance profile
1. **[MCP]** `build_kmer_database`(referenceGenomes, taxonomy, k?=31) → `entries[]` (canonical k-mer → LCA taxon).
2. **[MCP]** `classify_reads`(reads, kmerDatabase=entries, taxonomy, k?=31) → per-read `taxonId`, `rank`, `confidence`.
3. **[MCP]** `taxonomic_profile`(classifications=items) → `speciesAbundance[]`, `shannonDiversity`, `simpsonDiversity`, `totalReads`, `classifiedReads`.
- **[C# API]** `MetagenomicsAnalyzer.BuildKmerDatabase(refs,taxonomy,k)` → `.ClassifyReads(reads,db,taxonomy,k)` → `.GenerateTaxonomicProfile(classifications)`.
- **k must match** between build and classify (default 31); canonicalisation collapses shared k-mers to the taxonomy LCA.
```
Provenance
1) build_kmer_database(refs,taxonomy,k=31) → entries, count
2) classify_reads(reads,entries,taxonomy,k=31) → items[{taxonId,rank,confidence=C/Q}]
3) taxonomic_profile(items) → speciesAbundance, shannon, simpson, classifiedReads/totalReads
Cross-check: profile's species-level shannon/simpson vs alpha_diversity(speciesAbundance) — must agree.
Envelope: none guarded on this path. Caveat: alpha — validate before decision use.
```

### (b) Profile / abundance vector → alpha diversity (Shannon / Simpson / Chao1 / evenness)
1. **[MCP]** `alpha_diversity`(abundances=`[{name,fraction}]`) → `shannonIndex`, `simpsonIndex`, `inverseSimpson`, `chao1Estimate`, `observedSpecies`, `pielouEvenness`.
- **[C# API]** `MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances)`.
- Accepts raw counts **or** fractions (normalised internally); zero entries dropped. Chao1 needs integer counts, else falls back to observed richness.
- **Cross-check:** feed the `speciesAbundance` from pipeline (a) — its `shannonDiversity`/`simpsonDiversity` must equal `alpha_diversity`'s `shannonIndex`/`simpsonIndex`.

### (c) Two samples → beta diversity (between-sample distance)
1. **[MCP]** `beta_diversity`(sample1Name, sample1=`[{name,fraction}]`, sample2Name, sample2) → `brayCurtis`, `jaccardDistance`, `sharedSpecies`, `uniqueToSample1/2` (`uniFracDistance`=0, no tree).
- **[C# API]** `MetagenomicsAnalyzer.CalculateBetaDiversity(name1,s1,name2,s2)`.
- Bray-Curtis uses abundances; Jaccard uses presence/absence (`>0`). For a per-taxon test across replicates use `differential_abundance` instead.

### (d) Assembled contigs → binning into MAGs (⚠ guarded META-BIN-001)
1. **[MCP]** `bin_contigs`(contigs=`[{contigId,sequence,coverage}]`, numBins?=10, minBinSize?=500000, expectedGenomeSize?=4000000) → `items[]` with `contigIds`, `totalLength`, `gcContent`, `coverage`, `completeness`, `contamination`.
- **[C# API]** `MetagenomicsAnalyzer.BinContigs(contigs,numBins,minBinSize,expectedGenomeSize)` — requires `LimitationPolicy` ≥ `Moderate`; under `Strict` it throws `SeqeronLimitationException` (STOP rule above).
- Features: GC, normalised coverage, tetranucleotide (TETRA) frequency; k-means with deterministic seeding. `completeness`/`contamination` are **domain-level CheckM approximations** — report the caveat.

### (e) Multi-genome set → pan-genome (core / accessory / markers)
1. **[MCP]** `construct_pangenome`(genomes) → clusters. 2. `core_gene_clusters` → core set; `accessory_genes` → accessory summary. 3. `select_phylogenetic_markers` → single-copy core markers; `fit_heaps_law` → openness of the pan-genome.
- **[C# API]** `PanGenomeAnalyzer.ConstructPanGenome` → `.GetCoreGeneClusters` / `.AnalyzeAccessoryGenes` / `.SelectPhylogeneticMarkers` / `.FitHeapsLaw`.
- Full tool list + Method IDs: [`reference/tool-map.md`](reference/tool-map.md).

## End-to-end grounded example (extends the metagenomics tool docs, §7.6)

**Task.** From a tiny reference set and a handful of reads, (1) classify the reads, (2) build a
sample profile, (3) report within-sample diversity, then (4) compare two such samples.

Tool / `Method ID` chain (MCP names; C# path in parentheses):
1. `build_kmer_database`(referenceGenomes=[(100,"AGCTAAAA"),(101,"AGCTCCCC")], taxonomy, k=4)
   → `entries` (shared `AGCT` collapses to genus LCA 20). (`MetagenomicsAnalyzer.BuildKmerDatabase`)
2. `classify_reads`(reads, kmerDatabase=entries, taxonomy, k=4) → per-read `taxonId`+`confidence`
   (unanimous k-mers → species; split k-mers → genus LCA). (`MetagenomicsAnalyzer.ClassifyReads`)
3. `taxonomic_profile`(classifications=items) → `speciesAbundance`, `shannonDiversity`,
   `simpsonDiversity`, `classifiedReads`/`totalReads`. (`MetagenomicsAnalyzer.GenerateTaxonomicProfile`)
4. `alpha_diversity`(abundances=speciesAbundance) → must reproduce step-3 Shannon/Simpson
   (independent code path — cross-check). (`MetagenomicsAnalyzer.CalculateAlphaDiversity`)
5. Second sample the same way → `beta_diversity`(s1,s2) → `brayCurtis`, `jaccardDistance`,
   `sharedSpecies`. (`MetagenomicsAnalyzer.CalculateBetaDiversity`)

Expected-shape output (values illustrative; **compute them with the tools, do not eyeball**):
```
| sample | classified/total | shannon | simpson | observed_species |
|--------|-----------------:|--------:|--------:|-----------------:|
| S1     | 4/5              |   0.693 |    0.50 |                2 |

Provenance
1) build_kmer_database(refs,taxonomy,k=4) → entries
2) classify_reads(reads,entries,taxonomy,k=4) → items[{taxonId,rank,confidence=C/Q}]
3) taxonomic_profile(items) → speciesAbundance, shannon, simpson, classified/total
4) alpha_diversity(speciesAbundance) → shannonIndex, simpsonIndex (== step 3 — cross-check)
5) beta_diversity(S1,S2) → brayCurtis, jaccardDistance, sharedSpecies (uniFrac=0, no tree)
Envelope: none guarded on this path (binning would be META-BIN-001, Moderate).
Caveat: alpha software; not for clinical use — validate before relying on any call.
```

## Reference

- **Full domain tool index (all 19, generated — do NOT hand-edit):** [`_generated/tools.md`](_generated/tools.md)
  (produced by `scripts/skills/gen-catalog.py`; if absent, run `seqeron-discovery`).
- **Fuller recipes + parameter/gotcha guidance:** [`reference/pipelines.md`](reference/pipelines.md)
- **Tool map (all 19 grouped by sub-task, one-liners + Method ID):** [`reference/tool-map.md`](reference/tool-map.md)
- **Algorithm background (invariants/formulas — link, don't copy):**
  [`docs/algorithms/Metagenomics/Taxonomic_Classification.md`](../../../docs/algorithms/Metagenomics/Taxonomic_Classification.md) ·
  [`Taxonomic_Profile.md`](../../../docs/algorithms/Metagenomics/Taxonomic_Profile.md) ·
  [`Alpha_Diversity.md`](../../../docs/algorithms/Metagenomics/Alpha_Diversity.md) ·
  [`Beta_Diversity.md`](../../../docs/algorithms/Metagenomics/Beta_Diversity.md) ·
  [`Genome_Binning.md`](../../../docs/algorithms/Metagenomics/Genome_Binning.md) (⚠ META-BIN-001) ·
  [`Significant_Taxa_Detection.md`](../../../docs/algorithms/Metagenomics/Significant_Taxa_Detection.md) ·
  [`PanGenome_Core_Accessory.md`](../../../docs/algorithms/Metagenomics/PanGenome_Core_Accessory.md)
- **Envelope source of truth:** [`docs/Validation/LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) (META-BIN-001).
- **Cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor guardrail) · [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup).
