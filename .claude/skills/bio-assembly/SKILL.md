---
name: bio-assembly
description: >-
  Assemble sequencing reads into contigs and QC the result with Seqeron
  (MCP tools OR the C# API). Use to assemble reads via de Bruijn graph or
  overlap-layout-consensus (OLC), build/merge/scaffold contigs, compute read
  coverage over a reference, and report assembly statistics / N50; plus k-mer
  counting, k-mer spectrum, most-frequent / min-count / unique k-mers and
  clumps for assembly QC, and assembly-relevant repeat & low-complexity
  screening (direct/inverted/tandem repeats, palindromes, microsatellites,
  longest repeat). Triggers: "assemble these reads", "de Bruijn / OLC assembly",
  "compute N50 / assembly stats", "coverage of these reads over the reference",
  "k-mer spectrum of…", "most frequent k-mers", "find clumps", "screen for
  repeats before assembly". Servers: analysis + core (k-mer/repeat subset); the
  assembly-graph engine (de Bruijn/OLC/coverage/N50) is on the Alignment server.
allowed-tools: Read, Bash, Grep, Glob
---

# bio-assembly — read assembly, coverage, N50, k-mer & repeat QC

Routing + orchestration skill for the **assembly** workflow family. Its owned slice is the
**Analysis** + **Core** servers' assembly-relevant subset — **k-mer** counting/spectrum/clumps and
**repeat / low-complexity** screening. The de-Bruijn / OLC / coverage / N50 **assembly engine** lives
on the **Alignment** server (`SequenceAssembler.*`); this skill drives it end-to-end and cites its real
tool names. It gives a **dual-mode** recipe (MCP tool calls **and** the equivalent `Seqeron.Genomics`
C# `Method ID`s).

- **Rigor is delegated.** Parse-with-a-tool, envelope, provenance, cross-check, units / 0-based
  coordinates, and the alpha / not-for-clinical-use caveat are owned by **[`bio-rigor`](../bio-rigor/SKILL.md)** —
  it applies here by default; do not restate its rules.
- **Don't know the tool name?** Use **[`seqeron-discovery`](../seqeron-discovery/SKILL.md)**
  (`python3 scripts/skills/find-tool.py <kw> --server analysis|core`) — never guess.
- **Point, don't duplicate.** Full I/O schemas live in `docs/mcp/tools/{analysis,core,alignment}/*.md`;
  algorithm invariants in `docs/algorithms/{Assembly,Extended_Assembly,K-mer,K-mer_Analysis,Repeat_Analysis}/*.md`.
  This skill links, it does not copy.

## Scope split — read this first

| Belongs here (bio-assembly) | Goes elsewhere |
|---|---|
| Assemble reads → contigs (de Bruijn / OLC), merge/scaffold, coverage, N50/assembly stats | — |
| k-mer counting / spectrum / most-frequent / min-count / unique / clumps for **assembly QC** | — |
| assembly-relevant **repeats** (direct/inverted/tandem/palindrome/microsatellite/longest-repeat) & low-complexity screen | — |
| **motif discovery** (`discover_motifs`, PWM, PROSITE, degenerate), ORF/gene/promoter, **variant** calling/effect, RNA structure, protein features | → **[`bio-annotation`](../bio-annotation/SKILL.md)** (also uses the Analysis server) |
| pairwise/MSA alignment, alignment identity/similarity | → **[`bio-alignment`](../bio-alignment/SKILL.md)** (owns the Alignment server, incl. the `SequenceAssembler.*` engine tools) |

> **Analysis-server overlap:** the Analysis server is shared with `bio-annotation`. This skill claims
> only the **k-mer** and **repeat/low-complexity** tools for the assembly workflow; every
> motif/variant/ORF/RNA/protein Analysis tool is `bio-annotation`'s — defer there and link it, do not
> re-document it here.

## Decision guide — which tool for which question

| Question | Tool ([MCP] / `Method ID`) · server |
|---|---|
| **Assemble** reads, complex/short overlaps, error-tolerant | `assemble_de_bruijn` / `SequenceAssembler.AssembleDeBruijn` · Alignment |
| **Assemble** reads with clear long overlaps → consensus | `assemble_olc` / `SequenceAssembler.AssembleOLC` · Alignment |
| **Assembly stats / N50** of a contig set | `assembly_stats` / `SequenceAssembler.CalculateStats` · Alignment |
| **Coverage depth** of reads over a reference | `calculate_coverage` / `SequenceAssembler.CalculateCoverage` · Alignment |
| Pairwise / all-pairs **read overlaps** | `find_overlap`,`find_all_overlaps` / `SequenceAssembler.FindOverlap`,`.FindAllOverlaps` · Alignment |
| Merge / scaffold contigs | `merge_contigs`,`scaffold_contigs` / `SequenceAssembler.MergeContigs`,`.Scaffold` · Alignment |
| Pre-assembly read cleanup | `error_correct_reads`,`quality_trim_reads` / `SequenceAssembler.ErrorCorrectReads`,`.QualityTrimReads` · Alignment |
| **k-mer spectrum** (error vs genomic k-mers) | `kmer_spectrum` / `KmerAnalyzer.GetKmerSpectrum` · Analysis |
| **Count k-mers** (one / both strands) | `count_kmers`,`count_kmers_both_strands` / `KmerAnalyzer.CountKmers`,`.CountKmersBothStrands` · Analysis |
| **Most-frequent / min-count / unique** k-mers | `most_frequent_kmers`,`kmers_with_min_count`,`unique_kmers` / `KmerAnalyzer.FindMostFrequentKmers`,`.FindKmersWithMinCount`,`.FindUniqueKmers` · Analysis |
| k-mers **clumping** in a window | `find_clumps` / `KmerAnalyzer.FindClumps` · Analysis |
| Direct / inverted / tandem repeats, palindromes, microsatellites | `find_direct_repeats`,`find_inverted_repeats`,`find_tandem_repeats`,`find_palindromes`,`find_microsatellites` / `RepeatFinder.*`, `GenomicAnalyzer.FindTandemRepeats` · Analysis |
| **Longest repeated** region (suffix-tree, independent path) | `find_longest_repeat`,`suffix_tree_lrs` / `GenomicAnalyzer.FindLongestRepeat`,`SuffixTree.LongestRepeatedSubstring` · Core |
| Low-complexity screen / DUST | `find_low_complexity_regions`,`dust_score`,`mask_low_complexity` / `SequenceComplexity.*` · Analysis |

Rule of thumb: **de Bruijn** for many short reads / dense overlaps; **OLC** for fewer reads with clear
long overlaps. Always follow assembly with **`assembly_stats`** (N50) and screen inputs with the
**k-mer spectrum** and **repeat** tools before trusting a contig.

## Canonical dual-mode pipelines

de Bruijn / OLC share knobs: `kmerSize` (default 31; de-Bruijn operative), `minOverlap` (20),
`minIdentity` (0.9; OLC operative), `minContigLength` (100). Report them in provenance.

### (a) Reads → de Bruijn assembly → assembly stats (N50)
1. **[MCP]** (optional) `error_correct_reads`(reads, k) then `quality_trim_reads` → cleaned reads.
2. **[MCP]** `assemble_de_bruijn`(reads, kmerSize, minContigLength) → `contigs`, `n50`, `longestContig`, `totalLength`, `totalReads`.
3. **[MCP]** `assembly_stats`(contigs, totalReads) → N50 / longest / total (independent recompute; cross-check the engine's own `n50`).
- **[C# API]** `SequenceAssembler.ErrorCorrectReads` → `.AssembleDeBruijn(reads,kmerSize,minOverlap,minIdentity,minContigLength)` → `.CalculateStats(contigs,totalReads)`.
```
Provenance
1) assemble_de_bruijn(reads, kmerSize=31, minContigLength=100) → contigs, n50, totalLength
2) assembly_stats(contigs, totalReads) → n50 (recomputed) — must equal step-1 n50
Cross-check: engine n50 == assembly_stats n50; totalLength == Σ|contig|.
Envelope: none of the 9 guarded units apply. Caveat: alpha — validate before decision use.
```

### (b) Reads → OLC assembly → consensus contigs
1. **[MCP]** (optional) `find_all_overlaps`(reads, minOverlap) → inspect overlap graph before assembling.
2. **[MCP]** `assemble_olc`(reads, minOverlap, minIdentity, minContigLength) → `contigs`, N50 fields.
3. **[MCP]** `assembly_stats`(contigs, totalReads) → N50 report.
- **[C# API]** `SequenceAssembler.FindAllOverlaps` → `.AssembleOLC(reads,minOverlap,minIdentity,…)` → `.CalculateStats`.
- Choose OLC when reads are few with clear long overlaps; de Bruijn (pipeline a) when reads are many/short.

### (c) Coverage of reads over a reference
1. **[MCP]** `calculate_coverage`(reference, reads, minOverlap) → `coverage` (int array, length = |reference|; `coverage[i]` = reads spanning position *i*, **0-based**).
- **[C# API]** `SequenceAssembler.CalculateCoverage(reference, reads, minOverlap)`.
- Mean/min/zero-coverage gaps are derived from the array; report which positions are 0× (assembly gaps).

### (d) k-mer QC for assembly (spectrum / frequent / clumps)
1. **[MCP]** `kmer_spectrum`(sequence, k) → `spectrum` (count→#distinct-k-mers): low-count peak = errors, main peak ≈ coverage. Guides `kmerSize` and a coverage cutoff.
2. **[MCP]** `most_frequent_kmers`(sequence, k[, top]) and/or `kmers_with_min_count`(sequence, k, minCount) → over-represented k-mers (repeats/adapters).
3. **[MCP]** `find_clumps`(sequence, k, windowSize, minOccurrences) → locally clumped k-mers (e.g. ori, repeat arrays). `count_kmers_both_strands` for strand-aware totals.
- **[C# API]** `KmerAnalyzer.GetKmerSpectrum` → `.FindMostFrequentKmers` / `.FindKmersWithMinCount` → `.FindClumps` / `.CountKmersBothStrands`.

### (e) Repeat / low-complexity screen (cross-links `bio-annotation`)
1. **[MCP]** `find_tandem_repeats` / `find_direct_repeats` / `find_inverted_repeats` / `find_microsatellites` / `find_palindromes` → repeat structures that fragment assemblies.
2. **[MCP]** `find_low_complexity_regions` + `dust_score`; `mask_low_complexity` before feeding an assembler.
3. **[MCP]** cross-check the biggest repeat with `find_longest_repeat` / `suffix_tree_lrs` (Core, suffix-tree path).
- **[C# API]** `RepeatFinder.*` / `GenomicAnalyzer.FindTandemRepeats` → `SequenceComplexity.*` → `GenomicAnalyzer.FindLongestRepeat`.
- For **motif discovery, PWM/PROSITE scanning, ORFs, variant/effect** on these regions → use **[`bio-annotation`](../bio-annotation/SKILL.md)** (its Analysis-server subset), not this skill.

## End-to-end grounded example (extends `docs/mcp/README.md`)

**Task.** Given two overlapping reads `AAABBB` and `AABBBC`, (1) assemble them, (2) confirm N50
independently, (3) sanity-check the k-mer basis of the assembly.

Tool / `Method ID` chain (MCP names; C# path in parentheses):
1. `assemble_de_bruijn`(reads=["AAABBB","AABBBC"], kmerSize=3, minContigLength=3)
   → `contigs=["AAABBBBBBC"]`, `n50=10`, `totalLength=10`. (`SequenceAssembler.AssembleDeBruijn`)
2. `assembly_stats`(contigs=["AAABBBBBBC"], totalReads=2) → `n50=10`, `longestContig=10`, `totalLength=10`
   — matches step 1's engine `n50`. (`SequenceAssembler.CalculateStats`)
3. `kmer_spectrum`(sequence="AAABBBBBBC", k=3) → frequency-of-frequencies of the contig's 3-mers;
   a flat low-count spectrum confirms no spurious high-multiplicity k-mer. (`KmerAnalyzer.GetKmerSpectrum`)

Expected-shape output (values illustrative; **compute them with the tools, do not eyeball**):
```
| contig      | length | n50 (engine) | n50 (stats) | totalLength |
|-------------|-------:|-------------:|------------:|------------:|
| AAABBBBBBC  |     10 |           10 |          10 |          10 |

Provenance
1) assemble_de_bruijn(reads=[…], kmerSize=3, minContigLength=3) → contigs, n50=10, totalLength=10
2) assembly_stats(contigs, totalReads=2) → n50=10 (independent recompute) — agrees with step 1
3) kmer_spectrum(contig, k=3) → spectrum; no anomalous high-count k-mer
Cross-check: engine n50 == assembly_stats n50; totalLength == Σ|contig|.
Envelope: none of the 9 LimitationPolicy-guarded units apply to assembly/k-mer/repeat tools.
Caveat: alpha software; not for clinical use — independently validate before relying on a contig.
```

## Reference

- **Full domain tool index (analysis + core union, generated — do NOT hand-edit):** [`_generated/tools.md`](_generated/tools.md)
  (produced by `scripts/skills/gen-catalog.py`; if absent, run `seqeron-discovery`). Assembly-engine
  tools (`SequenceAssembler.*`) appear in the **Alignment** slice — see `bio-alignment`.
- **Assembly-relevant recipes + parameter guidance:** [`reference/pipelines.md`](reference/pipelines.md)
- **Tool map (assembly / coverage&stats / k-mer / repeats):** [`reference/tool-map.md`](reference/tool-map.md)
- **Algorithm background (invariants/formulas — link, don't copy):**
  [`docs/algorithms/Assembly/De_Bruijn_Graph_Assembly.md`](../../../docs/algorithms/Assembly/De_Bruijn_Graph_Assembly.md) ·
  [`Overlap_Layout_Consensus.md`](../../../docs/algorithms/Assembly/Overlap_Layout_Consensus.md) ·
  [`Assembly_Statistics.md`](../../../docs/algorithms/Assembly/Assembly_Statistics.md) ·
  [`Coverage_Calculation.md`](../../../docs/algorithms/Assembly/Coverage_Calculation.md) ·
  [`docs/algorithms/K-mer_Analysis/K-mer_Frequency_Analysis.md`](../../../docs/algorithms/K-mer_Analysis/K-mer_Frequency_Analysis.md) ·
  [`docs/algorithms/Repeat_Analysis/Repeat_Detection.md`](../../../docs/algorithms/Repeat_Analysis/Repeat_Detection.md)
- **Cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor guardrail) · [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup) · [`bio-annotation`](../bio-annotation/SKILL.md) (motif/variant/ORF Analysis-server subset) · [`bio-alignment`](../bio-alignment/SKILL.md) (owns the Alignment `SequenceAssembler.*` engine).
