---
name: bio-qc
description: >-
  Parse and QC biological sequences and files with Seqeron (MCP tools OR the
  C# API). Use to parse/validate FASTA, FASTQ, GenBank, GFF, VCF, BED, EMBL;
  validate DNA/RNA/protein; compute sequence composition (GC%, base/AA counts,
  Tm/melting temperature, molecular weight, hydrophobicity, pI); measure
  complexity/entropy (Shannon, DUST, linguistic, k-mer, compression); FASTQ
  quality stats (Q20/Q30, mean quality, encoding); and transcribe / translate /
  reverse-complement. Prompts like "parse this file", "extract features from
  GenBank/GFF", "validate this sequence", "GC content of…", "is this valid
  DNA/RNA/protein", "quality stats for this FASTQ", "translate this ORF",
  "reverse complement", "mask low-complexity". Servers: sequence + parsers.
allowed-tools: Read, Bash, Grep, Glob
---

# bio-qc — parse/validate files & sequences, composition, complexity, quality

Routing + orchestration skill for the **Sequence** (35 tools) and **Parsers** (41 tools) servers
= **76 tools**. It picks the right tool for a parsing / validation / composition / quality question
and gives a **dual-mode** recipe (MCP tool calls **and** the equivalent `Seqeron.Genomics` C# `Method ID`s).

- **Rigor is delegated.** Parse-with-a-tool (never interpret FASTA/FASTQ by hand), envelope,
  provenance, cross-check, units, and the alpha / not-for-clinical-use caveat are owned by
  **`bio-rigor`** — applies here by default; do not restate its rules.
- **Don't know the tool name?** Use **`seqeron-discovery`**
  (`python3 scripts/skills/find-tool.py <kw> --server sequence|parsers`) — never guess.
- **Point, don't duplicate.** Full I/O schemas live in `docs/mcp/tools/{sequence,parsers}/*.md`;
  algorithm invariants in `docs/algorithms/{Sequence_Composition,Complexity,Quality,Translation,FileIO,Statistics}/`.
  This skill links, it does not copy.

## Envelope STOP rule (one guarded unit in this domain)

- **PARSE-FASTQ-001** (`fastq_detect_encoding` / `FastqParser.DetectEncoding`): auto-disambiguation of
  Phred+33 vs Phred+64 is **guarded → MinimumMode `Permissive`**. When *every* read is confined to the
  ASCII 64–74 overlap the encoding is information-theoretically ambiguous; the call flags `Ambiguous`
  and defaults to Phred+33. **STOP** and surface the ambiguity — do not silently pick an encoding.
  In C# tests this branch needs the **Permissive bootstrap** (see `bio-rigor` → `reference/envelope.md`).

## Decision guide — which tool for which question

| Question | Tool ([MCP] / `Method ID`) |
|---|---|
| Parse a **FASTA / FASTQ / GenBank / GFF / VCF / BED / EMBL** file | `fasta_parse` / `fastq_parse` / `genbank_parse` / `gff_parse` / `vcf_parse` / `bed_parse` / `embl_parse` |
| **Is this valid** DNA / RNA / protein? | `dna_validate` / `rna_validate` / `protein_validate` (`*.TryCreate`); quick bool: `is_valid_dna` / `is_valid_rna` |
| **GC content** of a sequence | `gc_content` / `SequenceExtensions.CalculateGcContentFast` |
| **Base / amino-acid composition** | `nucleotide_composition` / `amino_acid_composition` |
| **One-shot summary** (length, GC, counts, Tm…) | `summarize_sequence` / `SequenceStatistics.SummarizeNucleotideSequence` |
| **Melting temperature (Tm)** | `melting_temperature` / `SequenceStatistics.CalculateMeltingTemperature` |
| **Molecular weight / hydrophobicity / pI** | `molecular_weight_nucleotide` · `molecular_weight_protein` · `hydrophobicity` · `isoelectric_point` |
| **Complexity / entropy** (screen low-complexity) | `complexity_shannon` · `complexity_dust_score` · `complexity_linguistic` · `complexity_kmer_entropy` · `complexity_mask_low` |
| **FASTQ quality stats** (Q20/Q30, mean Q, encoding) | `fastq_statistics` · `fastq_detect_encoding` (guarded — see STOP) |
| **Transcribe / translate / reverse-complement** | `rna_from_dna` · `translate_dna` / `translate_rna` · `dna_reverse_complement` |
| **Extract features** from GenBank / GFF / EMBL | `genbank_features` · `gff_parse`(+`gff_filter`) · `embl_features` |

Rule of thumb: **parse first** (tool, not by hand) → **validate** the alphabet → then **measure**
(composition / complexity / quality). See `reference/tool-map.md` for all 76 tools by sub-task.

## Canonical dual-mode pipelines

### (a) Parse FASTA → validate → composition summary
1. **[MCP]** `fasta_parse`(content) → records (id, sequence).
2. **[MCP]** `dna_validate`(sequence) → `valid`, `length` (stop if `valid=false`).
3. **[MCP]** `summarize_sequence`(sequence) or `gc_content` + `nucleotide_composition`.
- **[C# API]** `FastaParser.Parse` → `DnaSequence.TryCreate` → `SequenceStatistics.SummarizeNucleotideSequence` (or `SequenceExtensions.CalculateGcContentFast` + `SequenceStatistics.CalculateNucleotideComposition`).
```
Provenance: fasta_parse → dna_validate → summarize_sequence. Envelope: none guarded here. Caveat: alpha.
```

### (b) FASTQ → quality stats
1. **[MCP]** `fastq_detect_encoding`(reads) → encoding (**guarded**; if `Ambiguous`, STOP & report — do not guess).
2. **[MCP]** `fastq_statistics`(content, encoding) → `totalReads`, `meanQuality`, `q20Percentage`, `q30Percentage`, `gcContent`.
- **[C# API]** `FastqParser.DetectEncoding` (Permissive bootstrap) → `FastqParser.CalculateStatistics`.
```
Provenance: fastq_detect_encoding(encoding=…) → fastq_statistics. Envelope: PARSE-FASTQ-001 (Permissive if ambiguous). Caveat: alpha.
```

### (c) Parse GenBank / GFF → extract features
1. **[MCP]** `genbank_parse`(content) → record; `genbank_features` → feature table; `genbank_extract_sequence` → per-feature seq.
   For GFF: `gff_parse` → features; `gff_filter`(type=…) → subset; `gff_statistics` → counts.
- **[C# API]** `GenBankParser.Parse` → `GenBankParser.GetFeatures` → `GenBankParser.ExtractSequence`; `GffParser.Parse` → `GffParser.Filter*` → `GffParser.CalculateStatistics`.
```
Provenance: genbank_parse → genbank_features (→ genbank_extract_sequence). Envelope: none guarded. Caveat: alpha.
```

### (d) Transcribe / translate / reverse-complement
1. **[MCP]** `rna_from_dna`(dna) → RNA; `translate_dna`(dna) or `translate_rna`(rna) → protein; `dna_reverse_complement`(dna) → revcomp.
- **[C# API]** `RnaSequence.FromDna` → `Translator.Translate(DnaSequence|RnaSequence)`; `DnaSequence.GetReverseComplementString`.
- **Cross-check:** `translate_dna(dna)` should equal `translate_rna(rna_from_dna(dna))` (independent paths).
```
Provenance: (rna_from_dna →) translate_dna/translate_rna; dna_reverse_complement. Envelope: none guarded. Caveat: alpha.
```

### (e) Low-complexity / entropy screen
1. **[MCP]** `complexity_shannon`(seq) and/or `complexity_dust_score`(seq) → complexity metrics; `complexity_mask_low`(seq) → masked seq.
- **[C# API]** `SequenceComplexity.CalculateShannonEntropy` · `.CalculateDustScore` · `.MaskLowComplexity`.
- Use before k-mer/alignment work to flag repetitive / low-information regions.
```
Provenance: complexity_dust_score + complexity_shannon (→ complexity_mask_low). Envelope: none guarded. Caveat: alpha.
```

## End-to-end grounded example (extends `docs/mcp/README.md` cloning-QC)

**Task.** A cloning insert is provided in FASTA. (1) Parse it with a tool (never by hand),
(2) validate it is DNA, (3) report GC% and composition, (4) get its melting temperature,
(5) screen for low complexity — one QC report.

Tool / `Method ID` chain (MCP names; C# path in parentheses):
1. `fasta_parse`(content=">seq1\nGCGCGAATTCATGGATCCATAT") → `sequence`="GCGCGAATTCATGGATCCATAT". (`FastaParser.Parse`)
2. `dna_validate`("GCGCGAATTCATGGATCCATAT") → `valid`=true, `length`=22. (`DnaSequence.TryCreate`)
3. `gc_content`(sequence) → `gcContent`≈45.45, `gcCount`=10, `totalCount`=22. (`SequenceExtensions.CalculateGcContentFast`)
4. `melting_temperature`(sequence) → `tm` (report the method/params from the tool doc). (`SequenceStatistics.CalculateMeltingTemperature`)
5. `complexity_dust_score`(sequence) → DUST score (flag if low-complexity). (`SequenceComplexity.CalculateDustScore`)

Expected-shape output (**compute with the tools, do not eyeball**):
```
| id   | length | valid | gc_% | tm_C | dust |
|------|-------:|:-----:|-----:|-----:|-----:|
| seq1 |     22 |  yes  | 45.45|   …  |   …  |

Provenance
1) fasta_parse(content) → sequence            (parsed by tool, not manually)
2) dna_validate(sequence) → valid=true, length=22
3) gc_content(sequence) → gcContent=45.45, gcCount=10, totalCount=22
4) melting_temperature(sequence) → tm (see tool doc for method + salt/conc params)
5) complexity_dust_score(sequence) → dust
Envelope: none guarded in this chain (FASTQ encoding would be PARSE-FASTQ-001). 
Caveat: alpha software; not for clinical use — validate before relying on any construct decision.
```

## Reference

- **Full domain tool index (all 76, generated — do NOT hand-edit):** [`_generated/tools.md`](_generated/tools.md)
  (produced by `scripts/skills/gen-catalog.py`; if absent, run `seqeron-discovery`).
- **Fuller recipes + parameter guidance:** [`reference/pipelines.md`](reference/pipelines.md)
- **Tool map (76 tools by sub-task, one-liners + Method ID):** [`reference/tool-map.md`](reference/tool-map.md)
- **Algorithm background (invariants/formulas — link, don't copy):**
  [`docs/algorithms/FileIO/`](../../../docs/algorithms/FileIO/) (FASTA/FASTQ/GenBank/GFF/VCF/BED/EMBL parsing) ·
  [`Sequence_Composition/Sequence_Composition.md`](../../../docs/algorithms/Sequence_Composition/Sequence_Composition.md) ·
  [`Sequence_Composition/Sequence_Validation.md`](../../../docs/algorithms/Sequence_Composition/Sequence_Validation.md) ·
  [`Statistics/Melting_Temperature.md`](../../../docs/algorithms/Statistics/Melting_Temperature.md) ·
  [`Complexity/`](../../../docs/algorithms/Complexity/) (DUST/entropy/Lempel-Ziv) ·
  [`Quality/`](../../../docs/algorithms/Quality/) (Phred / quality stats) ·
  [`Translation/`](../../../docs/algorithms/Translation/) (codon/protein/six-frame).
- **Cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor guardrail) · [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup).
