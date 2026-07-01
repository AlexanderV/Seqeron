---
name: bio-rigor
version: 1.0.0
description: Enforce scientific rigor whenever you compute or derive a biological result from real data with the Seqeron library — via MCP tools OR the C# API. Use whenever a task involves parsing FASTA/FASTQ/GenBank/GFF/VCF/BED, computing GC%, Tm, distances, k-mers, alignments, variant calls, primers/CRISPR, phylogeny, or popgen metrics, or when asked "is this result reliable/valid?". Enforces tool-only computation (no manual parsing or mental math), a validated operating envelope (LimitationPolicy / LIMITATIONS.md), a provenance chain, cross-checking, unit/coordinate checks, and the alpha / not-for-clinical-use disclaimer.
allowed-tools: Read, Grep, Glob
triggers: [
  "GC content", "GC%", "melting temperature", "Tm", "compute", "calculate",
  "parse FASTA", "parse FASTQ", "parse VCF", "parse GFF", "parse GenBank", "parse BED",
  "call variants", "variant calling", "align", "alignment", "distance", "edit distance",
  "k-mer", "kmer", "motif", "ORF", "gene", "promoter", "primer", "probe", "CRISPR guide",
  "codon optimization", "restriction site", "phylogenetic tree", "popgen", "Fst",
  "taxonomic classification", "GC skew", "structural variant",
  "is this reliable", "is this valid", "is this result correct", "can I trust",
  "provenance", "reproducible", "operating envelope", "LimitationPolicy", "validated"
]
---

# bio-rigor — scientific-rigor guardrail for Seqeron

This is a **cross-cutting guardrail**. It applies whenever you turn real biological input into a
computed result using Seqeron — through **MCP tools** (LLM client) or the **C# `Seqeron.Genomics` API**.
It does not replace the domain skills; it constrains *how* they produce answers.

Seqeron is **PRE-RELEASE / ALPHA and NOT certified for clinical or diagnostic use**
(see the disclaimer atop [`README.md`](../../../README.md)). Enforce the six rules below.

## When this applies

- Parsing any bio format (FASTA/FASTQ/GenBank/GFF/VCF/BED/EMBL) into sequences/records.
- Computing any metric: GC%, Tm, composition, complexity, distances, counts, positions,
  alignment scores, variant calls, primer/CRISPR properties, tree topology, popgen stats.
- Any question of the form "is this result reliable / valid / trustworthy?".

## When this does NOT apply

- Pure library development (writing/refactoring/testing Seqeron code) — that is
  `clean-architecture` / `clean-code` / the testing campaigns.
- General discussion of biology with no computation over supplied data.

---

## The six rules

### 1. Tool-only computation
- **Parse with a parser, never by eye.** Read FASTA/FASTQ/VCF/GFF/GenBank via the parser tool
  (`fasta_parse`, `fastq_parse`, `vcf_parse`, …) or the C# parser API. Never hand-interpret headers,
  sequences, coordinates, or quality strings.
- **Compute every metric with a call.** GC%, Tm, distances, counts, positions, alignment scores —
  each comes from a tool/`Method ID`, never from mental math or ad-hoc scratch code.
- **No improvised calculation.** If no tool/API exists for a step, **say so and stop** — do not
  invent a formula, spreadsheet, or one-off script. Use `seqeron-discovery` to confirm a tool exists.
- This is the discipline the MCP examples repeat verbatim: *"Use tools only; no manual parsing or
  calculations."* See [`docs/mcp/README.md`](../../../docs/mcp/README.md).

### 2. Envelope awareness
- Before using a **guarded** algorithm unit, confirm the task is inside its **validated operating
  envelope**. The guard is `Seqeron.Genomics.Core.LimitationPolicy` (3 tiers:
  `Strict` < `Moderate` < `Permissive`, default **`Moderate`**).
- If the task pushes a unit **below its MinimumMode**, the guarded call throws
  `SeqeronLimitationException`. **STOP and report the limitation** (name it, say why, name the
  alternative) — do **not** raise the mode just to force an output.
- The 9 guarded units, their MinimumMode, and the honest scope boundaries live in
  [`docs/Validation/LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md).
  Details on checking the envelope + the **Permissive test bootstrap**: [`reference/envelope.md`](reference/envelope.md).

### 3. Provenance
- **Every result carries a provenance block**: the tools / `Method ID`s used, **in call order**,
  with the exact parameters. A result without provenance is not reportable.
- Exact block format + a worked example: [`reference/provenance-format.md`](reference/provenance-format.md).

### 4. Cross-checking
- **Corroborate any critical / decision-relevant result** via a second independent tool or an
  invariant, where one exists. Examples: verify a motif position by re-searching the
  reverse-complement; confirm a count via an independent counter; check an alignment score against
  its identity%; validate GC% against explicit A/C/G/T composition.
- If no independent check exists, **state that** in the provenance rather than implying certainty.

### 5. Units & coordinates
- Before reporting, verify: **coordinate base** (0-based vs 1-based — Seqeron tool outputs are
  0-based unless a tool doc states otherwise), **units** (°C for Tm, % for GC, bp for lengths),
  and **alphabet/validity** (validate DNA/RNA/protein with the validator tool before deriving
  anything from it).
- Report units and coordinate base explicitly in the output.

### 6. Disclaimer
- If a result could be read as **clinical, diagnostic, or decision-relevant** (e.g. variant
  pathogenicity, HLA typing, tumour purity, primer/guide for a real assay), surface the
  **alpha / not-for-clinical-use** caveat and the "independently validate before relying on it"
  note from [`README.md`](../../../README.md).

---

## Dual-mode quick reference

| Rule | MCP client | C# API caller |
|---|---|---|
| Parse | call `*_parse` tool | call the parser API, not string ops |
| Compute | one tool call per metric | one `Method ID` call per metric |
| Envelope | read tool doc; stop on guarded-error | check `LimitationPolicy`; catch `SeqeronLimitationException`; `Permissive` bootstrap in tests |
| Provenance | list tools + args in order | list `Method ID`s + args in order |
| Cross-check | second independent tool | second independent method/invariant |

**Bottom line:** parse with a tool, compute with a tool, stay inside the envelope, show your
provenance, cross-check what matters, and never let a computed number look more certain — or more
clinical — than it is.
