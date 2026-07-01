# bio-qc pipelines — fuller recipes, parameter guidance, gotchas

Detailed, on-demand companion to `SKILL.md`. Every recipe is **dual-mode** ([MCP] tool names +
`Seqeron.Genomics` C# `Method ID`s) and ends with a **Provenance** block. Rigor rules
(parse-with-a-tool, envelope, cross-check, units, alpha caveat) belong to **`bio-rigor`** — not repeated here.

Always: **parse the file with a tool** (never read FASTA/FASTQ/GenBank by hand), then **validate the
alphabet**, then **measure**. Report every tool + params in provenance so the result is reproducible.

---

## 1. FASTA → validate → composition summary

**Goal.** Given a FASTA blob or file, get a clean per-record QC row (length, valid?, GC%, composition).

1. **[MCP]** `fasta_parse`(content) → `[{ id, description, sequence }]`.
2. For each record: **[MCP]** `dna_validate`(sequence) → `valid`, `length`, `error`.
   - If `valid=false`: **STOP for that record**, report the `error`; do not compute downstream metrics on it.
3. **[MCP]** `summarize_sequence`(sequence) → one-shot (length, GC, composition, Tm…), **or** compose
   `gc_content` + `nucleotide_composition` if you only need those.
- **[C# API]** `FastaParser.Parse` → `DnaSequence.TryCreate(out seq, out error)` →
  `SequenceStatistics.SummarizeNucleotideSequence(seq)` (or `SequenceExtensions.CalculateGcContentFast` +
  `SequenceStatistics.CalculateNucleotideComposition`).

**Gotchas.** `gc_content.gcContent` is a **percentage** (e.g. 45.45), while `fastq_statistics.gcContent`
is a **fraction** (0.0–1.0) — do not mix them. Use `is_valid_dna` for a fast boolean gate;
use `dna_validate`/`TryCreate` when you need the failure reason and length.

```
Provenance
1) fasta_parse(content) → records                (parsed by tool, not manually)
2) dna_validate(seq) → valid, length, error
3) summarize_sequence(seq) → length, gc, composition, (tm)
Envelope: none guarded. Caveat: alpha — validate before decisions.
```

---

## 2. FASTQ → encoding → quality stats (guarded)

**Goal.** QC a sequencing run: mean quality, Q20/Q30, GC, read-length distribution.

1. **[MCP]** `fastq_detect_encoding`(reads) → encoding. **Guarded: PARSE-FASTQ-001.**
   - If the result is **`Ambiguous`** (every read confined to ASCII 64–74, the Phred+33/+64 overlap):
     **STOP and surface it.** The tool defaults to Phred+33, but this is information-theoretically
     unresolvable — ask the user or use externally-known encoding. In C# this branch needs the
     **Permissive bootstrap** (`bio-rigor` → `reference/envelope.md`).
2. **[MCP]** `fastq_statistics`(content, encoding) → `totalReads`, `totalBases`, `meanReadLength`,
   `meanQuality`, `minReadLength`, `maxReadLength`, `q20Percentage`, `q30Percentage`, `gcContent` (fraction).
3. Optional cleanup before re-stats: `fastq_filter` (drop low-mean-Q reads), `fastq_trim_quality`
   (trim low-Q ends), `fastq_trim_adapter` (remove adapters).
- **[C# API]** `FastqParser.DetectEncoding(IEnumerable<string>)` → `FastqParser.CalculateStatistics(content, encoding)`
  → `FastqParser.FilterByQuality` / `TrimByQuality` / `TrimAdapter`.

**Gotchas.** Pass an **explicit** `encoding` once known — don't re-detect per call. `meanQuality` is a
Phred score; convert to/from error probabilities with `fastq_phred_to_error` / `fastq_error_to_phred`.

```
Provenance
1) fastq_detect_encoding(reads) → encoding=phred33|phred64|Ambiguous(→STOP)
2) fastq_statistics(content, encoding) → meanQuality, q20%, q30%, gcContent(fraction)
Envelope: PARSE-FASTQ-001 (MinimumMode Permissive if Ambiguous). Caveat: alpha.
```

---

## 3. GenBank / GFF / EMBL → extract & count features

**Goal.** Pull a feature table (CDS, gene, exon…) and the sequence behind each feature.

- **GenBank:** `genbank_parse`(content) → record → `genbank_features` → feature list
  (type, location, qualifiers). For a feature's bases: `genbank_extract_sequence`(record, location);
  parse a raw location string with `genbank_parse_location`; overall counts with `genbank_statistics`.
- **GFF:** `gff_parse`(content) → features → `gff_filter`(type="CDS"|…) → subset → `gff_statistics` → counts.
- **EMBL:** `embl_parse` → `embl_features` → `embl_statistics`.
- **[C# API]** `GenBankParser.Parse` → `GetFeatures` → `ExtractSequence` / `ParseLocation` / `Statistics`;
  `GffParser.Parse` → `Filter*` → `CalculateStatistics`; `EmblParser.Parse` → `GetFeatures` → `Statistics`.

**Gotchas.** GenBank locations can be `join(...)`, `complement(...)`, and are **1-based inclusive** in the
file — let `genbank_parse_location` / `ExtractSequence` handle the arithmetic; do not slice by hand.
**BED is 0-based half-open** — do not compare BED coords directly against GenBank/GFF 1-based coords
without conversion (see `bio-rigor` on coordinate systems).

```
Provenance
1) genbank_parse(content) → record
2) genbank_features(record) → features[]
3) genbank_extract_sequence(record, feature.location) → sub-seq  (coords handled by tool)
Envelope: none guarded. Caveat: alpha.
```

---

## 4. Transcribe / translate / reverse-complement

**Goal.** Move between strands and alphabets; get the protein of an ORF.

1. **[MCP]** `rna_from_dna`(dna) → RNA (T→U).
2. **[MCP]** `translate_dna`(dna) **or** `translate_rna`(rna) → protein (report the genetic-code table used;
   see the tool doc / `docs/algorithms/Translation/Codon_Translation.md`).
3. **[MCP]** `dna_reverse_complement`(dna) → reverse complement (for the antisense frame / primer design).
- **[C# API]** `RnaSequence.FromDna` → `Translator.Translate(DnaSequence|RnaSequence)`;
  `DnaSequence.GetReverseComplementString`; single base via `SequenceExtensions.GetComplementBase`.
- **Cross-check (two independent paths):** `translate_dna(dna)` should equal
  `translate_rna(rna_from_dna(dna))`. Disagreement ⇒ investigate before trusting either.

**Gotchas.** Translation is **frame-0** on the given strand — supply the correct reading frame / strand
yourself (reverse-complement first for antisense). For ORF finding / six-frame scanning use **`bio-annotation`**.

```
Provenance
1) rna_from_dna(dna) → rna
2) translate_dna(dna) → protein   (genetic code per tool doc)
Cross-check: translate_rna(rna) == translate_dna(dna)
Envelope: none guarded. Caveat: alpha.
```

---

## 5. Low-complexity / entropy screen

**Goal.** Flag repetitive or low-information regions before k-mer/alignment/primer work.

1. **[MCP]** `complexity_dust_score`(seq) → DUST score (higher = lower complexity).
2. **[MCP]** `complexity_shannon`(seq) and/or `complexity_linguistic`(seq) / `complexity_kmer_entropy`(seq)
   / `complexity_compression_ratio`(seq) → corroborating metrics.
3. **[MCP]** `complexity_mask_low`(seq) → sequence with low-complexity stretches masked (soft/hard).
- **[C# API]** `SequenceComplexity.CalculateDustScore` · `CalculateShannonEntropy` ·
  `CalculateLinguisticComplexity` · `CalculateKmerEntropy` · `EstimateCompressionRatio` · `MaskLowComplexity`.

**Gotchas.** There are **two** Shannon/linguistic implementations — `SequenceComplexity.*` (complexity module)
and `SequenceStatistics.*` (stats module). They can differ in windowing/normalization; pick one path and
keep it consistent across a comparison. Screen with a couple of metrics (DUST + one entropy) rather than one.

```
Provenance
1) complexity_dust_score(seq) → dust
2) complexity_shannon(seq) → entropy
3) complexity_mask_low(seq) → masked_seq
Envelope: none guarded. Caveat: alpha.
```

---

## Cross-references

- **Rigor guardrail (always on):** [`bio-rigor`](../../bio-rigor/SKILL.md) — parse-with-a-tool, envelope,
  provenance, cross-check, units/coords, alpha caveat, Permissive test bootstrap.
- **Tool lookup:** [`seqeron-discovery`](../../seqeron-discovery/SKILL.md) —
  `python3 scripts/skills/find-tool.py <kw> --server sequence|parsers`.
- **All 76 tools by sub-task:** [`tool-map.md`](tool-map.md).
- **ORF finding / six-frame / variant effect:** `bio-annotation`. **k-mer assembly:** `bio-assembly`.
- **Algorithm background (link, don't copy):**
  [`docs/algorithms/FileIO/`](../../../../docs/algorithms/FileIO) ·
  [`Sequence_Composition/`](../../../../docs/algorithms/Sequence_Composition) ·
  [`Complexity/`](../../../../docs/algorithms/Complexity) ·
  [`Quality/`](../../../../docs/algorithms/Quality) ·
  [`Translation/`](../../../../docs/algorithms/Translation) ·
  [`Statistics/`](../../../../docs/algorithms/Statistics).
```
