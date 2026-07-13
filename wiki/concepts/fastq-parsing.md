---
type: concept
title: "FASTQ parsing (4-line records; tolerant assembly + quality/filter/trim/paired-end helpers)"
tags: [file-io, algorithm]
sources:
  - docs/algorithms/FileIO/FASTQ_Parsing.md
source_commit: bfe0725d505fc6e6ee78ee0fe2c4ed25704ca9cd
created: 2026-07-13
updated: 2026-07-13
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: parse-fastq-001-evidence
      evidence: "Test Unit ID: PARSE-FASTQ-001 ... Algorithm Group: FileIO (FastqParser)"
      confidence: high
      status: current
---

# FASTQ parsing (4-line records; tolerant assembly + quality/filter/trim/paired-end helpers)

**FASTQ** pairs each nucleotide sequence with a per-base **quality string**. A record is the
canonical **four logical lines**: an `@`-**header**, the **sequence**, a `+` **separator**, and a
**quality string of the same length as the sequence**. Seqeron's `FastqParser` (test unit
**PARSE-FASTQ-001**, status *Simplified*) parses such content from strings, files, and readers into
typed `FastqRecord` values, decodes Phred quality, and layers filtering, trimming, statistics,
writing, and paired-end helpers on top.

This concept is the **parser-surface** synthesis of the primary algorithm spec
`docs/algorithms/FileIO/FASTQ_Parsing.md` — it owns the record state machine and the `FastqParser`
contract. Two aspects with real algorithmic content are factored into sibling concepts and are
**consumed, not re-explained** here: the quality-score **encoding / offset math** lives in
[[phred-quality-encoding]] and the run-QC **summary statistics** live in
[[fastq-quality-statistics]]; quality-based end trimming is elaborated by
[[quality-trimming-running-sum]]. This is a member of the FileIO file-parsing **`PARSE-*`** family
anchored by [[bed-format-parsing]]; its simpler format sibling (no quality line) is
[[fasta-parsing]]. The literature-traced format facts and test oracles live in the Evidence source
page [[parse-fastq-001-evidence]]; [[test-unit-registry]] tracks the unit,
[[algorithm-validation-evidence]] describes the evidence-artifact pattern, and [[fuzzing]] explains
why parsers are the campaign's highest-priority malformed-input target.

## The record (what the parser consumes)

```text
@<identifier> <optional description>
<sequence>
+
<quality string>
```

- **Header:** begins with `@`; split at the **first space** into `Id` (identifier) and the optional
  `Description` (remainder). The `+` separator may optionally repeat the identifier.
- **Core invariant (INV-01):** in a valid record, **quality-string length == sequence length** —
  FASTQ encodes exactly one quality score per base.
- **Quality (INV-02/INV-03):** Phred scores in printable ASCII, either **Phred+33** (Q 0–93) or
  **Phred+64** (Q 0–62); the score math `Q = −10·log₁₀(p)` / `p = 10^(−Q/10)`, ASCII ranges,
  boundary characters, and offset auto-detection are the contract of [[phred-quality-encoding]].

## The parser as a (tolerant) state machine

`Parse` is a **line-oriented state machine** that skips blank lines and only **starts a record when
a line begins with `@`**. It then:

1. Splits the header into `Id` and optional `Description` at the first space.
2. **Accumulates sequence lines** until a line beginning with `+` is encountered (multi-line
   sequences are accepted).
3. **Accumulates quality lines** until the accumulated quality length **reaches the sequence
   length** (this length-driven stop is how multi-line quality is reassembled).
4. Detects or applies the selected encoding, decodes the quality string to Phred scores, and yields
   the `FastqRecord`.

Parsing is deliberately **tolerant**: it favors record assembly over strict validation and performs
**no malformed-record rejection pass** beyond these rules — so a malformed record may be skipped or
partially assembled rather than raising an error. Null/empty input returns **no records**.

## Contract & surface (`FastqParser`, `FastqRecord`, `FastqStatistics`)

Implementation: `src/Seqeron/Algorithms/Seqeron.Genomics.IO/FastqParser.cs`.

| Group | Entry points | Notes |
|-------|--------------|-------|
| **Parse** | `Parse(string, QualityEncoding)`, `ParseFile(string, …)`, `Parse(TextReader, …)` | String / file / reader; `encoding` defaults to **`Auto`** (heuristic detector). |
| **Encoding** | `DetectEncoding`, `DecodeQualityScores`, `EncodeQualityScores` | `Decode` subtracts the ASCII offset and **clamps negatives to 0**; `Encode` **clamps to the encoding's representable range**. Contract: [[phred-quality-encoding]]. |
| **Phred math** | `PhredToErrorProbability`, `ErrorProbabilityToPhred` | `ErrorProbabilityToPhred` returns **93** when `p ≤ 0` (Sanger Phred+33 max — a representability bound). |
| **Filter** | `FilterByQuality(minAverageQuality)`, `FilterByLength` | Average-quality / length-range filters over decoded scores. |
| **Trim** | `TrimByQuality(minQuality=20)`, `TrimAdapter(adapter, minOverlap)` | End trimming (details: [[quality-trimming-running-sum]]) and adapter removal. |
| **Statistics** | `CalculateStatistics`, `CalculatePositionQuality` | Aggregate QC scalars; surface: [[fastq-quality-statistics]]. |
| **Write / pair** | `WriteToStream`, `ToFastqString`, `InterleavePairedReads`, `SplitInterleavedReads` | Serialization + paired-end interleave/split. |

`FastqRecord` exposes `Id`, `Description`, `Sequence`, `QualityString`, and decoded
`QualityScores` (`IReadOnlyList<int>`). `FastqStatistics` exposes `TotalReads`, `TotalBases`,
`MeanReadLength`, `MeanQuality`, `Q20Percentage` / `Q30Percentage`, and `GcContent`.

### Complexity

| Operation | Time | Space |
|-----------|------|-------|
| `Parse` / `ParseFile` | `O(n)` | `O(1)` aux |
| `DetectEncoding` | `O(m)` (m = quality length) | `O(1)` |
| `DecodeQualityScores` / `EncodeQualityScores` | `O(m)` | `O(m)` |
| `FilterByQuality` / `CalculateStatistics` | `O(r·m)` (r = records) | `O(1)` aux |
| `TrimByQuality` | `O(m)` | `O(m)` |
| `TrimAdapter` | `O(m·a)` worst case (a = adapter length) | `O(m)` |

## Paired-end support

Two forms are modeled: **separate files** (R1 and R2 in matched order) and **interleaved** (R1/R2
alternate within one stream). `InterleavePairedReads` / `SplitInterleavedReads` convert between the
two. The splitter **alternates records without pair validation**, so an odd (unmatched) final record
simply stays on whichever alternating side the splitter reaches.

## Edge cases & intentional simplifications

- **Null / empty input** → no records (explicit early-return guards); `DecodeQualityScores` returns
  an **empty array** for null/empty quality strings.
- **Empty quality string** → detected as Phred+33, decodes to an empty score list.
- **All-low-quality read** → `TrimByQuality` can return an **empty-sequence** record (end trimming
  may remove the whole read); it trims **ends only**, not internal low-quality segments.
- **No adapter match / short adapter** → `TrimAdapter` returns the record unchanged (no trimming
  when the adapter is null, empty, or shorter than `minOverlap`). It first looks for an adapter
  **overlap at the 3′ end**, then a **full internal match starting after position 0** — a full
  adapter already at the first base is left unless the end-overlap path trims it.
- **Auto-detection is heuristic** (see [[phred-quality-encoding]]): the ambiguous `@`–`I` window
  defaults to Phred+33, and a high-quality Phred+33 string containing `J`–`~` can be **misclassified
  as Phred+64** unless the caller passes the encoding explicitly.
- **`GcContent` unit mismatch:** `CalculateStatistics` reports `Q20Percentage`/`Q30Percentage` as
  percentages but **`GcContent` as a 0..1 fraction** — a documented API inconsistency, not a spec
  rule.
- **No strict validator:** the repository provides practical parsing + quality utilities, not a
  conformance validator; encoding detection cannot reliably disambiguate the full Phred+33/Phred+64
  printable-range overlap.
