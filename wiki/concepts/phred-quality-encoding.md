---
type: concept
title: "Phred quality-score encoding (Phred+33 / Phred+64)"
tags: [file-io, algorithm]
sources:
  - docs/Evidence/PARSE-FASTQ-001-Evidence.md
  - docs/Evidence/QUALITY-PHRED-001-Evidence.md
  - docs/algorithms/Quality/Phred_Score_Handling.md
source_commit: 0677fe25ab89d866e2e61ba7ad892fafa8708fec
created: 2026-07-10
updated: 2026-07-16
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: parse-fastq-001-evidence
      evidence: "Test Unit ID: PARSE-FASTQ-001 ... Area: FileIO ... FASTQ Parsing ... Quality encoding: Phred+33 and Phred+64"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: quality-phred-001-evidence
      evidence: "Test Unit ID: QUALITY-PHRED-001 ... Algorithm: Phred Score Handling (FASTQ quality string parsing, encoding, and Phred+33 ↔ Phred+64 conversion)"
      confidence: high
      status: current
---

# Phred quality-score encoding (Phred+33 / Phred+64)

A **Phred quality score** attaches a per-base error probability to a sequencing read. It is the
numeric payload of the FASTQ quality line and the input to quality trimming, filtering, and QC
statistics — so its ASCII encoding is a cross-cutting concern rather than a detail of any one
format. This page is the shared reference for that encoding; it is the quality-decoding step consumed by
the FASTQ parser concept [[fastq-parsing]], first traced in the FASTQ parsing evidence
[[parse-fastq-001-evidence]], and reused by [[quality-trimming-running-sum]]. The dedicated **Phred score handling** unit
[[quality-phred-001-evidence]] (parse / encode / cross-variant convert) is anchored on this page's
formulae and pins them to the primary FASTQ specification, **Cock et al. (2010)**.
[[test-unit-registry]] tracks the owning units.

## The score: a log-scaled error probability

- **Definition:** `Q = −10 · log₁₀(p)`, where `p` is the probability the base call is wrong.
- **Inverse:** `p = 10^(−Q/10)`.
- Every **10 Q points = a 10× drop in error rate**:

| Q | Error probability `p` | Accuracy |
|---|-----------------------|----------|
| 0 | 1.0 (100%) | — |
| 10 | 0.1 (1 in 10) | 90% |
| 20 | 0.01 (1 in 100) | 99% |
| 30 | 0.001 (1 in 1,000) | 99.9% |
| 40 | 0.0001 (1 in 10,000) | 99.99% |

**Q20** and **Q30** are the conventional QC thresholds — the "% of bases ≥ Q20 / ≥ Q30"
reported by FASTQ statistics. Those run-quality summaries (Q20/Q30 fractions, mean, median,
min/max, std dev) are a **separate surface** that *consumes* these decoded scores; they live on
their own concept page [[fastq-quality-statistics]] (unit [[quality-stats-001-evidence]]). The
inverse conversion caps at **Q93** for `p ≤ 0` (the maximum value representable in a Sanger
Phred+33 string).

## The two live encodings (ASCII offset)

A quality character encodes `Q` as `chr(Q + offset)`; decoding is `Q = ord(char) − offset`.
Two offsets are in use, and the difference is the single most common FASTQ parsing hazard:

| Encoding | Offset | ASCII range | Q range | Era |
|----------|--------|-------------|---------|-----|
| **Sanger / Illumina 1.8+** | **33** | 33–126 (`!`–`~`) | 0–93 | Current standard |
| **Illumina 1.3–1.7** | **64** | 64–126 (`@`–`~`) | 0–62 | Legacy |
| Solexa (obsolete) | 64 | 59–126 | −5 to 62 | Uses a different odds-ratio score, not straight Phred |

Boundary characters worth memorizing: in **Phred+33**, `!`=Q0, `I`=Q40, `~`=Q93 (max); in
**Phred+64**, `@`=Q0, `h`=Q40, `~`=Q62 (max).

## Auto-detecting the offset

The encoding is not written in the file, so it is inferred from the character range observed
across quality strings:

- A character **below `@` (ASCII 64)** can only be Phred+33 → **Phred+33**.
- A character **above `I` (ASCII 73)** in the low band that would imply Q > 40 under +33 without
  any low chars indicates **Phred+64**.
- The ambiguous window `@`–`I` (both encodings overlap there) is undecidable from range alone
  and **defaults to Phred+33** (the current standard).

Because detection reads a bounded alphabet, it is **deterministic** for a given input; Seqeron
runs it **per record**.

## Converting between the two offsets

Because the **Phred score is invariant** across the Sanger (offset 33) and Illumina 1.3+
(offset 64) variants, converting a quality string between them is a **pure re-offset** — shift
every byte by ±31, with no numeric rescaling (Cock et al. 2010). Two representability rules bound
the conversion:

- **Phred+64 → Phred+33 is always safe:** Phred+64 holds `Q ∈ [0, 62] ⊆ [0, 93]`, so every valid
  Phred+64 score fits Phred+33. Worked: `@h~` → `!I_` (Q 0/40/62).
- **Phred+33 → Phred+64 can overflow:** a Phred+33 score in `(62, 93]` exceeds the Phred+64
  maximum and cannot be encoded (Seqeron raises `ArgumentOutOfRangeException`; e.g. `~`=Q93).
  Worked in-range: `!I` → `@h` (Q 0/40).

A byte **below the variant's offset** decodes to a negative Phred score — invalid for either
variant (`Q ≥ 0`) — and is treated as malformed input (`ArgumentOutOfRangeException`). Contrast
this with the **Solexa** obsolete variant, whose odds-ratio score `Q = −10·log₁₀(p/(1−p))`
requires a *lossy numeric* conversion (not a re-offset) and is out of scope for the Phred path.

## Why the distinction bites

Decoding a Phred+64 file as Phred+33 (or vice versa) shifts every score by 31 and silently
corrupts all downstream Q20/Q30 stats, trimming cut-points, and quality filters — with no parse
error, because both offsets produce printable ASCII. Auto-detection exists precisely to stop that
class of silent error; when it cannot decide, defaulting to the modern Phred+33 is the
lowest-risk choice for contemporary data.

## Implementation (QualityScoreAnalyzer)

The **Phred score handling** unit [[quality-phred-001-evidence]] (unit `QUALITY-PHRED-001`) exposes
three canonical methods on `QualityScoreAnalyzer` in `Seqeron.Genomics.IO`
(`src/Seqeron/Algorithms/Seqeron.Genomics.IO/QualityScoreAnalyzer.cs`), each a single **O(n) linear
pass** over the n characters/scores (O(n) space):

- `ParseQualityString(qualityString, encoding) → int[]` — decodes `Q = ord(char) − offset`, validating
  each `Q` within the encoding's range; `Auto` is resolved through the existing `DetectEncoding`
  heuristic before decoding.
- `ToQualityString(scores, encoding) → string` — encodes `char = chr(Q + offset)` with the same
  range validation; `Auto` here is treated as the modern **Phred+33** default (there is nothing to
  detect on the encode path).
- `ConvertEncoding(qualityString, fromEncoding, toEncoding) → string` — parse under the source offset,
  re-encode under the target offset (the ±31 re-offset).

The offsets (33/64) and valid ranges (`[0,93]` / `[0,62]`) are **named constants** citing Cock et al.
(2010). Parse/encode form mutual inverses on the valid range — the round-trip invariant
`ToQualityString(ParseQualityString(s, e), e) == s` (INV-03). Contract edges: **null input** raises
`ArgumentNullException`; a character or score **outside the encoding's valid range** (including the
Phred+33→Phred+64 overflow for `Q > 62`) raises `ArgumentOutOfRangeException`; empty input yields empty
output. This is a specification-driven exact re-offset, not a search/matching unit, so the repository
suffix tree does not apply. Pre-existing `QualityStringToPhred` / `PhredToQualityString` helpers remain
for other callers; the canonical methods above add the explicit range validation those helpers lack.
