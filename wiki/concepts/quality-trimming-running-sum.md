---
type: concept
title: "Quality trimming (BWA/cutadapt running-sum)"
tags: [assembly, algorithm]
sources:
  - docs/Evidence/ASSEMBLY-TRIM-001-Evidence.md
source_commit: 6b24d624caf0d1b12aba8c6569fd25efe0c496ee
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: assembly-trim-001-evidence
      evidence: "Test Unit ID: ASSEMBLY-TRIM-001 ... Quality Trimming (BWA / cutadapt running-sum)"
      confidence: high
      status: current
---

# Quality trimming (BWA/cutadapt running-sum)

Removing low-quality bases from the ends of a sequencing read using the **running-sum** algorithm
shared by **BWA** (`bwa_trim_read`) and **cutadapt**. This is a read-preprocessing/QC step that runs
*before* assembly or alignment; it is the anchor for the assembly **TRIM** family. Validated under
test unit **ASSEMBLY-TRIM-001**; the validation record is [[assembly-trim-001-evidence]], and
[[test-unit-registry]] tracks the unit. See [[algorithm-validation-evidence]] for the artifact
pattern.

Unlike the other assembly-family units, trimming reconstructs no sequence — it operates on a single
read's Phred quality string. It sits upstream of [[kmer-spectrum-error-correction|error correction]]
and the [[de-bruijn-graph-assembly|DBG]] / [[overlap-layout-consensus-assembly|OLC]] build steps.

## The running-sum core

cutadapt states its algorithm "is the same as the one used by BWA". Given a per-base Phred quality
array and a cutoff threshold, the 3'-end trim is:

1. **Subtract the threshold** from every quality: `q - cutoff`.
2. **Compute partial sums from each index to the 3' end** of the read.
3. **Cut at the index where that partial sum is minimal** — keep everything 5' of it.

The intent (cutadapt, verbatim): "remove all bases starting from the end of the read whose quality is
smaller than the given threshold ... refined a bit by allowing some good-quality bases among the
bad-quality ones." A lone high-quality base inside a bad tail survives only if the cumulative sum
does not hit a new minimum before it.

**Both ends:** if requested, "repeat this for the other end" — a 5' pass on the surviving window.

## BWA's equivalent argmax form

BWA's `bwa_trim_read` (bwaseqio.c, Heng Li) accumulates `s += trim_qual - (q - 33)` walking from the
3' end and tracks the **argmax** position `max_l`. This is algebraically the same optimum: the argmax
of accumulated `(threshold − q)` from the end equals cutadapt's argmin of partial sums of
`(q − threshold)`. Two BWA-specific implementation details:

- **Early break** when `s < 0` (the running sum can never recover to a new max after that).
- **Hard floor** `BWA_MIN_RDLEN = 35`: BWA never trims a read below 35 bp, independent of the
  running-sum optimum.

## Phred+33 decoding

Qualities are decoded from ASCII with the **Sanger/Phred+33** offset: `q = ASCII − 33` (Cock et al.
2010: ASCII 33–126 encode Phred 0–93). Phred `Q = −10·log₁₀(P)`, but the trimming operates on the
integer Phred values directly.

## Worked example (cutadapt docs, the test oracle)

Qualities `42, 40, 26, 27, 8, 7, 11, 4, 2, 3` (Phred+33 string `KI;<)(,%#$`), threshold `10`. After
subtracting: `32, 30, 16, 17, -2, -3, 1, -6, -8, -7`. Partial sums from the end:
`(70), (38), 8, -8, -25, -23, -20, -21, -15, -7`. The minimum `-25` is at index 4 → the read is
trimmed to the **first four bases** (qualities 42, 40, 26, 27).

## Corner cases / failure modes

- **Threshold disables trimming:** BWA returns 0 (no trim) when `trim_qual < 1`; a cutoff of 0 leaves
  all partial sums non-negative, minimum at the last index → nothing removed.
- **All-high-quality read:** every `q − threshold ≥ 0` → no bases removed.
- **All-low-quality read:** every `q − threshold < 0` → the whole 3' (and, on the 5' pass, the whole
  5') span is removed; the read may be fully consumed.
- **Good base among bad ones:** retained iff the cumulative sum has not already reached a new minimum.

## Downstream min-length filter (assumption)

The running-sum core defines **no** minimum-length drop. The `minLength` parameter in the test-unit
signature is the standard post-trim length filter (cutadapt `--minimum-length`): reads whose trimmed
length `< minLength` are dropped. This is a documented downstream filter, not part of the running-sum
optimum. Similarly, the **both-end pass order** (this unit does 3' then 5' on the surviving window) is
not numerically significant because the two passes act on disjoint ends of the surviving window.
