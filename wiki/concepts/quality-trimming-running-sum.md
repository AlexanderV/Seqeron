---
type: concept
title: "Quality trimming (BWA/cutadapt running-sum)"
tags: [assembly, algorithm]
sources:
  - docs/Evidence/ASSEMBLY-TRIM-001-Evidence.md
  - docs/algorithms/Assembly/Quality_Trimming.md
  - docs/Validation/reports/ASSEMBLY-TRIM-001.md
source_commit: 9e3afa723f48771feef69632da397e2992f74114
created: 2026-07-09
updated: 2026-07-10
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
test unit **ASSEMBLY-TRIM-001**; the pre-implementation record is [[assembly-trim-001-evidence]], the
independent two-stage verdict is [[assembly-trim-001-report]] (Stage B **FAIL → FIXED** — see below),
and [[test-unit-registry]] tracks the unit. See [[algorithm-validation-evidence]] for the artifact
pattern.

Unlike the other assembly-family units, trimming reconstructs no sequence — it operates on a single
read's Phred quality string. It sits upstream of [[kmer-spectrum-error-correction|error correction]]
and the [[de-bruijn-graph-assembly|DBG]] / [[overlap-layout-consensus-assembly|OLC]] build steps.

## The running-sum core

cutadapt states its algorithm "is the same as the one used by BWA". Given a per-base Phred quality
array and a cutoff threshold, the 3'-end trim is:

1. **Subtract the threshold** from every quality: `q - cutoff`.
2. **Accumulate from the end, tracking the running max** — and **break as soon as the sum goes
   negative** (`if s < 0: break`), i.e. stop at the argmin/argmax before the sum could recover.
3. **Cut at the tracked boundary** — keep everything 5' of it.

The intent (cutadapt, verbatim): "remove all bases starting from the end of the read whose quality is
smaller than the given threshold ... refined a bit by allowing some good-quality bases among the
bad-quality ones." The `s < 0` early break **is** that refinement: a lone high-quality base inside a
bad tail survives only if the cumulative sum has not already gone negative before it (validated by
[[assembly-trim-001-report]] — omitting the break was the shipped defect).

**Both ends** (`quality_trim_index`, the authoritative reference): the 5' and 3' passes run
**independently over the full read `[0, n)`** — the 5' pass is *not* chained onto the 3'-surviving
window — yielding `(start, stop)`. Then the **drop rule**: `if start >= stop: (start, stop) = (0, 0)`,
i.e. when the two passes cross, the whole read is discarded.

## BWA's equivalent argmax form

BWA's `bwa_trim_read` (bwaseqio.c, Heng Li) accumulates `s += trim_qual - (q - 33)` walking from the
3' end and tracks the **argmax** position `max_l`. This is algebraically the same optimum: the argmax
of accumulated `(threshold − q)` from the end equals cutadapt's argmin of partial sums of
`(q − threshold)`. Two details:

- **Early break** when `s < 0` (the running sum can never recover to a new max after that). This is
  **shared** with cutadapt's `quality_trim_index` reference — it is *not* a BWA-only quirk, and its
  absence was the ASSEMBLY-TRIM-001 defect ([[assembly-trim-001-report]]).
- **Hard floor** `BWA_MIN_RDLEN = 35` (BWA-specific): BWA never trims a read below 35 bp, independent
  of the running-sum optimum.

## Phred+33 decoding

Qualities are decoded from ASCII with the **Sanger/Phred+33** offset: `q = ASCII − 33` (Cock et al.
2010: ASCII 33–126 encode Phred 0–93). Phred `Q = −10·log₁₀(P)`, but the trimming operates on the
integer Phred values directly. See [[phred-quality-encoding]] for the full encoding scheme
(Phred+33 vs Phred+64 offsets and offset auto-detection).

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
optimum. Similarly, the **both-end pass order** (5' vs 3' first) is not numerically significant: after
the fix ([[assembly-trim-001-report]]) both passes run **independently over the full read**, so they
act on disjoint ends and the order cannot change the `(start, stop)` result.
