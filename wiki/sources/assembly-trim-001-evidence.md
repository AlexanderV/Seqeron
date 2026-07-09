---
type: source
title: "Evidence: ASSEMBLY-TRIM-001 (Quality trimming — BWA/cutadapt running-sum)"
tags: [validation, assembly]
doc_path: docs/Evidence/ASSEMBLY-TRIM-001-Evidence.md
sources:
  - docs/Evidence/ASSEMBLY-TRIM-001-Evidence.md
source_commit: 6b24d624caf0d1b12aba8c6569fd25efe0c496ee
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ASSEMBLY-TRIM-001

The validation-evidence artifact for test unit **ASSEMBLY-TRIM-001** — quality trimming: the
BWA/cutadapt running-sum algorithm that removes low-quality bases from read ends. One instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the running-sum
core, the BWA argmax equivalence, Phred+33 decoding, the published oracle and the two assumptions are
summarized in [[quality-trimming-running-sum]], the anchor for the assembly TRIM family. See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** (all accessed 2026-06-13):
  - **cutadapt — "Algorithm details (Quality trimming)"** (rank 3; readthedocs, WebSearch then
    WebFetch, cross-checked v1.18) — states the algorithm "is the same as the one used by BWA" and
    gives the verbatim procedure: subtract the cutoff from all qualities, compute partial sums from
    each index to the end, cut at the index where the sum is minimal; "repeat this for the other end"
    if both ends are trimmed; plus the intent (remove low-quality end bases, allowing some good bases
    among the bad).
  - **BWA source — `bwa_trim_read` in `bwaseqio.c`** (rank 3; canonical reference implementation,
    Heng Li) — the verbatim C accumulating `s += trim_qual - (p->qual[l] - 33)` from the 3' end,
    tracking the argmax `max_l`, with the `s < 0` early break and the `BWA_MIN_RDLEN` floor; Phred+33
    decode `qual - 33`; the `trim_qual < 1 || qual == 0` no-trim guard.
  - **BWA source — `BWA_MIN_RDLEN` in `bwtaln.h`** (rank 3) — `#define BWA_MIN_RDLEN 35`, BWA's hard
    floor on trimmed read length, separate from the running-sum optimum.
  - **Cock et al. (2010) — "The Sanger FASTQ file format ..."** (rank 1; *Nucleic Acids Research*
    38(6):1767–1771) — Sanger FASTQ uses ASCII 33–126 for Phred 0–93 (offset 33); Phred
    `Q = −10·log₁₀(P)` (background; trimming uses integer Phred directly).
- **Dataset (published oracle):**
  - **cutadapt worked example** — qualities `42,40,26,27,8,7,11,4,2,3` (Phred+33 `KI;<)(,%#$`),
    threshold 10, partial sums from end `(70),(38),8,-8,-25,-23,-20,-21,-15,-7`, minimum `-25` at
    index 4 → **first 4 bases kept**. The evidence file also gives the full Phred+33 ASCII derivation
    (42→`K`, 40→`I`, 26→`;`, 27→`<`, 8→`)`, 7→`(`, 11→`,`, 4→`%`, 2→`#`, 3→`$`).
- **Corner cases / failure modes** — threshold `< 1` disables trimming (BWA guard; cutoff 0 → all
  sums non-negative → nothing trimmed); all-high-quality → unchanged; all-low-quality → read fully
  removed; a good base among bad ones retained only if the cumulative sum has not hit a new minimum
  (cutadapt's "refinement").
- **Recommended coverage** — MUST: the cutadapt worked example (→ first 4 bases); Phred+33 decode
  `q = ASCII − 33`; all-high-quality unchanged; all-low-quality dropped (length 0 < minLength).
  SHOULD: 5'-end trimming; `minLength` filter drops short survivors / keeps ≥ minLength; threshold ≤ 0
  disables trimming. COULD: good-base-among-bad retention; empty/null/empty-quality inputs.

## Assumptions (from the artifact)

Two assumption records, both outside the source-backed running-sum optimum:

1. **Both-end pass order.** cutadapt trims both ends "in turn"; the order (5' vs 3' first) is not
   numerically significant because the two passes operate on disjoint ends of the surviving window.
   The repository does 3' then 5' on the remaining window — consistent with cutadapt's "repeat for the
   other end" and unchanged for the published (pure-3') example.
2. **`minLength` filter semantics.** BWA/cutadapt's running-sum core defines no minimum-length drop;
   the standard post-trim min-length filter (cutadapt `--minimum-length`) drops reads whose trimmed
   length `< minLength`. A documented downstream filter, not part of the running-sum optimum.

No contradictions among the sources — cutadapt explicitly identifies its algorithm with BWA's, the
BWA argmax of accumulated `(threshold − q)` is algebraically the argmin of cutadapt's partial sums of
`(q − threshold)`, and Cock et al. supply the Phred+33 encoding both rely on.
