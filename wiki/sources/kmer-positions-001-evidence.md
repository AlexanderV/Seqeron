---
type: source
title: "Evidence: KMER-POSITIONS-001 (K-mer positions / occurrence index)"
tags: [validation, analysis]
doc_path: docs/Evidence/KMER-POSITIONS-001-Evidence.md
sources:
  - docs/Evidence/KMER-POSITIONS-001-Evidence.md
source_commit: f43e0daf5dd06eba936e6e14939946e4cd980b67
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: KMER-POSITIONS-001

The validation-evidence artifact for test unit **KMER-POSITIONS-001** — **k-mer positions**
(`KmerAnalyzer.FindKmerPositions`): the ascending, 0-based list of every starting index at which a
given k-mer (a fixed pattern) occurs in a sequence — a *position / occurrence index*, distinct from
counting (how many times). This is the fifth **K-mer** family Evidence file (after
[[kmer-async-001-evidence|KMER-ASYNC-001]], [[kmer-both-001-evidence|KMER-BOTH-001]],
[[kmer-dist-001-evidence|KMER-DIST-001]], [[kmer-generate-001-evidence|KMER-GENERATE-001]]) and one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
the operation itself is synthesized in [[k-mer-positions]]. See [[test-unit-registry]] for how units
are tracked.

It validates the classical **Pattern Matching Problem** — "find all occurrences of a pattern in a
string", reporting all 0-based starting positions including overlapping occurrences.

## What this file records

- **Online sources:**
  - **Rosalind — "Find All Occurrences of a Pattern in a String" (Problem BA1D)** (authority rank 4)
    — the binding worked example and conventions: output "all starting positions in *Genome* where
    *Pattern* appears as a substring", **0-based**; overlapping occurrences all reported; verbatim
    example Pattern `ATAT` / Genome `GATATATGCATATACTT` → `1 3 9`.
  - **Wikipedia — "k-mer"** (authority rank 4) — the formal definition (substring of length k), the
    **L − k + 1** candidate-position count, and the overlapping-consecutive-k-mer property
    (`AGAT` → 2-mers `AG`@0, `GA`@1, `AT`@2).
  - **Compeau & Pevzner — *Bioinformatics Algorithms: An Active Learning Approach*** (authority rank
    1, textbook) — the Pattern Matching Problem statement. The textbook prose narrates 1-based
    positions, but the canonical machine-checked Rosalind BA1D exercise specifies **0-based** output,
    recorded as the binding convention for this unit.

- **Documented corner cases / failure modes:** self-overlap → every overlapping start reported
  (`AA` in `AAAA` → `0 1 2`); pattern longer than text → 0 positions (`L − k + 1 ≤ 0`) → empty;
  pattern absent → empty; pattern equals whole sequence → exactly `[0]`.

- **Datasets (documented oracles):**
  - Rosalind BA1D: `ATAT` in `GATATATGCATATACTT` → `[1, 3, 9]` (in order).
  - Wikipedia `AGAT`, k=2 → 2-mers in order `AG`@0, `GA`@1, `AT`@2 (each start ascending).
  - Self-overlap derivation: `AA` in `AAAA` → `[0, 1, 2]` (`L − k + 1 = 3`).

- **Test-coverage recommendations:** MUST — Rosalind BA1D `[1,3,9]`; overlapping `AA`/`AAAA`→`[0,1,2]`;
  `ana`/`banana`-style ascending starts; pattern absent → empty. SHOULD — pattern longer than text →
  empty; pattern equals whole sequence → `[0]`; case-insensitive match (`atat` ≡ `ATAT`). COULD —
  null/empty sequence or k-mer → empty.

## Deviations and assumptions

**Three ASSUMPTIONS**, all API-shape / repository-interoperability and non-correctness-affecting on
the all-uppercase evidence examples:

- **0-based indexing** — adopted per the canonical machine-checked Rosalind BA1D exercise (over the
  textbook's 1-based prose), consistent with C# string indexing and the repository `SuffixTree`
  (`FindAllOccurrences("ana")` → `[1,3]` for `banana`).
- **Case-insensitive matching** — no source mandates case-folding; both arguments are upper-cased to
  match the sibling `KmerAnalyzer.CountKmers` methods (so `atat` matches `ATAT`).
- **Null/empty input returns empty** — no source defines null/empty behaviour; resolved to an empty
  result (no exception), matching the sibling K-mer methods.

No source contradictions — Rosalind BA1D and Wikipedia agree on the occurrence-set definition,
overlapping rule, and the `L − k + 1` bound; the only 1-based/0-based tension is resolved to the
machine-checked BA1D convention.
