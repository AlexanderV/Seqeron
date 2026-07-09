---
type: source
title: "Evidence: KMER-GENERATE-001 (K-mer generation / n^k universe enumeration)"
tags: [validation, analysis]
doc_path: docs/Evidence/KMER-GENERATE-001-Evidence.md
sources:
  - docs/Evidence/KMER-GENERATE-001-Evidence.md
source_commit: 07cd59444e1d3f85403bf468c40c6de97216c385
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: KMER-GENERATE-001

The validation-evidence artifact for test unit **KMER-GENERATE-001** — **k-mer generation**
(`KmerAnalyzer.GenerateAllKmers`): exhaustive enumeration of *every* possible k-mer of length `k`
over an alphabet (the complete `n^k` k-mer universe), independent of any particular sequence. This
is the fourth **K-mer** family Evidence file (after [[kmer-async-001-evidence|KMER-ASYNC-001]],
[[kmer-both-001-evidence|KMER-BOTH-001]], and [[kmer-dist-001-evidence|KMER-DIST-001]]) and one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
the operation itself is synthesized in [[k-mer-generation]]. See [[test-unit-registry]] for how units
are tracked.

It validates that generation is **distinct from k-mer counting** — it produces the full Cartesian
product `Σ^k` of all possible words (sequence-independent), not the multiset of substrings observed
in a given sequence.

## What this file records

- **Online sources:**
  - **Wikipedia — K-mer** (authority rank 4) — the foundational definition ("substrings of length k
    contained within a biological sequence") and the universe-size formula: "there exist n^k total
    possible k-mers, where n is number of possible monomers (e.g. four in the case of DNA)" ⇒ **4^k**
    for DNA {A,C,G,T}. The `AGAT` worked example fixes a k-mer as a contiguous length-k string.
  - **BioInfoLogics — k-mer counting part I (Clavijo 2018)** (rank 3) — independently confirms and
    gives the constructive reason: "Since each of the k nucleotides in a k-mer can take any of the A,
    C, G or T values, the possible combinations of k positions are computed as 4^k" (i.e. the k-fold
    Cartesian product of the alphabet).
  - **Python stdlib — `itertools.product`** (rank 3) — the canonical reference-implementation
    behaviour the generation algorithm realises: `product(A, repeat=k)` = the k-fold Cartesian
    product = all length-k tuples; **odometer / lexicographic emission** ("if the input's iterables
    are sorted, the product tuples are emitted in sorted order", rightmost position advances fastest),
    with the `product(range(2), repeat=3) → 000 001 010 011 100 101 110 111` enumeration structurally
    identical to k-mer generation over sorted {A,C,G,T}.

- **Documented corner cases / failure modes:** k-mer length must be positive (`k ≤ 0` has no meaning
  as a length); single-letter alphabet → exactly one k-mer (`1^k = 1`, the homopolymer); the
  lexicographic-emission guarantee holds **only if the input alphabet is sorted** (an unsorted
  alphabet still yields all `n^k` strings, but in the alphabet's own positional order).

- **Datasets (documented oracles):**
  - DNA `k=1` → `{A, C, G, T}` (4^1 = 4).
  - DNA `k=2` → all 16 (4^2) 2-mers in lexicographic order `AA, AC, AG, AT, CA, …, TT` (first `AA`,
    last `TT`).
  - DNA `k=3` → 64 (4^3) distinct k-mers, first `AAA`, second `AAC`, last `TTT` (`|set| = 64`).
  - Protein alphabet `ACDEFGHIKLMNPQRSTVWY` (20 letters), `k=2` → 400 (20^2).
  - Single-letter alphabet `A`, `k=4` → exactly `{AAAA}` (1^4 = 1).

- **Test-coverage recommendations:** MUST — `GenerateAllKmers(1)/(2)/(3)` exact sets and ordering;
  count = `4^k` for `k = 1..6` on the default DNA alphabet; custom-alphabet count = `(len)^k` (protein
  20^2 = 400). SHOULD — single-letter alphabet `1^k = 1`; no duplicates (output = the Cartesian-product
  set); `k ≤ 0` → `ArgumentOutOfRangeException`; empty/null alphabet → `ArgumentException`.

## Deviations and assumptions

**One ASSUMPTION** (documented property, not an invented value): the **default alphabet is `"ACGT"`
in upper case**, matching the sibling `KmerAnalyzer` methods; the alphabet is used verbatim, so the
lexicographic-ordering guarantee holds only when the supplied alphabet is itself sorted (`{A,C,G,T}`
already is, so default DNA output is lexicographic). No source contradictions.
</content>
</invoke>
