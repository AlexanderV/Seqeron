---
type: source
title: "Evidence: GENOMIC-SIMILARITY-001 (Sequence similarity, k-mer Jaccard index)"
tags: [validation, analysis]
doc_path: docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md
sources:
  - docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md
source_commit: f2b9ce29b93a0977bf8cc2d4d003a59711a6534b
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: GENOMIC-SIMILARITY-001

The validation-evidence artifact for test unit **GENOMIC-SIMILARITY-001** — **sequence similarity by
the k-mer Jaccard index** (`GenomicAnalyzer.CalculateSimilarity`). One instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm itself is
synthesized in the concept [[kmer-jaccard-similarity]]. See [[test-unit-registry]] for how units are
tracked.

## What this file records

- **Online sources:**
  - **Jaccard index** (Wikipedia citing the primary **Jaccard 1901**; authority rank 4/primary rank 1) — the verbatim formula `J(A,B) = |A∩B|/|A∪B| = |A∩B|/(|A|+|B|−|A∩B|)`, range `0 ≤ J ≤ 1`, the scope note that it "measures similarity between finite **non-empty** sample sets" and is "not well-defined when μ(A ∪ B) = 0", and Jaccard distance `d_J = 1 − J`.
  - **Ondov et al. 2016, *Mash*** (peer-reviewed, *Genome Biology* 17:132; authority rank 1) — the k-mer-set application: the Jaccard index is "the fraction of shared hashes … out of all distinct hashes in A and B", `J(A,B)=|A∩B|/|A∪B|` over k-mer sets, and the MinHash *sketch* estimate `j(A_s,B_s)=|A_s∩B_s|/s`. Mash uses k=21 for whole genomes.
  - **Mash documentation — Distance Estimation** (reference-implementation docs, rank 3) — confirms the sketch Jaccard `j(A_s,B_s)=|A_s∩B_s|/s` as shared-k-mers / distinct-k-mers.

- **Documented corner cases / failure modes:** empty union (both sets empty) = undefined, no authoritative value (impl returns 0); identical sets → J=1; disjoint sets → J=0; distinct-k-mer (set, not multiset) decomposition so repeated k-mers count once; choice of k changes *resolution* only, not the formula.

- **Datasets (hand-derived oracles, k=3):** `ACGTACGT`/`ACGTACGA` → J×100 = **80.0**; `ACGT`/`ACGA` → **100/3 ≈ 33.33**; `ACGTACGT`/`ACGTACGT` → **100.0**; `AAAAA`/`CCCCC` → **0.0**. Distinct-set semantics: `AAAAAA`/`AAAA` at k=3 → **100** (both sets `{AAA}`).

- **Test-coverage recommendations:** MUST — partial-overlap exact fraction (80.0); identical → 100; disjoint → 0; non-integer fraction (100/3); distinct-k-mer set semantics (repeats once); symmetry `J(A,B)=J(B,A)`. SHOULD — empty-union convention → 0 (both empty / both shorter than k); input validation (null → `ArgumentNullException`, `kmerSize < 1` → `ArgumentOutOfRangeException`). COULD — range invariant `0 ≤ result ≤ 100`.

## Deviations and assumptions

**Three ASSUMPTIONS**, all source-backed, none an open correctness gap:

- **Empty-union return value = 0.0** — Jaccard is undefined when the union is empty (Jaccard 1901; "not well-defined when μ(A ∪ B) = 0"). The implementation returns `0.0` ("no shared content"), an implementation convention (either 0 or 1 appears in practice), documented and tested as the contract, not asserted as a literature value.
- **Percentage scaling (×100)** — the formal index is in `[0,1]`; ×100 is a presentation convention that leaves relative ordering unchanged.
- **Default k = 5** — no source mandates a default k for short-DNA similarity (Mash uses k=21 for genomes); k=5 is a project default setting resolution only. All evidence-based tests pass k explicitly.

No source contradictions — Jaccard's non-empty-set definition and Mash's k-mer-set application are mutually consistent; the empty-union value, ×100 scaling, and default k are the only (documented) implementation choices. The algorithm doc `docs/algorithms/Analysis/Sequence_Similarity.md` additionally notes the **suffix tree was evaluated but not used** (a `HashSet<string>` gives optimal O(n+m) set-overlap; a suffix tree solves occurrence-search, not set resemblance).
