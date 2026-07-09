---
type: source
title: "Evidence: KMER-DIST-001 (K-mer Euclidean distance, alignment-free)"
tags: [validation, analysis]
doc_path: docs/Evidence/KMER-DIST-001-Evidence.md
sources:
  - docs/Evidence/KMER-DIST-001-Evidence.md
source_commit: e59aee7693f6331696abe9d68f59844f570fcf17
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: KMER-DIST-001

The validation-evidence artifact for test unit **KMER-DIST-001** — **K-mer Euclidean distance**
(`KmerAnalyzer.KmerDistance`): the alignment-free L2 distance between two sequences' normalized
k-mer **frequency** vectors. This is the third **K-mer** family Evidence file (after
[[kmer-async-001-evidence|KMER-ASYNC-001]] and [[kmer-both-001-evidence|KMER-BOTH-001]]) and one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the operation itself is synthesized in [[k-mer-euclidean-distance]]. See
[[test-unit-registry]] for how units are tracked.

It validates the **frequency-vector** variant — counts normalized by the number of k-mer windows
`L − k + 1` (Lau et al. 2022), then Euclidean distance over the union of observed k-mers — and is
careful to distinguish it from the count-based (un-normalized) form used only for an intermediate
derivation.

## What this file records

- **Online sources (all authority rank 1):**
  - **Zielezinski, Vinga, Almeida & Karlowski 2017** — *Alignment-free sequence comparison*
    (Genome Biology 18:186). The word-vector representation and the Figure 1 worked example
    (x=`ATGTGTG`, y=`CATGTG`, k=3; union W₃={ATG,CAT,GTG,TGT}; counts c_x=(1,0,2,2),
    c_y=(1,1,1,1)); "this difference is very commonly computed by the Euclidean distance";
    identical sequences → distance 0; absent word → 0 component.
  - **Lau et al. 2022** — *Interpreting alignment-free sequence comparison* (NAR Genomics and
    Bioinformatics). The **frequency definition** (verbatim): counts divided by the total number
    of k-mers = `L − k + 1`; lists Euclidean/Manhattan/Canberra/Chebyshev as usable metrics; Table 1
    marks Euclidean as suitable for variable-length sequences.
  - **Vinga & Almeida 2003** — *Alignment-free sequence comparison—a review* (Bioinformatics
    19(4):513). The 4^k-dimensional word-frequency vector mapping and the family of scores
    (Euclidean, Pearson, Kullback–Leibler, cosine). Abstract/metadata only (full text paywalled).
  - **Boden et al. 2014** — *Fast alignment-free sequence comparison using spaced-word frequencies*
    (Bioinformatics 30(14):1991). Confirms the standard Euclidean distance is taken over
    **relative (normalized) frequency** vectors, not raw counts.

- **Documented corner cases / failure modes:** identical sequences → 0; the comparison vector spans
  the **union** of words in either sequence (absent word contributes a 0 component); frequency
  normalization requires `L ≥ k` (at least one window), else an empty vector.

- **Datasets (documented oracles):**
  - **Zielezinski 2017 Fig. 1** — x=`ATGTGTG`, y=`CATGTG`, k=3. Count-based Euclidean (intermediate
    derivation): differences (0,−1,1,1) ⇒ √3 ≈ 1.7320508076. **Frequency-based** (what
    `KmerDistance` computes): f_x=(0.2,0,0.4,0.4), f_y=(0.25,0.25,0.25,0.25), squared diffs sum
    0.11 ⇒ **√0.11 ≈ 0.33166247903553997**.
  - **Single-substitution k=1** (direct derivation) — `AAAA` vs `AAAT`: f₁=(A 1.0, T 0),
    f₂=(A 0.75, T 0.25) ⇒ √0.125 ≈ **0.3535533906**.

- **Test-coverage recommendations:** MUST — Fig. 1 example (√0.11), identical → exactly 0,
  single-substitution k=1 (√0.125), symmetry d(x,y)=d(y,x). SHOULD — non-overlapping single-k-mer
  sets → √2, k>min(length) → empty vector so distance = norm of the other's frequency vector.
  COULD — k≤0 throws `ArgumentOutOfRangeException`, case-insensitivity (lower==upper).

## Deviations and assumptions

**Two ASSUMPTIONS**, both benign / non-boundary-affecting on canonical input:

- **Case folding (ASM-01)** — inputs are upper-cased before counting (via `CountKmers`), so
  mixed-case inputs yield identical k-mers. No source specifies case handling; benign for canonical
  upper-case input.
- **Empty / too-short input (ASM-02)** — a sequence with no k-mer windows (`L < k`, empty, or null)
  yields an empty frequency vector treated as the zero vector; the distance then equals the L2 norm
  of the other sequence's frequency vector, and 0 when both are empty. The natural extension where
  sources define frequencies only for `L ≥ k`.

No source contradictions — Zielezinski's word-vector model, Lau's `L − k + 1` normalization, and
Boden's "Euclidean over relative-frequency vectors" are mutually consistent; the only implementation
choices (case folding, empty→zero-vector) are documented assumptions.
</content>
</invoke>
