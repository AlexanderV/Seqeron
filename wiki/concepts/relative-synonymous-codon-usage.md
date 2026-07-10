---
type: concept
title: "Relative Synonymous Codon Usage (RSCU)"
tags: [annotation, algorithm]
sources:
  - docs/Evidence/ANNOT-CODONUSAGE-001-Evidence.md
  - docs/Evidence/CODON-RSCU-001-Evidence.md
  - docs/algorithms/Codon/Relative_Synonymous_Codon_Usage.md
  - docs/algorithms/Annotation/Relative_Synonymous_Codon_Usage.md
  - docs/Validation/reports/ANNOT-CODONUSAGE-001.md
source_commit: 987ea6c1cf04c61c6257f0034ea4d51e00e0fffc
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: annot-codonusage-001-evidence
      evidence: "Test Unit ID: ANNOT-CODONUSAGE-001 ... Algorithm: Relative Synonymous Codon Usage (RSCU)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: codon-rscu-001-evidence
      evidence: "Test Unit ID: CODON-RSCU-001 ... Algorithm: Relative Synonymous Codon Usage (RSCU) and codon counting"
      confidence: high
      status: current
---

# Relative Synonymous Codon Usage (RSCU)

A per-codon **codon-usage-bias** measure: within each amino acid's set of synonymous codons,
how over- or under-represented one codon is relative to uniform (unbiased) usage. Introduced by
**Sharp & Li (1986)**, it is the canonical normalization for comparing codon preference across
genes, genomes, and organisms independent of amino-acid composition. Validated by **two**
evidence units: [[annot-codonusage-001-evidence|ANNOT-CODONUSAGE-001]] (against the LIRMM
methods page, PMC2528880, and the CodonU reference implementation) and
[[codon-rscu-001-evidence|CODON-RSCU-001]] (same measure plus the supporting `CountCodons`
operation, cross-checked against LIRMM, GenomicSig, seqinr `uco`, and cubar `est_rscu`). See
[[test-unit-registry]] for how the units are tracked. The two-stage validation write-up
([[annot-codonusage-001-report|ANNOT-CODONUSAGE-001 report]]) verdict is **✅ CLEAN** — no code
defect; the formula was confirmed verbatim and only three test-coverage gaps (Trp single-codon,
the 61-sense-codon set, the empty-enumerable branch) were closed in-session.

**Primary attribution.** RSCU was introduced by **Sharp, Tuohy & Mosurski (1986)**, *Nucleic
Acids Res.* 14(13):5125-5143 (the precise citation given by CODON-RSCU-001 and seqinr); the
looser "Sharp & Li (1986)" attribution also circulates (used by the begomovirus restatement and
ANNOT-CODONUSAGE-001) and refers to the same 1986 origin.

## The measure

For codon `j` of amino acid `i`, with `n_i` codons in that amino acid's synonymous family and
`x_{i,j}` occurrences of codon `j`:

    RSCU_{i,j} = n_i · x_{i,j} / Σ_j x_{i,j}

The denominator sums occurrences over **all** synonymous codons of amino acid `i`.
Equivalently, RSCU is the observed codon frequency divided by the frequency expected if all
synonymous codons were used equally.

**Reading the value:** RSCU = 1.0 ⇒ no bias (used exactly as often as expected); > 1.0 ⇒
preferred / over-represented; < 1.0 ⇒ under-represented. Values are bounded in **[0, n_i]**.

**Invariant:** the RSCU values over an observed synonymous family sum to `n_i` (an algebraic
identity of the formula, and a good differential-test oracle).

## Edge behaviour

- **Aggregation.** Codon counts are **pooled across all input reference sequences** before the
  ratio is computed — RSCU is a property of the aggregate codon pool, not an average of
  per-sequence values.
- **Sense codons only.** Stop codons (TAA, TAG, TGA in the standard code) are excluded from the
  families (Biopython `forward_table` convention). *Cross-page nuance:*
  [[codon-rscu-001-evidence|CODON-RSCU-001]] instead describes the repository as treating the three
  stops as a **degeneracy-3 synonymous family** and computing their RSCU like any other; both
  framings agree that stop handling has **no effect on the RSCU of any amino-acid codon**.
- **Single-codon amino acids.** Met (ATG) and Trp (TGG) have `n_i = 1`, so their RSCU is always
  exactly **1.0** regardless of count.
- **Unobserved family.** If a whole synonymous family has zero occurrences the denominator is
  0 and the base definition is undefined; Seqeron reports **0.0** for each member.
- **Genetic code.** Defaults to NCBI Standard table 1; an overload accepts a non-standard
  `GeneticCode` for alternate tables (the six-codon families are Leu, Arg, Ser).
- **Codon counting (`CountCodons`, per CODON-RSCU-001).** The upstream count step reads
  **non-overlapping triplets from offset 0**, ignores trailing 1–2 bases, and **excludes any triplet
  containing a non-`{A,C,G,T}` character** (Kazusa CUTG convention; string overload uppercases first).

## Relation to other codon-usage statistics

RSCU is the shared normalization at the base of the codon-usage-bias family. The **[[codon-adaptation-index|CAI]]**
(Codon Adaptation Index) reuses RSCU-style relative adaptiveness — but referenced to the family
**maximum** codon and reduced to a single geometric-mean gene score — and guards `log(0)` for
unobserved codons with a pseudocount-style adjustment that plain RSCU does **not** apply (Sharp &
Li 1987 use a **0.5 pseudocount** at reference-table build time; Seqeron's CAI implementation uses
a **`1e-6` clamp** at score time — see [[codon-adaptation-index]]). The
**[[effective-number-of-codons|ENC / Nc]]** (Wright 1990) is the reference-free sibling: it reduces a
gene's codon bias to a single number in [20, 61] using codon homozygosity `F̂` built from the same
synonymous-codon frequencies `p_i`. **[[codon-optimization]]** is the family's *rewriting* operation:
it consumes a codon-usage table and substitutes synonymous codons to improve host expression. The
**aggregation / reporting view** is [[codon-stats-001-evidence|CODON-STATS-001]] (`GetStatistics`),
which bundles RSCU, ENC and CAI over one input together with a positional-GC block
(GC1/GC2/GC3, plus **GC3s** = GC of synonymous third positions, Met/Trp/stop excluded, Peden 1999).
Other siblings in `docs/Evidence/` include rare-codon analysis and the raw codon-count table +
TVD distribution comparison ([[codon-usage-comparison]], `CalculateCodonUsage`/`CompareCodonUsage`);
each builds on the same synonymous-family counting that RSCU formalizes.
