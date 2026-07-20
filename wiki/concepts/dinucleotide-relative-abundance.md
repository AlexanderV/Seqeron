---
type: concept
title: "Dinucleotide composition & relative abundance (Karlin genomic signature)"
tags: [sequence-statistics, composition, algorithm]
sources:
  - docs/algorithms/Statistics/Dinucleotide_Analysis.md
  - docs/Evidence/SEQ-DINUC-001-Evidence.md
source_commit: 42acb6214cd322fa4b68345a3d0292cdac427fd7
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-dinuc-001-evidence
      evidence: "Test Unit ID: SEQ-DINUC-001 ... Algorithm: Dinucleotide Analysis (frequencies, observed/expected relative abundance, codon frequencies)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:cpg-island-detection
      source: seq-dinuc-001-evidence
      evidence: "Gardiner-Garden CpG O/E is the same odds-ratio shape as Karlin's ρ_XY = f_XY/(f_X·f_Y), specialized to the CG dinucleotide; differs only by the N/(N−1) normalization factor"
      confidence: high
      status: current
---

# Dinucleotide composition & relative abundance (Karlin genomic signature)

The **dinucleotide** layer of sequence composition — counting adjacent base pairs and scoring how
over- or under-represented each pair is versus what single-base composition alone would predict.
This is the *2-mer* generalization of [[base-composition|single-base composition]] and the general
case of the **CpG observed/expected ratio** ([[cpg-island-detection]] is the same odds ratio
specialized to the `CG` dinucleotide). Validated under test unit **SEQ-DINUC-001**; the record is
[[seq-dinuc-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

## The two dinucleotide outputs

For a sequence of length `N` there are `N − 1` adjacent dinucleotide (base-step) positions.

1. **Dinucleotide frequency** — normalized frequency `f_XY = count(XY) / (N − 1)` for each observed
   base step. (Karlin's convention: `f_xy` is a *normalized* relative frequency, count over the
   `N − 1` dinucleotide positions, not a raw count.)
2. **Relative abundance / odds ratio (the genomic signature)** —
   `ρ_XY = f_XY / (f_X · f_Y)`, where `f_X`, `f_Y` are the normalized single-base frequencies. This
   is the **Karlin & Burge (1995) dinucleotide relative-abundance** statistic — the "genomic
   signature". `ρ_XY = 1.0` is the no-bias baseline (the frequency expected under statistical
   independence of neighbouring bases); `ρ_XY > 1` is over-represented, `ρ_XY < 1` under-represented.

The method under test is `CalculateDinucleotideRatios`; it returns the **raw ratio**, not a
classification. For *interpretation only*, Karlin & Burge's thresholds (restated by Tsirigos &
Rigoutsos, MBE 19(6):964) call a step **under-represented** at `ρ ≤ 0.78` and **over-represented**
at `ρ ≥ 1.23`.

## Worked oracles (sequence `ATGCGCGT`, length 8)

Mononucleotide counts A=1, T=2, G=3, C=2; dinucleotide counts (7 positions) AT=1, TG=1, GC=2, CG=2,
GT=1. All exact rationals:

| Base step | Frequency `f_XY` | Odds ratio `ρ_XY` |
|-----------|------------------|-------------------|
| GC | 2/7 ≈ 0.28571 | (2/7)/((3/8)(2/8)) = 64/21 ≈ 3.04762 |
| CG | 2/7 ≈ 0.28571 | (2/7)/((2/8)(3/8)) = 64/21 ≈ 3.04762 |
| AT | 1/7 ≈ 0.14286 | (1/7)/((1/8)(2/8)) = 32/7 ≈ 4.57143 |
| TG | 1/7 ≈ 0.14286 | (1/7)/((2/8)(3/8)) = 32/21 ≈ 1.52381 |
| GT | 1/7 ≈ 0.14286 | (1/7)/((3/8)(2/8)) = 32/21 ≈ 1.52381 |

Note `ρ_GC = ρ_CG` and `ρ_TG = ρ_GT` here — the odds ratio is symmetric in the two constituent
base frequencies (`f_X·f_Y = f_Y·f_X`), though `f_XY` itself is directional.

## Relationship to CpG O/E (same shape, N vs N−1)

The Gardiner-Garden & Frommer **CpG O/E** ([[cpg-island-detection]]) is the *same odds-ratio shape*
specialized to `CG`, but normalizes the dinucleotide count by `N` rather than `N − 1`:
`O/E = (#CpG/N) / ((#C/N)·(#G/N))`. The two conventions differ only by the factor `N/(N − 1)`,
which → 1 for long sequences. The repository uses the **Karlin `(N − 1)` normalization** for
`CalculateDinucleotideRatios` and the **Gardiner-Garden `N` normalization** for the dedicated CpG
O/E method — a documented modeling choice, *not* an error, since each convention is authoritative in
its own literature. Both are the dinucleotide-frequency generalization of single-base composition.

## Codon frequency (same source, cross-linked not duplicated)

SEQ-DINUC-001 also validates a **codon-frequency** output: a coding sequence is read as consecutive
**non-overlapping** triplets within a reading frame, and a codon's frequency is `count / total`
counted codons (Kazusa CUTG convention; ambiguous non-ACGT triplets excluded; trailing 1–2 leftover
bases ignored). Oracles: `ATGATGAAA` frame 0 → ATG=2/3, AAA=1/3; frame 1 → TGA=1.0. This is the
**same metric** as `SequenceStatistics.CalculateCodonFrequencies`
([[seq-codon-freq-001-evidence|SEQ-CODON-FREQ-001]], concept [[codon-usage-comparison]]) — see
those pages for the full write-up; it is not re-derived here.

## Method contract (algorithm spec)

The three outputs are three public methods on `SequenceStatistics`
(`Seqeron.Genomics.Analysis/SequenceStatistics.cs`), each returning
`IReadOnlyDictionary<string,double>`:

- `CalculateDinucleotideFrequencies(string)` — normalized dinucleotide frequencies `count / (N − 1)`.
- `CalculateDinucleotideRatios(string)` — the odds ratio `ρ_XY = f_XY / (f_X · f_Y)`; it derives the
  single-base frequencies `f_X` from [[base-composition|`CalculateNucleotideComposition`]] over
  `{A,T,G,C,U}` and the dinucleotide frequencies from `CalculateDinucleotideFrequencies`, so RNA (`U`)
  is supported.
- `CalculateCodonFrequencies(string, int readingFrame = 0)` — non-overlapping triplet frequencies for
  `readingFrame ∈ {0, 1, 2}`.

**Alphabet filters:** dinucleotides are counted only over `{A,T,G,C,U}`; codons only over `{A,T,G,C}`
(ambiguous / non-alphabet units are excluded from both the numerator and the total).

**Complexity:** each method is a single linear counting pass — `O(n)` time, `O(k)` space
(`k` = distinct units observed, `≤ 25` dinucleotides, `≤ 64` codons). These are **tabulation, not
substring search** (no occurrence enumeration), so the repository suffix tree does not apply.

**Frequency invariants (beyond the per-pair `ρ` invariant below):** INV-01 — dinucleotide frequencies
sum to `1.0` (given ≥ 1 valid dinucleotide); INV-02 — codon frequencies sum to `1.0` (given ≥ 1 valid
codon); all returned values are finite and `≥ 0`.

**Single-strand vs strand-symmetrized ρ\* (deviation).** `CalculateDinucleotideRatios` uses
**single-strand** base and dinucleotide frequencies — it is *not* the strand-symmetrized `ρ*` that
Karlin computes by concatenating the sequence with its reverse complement for a true genomic
signature. Consequence: the returned values are single-strand relative abundances, exactly right for
per-sequence O/E (e.g. CpG) but not the double-strand genomic signature. The over-/under-representation
classification (the `0.78` / `1.23` thresholds) is likewise **not implemented** — the caller compares
the raw `ρ` against them.

## Invariants and edge cases

- **INV:** `ρ_XY` is uniquely determined by the three counts (`count(XY)`, `count(X)`, `count(Y)`)
  and `N`; deterministic for a given input.
- **Division-by-zero guard:** when a constituent base is absent (`f_X = 0` ⇒ expected `f_X·f_Y = 0`),
  the ratio is undefined; the method **returns `0`** for that dinucleotide (expected = 0 guard).
- **Case-insensitive:** input is upper-cased before counting.
- **RNA / `U` as a fifth base** — the single-base frequency denominator includes `U` for RNA inputs
  (matching [[base-composition|nucleotide composition]]); RNA dinucleotide signatures are defined
  the same way with `U` replacing `T`. This is a documented assumption, not contradicted by any
  source. See [[seq-dinuc-001-evidence]].
- **Guards:** null / empty / length-<2 input → empty ratios and frequencies; length-<3 → empty
  codon frequencies.

## Relationship to neighbouring composition statistics

- **Single-base layer** — [[base-composition]] is the `f_X` input this page builds the odds ratio
  on; dinucleotide relative abundance asks whether adjacent bases are *independent* of that base
  composition.
- **CpG special case** — [[cpg-island-detection]] (the `CG` odds ratio + island calling).
- **Strand asymmetry** — [[nucleotide-composition-skew]] (`(G−C)/(G+C)`, `(A−T)/(A+T)`) is the
  single-base *asymmetry* cousin; relative abundance is the *neighbour-dependence* view.
- **k-mer diversity** — a dinucleotide is a 2-mer; [[k-mer-statistics]] summarizes the *diversity*
  (Shannon entropy) of a k-mer count profile, whereas this unit scores per-pair *over/under-representation*
  against the independence model — different questions over the same 2-mer counts.
- **Codon layer** — the codon-usage family ([[codon-usage-comparison]],
  [[seq-codon-freq-001-evidence]]) is the *non-overlapping triplet* analog.

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for dinucleotide frequencies,
the Karlin odds ratio, and codon frequencies. The method returns raw ratios, not the
over/under-representation classification (that is documentation-only). No source contradictions —
Karlin (PMC126251), Karlin & Burge 1995 (via MBE 19(6):964), Gardiner-Garden & Frommer 1987, and
the Kazusa CUTG readme are mutually consistent; the two recorded items are documented modeling
choices (`(N−1)` vs `N` normalization; `U` counted for RNA), not contradictions.
