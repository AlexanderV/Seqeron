---
type: concept
title: "Primer structure QC screens (hairpin / dimer / 3'-stability / run heuristics)"
tags: [primer, algorithm]
sources:
  - docs/algorithms/MolTools/Primer_Structure_Analysis.md
source_commit: 208c7e40e603985d2f569e45000b95b3201a934d
created: 2026-07-13
updated: 2026-07-15
graph:
  relationships:
    - predicate: alternative_to
      object: concept:primer-dimer-thermodynamics-tm
      source: primer-struct-001-report
      evidence: "docs/algorithms/MolTools/Primer_Structure_Analysis.md §1/§5.3: 'exposes discrete boolean or scalar quality signals rather than a full thermodynamic folding model'; hairpin/dimer are boolean structural screens, explicitly NOT the full free-energy folding / ntthal Tm of PRIMER-TM-001."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:primer-design
      source: primer-struct-001-report
      evidence: "docs/algorithms/MolTools/Primer_Structure_Analysis.md §5.1 methods (HasHairpinPotential, HasPrimerDimer, Calculate3PrimeStability, FindLongestHomopolymer, FindLongestDinucleotideRepeat) are the per-candidate screens consumed by EvaluatePrimer / DesignPrimers (PRIMER-DESIGN-001)."
      confidence: high
      status: current
---

# Primer structure QC screens (hairpin / dimer / 3'-stability / run heuristics)

The **boolean/scalar secondary-structure QC screens** a single PCR primer is checked against
(test unit [[primer-struct-001-report|PRIMER-STRUCT-001]], `Implementation Status: Simplified`): five fast heuristics on
`PrimerDesigner` — hairpin potential, self/hetero primer-dimer, 3′-end stability, longest
homopolymer run, and longest dinucleotide repeat. The unit deliberately **exposes discrete
boolean or scalar quality signals rather than a full thermodynamic folding model** — it is the
**screening** counterpart, not the thermodynamic engine. Primary spec:
`docs/algorithms/MolTools/Primer_Structure_Analysis.md`.

These are exactly the per-candidate screens that the [[primer-design]] pipeline's `EvaluatePrimer`
/ `CalculatePrimerScore` consume (and that its pair-compatibility gate reads `HasPrimerDimer`
from). They are the **low-fidelity alternative** to the full thermodynamic path
[[primer-dimer-thermodynamics-tm]] (PRIMER-TM-001): where that engine runs the ntthal DP over the
SantaLucia NN model to return a dimer/hairpin **Tm**, these methods return a **yes/no flag** (or a
scalar ΔG / run length). [[algorithm-validation-evidence]] describes the artifact pattern.

## The five screens

| Method | Returns | Core rule | Default thresholds |
|--------|---------|-----------|--------------------|
| `HasHairpinPotential(seq, minStemLength, minLoopLength)` | `bool` | self-complementary stem separated by a loop | `minStemLength = 4`, `minLoopLength = 3` |
| `HasPrimerDimer(p1, p2, minComplementarity)` | `bool` | terminal complementarity of p1's 3′ window vs `revcomp(p2)` | `minComplementarity = 4` |
| `Calculate3PrimeStability(seq)` | `double` (kcal/mol ΔG) | NN ΔG°37 of the terminal 5-mer + initiation | last 5 bases only |
| `FindLongestHomopolymer(seq)` | `int` | longest mononucleotide run | — |
| `FindLongestDinucleotideRepeat(seq)` | `int` | longest repeated-dinucleotide count | — |

All string methods **normalise to uppercase** before comparison. `HasHairpinPotential` /
`HasPrimerDimer` return `false` on null/empty input; `Calculate3PrimeStability` returns `0` for
null/empty/`< 5` nt; `FindLongestHomopolymer` returns `0` for empty; `FindLongestDinucleotideRepeat`
returns `0` for `< 4` nt.

## Hairpin: boolean stem-loop screen, length-dispatched

A hairpin is a self-complementary **stem** of `≥ minStemLength` bp separated by a **loop** of
`≥ minLoopLength` nt. The screen is **length-dispatched** on the primer:

| Sequence length | Strategy | Rationale |
|-----------------|----------|-----------|
| `< 100 bp` | nested-loop search, `O(n²)` time / `O(1)` aux | low overhead for typical 18–25 bp primers (the intended path) |
| `≥ 100 bp` | suffix-tree-assisted search | avoids the quadratic scan on long inputs |

`HasHairpinPotential` returns `false` whenever `len < 2·minStemLength + minLoopLength` (**INV-01** —
a valid stem-loop cannot fit). A perfect palindrome *may or may not* form a hairpin: loop
feasibility still governs.

## Primer-dimer: terminal 3′ complementarity

`HasPrimerDimer` reverse-complements the second primer and compares the **last `min(8, len1, len2)`
bases** of the first primer to the start of `revcomp(p2)`, flagging a dimer when the complementary
run meets `minComplementarity` (default 4). This models the classic mechanism — polymerase-extensible
complementary **3′ ends** — as an **`O(n)`** window comparison, and deliberately does **not** model
non-terminal or context-dependent duplex interactions (that is the ntthal path in
[[primer-dimer-thermodynamics-tm]]).

## 3′-end stability: SantaLucia-1998 NN ΔG of the terminal 5-mer

`Calculate3PrimeStability` scores the **final 5 bases** as a nearest-neighbour ΔG°37 sum — the
Primer3 `PRIMER_MAX_END_STABILITY` convention (SantaLucia 1998). It sums the four dinucleotide
steps of the terminal 5-mer from this **distinct SantaLucia-1998 ΔG table** (note: **not** the 2004
`NnUnifiedParams` ΔH°/ΔS° set used by the thermodynamic Tm engine):

| Step | ΔG | Step | ΔG | Step | ΔG |
|------|----|------|----|------|----|
| AA/TT | −1.00 | CA/TG | −1.45 | GA/TC | −1.30 |
| AT | −0.88 | GT/AC | −1.44 | CG | −2.17 |
| TA | −0.58 | CT/AG | −1.28 | GC | −2.24 |
| | | | | GG/CC | −1.84 |

plus **terminal initiation** `+0.98` kcal/mol (terminal G·C) or `+1.03` (terminal A·T). More
negative ⇒ more stable 3′ end (raises mispriming risk). The doc records `GCGCG` as the most stable
5-mer (`−6.86`) and `TATAT` the least (`−0.86`). `O(1)` — only the last 5 bases are read;
`< 5` nt returns `0` (**INV-02**). The [[primer-design]] evaluator surfaces this as a `ΔG < −9`
issue flag when `Check3PrimeStability` is on.

## Run heuristics: homopolymer and dinucleotide repeats

`FindLongestHomopolymer` is a single left-to-right scan returning the longest mononucleotide run:
`0` for empty, `1` when all bases are unique, `len` when all bases are identical (**INV-03**).
`FindLongestDinucleotideRepeat` scans for the longest repeated dinucleotide, returning `0` for
`< 4` nt (**INV-04**). Both are `O(n)` / `O(1)` run-length heuristics — they flag slippage /
mispriming liabilities without estimating PCR-yield impact.

## Scope and simplifications (§5.3, §6.2)

- **Boolean/scalar, not thermodynamic** — hairpin and dimer are structural *screens*; they report
  a potential stem-loop or terminal complementarity without ranking a free-energy ensemble or a
  full duplex landscape. Consequence: no ΔG ranking, no full-fold enumeration.
- **Run checks are heuristics** — homopolymer/dinucleotide detection flags repetitive structure
  only.
- **Not implemented** — full secondary-structure thermodynamics for primer hairpins and dimers;
  the spec documents *no in-unit alternative* (the thermodynamic sibling
  [[primer-dimer-thermodynamics-tm]] is the separate PRIMER-TM-001 unit). For typical 18–25 bp
  primers the short-sequence branch is the intended path; the suffix-tree branch is a
  long-sequence optimisation.

No source contradictions.

**Implementation:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs`
(`HasHairpinPotential`, `HasPrimerDimer`, `Calculate3PrimeStability`, `FindLongestHomopolymer`,
`FindLongestDinucleotideRepeat`). References: SantaLucia (1998); Primer3 manual
(`PRIMER_MAX_HAIRPIN_TH`, `PRIMER_MAX_SELF_END`, `PRIMER_MAX_END_STABILITY`, `PRIMER_MAX_POLY_X`);
Wikipedia (Primer dimer, Stem-loop, Nucleic acid thermodynamics).
