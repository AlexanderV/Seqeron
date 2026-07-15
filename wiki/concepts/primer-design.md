---
type: concept
title: "PCR primer-pair design pipeline (heuristic candidate selection + pair compatibility)"
tags: [primer, algorithm]
sources:
  - docs/algorithms/MolTools/Primer_Design.md
source_commit: a51628fa98adbfe52e3789aa7a8fc13d0474ab9b
created: 2026-07-13
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:primer3-weighted-penalty-objective
      source: primer-design-001-report
      evidence: "docs/algorithms/MolTools/Primer_Design.md §5.3: pair selection is greedy and driven by the legacy additive Score (CalculatePrimerScore), explicitly 'rather than a full pairwise optimization' and 'Primer3-style pair penalties … are not represented' — the pipeline that the CalculatePrimer3Penalty objective is the future replacement for."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:primer-dimer-thermodynamics-tm
      source: primer-design-001-report
      evidence: "docs/algorithms/MolTools/Primer_Design.md §2.4 INV-02 / §4.1 step 6 / §7.2: pair compatibility requires |Tm_f − Tm_r| ≤ 5 °C AND !HasPrimerDimer — consuming the PRIMER-TM-001 melting-temperature and dimer signals as prerequisites (PRIMER-TM-001, PRIMER-STRUCT-001)."
      confidence: high
      status: current
---

# PCR primer-pair design pipeline (heuristic candidate selection + pair compatibility)

The **end-to-end PCR primer-pair designer** (test unit
[[primer-design-001-report|PRIMER-DESIGN-001]],
`PrimerDesigner.DesignPrimers`): given a template and a target region, enumerate candidate
oligos in the flanking regions, score each with a per-primer heuristic, **greedily** pick the
best forward and best reverse candidate independently, then accept the pair only if the two
primers are **compatible** (Tm agreement + no primer-dimer). It is a fast heuristic selector,
**not** a Primer3-style combinatorial pair optimizer (`Implementation Status: Simplified`).
Primary spec: `docs/algorithms/MolTools/Primer_Design.md`.

This is the **orchestrator** of the PCR primer-design (`PRIMER-*` / MolTools) family: it
*consumes* the per-primer melting-temperature and primer-dimer signals of
[[primer-dimer-thermodynamics-tm]] (PRIMER-TM-001 / PRIMER-STRUCT-001) as its compatibility
gate, and it is the pair-level pipeline that the validated single-oligo
[[primer3-weighted-penalty-objective]] penalty is designed to eventually replace (see
[[#Relation to the Primer3 penalty]] below). It is the primer analogue of the single-oligo
[[hybridization-probe-design|hybridization-probe designer]] (a *probe* detects, a *primer pair*
amplifies — different geometry, different scorer). [[algorithm-validation-evidence]] describes
the artifact pattern.

## The four-stage pipeline

`DesignPrimers(DnaSequence template, int targetStart, int targetEnd, PrimerParameters?)`
runs four stages (§4.1):

1. **Search-region definition** — a forward region up to **200 bp upstream** of `targetStart`
   and a reverse region up to **200 bp downstream** of `targetEnd`.
2. **Candidate generation** — enumerate every window in the configured length range
   (default **18–25 bp**) across each search region. Reverse-region candidates are
   **reverse-complemented before evaluation** so they are scored in primer (5′→3′) orientation.
3. **Greedy per-side selection** — score every candidate with `EvaluatePrimer` and keep the
   single highest-scoring forward and highest-scoring reverse candidate, chosen **independently**.
4. **Pair compatibility check** — accept the pair only if both conditions hold (INV-02):
   `|Tm_f − Tm_r| ≤ 5 °C` **and** `!HasPrimerDimer(forward, reverse)`.

`ProductSize = reverse.Position + reverse.Sequence.Length − forward.Position` (INV-03),
computed directly from the two chosen candidates.

## Per-primer heuristic score

`EvaluatePrimer` screens a single candidate on GC content, Tm, homopolymers, dinucleotide
repeats, hairpin potential, and 3′ stability — the boolean/scalar
[[primer-structure-qc-screens]] (PRIMER-STRUCT-001) — and `CalculatePrimerScore` collapses those into a
**fixed additive heuristic** (§2.2, **higher is better**, base 100):

```
score = 100 − 2·|length − optimalLength|
            − 2·|Tm − optimalTm|
            − 0.5·|GC% − 50|
            − 5·homopolymerLength
            + bonus_GC-clamp        (+5 when the final base is G or C)
```

3′ stability is calculated and exposed as an additional diagnostic; when
`Check3PrimeStability` is enabled the evaluator can emit a threshold-based `ΔG < −9` issue, but
final ranking is still dominated by this additive score plus the later compatibility checks.
`Avoid3PrimeGC` is **off by default**; when enabled it applies a GC-clamp-style check flagging
primers whose final **two** bases contain no G or C.

## Relation to the Primer3 penalty

This pipeline's ranking is the **legacy** additive `Score` above
(`CalculatePrimerScore`) — **not** the validated Primer3 objective
[[primer3-weighted-penalty-objective]] (`CalculatePrimer3Penalty`). They differ in direction
and shape: the Primer3 penalty is a *lower-is-better* weighted sum of sign-gated one-sided
deviations that collapses to `|Tm−60| + |len−20|` under defaults, whereas this heuristic is a
*higher-is-better* 100-based score with a GC-distance term, a homopolymer penalty, and a
GC-clamp bonus. Per that concept's own scope note, pair ranking is "a separate (future)
concern" and `CalculatePrimerScore` is left unchanged for backward compatibility — so
**PRIMER-DESIGN-001 is exactly that separate pair-level pipeline**, still driven by the legacy
score rather than the Primer3-anchored objective.

## Contract and defaults

`PrimerDesigner.DefaultParameters` (§3.1): length **18–25 bp** (`OptimalLength = 20`), GC
**40–60 %**, Tm **57–63 °C** (`OptimalTm = 60`), `MaxHomopolymer = 4`,
`MaxDinucleotideRepeats = 4`, `Avoid3PrimeGC = false`, `Check3PrimeStability = true`.

- **Preconditions** — `DesignPrimers` throws `ArgumentException` on an invalid target region
  (`targetStart ≥ 0`, `targetEnd < template.Length`, `targetStart < targetEnd`).
- **INV-01** — returns `IsValid = false` (null candidates) when either side has no valid
  candidate.
- Large Tm difference (> 5 °C) or a detected primer-dimer ⇒ `IsValid = false`.

Return `PrimerPairResult`: `Forward` / `Reverse` (`PrimerCandidate?`), `IsValid`, `Message`,
`ProductSize`.

**Complexity:** `DesignPrimers` is **O(n²)** time / **O(k)** space (candidate enumeration over
positions × lengths); `EvaluatePrimer` is **O(n)** / **O(1)**.

## Scope and simplifications

Intentionally simplified (§5.3, §6.2):

- **Greedy, not global** — the best individual forward and reverse primers are *not* guaranteed
  to be the globally best pair (no pairwise optimization over the candidate cross-product).
- **Fixed additive heuristic** — Primer3-style pair penalties and richer thermodynamic pair
  interactions are not represented; 3′ stability is a simple issue flag, not a pair model.
- **Not implemented** — full Primer3 combinatorial pair optimization, richer laboratory
  constraints, and genome-wide specificity analysis; the spec directs users to external
  primer-design tools when those are required. It inherits the simplified Tm and structure
  screens defined elsewhere in the repository.

**Deviation (accepted).** The current source defaults narrow the acceptable Tm to
**57–63 °C**, versus the broader **55–65 °C** summary in the original narrative document —
confirmed from `PrimerDesigner.DefaultParameters`; default filtering is slightly narrower than
the prose. No source contradictions.

**Implementation:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs`
(`DesignPrimers`, `EvaluatePrimer`, `CalculatePrimerScore`). References: Primer3 manual, Addgene
primer-design guidance, Wikipedia (Primer), SantaLucia (1998).
