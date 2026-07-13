---
type: concept
title: "Hybridization probe design (heuristic candidate generation + additive-penalty ranking)"
tags: [primer, algorithm, validation]
mcp_tools:
  - design_probes
  - design_tiling_probes
  - design_antisense_probes
sources:
  - docs/algorithms/MolTools/Hybridization_Probe_Design.md
source_commit: acf2ff3439b667439bd5b830c84edb1f946b9f4a
created: 2026-07-13
updated: 2026-07-13
graph:
  relationships:
    - predicate: relates_to
      object: concept:taqman-probe-design-rules
      source: hybridization-probe-design
      evidence: "docs/algorithms/MolTools/Hybridization_Probe_Design.md §2.2.1: the TaqMan (5'-nuclease hydrolysis probe) rules are an opt-in extension that 'is separate from and does not alter the generic designer, which remains the default'"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:probe-offtarget-specificity-scan
      source: hybridization-probe-design
      evidence: "The suffix-tree overload calls CheckSpecificity (docs §5.1/§4.2) to filter/rescale candidates for genome-wide uniqueness; the standalone gapped-SW off-target scan is the specificity-checking sibling"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:melting-temperature
      source: hybridization-probe-design
      evidence: "docs §2.2: Tm uses the Wallace rule for lengths < 14 and ThermoConstants.CalculateSaltAdjustedTm for longer probes"
      confidence: high
      status: current
---

# Hybridization probe design (heuristic candidate generation + additive-penalty ranking)

The **generic hybridization-probe designer** (test unit **PROBE-DESIGN-001**,
`ProbeDesigner.DesignProbes`): enumerate every oligonucleotide window in a target within an
application-specific length range, score each with a **fixed additive-penalty heuristic**, keep
the positive-scoring candidates, and return the top `maxProbes` by score. It supports FISH, DNA
microarrays, Northern/Southern blots, and qPCR by swapping **application-specific parameter sets**.
Primary spec: `docs/algorithms/MolTools/Hybridization_Probe_Design.md`.

This is the **default** designer. The [[taqman-probe-design-rules|TaqMan (5'-nuclease
hydrolysis-probe) rules]] are a separate **opt-in** layer that does not alter it, and
[[probe-offtarget-specificity-scan]] is the standalone specificity/off-target scan. It sits in the
MolTools reagent-design family alongside the primer units [[primer-dimer-thermodynamics-tm]] and
[[primer3-weighted-penalty-objective]]; [[algorithm-validation-evidence]] describes the artifact
pattern.

## Distinction from primer design

A **primer** is short (~18–25 nt), selected in a *pair* that flanks an amplicon, and scored by the
Primer3 weighted objective ([[primer3-weighted-penalty-objective]]) with dimer/hairpin
thermodynamics ([[primer-dimer-thermodynamics-tm]]). A **probe** is a *single* oligo (20–500 nt
depending on assay) that must *detect* — not prime — a target, so the design goal is a Tm/GC window,
low self-structure, low repetitiveness, and genome-wide uniqueness, not amplicon geometry. The
scoring machinery is therefore its own: a fixed additive-penalty heuristic rather than the Primer3
sign-gated weighted deviation sum.

## The additive-penalty score

Each candidate window starts at a raw score of **1.0**, reduced by a fixed set of penalties
(docs §2.2):

| Factor | Penalty |
|--------|---------|
| GC content outside the configured range | −0.30 |
| Tm outside the configured range | −0.30 |
| Homopolymer run above the configured max | −0.20 |
| Self-complementarity above the configured max | −0.20 |
| Secondary-structure potential | −0.15 |
| Simple repeats | −0.10 |
| Starts with G/C | −0.02 |
| Ends with G/C | −0.02 |

The score **ranks candidates heuristically; it is not a hybridization probability** (an
intentional simplification — probe quality is a fixed additive penalty, not a thermodynamic
binding model). Candidates with a non-positive raw score are rejected (`Evaluate` returns `null`).

## Melting temperature

Probe Tm follows the repository's two-regime rule (docs §2.2): the **Wallace rule** for oligos
shorter than 14 nt, and a **salt-adjusted longer-probe formula**
`81.5 + 16.6·log₁₀[Na⁺] + 41·GC − 600/N` via `ThermoConstants.CalculateSaltAdjustedTm(...)` for
longer probes — the same shared helper family as [[melting-temperature]], **not** a full
nearest-neighbour duplex calculation (intentional simplification: fine-grained context effects are
not captured).

## Application-specific defaults

`ProbeDesigner.Defaults` supplies length / Tm / GC windows per assay (docs §4.2):

| Application | Length (bp) | Tm (°C) | GC | Notes |
|-------------|-------------|---------|------|-------|
| Microarray | 50–60 | 75–85 | 0.40–0.60 | **Default** parameter set |
| FISH | 200–500 | 70–90 | 0.35–0.65 | Higher self-complementarity tolerance |
| Northern blot | 100–300 | 65–80 | 0.40–0.60 | Intermediate sizes |
| qPCR | 20–30 | 68–72 | 0.40–0.60 | Shortest default probes |
| Southern blot | 150–500 | 65–75 | 0.35–0.65 | Long-probe setting |

## Algorithm and complexity

1. Uppercase the target; **precompute GC prefix sums** for O(1) per-window GC lookup.
2. Enumerate all candidate windows in the configured length range.
3. Evaluate GC, Tm, homopolymers, self-complementarity, secondary structure, repeats, terminal G/C.
4. Keep raw-score-positive candidates, sort by score, return the top `maxProbes` (default 10).
5. **Genome-index overload** (`DesignProbes(string, ISuffixTree, …)`): build a larger raw-score
   shortlist (currently `maxProbes × 5`), then either **filter** it for uniqueness or **rescale**
   shortlisted scores by the specificity value — see below.

`DesignProbes` runs in **O(n × m)** time (n = sequence length, m = scanned length range),
`CheckSpecificity` in **O(m)** per probe via suffix-tree lookups, and `DesignTilingProbes` in
**O(n)** over fixed-length windows (it deliberately includes suboptimal probes when needed for
coverage, reporting coverage, mean Tm, and Tm range). Probe types: `Standard`, `Tiling`,
`Antisense`, `LNA`, `MolecularBeacon`.

## Specificity is a post-shortlist filter, not a full rerank

`CheckSpecificity(probe, index)` maps suffix-tree hit counts to a specificity score: **0** for no
hits, **1** for a unique hit, else **1 / hits** (INV-03). The overload applies this *after* the
raw-score shortlist is already formed (an intentional simplification), so uniqueness-aware results
are a specificity-filtered or specificity-scaled **subset of the top raw-score candidates**, not a
full-candidate rerank. Consequence: when `requireUnique` is `false` and specificity is `0`, a
shortlisted probe can remain in the output with **final score 0**. This is the gap the standalone
gapped-SW scan [[probe-offtarget-specificity-scan]] (PROBE-VALID-001) addresses independently.

## Invariants

- **INV-01** — the base ranking pass retains only raw-score-positive candidates before optional
  specificity rescaling (`DesignProbesOptimized` rejects `score <= 0`).
- **INV-02** — probe GC content is a fraction of length (`gcCount / length`).
- **INV-03** — `CheckSpecificity` returns 0 (no hits) / 1 (unique) / `1/hits` (otherwise).

## Edge cases

Null, empty, or shorter-than-`MinLength` target → **no probes** (explicit early return). Candidate
with score ≤ 0 → rejected. Specificity check with no genome hits → 0.

## Scope and simplifications

Suitable for **fast candidate generation and filtering**, not high-confidence experimental
validation on its own (docs §6.2). Self-complementarity and secondary structure use simple sequence
rules (no explicit folding model). **Not implemented:** database-style alignment / calibrated
hybridization prediction, and MGB / LNA / dual-quencher probe *chemistries* — the LNA Tm-adjustment
that does exist is the [[taqman-probe-design-rules|separate PROBE-DESIGN-001 LNA variant]]; for the
other chemistries the spec directs users to the relevant chemistry's own design tool. One accepted
assumption: the longer-probe Tm is delegated to the shared `CalculateSaltAdjustedTm` helper rather
than an inline formula.
