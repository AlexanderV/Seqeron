---
type: concept
title: "Protein secondary structure — Chou-Fasman Pα/Pβ/Pt sliding-window propensity profile"
tags: [sequence-statistics, protein, algorithm]
mcp_tools:
  - predict_chou_fasman
sources:
  - docs/algorithms/Statistics/Secondary_Structure_Prediction.md
  - docs/Evidence/SEQ-SECSTRUCT-001-Evidence.md
source_commit: 016f1502327c2e78f7769d623de5278dbb2f8c40
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-secstruct-001-evidence
      evidence: "Test Unit ID: SEQ-SECSTRUCT-001 ... Algorithm: Protein Secondary Structure Prediction — Chou-Fasman conformational propensities (sliding-window profile)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:hydrophobicity-gravy-and-profile
      source: seq-secstruct-001-evidence
      evidence: "SEQ-SECSTRUCT-001 is a generic protein sliding-window mean-propensity profile sharing the exact mechanism SEQ-HYDRO-001 validates: N−W+1 unweighted window means over a per-residue value table, W>N ⇒ empty, unknown residues skip-and-excluded from the average, case-insensitive. Chou-Fasman propensity (Pα/Pβ/Pt) is the conformational analogue of the Kyte-Doolittle hydropathy value averaged over the same window."
      confidence: medium
      status: current
---

# Protein secondary structure — Chou-Fasman Pα/Pβ/Pt sliding-window propensity profile

**Chou-Fasman secondary-structure prediction** scores *how α-helix-, β-sheet- or β-turn-prone*
each residue of a protein is, using the classic **Chou & Fasman (1974, 1978) conformational
propensities** Pα, Pβ, Pt — each the observed/expected frequency of that amino acid in the
corresponding structure in proteins of known 3-D structure (propensity = observed/expected, ×100
in the integer convention). The **SEQ-SECSTRUCT-001** unit
([[seq-secstruct-001-evidence]]) validates the **sliding-window mean-propensity profile**: for a
window length `W`, emit the arithmetic mean of each of Pα, Pβ, Pt over each window.
[[test-unit-registry]] tracks the unit; [[algorithm-validation-evidence]] describes the artifact
pattern.

This is a **protein-property member of the SEQ-\* sequence-statistics family** — the
conformational-propensity analogue of the [[hydrophobicity-gravy-and-profile]] Kyte-Doolittle
hydropathy profile: **same sliding-window-mean machinery** (exactly `N − W + 1` unweighted window
means, `W > N` ⇒ empty, unknown residues skip-and-excluded, case-insensitive), different
per-residue value table. It is a **sequence-only propensity profile**, distinct from the other
protein-feature predictors it sits beside — the segment-calling
[[transmembrane-helix-prediction]] (thresholds a hydropathy window to emit membrane spans), the
charge–hydropathy [[intrinsic-disorder-prediction-top-idp]], and the SEG composition scan
[[protein-low-complexity-seg]] — all of which read the same amino-acid sequence but score a
different structural property.

## Scope — profile, not the classic state machine

The **method under test is the generic mean-propensity profile, not the full Chou-Fasman
nucleation/extension state machine.** The published algorithm nucleates a helix where **4 of any
6 contiguous** residues are formers (**3 of 5** for a sheet), extends the region in both
directions until the average 4-peptide propensity drops below 1, and calls a helix iff the region
is long enough and ⟨Pα⟩ > ⟨Pβ⟩ (helix former-threshold 1.03, strand 1.05; turn
p(t)=∏p_t over 4 positions, cutoff 7.5e-3). This unit does **not** implement that state machine —
it computes the raw windowed means that such a caller would threshold. The **default window of 7**
is therefore an API convenience, **not** a Chou-Fasman constant (the method's own windows are 6
for helix, 5 for sheet); the window length is caller-supplied and tests pass it explicitly.

## The Chou-Fasman propensity table (1978, all 20 residues)

Verbatim Pα / Pβ / Pt (decimal convention; the reference implementation uses the ×100 integers):

| AA | Pα | Pβ | Pt | AA | Pα | Pβ | Pt |
|----|----|----|----|----|----|----|----|
| A | 1.42 | 0.83 | 0.66 | L | 1.21 | 1.30 | 0.59 |
| R | 0.98 | 0.93 | 0.95 | K | **1.14** | 0.74 | 1.01 |
| D | 1.01 | 0.54 | 1.46 | M | 1.45 | 1.05 | 0.60 |
| N | 0.67 | 0.89 | 1.56 | F | 1.13 | 1.38 | 0.60 |
| C | 0.70 | 1.19 | 1.19 | P | 0.57 | 0.55 | 1.52 |
| E | 1.51 | 0.37 | 0.74 | Q | 1.11 | 1.10 | 0.98 |
| G | 0.57 | 0.75 | 1.56 | S | 0.77 | 0.75 | 1.43 |
| H | 1.00 | 0.87 | 0.95 | T | 0.83 | 1.19 | 0.96 |
| I | 1.08 | 1.60 | 0.47 | W | 1.08 | 1.37 | 0.96 |
| V | 1.06 | 1.70 | 0.50 | Y | 0.69 | 1.47 | 1.14 |

Turn *position* frequencies f(i)..f(i+3) exist in the reference source but are **not** consumed by
this mean profile — only Pα/Pβ/Pt are averaged.

## The one contested value — Lys Pα = 1.14

Two academic sources disagree on the lysine α-helix parameter: CSB|SJU lists **1.16**, while the
Przytycka NCBI lecture and the ravihansa3000 reference implementation (integer **114**) list
**1.14**. The unit adopts **1.14** — supported by two independent retrieved sources (one a
reference implementation) versus one for 1.16, and consistent with the integer convention used for
every other residue. The remaining 19 residues are identical across all retrieved sources.

## Canonical oracles

Closed-form over the tabulated propensities (no library run needed):

| Input | W | Helix mean (Pα) | Sheet mean (Pβ) | Turn mean (Pt) |
|-------|---|------|------|------|
| `A` | 1 | **1.42** | 0.83 | 0.66 |
| `K` | 1 | **1.14** | 0.74 | 1.01 |
| `AE` | 2 | **1.465** | 0.60 | 0.70 |
| `AEV` | 3 | **1.330** | 0.9666… | 0.6333… |

Unknown residues excluded: `AXE` window 3 averages **only A and E**. Sliding step 1 yields
**N − W + 1** windows in N-terminus order.

## Contract and assumptions

- **Case-insensitive** — the implementation uppercases; the table is keyed on uppercase.
- **`W > N` / empty / null / non-positive `W`** ⇒ empty result.
- **Unknown-residue handling** — X, B, Z, `*`, gaps carry no propensity and are **excluded** from
  the per-window count/mean; a window of only unknown residues emits nothing. No source specifies
  non-standard-residue behaviour, so this is the documented deterministic contract (not a
  deviation from a mandated rule).
- **No scoring-constant deviation** — every value produced for in-alphabet residues is exactly the
  published 1978 propensity (Lys resolved to the two-source 1.14).

## Implementation binding

The canonical **Statistics-domain spec** (`docs/algorithms/Statistics/Secondary_Structure_Prediction.md`,
unit SEQ-SECSTRUCT-001) is realised by a single static method in the **Analysis assembly**:

- **File:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs`
- **Signature:** `IEnumerable<(double Helix, double Sheet, double Turn)> PredictSecondaryStructure(string proteinSequence, int windowSize = 7)`
- **Default window:** `DefaultWindowSize = 7` (helix hexapeptide + 1 residue of context; a caller
  parameter, **not** a Chou-Fasman constant — see Scope above).
- **Output contract:** a **lazy `yield return` iterator** streaming **one `(Helix, Sheet, Turn)`
  `double`-tuple per window position**, N-terminus → C-terminus; empty enumerable (no exception)
  when the sequence is null/empty, `windowSize < 1`, or `windowSize > length`. The propensity table
  is a private `Dictionary<char,(double Helix,double Sheet,double Turn)>` keyed on the 20 uppercase
  standard residues; lookup miss ⇒ residue skipped and excluded from the window `count`.

**Named invariants** (validated by
`tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceStatistics_PredictSecondaryStructure_Tests.cs`):

| ID | Invariant |
|----|-----------|
| INV-01 | A single-residue window returns that residue's `(Pα, Pβ, Pt)` tuple (mean of one value). |
| INV-02 | Each emitted tuple is the per-component arithmetic mean of its window's known residues. |
| INV-03 | For an all-known length-`n` sequence, `max(0, n − w + 1)` windows are emitted. |
| INV-04 | Case-insensitive (input upper-cased before lookup). |
| INV-05 | Unknown residues excluded from the window `count`/mean; an all-unknown window emits nothing. |
| INV-06 | Null/empty input, `w > n`, or `w < 1` ⇒ empty result, no exception. |

**Complexity:** `O(n·w)` time (`≤ O(n²)`), `O(1)` extra space (streamed, one tuple per position);
table lookup is `O(1)`. No substring search is performed, so the repository suffix tree is not
applicable.

## Scope and reliability

A **sequence-only propensity profile**, not a structure caller — it reports windowed Pα/Pβ/Pt, it
does not emit helix/sheet/turn *segments* (that would threshold this profile through the
nucleation/extension state machine). The original propensities were derived from a small
non-representative sample (**29 proteins**) and the method has limited accuracy (~50–60% Q3,
Wikipedia); the values are nonetheless the formally defined Chou-Fasman parameters. A
[[research-grade-limitations|research-grade]] implementation of the Chou-Fasman propensity method.

## References

Chou P.Y. & Fasman G.D. (1978) *Ann. Rev. Biochem.* 47:251–276 (empirical predictions of protein
conformation, the 1978 propensity table); Chou P.Y. & Fasman G.D. (1974) *Biochemistry*
13(2):222–245 (prediction of protein conformation, nucleation rules). Reproduced by Wikipedia
(Chou–Fasman method), the Kelley bioinfo lecture, CSB|SJU (Jakubowski), the Przytycka NCBI
lecture, the ravihansa3000/ChouFasman reference implementation, and BMC Bioinformatics PMC1780123.
Full citations in [[seq-secstruct-001-evidence]] (not duplicated here).
