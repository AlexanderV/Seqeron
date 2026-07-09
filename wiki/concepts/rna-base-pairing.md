---
type: concept
title: "RNA base pairing (Watson-Crick + G-U wobble) and the miRNA-target duplex"
tags: [mirna, algorithm]
sources:
  - docs/Evidence/MIRNA-PAIR-001-Evidence.md
  - docs/algorithms/MiRNA/MiRNA_Target_Pairing.md
source_commit: da06ef5590c4a92834125cfd650b942f61f9b586
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: mirna-pair-001-evidence
      evidence: "Test Unit ID: MIRNA-PAIR-001 ... Algorithm: MiRNA-Target Pairing Analysis (miRNA-mRNA duplex base pairing) ... Algorithm Group: MiRNA"
      confidence: high
      status: current
---

# RNA base pairing (Watson-Crick + G-U wobble) and the miRNA-target duplex

Deciding **whether two RNA bases can pair**, and building the **antiparallel miRNA-mRNA duplex** from
that rule. This is the **first ingested unit of the MiRNA family** (test unit **MIRNA-PAIR-001**) and
the home of the **shared RNA base-pairing primitive** — the Watson-Crick set {A-U, G-C} plus the
single standard RNA **G-U wobble** — a rule that RNA secondary-structure folding and miRNA targeting
both build on. The record is [[mirna-pair-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

A mature miRNA (~22 nt) guides the silencing complex to mRNA targets primarily through **Watson-Crick
pairing** between the miRNA 5' **seed** (nucleotides 2-8) and complementary sites in the target 3'UTR
(Bartel 2009; Agarwal et al. 2015). The target site is the **reverse complement of the seed read
antiparallel** (Lewis et al. 2005). Beyond canonical A-U / G-C pairing, RNA tolerates the **G-U
"wobble"** pair (Crick 1966), which is weaker and counted **separately** from canonical pairs
(PMC4870184). The unit exposes four operations on `MiRnaAnalyzer`, all case-insensitive and
normalising DNA `T` → RNA `U`:

1. `CanPair(base1, base2)` — pairing predicate.
2. `IsWobblePair(base1, base2)` — G-U wobble predicate.
3. `GetReverseComplement(rnaSequence)` — antiparallel RNA reverse complement (for seed → target motif).
4. `AlignMiRnaToTarget(miRna, target)` — the antiparallel duplex alignment + estimated ΔG.

## 1. The pairing rules (shared primitive)

For two nucleotides, classified after `T` → `U` normalisation:

```
CanPair       ⟺  {A-U, U-A, G-C, C-G, G-U, U-G}
IsWobblePair  ⟺  {G-U, U-G}
```

- **Watson-Crick** = {A-U, U-A, G-C, C-G} — "A pairs with U, C pairs with G" (Agarwal et al. 2015).
- **G-U wobble** = {G-U, U-G} — the only standard RNA wobble (Crick 1966). A-C, C-U, etc. **do not
  pair**; they are mismatches.
- **Invariant:** wobble pairs ⊆ pairable pairs (G-U is both a wobble *and* a valid pair); a
  Watson-Crick pair is never a wobble. Both predicates are `O(1)`.

This {A-U, G-C} + G-U rule is the **same base-pairing primitive** used by RNA secondary-structure
classification (Watson-Crick vs G-U wobble in a fold); this concept is the reusable anchor for it.
The sibling MiRNA-family unit [[pre-mirna-hairpin-detection]] builds directly on this rule — a
precursor stem is scored as consecutive {A-U, G-C} + G-U wobble pairs between the two arms.

## 2. Reverse complement for seed → target matching

`GetReverseComplement` reverses and complements a sequence (A↔U, G↔C, `T`→A complement, output in the
RNA alphabet), length-preserving. It turns a miRNA seed — extracted by the sibling unit
[[seed-sequence-analysis]] (positions 2-8) — into the target motif that pairs with it antiparallel
(Lewis et al. 2005). Oracle: the **hsa-let-7a-5p** seed (pos 2-8) `GAGGUAG` →
complement `CUCCAUC` → reversed `CUACCUC`, trivially checkable from A-U/G-C. Unrecognised bases
complement to `N`; empty/null input → `""`.

## 3. The antiparallel duplex (`AlignMiRnaToTarget`)

For miRNA `m` (5'→3') and target `t` (5'→3'), the ungapped antiparallel duplex pairs miRNA index `i`
with target index `len(t)−1−i`, over the overlap `0..min(len)−1`. Each position is classified and an
alignment symbol emitted:

- `|` **Watson-Crick match**, `:` **G:U wobble**, space **mismatch**.
- **Counts invariant:** `Matches + Mismatches + GUWobbles = min(len(m), len(t))`, `Gaps = 0`
  (ungapped — every overlap position is classified into exactly one class).

Duplex folding free energy `FreeEnergy` (ΔG) is a **simplified** Turner 2004 nearest-neighbor estimate:
`CalculateDuplexEnergy` sums stacking energies **only over runs of consecutive paired positions**
(mismatches break stacking but carry no explicit penalty; loop/bulge/initiation/terminal-mismatch and
AU/GU-end terms are omitted). Because Turner 2004 paired-stack energies at 37 °C are all negative, a
**fully paired duplex has ΔG ≤ 0** and an **all-mismatch duplex has ΔG ≥ 0** (sum of zero stacks) —
the **sign and relative ordering are reliable, not the absolute kcal/mol** (Assumption 1).

Worked oracles:

| miRNA | target | Alignment | Matches | Wobbles | Mismatches |
|-------|--------|-----------|---------|---------|------------|
| `AAAA` | `UUUU` | `\|\|\|\|`  | 4 | 0 | 0 |
| `GGGG` | `UUUU` | `::::`     | 0 | 4 | 0 |
| `AAAA` | `AAAA` | (spaces)  | 0 | 0 | 4 |

## Invariants and edge cases

- **INV:** `CanPair` true ⇔ {A-U,U-A,G-C,C-G,G-U,U-G}; `IsWobblePair` true ⇔ {G-U,U-G}; wobble ⊆ pairable.
- **INV:** `GetReverseComplement` reverses + complements, length-preserving; RNA output.
- **INV:** duplex counts sum to `min(len(m),len(t))`; `Gaps = 0`; fully-WC ΔG ≤ 0, all-mismatch ΔG ≥ 0.
- Empty/null miRNA or target → empty `MiRnaDuplex` (all counts 0, empty strings); empty/null to
  `GetReverseComplement` → `""`; no exceptions for these.
- Unequal lengths → align over the shorter overlap. DNA `T` normalised to `U` before pairing.
- Unknown base → complements to `N`, never pairs.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for the base-pairing rule, seed reverse
complement, and the ungapped antiparallel duplex. **It is a pairing primitive, not target-site
prediction:** seed-context efficacy scoring and site-type (8mer / 7mer / 6mer) classification are out
of scope — those belong to [[mirna-target-site-prediction]] (`MiRnaAnalyzer.FindTargetSites`,
MIRNA-TARGET-001), which **depends on** this concept's `GetReverseComplement` + `AlignMiRnaToTarget`. The **A-opposite-
position-1** preference (Agarwal et al. 2015) is an **Argonaute recognition event, not base pairing**,
so it is deliberately **not** represented here. The **free-energy magnitude is approximate**
(stacking-only) — do not treat it as a calibrated thermodynamic ΔG; the aligner models no bulges,
internal loops, or 3' supplementary/compensatory pairing. This is a single positional alignment, not a
pattern search, so the repository suffix tree is **not** applicable. No source contradictions —
Bartel 2009, Agarwal 2015, PMC4870184, Crick 1966, Lewis 2005, and Turner 2004 are mutually
consistent; the pairing classification is exact and only the ΔG magnitude is an intentional
simplification.
