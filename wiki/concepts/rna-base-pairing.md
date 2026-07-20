---
type: concept
title: "RNA base pairing (Watson-Crick + G-U wobble) and the miRNA-target duplex"
tags: [mirna, algorithm]
mcp_tools:
  - can_pair
  - is_wobble_pair
sources:
  - docs/Evidence/MIRNA-PAIR-001-Evidence.md
  - docs/Evidence/RNA-PAIR-001-Evidence.md
  - docs/Evidence/SEQ-RNACOMP-001-Evidence.md
  - docs/algorithms/MiRNA/MiRNA_Target_Pairing.md
  - docs/algorithms/RnaStructure/RNA_Base_Pairing.md
  - docs/algorithms/Sequence_Composition/RNA_Complement.md
  - docs/Validation/reports/MIRNA-PAIR-001.md
source_commit: 19b199ec6f2ace46231d61652fbd7b01e44caa6d
created: 2026-07-09
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: mirna-pair-001-evidence
      evidence: "Test Unit ID: MIRNA-PAIR-001 ... Algorithm: MiRNA-Target Pairing Analysis (miRNA-mRNA duplex base pairing) ... Algorithm Group: MiRNA"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: rna-pair-001-evidence
      evidence: "Test Unit ID: RNA-PAIR-001 ... Algorithm: RNA Base Pairing (CanPair / GetBasePairType / GetComplement) — the RNA-secondary-structure family's base-pairing primitive"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-rnacomp-001-evidence
      evidence: "Test Unit ID: SEQ-RNACOMP-001 ... Algorithm: RNA-specific Complement (per-base, IUPAC-complete)"
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
The **notation** in which a fold's pairs are written and validated — dot-bracket / extended WUSS — is
a distinct, purely syntactic layer covered by [[rna-dot-bracket-notation]].
The sibling MiRNA-family unit [[pre-mirna-hairpin-detection]] builds directly on this rule — a
precursor stem is scored as consecutive {A-U, G-C} + G-U wobble pairs between the two arms.

### The RNA-secondary-structure family's own copy (RNA-PAIR-001)

The RNA-secondary-structure family exposes the **same primitive** on `RnaSecondaryStructure` as its
own test unit **RNA-PAIR-001** ([[rna-pair-001-evidence]]) — same {A-U, G-C} + G-U wobble rule, same
`T`→`U` normalisation, same case-insensitivity — with two shape differences from the `MiRnaAnalyzer`
surface:

- `CanPair(base1, base2)` — the pairing predicate (identical semantics).
- `GetBasePairType(base1, base2)` — a **typed classifier** returning `WatsonCrick` for {A-U, U-A,
  G-C, C-G}, `Wobble` for {G-U, U-G}, and `null` for every non-pair. This makes the
  Watson-Crick-vs-wobble distinction a first-class return value (where the miRNA surface splits it
  across `CanPair` + `IsWobblePair`). `GetBasePairType('G','U')` must be `Wobble`, never
  `WatsonCrick` (Crick 1966).
- `GetComplement(base)` — the single-base RNA complement (A→U, U→A, G→C, C→G, **T→A** as a DNA input,
  degenerate IUPAC preserved: N→N, R→Y, Y→R), the base-level counterpart of the sequence-level
  `GetReverseComplement`.

Both surfaces are **symmetric** (`f(x,y) == f(y,x)`) and both treat non-alphabet input as a non-pair
(false / null, no exception). This typed classifier is the pairing chemistry that
[[rna-free-energy-turner-model]] assigns stacking energies to and that [[rna-dot-bracket-notation]]
folds encode; no source contradictions between the miRNA and RNA-structure copies of the rule.

**Implementation (per the RnaStructure spec `RNA_Base_Pairing.md`).** The two predicates live on
`RnaSecondaryStructure` (`Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs`) and share a **single
precomputed `byte[128*128]` lookup** (`PairLookup`) indexed by `(b1, b2)` after upper-casing —
value `0` = no pair, `1` = WatsonCrick, `2` = Wobble. This makes classification **branch-free O(1)**
(a fixed ~16 KB byte array) and guarantees the symmetry invariants (`INV-01`/`INV-02`:
`f(x,y) == f(y,x)`) and `INV-03` (`CanPair(x,y) == (GetBasePairType(x,y) != null)`) **by construction**
from one symmetrically-seeded table. The bounds check is the source of the "non-pair for junk input"
behaviour: **any char outside the 0–127 ASCII range returns `false`/`null`** with no exception (the
index is range-checked before the lookup). Note the RnaStructure surface treats a **DNA `T` in
`CanPair`/`GetBasePairType` as a non-RNA base that does NOT pair** (returns `false`/`null`) — this is
the deliberate divergence from the `MiRnaAnalyzer` predicates, which normalise `T`→`U` before pairing
(see the MIRNA-PAIR-001 fix note below); only `GetComplement` on this surface treats `T` as `U`
(complement `A`). `GetComplement(char)` is a thin **delegate to
`SequenceExtensions.GetRnaComplementBase`** (Core), the IUPAC-complete complement cross-verified
against Biopython `complement_rna` (the SEQ-RNACOMP-001 copy described above).

**Not implemented (out of scope by spec).** Only the standard four-letter RNA alphabet is modelled:
non-canonical pairs beyond G-U wobble — **Hoogsteen pairs, sheared G•A, and the inosine wobble pairs
I•U / I•A / I•C** — are deliberately excluded, consistent with the nearest-neighbour folding model the
sibling RNA-secondary-structure methods use. `GetBasePairType` has **no "intentionally simplified"
terms** — the pairing classification is exact; only downstream ΔG magnitude (below) is approximate.

### The SEQ-family full-IUPAC RNA complement (SEQ-RNACOMP-001)

A **third** copy of the RNA base complement lives on the SEQ-\* sequence-utility surface as
`SequenceExtensions.GetRnaComplementBase(char)` (test unit **SEQ-RNACOMP-001**,
[[seq-rnacomp-001-evidence]]; primary spec `docs/algorithms/Sequence_Composition/RNA_Complement.md`) —
the RNA sibling of the DNA per-base `GetComplementBase` (SEQ-COMP-001, not yet ingested). It is a
single-`char` → single-`char` operation, implemented in
`Seqeron.Genomics.Core/SequenceExtensions.cs` as an `[AggressiveInlining]` `switch` expression over a
fixed alphabet — **O(1) time / O(1) space**, no lookup table, no search (the repository suffix tree is
N/A). It is the same base chemistry (A→U, U→A, G→C, C→G, **T→A** with DNA T treated as U) but
**IUPAC-complete**: it maps *every* ambiguity code, not just the canonical bases and the handful
(N/R/Y) the `RnaSecondaryStructure` `GetComplement` above documents. The complement table is
Biopython's `ambiguous_rna_complement`, which is **identical to `ambiguous_dna_complement` except the
alphabet swaps T→U**:

- **Reciprocal codes:** A↔U, C↔G, R↔Y, M↔K, D↔H, B↔V.
- **Self-complementary:** W→W, S→S, N→N.
- **Pass-through:** any character outside the IUPAC nucleotide set (gaps `-`/`.`, digits, `Z`,
  whitespace, and `X`) is returned **unchanged, including its original case** — never an exception for
  any `char` input.
- **Casing (only divergence from Biopython):** recognized bases/codes return **uppercase** (repo
  normalize-to-uppercase convention, mirroring SEQ-COMP-001 MUST-02); Biopython preserves input case.
  Casing-only — the complement identity is unchanged; non-IUPAC pass-through characters keep their case.

Complement is an **involution** on the canonical bases and ambiguity codes (within the U-alphabet;
`T` is absorbed into U, so `complement(complement(T)) = U`, not T). It is distinct from the DNA path:
`GetRnaComplementBase('A') = 'U'` vs `GetComplementBase('A') = 'T'`. Worked oracle — Biopython's
`"ACGTUacgtuXYZxyz"` forward-complements per base to `"UGCAAugcaaXRZxrz"` in Biopython (case
preserved); under the repo uppercase convention recognized bases uppercase while non-IUPAC `X`/`Z`/`x`/`z`
pass through verbatim, giving `"UGCAAUGCAAXRZxRz"`. Sources: Biopython `IUPACData.py` / `Seq.py` / API
docs, bioinformatics.org SMS, and the NC-IUB 1984 standard (Cornish-Bowden 1985) — all mutually
consistent.

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

The two-stage validation verdict is recorded in [[mirna-pair-001-report]] — **Stage A PASS,
Stage B PASS-WITH-NOTES, End-state CLEAN** (full suite 6543/0). The review found and fixed a real
gap: `CanPair`/`IsWobblePair` did not honour the documented DNA-`T`→`U` contract
(`CanPair('A','T')` returned false while `AlignMiRnaToTarget`/`GetReverseComplement` normalise T→U),
and the canonical test was green-washed by silently omitting the `A-T` assertion. A private
`NormalizeBase` was added so the predicates (and the MCP tools that delegate to them) treat T as U,
with the tests rewritten to the sourced contract values.
