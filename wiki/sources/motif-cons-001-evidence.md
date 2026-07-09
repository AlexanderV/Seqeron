---
type: source
title: "Evidence: MOTIF-CONS-001 (Consensus from a multiple alignment / most-frequent residue)"
tags: [validation, motif]
doc_path: docs/Evidence/MOTIF-CONS-001-Evidence.md
sources:
  - docs/Evidence/MOTIF-CONS-001-Evidence.md
source_commit: de59ece45cd0b9e5969d6589c1c935e8522d4e4c
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: MOTIF-CONS-001

The validation-evidence artifact for test unit **MOTIF-CONS-001** — **Consensus
Sequence from a Multiple Alignment**, the plurality / most-frequent-residue rule
(`MotifFinder.CreateConsensusFromAlignment(alignedSequences)`). It is one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the
algorithm, its decision rule, invariants, worked oracles, and corner cases are synthesized in
[[consensus-from-alignment]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Wikipedia "Consensus sequence"** (rank 4, citing Schneider & Stephens 1990; Pierce
    2002) — a consensus is *"the calculated sequence of most frequent residues, either
    nucleotide or amino acid, found at each position in a sequence alignment"*; it reduces
    variability to a single residue per position and does **not** specify a tie-break
    (deferring richer representation to sequence logos).
  - **Rosalind "Consensus and Profile" (CONS)** (rank 5, curated datasets) — the 4×n
    **profile matrix** `P[base, j]` = count of that base in column `j`, and the consensus
    rule *"the most common symbol at each position … the symbol with the maximum value in
    the j-th column."* Explicitly notes **ties → multiple valid consensus strings** (so a
    deterministic implementation must fix one rule). Defined only for **equal-length** strings.
  - **EMBOSS `cons`** (rank 3, reference implementation) — parameterised **plurality**
    variant: a residue is called only if its weighted match count exceeds a plurality
    threshold (default = half the total sequence weight), else `n` (nucleotide) / `x`
    (protein) is written. This is the *threshold* alternative that MOTIF-CONS-001 does **not** adopt.
  - **Geneious / LANL HIV database** (rank 5) — documents the tie-break options
    (IUPAC ambiguity code, specified residue order, `?` symbol) and, for the Geneious
    family, an explicit **alphabetical** tie-break: *"the residue letter occurring earlier
    in the alphabet was chosen."*
- **Datasets (documented oracles):**
  - *Rosalind CONS sample*: 7 strings of length 8 (`ATCCAGCT`, `GGGCAACT`, `ATGGATCT`,
    `AAGCAACC`, `TTGGAACT`, `ATGCCATT`, `ATGGCACT`) → profile A=`5 1 0 0 5 5 0 0`,
    C=`0 0 1 4 2 0 6 1`, G=`1 1 6 3 0 1 0 0`, T=`1 5 0 0 0 1 1 6` → consensus **`ATGCAACT`**.
  - *Alphabetical tie-break* (derived): `AT`, `GT` → column 1 has A and G tied at 1;
    alphabetical rule picks **A**; column 2 unanimous T → consensus **`AT`**.

## Deviations and assumptions

Two documented **assumptions** (both determinism choices, correctness-affecting only on
tied columns / out-of-scope thresholds — the rank-5 Rosalind oracle is unaffected):

1. **Alphabetical tie-break (A<C<G<T).** Rosalind permits any most-common symbol on a tie;
   EMBOSS uses scoring/plurality; the Geneious/LANL family documents an explicit
   alphabetical tie-break, adopted here to make the method deterministic (a library
   requirement).
2. **Pure most-frequent, no plurality threshold.** The canonical signature
   `CreateConsensusFromAlignment(alignedSequences)` takes no threshold parameter, matching
   the Rosalind/Wikipedia "most common symbol" definition rather than EMBOSS's plurality;
   threshold-based no-consensus (`n`/`x`) output is out of scope (the area separately exposes
   IUPAC-degenerate consensus via `GenerateConsensus` and a PWM via `CreatePwm`).

Recommended coverage (MUST): Rosalind CONS → `ATGCAACT`; identical sequences → that exact
sequence; alphabetical tie-break A,G → A; single sequence → returned unchanged;
case-insensitive (lowercase normalised via `ToUpperInvariant`). SHOULD: null →
`ArgumentNullException`, empty collection → empty string, unequal-length →
`ArgumentException` (Rosalind equal-length precondition). COULD: non-ACGT → `ArgumentException`
(alphabet validation consistent with `CreatePwm`). No source contradictions — Wikipedia,
Rosalind, EMBOSS, and Geneious/LANL agree on the most-frequent rule and differ only on the
(implementation-defined) tie-break and the optional plurality threshold.
