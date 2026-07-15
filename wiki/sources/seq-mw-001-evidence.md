---
type: source
title: "Evidence: SEQ-MW-001 (molecular weight — protein / DNA / RNA molecular mass in Daltons)"
tags: [validation, sequence-statistics, protein]
doc_path: docs/Evidence/SEQ-MW-001-Evidence.md
sources:
  - docs/Evidence/SEQ-MW-001-Evidence.md
source_commit: e058738ff312bb90e5022081cf85e0b9da5b67cb
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-MW-001

The validation-evidence artifact for test unit **SEQ-MW-001** — **molecular weight
calculation**: the average-isotopic **molecular mass (Daltons)** of a **protein, DNA, or RNA**
sequence via a single shared polymerization formula. It is one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; [[test-unit-registry]]
tracks the unit. The formula, mass tables, contract, oracles and the one deviation are
synthesized on the concept [[molecular-weight]].

## What this file records

- **Online sources:**
  - **Expasy Compute pI/Mw documentation** (rank 2, SIB) — the protein Mw definition, verbatim:
    "Protein Mw is calculated by the addition of average isotopic masses of amino acids in the
    protein and the average isotopic mass of one water molecule", expressed in **Daltons (Da)**.
    Uses *average* (not monoisotopic) residue masses.
  - **Expasy ProtParam documentation** (rank 2, SIB) — delegation, verbatim: "Molecular weight
    and theoretical pI are calculated as in Compute pI/Mw" (ProtParam Mw = the Compute pI/Mw
    formula).
  - **Expasy FindMod — average masses of amino-acid residues** (rank 2, SIB reference table) —
    *residue* masses in Da (free amino acid minus one water lost on peptide-bond formation):
    Ala 71.0788, Gly 57.0519, Trp 186.2132, Arg 156.1875, Cys 103.1388; **water H₂O = 18.01524 Da**;
    free-amino-acid mass = residue mass + 18.0153 (A: 71.0788 + 18.0153 = 89.0941).
  - **Biopython `Bio/Data/IUPACData.py`** (rank 3, tables "from PubChem") — the average
    **free-amino-acid** `protein_weights` (A 89.0932 … W 204.2252, G 75.0666), the average
    **monophosphate** `unambiguous_dna_weights` (A 331.2218, C 307.1971, G 347.2212, T 322.2085)
    and `unambiguous_rna_weights` (A 347.2212, C 323.1965, G 363.2206, U 324.1813).
  - **Biopython `Bio/SeqUtils/__init__.py` `molecular_weight`** (rank 3) — the shared single-strand
    formula `weight = sum(weight_table[x] for x in seq) − (len(seq) − 1) * water`; `water = 18.0153`
    (average) or `18.010565` (monoisotopic); subtracting `(len−1)·water` removes one water per bond
    (peptide bonds for protein, phosphodiester bonds for nucleic acids); circular ⇒ `weight −= water`
    (one extra ring-closing bond); double-stranded adds the complement strand's single-strand weight
    (error for protein); "Only unambiguous letters are allowed"; nucleotide sequences "are assumed to
    have a 5' phosphate" (built into the monophosphate masses). Docstring worked values
    `molecular_weight("AGC", "DNA")` → 949.61, `("AGC", "RNA")` → 997.61, `("AGC", "protein")` → 249.29.

- **Datasets (re-derived in-session, average tables, water = 18.0153):**
  - `AGC` protein = 89.0932 + 75.0666 + 121.1582 − 2·18.0153 = **249.2874 ≈ 249.29**
  - `AGC` DNA = 331.2218 + 347.2212 + 307.1971 − 2·18.0153 = **949.6095 ≈ 949.61**
  - `AGC` RNA = 347.2212 + 363.2206 + 323.1965 − 2·18.0153 = **997.6077 ≈ 997.61**
  - Single monomer (zero bonds ⇒ free-monomer mass): `G` protein = **75.0666**; `A` DNA = **331.2218**;
    `A` RNA = **347.2212**.

- **Corner cases / failure modes:** **single monomer** ⇒ `(len−1)·water = 0`, result is the free
  monomer mass (free amino acid for protein, monophosphate for nucleotide); **empty/null** ⇒ 0 (no
  monomers, degenerate — sources define ≥1 monomer only); **unknown/ambiguous letters** — Biopython
  raises a lookup error ("Only unambiguous letters are allowed"); **double-stranded protein** is an
  error; the **5' phosphate** is already in the monophosphate masses (not added/removed separately).

## Deviations and assumptions

**One documented deviation — unknown-symbol handling.** Biopython *rejects* letters outside the
weight table (KeyError). The **repository implementation instead *skips* unknown amino-acid /
nucleotide symbols** — they contribute no mass and no bond, so the reported mass reflects only
recognized monomers and every value stays source-backed (no invented "average" mass). Consistent
with the sibling SEQ-\* methods' skip-unknown / non-throwing convention (cf. the same deviation on
[[hydrophobicity-gravy-and-profile]] and the CountOther routing on [[base-composition]]).

**Two API-shape assumptions:** (1) **input normalization** — free-form `string` is upper-cased
(`ToUpperInvariant`) before table lookup, matching the sibling methods; case folding does not change
any cited numeric value. (2) **unknown-symbol resolution** as above. Both are shape/robustness
choices, not mass-constant changes.

Recommended coverage (from the artifact): MUST — protein `AGC` = 249.29, DNA `AGC` = 949.6095,
RNA `AGC` = 997.6177, single amino acid `G` = 75.0666, single nucleotide `A` DNA 331.2218 / RNA
347.2212, empty/null → 0. SHOULD — case-insensitivity (`agc` == `AGC`); unknown symbol contributes
no mass (the Biopython-reject deviation). COULD — a two-monomer input subtracts exactly one water
(the bond-count invariant). **No source contradictions** — Expasy and Biopython agree on the formula
(protein `Σ residue + water` is algebraically the same `Σ free-amino-acid − (len−1)·water`) and on
average-mass tables to rounding.
