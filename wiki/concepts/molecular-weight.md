---
type: concept
title: "Molecular weight — average molecular mass of a protein / DNA / RNA sequence (Daltons)"
tags: [sequence-statistics, protein, algorithm]
sources:
  - docs/Evidence/SEQ-MW-001-Evidence.md
source_commit: e058738ff312bb90e5022081cf85e0b9da5b67cb
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-mw-001-evidence
      evidence: "Test Unit ID: SEQ-MW-001 ... Algorithm: Molecular Weight Calculation (protein and nucleotide)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:base-composition
      source: seq-mw-001-evidence
      evidence: "The single-strand formula `sum(weight_table[x] for x in seq) − (len−1)·water` is the per-monomer mass sum over the same {A,T,G,C,U}+amino-acid alphabets that base-composition counts — molecular weight is the mass-weighted composition, sharing the case-fold + skip-unknown contract."
      confidence: high
      status: current
---

# Molecular weight — average molecular mass of a protein / DNA / RNA sequence (Daltons)

**Molecular weight (MW)** is the average-isotopic **molecular mass, in Daltons (Da)**, of a
biological polymer — a **protein, DNA, or RNA** sequence. The **SEQ-MW-001** unit
([[seq-mw-001-evidence]]) validates one method that handles all three molecule types through a
single shared **polymerization formula**. [[test-unit-registry]] tracks the unit;
[[algorithm-validation-evidence]] describes the artifact pattern.

This is another **whole-sequence scalar member of the SEQ-\* sequence-statistics family** — the
mass counterpart of the nucleotide [[base-composition]] tally and a sibling of the protein-only
[[hydrophobicity-gravy-and-profile]] GRAVY scalar. Where composition *counts* monomers and GRAVY
*averages* a per-residue property, MW is the **mass-weighted sum** of the same per-monomer tally,
minus the water lost when the monomers polymerize.

## The single shared formula

For any monomer sequence the mass is (Biopython `Bio.SeqUtils.molecular_weight`):

```
weight = Σ (mass of each monomer) − (len − 1) · water
```

Subtracting `(len − 1) · water` removes **one water molecule per bond formed** — peptide bonds for
protein, phosphodiester bonds for nucleic acids (`water = 18.0153` Da average, `18.010565`
monoisotopic). For **protein** this is algebraically identical to the **Expasy definition**
`Σ (residue masses) + water` (Compute pI/Mw, delegated to by ProtParam), because a residue mass is
already the free amino acid minus one water.

- **Single monomer** ⇒ zero bonds ⇒ `(len−1)·water = 0` ⇒ the result is the **free monomer mass**
  (free amino acid for protein; monophosphate for a nucleotide).
- **Circular** nucleic acid ⇒ `weight −= water` once more (the extra bond that closes the ring).
- **Double-stranded** ⇒ add the complement strand's single-strand weight; **error for protein**.

## Mass tables (average, Da)

- **Protein** — *free-amino-acid* masses (Biopython `IUPACData.protein_weights`, "from PubChem"):
  A 89.0932, C 121.1582, D 133.1027, E 147.1293, F 165.1891, **G 75.0666**, H 155.1546, I/L 131.1729,
  K 146.1876, M 149.2113, N 132.1179, P 115.1305, Q 146.1445, R 174.201, S 105.0926, T 119.1192,
  V 117.1463, W 204.2252, Y 181.1885. (Expasy FindMod gives the corresponding *residue* masses —
  free-amino-acid − 18.0153; e.g. Ala residue 71.0788 + 18.0153 = 89.0941.)
- **DNA** — *5'-monophosphate deoxynucleotide* masses: A 331.2218, C 307.1971, G 347.2212, T 322.2085.
- **RNA** — *5'-monophosphate nucleotide* masses: A 347.2212, C 323.1965, G 363.2206, U 324.1813.
- Nucleotide sequences are **assumed to carry a 5' phosphate** (built into the monophosphate masses;
  no terminal phosphate is added or removed beyond this).

## Canonical oracles

Closed-form over the tables (average, water = 18.0153) — Biopython docstring worked values:

| Input | Type | MW (Da) | Derivation |
|-------|------|---------|-----------|
| `AGC` | protein | **249.29** | 89.0932 + 75.0666 + 121.1582 − 2·18.0153 |
| `AGC` | DNA | **949.61** | 331.2218 + 347.2212 + 307.1971 − 2·18.0153 |
| `AGC` | RNA | **997.61** | 347.2212 + 363.2206 + 323.1965 − 2·18.0153 |
| `G`   | protein | **75.0666** | free amino acid (zero bonds) |
| `A`   | DNA | **331.2218** | monophosphate (zero bonds) |
| `A`   | RNA | **347.2212** | monophosphate (zero bonds) |

## Contract and the one deviation

- **Case-insensitive** — the implementation uppercases (`ToUpperInvariant`); tables are keyed on
  uppercase, so `agc` == `AGC` (no numeric change).
- **Empty / null** ⇒ **0** (no monomers; the sources define ≥1 monomer only).
- **Unknown-symbol handling — the sole deviation.** Biopython *rejects* letters outside the weight
  table ("Only unambiguous letters are allowed" → lookup error). This repository instead **skips**
  unrecognized amino-acid / nucleotide symbols — they contribute **no mass and no bond**, so the
  reported mass reflects only recognized monomers and stays entirely source-backed (no invented
  "average" mass). This is the **same skip-unknown / non-throwing convention** as its SEQ-\* siblings
  (cf. GRAVY skipping unknown residues on [[hydrophobicity-gravy-and-profile]], and `CountOther`
  routing on [[base-composition]]); it is an API-shape/robustness choice, not a mass-constant change —
  every value produced for in-alphabet monomers is **exactly source-conformant**.

## Scope

A **sequence-only mass scalar** — one number per sequence. It reports mass; it does not derive
charge or isoelectric point (pI), which the Expasy Compute pI/Mw source pairs with Mw but which is a
separate calculation not in this unit. A [[research-grade-limitations|research-grade]] implementation
of the Expasy/Biopython average-mass method.

## References

Gasteiger et al. 2005 (Expasy Compute pI/Mw — protein Mw = Σ average residue masses + one water);
Expasy ProtParam (delegates Mw to Compute pI/Mw); Expasy FindMod (average residue-mass table +
H₂O 18.01524); Biopython `Bio.SeqUtils.molecular_weight` and `Bio.Data.IUPACData` (the shared
`Σ mass − (len−1)·water` formula and the average protein/DNA/RNA mass tables, "from PubChem"). Full
citations in [[seq-mw-001-evidence]] (not duplicated here).
