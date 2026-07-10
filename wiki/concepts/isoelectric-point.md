---
type: concept
title: "Isoelectric point (pI) — pH at which a protein's net charge is zero"
tags: [sequence-statistics, protein, algorithm]
mcp_tools:
  - isoelectric_point
sources:
  - docs/Evidence/SEQ-PI-001-Evidence.md
source_commit: 8a4f33ace0d47f5ad2116aa8e775cab5608ccfc2
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-pi-001-evidence
      evidence: "Test Unit ID: SEQ-PI-001 ... Algorithm: Isoelectric Point (pI) Calculation"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:molecular-weight
      source: seq-pi-001-evidence
      evidence: "ExPASy Compute pI/Mw computes pI and molecular weight together in one tool; SEQ-PI-001 and SEQ-MW-001 are the two protein-property scalars that source pairs — pI is the charge/pH calculation MW omits."
      confidence: high
      status: current
---

# Isoelectric point (pI) — pH at which a protein's net charge is zero

The **isoelectric point (pI)** is the **pH at which a protein carries no net electric charge** —
its positive (basic) and negative (acidic) group charges exactly cancel. The **SEQ-PI-001** unit
([[seq-pi-001-evidence]]) validates its calculation by evaluating a per-pH net-charge function and
locating the zero crossing. [[test-unit-registry]] tracks the unit;
[[algorithm-validation-evidence]] describes the artifact pattern.

This is a **protein-property member of the SEQ-\* sequence-statistics family** — a whole-sequence
scalar sibling of [[molecular-weight]] (both are produced together by **ExPASy Compute pI/Mw**;
pI is the charge/pH calculation MW omits) and of the [[hydrophobicity-gravy-and-profile]] GRAVY
index. Like GRAVY and MW it depends on **amino-acid composition only, not sequence order** — the
EMBOSS "no electrostatic interactions" assumption means each ionizable group titrates independently,
so any permutation of a sequence yields the identical pI (the composition-only trait it shares with
[[base-composition]]).

## The net-charge function

At a given pH the net charge is the sum over all ionizable groups (Henderson–Hasselbalch, Moore 1985;
Peptides `charge_pI.cpp`):

```
charge(pH) = Σ_basic  +1 / (1 + 10^(pH − pKa))     +  Σ_acidic  −1 / (1 + 10^(pKa − pH))
```

- **Basic / positive groups** (charge → +1 at low pH): **R, K, H**, and the **N-terminus**.
- **Acidic / negative groups** (charge → −1 at high pH): **D, E, C, Y, the C-terminus**.
- Both **termini are added exactly once per chain**. `charge(pH)` is **monotonically decreasing**
  in pH, so it has a single zero — the pI.

## Finding the zero: bisection over [0, 14]

pI is the root of `charge(pH) = 0`. The implementation brackets it on the **[0, 14] pH interval**
(EMBOSS `iep`) and **bisects to ±0.01**. Because charge is monotonic, bisection converges to the
unique crossing, and the reported pI is always within the **[0, 14] bounds** (an invariant for any
input).

## The EMBOSS pKa scale (`Epk.dat`)

The single-pKa-per-residue constants (EMBOSS `iep`; **the scale this implementation adopts**):

| Group | pKa | Sign | | Group | pKa | Sign |
|-------|-----|------|-|-------|-----|------|
| N-terminus | 8.6 | + | | E (Glu) | 4.1 | − |
| C-terminus | 3.6 | − | | H (His) | 6.5 | + |
| C (Cys) | 8.5 | − | | K (Lys) | 10.8 | + |
| D (Asp) | 3.9 | − | | R (Arg) | 12.5 | + |
| Y (Tyr) | 10.1 | − | | | | |

**Scale dependence is real.** Peptides lists nine scales (Bjellqvist, Dawson, EMBOSS, Lehninger,
Murray, Rodwell, Sillero, Solomon, Stryer); each gives a slightly different pI. This repository's
single-pKa-per-residue model matches **EMBOSS**, *not* the position-dependent Bjellqvist model used
by ExPASy/seqinr — so the seqinr Bjellqvist worked value (pI of `ACDEFGHIKLMNPQRSTVWY` = 6.78454) is
**not** an expected value here; the EMBOSS scale gives **7.36** for the same sequence.

## Canonical oracles

Net charge (EMBOSS scale, Peptides worked example): `FLPVLAGLTPSIVPKLVCLLTKKC` →
**3.037398 / 2.914112 / 0.7184524** at pH **5 / 7 / 9** (reproduced to 6 dp — validates the formula +
pKa table). Derived pI (EMBOSS, bisection ±0.01):

| Sequence | pI | Note |
|----------|-----|------|
| `A`, `AG` | **6.10** | termini only ⇒ midpoint of N 8.6 / C 3.6 |
| `D` | **3.75** | one acidic side chain + termini |
| `K` | **9.70** | one basic side chain + termini |
| `DDDD` | **3.23** | acidic-dominated (low pI) |
| `KKKK` | **11.27** | basic-dominated (high pI) |
| `FLPVLAGLTPSIVPKLVCLLTKKC` | **9.67** | basic (charge still +0.72 at pH 9) |
| `ACDEFGHIKLMNPQRSTVWY` | **7.36** | one of each residue |

## Contract and assumptions

- **Composition-only / order-independent** — permutations give identical pI (no electrostatic coupling).
- **pI ∈ [0, 14]** for any input (the bisection interval).
- **Empty / null ⇒ 7.0** — an input-guard convention (a real protein always has both termini; no source
  defines pI for a zero-length chain), following the neutral/zero sentinel of its SEQ-\* siblings. A
  non-correctness-affecting guard, not an algorithm output.
- **pKa scale = EMBOSS** — a documented parameter choice, not a deviation; the charge model itself matches
  EMBOSS/Peptides/seqinr/ExPASy exactly (they differ only in pKa constants). **No source contradictions.**

## Scope

A **sequence-only charge/pH scalar** — one number per protein. It reports the neutral-charge pH; it does
not emit a full titration curve or per-residue charge assignment. Accuracy caveat (ExPASy): predictions
for **small and highly basic proteins** can be unreliable. A
[[research-grade-limitations|research-grade]] implementation of the EMBOSS/Henderson–Hasselbalch method.

## References

EMBOSS `iep` (pI = net-charge-zero pH, `Epk.dat` pKa table, [0,14] interval); Osorio et al. 2015 Peptides
`charge_pI.cpp` / `charge()` (Henderson–Hasselbalch net-charge formula, Moore 1985; nine pKa scales;
EMBOSS worked example); seqinr `computePI` (Bjellqvist scale, ExPASy Compute pI algorithm); Bjellqvist
et al. 1993 *Electrophoresis* 14:1023; ExPASy Compute pI/Mw (pairs pI with Mw; small/basic-protein caveat).
Full citations in [[seq-pi-001-evidence]] (not duplicated here).
