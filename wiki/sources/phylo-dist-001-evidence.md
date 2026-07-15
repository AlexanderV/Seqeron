---
type: source
title: "Evidence: PHYLO-DIST-001 (Phylogenetic Distance Matrix — p-distance / JC69 / K2P / Hamming)"
tags: [validation, phylogenetics]
doc_path: docs/Evidence/PHYLO-DIST-001-Evidence.md
sources:
  - docs/Evidence/PHYLO-DIST-001-Evidence.md
source_commit: 3a53115ec5fbdbc54448d69550c3b961c40a320a
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PHYLO-DIST-001

The validation-evidence artifact for test unit **PHYLO-DIST-001** — **Phylogenetic Distance Matrix
Calculation**: the pairwise evolutionary-distance primitive that emits the symmetric *n×n* matrix
distance-based tree construction (NJ/UPGMA) consumes. Four methods — **Hamming** (raw count),
**p-distance** (`differences/comparableSites`), **Jukes–Cantor JC69** (`-3/4·ln(1−4p/3)`), and
**Kimura-2-parameter K2P** (`-1/2·ln((1−2S−V)·√(1−2V))`). This is the **third phylogenetics-family
(`PHYLO-*`) Evidence file** (after PHYLO-BOOT-001 and PHYLO-COMP-001) and one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the definitions, formulas,
invariants, corner cases, and worked oracles are synthesized in the dedicated concept
[[evolutionary-distance-matrix]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources (§1):** Wikipedia **Models of DNA evolution** (p-distance; JC69 formula + equal-base
  / equal-rate assumptions; K2P S/V transition-vs-transversion definition; gaps ignored) + Wikipedia
  **Substitution model** (symmetric matrix `d(i,j)=d(j,i)`, zero diagonal, time-reversibility;
  JC69/Kimura primary refs) + Wikipedia **Distance matrices in phylogeny** (distance = fraction of
  aligned mismatches, gaps ignored *or* counted as mismatch by implementation choice, raw Hamming count).
- **Formulas & invariants (§4, §6):** matrix symmetry / zero diagonal / non-negativity / *n×n* shape;
  **correction ordering** `JC69 ≥ p` and `K2P ≥ p` (corrections only inflate the raw proportion);
  Hamming = integer count not a proportion; gaps excluded from `comparableSites`; case-insensitive
  (`ToUpperInvariant`); throws on unequal-length inputs. Triangle inequality *expected in most cases* but
  **not guaranteed** for the corrected distances (flagged, not asserted).
- **Corner cases (§3):** identical → 0; single diff → small positive; **all-gaps / empty → 0** (no
  comparable sites, no differences, `0/n→0`); unequal length → `ArgumentException`; null →
  `ArgumentNullException`; **JC69 saturates at p ≥ 3/4 → +∞**, **K2P saturates at V ≥ 1/2 → +∞**;
  ambiguous IUPAC bases (N/R/Y…) skipped like gaps (pairwise deletion — only A,C,G,T compared).
- **Worked datasets (§5):** `ACGTACGT` vs `TCGTACGT` → Hamming 1 / p 0.125 / JC69 ≈ 0.137; pure
  transition `ACGT`→`GCGT` → K2P ≈ 0.34657; pure transversion `ACGT`→`CCGT` → K2P ≈ 0.31713; mixed
  1-transition+1-transversion → K2P ≈ 0.30679; gap case `ACG-ACGT` vs `ACGTACGT` → 7 comparable sites;
  JC69 `p=0.25 → ≈0.304`.

## Deviations and assumptions

No deviations from the literature — the JC69/K2P formulas, symmetric zero-diagonal matrix, and
saturation limits are the standard textbook definitions. Two documented **API-contract assumptions**:
(1) **empty / all-gap input → distance 0** (mathematical limit `0/n→0`, not NaN or throw); (2)
**pairwise deletion** for gaps *and* ambiguous IUPAC codes (skipped per-position, so `comparableSites`
may be fewer than the alignment length — complete-deletion and gap-as-mismatch are literature
alternatives this unit does not take). No source contradictions. References: Jukes & Cantor 1969,
Kimura 1980, Felsenstein 2004 *Inferring Phylogenies*.
