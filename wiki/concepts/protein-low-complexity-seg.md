---
type: concept
title: "Protein low-complexity region detection (SEG algorithm)"
tags: [analysis, algorithm]
sources:
  - docs/Evidence/DISORDER-LC-001-Evidence.md
source_commit: 35bcacb8fa9a2080358233de559d2678d7600b14
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: disorder-lc-001-evidence
      evidence: "Test Unit ID: DISORDER-LC-001 ... Algorithm: Low-Complexity Region Detection in Protein Sequences (SEG algorithm; Wootton & Federhen)"
      confidence: high
      status: current
---

# Protein low-complexity region detection (SEG algorithm)

Partitioning a protein into **low-complexity** and **high-complexity** segments — the
compositionally biased runs (homopolymers, dipeptide repeats, `X/Y`-rich stretches) that confound
alignment and database search. Seqeron implements the **SEG algorithm** of Wootton & Federhen
(1993/1996), validated under test unit **DISORDER-LC-001**; the validation record is
[[disorder-lc-001-evidence]] and [[test-unit-registry]] tracks the unit. See
[[algorithm-validation-evidence]] for the artifact pattern.

This is the **first ingested unit of the protein disorder / features family** (DISORDER-LC / MORF /
PRED / PROPENSITY / REGION). SEG low-complexity detection is a *distinct algorithm* from intrinsic
disorder prediction (TOP-IDP) and MoRF prediction, so those units are expected to warrant their own
concepts when ingested — low-complexity regions overlap with but are not identical to intrinsically
disordered regions. It is also the **protein counterpart** of the genomic/DNA low-complexity handled
under [[repetitive-element-detection]] — a different alphabet and complexity measure (SEG's Shannon
entropy over 20 amino acids vs the DNA repeats/masking family).

## Complexity measure — Shannon entropy per window

Complexity is the **Shannon entropy of the residue composition** of a fixed-length window, in
**bits per residue**:

```
H = −Σ pᵢ · log₂(pᵢ)
```

where `pᵢ` is the fraction of the window occupied by residue type `i` (matches `s_Entropy` in NCBI
`blast_seg.c`, which normalizes composition counts by window length and converts to base-2 via
`NCBIMATH_LN2`). For a 20-letter amino-acid alphabet the maximum is `log₂(20) ≈ 4.322` bits/residue
(a random equiprobable sequence); a homopolymer window has `H = 0`.

## Two-stage scan and the three parameters

SEG has three user parameters (Wootton & Federhen; NCBI/GCG defaults given):

| Parameter | Symbol | NCBI/GCG default | Role |
|-----------|--------|------------------|------|
| Trigger window length | W | `kSegWindow` / `-WINdow` = **12** | minimum first-stage segment size |
| Trigger complexity | K1 | `kSegLocut` / `-LOWcut` = **2.2** bits | stage-1 low-complexity cutoff |
| Extension complexity | K2 | `kSegHicut` / `-HIGhcut` = **2.5** bits | stage-2 extension cutoff |

- **Stage 1 (trigger):** scan length-`W` windows; a window with entropy **≤ K1** triggers a
  low-complexity segment.
- **Stage 2 (extension):** each trigger is extended in both directions by merging length-`W`
  windows whose complexity is **≤ K2** (K2 ≥ K1), growing the segment into a contig.

Because `K1 ≤ K2`, the trigger is a stricter test than the extension.

## Reference oracle window entropies (W = 12)

Hand-derived `H` for canonical compositions (independently computed, matching `s_Entropy`):

| Window composition (L = 12) | Distinct residues | H (bits) | Behaviour at K1=2.2 / K2=2.5 |
|-----------------------------|-------------------|----------|------------------------------|
| `QQQQQQQQQQQQ` (homopolymer) | 1 | 0.000000 | ≤ K1 → **triggers** |
| `AAAAAALLLLLL` (two residues 6+6) | 2 | 1.000000 | ≤ K1 → **triggers** |
| `AAABBBCCCDDD` (four residues 3×4) | 4 | 2.000000 | ≤ K1 → triggers (but > strict K1=0.5 → no trigger) |
| `ACDEFGHIKLMN` (12 distinct) | 12 | 3.584963 | > K2 → **no trigger, no extension** |

A window of `W` distinct residues has entropy `log₂(W) = log₂(12) ≈ 3.585 > K2`, so maximal-complexity
windows neither trigger nor extend.

## Corner cases

- **Sequence shorter than W** → no full trigger window exists → **no segments** (empty result).
- **Homopolymer of length ≥ W** (e.g. 26×Q) → exactly one segment spanning the whole sequence
  (`H = 0 ≤ K1`).
- **Maximal-complexity sequence** (every window `W` distinct) → no segments.

## Deviations from the reference (repository extensions)

Two documented assumptions, neither altering source-defined segment boundaries on the canonical
cases:

1. **Region-type label (`"X-rich"` / `"X/Y-rich"`)** — a *presentation extension* on top of SEG:
   after a segment is found, it is labelled by its dominant residue (single most frequent residue
   when its fraction > 0.5 → `"X-rich"`, else the top two → `"X/Y-rich"`). SEG itself defines only
   *where* low-complexity segments are, not a textual composition label; only the 50% dominance
   threshold affects which label is emitted, and it does not move segment boundaries.
2. **Greedy single-residue extension** — the reference extends by merging length-`W` windows with
   complexity ≤ K2; the repository grows the contig one residue at a time while the *whole growing
   segment's* entropy stays ≤ K2. For the homopolymer / dipeptide oracle cases the boundaries are
   fixed by the trigger spans and are identical.

A `minLength` post-filter (drop segments shorter than a threshold) and case-insensitive input
(uppercasing) are additional API-contract behaviours.

## References

Wootton J.C. & Federhen S. (1993) *Computers & Chemistry* 17(2):149–163; (1996) *Methods in
Enzymology* 266:554–571. Reference implementation: NCBI C++ Toolkit `blast_seg.c` (`s_Entropy`,
defaults `kSegWindow=12` / `kSegLocut=2.2` / `kSegHicut=2.5`); program documentation: GCG/Weizmann
SEG help + `ncbi-seg` manpage. Full citations in [[disorder-lc-001-evidence]] (do not duplicate
here). A [[research-grade-limitations|research-grade]] implementation of the standard SEG method.
