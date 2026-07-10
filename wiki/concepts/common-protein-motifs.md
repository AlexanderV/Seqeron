---
type: concept
title: "Common protein motifs (built-in PROSITE-pattern dictionary scan)"
tags: [analysis, algorithm, protein, motif]
sources:
  - docs/Evidence/PROTMOTIF-COMMON-001-Evidence.md
source_commit: 12b13e4ecc31636e0c27a2c4b0098bf11d6cc054
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: protmotif-common-001-evidence
      evidence: "Test Unit ID: PROTMOTIF-COMMON-001 ... Algorithm: Common Motif Finding (ProteinMotifFinder.FindCommonMotifs)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:coiled-coil-prediction
      source: protmotif-common-001-evidence
      evidence: "Both are ProteinMotif-family units on ProteinMotifFinder sharing the 0-based-inclusive MotifMatch.Start/End shape (Evidence cites sibling units PROTMOTIF-FIND-001 / PROTMOTIF-PATTERN-001); coiled-coil scores a/d hydrophobic periodicity while common-motif matches a fixed PROSITE-pattern dictionary — distinct algorithms in one family"
      confidence: medium
      status: current
---

# Common protein motifs (built-in PROSITE-pattern dictionary scan)

**Common motif finding** scans an amino-acid sequence against a small **built-in
dictionary of canonical, well-known PROSITE patterns** and reports every match of every
entry. Seqeron exposes it as `ProteinMotifFinder.FindCommonMotifs`: it whole-dictionary
scans the fixed `CommonMotifs` catalog (each entry a PROSITE regular expression) and
aggregates all hits, each carrying the matching entry's identity. Validated under test unit
**PROTMOTIF-COMMON-001**; the validation record is [[protmotif-common-001-evidence]] and
[[test-unit-registry]] tracks the unit. See [[algorithm-validation-evidence]] for the
artifact pattern.

This is a **protein** motif finder — distinct from the DNA-side family. It is the
protein analogue of the fixed-catalog DNA scanner [[regulatory-element-detection]] (which
matches a cited catalog of regulatory consensus strings), and unlike the caller-supplies-the-set
exact matcher [[known-motif-search]] (DNA), the query set here is a **built-in, curated
protein-motif dictionary**. Within the ProteinMotif family it is a sibling of the windowed
[[coiled-coil-prediction]] and of
[[protein-domain-and-signal-peptide-prediction]] (which matches longer PROSITE **domain**
patterns and adds an opt-in Plan7 profile-HMM engine for profile-only families): all live on
`ProteinMotifFinder`; this unit scans a fixed short-motif dictionary while domain prediction
matches domain-length signatures, and coiled-coil scores a/d hydrophobic periodicity.

## PROSITE pattern syntax

Each catalog entry is a PROSITE pattern; matching follows the ScanProsite rules:

| Element | Meaning |
|---------|---------|
| `A` (single letter) | required IUPAC amino acid at that position |
| `x` | any amino acid (wildcard) |
| `[ST]` | ambiguity — any one of the bracketed residues (Ser **or** Thr) |
| `{P}` | negated set — any residue **except** those braced (any but Proline) |
| `x(2)` | fixed repetition — `x-x` |
| `x(2,4)` | variable repetition — `x-x`, `x-x-x`, or `x-x-x-x` |
| `-` | element separator |

## The built-in `CommonMotifs` catalog

The dictionary is a fixed set of well-known short functional motifs. The regex column is the
translation used by the implementation:

| Accession | Name | PROSITE pattern | Regex form | Meaning |
|-----------|------|-----------------|------------|---------|
| PS00001 | `ASN_GLYCOSYLATION` | `N-{P}-[ST]-{P}` | `N[^P][ST][^P]` | N-glycosylation site |
| PS00005 | `PKC_PHOSPHO_SITE` | `[ST]-x-[RK]` | `[ST].[RK]` | Protein kinase C phosphorylation site |
| PS00006 | `CK2_PHOSPHO_SITE` | `[ST]-x(2)-[DE]` | `[ST].{2}[DE]` | Casein kinase II phosphorylation site |
| PS00016 | `RGD` | `R-G-D` | `RGD` | Cell attachment (RGD) sequence |
| PS00017 | `ATP_GTP_A` | `[AG]-x(4)-G-K-[ST]` | `[AG].{4}GK[ST]` | ATP/GTP-binding P-loop (motif A) |

## API contract and invariants

| Aspect | Behaviour |
|--------|-----------|
| Coordinates | **0-based, inclusive** `MotifMatch.Start`/`End` (repository convention; PROSITE/ScanProsite reports **1-based** — see below) |
| Substring invariant | every match's `Sequence == protein.Substring(Start, End - Start + 1)` |
| Identity | each match's `MotifName` / `Pattern` (accession) come from the matching `CommonMotifs` entry — aggregation carries the entry's identity |
| Aggregation | one whole-dictionary scan aggregates hits of **all** entries from a single sequence |
| Overlaps | overlapping occurrences of a pattern are **all reported** (ScanProsite default "overlaps, no includes"); only a match fully contained within another is suppressed |
| null / empty input | returns empty |

## Worked oracles

From synthetic windows constructed to satisfy/violate each pattern (0-based inclusive):

- `AAAANFTAAAA` / PS00001 → `4..7` = `NFTA`.
- `AAAANPSAAAAANPTAAA` / PS00001 → **no match** (Proline at the excluded `{P}` position rejects the site).
- `AAAAASARKAAA` / PS00005 → `5..7` = `SAR`.
- `AAAASAAEASDEDAAA` / PS00006 → `4..7` = `SAAE` and `9..12` = `SDED` (two hits).
- `AAAAAGXXXXGKSAAAA` / PS00017 → `5..12` = `GXXXXGKS`.
- `AARGDKK` / PS00016 → `2..4` = `RGD`.
- `RGDRGD` / PS00016 (overlap test) → `0..2` = `RGD` and `3..5` = `RGD` (two non-overlapping hits).

## Deviations and assumptions

Only **one** assumption, an API-shape convention with no correctness effect: PROSITE/ScanProsite
report **1-based** inclusive coordinates, whereas `MotifMatch` records **0-based** `Start`/`End`
(matching the sibling units PROTMOTIF-FIND-001 and PROTMOTIF-PATTERN-001). The matched substring
content and relative positions are identical; only the coordinate origin differs. No algorithm
deviations. The catalog is a small curated subset of PROSITE, not the full database — the general
engine that scans an **arbitrary caller-supplied** PROSITE pattern (with PROSITE→regex
conversion, overlapping-match lookahead and information-content scoring) is
[[protein-motif-pattern-search]] (`FindMotifByPattern`, unit PROTMOTIF-FIND-001), of which
this fixed-dictionary scan is one application.

## References

ExPASy PROSITE entries PS00001 / PS00005 / PS00006 / PS00016 / PS00017 and the ScanProsite
pattern-syntax documentation (accessed 2026-06-14); Sigrist et al. 2013, *Nucleic Acids Research*
41(D1):D344–D347. Full citations in [[protmotif-common-001-evidence]] (do not duplicate here).
