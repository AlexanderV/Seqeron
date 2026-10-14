# restriction_digest

Simulate a restriction digest of a linear or circular molecule.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `restriction_digest` |
| **Method ID** | `RestrictionAnalyzer.Digest` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Cuts the molecule at the top-strand (forward-strand) cut positions of the named enzymes and returns the resulting fragments in ascending `startPosition`. Each fragment carries its sequence, start position, length, flanking enzymes (`leftEnzyme`/`rightEnzyme`), fragment number, and the chemistry of each end created by a cut (`leftOverhangType`/`leftOverhang`, `rightOverhangType`/`rightOverhang`).

**Topology.** With `topology=Linear` (the default) `k` distinct cut positions yield `k+1` fragments; with zero cuts the whole molecule is returned as one fragment, and its two ends are the original termini (overhang fields `null`). With `topology=Circular` the site search wraps past the origin — so an origin-spanning recognition site is found and cut — and `k` cut positions yield exactly `k` fragments. The origin-spanning fragment is reported **last** (its `startPosition` is the largest cut), and every circular fragment end is a cut, so both flanks carry an enzyme and an overhang. A circular molecule with zero cuts stays whole (one full-length uncut fragment, no ends); a single cutter linearizes it into one full-length fragment.

**Fragment sequences are top-strand bases.** A digest partitions the top strand at the top-strand cuts, the same convention as Biopython `Bio.Restriction` `enzyme.catalyse(seq[, linear=False])`, so a sticky end's bases belong to the fragment downstream of the top-strand cut and the single-stranded extension is reported separately in the overhang fields (read 5′→3′ on the strand that carries it, taken from the molecule). Blunt cutters report `Blunt` with an empty overhang. Biopython emits the origin-spanning fragment first rather than last, so the circular fragment order differs from Biopython's by a rotation of one.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L761](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L761)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to digest (non-empty). |
| `enzyme_names` | string[] | Yes | Enzyme names to use (≥1). |
| `topology` | string | No | `Linear` (default) or `Circular`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `fragments` | array | Fragments, each `{sequence, startPosition, length, leftEnzyme, rightEnzyme, fragmentNumber, leftOverhangType, leftOverhang, rightOverhangType, rightOverhang}`. |

`leftOverhangType`/`rightOverhangType` are `Blunt`, `FivePrime` or `ThreePrime`, or `null` when that end was not created by a cut. `leftOverhang`/`rightOverhang` give the extension's bases (empty string for a blunt end, `null` for an uncut end).

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | At least one enzyme name is required |

## Examples

### Example 1: EcoRI on linear `AAAGAATTCAAA`

Cut at position 4 → `AAAG` (0–4) and `AATTCAAA` (4–12); the cut ends are 5′ overhangs `AATT`.

### Example 2: EcoRI on circular `TTCAAAAAAAGAA`

The `GAATTC` site spans the origin (starting at position 10), so the plasmid is linearized into one 13 bp fragment `AATTCAAAAAAAG` starting at the cut, position 11 — matching `Bio.Restriction.EcoRI.catalyse(Seq("TTCAAAAAAAGAA"), linear=False) == (Seq("AATTCAAAAAAAG"),)`. With `topology=Linear` the same sequence is uncut.

### Example 3: No cut site → single whole-sequence fragment (linear or circular).

## See Also

- [digest_summary](digest_summary.md), [find_restriction_sites](find_restriction_sites.md), [restriction_map](restriction_map.md)
