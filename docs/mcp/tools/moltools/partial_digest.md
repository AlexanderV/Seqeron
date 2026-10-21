# partial_digest

Simulate a partial restriction digest: every fragment that can arise when any subset of the cut sites is left uncut.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `partial_digest` |
| **Method ID** | `RestrictionAnalyzer.PartialDigest` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

A partial (incomplete) digest leaves some recognition sites uncut, so besides the complete-digest fragments the gel also shows every fragment that spans one or more uncut sites. This tool lists all of them. It needs no data beyond the complete digest: it uses exactly the cut positions, the per-cut enzyme lists and the end overhangs that [restriction_digest](restriction_digest.md) reports, and every fragment carries the same members (`sequence`, `startPosition`, `length`, `leftEnzyme(s)`/`rightEnzyme(s)`, overhang type and sequence of each end). `missedCuts` is the number of cut positions strictly inside the fragment.

**Definition.** For a linear molecule of length `L` with cut positions `c_1 < … < c_k`, the fragments are all intervals `[a, b)`, `a < b`, over `X = {0, c_1, …, c_k, L}` — the difference multiset ΔX of the partial digest problem, `(k+1)(k+2)/2` fragments (Rosalind PDPL "Creating a Restriction Map", after Jones & Pevzner). For a circular molecule with `k ≥ 1` cuts, every arc from one cut to another (wrapping past the origin; `k−1` missed cuts is the molecule linearised at that cut) plus the intact circle, reported with `missedCuts = k`, `startPosition 0` and no ends — `k² + 1` fragments. With no site the single uncut molecule is returned (`missedCuts = 0`), as in `restriction_digest`.

**Order.** By `missedCuts`, then `startPosition`; `fragmentNumber` is the 1-based position in that order. The `missedCuts = 0` fragments come first and are exactly the `restriction_digest` output for the same input.

**Methylation (opt-in).** `methylation` (`none`|`dam`|`dcm`|`dam+dcm`|`emboss`) masks the bases the chosen methylase set methylates and drops every enzyme window containing one before the cut list is built — the semantics of EMBOSS `restrict -methylation -mfile Emethylsites.dat`, as in [restriction_digest](restriction_digest.md) (RESTR-FIND-001 Evidence §7). A blocked site is neither a cut nor a missed cut: `k` counts the surviving cuts only, so the fragment set is the constructed set above over the EMBOSS-surviving cuts (B08 F83). Default `none` leaves the output unchanged. Example: `GAATTCAAAAGATCAAAAGAATTCAAACCCGGGAAAGATCGAAATCGAAAA` with `EcoRI, MboI, SmaI, TaqI` has 7 cuts and 36 fragments; with `methylation=dam` Dam blocks MboI (10, 36), leaving cuts `{1, 19, 30, 39, 45}` and the 21 fragments of `X = {0, 1, 19, 30, 39, 45, 51}` (lengths `1, 6, 6, 9, 11, 12, 15, 18, 19, 20, 21, 26, 29, 30, 32, 38, 39, 44, 45, 50, 51`).

**Limits.** `max_missed_cuts` (optional, ≥ 0) keeps only fragments with at most that many uncut sites (`0` = complete digest). A request that would produce more than 100,000 fragments is refused; lower `max_missed_cuts`. Only which fragments **can** arise is computed — the yield or relative band intensity of each partial product (which depends on enzyme amount, time and per-site rates) is not modelled.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.PartialDigest.cs#L112](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.PartialDigest.cs#L112)
- Algorithm: [Restriction_Digest_Simulation.md](../../../algorithms/MolTools/Restriction_Digest_Simulation.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to digest (non-empty). |
| `enzyme_names` | string[] | Yes | Enzyme names to use (≥1). |
| `topology` | string | No | `Linear` (default) or `Circular`. |
| `max_missed_cuts` | integer \| null | No | Upper bound on uncut sites per fragment; omit for no bound. |
| `methylation` | string | No | `none` (default), `dam` (`GmATC`), `dcm` (`mCmCAGG`, `mCmCTGG`), `dam+dcm` or `emboss` (all five `Emethylsites.dat` entries incl. EcoKI). Case-insensitive; an unknown value is an error. A recognition window containing a methylated base is not a site (EMBOSS `restrict -methylation`; a window rule, not an enzyme-sensitivity table — BamHI is blocked by Dam). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `fragments` | array | Each `{fragment, missedCuts}`; `fragment` has the `restriction_digest` fragment fields. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | At least one enzyme name is required |
| 1003 | The missed-cut limit cannot be negative. |
| 1004 | A partial digest with {k} cut positions would produce {n} fragments, above the limit of 100000; pass a smaller maxMissedCuts. |

## Examples

### Example 1: Rosalind's sample, linear, `max_missed_cuts = 1`

`AGCTGATCGA` with AluI, Sau3AI and TaqI: Biopython 1.88 `search` gives AluI `[3]`, Sau3AI `[5]`, TaqI `[8]` (1-based), i.e. cuts `{2, 4, 7}` on a 10 bp molecule — Rosalind's `X = {0, 2, 4, 7, 10}`. Without a limit the ten fragment lengths are Rosalind's sample `ΔX = {2, 2, 3, 3, 4, 5, 6, 7, 8, 10}`; with `max_missed_cuts = 1` the result is the seven fragments `AG`, `CT`, `GAT`, `CGA` (complete digest) and `AGCT`, `CTGAT`, `GATCGA` (one site uncut).

### Example 2: circular, AluI + TaqI

Two cuts (2, 7) give `k² + 1 = 5` fragments: `CTGAT`, `CGAAG` (complete digest, the second wraps past the origin), `CTGATCGAAG`, `CGAAGCTGAT` (linearised at the AluI / TaqI cut), and the intact circle `AGCTGATCGA` with `missedCuts = 2`.

### Example 3: a negative `max_missed_cuts` is rejected.

## See Also

- [restriction_digest](restriction_digest.md), [digest_summary](digest_summary.md), [restriction_map](restriction_map.md)
