---
type: concept
title: "Restriction-enzyme filtering (cutter classification + recognition-length selection)"
tags: [moldesign, restriction]
mcp_tools:
  - find_all_restriction_sites
  - find_restriction_sites
  - get_enzyme
sources:
  - docs/algorithms/MolTools/Restriction_Enzyme_Filtering.md
  - docs/Evidence/RESTR-FILTER-001-Evidence.md
source_commit: 6a7651515bbf8015d31e23697d8252f85ba10258
created: 2026-07-10
updated: 2026-07-14
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: restr-filter-001-evidence
      evidence: "Test Unit ID: RESTR-FILTER-001 ... Algorithm: Restriction Enzyme Filtering (by recognition-site length, blunt cutters, sticky cutters)"
      confidence: high
      status: current
---

# Restriction-enzyme filtering

The anchor for the **RESTR** family — selecting/subsetting a library of restriction
enzymes by properties of their cut, rather than searching a sequence for cut sites.
The first unit, `RESTR-FILTER-001`, covers three query helpers over the enzyme
library: `GetBluntCutters()`, `GetStickyCutters()`, and `GetEnzymesByCutLength(min, max)`
(plus the single-length overload `GetEnzymesByCutLength(length)`). This is a
**selection** operation in the same MolTools/reagent-design surface as
[[codon-optimization]] and the PCR primer units — it operates on enzyme metadata,
not on a target sequence. Validated per the templated
[[algorithm-validation-evidence|evidence artifact]] pattern (see
[[restr-filter-001-evidence]]); tracked in [[test-unit-registry]].

## Two axes of the filter

**1. End type — blunt vs sticky (a total, disjoint partition).**
A Type II enzyme cuts either symmetrically at the centre of its palindrome, leaving a
**blunt** end (both strands terminate in a base pair, no unpaired overhang), or at a
staggered position, leaving a **sticky** end (a stretch of unpaired nucleotides — a
5' or 3' **overhang**, also called a cohesive/cohesive end). There is no third
category: every end is blunt or an overhang, so the blunt-cutter set and the
sticky-cutter set are **complementary and disjoint**, and their union is the whole
library. A sticky end is simply a non-blunt end regardless of overhang polarity.

- Blunt cutters: SmaI (`CCC^GGG`), EcoRV (`GAT^ATC`), AluI (`AG^CT`), HaeIII (`GG^CC`).
- Sticky cutters: EcoRI (`G^AATTC`, 5' overhang), KpnI (`GGTAC^C`, 3' overhang),
  PstI (`CTGCA^G`, 3' overhang), NotI (`GCGGCCGC`, 8-bp), TaqI (`TCGA`).

**2. Recognition-site length.** `GetEnzymesByCutLength(min, max)` returns exactly the
enzymes whose recognition-string length lies in the **inclusive** interval `[min, max]`.
Type II undivided recognition sites are canonically **4–8 nt** (4-cutters AluI/HaeIII/TaqI;
6-cutters EcoRI/KpnI/PstI/SmaI/EcoRV; 8-cutter NotI), so `[4, 8]` returns the whole
undivided library.

## Key invariants and contract

- **Total partition:** `GetBluntCutters()` ∪ `GetStickyCutters()` = full library, and the
  two sets are disjoint. A blunt query excludes EcoRI/KpnI; a sticky query excludes SmaI.
- **Inclusive range bounds:** `[min, max]` includes both endpoints — this is an
  **API-shape assumption** (the conventional min/max meaning), not a biological
  parameter; the recognition-length values themselves are source-backed. The
  single-length overload agrees with the range overload at `min == max == length`.
- **Empty-interval / boundary behaviour:** `min > max` returns empty; a range disjoint
  from all lengths (e.g. `9..10`) returns none; negative/zero bounds return empty (no
  recognition site has length ≤ 0).
- **Interrupted (divided) palindromes fall outside 4–8:** SfiI recognizes the interrupted
  palindrome `GGCCNNNN^NGGCC` (a 13-nt string with a degenerate `NNNNN` spacer), so it is
  a sticky cutter but is **correctly excluded** by the `[4, 8]` length filter — the 4–8 nt
  range is for **undivided** sites only.

## Implementation surface

The primary spec (`docs/algorithms/MolTools/Restriction_Enzyme_Filtering.md`, unit
`RESTR-FILTER-001`, status *Production*) pins the four entry points on
`RestrictionAnalyzer` (`Seqeron.Genomics.MolTools`, `RestrictionAnalyzer.cs`):
`GetEnzymesByCutLength(int length)`, `GetEnzymesByCutLength(int minLength, int maxLength)`,
`GetBluntCutters()`, and `GetStickyCutters()`. All four are **total** — they never throw and
never return null; an empty/degenerate range (e.g. `min > max`, non-positive bounds) yields an
empty sequence rather than an error. There is no sequence input, so no alphabet/normalization
concerns apply.

- **Data structure:** the library is a fixed static `Dictionary<string, RestrictionEnzyme>` of
  Type II enzymes; each `RestrictionEnzyme` record stores its recognition string and per-strand
  cut positions, from which `RecognitionLength` and the record-derived `IsBluntEnd`
  (`CutPositionForward == CutPositionReverse`) are computed. Blunt/sticky uses only the
  cut-position equality — the center-vs-staggered criterion — and does **not** re-derive
  cleavage from sequence.
- **Evaluation:** filters are lazy LINQ `Where` over the dictionary values, so the return is a
  deferred `IEnumerable<RestrictionEnzyme>` (order follows dictionary insertion); callers
  materialize with `ToList()` when a snapshot is needed.
- **Complexity:** O(e) time (e = library size, single linear pass), O(1) extra space, for any
  filter. This is a metadata selection, not a search — the repository suffix tree was evaluated
  and found inapplicable (no text to search).
- **Not implemented:** filtering by overhang direction (5' vs 3') or by overhang sequence — for
  end-compatibility use `AreCompatible` / `FindCompatibleEnzymes` on
  [[restriction-digest-simulation]].

## Relation to the rest of MolTools

Filtering selects candidate enzymes; the complementary operations — finding where an
enzyme actually cuts a target and simulating a digest — are separate RESTR units.
[[restriction-digest-simulation]] covers the digest surface (sequence → fragments,
restriction map, and the `AreCompatible` end-compatibility test). End-type compatibility
(blunt ends are always mutually compatible; matching overhangs anneal) is the downstream
reason a user filters by cutter class when planning a ligation/cloning step.

## Sources

Wikipedia *Sticky and blunt ends* (blunt/overhang definitions, total partition, blunt-blunt
compatibility), Wikipedia *Restriction enzyme* (Type II 4–8 nt undivided palindromes;
center-cut→blunt vs staggered→sticky; EcoRI/SmaI/KpnI/PstI worked examples), Wikipedia
*List of restriction enzyme cutting sites* (4/6/8-bp length categories), and NEB/REBASE +
PMC for the KpnI 3'-overhang, EcoRI 5'-overhang, and SfiI interrupted-palindrome facts. See
[[restr-filter-001-evidence]] for the full source list and worked oracles.
