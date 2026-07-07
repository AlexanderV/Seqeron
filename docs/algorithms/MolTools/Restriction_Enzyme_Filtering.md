# Restriction Enzyme Filtering

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | RESTR-FILTER-001 |
| Related Projects | Seqeron.Genomics.MolTools |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Restriction enzyme filtering selects subsets of a restriction-enzyme library by two properties that matter for cloning design: the length of the recognition sequence and the type of DNA end the enzyme produces (blunt vs sticky). These are pure, specification-driven set operations over a curated enzyme table — no sequence input is involved. They let a user, for example, list all 6-cutters, or all blunt-end producers for a particular ligation strategy.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Type II restriction endonucleases recognize short, usually palindromic DNA sequences and cleave at or near them. Their undivided recognition sites are "4–8 nucleotides in length" [2] (a few enzymes, e.g. SfiI, use a longer interrupted palindrome with a degenerate spacer [6]). Cleavage either happens "at the center of both strands to yield a blunt end, or at a staggered position leaving overhangs called sticky ends" [2]. A blunt-ended molecule is one where "both strands terminate in a base pair" [1]; a sticky (cohesive) end carries "a stretch of unpaired nucleotides" — a 5' or 3' overhang [1].

### 2.2 Core Model

For an enzyme `e` with recognition sequence `r(e)` and per-strand cut positions `cf(e)` (forward) and `cr(e)` (reverse, measured from the 5' end of the recognition site):

- **Recognition length:** `len(e) = |r(e)|`.
- **End type:** `e` is **blunt** iff `cf(e) = cr(e)` (the two strands are cut at the same offset, i.e. a center/symmetric cut [2]); otherwise `e` is **sticky** (a staggered cut leaving a 5' or 3' overhang [1][2]).
- **Length filter (range):** `ByLength(min, max) = { e : min ≤ len(e) ≤ max }` (both bounds inclusive).
- **Blunt set:** `Blunt = { e : cf(e) = cr(e) }`; **Sticky set:** `Sticky = { e : cf(e) ≠ cr(e) }`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `Blunt ∩ Sticky = ∅` and `Blunt ∪ Sticky = Library` (total partition). | Every DNA end is either blunt or an overhang [1]; the predicate `cf = cr` is a Boolean split of the library. |
| INV-02 | `GetBluntCutters()` returns exactly enzymes with `cf = cr`. | Center cut → blunt [2]; implemented via `IsBluntEnd`. |
| INV-03 | `ByLength(min,max)` returns exactly enzymes with `min ≤ len ≤ max`. | Direct inclusive-range predicate; recognition lengths are 4–8 nt [2][3]. |
| INV-04 | `ByLength(L) = ByLength(L,L)` for any L. | Single-length filter is the degenerate equal-bounds range. |
| INV-05 | `min > max ⇒ ByLength(min,max) = ∅`. | No integer length satisfies an empty interval. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `length` | `int` | required | Exact recognition-sequence length to match (bp) | any int; meaningful 4–8 |
| `minLength` | `int` | required | Inclusive lower bound on recognition length (bp) | any int |
| `maxLength` | `int` | required | Inclusive upper bound on recognition length (bp) | any int |

`GetBluntCutters()` and `GetStickyCutters()` take no parameters.

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `IEnumerable<RestrictionEnzyme>` | The matching enzymes from the built-in library (deferred LINQ enumeration; order follows dictionary insertion). |

### 3.3 Preconditions and Validation

All four methods are total: they never throw and never return null. They enumerate the static built-in enzyme dictionary and apply a predicate. A range with `minLength > maxLength`, or non-positive bounds, yields an empty sequence rather than an error. There is no sequence input, so no alphabet/normalization concerns apply.

## 4. Algorithm

### 4.1 High-Level Steps

1. Enumerate the static enzyme library values.
2. Apply the predicate: recognition-length equality, inclusive-range membership, or `IsBluntEnd` / `!IsBluntEnd`.
3. Yield the matching enzymes (lazy LINQ `Where`).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The library is a fixed `Dictionary<string, RestrictionEnzyme>` of Type II enzymes; each `RestrictionEnzyme` stores its recognition sequence and per-strand cut positions, from which `RecognitionLength` and `IsBluntEnd` (`cf == cr`) are derived. The blunt/sticky predicate uses only the cut-position equality, matching the center-vs-staggered criterion [2]. End-type classification of representative enzymes (SmaI/EcoRV/AluI/HaeIII = blunt; EcoRI = 5' overhang; KpnI/PstI = 3' overhang) is corroborated against [2], [4], [5].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Any filter | O(e) | O(1) extra | `e` = library size; single linear pass, deferred enumeration. No search over sequences, so the repository suffix tree is not applicable. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RestrictionAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs)

- `RestrictionAnalyzer.GetEnzymesByCutLength(int length)`: enzymes whose recognition length equals `length`.
- `RestrictionAnalyzer.GetEnzymesByCutLength(int minLength, int maxLength)`: enzymes whose recognition length is in the inclusive range.
- `RestrictionAnalyzer.GetBluntCutters()`: enzymes with `IsBluntEnd == true`.
- `RestrictionAnalyzer.GetStickyCutters()`: enzymes with `IsBluntEnd == false`.

### 5.2 Current Behavior

Filters operate over the static built-in library only (not over any user sequence). `IsBluntEnd` is a record-derived property: `CutPositionForward == CutPositionReverse`. The methods return deferred `IEnumerable<>`; callers materialize with `ToList()` when a snapshot is needed. The unit is not a search/matching operation, so the repository suffix tree was evaluated and found inapplicable (no text to search).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Blunt vs sticky split by center (`cf = cr`) vs staggered cut [2]; blunt = "both strands terminate in a base pair", sticky = an overhang [1].
- Inclusive recognition-length range filter over Type II 4–8 nt sites [2][3].
- Total partition of the library into blunt and sticky (INV-01) [1].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Filtering by overhang direction (5' vs 3') or overhang sequence; for end-compatibility use `RestrictionAnalyzer.AreCompatible` / `FindCompatibleEnzymes`.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Range bounds inclusive | Assumption | Determines whether boundary lengths are returned | accepted | API-shape convention; recognition-length values are source-backed [2][3]. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `ByLength(min,max)` with `min > max` | empty | Empty interval (INV-05) |
| `ByLength(4,8)` | full library except SfiI | Undivided sites are 4–8 nt [2]; SfiI is a 13-nt interrupted palindrome `GGCCNNNN^NGGCC` [6], excluded |
| `ByLength(9,10)` | empty | Above the 8-nt maximum for undivided sites [2] |
| Non-positive bounds | empty | No recognition site has length ≤ 0 |

### 6.2 Limitations

The library is a curated subset of common Type II enzymes, not the full REBASE catalog; filters reflect only the enzymes present. Blunt/sticky classification depends on the stored cut positions being correct; it does not re-derive cleavage from sequence.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// All blunt-end producers in the library
var blunt = RestrictionAnalyzer.GetBluntCutters().Select(e => e.Name).ToList();

// All 6-cutters
var sixCutters = RestrictionAnalyzer.GetEnzymesByCutLength(6, 6).ToList();

// Every undivided 4-8 nt site (excludes the 13-nt interrupted palindrome SfiI)
var undivided = RestrictionAnalyzer.GetEnzymesByCutLength(4, 8).ToList();
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RestrictionAnalyzer_Filter_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/RestrictionAnalyzer_Filter_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [RESTR-FILTER-001-Evidence.md](../../../docs/Evidence/RESTR-FILTER-001-Evidence.md)
- Related algorithms: [Restriction_Digest_Simulation](../MolTools/Restriction_Digest_Simulation.md)

## 8. References

1. Wikipedia. 2026. Sticky and blunt ends. https://en.wikipedia.org/wiki/Sticky_and_blunt_ends
2. Wikipedia. 2026. Restriction enzyme. https://en.wikipedia.org/wiki/Restriction_enzyme
3. Wikipedia. 2026. List of restriction enzyme cutting sites. https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites
4. New England Biolabs / REBASE. KpnI (R0142). https://www.neb.com/en/products/r0142-kpni
5. New England Biolabs / REBASE. EcoRI-HF (R3101). https://www.neb.com/en/products/r3101-ecori-hf
6. SfiI interrupted palindrome `5'-GGCCNNNN^NGGCC-3'` (homology model article). PMC. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC548270/
