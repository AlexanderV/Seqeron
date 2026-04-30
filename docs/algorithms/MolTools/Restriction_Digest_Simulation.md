# Restriction Digest Simulation

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | N/A |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Restriction digest simulation models DNA fragmentation after cleavage by one or more restriction enzymes. In this repository, digestion is performed by locating cut sites, taking forward-strand cut positions to avoid double-counting palindromic sites, and emitting fragments bounded by those cuts. The same source surface also provides digest summaries, restriction maps, and overhang-compatibility checks.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A restriction digest cleaves DNA at enzyme recognition sequences and produces fragments whose sizes can be checked experimentally by gel electrophoresis. The original document emphasizes the fragment-sum principle: the sum of fragment lengths must equal the original sequence length. It also states that two enzymes produce compatible ligatable ends when both are blunt cutters or when they generate identical sticky-end overhangs. Sources: Wikipedia (Restriction digest, Restriction enzyme, Restriction map), Addgene Restriction Digest Protocol, REBASE.

### 2.2 Core Model

Given a sequence of length `L` and forward-strand cut positions `c1 < c2 < ... < ck`, the digest fragments are the half-open intervals:

$$
[0, c_1), [c_1, c_2), \ldots, [c_k, L)
$$

This yields `k + 1` fragments when `k` cuts are present. `GetDigestSummary(...)` sorts fragment sizes in descending order, and `CreateMap(...)` groups sites by enzyme name while treating unique cutters as enzymes with exactly one forward-strand site.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `k` forward-strand cut positions produce `k + 1` fragments | `Digest(...)` inserts boundaries `[0, cuts..., sequence.Length]` |
| INV-02 | The sum of fragment lengths equals the original sequence length | Adjacent boundaries partition the sequence |
| INV-03 | `FragmentSizes` in `DigestSummary` are sorted descending | `GetDigestSummary(...)` orders sizes descending before constructing the record |
| INV-04 | `AreCompatible(A, B) == AreCompatible(B, A)` | Compatibility is determined from bluntness or identical overhang type and sequence |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` | required | DNA sequence to digest or map | Null input throws `ArgumentNullException` |
| `enzymeNames` | `string[]` | required for digestion | Enzymes used for site discovery and fragment generation | `Digest(...)` throws `ArgumentException` when none are supplied |
| `enzyme1Name`, `enzyme2Name` | `string` | required | Enzyme names for compatibility analysis | Unknown names produce `false` in `AreCompatible(...)` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `DigestFragment` | record | Sequence, start position, length, one representative left/right enzyme name per fragment boundary, and fragment number |
| `DigestSummary` | record | Total fragment count, sorted sizes, largest/smallest fragment, average size, and enzyme list |
| `RestrictionMap` | record | Sequence length, sites, grouped positions, forward-strand total site count, unique cutters, and non-cutters |
| `compatible` | `bool` | Whether two enzymes generate compatible ends |

### 3.3 Preconditions and Validation

`Digest(...)` requires a non-null sequence and at least one enzyme name. When no cuts are found, it returns a single fragment equal to the full sequence. `CreateMap(...)` accepts zero enzyme names and then scans the full built-in enzyme catalog. `AreCompatible(...)` returns `false` if either enzyme name is not present in the catalog.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate the input sequence and enzyme list.
2. Find restriction sites for each requested enzyme.
3. Keep only forward-strand cut positions to avoid double-counting palindromic sites.
4. Sort and deduplicate cut positions.
5. Emit fragments between adjacent boundaries.
6. Build summaries or maps by sorting fragment sizes and grouping site positions by enzyme.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Compatibility rules preserved from the original document and confirmed in source:

| Rule | Result |
|------|--------|
| Both enzymes produce blunt ends | Compatible |
| Overhang types differ | Not compatible |
| Overhang types match and overhang sequences match | Compatible |
| Overhang sequences differ | Not compatible |

Examples called out in the original document:

| Pair | Outcome |
|------|---------|
| BamHI and BglII | Compatible (`GATC` overhang) |
| EcoRV and SmaI | Compatible (both blunt) |
| EcoRI and PstI | Not compatible |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Digest` | `O(n + k log k)` | `O(k)` | Site discovery plus sorted cut handling |
| `GetDigestSummary` | `O(n + k log k)` | `O(k)` | Materializes digest fragments and sorts sizes |
| `CreateMap` | `O(n + k log k)` | `O(k)` | Groups positions by enzyme and identifies unique/non-cutters |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RestrictionAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs)

- `RestrictionAnalyzer.Digest(DnaSequence, params string[])`: Generates digest fragments.
- `RestrictionAnalyzer.GetDigestSummary(DnaSequence, params string[])`: Produces sorted fragment statistics.
- `RestrictionAnalyzer.CreateMap(DnaSequence, params string[])`: Builds a restriction map from requested or all enzymes.
- `RestrictionAnalyzer.AreCompatible(string, string)`: Tests end compatibility.
- `RestrictionAnalyzer.FindCompatibleEnzymes()`: Enumerates compatible enzyme pairs from the built-in catalog.

### 5.2 Current Behavior

The current implementation records forward-strand cut positions only when digesting so that palindromic sites are not double-counted. If no cut positions are found, the digest yields a single fragment equal to the original sequence. When multiple enzymes cut at the same coordinate, the cut position is preserved but the digest stores only one representative enzyme name for that fragment boundary. `GetDigestSummary(...)` sorts fragment sizes descending before constructing the record. `CreateMap(...)` counts `TotalSites` on the forward strand only and defines `UniqueCutters` from distinct forward-strand positions, while `NonCutters` are requested enzymes that do not appear in the grouped map.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Fragment generation from ordered restriction cut positions.
- Restriction-map construction by grouping sites per enzyme.
- Compatibility rules for blunt ends and identical sticky-end overhangs.

**Intentionally simplified:**

- Digest simulation uses forward-strand cut positions only; **consequence:** strand-paired palindromic site reports are collapsed to a single cut for digestion purposes.
- The model reports virtual fragments and statistics only; **consequence:** gel-migration behavior is not simulated.

**Not implemented:**

- Gel electrophoresis simulation and other downstream laboratory effects; **users should rely on:** external analysis or visualization tools for those steps.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| No cut sites found | Returns one fragment equal to the original sequence | Explicit special case in `Digest(...)` |
| First fragment | `LeftEnzyme = null` | No enzyme cuts before position 0 |
| Last fragment | `RightEnzyme = null` | No enzyme cuts after the final boundary |
| Zero-length fragments | Not generated | The source yields only when `length > 0` |

### 6.2 Limitations

The digest model is a sequence-partitioning simulation. It does not simulate gel electrophoresis, incomplete digestion, methylation effects, or circular-DNA behavior, and it relies on the enzyme definitions available in the built-in restriction catalog. When different enzymes cut at the same coordinate, fragment-boundary provenance is reduced to a single representative enzyme name.

## 8. References

1. Wikipedia: Restriction digest - https://en.wikipedia.org/wiki/Restriction_digest
2. Wikipedia: Restriction enzyme - https://en.wikipedia.org/wiki/Restriction_enzyme
3. Wikipedia: Restriction map - https://en.wikipedia.org/wiki/Restriction_map
4. Addgene: Restriction Digest Protocol - https://www.addgene.org/protocols/restriction-digest/
5. Roberts RJ (1976) - Restriction endonucleases, CRC Critical Reviews in Biochemistry.
6. REBASE (Restriction Enzyme Database) - http://rebase.neb.com/
