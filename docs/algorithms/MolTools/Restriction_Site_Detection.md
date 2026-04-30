# Restriction Site Detection

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | N/A |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Restriction site detection identifies DNA sequences recognized by restriction endonucleases and computes the corresponding cut positions. In this repository, the implementation scans both strands for recognition-sequence matches, supports IUPAC ambiguity codes in enzyme patterns, and returns site records that include the enzyme, strand, cut position, and matched sequence. The current implementation is sequence-based and centered on a built-in enzyme catalog.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Restriction enzymes cleave DNA at specific recognition sequences, and Type II enzymes are the standard molecular-biology workhorses because they cut at or near their recognition sites. Recognition sequences are commonly 4-8 bases long, and many are palindromic so the same site appears on both strands. The original document also classifies cleavage products into 5' overhangs, 3' overhangs, and blunt ends. Sources: Wikipedia (Restriction enzyme, Restriction site, EcoRI), Roberts (1976), IUPAC-IUB (1970).

### 2.2 Core Model

For each enzyme, the algorithm slides a window of the recognition-sequence length across the input sequence and checks character-by-character matches under IUPAC ambiguity rules. A forward-strand match at position `i` yields a cut position of:

$$
cutPosition = i + CutPositionForward
$$

Reverse-strand matches are found by scanning the reverse complement, converting the recognition-site start back to forward coordinates, and using:

$$
forwardPos = sequenceLength - i - patternLength
$$

with a reverse-strand cut position of `forwardPos + CutPositionReverse`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 <= Position <= sequence.Length - recognitionLength` for every reported site | Sites are yielded only from valid window starts |
| INV-02 | `RecognizedSequence.Length == Enzyme.RecognitionSequence.Length` | The source slices exactly the pattern length |
| INV-03 | Enzyme-name lookup is case-insensitive | The built-in enzyme dictionary uses `StringComparer.OrdinalIgnoreCase` |
| INV-04 | Empty raw-string input yields no sites | The raw-string overload short-circuits on null or empty input |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA sequence to scan | Null `DnaSequence` input throws; empty raw-string input yields no results |
| `enzymeName` | `string` | required | Name of the restriction enzyme to search | Unknown names throw `ArgumentException` |
| `enzymeNames` | `string[]` | required for multi-enzyme overload | Set of enzyme names to scan | Multi-enzyme search yields the union of per-enzyme results |
| `enzyme` | `RestrictionEnzyme` | required for custom-enzyme overload | Enzyme definition to apply directly | Null input throws `ArgumentNullException` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | Start position of the recognition sequence |
| `Enzyme` | `RestrictionEnzyme` | Enzyme definition that matched |
| `IsForwardStrand` | `bool` | `true` for forward-strand matches and `false` for reverse-strand matches |
| `CutPosition` | `int` | Computed cut position |
| `RecognizedSequence` | `string` | Actual sequence segment that matched the enzyme pattern |

### 3.3 Preconditions and Validation

The `DnaSequence` overloads require non-null sequence input. The raw-string overload returns no results for null or empty strings. Unknown enzyme names raise `ArgumentException`, while enzyme lookup itself is case-insensitive. Pattern matching uses IUPAC ambiguity codes through `IupacHelper.MatchesIupac(...)`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Resolve the enzyme definition from the built-in catalog or use the supplied custom enzyme.
2. Scan the forward strand for IUPAC-aware recognition-sequence matches.
3. Yield a forward-strand `RestrictionSite` for each match.
4. Reverse-complement the sequence and repeat the recognition search.
5. Convert reverse-strand positions back to forward coordinates and yield reverse-strand `RestrictionSite` records.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Recognition-sequence categories called out in the original document:

| Category | Examples | Recognition Length |
|----------|----------|-------------------|
| 4-cutters | AluI, MspI, TaqI | 4 bp |
| 6-cutters | EcoRI, BamHI, HindIII | 6 bp |
| 8-cutters | NotI, PacI, AscI | 8 bp |

Overhang types:

| Type | Description | Example |
|------|-------------|---------|
| 5' overhang | `CutPositionForward < CutPositionReverse` | EcoRI |
| 3' overhang | `CutPositionForward > CutPositionReverse` | PstI |
| Blunt end | `CutPositionForward == CutPositionReverse` | EcoRV |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Single-enzyme `FindSites` | `O(n × m)` | `O(1)` streaming | `n` is sequence length and `m` is recognition-sequence length |
| Multi-enzyme `FindSites` / `FindAllSites` | `O(n × k × m)` | `O(1)` streaming | `k` is the number of enzymes scanned |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RestrictionAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs), [IupacHelper.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/IupacHelper.cs)

- `RestrictionAnalyzer.FindSites(DnaSequence, string)`: Searches a validated DNA sequence for one enzyme.
- `RestrictionAnalyzer.FindSites(DnaSequence, params string[])`: Aggregates per-enzyme searches.
- `RestrictionAnalyzer.FindSites(string, string)`: Raw-string overload that uppercases input and yields no results for empty input.
- `RestrictionAnalyzer.FindAllSites(DnaSequence)`: Scans the full built-in enzyme set.
- `RestrictionAnalyzer.GetEnzyme(string)`: Case-insensitive enzyme lookup.

### 5.2 Current Behavior

The implementation ships with an in-memory catalog of common enzymes, including 4-cutters, 6-cutters, 8-cutters, blunt cutters, and enzymes with degenerate recognition motifs such as `SfiI` and `HincII`. Both strands are reported, so palindromic sites can appear twice at the same genomic position with different `IsForwardStrand` values. Reverse-strand matches compute `Position` in forward coordinates and compute the cut site with `CutPositionReverse`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Recognition-sequence scanning for a catalog of Type II-style enzymes.
- IUPAC ambiguity handling in enzyme patterns.
- Strand-aware reporting of recognition and cut positions.

**Intentionally simplified:**

- The implementation relies on a built-in enzyme catalog; **consequence:** only enzymes present in that catalog are available through name-based lookup.
- Site finding is purely sequence-based; **consequence:** methylation sensitivity and experimental context are not modeled.

**Not implemented:**

- Contextual enzyme activity effects such as methylation-state filtering; **users should rely on:** external enzymology references or downstream validation when those matter.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty raw-string input | Returns no sites | Explicit short-circuit in source |
| Unknown enzyme name | Throws `ArgumentException` | Name lookup must succeed before scanning |
| Palindromic recognition site | May appear on both strands at the same position | Both strands are searched independently |
| Degenerate recognition sequence | Matched with IUPAC ambiguity rules | `IupacHelper.MatchesIupac(...)` handles the pattern |

### 6.2 Limitations

The current implementation is limited to the built-in enzyme definitions and sequence-only matching. It does not model methylation dependence, assay conditions, or experimental digestion efficiency, and double reporting of palindromic sites must be interpreted downstream when unique cut counts are required.

## 8. References

1. Wikipedia: Restriction enzyme - https://en.wikipedia.org/wiki/Restriction_enzyme
2. Wikipedia: Restriction site - https://en.wikipedia.org/wiki/Restriction_site
3. Wikipedia: EcoRI - https://en.wikipedia.org/wiki/EcoRI
4. Roberts RJ (1976) - Restriction endonucleases, CRC Critical Reviews in Biochemistry.
5. REBASE (Restriction Enzyme Database) - http://rebase.neb.com/
6. IUPAC-IUB (1970) - Nucleic acid notation.
