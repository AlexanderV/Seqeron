# Restriction Digest Simulation

## Overview

Restriction digest simulation is the computational modeling of DNA fragmentation by restriction endonucleases. Given a DNA sequence and one or more restriction enzymes, the algorithm identifies cut sites, calculates fragment positions and sizes, and produces a virtual representation of the resulting fragments—equivalent to what would be observed on a gel electrophoresis.

## Algorithm Description

### Core Algorithm: `Digest`

The `Digest` method simulates restriction enzyme digestion:

1. **Validate Input**
   - Check sequence is not null
   - Verify at least one enzyme is specified (throw `ArgumentException` if none)

2. **Find Cut Positions**
   - For each enzyme, call `FindSites` to locate recognition sites
   - Filter to forward-strand sites only (to avoid double-counting palindromic sites)
   - Extract cut positions and store in a sorted set
   - Maintain mapping from cut position to enzyme

3. **Handle No Cuts Case**
   - If no cut sites found, return entire sequence as single fragment

4. **Generate Fragments**
   - Create ordered list: [0, cut1, cut2, ..., cutN, sequenceLength]
   - For each adjacent pair of positions, create a fragment:
     - Extract subsequence between positions
     - Record start position, length, flanking enzyme names
     - Assign sequential fragment numbers

5. **Return Results**
   - Yield `DigestFragment` records

### Complexity

- **Time**: O(n + k log k) where n = sequence length, k = number of cut sites
  - O(n) for finding sites
  - O(k log k) for sorting cut positions
- **Space**: O(k) for storing cut positions

### GetDigestSummary Algorithm

Produces a statistical summary of the digest:

```csharp
var fragments = Digest(sequence, enzymeNames).ToList();
var sizes = fragments.Select(f => f.Length).OrderByDescending(x => x).ToList();
return new DigestSummary(
    TotalFragments: fragments.Count,
    FragmentSizes: sizes,  // Sorted descending
    LargestFragment: sizes.First(),
    SmallestFragment: sizes.Last(),
    AverageFragmentSize: sizes.Average(),
    EnzymesUsed: enzymeNames);
```

### CreateMap Algorithm

Creates a comprehensive restriction map showing all enzyme sites:

1. Find all sites for specified enzymes (or all enzymes if none specified)
2. Group sites by enzyme name
3. Identify unique cutters (exactly one forward-strand site)
4. Identify non-cutters (enzymes with zero sites)
5. Return `RestrictionMap` with all data

## Scientific Background

### Restriction Digest

A restriction digest is a laboratory procedure where DNA is cleaved by restriction enzymes at specific recognition sequences. The resulting fragments can be separated by gel electrophoresis to verify size and count.

**Key validation principle**: The sum of all fragment sizes must equal the original sequence length.

Source: Wikipedia (Restriction digest), Addgene Protocol

### Fragment Generation

When an enzyme cuts at position p in a sequence of length L:
- Single cut produces 2 fragments: [0, p) and [p, L)
- n cuts produce n+1 fragments

Source: Wikipedia (Restriction digest)

### Gel Electrophoresis Validation

In laboratory practice, restriction digest results are verified by:
1. Running digest products on an agarose gel
2. Comparing band sizes to a DNA ladder
3. Verifying that fragment sizes sum to the expected total

Source: Addgene Protocol, Wikipedia (Gel electrophoresis)

### Enzyme Compatibility

Two enzymes produce compatible (ligatable) ends if:
- Both produce blunt ends, OR
- Both produce identical sticky-end overhangs

Examples:
- BamHI (GATC overhang) and BglII (GATC overhang) are compatible
- EcoRV (blunt) and SmaI (blunt) are compatible
- EcoRI (AATT overhang) and PstI (TGCA overhang) are NOT compatible

Source: Wikipedia (Restriction enzyme), NEB

## Implementation Notes

### Key Methods

```csharp
// Simulate restriction digest
public static IEnumerable<DigestFragment> Digest(DnaSequence sequence, params string[] enzymeNames)

// Get digest summary with fragment sizes
public static DigestSummary GetDigestSummary(DnaSequence sequence, params string[] enzymeNames)

// Create comprehensive restriction map
public static RestrictionMap CreateMap(DnaSequence sequence, params string[] enzymeNames)

// Check enzyme compatibility
public static bool AreCompatible(string enzyme1Name, string enzyme2Name)

// Find all compatible enzyme pairs
public static IEnumerable<(string Enzyme1, string Enzyme2, string CompatibleEnd)> FindCompatibleEnzymes()
```

### DigestFragment Record

```csharp
public sealed record DigestFragment(
    string Sequence,        // The DNA sequence of this fragment
    int StartPosition,      // Start position in original sequence
    int Length,             // Fragment length in base pairs
    string? LeftEnzyme,     // Enzyme that cut left end (null for first fragment)
    string? RightEnzyme,    // Enzyme that cut right end (null for last fragment)
    int FragmentNumber);    // Sequential number (1-based)
```

### DigestSummary Record

```csharp
public sealed record DigestSummary(
    int TotalFragments,
    IReadOnlyList<int> FragmentSizes,    // Sorted descending
    int LargestFragment,
    int SmallestFragment,
    double AverageFragmentSize,
    IReadOnlyList<string> EnzymesUsed);
```

### RestrictionMap Record

```csharp
public sealed record RestrictionMap(
    int SequenceLength,
    IReadOnlyList<RestrictionSite> Sites,
    IReadOnlyDictionary<string, List<int>> SitesByEnzyme,
    int TotalSites,                        // Forward-strand sites only
    IReadOnlyList<string> UniqueCutters,   // Exactly one site
    IReadOnlyList<string> NonCutters);     // Zero sites
```

## Invariants

### Digest Invariants

1. **Fragment Count**: k cut positions → k+1 fragments
2. **Fragment Sum**: Σ(fragment.Length) = sequence.Length
3. **Fragment Numbering**: Fragments numbered 1 to n sequentially
4. **Position Order**: Fragments ordered by ascending start position
5. **Boundary Enzymes**: First fragment LeftEnzyme=null, last fragment RightEnzyme=null
6. **Positive Length**: All fragments have Length > 0 (zero-length fragments not generated)
7. **No Cuts**: Zero cut sites returns single fragment equal to original sequence

### Summary Invariants

8. **Sorted Descending**: FragmentSizes[i] ≥ FragmentSizes[i+1]
9. **Size Bounds**: SmallestFragment ≤ AverageFragmentSize ≤ LargestFragment

### Compatibility Invariants

10. **Blunt-Blunt Compatible**: All blunt enzymes compatible with each other
11. **Matching Overhangs Compatible**: Same overhang sequence = compatible
12. **Symmetry**: AreCompatible(A, B) = AreCompatible(B, A)

## Related Algorithms

- `FindSites`: Identifies restriction sites (RESTR-FIND-001)
- Gel electrophoresis simulation (not implemented)

## References

1. Wikipedia: Restriction digest - https://en.wikipedia.org/wiki/Restriction_digest
2. Wikipedia: Restriction enzyme - https://en.wikipedia.org/wiki/Restriction_enzyme
3. Wikipedia: Restriction map - https://en.wikipedia.org/wiki/Restriction_map
4. Addgene: Restriction Digest Protocol - https://www.addgene.org/protocols/restriction-digest/
5. Roberts RJ (1976) - Restriction endonucleases, CRC Critical Reviews in Biochemistry
6. REBASE (Restriction Enzyme Database) - http://rebase.neb.com/
