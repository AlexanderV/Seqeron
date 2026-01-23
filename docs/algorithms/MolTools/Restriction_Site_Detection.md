# Restriction Site Detection

## Overview

Restriction site detection is the computational identification of DNA sequences recognized by restriction endonucleases (restriction enzymes). These enzymes cleave DNA at specific recognition sequences, producing defined fragments used in molecular cloning, DNA mapping, and genetic engineering.

## Algorithm Description

### Core Algorithm: `FindSites`

The `FindSites` method scans a DNA sequence for recognition sites of specified restriction enzymes:

1. **Validate Input**
   - Check sequence and enzyme name are not null/empty
   - Retrieve enzyme from database (case-insensitive lookup)
   - Throw `ArgumentException` for unknown enzymes

2. **Forward Strand Search**
   - Slide a window of recognition sequence length across the sequence
   - At each position, check if the substring matches the recognition pattern
   - Support IUPAC ambiguity codes in recognition sequences (e.g., N = any base)
   - Record position, enzyme, and cut position for each match

3. **Reverse Strand Search**
   - Compute reverse complement of the sequence
   - Search for recognition sites on the reverse complement
   - Convert positions back to forward strand coordinates
   - Record matches with `IsForwardStrand = false`

4. **Return Results**
   - Yield `RestrictionSite` records with position, enzyme, strand, cut position, and recognized sequence

### Complexity

- **Time**: O(n × k × m) where n = sequence length, k = number of enzymes, m = recognition sequence length
- **Space**: O(1) for streaming results (excluding input)

### FindAllSites Algorithm

Searches for all 40+ enzymes in the database against the sequence:

```csharp
foreach (var enzyme in Enzymes.Values)
    foreach (var site in FindSitesCore(sequence, enzyme))
        yield return site;
```

## Scientific Background

### Restriction Enzymes

Restriction enzymes (restriction endonucleases) are enzymes that cleave DNA at specific recognition sequences. Type II restriction enzymes are most commonly used in molecular biology because they cleave at or near their recognition site.

Source: Wikipedia (Restriction enzyme), Roberts RJ (1976)

### Recognition Sequences

- **Palindromic**: Most Type II enzymes recognize palindromic sequences (reads same on reverse strand)
- **Length**: Typically 4-8 base pairs
  - 4-cutter: cuts ~every 256 bp (4^4)
  - 6-cutter: cuts ~every 4,096 bp (4^6)
  - 8-cutter: cuts ~every 65,536 bp (4^8)

Source: Wikipedia (Restriction site), Cooper S (2003)

### Overhang Types

| Type | Description | Example |
|------|-------------|---------|
| **5' overhang (sticky)** | Cut position forward < reverse | EcoRI: G↓AATTC / CTTAA↓G |
| **3' overhang (sticky)** | Cut position forward > reverse | PstI: CTGCA↓G / G↓ACGTC |
| **Blunt end** | Cut position forward = reverse | EcoRV: GAT↓ATC / CTA↓TAG |

Source: Wikipedia (Restriction enzyme), Wikipedia (EcoRI)

### IUPAC Ambiguity Codes

Some enzymes have degenerate recognition sequences using IUPAC codes:
- N = any nucleotide (A, C, G, T)
- R = purine (A, G)
- Y = pyrimidine (C, T)
- W = weak (A, T)
- S = strong (C, G)

Source: IUPAC-IUB (1970)

## Implementation Notes

### Enzyme Database

The implementation includes 40+ common restriction enzymes:

| Category | Examples | Recognition Length |
|----------|----------|-------------------|
| 4-cutters | AluI, MspI, TaqI | 4 bp |
| 6-cutters | EcoRI, BamHI, HindIII | 6 bp |
| 8-cutters | NotI, PacI, AscI | 8 bp |

### Key Methods

```csharp
// Find sites for specific enzyme(s)
public static IEnumerable<RestrictionSite> FindSites(DnaSequence sequence, string enzymeName)
public static IEnumerable<RestrictionSite> FindSites(DnaSequence sequence, params string[] enzymeNames)
public static IEnumerable<RestrictionSite> FindSites(string sequence, string enzymeName)

// Find sites for all known enzymes
public static IEnumerable<RestrictionSite> FindAllSites(DnaSequence sequence)

// Enzyme lookup
public static RestrictionEnzyme? GetEnzyme(string name)
```

### RestrictionSite Record

```csharp
public sealed record RestrictionSite(
    int Position,           // Start position of recognition sequence
    RestrictionEnzyme Enzyme,
    bool IsForwardStrand,   // True for forward, false for reverse complement
    int CutPosition,        // Position where enzyme cuts
    string RecognizedSequence);
```

## Invariants

1. **Position Range**: 0 ≤ Position ≤ sequence.Length - recognitionLength
2. **Cut Position**: Position ≤ CutPosition ≤ Position + recognitionLength
3. **Recognition Sequence Length**: RecognizedSequence.Length == Enzyme.RecognitionSequence.Length
4. **Case Insensitivity**: Enzyme lookup is case-insensitive
5. **Empty Sequence**: Returns no sites for empty/null sequence
6. **Unknown Enzyme**: Throws `ArgumentException` for unknown enzyme names
7. **Both Strands**: Palindromic sites may appear on both strands at same position

## Related Algorithms

- `Digest`: Simulates restriction digestion to produce fragments (RESTR-DIGEST-001)
- `CreateMap`: Creates comprehensive restriction map
- `AreCompatible`: Checks if two enzymes produce compatible sticky ends

## References

1. Wikipedia: Restriction enzyme - https://en.wikipedia.org/wiki/Restriction_enzyme
2. Wikipedia: Restriction site - https://en.wikipedia.org/wiki/Restriction_site
3. Wikipedia: EcoRI - https://en.wikipedia.org/wiki/EcoRI
4. Roberts RJ (1976) - Restriction endonucleases, CRC Critical Reviews in Biochemistry
5. REBASE (Restriction Enzyme Database) - http://rebase.neb.com/
6. IUPAC-IUB (1970) - Nucleic acid notation
