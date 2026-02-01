# IUPAC Degenerate Motif Matching

**Test Unit ID:** PAT-IUPAC-001
**Algorithm Group:** Pattern Matching
**Implementation:** `MotifFinder.FindDegenerateMotif()`, `IupacHelper.MatchesIupac()`

---

## 1. Definition

IUPAC degenerate motif matching extends exact pattern matching to handle ambiguity codes in DNA/RNA patterns. Each position in the pattern can specify one or more acceptable nucleotides using the IUPAC notation system.

### Formal Problem Statement

**Input:**
- Sequence S of length n (DNA: A, C, G, T only)
- Pattern P of length m (may contain IUPAC ambiguity codes)

**Output:**
- All positions i where S[i..i+m-1] matches P according to IUPAC rules

### IUPAC Nucleotide Code Table

| Code | Represents | Mnemonic | Count |
|------|------------|----------|-------|
| A | A | Adenine | 1 |
| C | C | Cytosine | 1 |
| G | G | Guanine | 1 |
| T | T | Thymine | 1 |
| R | A, G | puRine | 2 |
| Y | C, T | pYrimidine | 2 |
| S | G, C | Strong (3 H-bonds) | 2 |
| W | A, T | Weak (2 H-bonds) | 2 |
| K | G, T | Keto | 2 |
| M | A, C | aMino | 2 |
| B | C, G, T | not A | 3 |
| D | A, G, T | not C | 3 |
| H | A, C, T | not G | 3 |
| V | A, C, G | not T | 3 |
| N | A, C, G, T | aNy | 4 |

**Source:** IUPAC-IUB Commission on Biochemical Nomenclature (1970), NC-IUB (1984)

---

## 2. Algorithm

### Matching Procedure

For each position i in sequence S where pattern P fits:
```
match = true
for j = 0 to m-1:
    if not MatchesIupac(S[i+j], P[j]):
        match = false
        break
if match:
    report position i
```

### MatchesIupac Function

```
MatchesIupac(nucleotide n, IUPAC code c):
    case c of:
        'A': return n == 'A'
        'C': return n == 'C'
        'G': return n == 'G'
        'T': return n == 'T'
        'N': return n ∈ {A, C, G, T}
        'R': return n ∈ {A, G}
        'Y': return n ∈ {C, T}
        'S': return n ∈ {G, C}
        'W': return n ∈ {A, T}
        'K': return n ∈ {G, T}
        'M': return n ∈ {A, C}
        'B': return n ∈ {C, G, T}  // not A
        'D': return n ∈ {A, G, T}  // not C
        'H': return n ∈ {A, C, T}  // not G
        'V': return n ∈ {A, C, G}  // not T
        default: return n == c     // exact match for unknown codes
```

### Complexity

| Operation | Complexity | Notes |
|-----------|------------|-------|
| Single position match | O(m) | Pattern length |
| Full sequence scan | O(n × m) | Brute force |
| Space | O(1) | No auxiliary data structures |

---

## 3. Implementation Details

### IupacHelper Class

```csharp
public static class IupacHelper
{
    public static bool MatchesIupac(char nucleotide, char iupacCode) => iupacCode switch
    {
        'A' => nucleotide == 'A',
        'C' => nucleotide == 'C',
        'G' => nucleotide == 'G',
        'T' => nucleotide == 'T',
        'N' => nucleotide is 'A' or 'C' or 'G' or 'T',
        'R' => nucleotide is 'A' or 'G',          // puRine
        'Y' => nucleotide is 'C' or 'T',          // pYrimidine
        'S' => nucleotide is 'G' or 'C',          // Strong
        'W' => nucleotide is 'A' or 'T',          // Weak
        'K' => nucleotide is 'G' or 'T',          // Keto
        'M' => nucleotide is 'A' or 'C',          // aMino
        'B' => nucleotide is 'C' or 'G' or 'T',   // not A
        'D' => nucleotide is 'A' or 'G' or 'T',   // not C
        'H' => nucleotide is 'A' or 'C' or 'T',   // not G
        'V' => nucleotide is 'A' or 'C' or 'G',   // not T
        _ => nucleotide == iupacCode              // fallback: exact match
    };
}
```

### MotifFinder.FindDegenerateMotif

```csharp
public static IEnumerable<MotifMatch> FindDegenerateMotif(DnaSequence sequence, string motif)
{
    ArgumentNullException.ThrowIfNull(sequence);
    if (string.IsNullOrEmpty(motif)) yield break;

    string seq = sequence.Sequence;
    string motifUpper = motif.ToUpperInvariant();

    for (int i = 0; i <= seq.Length - motifUpper.Length; i++)
    {
        bool matches = true;
        for (int j = 0; j < motifUpper.Length && matches; j++)
        {
            char motifChar = motifUpper[j];
            char seqChar = seq[i + j];

            if (IupacCodes.TryGetValue(motifChar, out string? allowed))
                matches = allowed.Contains(seqChar);
            else
                matches = motifChar == seqChar;
        }

        if (matches)
        {
            yield return new MotifMatch(
                Position: i,
                MatchedSequence: seq.Substring(i, motifUpper.Length),
                Pattern: motifUpper,
                Score: 1.0);
        }
    }
}
```

**Implementation Note:** The MotifFinder uses an internal dictionary lookup rather than calling IupacHelper directly. Both implementations are equivalent and follow the IUPAC specification.

### Edge Case Handling

| Condition | Behavior | Rationale |
|-----------|----------|-----------|
| Null sequence | ArgumentNullException | Standard .NET convention |
| Empty pattern | Returns empty | No pattern to match |
| Empty sequence | Returns empty | No content to match |
| Pattern longer than sequence | Returns empty | Cannot match |
| Lowercase input | Normalized to uppercase | Case-insensitive matching |
| Unknown IUPAC code | Exact match | Fallback behavior |

---

## 4. Biological Applications

### Common Degenerate Motifs

| Motif | Pattern | Description |
|-------|---------|-------------|
| E-box | CANNTG | Transcription factor binding site |
| TATA box | TATAAA | Promoter element |
| Kozak sequence | GCCGCCRCCATG | Translation initiation |
| Restriction sites | Various | Many use degenerate codes |

### Use Cases

1. **Motif discovery:** Searching for consensus binding sites
2. **Primer design:** Degenerate primers for related sequences
3. **SNP tolerance:** Patterns that match multiple alleles
4. **Consensus matching:** Finding sequences matching a profile

---

## 5. References

1. IUPAC-IUB Commission on Biochemical Nomenclature (1970). "Abbreviations and symbols for nucleic acids, polynucleotides, and their constituents." Biochemistry 9(20):4022–4027. doi:10.1021/bi00822a023

2. NC-IUB (1984). "Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences." Nucleic Acids Research 13(9):3021–3030. doi:10.1093/nar/13.9.3021

3. Wikipedia contributors. "Nucleic acid notation." Wikipedia, The Free Encyclopedia. https://en.wikipedia.org/wiki/Nucleic_acid_notation

4. Bioinformatics.org. "IUPAC Codes." https://www.bioinformatics.org/sms/iupac.html

---

## 6. Implementation-Specific Notes

### Deviations from IUPAC Standard

- None identified. Implementation follows IUPAC-IUB 1970/1984 specifications exactly.

### Assumptions

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | Sequence contains only uppercase A, C, G, T | DnaSequence class guarantees this |
| A2 | Unknown pattern codes fall back to exact match | Safe default behavior |
