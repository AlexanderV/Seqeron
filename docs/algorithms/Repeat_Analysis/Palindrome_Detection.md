# Palindrome Detection

## Algorithm Overview

| Property | Value |
|----------|-------|
| **Algorithm** | DNA Palindrome Detection |
| **Implementation** | `RepeatFinder.FindPalindromes` |
| **Complexity** | O(n² × L) where n = sequence length, L = length range |
| **Category** | Repeat Analysis |

---

## Definition

A **DNA palindrome** (also called "reverse palindrome" or "biological palindrome") is a nucleotide sequence that equals its reverse complement. This differs from textual palindromes, which read the same forwards and backwards on a single strand.

**Key principle:** Reading 5'→3' on the forward strand gives the same sequence as reading 5'→3' on the reverse strand.

**Example: GAATTC (EcoRI recognition site)**
```
5'-GAATTC-3'   (forward strand, read 5'→3')
3'-CTTAAG-5'   (reverse strand, read 3'→5')
     ↓
5'-GAATTC-3'   (reverse strand, read 5'→3')
```

Mathematically:
```
Palindrome(S) ⟺ S == ReverseComplement(S)
```

---

## Sources

| Source | Type | Reference |
|--------|------|-----------|
| Wikipedia - Palindromic sequence | Definition | https://en.wikipedia.org/wiki/Palindromic_sequence |
| Wikipedia - Restriction enzyme | Application | https://en.wikipedia.org/wiki/Restriction_enzyme |
| Rosalind - REVP | Test data | https://rosalind.info/problems/revp/ |
| REBASE | Database | https://rebase.neb.com/ |

---

## Theory

### Structural Characteristics

1. **Even length required**: DNA palindromes must have even length because each base pairs with a complementary base in the opposite strand. An odd-length sequence would have a central base that cannot complement itself.

2. **Symmetry axis**: The symmetry axis lies between the two central nucleotides.

3. **Self-complementary**: The sequence can form a perfect duplex with itself.

### Mathematical Verification

For a sequence S of even length 2n:
- S = s₁s₂...s₂ₙ
- S is a palindrome iff sᵢ = complement(s₂ₙ₊₁₋ᵢ) for all i ∈ [1, 2n]

Where complement mapping is: A↔T, G↔C

### Biological Significance

DNA palindromes serve important biological functions:

1. **Restriction enzyme recognition**: Type II restriction endonucleases recognize palindromic sequences (4-8 bp typically)
2. **Cruciform structures**: Long palindromes can form cruciform DNA under supercoiling
3. **Hairpin formation**: Single-stranded palindromes can fold into hairpin structures
4. **Genome instability**: Long palindromes are associated with chromosomal fragility

### Common Restriction Enzyme Sites

| Enzyme | Sequence | Length | Source Organism |
|--------|----------|--------|-----------------|
| EcoRI | GAATTC | 6 | Escherichia coli |
| BamHI | GGATCC | 6 | Bacillus amyloliquefaciens |
| HindIII | AAGCTT | 6 | Haemophilus influenzae |
| TaqI | TCGA | 4 | Thermus aquaticus |
| AluI | AGCT | 4 | Arthrobacter luteus |
| NotI | GCGGCCGC | 8 | Nocardia otitidis |
| SmaI | CCCGGG | 6 | Serratia marcescens |
| EcoRV | GATATC | 6 | Escherichia coli |

---

## Implementation

### Method Signature

```csharp
public static IEnumerable<PalindromeResult> FindPalindromes(
    DnaSequence sequence,
    int minLength = 4,
    int maxLength = 12)
```

### Parameters

| Parameter | Default | Description | Constraints |
|-----------|---------|-------------|-------------|
| `sequence` | required | DNA sequence to search | Not null |
| `minLength` | 4 | Minimum palindrome length | Must be ≥ 4 and even |
| `maxLength` | 12 | Maximum palindrome length | Must be ≥ minLength |

### Return Value

```csharp
public readonly record struct PalindromeResult(
    int Position,    // 0-based position in sequence
    string Sequence, // The palindromic sequence
    int Length);     // Length of the palindrome
```

### Algorithm

```
FOR each length L from minLength to maxLength, step 2:
    FOR each position i from 0 to (sequenceLength - L):
        candidate = sequence[i..i+L]
        IF candidate == ReverseComplement(candidate):
            YIELD PalindromeResult(i, candidate, L)
```

### Complexity Analysis

- **Time**: O(n × L × k) where n = sequence length, L = length range / 2, k = cost of reverse complement comparison (O(length))
- **Space**: O(L_max) for temporary strings

### Constraints

1. **minLength ≥ 4**: Prevents trivial 2bp matches which are extremely common
2. **Even lengths only**: Biological palindromes require even length
3. **Length step of 2**: Only even lengths checked

---

## Implementation Notes

### String Overload

The string overload converts input to uppercase before processing:

```csharp
public static IEnumerable<PalindromeResult> FindPalindromes(
    string sequence,
    int minLength = 4,
    int maxLength = 12)
{
    if (string.IsNullOrEmpty(sequence))
        yield break;

    foreach (var result in FindPalindromesCore(
        sequence.ToUpperInvariant(), minLength, maxLength))
        yield return result;
}
```

### GenomicAnalyzer Alternate Implementation

`GenomicAnalyzer.FindPalindromes` provides an alternate implementation with different return type (`PalindromeInfo` vs `PalindromeResult`). Both implementations use the same algorithm but have slightly different APIs:

| Aspect | RepeatFinder | GenomicAnalyzer |
|--------|--------------|-----------------|
| Return type | `PalindromeResult` (record struct) | `PalindromeInfo` (struct) |
| Null handling | Throws `ArgumentNullException` | Not documented |
| Length validation | Explicit exception for odd minLength | No explicit validation |

---

## Edge Cases

| Case | Behavior | Rationale |
|------|----------|-----------|
| Empty sequence | Returns empty | No positions to check |
| Sequence shorter than minLength | Returns empty | Cannot contain palindrome |
| No palindromes in sequence | Returns empty | Expected behavior |
| Entire sequence is palindrome | Reports at valid positions | Standard behavior |
| Overlapping palindromes (different lengths) | Both reported | 4bp and 6bp palindrome at same start |
| Odd minLength | Throws ArgumentOutOfRangeException | Biological constraint |
| minLength < 4 | Throws ArgumentOutOfRangeException | Implementation constraint |
| maxLength < minLength | Throws ArgumentOutOfRangeException | Invalid parameter |

---

## Relation to Other Algorithms

### Inverted Repeat Detection (REP-INV-001)

Inverted repeats are related but distinct:
- **Palindrome**: Sequence equals its reverse complement (no gap)
- **Inverted repeat**: Two sequences that are reverse complements, separated by a loop

A palindrome is essentially an inverted repeat with loop length = 0 and overlapping arms.

### Restriction Site Analysis (RESTR-FIND-001)

Palindrome detection is the foundation for restriction site analysis. Most Type II restriction enzymes recognize specific palindromic sequences.

---

## Test Coverage

See [TestSpec REP-PALIN-001](../../tests/TestSpecs/REP-PALIN-001.md) for complete test specification.

Key test categories:
1. Known restriction enzyme sites (EcoRI, BamHI, HindIII, etc.)
2. Parameter validation (even length, range checks)
3. Edge cases (empty, too short, no matches)
4. Rosalind REVP sample dataset validation
