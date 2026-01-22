# Direct Repeat Detection

## Algorithm Overview

| Property | Value |
|----------|-------|
| **Algorithm** | Direct Repeat Detection |
| **Implementation** | `RepeatFinder.FindDirectRepeats` |
| **Complexity** | O(n²) where n = sequence length |
| **Category** | Repeat Analysis |

---

## Definition

A **direct repeat** is a nucleotide sequence that appears two or more times in a genome, with copies oriented in the same 5' to 3' direction. Unlike inverted repeats (which are reverse complements), direct repeats preserve the original sequence orientation.

**Canonical format:**
```
5' TTACG------TTACG 3'
3' AATGC------AATGC 5'
```

Where `------` represents intervening nucleotides (spacing). Adjacent direct repeats (spacing = 0) are also valid.

---

## Sources

| Source | Type | Reference |
|--------|------|-----------|
| Wikipedia - Direct repeat | Definition | https://en.wikipedia.org/wiki/Direct_repeat |
| Wikipedia - Repeated sequence (DNA) | Context | https://en.wikipedia.org/wiki/Repeated_sequence_(DNA) |
| Ussery et al. (2009) | Technical | Computing for Comparative Microbial Genomics, Chapter 8 |
| Richard (2021) | Clinical | PMC8145212 - Trinucleotide repeat expansions |

---

## Theory

### Structural Characteristics

1. **Same directionality**: Both copies read identically in the 5' to 3' direction
2. **Variable spacing**: Copies may be adjacent (tandem) or separated by intervening nucleotides
3. **Length variability**: From short (2-6 bp) to long (hundreds of bp)

### Biological Significance

Direct repeats serve important biological functions:

1. **Transposable elements**: Long terminal repeats (LTRs) flank retrotransposons
2. **DNA recombination**: Direct repeats are hotspots for homologous recombination
3. **Genome instability**: Associated with deletions via replication slippage
4. **Gene regulation**: Some regulatory elements contain direct repeats

### Clinical Relevance

Trinucleotide repeat expansions (a form of tandem direct repeat) underlie several human diseases:

| Disease | Repeat | Gene |
|---------|--------|------|
| Huntington's disease | CAG | HTT |
| Fragile X syndrome | CGG | FMR1 |
| Spinocerebellar ataxias | CAG | Various |
| Friedreich's ataxia | GAA | FXN |
| Myotonic dystrophy | CTG/CCTG | DMPK/ZNF9 |

---

## Implementation

### Method Signature

```csharp
public static IEnumerable<DirectRepeatResult> FindDirectRepeats(
    DnaSequence sequence,
    int minLength = 5,
    int maxLength = 50,
    int minSpacing = 1)
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sequence` | required | DNA sequence to search |
| `minLength` | 5 | Minimum repeat unit length |
| `maxLength` | 50 | Maximum repeat unit length |
| `minSpacing` | 1 | Minimum spacing between repeat copies |

### Return Value

```csharp
public readonly record struct DirectRepeatResult(
    int FirstPosition,    // 0-based position of first copy
    int SecondPosition,   // 0-based position of second copy
    string RepeatSequence,// The repeated sequence
    int Length,           // Length of the repeat
    int Spacing);         // Gap between copies
```

### Algorithm Strategy

The implementation uses a **suffix tree** for efficient pattern matching:

1. Build suffix tree from the input sequence (O(n))
2. For each position and length, extract candidate pattern
3. Use suffix tree to find all occurrences in O(m + k) time
4. Filter occurrences by spacing constraints
5. Report non-duplicate results

### Invariants

1. **Spacing calculation**: `Spacing = SecondPosition - FirstPosition - Length`
2. **Non-overlap** (when minSpacing > 0): `SecondPosition >= FirstPosition + Length + minSpacing`
3. **Result uniqueness**: No duplicate (FirstPosition, SecondPosition, Length) tuples

---

## Implementation-Specific Notes

### Case Handling

The string overload normalizes input to uppercase via `ToUpperInvariant()` before processing.

### Empty/Null Handling

- Null `DnaSequence` input throws `ArgumentNullException`
- Empty string returns empty enumerable (no exception)

### Parameter Validation

- `minLength < 2` throws `ArgumentOutOfRangeException`
- `maxLength < minLength` throws `ArgumentOutOfRangeException`
- `minSpacing` can be 0 (allows adjacent repeats)

---

## Edge Cases

| Case | Behavior |
|------|----------|
| Empty sequence | Returns empty enumerable |
| Sequence too short | Returns empty enumerable |
| No repeats found | Returns empty enumerable |
| Adjacent repeats (spacing=0) | Found when minSpacing=0 |
| Multiple occurrences | Reports pairwise matches |
| Case variation | Normalized to uppercase |

---

## Deviations / Assumptions

| ID | Note | Type |
|----|------|------|
| D1 | Only reports first-second pairs, not all pairwise combinations | Implementation |
| A1 | minLength ≥ 2 prevents single-nucleotide "repeats" | ASSUMPTION |
