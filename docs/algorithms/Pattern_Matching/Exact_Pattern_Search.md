# Exact Pattern Search (Suffix Tree)

**Test Unit ID:** PAT-EXACT-001
**Algorithm Group:** Pattern Matching
**Implementation:** `SuffixTree` class (`SuffixTree.Search.cs`)

---

## 1. Definition

Exact pattern matching is the problem of finding all occurrences of a pattern string P in a text string T. Using a suffix tree, this problem can be solved in optimal time.

### Formal Problem Statement

**Input:**
- Text T of length n (pre-indexed in suffix tree)
- Pattern P of length m

**Output:**
- All positions i where T[i..i+m-1] = P

---

## 2. Algorithm

### Suffix Tree Search (Gusfield, 1997)

The suffix tree for text T contains all suffixes of T as paths from root to leaves. Pattern P occurs at position i if and only if P is a prefix of the suffix starting at position i.

**Search procedure:**
1. Start at the root of the suffix tree
2. Match P character-by-character along tree edges
3. If mismatch occurs before P is exhausted → P not found
4. If P is fully matched → all leaves in subtree below match point correspond to occurrences
5. Collect leaf positions as occurrence positions

### Complexity (from Wikipedia)

| Operation | Complexity | Notes |
|-----------|------------|-------|
| Build suffix tree | O(n) | Ukkonen's algorithm |
| Pattern search | O(m + z) | m = pattern length, z = occurrences |
| Contains check | O(m) | Stop at first match |
| Count occurrences | O(m) | Pre-computed leaf counts |

**Source:** Gusfield (1997), p.92, 123; Wikipedia "Suffix tree" article.

---

## 3. Implementation Details

### SuffixTree Class Methods

```csharp
// Core pattern matching
public IReadOnlyList<int> FindAllOccurrences(string pattern);
public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern);

// Existence check (optimized - stops early)
public bool Contains(string value);
public bool Contains(ReadOnlySpan<char> value);

// Count (uses pre-computed LeafCount - no leaf enumeration)
public int CountOccurrences(string pattern);
public int CountOccurrences(ReadOnlySpan<char> pattern);
```

### Internal Implementation

The implementation uses `MatchPatternCore` for edge-by-edge traversal:

1. **Hybrid SIMD optimization:** Uses `SequenceEqual` for comparisons ≥8 chars, scalar loop for shorter
2. **LeafCount pre-computation:** Each node stores count of leaves in subtree for O(m) counting
3. **Thread-static buffers:** Reuses buffers to minimize allocations

### Edge Case Handling

| Condition | Behavior | Rationale |
|-----------|----------|-----------|
| Null pattern | ArgumentNullException | Standard .NET convention |
| Empty pattern (string) | Returns all positions [0..n-1] | Matches regex "" semantics |
| Empty pattern (Span) | Returns empty list | Span API difference (documented) |
| Empty text | Returns empty for non-empty pattern | No content to match |
| Pattern longer than text | Returns empty | Cannot match |

---

### Related Documentation

- [Suffix Tree (Ukkonen)](Suffix_Tree.md)

---

## 4. Genomics Wrappers

### GenomicAnalyzer.FindMotif

```csharp
public static IReadOnlyList<int> FindMotif(DnaSequence sequence, string motif)
{
    if (string.IsNullOrEmpty(motif)) return Array.Empty<int>();
    string normalizedMotif = motif.ToUpperInvariant();
    return sequence.SuffixTree.FindAllOccurrences(normalizedMotif);
}
```

**Key difference:** Normalizes motif to uppercase for case-insensitive DNA matching.

### MotifFinder.FindExactMotif

```csharp
public static IEnumerable<int> FindExactMotif(DnaSequence sequence, string motif)
{
    // ... validation ...
    string motifUpper = motif.ToUpperInvariant();
    var positions = sequence.SuffixTree.FindAllOccurrences(motifUpper);
    foreach (int pos in positions.OrderBy(p => p))
        yield return pos;
}
```

**Key differences:**
- Returns IEnumerable (lazy enumeration)
- Orders results by position
- Normalizes to uppercase

---

## 5. Test Strings from Literature

### "banana" (Wikipedia Suffix Tree Article)

The string "banana" is the canonical example used in Wikipedia's suffix tree article:

| Pattern | Occurrences (0-indexed) |
|---------|-------------------------|
| "a" | [1, 3, 5] |
| "na" | [2, 4] |
| "ana" | [1, 3] |
| "ban" | [0] |
| "banana" | [0] |
| "nan" | [2] |

### "mississippi" (Gusfield, 1997)

Classic example from Gusfield's textbook:

| Pattern | Occurrences (0-indexed) |
|---------|-------------------------|
| "i" | [1, 4, 7, 10] |
| "s" | [2, 3, 5, 6] |
| "issi" | [1, 4] |
| "pp" | [8] |
| "ss" | [2, 5] |

### Rosalind SUBS Problem

Bioinformatics problem demonstrating overlapping occurrences:

**Text:** "GATATATGCATATACTT"
**Pattern:** "ATAT"
**Occurrences (0-indexed):** [1, 3, 9]

Note: Rosalind uses 1-indexed positions (2, 4, 10).

---

## 6. Invariants

| Property | Description |
|----------|-------------|
| Correctness | Every returned position i satisfies text[i..i+m-1] = pattern |
| Completeness | All positions where pattern occurs are returned |
| Consistency | `CountOccurrences(P)` = `FindAllOccurrences(P).Count` |
| Existence | `Contains(P)` = `FindAllOccurrences(P).Count > 0` |
| Substring property | Every substring of text can be found |

---

## 7. References

1. **Gusfield, D.** (1997). *Algorithms on Strings, Trees and Sequences: Computer Science and Computational Biology*. Cambridge University Press. ISBN 0-521-58519-8.

2. **Ukkonen, E.** (1995). "On-line construction of suffix trees." *Algorithmica*, 14(3), 249-260.

3. **Wikipedia contributors.** "Suffix tree." *Wikipedia, The Free Encyclopedia*. https://en.wikipedia.org/wiki/Suffix_tree

4. **Rosalind Team.** "Finding a Motif in DNA." *Rosalind Bioinformatics*. https://rosalind.info/problems/subs/
