# RNA Stem-Loop Detection

**Algorithm Group:** RNA Structure  
**Algorithm Name:** Stem-Loop Detection  
**Implementation:** `Seqeron.Genomics.Analysis.RnaSecondaryStructure.FindStemLoops`  
**Test Unit:** RNA-STEMLOOP-001  
**Version:** 1.0  
**Last Updated:** 2026-02-05

---

## 1. Overview

Stem-loop (hairpin) structures are fundamental RNA secondary structure motifs formed when a single-stranded RNA folds back on itself through Watson-Crick and wobble base pairing. The `FindStemLoops` method detects these structures in RNA sequences.

---

## 2. Biological Background

### 2.1 Structure Definition

A stem-loop consists of two components:

| Component | Description | Constraints |
|-----------|-------------|-------------|
| **Stem** | Double-stranded helical region | Antiparallel base pairing |
| **Loop** | Unpaired nucleotides at stem end | Minimum 3 nucleotides (steric) |

**Diagram:**
```
      L O O P
      A-A-A-A
     /       \
    G         C    <- Loop closing pair
    |         |
    G---------C    <- Stem (base pairs)
    |         |
    G---------C
    |         |
    C---------G
   5'         3'
```

**Source:** Wikipedia (Stem-loop)

### 2.2 Valid Base Pairs

| Pair Type | Bases | Stability |
|-----------|-------|-----------|
| Watson-Crick | A-U, U-A | Standard |
| Watson-Crick | G-C, C-G | Strong (3 H-bonds) |
| Wobble | G-U, U-G | Valid in RNA |

**Source:** IUPAC conventions, Wikipedia

### 2.3 Loop Size Constraints

| Constraint | Value | Reason |
|------------|-------|--------|
| Minimum | 3 nucleotides | Steric impossibility for <3 |
| Optimal | 4-8 nucleotides | Balance of stability |
| Maximum | No hard limit | Larger loops are less stable |

**Source:** Wikipedia (Stem-loop) - "loops that are fewer than three bases long are sterically impossible"

### 2.4 Tetraloops

Special 4-nucleotide loops with enhanced stability:

| Family | Pattern | Examples | Bonus |
|--------|---------|----------|-------|
| GNRA | G-N-R-A | GAAA, GCAA, GGAA, GUAA | ~3 kcal/mol |
| UNCG | U-N-C-G | UUCG, UACG, UGCG | ~3 kcal/mol |
| CUUG | C-U-U-G | CUUG | ~2 kcal/mol |

**Source:** Wikipedia (Tetraloop), Woese et al. (1990), Heus & Pardi (1991)

---

## 3. Algorithm

### 3.1 Method Signature

```csharp
public static IEnumerable<StemLoop> FindStemLoops(
    string rnaSequence,
    int minStemLength = 3,
    int minLoopSize = 3,
    int maxLoopSize = 10,
    bool allowWobble = true)
```

### 3.2 Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `rnaSequence` | - | RNA sequence to analyze |
| `minStemLength` | 3 | Minimum base pairs in stem |
| `minLoopSize` | 3 | Minimum loop nucleotides |
| `maxLoopSize` | 10 | Maximum loop nucleotides |
| `allowWobble` | true | Allow G-U wobble pairs |

### 3.3 Algorithm Steps

1. **Validation:** Check sequence length is sufficient
2. **Loop Scanning:** Iterate potential loop start positions
3. **Loop Sizing:** For each position, try loop sizes in range
4. **Stem Extension:** Extend stem by checking base pair validity
5. **Collection:** Return structures meeting minimum requirements

### 3.4 Complexity

- **Time:** O(n² × L) where n = sequence length, L = max loop size
- **Space:** O(k) where k = number of stem-loops found

---

## 4. Return Type

### 4.1 StemLoop Record

```csharp
public readonly record struct StemLoop(
    int Start,                    // Start position
    int End,                      // End position
    Stem Stem,                    // Stem details
    Loop Loop,                    // Loop details
    double TotalFreeEnergy,       // Calculated energy
    string DotBracketNotation);   // Structure notation
```

### 4.2 Stem Record

```csharp
public readonly record struct Stem(
    int Start5Prime,
    int End5Prime,
    int Start3Prime,
    int End3Prime,
    int Length,
    IReadOnlyList<BasePair> BasePairs,
    double FreeEnergy);
```

### 4.3 Loop Record

```csharp
public readonly record struct Loop(
    LoopType Type,   // Always Hairpin for stem-loops
    int Start,
    int End,
    int Size,
    string Sequence);
```

---

## 5. Related Methods

### 5.1 Pseudoknot Detection

```csharp
public static IEnumerable<Pseudoknot> DetectPseudoknots(
    IReadOnlyList<BasePair> basePairs)
```

Detects crossing base pairs that indicate pseudoknots. A pseudoknot occurs when pairs (i,j) and (k,l) satisfy: i < k < j < l.

**Source:** Wikipedia (Pseudoknot)

### 5.2 Inverted Repeat Finding

```csharp
public static IEnumerable<(int, int, int, int, int)> FindInvertedRepeats(
    string sequence,
    int minLength = 4,
    int minSpacing = 3,
    int maxSpacing = 100)
```

Finds potential stem regions based on antiparallel complementarity.

---

## 6. Usage Examples

### 6.1 Basic Usage

```csharp
var rna = "GGGAAAACCC";
var stemLoops = RnaSecondaryStructure.FindStemLoops(rna);

foreach (var sl in stemLoops)
{
    Console.WriteLine($"Stem: {sl.Stem.Length} bp");
    Console.WriteLine($"Loop: {sl.Loop.Sequence}");
    Console.WriteLine($"Structure: {sl.DotBracketNotation}");
}
```

### 6.2 Custom Parameters

```csharp
// Find only strong hairpins with 4-nucleotide tetraloops
var stemLoops = RnaSecondaryStructure.FindStemLoops(
    rnaSequence,
    minStemLength: 4,
    minLoopSize: 4,
    maxLoopSize: 4,
    allowWobble: false);
```

---

## 7. Limitations

1. **No pseudoknot prediction:** Cannot predict pseudoknots from sequence alone
2. **Simplified energy model:** Uses Turner 2004 parameters, less accurate than ViennaRNA
3. **No internal loops:** Focus is on hairpin loops only
4. **Overlapping structures:** Returns all candidates, may need filtering

---

## 8. References

1. Wikipedia. "Stem-loop." https://en.wikipedia.org/wiki/Stem-loop
2. Wikipedia. "Tetraloop." https://en.wikipedia.org/wiki/Tetraloop
3. Wikipedia. "Pseudoknot." https://en.wikipedia.org/wiki/Pseudoknot
4. Woese CR, Winker S, Gutell RR. (1990). "Architecture of ribosomal RNA: constraints on the sequence of 'tetra-loops'." PNAS 87(21):8467-8471.
5. Heus HA, Pardi A. (1991). "Structural features that give rise to the unusual stability of RNA hairpins containing GNRA loops." Science 253(5016):191-194.
6. Svoboda P, Cara A. (2006). "Hairpin RNA: A secondary structure of primary importance." Cell Mol Life Sci 63(7):901-908.

---

## 9. See Also

- [RNA Secondary Structure](./RNA_Secondary_Structure.md) - Full structure prediction
- [RNA-STRUCT-001](../../tests/TestSpecs/RNA-STRUCT-001.md) - Structure prediction tests
- [RNA-STEMLOOP-001](../../tests/TestSpecs/RNA-STEMLOOP-001.md) - Stem-loop tests
