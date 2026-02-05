# RNA Secondary Structure Prediction

**Algorithm Group:** RNA Structure  
**Implementation:** `Seqeron.Genomics.Analysis.RnaSecondaryStructure`  
**Version:** 1.0  
**Last Updated:** 2026-02-05

---

## 1. Overview

RNA secondary structure prediction determines the base pairing interactions within an RNA molecule from its primary sequence. The secondary structure is critical for RNA function, as it determines binding sites, catalytic activity, and overall three-dimensional folding.

---

## 2. Theoretical Background

### 2.1 Base Pairing

RNA secondary structure is primarily determined by Watson-Crick and wobble base pairing:

| Pair Type | Bases | Hydrogen Bonds | Stability |
|-----------|-------|----------------|-----------|
| Watson-Crick | A-U | 2 | Standard |
| Watson-Crick | G-C | 3 | Strong |
| Wobble | G-U | 2 | Valid in RNA |

**Source:** Wikipedia (Nucleic acid secondary structure), IUPAC conventions

### 2.2 Structural Motifs

| Motif | Description |
|-------|-------------|
| Stem | Double-stranded helical region |
| Hairpin Loop | Unpaired bases at end of stem |
| Internal Loop | Unpaired bases within helix |
| Bulge | Asymmetric unpaired region |
| Multi-loop | Junction of 3+ stems |
| Pseudoknot | Crossing base pairs (not nested) |

### 2.3 Nussinov Algorithm

The Nussinov algorithm uses dynamic programming to maximize the number of base pairs:

**Recurrence:**
```
M(i,j) = max{
  M(i, j-1),                           // j unpaired
  max_{i≤k<j} [M(i,k-1) + M(k+1,j-1) + score(k,j)]  // k pairs with j
}
```

**Complexity:**
- Time: O(n³)
- Space: O(n²)

**Source:** Nussinov & Jacobson (1980), PNAS 77(11):6309-6313

### 2.4 Minimum Free Energy (MFE) Approach

The MFE method uses thermodynamic parameters to find the most stable structure:

**Principle:** Structure with lowest free energy is most probable.

**Energy Model:** Turner 2004 nearest-neighbor parameters
- Stacking energies between adjacent base pairs
- Loop initiation costs
- Special bonuses (tetraloops like GNRA, UNCG)

**Source:** Mathews et al. (2004), PNAS 101(19):7287-7292

---

## 3. Implementation Details

### 3.1 Core Methods

```csharp
// Structure prediction
SecondaryStructure PredictStructure(string rnaSequence, 
    int minStemLength = 3, int minLoopSize = 3, int maxLoopSize = 10)

// Stem-loop finding
IEnumerable<StemLoop> FindStemLoops(string rnaSequence,
    int minStemLength = 3, int minLoopSize = 3, int maxLoopSize = 10,
    bool allowWobble = true)

// MFE calculation
double CalculateMinimumFreeEnergy(string rnaSequence, int minLoopSize = 3)

// Base pairing utilities
bool CanPair(char base1, char base2)
BasePairType? GetBasePairType(char base1, char base2)
```

### 3.2 Data Structures

```csharp
// Base pair representation
record struct BasePair(int Position1, int Position2, char Base1, char Base2, BasePairType Type);

// Stem structure
record struct Stem(int Start5Prime, int End5Prime, int Start3Prime, int End3Prime,
    int Length, IReadOnlyList<BasePair> BasePairs, double FreeEnergy);

// Loop structure
record struct Loop(LoopType Type, int Start, int End, int Size, string Sequence);

// Complete structure
record struct SecondaryStructure(string Sequence, string DotBracket,
    IReadOnlyList<BasePair> BasePairs, IReadOnlyList<StemLoop> StemLoops,
    IReadOnlyList<Pseudoknot> Pseudoknots, double MinimumFreeEnergy);
```

### 3.3 Dot-Bracket Notation

Standard representation of secondary structure:
- `.` = unpaired base
- `(` = base paired with a downstream base
- `)` = base paired with an upstream base
- `[`, `]`, `{`, `}`, `<`, `>` = for pseudoknots and complex structures

**Example:**
```
Sequence:    GGGGAAAACCCC
Dot-bracket: ((((....))))
```

**Source:** Wikipedia (Nucleic acid secondary structure)

---

## 4. Energy Parameters

### 4.1 Stacking Energies (Turner 2004)

Selected nearest-neighbor stacking energies at 37°C (kcal/mol):

| Stack | Energy |
|-------|--------|
| GC/CG | -3.4 |
| CG/GC | -2.4 |
| GG/CC | -3.3 |
| AU/UA | -1.1 |
| UA/AU | -1.3 |
| AA/UU | -0.9 |
| GU/UG | -1.3 |

### 4.2 Hairpin Loop Energies

| Loop Size | Energy (kcal/mol) |
|-----------|-------------------|
| 3 | 5.4 |
| 4 | 5.6 |
| 5 | 5.7 |
| 6 | 5.4 |
| 7 | 6.0 |

### 4.3 Tetraloop Bonuses

| Loop | Bonus |
|------|-------|
| GAAA (GNRA) | -3.0 |
| UUCG (UNCG) | -3.0 |
| CUUG (CUYG) | -2.0 |

---

## 5. Limitations

1. **Pseudoknot Prediction:** Standard DP cannot efficiently predict pseudoknots (O(n⁶) complexity for general pseudoknots)

2. **Simplified Energy Model:** Implementation uses a subset of Turner parameters; full model is more complex

3. **Single Sequence:** No comparative/consensus structure prediction

4. **Modified Bases:** Does not account for post-transcriptional modifications

---

## 6. Usage Example

```csharp
// Predict structure for a hairpin sequence
string rna = "GCGCGAAACGCGC";
var structure = RnaSecondaryStructure.PredictStructure(rna);

Console.WriteLine($"Sequence: {structure.Sequence}");
Console.WriteLine($"Structure: {structure.DotBracket}");
Console.WriteLine($"MFE: {structure.MinimumFreeEnergy} kcal/mol");
Console.WriteLine($"Base pairs: {structure.BasePairs.Count}");

// Find stem-loops specifically
var stemLoops = RnaSecondaryStructure.FindStemLoops(rna, minStemLength: 3).ToList();
foreach (var sl in stemLoops)
{
    Console.WriteLine($"Stem-loop at {sl.Start}-{sl.End}: {sl.DotBracketNotation}");
}

// Calculate MFE directly
double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna);
```

---

## 7. References

1. Nussinov R, Jacobson AB (1980). "Fast algorithm for predicting the secondary structure of single-stranded RNA." PNAS 77(11):6309-6313.

2. Zuker M, Stiegler P (1981). "Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information." Nucleic Acids Res 9(1):133-148.

3. Mathews DH, Disney MD, Childs JL, Schroeder SJ, Zuker M, Turner DH (2004). "Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure." PNAS 101(19):7287-7292.

4. Wikipedia contributors. "Nucleic acid structure prediction." Wikipedia.

5. Wikipedia contributors. "Nussinov algorithm." Wikipedia.

---

## 8. See Also

- [RNA-STEMLOOP-001](./RNA_Stemloop.md) - Stem-loop detection
- [RNA-ENERGY-001](./RNA_Energy.md) - Free energy calculation
