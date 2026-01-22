# Inverted Repeat Detection

**Test Unit:** REP-INV-001  
**Last Updated:** 2026-01-22  
**Status:** Active

---

## 1. Definition

An **inverted repeat** (IR) is a single-stranded sequence of nucleotides followed downstream by its **reverse complement**. The intervening sequence between the initial sequence and its reverse complement can be any length, including zero.

### Structure

```
5'---TTACG------nnnnnn------CGTAA---3'
     └─Left Arm─┘   └─Loop─┘  └─Right Arm─┘
                              (reverse complement of left arm)
```

Where:
- **Left Arm:** Initial sequence
- **Loop:** Intervening nucleotides (n ≥ 0)
- **Right Arm:** Reverse complement of the left arm

### Terminology
- **Stem-loop (Hairpin):** When an inverted repeat folds back on itself to form a double helix ending in a loop
- **Palindrome:** An inverted repeat with no intervening nucleotides (loop = 0)
- **Cruciform:** Extruded structure from inverted repeat regions

### Sources
- Wikipedia: [Inverted repeat](https://en.wikipedia.org/wiki/Inverted_repeat)
- Wikipedia: [Stem-loop](https://en.wikipedia.org/wiki/Stem-loop)
- Wikipedia: [Palindromic sequence](https://en.wikipedia.org/wiki/Palindromic_sequence)
- EMBOSS: [einverted](https://emboss.sourceforge.net/apps/cvs/emboss/apps/einverted.html)

---

## 2. Biological Significance

### Structural Role
Inverted repeats can form secondary structures:
- **Hairpins/Stem-loops:** Optimal loop length is 4–8 bases; loops < 3 bases are sterically impossible
- **Cruciforms:** Four-way junctions formed by extrusion of inverted repeats in double-stranded DNA
- **Pseudoknots:** Nested stem-loops in RNA

### Functional Roles
| Function | Description |
|----------|-------------|
| Replication origins | Found at origins of replication in phages, plasmids, mitochondria, eukaryotic viruses |
| Transposon boundaries | Terminal inverted repeats (TIRs) define transposable element boundaries |
| Transcription termination | Rho-independent terminators use stem-loop structures |
| Riboswitches | RNA regulatory elements using alternative stem-loop conformations |
| tRNA structure | Three stem-loops in cloverleaf structure |

### Disease Associations
Inverted repeats are "hotspots" of genomic instability:
- Long IRs are often deleted in E. coli and yeast
- Associated with recombination and mutations in mammalian cells
- Can lead to frameshift mutations and point mutations

**Sources:** Bissler (1998), Mirkin et al. (2008), Pearson et al. (1996)

---

## 3. Algorithm Description

### Input Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sequence` | DnaSequence or string | - | DNA sequence to search |
| `minArmLength` | int | 4 | Minimum length of each arm |
| `maxLoopLength` | int | 50 | Maximum loop length between arms |
| `minLoopLength` | int | 3 | Minimum loop length |

### Output

Collection of `InvertedRepeatResult` containing:
- `LeftArmStart`: 0-based start position of left arm
- `RightArmStart`: 0-based start position of right arm
- `ArmLength`: Length of each arm
- `LoopLength`: Length of intervening sequence
- `LeftArm`: Left arm sequence
- `RightArm`: Right arm sequence (reverse complement of left)
- `Loop`: Intervening sequence
- `CanFormHairpin`: True if loop ≥ 3 (biologically viable)
- `TotalLength`: 2 × ArmLength + LoopLength

### Algorithm (Implementation)

```
1. For each position i from 0 to (length - 2×minArmLength - minLoopLength):
   a. For each armLength from minArmLength upward:
      i.   Extract leftArm of armLength at position i
      ii.  Compute leftArmRevComp = ReverseComplement(leftArm)
      iii. For each loop position j from (i + armLength + minLoopLength) to min(i + armLength + maxLoopLength, length - armLength):
           - Extract rightArm at position j
           - If rightArm == leftArmRevComp:
             • Calculate loopLength = j - (i + armLength)
             • If not already reported:
               → Yield InvertedRepeatResult
```

### Complexity
- **Time:** O(n² × L) where n = sequence length, L = maxLoopLength
- **Space:** O(k) where k = number of results (plus deduplication set)

### Comparison with EMBOSS einverted

| Feature | This Implementation | EMBOSS einverted |
|---------|---------------------|------------------|
| Algorithm | Exact matching | Dynamic programming |
| Mismatches | Not allowed | Allowed (with penalty) |
| Gaps | Not allowed | Allowed (with penalty) |
| Threshold | Arm length based | Score based |
| Overlapping | Not reported | Not reported |

---

## 4. Invariants

| Invariant | Formula | Verification |
|-----------|---------|--------------|
| Reverse complement match | ReverseComplement(LeftArm) = RightArm | String equality |
| Total length | TotalLength = 2 × ArmLength + LoopLength | Arithmetic |
| Position ordering | LeftArmStart < RightArmStart | Always true |
| Loop calculation | LoopLength = RightArmStart - (LeftArmStart + ArmLength) | Arithmetic |
| Hairpin viability | CanFormHairpin ⟺ LoopLength ≥ 3 | Boolean |

---

## 5. Edge Cases

| Case | Behavior | Rationale |
|------|----------|-----------|
| Empty sequence | Returns empty | No structure possible |
| Sequence shorter than 2×minArm + minLoop | Returns empty | Cannot form structure |
| No complementary regions | Returns empty | No inverted repeats exist |
| All same nucleotide (e.g., AAAA) | Returns empty | RevComp would be TTTT, no match |
| Self-complementary sequence (GCGC) | May find palindromic IR | GCGC revcomp = GCGC |
| Loop length = 0 | Not returned by default | minLoopLength enforces this |

---

## 6. Implementation Notes

### Current Implementation
- Located in `RepeatFinder.cs`
- Uses `DnaSequence.GetReverseComplementString()` for complement calculation
- Deduplicates using HashSet with (leftStart, rightStart, armLength) key
- Case-insensitive: input converted to uppercase

### RnaSecondaryStructure Alternative
- `RnaSecondaryStructure.FindInvertedRepeats()` provides similar functionality for RNA
- Uses RNA complement rules (U↔A instead of T↔A)
- Returns tuple-based results instead of structured record

---

## 7. References

1. Ussery DW, Wassenaar T, Borini S (2008). "Word Frequencies, Repeats, and Repeat-related Structures in Bacterial Genomes." Computing for Comparative Microbial Genomics.
2. Pearson CE, Zorbas H, Price GB, Zannis-Hadjopoulos M (1996). "Inverted repeats, stem-loops, and cruciforms: significance for initiation of DNA replication." J Cell Biochem 63(1):1-22.
3. Bissler JJ (1998). "DNA inverted repeats and human disease." Frontiers in Bioscience 3:d408-18.
4. Rice P, Longden I, Bleasby A (2000). "EMBOSS: the European Molecular Biology Open Software Suite." Trends Genet 16(6):276-7.
