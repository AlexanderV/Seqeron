# RNA-STEMLOOP-001 Evidence Document

**Test Unit:** RNA-STEMLOOP-001  
**Title:** Stem-Loop Detection  
**Area:** RnaStructure  
**Created:** 2026-02-05  
**Status:** Complete

---

## 1. Overview

This document provides evidence for testing the stem-loop (hairpin) detection algorithms in RNA sequences. Stem-loops are fundamental secondary structure motifs that form when a single-stranded RNA folds back on itself via Watson-Crick and wobble base pairing.

---

## 2. Authoritative Sources

### 2.1 Primary References

| Source | Type | Key Information | URL/Citation |
|--------|------|-----------------|--------------|
| Wikipedia - Stem-loop | Encyclopedia | Formation, stability, structural contexts | https://en.wikipedia.org/wiki/Stem-loop |
| Wikipedia - Tetraloop | Encyclopedia | GNRA, UNCG, CUUG special loops | https://en.wikipedia.org/wiki/Tetraloop |
| Wikipedia - Pseudoknot | Encyclopedia | Crossing base pairs, detection difficulty | https://en.wikipedia.org/wiki/Pseudoknot |
| Woese et al. (1990) | Peer-reviewed | GNRA/UNCG tetraloops make up 70% in 16S rRNA | PNAS 87(21):8467-8471 |
| Heus & Pardi (1991) | Peer-reviewed | GNRA tetraloop stability mechanisms | Science 253(5016):191-194 |
| Antao et al. (1991) | Peer-reviewed | UUCG is most stable tetraloop | Nucleic Acids Res 19(21):5901-5905 |
| Svoboda & Cara (2006) | Review | Hairpin RNA biology and importance | Cell Mol Life Sci 63(7):901-908 |
| Rivas & Eddy (1999) | Peer-reviewed | Pseudoknot detection algorithm | J Mol Biol 285(5):2053-2068 |

### 2.2 Key Evidence from Sources

**From Wikipedia - Stem-loop:**
> "Optimal loop length tends to be about 4-8 bases long; loops that are fewer than three bases long are sterically impossible and thus do not form"

**From Wikipedia - Tetraloop:**
> "UUCG and GNRA tetraloops make up 70% of all tetraloops in 16S-rRNA"
> "The UUCG tetraloop is the most stable tetraloop"

**From Wikipedia - Pseudoknot:**
> "A pseudoknot occurs when base pairs cross each other, i.e., for pairs (i,j) and (k,l), we have i < k < j < l"

---

## 3. Algorithm Specification

### 3.1 Stem-Loop (Hairpin) Definition

A stem-loop consists of:
1. **Stem:** Double-stranded helical region formed by antiparallel Watson-Crick/wobble base pairing
2. **Loop:** Single-stranded region at the end of the stem (minimum 3 nucleotides due to steric constraints)

**Structure:**
```
      Loop (≥3 nt)
    5'-G-A-A-A-3'
       |     |
  5'-G-C     G-C-3'
     | |     | |
     G-C     G-C
     | |     | |
     C-G     C-G
     Stem    Stem
```

### 3.2 Detection Algorithm

**Input:** RNA sequence, minStemLength, minLoopSize, maxLoopSize, allowWobble

**Approach:**
1. Scan sequence for potential loop positions
2. For each potential loop, extend stem on both sides
3. Check for valid base pairing (WC or wobble)
4. Return stem-loops meeting minimum requirements

**Complexity:** O(n² × L) where n = sequence length, L = max loop size

### 3.3 Special Cases

**Tetraloops (4-nucleotide loops):**
- GNRA (N=any, R=purine): GAAA, GCAA, GGAA, GUAA
- UNCG: UACG, UCCG, UGCG, UUCG
- CUUG: CUUG, CCUG

**Stability bonus:** ~3.0 kcal/mol for GNRA/UNCG tetraloops

### 3.4 Pseudoknot Detection

Pseudoknots occur when base pairs "cross" each other:
- For pairs (i,j) and (k,l): pseudoknot if i < k < j < l

**Note:** Standard secondary structure prediction cannot detect pseudoknots efficiently (NP-complete for general case).

---

## 4. Constraints and Invariants

### 4.1 Biological Constraints

| Constraint | Value | Source |
|------------|-------|--------|
| Minimum loop size | 3 nucleotides | Wikipedia (steric impossibility for <3) |
| Optimal loop size | 4-8 nucleotides | Wikipedia |
| Minimum stem length | 2-3 bp typically | Implementation choice |
| Valid base pairs | A-U, U-A, G-C, C-G, G-U, U-G | IUPAC, Wikipedia |

### 4.2 Algorithm Invariants

1. Stem-loops must not overlap (for simple structure)
2. Each base can participate in at most one base pair
3. Loop must be contiguous
4. Stem regions must be antiparallel (5'→3' pairs with 3'→5')

---

## 5. Test Methodology

### 5.1 Methods Under Test

| Method | Type | Priority |
|--------|------|----------|
| `FindStemLoops(sequence, minStem, minLoop, maxLoop, allowWobble)` | Canonical | Must |
| `FindHairpins(sequence, params)` | Variant | Should |
| `FindPseudoknots(sequence)` | Structural | Should |

### 5.2 Test Categories

**Must Tests:**
1. Basic hairpin detection (GGGAAAACCC)
2. Minimum loop size constraint (3 nt)
3. Wobble pair inclusion/exclusion
4. Empty/invalid input handling
5. Case insensitivity
6. Tetraloop detection
7. Dot-bracket notation generation

**Should Tests:**
1. Multiple stem-loops in sequence
2. Pseudoknot detection from base pairs
3. Energy calculation for stems
4. tRNA-like complex structures

---

## 6. Reference Test Data

### 6.1 Simple Hairpin

**Sequence:** `GGGAAAACCC`
**Expected:**
- Stem: 3 bp (GGG:CCC)
- Loop: 4 nt (AAAA)
- Structure: `(((...)))` 

**Source:** Basic Watson-Crick complementarity

### 6.2 GNRA Tetraloop

**Sequence:** `GGGGCGAACCCC` 
**Expected:**
- Stem: 4 bp
- Loop: CGAA (follows GN-RA pattern where first loop base is closing pair)
- Tetraloop bonus applicable

**Source:** Wikipedia Tetraloop, Heus & Pardi (1991)

### 6.3 No Structure Possible

**Sequence:** `AAAAAAAAAAAAAAA`
**Expected:** No stem-loops (no complementary bases)

### 6.4 Minimum Loop Size

**Sequence with <3 loop:** Cannot form (steric constraint)
**Sequence:** `GCAUC` → Too short

### 6.5 Pseudoknot Example

**Base pairs:** (0,6) and (3,9)
- Positions: 0 < 3 < 6 < 9 → Crossing detected
- Result: Pseudoknot identified

**Source:** Wikipedia Pseudoknot

---

## 7. Edge Cases

| Case | Input | Expected | Source |
|------|-------|----------|--------|
| Empty string | `""` | Empty result | Implementation |
| Null | `null` | Empty result or exception | Implementation |
| Too short | `"GC"` | No stem-loops | Constraint |
| No complement | `"AAAA"` | Empty | Biology |
| Lowercase | `"gggaaaaccc"` | Same as uppercase | Convention |
| Minimum stem | 2 bp | Depends on minStemLength param | Param |
| Maximum loop | >30 nt | Higher energy penalty | Turner 2004 |
| All-C loop | `"CCCC"` | Energy penalty | Turner 2004 |
| Wobble only | Stem with only G-U | Valid if allowWobble=true | Biology |

---

## 8. Known Limitations

1. **No pseudoknot prediction from sequence:** Only detection from known base pairs
2. **Simplified energy model:** Not as accurate as ViennaRNA
3. **No internal loops/bulges:** Focus is on hairpin loops
4. **Single structure:** Does not enumerate suboptimal structures

---

## 9. Implementation Notes

The implementation uses a scanning approach:
1. Iterate over potential loop start positions
2. For each loop size in [minLoopSize, maxLoopSize]
3. Extend stem by checking base pair validity
4. Collect results meeting minStemLength requirement

**Key design decisions:**
- Allow wobble pairs by default (biologically valid)
- Default minimum loop size = 3 (steric constraint)
- Default minimum stem length = 3 (stability consideration)

---

## 10. Sign-off Checklist

- [x] Primary sources identified and cited
- [x] Algorithm specification documented
- [x] Test data with expected values
- [x] Edge cases enumerated
- [x] Invariants defined
- [x] Limitations acknowledged
