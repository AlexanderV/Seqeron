# Primer Structure Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Molecular Tools |
| Test Unit ID | PRIMER-STRUCT-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Primer structure analysis evaluates PCR primers for secondary-structure formation and self-complementarity issues that can reduce amplification efficiency. In this repository, the documented surface covers hairpin detection, primer-dimer detection, 3' end stability estimation, homopolymer detection, and dinucleotide-repeat detection. The implementation combines exact small-string heuristics with a suffix-tree-assisted branch for long sequences, and it exposes discrete boolean or scalar quality signals rather than a full thermodynamic folding model.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Hairpins form when a primer contains self-complementary regions separated by a loop; the original document identifies a minimum stem length of 4 bp and a minimum loop length of 3 nt as the practical threshold for the documented workflow. Primer-dimers arise when primers have complementary 3' ends that can be extended by polymerase, and high GC content at the 3' end increases their stability. The 3' terminal stability metric is documented as a nearest-neighbor $\Delta G$ calculation over the last 5 bases, following SantaLucia (1998) and the Primer3 `PRIMER_MAX_END_STABILITY` convention. Homopolymer runs and dinucleotide repeats are treated as primer-design liabilities because they can promote slippage, mispriming, and secondary structure. Sources: Wikipedia (Stem-loop, Primer dimer, Nucleic acid thermodynamics), SantaLucia (1998), Primer3 Manual.

### 2.2 Core Model

The repository models hairpins as complementary stems separated by at least `minLoopLength`, using a minimum stem length parameter to define candidate stems. Primer-dimer detection focuses on complementarity between the last bases of one primer and the reverse complement of the other primer. The 3' stability score is the sum of nearest-neighbor $\Delta G^\circ_{37}$ values for the four dinucleotide steps in the terminal 5-mer, plus terminal initiation parameters of `+0.98` kcal/mol for terminal G·C and `+1.03` kcal/mol for terminal A·T, as documented in the original file. The source comments identify `GCGCG` as the most stable 5-mer (`-6.86 kcal/mol`) and `TATAT` as the least stable (`-0.86 kcal/mol`) under the cited parameterization.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `HasHairpinPotential(...)` returns `false` when the sequence is shorter than `2 * minStemLength + minLoopLength` | The source guards on that minimum structure size |
| INV-02 | `Calculate3PrimeStability(...)` returns `0` for sequences shorter than 5 nt | The method explicitly short-circuits for short inputs |
| INV-03 | `FindLongestHomopolymer(...)` returns `0` for empty input and at least `1` for any non-empty input | The implementation tracks consecutive identical bases |
| INV-04 | `FindLongestDinucleotideRepeat(...)` returns `0` for inputs shorter than 4 nt | The implementation short-circuits when no repeated dinucleotide can exist |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[HasHairpinPotential] sequence` | `string` | required | Primer sequence to inspect for self-complementary stems | Case-insensitive in implementation |
| `[HasHairpinPotential] minStemLength` | `int` | `4` | Minimum complementary stem length | Used in both simple and suffix-tree branches |
| `[HasHairpinPotential] minLoopLength` | `int` | `3` | Minimum loop length between complementary stems | Enforced as a positional separation constraint |
| `[HasPrimerDimer] primer1, primer2` | `string` | required | Primer sequences to compare for 3' complementarity | Empty input returns `false` |
| `[HasPrimerDimer] minComplementarity` | `int` | `4` | Minimum complementary pairs to flag a dimer | Applied to the compared terminal window |
| `[Calculate3PrimeStability] sequence` | `string` | required | Primer sequence whose 3' 5-mer is scored | Inputs shorter than 5 return `0` |
| `[FindLongestHomopolymer/FindLongestDinucleotideRepeat] sequence` | `string` | required | Sequence to scan for runs and repeats | Case-insensitive in implementation |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `[HasHairpinPotential] hasHairpin` | `bool` | `true` when the implementation finds a valid self-complementary stem-loop candidate |
| `[HasPrimerDimer] hasPrimerDimer` | `bool` | `true` when terminal complementarity meets or exceeds the threshold |
| `[Calculate3PrimeStability] deltaG` | `double` | 3' end stability in kcal/mol; more negative values indicate greater stability |
| `[FindLongestHomopolymer] maxRun` | `int` | Length of the longest mononucleotide run |
| `[FindLongestDinucleotideRepeat] maxRepeatCount` | `int` | Longest repeated dinucleotide count |

### 3.3 Preconditions and Validation

All string-based methods normalize to uppercase before character comparisons. `HasHairpinPotential(...)` and `HasPrimerDimer(...)` return `false` for null or empty input. `Calculate3PrimeStability(...)` returns `0` for null, empty, or shorter-than-5 input. `FindLongestHomopolymer(...)` returns `0` for empty input, and `FindLongestDinucleotideRepeat(...)` returns `0` for inputs shorter than 4 nt.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the primer sequence to uppercase.
2. For hairpin detection, choose the simple or suffix-tree-assisted branch based on sequence length.
3. Search for complementary stems that satisfy the minimum stem and loop constraints.
4. For primer-dimer detection, reverse-complement the second primer and compare the terminal windows.
5. For 3' stability, evaluate the final 5-mer with the nearest-neighbor table and terminal initiation terms.
6. For run and repeat metrics, scan the sequence for the longest mononucleotide or dinucleotide repetition.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The source uses two hairpin-detection strategies:

| Sequence Length | Strategy | Documented Rationale |
|-----------------|----------|----------------------|
| `< 100 bp` | Nested-loop search | Lower overhead for typical PCR primer lengths |
| `>= 100 bp` | Suffix-tree-assisted search | Avoids the short-sequence quadratic scan for longer inputs |

Nearest-neighbor $\Delta G$ values documented for the 3' stability calculation:

| Dinucleotide | ΔG (kcal/mol) |
|--------------|---------------|
| AA/TT | -1.00 |
| AT | -0.88 |
| TA | -0.58 |
| CA/TG | -1.45 |
| GT/AC | -1.44 |
| CT/AG | -1.28 |
| GA/TC | -1.30 |
| CG | -2.17 |
| GC | -2.24 |
| GG/CC | -1.84 |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `HasHairpinPotential` | `O(n²)` for `< 100 bp`; suffix-tree-assisted for `>= 100 bp` | `O(1)` auxiliary for the simple branch | The long-sequence branch builds and queries a suffix tree |
| `HasPrimerDimer` | `O(n)` | `O(n)` | Uses the reverse complement of the second primer and compares a terminal window |
| `Calculate3PrimeStability` | `O(1)` | `O(1)` | Only the last 5 bases are evaluated |
| `FindLongestHomopolymer` | `O(n)` | `O(1)` | Single left-to-right scan |
| `FindLongestDinucleotideRepeat` | `O(n)` | `O(1)` | Repeated dinucleotide scan as documented in the original file |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PrimerDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs)

- `PrimerDesigner.HasHairpinPotential(string, int, int)`: Detects self-complementary stem-loop candidates and switches between simple and suffix-tree-assisted branches.
- `PrimerDesigner.HasPrimerDimer(string, string, int)`: Checks 3' complementarity between primers.
- `PrimerDesigner.Calculate3PrimeStability(string)`: Computes the last-5-base nearest-neighbor stability with initiation terms.
- `PrimerDesigner.FindLongestHomopolymer(string)`: Returns the longest mononucleotide run.
- `PrimerDesigner.FindLongestDinucleotideRepeat(string)`: Returns the longest repeated dinucleotide count.

### 5.2 Current Behavior

The current implementation uses a simple nested-loop hairpin search below 100 bp and a suffix-tree-assisted branch at 100 bp or above. Primer-dimer detection compares the last `min(8, len1, len2)` bases of the first primer to the start of the reverse complement of the second primer. The 3' stability calculation evaluates the final 5-mer only, sums the four dinucleotide contributions, and adds terminal initiation parameters exactly as described in the original document and in the source comments. Homopolymer and dinucleotide-repeat detection are direct sequence scans without thermodynamic weighting.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Hairpin detection based on self-complementary stems separated by a minimum loop.
- Primer-dimer detection focused on 3' end complementarity.
- 3' terminal stability scoring from nearest-neighbor $\Delta G$ values with terminal initiation parameters.

**Intentionally simplified:**

- Hairpin detection is a boolean structural screen rather than a full free-energy folding model; **consequence:** the method reports potential stem-loops without ranking complete secondary-structure ensembles.
- Primer-dimer detection inspects terminal complementarity instead of a full duplex thermodynamic landscape; **consequence:** non-terminal or context-dependent dimer interactions are not separately modeled.
- Homopolymer and dinucleotide-repeat checks use run-length heuristics; **consequence:** they flag simple repetitive structure without estimating PCR yield impact.

**Not implemented:**

- Full secondary-structure thermodynamics for primer hairpins and dimers; **users should rely on:** no current alternative documented in this test unit.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Sequence too short for a hairpin | Returns `false` | A valid stem-loop cannot satisfy the minimum stem and loop constraints |
| No self-complementary regions | Returns `false` | No candidate stem-loop is found |
| Perfect palindrome | May or may not form a hairpin | Loop feasibility still matters |
| Empty primer in dimer check | Returns `false` | No terminal complementarity can be evaluated |
| Empty sequence in homopolymer detection | Returns `0` | No run exists |
| All bases unique | Homopolymer result is `1` | The longest run is a single base |
| All bases identical | Homopolymer result is sequence length | Every position extends the same run |
| Sequence shorter than 4 in dinucleotide-repeat detection | Returns `0` | No repeated dinucleotide can exist |

### 6.2 Limitations

The documented workflow is a screening-oriented implementation. It does not model full RNA/DNA folding thermodynamics, full primer-pair interaction landscapes, or polymerase- and buffer-specific effects. For typical PCR primers in the 18-25 bp range, the short-sequence branch is the intended path and the richer suffix-tree branch is primarily a long-sequence optimization.

## 8. References

1. Wikipedia - Primer (molecular biology): PCR primer design section.
2. Wikipedia - Primer dimer: Mechanism of formation.
3. Wikipedia - Stem-loop (Hairpin loop): Formation and stability.
4. Wikipedia - Nucleic acid thermodynamics: Nearest-neighbor method.
5. SantaLucia JR (1998) "A unified view of polymer, dumbbell and oligonucleotide DNA nearest-neighbor thermodynamics", PNAS 95:1460-65.
6. Primer3 Manual (primer3.org): `PRIMER_MAX_HAIRPIN_TH`, `PRIMER_MAX_SELF_END`, `PRIMER_MAX_END_STABILITY`, `PRIMER_MAX_POLY_X`.
