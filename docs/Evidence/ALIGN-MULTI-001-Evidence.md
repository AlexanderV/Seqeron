# ALIGN-MULTI-001: Multiple Sequence Alignment - Evidence Document

## Test Unit
**ID:** ALIGN-MULTI-001
**Area:** Alignment
**Methods:** `SequenceAligner.MultipleAlign(IEnumerable<DnaSequence>, ScoringMatrix?)`

---

## Authoritative Sources

### Source 1: Wikipedia - Multiple Sequence Alignment
**URL:** https://en.wikipedia.org/wiki/Multiple_sequence_alignment

**Key Information Extracted:**
1. **Definition:** MSA is the process of aligning three or more biological sequences (protein, DNA, or RNA)
2. **Objective:** Insert gaps into sequences until all sequences conform to length L ≥ max{n_i}
3. **Constraint:** No column may consist only of gaps
4. **Complexity:** Finding optimal MSA is NP-complete (Wang & Jiang, 1994; Just, 2001; Elias, 2006)
5. **Progressive Alignment:** Most widely used heuristic approach developed by Feng & Doolittle (1987)
   - Stage 1: Create guide tree representing sequence relationships
   - Stage 2: Build MSA by adding sequences according to guide tree
6. **Sum-of-Pairs Score:** Common objective function for MSA optimization
7. **Consensus Sequence:** Can be derived from aligned columns using majority voting
8. **Error Propagation:** Progressive methods propagate early alignment errors to final result
9. **Benchmark:** BAliBase evaluation found at least 24% of aligned amino acid pairs were incorrect (Nuin et al., 2006)

### Source 2: Wikipedia - Clustal
**URL:** https://en.wikipedia.org/wiki/Clustal (ClustalW redirect)

**Key Information Extracted:**
1. **Star Alignment:** Simple variant where one sequence serves as reference (center)
2. **Progressive Method:** Aligns sequences in most-to-least similarity order
3. **Guide Tree:** Can be generated via UPGMA or neighbor-joining
4. **Pairwise Alignment:** First step computes distance matrix between all pairs
5. **Gap Penalties:** Gap opening and gap extension penalties are configurable
6. **Output Symbols:**
   - `*` asterisk: fully conserved position
   - `:` colon: conservation between strongly similar groups
   - `.` period: conservation between weakly similar groups
   - ` ` blank: non-conserved

---

## Implementation Notes (from source code analysis)

The `SequenceAligner.MultipleAlign()` method implements a **simple star alignment** approach:

1. Uses **first sequence as reference** (star center)
2. Aligns all other sequences to the reference using `GlobalAlign()`
3. **Pads sequences** to equal length with gaps (`-`)
4. Generates **consensus sequence** via majority voting at each position

### Implementation Characteristics:
- **Complexity:** O(n² × m) where n = number of sequences, m = sequence length
- **No guide tree construction** (simplified star alignment)
- **Uses global alignment** (Needleman-Wunsch) for pairwise alignments
- **Returns:** `MultipleAlignmentResult` containing:
  - `AlignedSequences`: Array of aligned sequence strings
  - `Consensus`: Majority-voted consensus sequence
  - `TotalScore`: Sum of pairwise alignment scores

---

## Test Datasets from Sources

### From Wikipedia MSA Definition:
1. **Minimum sequences:** MSA requires 3+ sequences (pairwise alignment handles 2)
2. **Empty set:** Edge case returning empty result
3. **Single sequence:** Trivial case returning the sequence unchanged
4. **Identical sequences:** Should produce perfect alignment
5. **Variable length sequences:** Requires gap insertion

### From ClustalW/Progressive Alignment:
1. **Highly similar sequences:** Expected to align well
2. **Divergent sequences:** Progressive methods may produce suboptimal results
3. **Gap-heavy sequences:** Tests gap padding behavior

---

## Edge Cases and Corner Cases

| Case | Expected Behavior | Source |
|------|-------------------|--------|
| Empty input | Return `MultipleAlignmentResult.Empty` | Implementation contract |
| Single sequence | Return alignment containing just that sequence | Implementation contract |
| Null input | Throw `ArgumentNullException` | .NET conventions |
| Two identical sequences | Perfect alignment, consensus equals input | MSA definition |
| All identical sequences | All outputs equal input, consensus equals input | MSA definition |
| Sequences of different lengths | All aligned sequences equal length (padded) | Wikipedia MSA |
| All-gap column | Should not occur per MSA invariant | Wikipedia MSA |

---

## Invariants (from sources)

1. **Length invariant:** All aligned sequences have equal length
2. **Gap-only column prohibition:** No column consists only of gaps
3. **Reversibility:** Removing gaps from each aligned sequence recovers the original
4. **Count preservation:** Number of aligned sequences equals number of input sequences
5. **Score non-negativity for identical sequences:** Identical sequences should yield positive/zero score
6. **Consensus validity:** Consensus characters come from {A, C, G, T, -}

---

## Documented Failure Modes (from sources)

1. **Error propagation:** Early alignment errors propagate through progressive alignment
2. **Order dependence:** Result depends on alignment order (star alignment uses first sequence)
3. **Distantly related sequences:** Performance degrades with evolutionary distance
4. **Suboptimality:** Heuristic methods do not guarantee optimal solution

---

## Scoring Considerations

- **Sum-of-pairs score:** Total score is sum of all pairwise alignment scores
- Implementation calculates: sum of scores from aligning each sequence to reference

---

## References

1. Feng DF, Doolittle RF (1987). "Progressive sequence alignment as a prerequisite to correct phylogenetic trees." J Mol Evol. 25(4):351-360.
2. Wang L, Jiang T (1994). "On the complexity of multiple sequence alignment." J Comput Biol. 1(4):337-348.
3. Nuin PA, Wang Z, Tillier ER (2006). "The accuracy of several multiple sequence alignment programs for proteins." BMC Bioinformatics. 7:471.
4. Thompson JD, Higgins DG, Gibson TJ (1994). "CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment." Nucleic Acids Res. 22(22):4673-80.
