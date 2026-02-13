# CpG Site Detection

## Documented Theory

### Purpose

CpG site detection identifies all positions in a DNA sequence where a cytosine nucleotide is immediately followed by a guanine nucleotide in the 5'→3' direction (CpG dinucleotides). CpG island detection identifies genomic regions with high CpG density that meet defined criteria for length, GC content, and CpG observed/expected ratio.

**Sources:** Gardiner-Garden & Frommer (1987), Wikipedia "CpG site" citing [13].

### Core Mechanism

**CpG site scanning:** Linear scan of a DNA sequence checking each consecutive pair of nucleotides for the pattern C followed by G. Returns 0-based positions of the C in each CpG dinucleotide.

**CpG O/E ratio** (Gardiner-Garden & Frommer, 1987):

$$\text{O/E} = \frac{\text{CpG count}}{\frac{C_{\text{count}} \times G_{\text{count}}}{L}}$$

Where:
- CpG count = number of CG dinucleotides in the sequence
- C_count = number of cytosines in the sequence
- G_count = number of guanines in the sequence
- L = total length of the sequence

**CpG island detection:** Sliding window approach evaluating regions against three criteria (Gardiner-Garden & Frommer, 1987):
1. Length ≥ 200 bp (default minimum)
2. GC content > 50%
3. CpG O/E ratio > 0.6

### Properties

- CpG site scanning is deterministic and exact — no heuristics or approximations.
- CpG O/E ratio is a well-defined mathematical formula with a single correct output for any given input sequence.
- CpG island detection uses a sliding window heuristic; boundary positions may vary with step size but the classification criteria are fixed and formal.

### Complexity

| Aspect | Value | Source |
|--------|-------|--------|
| Time — FindCpGSites | O(n) | Single pass over sequence |
| Time — CalculateCpGObservedExpected | O(n) | Single pass for counts + scan |
| Time — FindCpGIslands | O(n) | Sliding window with fixed step |
| Space | O(1) per site | Yields results incrementally |

---

## Implementation Notes

**Implementation location:** [EpigeneticsAnalyzer.cs](src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs)

- `FindCpGSites(sequence)`: Scans for 'C' followed by 'G' after converting to uppercase. Returns 0-based positions via `IEnumerable<int>`. Handles null/empty input by yielding nothing.
- `CalculateCpGObservedExpected(sequence)`: Computes CpG O/E ratio using the Gardiner-Garden & Frommer formula. Returns 0.0 for null/empty/short sequences or when expected = 0.
- `FindCpGIslands(sequence, minLength, minGc, minCpGRatio)`: Sliding window scan (step=1, per EMBOSS cpgplot canonical algorithm) identifying contiguous regions meeting all three CpG island criteria. Default parameters match Gardiner-Garden & Frommer (1987): minLength=200, minGc=0.5, minCpGRatio=0.6.

---

## Deviations and Assumptions

None. Implementation fully conforms to Gardiner-Garden & Frommer (1987) criteria and EMBOSS cpgplot canonical algorithm.

---

## Sources

- Gardiner-Garden M, Frommer M (1987). "CpG islands in vertebrate genomes." J Mol Biol. 196(2):261–282. https://doi.org/10.1016/0022-2836(87)90689-9
- Takai D, Jones PA (2002). "Comprehensive analysis of CpG islands in human chromosomes 21 and 22." Proc Natl Acad Sci USA. 99(6):3740–5. https://doi.org/10.1073/pnas.052410099
- Saxonov S, Berg P, Brutlag DL (2006). "A genome-wide analysis of CpG dinucleotides in the human genome." Proc Natl Acad Sci USA. 103(5):1412–1417. https://doi.org/10.1073/pnas.0510310103
- Wikipedia. "CpG site." https://en.wikipedia.org/wiki/CpG_site (accessed 2026-02-13).
