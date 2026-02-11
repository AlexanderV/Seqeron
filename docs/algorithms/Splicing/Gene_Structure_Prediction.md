# Gene Structure Prediction

## Documented Theory

### Purpose

Gene structure prediction identifies the boundaries between exons (expressed regions)
and introns (intervening regions) in a pre-mRNA sequence, predicting which segments
will be retained in the mature mRNA after splicing (Gilbert, 1978; Breathnach & Chambon, 1981).

### Core Mechanism

1. **Splice site identification:** Scan for canonical donor (GT/GU) and acceptor (AG)
   dinucleotides using position weight matrices (PWMs) derived from consensus sequences
   (Shapiro & Senapathy, 1987).

2. **Intron prediction:** Pair each donor with each acceptor where:
   - acceptor position > donor position + minimum intron length
   - intron length ≤ maximum intron length
   - Combined score (donor + acceptor + branch point) / 3 ≥ minimum score threshold

3. **Non-overlapping selection:** Select highest-scoring introns greedily, rejecting
   any intron that overlaps with a previously selected one.

4. **Exon derivation:** Regions between (and flanking) selected introns are classified
   as exons if they meet minimum exon length. Exon types: Initial (before first intron),
   Internal (between introns), Terminal (after last intron), Single (no introns).

5. **Spliced sequence generation:** Concatenate exonic regions, removing intronic
   sequences, to produce the mature mRNA.

### Properties

- Deterministic for fixed input and parameters.
- Greedy non-overlapping selection is a heuristic (not optimal in the dynamic
  programming sense used by GenScan or Augustus).
- Exon phase tracks reading frame: phase = (sum of previous exon lengths) mod 3.

### Complexity

| Aspect | Value | Source |
|--------|-------|--------|
| Time | O(D × A) where D=donors, A=acceptors | Pairwise pairing |
| Space | O(n) for selected introns | Linear in sequence length |

---

## Implementation Notes

**Implementation location:** [SpliceSitePredictor.cs](src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs)

- `PredictGeneStructure(sequence, minExonLength=30, minIntronLength=60, minScore=0.5)`:
  Orchestrates intron prediction, non-overlapping selection, exon derivation, and
  spliced sequence generation. Returns a `GeneStructure` record.

- `PredictIntrons(sequence, minIntronLength=60, maxIntronLength=100000, minScore=0.5)`:
  Pairs donor and acceptor sites. Uses 0.8 × minScore for individual site thresholds.
  Finds branch points in [acceptor-50, acceptor-18]. Combined score = average of
  donor, acceptor, and branch point scores.

- `SelectNonOverlappingIntrons(introns)`: Private. Greedy selection by descending score.
  Tracks used positions via HashSet.

- `DeriveExons(sequence, introns, minExonLength)`: Private. Maps gaps between introns
  to exons with type classification and phase calculation.

- `GenerateSplicedSequence(sequence, introns)`: Private. Concatenates non-intronic segments.

---

## Deviations and Assumptions

1. **Greedy vs. DP selection:** The implementation uses greedy non-overlapping intron
   selection (highest score first), whereas formal gene finders like GenScan (Burge & Karlin, 1997)
   use generalized HMMs. This is a simplification.

2. **Intron score as average:** Combined intron score = (donor + acceptor + branchPoint) / 3.
   No formal standard defines this averaging; most tools use log-likelihood models.

3. **Default branch point score:** When no branch point is found, a default score of 0.3
   is used. This is implementation-specific.

4. **Exon phase starts at 0:** The first exon always has phase 0, assuming the reading
   frame begins at the start of the sequence. Real genes may have upstream UTR.

5. **Single exon on no introns:** When no introns pass the threshold, the entire sequence
   is reported as a single exon (ExonType.Single). This is correct for intronless genes
   but may be misleading for sequences where splice sites simply weren't detected.

---

## Sources

- Gilbert W (1978). Why genes in pieces? Nature 271:501.
- Breathnach R, Chambon P (1981). Organization and expression of eucaryotic split genes.
  Annu Rev Biochem 50:349-383.
- Shapiro MB, Senapathy P (1987). RNA splice junctions of different classes of eukaryotes.
  Nucleic Acids Res 15:7155-7174.
- Burge CB, Tuschl T, Sharp PA (1999). Splicing of precursors to mRNAs by the spliceosomes.
  In: The RNA World, 2nd ed. CSHL Press.
- Sakharkar MK et al. (2002). ExInt: an Exon Intron Database. Nucleic Acids Res 30:191-194.
- Alberts B et al. (2002). Molecular Biology of the Cell, 4th ed. Garland Science.
