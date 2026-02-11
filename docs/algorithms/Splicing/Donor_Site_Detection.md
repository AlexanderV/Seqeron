# Donor (5') Splice Site Detection

## Documented Theory

### Purpose

The donor splice site (also called the 5' splice site) is the boundary at the 5' end
of an intron in pre-mRNA. Detection of donor splice sites is fundamental to gene
structure prediction, RNA splicing analysis, and understanding alternative splicing
events. The canonical donor site contains the almost invariant GU (GT in DNA)
dinucleotide at the start of the intron (Shapiro & Senapathy, 1987; Burge et al., 1999).

### Core Mechanism

Donor splice site detection uses a position weight matrix (PWM) model to score
candidate sites. The PWM encodes nucleotide frequency distributions at positions
-3 to +6 relative to the exon-intron junction, compiled from large datasets of
verified splice sites (Shapiro & Senapathy, 1987).

**Extended consensus:** MAG|GURAGU

Where:
- M = A or C (position -3)
- A = adenine (position -2, ~60% conserved)
- G = guanine (position -1, ~80% conserved)
- | = exon-intron boundary
- G = guanine (position 0, ~100% invariant)
- U = uracil (position +1, ~100% invariant)
- R = purine A/G (position +2, ~60% A)
- A = adenine (position +3, ~70%)
- G = guanine (position +4, ~80%)
- U = uracil (position +5, moderately conserved)

**Scoring:** Log-odds score at each position = log₂(observed_frequency / background_frequency),
where background is uniform 0.25 for each nucleotide (Yeo & Burge, 2004).

**Non-canonical donors:**
- GC-AG introns: ~0.5–1% of U2-type introns use GC at the donor site (Burge et al., 1999).
- U12-type AT-AC introns: ~0.3% of introns use AT donor / AC acceptor, processed by the minor spliceosome.

### Properties

- **Deterministic:** Given the same sequence and parameters, results are identical.
- **Position-invariant:** Each candidate site scored independently.
- **Sensitivity/specificity tradeoff:** Controlled by the `minScore` threshold parameter.

### Complexity

| Aspect | Value | Source |
|--------|-------|--------|
| Time | O(n) where n = sequence length | Single-pass scan |
| Space | O(1) per site (streaming via IEnumerable) | Implementation |

---

## Implementation Notes

**Implementation location:** [SpliceSitePredictor.cs](src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs)

- `FindDonorSites(sequence, minScore, includeNonCanonical)`: Scans sequence for GU dinucleotides
  (and optionally GC and AU), scores each candidate using the donor PWM, yields sites
  above `minScore`. Converts T→U and uppercase internally.
- `ScoreDonorSite(sequence, position)`: Private. Computes log-odds PWM score at positions
  -3 to +5 relative to the candidate GU, normalizes to [0, 1] range.
- `ScoreU12DonorSite(sequence, position)`: Private. Scores U12-type AT donors by matching
  against /AUAUCC/ consensus, returns fraction of matching positions.

---

## Deviations and Assumptions

1. **Score normalization:** The implementation normalizes the raw log-odds sum to [0, 1]
   using `(score/count + 2) / 4`, an approximate linear mapping. This differs from
   MaxEntScan's exact scoring but preserves relative ordering.

2. **GC donor penalty:** Non-canonical GC donors receive a 0.7 multiplier on their score.
   This is an implementation heuristic; no specific penalty factor was found in
   authoritative sources.

3. **Minimum sequence length:** Implementation requires sequence length ≥ 6 to produce
   any results, since the PWM window extends from position -3 to +5.

4. **T/U equivalence:** Input T is converted to U internally, so both DNA and RNA
   representations are accepted.

---

## Sources

- Shapiro MB, Senapathy P (1987). RNA splice junctions of different classes of eukaryotes. Nucleic Acids Res 15(17):7155-7174.
- Burge CB, Tuschl T, Sharp PA (1999). Splicing precursors to mRNAs by the spliceosomes. The RNA World, 2nd ed. CSHL Press, pp. 525-560.
- Yeo G, Burge CB (2004). Maximum entropy modeling of short sequence motifs. J Comput Biol 11(2-3):377-394.
- Wikipedia: RNA splicing. https://en.wikipedia.org/wiki/RNA_splicing
- Wikipedia: Spliceosome. https://en.wikipedia.org/wiki/Spliceosome
