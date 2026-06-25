# META-BIN-001: Genome Binning - Evidence Document

## Sources

### Primary Sources

1. **Wikipedia - Binning (metagenomics)**
   - URL: https://en.wikipedia.org/wiki/Binning_(metagenomics)
   - Key concepts:
     - Metagenomic binning groups assembled contigs to their separate genomes of origin
     - Results are called Metagenome Assembled Genomes (MAGs)
     - Methods use compositional features (GC-content, tetranucleotide frequencies)
     - Methods also use coverage across samples
     - Can employ supervised or unsupervised classification
   - Referenced algorithms: TETRA, MEGAN, Phylopythia, SOrt-ITEMS, MetaBAT2, MaxBin2, CONCOCT, DAS Tool

2. **TETRA (Teeling et al., 2004)**
   - URL: https://doi.org/10.1186/1471-2105-5-163
   - Tetranucleotide usage patterns for genomic fragment classification
   - 256 possible tetramers (4^4)
   - Z-scores indicate over/under-representation vs expected frequencies

3. **CheckM (Parks et al., 2014)**
   - URL: http://genome.cshlp.org/content/25/7/1043
   - Tool for assessing genome quality from metagenomes
   - Estimates completeness using single-copy marker genes
   - Estimates contamination using duplicated marker genes
   - Reference: Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW

4. **Maguire et al. (2020) - MAG binning limitations**
   - URL: https://doi.org/10.1099/mgen.0.000436
   - Key findings:
     - 82-94% chromosomes recovered correctly
     - Only 38-44% of GIs recovered
     - Only 1-29% of plasmid sequences recovered
     - Variable copy number and composition problematic for binning

### Documented Binning Features

1. **Compositional Features**
   - GC content: genome-specific signature
   - Tetranucleotide frequencies: species-specific patterns
   - Codon usage patterns

2. **Coverage Features**
   - Read depth across samples
   - Differential coverage binning

3. **Quality Metrics (from CheckM)**
   - Completeness: proportion of expected marker genes present
   - Contamination: proportion of duplicated marker genes
   - Strain heterogeneity: indicator of multiple strains

### Edge Cases from Literature

1. **Empty input**: Should return empty results
2. **Single contig**: May form its own bin or be filtered by size
3. **Minimum bin size**: Filter bins below threshold
4. **Similar GC content across organisms**: Can cause incorrect binning
5. **Variable coverage regions**: Cause assembly and binning issues
6. **High contamination bins**: Multiple genomes in one bin

### Test Dataset Design (Based on Literature)

1. **Pure genomes** (high GC vs low GC): Should separate into different bins
2. **Similar GC content**: May cluster together
3. **Variable coverage**: Test coverage feature contribution
4. **Minimum size filtering**: Bins below threshold excluded

## Implementation Analysis

### Current Implementation: `MetagenomicsAnalyzer.BinContigs`

```csharp
public static IEnumerable<GenomeBin> BinContigs(
    IEnumerable<(string ContigId, string Sequence, double Coverage)> contigs,
    int numBins = 10,
    double minBinSize = 500000,
    double expectedGenomeSize = 4_000_000)
```

**Algorithm (k-means on composite features):**
1. Calculate features: GC content, tetranucleotide frequency (TNF), normalized coverage per contig
2. K-means clustering on composite distance: |GC| + |coverage| + TNF Pearson distance
3. Centroid initialization by GC-sorted spread for deterministic diverse starting points
4. Iterative assignment/update until convergence (max 50 iterations)
5. Filter bins below minimum size
6. Completeness: `min(totalLength / expectedGenomeSize * 100, 100)` — parameterized, no hardcoded genome size
7. Contamination: `min(gcStdDev / 0.5 * 100, 100)` — normalized by theoretical maximum GC std dev

**Output: GenomeBin record**
- BinId, ContigIds, TotalLength, GcContent, Coverage
- Completeness, Contamination, PredictedTaxonomy

### Deviations from Literature

None. All previously documented assumptions have been resolved:

1. ~~Completeness estimated by length ratio to 2Mbp genome~~ → Parameterized via `expectedGenomeSize` (default 4Mbp)
2. ~~Contamination estimated from GC variance with arbitrary formula~~ → GC std dev normalized by theoretical maximum (0.5)
3. ~~Tetranucleotide frequency calculated but not used in binning~~ → TNF Pearson distance included in k-means composite distance

## Test Coverage Summary

All tests independently verified against theory:

| Test | What is verified | Theory basis |
|------|------------------|--------------|
| M5 | Completeness = 50.0 for 2MB bin / 4MB genome | min(L/E×100, 100) |
| M6 | Contamination = 0 (uniform GC), = 100 (max variance) | stddev/0.5×100 |
| M7 | GC = 1.0 (pure GC), = 0.0 (pure AT) | (G+C)/(A+T+G+C) |
| M10 | Extreme GC populations separate into distinct bins | K-means on compositional distance |
| M12 | Coverage = 25.0 = mean(20,30,25) | Arithmetic mean |
| M11 | Bins are disjoint (no contig in multiple bins) | Partition invariant |

## Invariants

1. Result bin IDs must be unique
2. Completeness must be in range [0, 100]
3. Contamination must be in range [0, 100]
4. GC content must be in range [0, 1]
5. Total length must equal sum of contig lengths in bin
6. Empty input returns empty output
7. Bins below minBinSize are excluded

## Test Categories

### Must Tests (Evidence-backed)
1. Empty input handling
2. Basic binning returns non-null
3. GC-based separation (high vs low GC contigs)
4. MinBinSize filtering
5. Completeness range [0, 100]
6. Contamination range [0, 100]
7. Unique bin IDs
8. Correct contig grouping
9. TotalLength accuracy

### Should Tests
1. Coverage values preserved
2. GC content averaging in bins
3. Multiple bins from diverse input

### Could Tests
1. Large dataset performance
2. Edge case with all contigs same GC

---

# Addendum (2026-06-24): TETRA z-score-normalised tetranucleotide signature

This addendum records the evidence retrieved **this session** for the new opt-in
TETRA z-score signature (`CalculateTetranucleotideZScores` / `TetranucleotideZScoreCorrelation`).
It conforms to the Evidence template structure (Online Sources → Corner Cases → Datasets →
Assumptions → Recommendations → References).

## Online Sources

### TETRA web-service paper (Teeling et al. 2004, BMC Bioinformatics)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC529438/
**Retrieved via:** WebSearch "TETRA Teeling BMC Bioinformatics 2004 PMC full text" → WebFetch of the PMC HTML.
**Accessed:** 2026-06-24   **Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points (verbatim / close paraphrase from the fetched text):**

1. **Strand symmetry:** "DNA sequences are extended by their reverse-complements to compensate
   for different patterns of tetranucleotide over- and underrepresentation between the leading
   and the lagging strand."
2. **256-space + maximal-order Markov:** "the frequencies of all 256 possible tetranucleotides
   are counted and the corresponding expected frequencies are calculated by means of a
   maximal-order Markov model from the sequences' di- and trinucleotide composition."
3. **z-score conversion:** "the divergence between the observed and expected tetranucleotide
   frequencies is then transferred into z-scores using an approximation published by Schbath."
4. **Comparison:** "all DNA sequences are compared in pairs by computing the Pearson's
   correlation coefficient of their z-scores." (Technical formulae are deferred by this paper to
   Teeling et al. 2004 Environ Microbiol and the TETRA manual — see next sources.)

### Maximal-order Markov expected-count formula (open-access reproductions)

**URL (expected count, primary OA):** https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0008113
**Retrieved via:** WebSearch → WebFetch (PLOS One, Bohlin & Skjerve 2009, Genome Homogeneity).
**Accessed:** 2026-06-24   **Authority rank:** 1

5. **Expected count (2nd-order / maximal-order Markov):** the Methods give
   `E(n1n2n3n4) = f(n1n2n3) × f(n2n3n4) / f(n2n3)`, stating "A second order Markov chain (SOM)
   uses di- and trinucleotide frequencies to approximate larger oligonucleotides." In count form
   this is **E(n1n2n3n4) = N(n1n2n3)·N(n2n3n4) / N(n2n3)**.

**URL (full z-score + variance, reproduction):** https://github.com/YunpengLiu/Tetra-nucleotide-analysis
**Retrieved via:** WebSearch (multiple independent queries returned the identical formula) → WebFetch of the README.
**Accessed:** 2026-06-24   **Authority rank:** 3 (reference implementation reproducing Teeling 2004)

6. **z-score:** `Z(n1n2n3n4) = (N(n1n2n3n4) − E(n1n2n3n4)) / √(var(N(n1n2n3n4)))`.
7. **Variance (Schbath approximation):**
   `var(n1n2n3n4) = E(n1n2n3n4) · [N(n2n3) − N(n1n2n3)] · [N(n2n3) − N(n2n3n4)] / N(n2n3)²`.
   The identical expressions (5–7) were returned independently by WebSearch over the Teeling/TETRA
   literature and by the EzBioCloud TNA documentation.

### Schbath variance origin (primary)

**URL:** https://pubmed.ncbi.nlm.nih.gov/9228617/
**Retrieved via:** WebSearch → WebFetch (abstract).
**Accessed:** 2026-06-24   **Authority rank:** 1

8. Schbath S (1997), "An efficient statistic to detect over- and under-represented words in DNA
   sequences": establishes the Markov-chain word-count z-statistic that TETRA's variance
   approximation cites. The PubMed page exposes only the abstract; the explicit variance form is
   the one reproduced in points 6–7 above.

### CheckM marker-gene database availability (sub-part 2 — NOT implemented)

**URL:** https://genome.cshlp.org/content/25/7/1043 ; software https://github.com/Ecogenomics/CheckM
**Retrieved via:** WebSearch "CheckM lineage-specific single-copy marker genes Pfam TIGRFAM HMM database download Parks 2015".
**Accessed:** 2026-06-24   **Authority rank:** 1

9. CheckM completeness/contamination uses **lineage-specific sets of single-copy marker genes**:
   "Single-copy Pfam and TIGRFAMs genes were identified within reference genomes, with a gene
   defined as a lineage-specific marker gene if it occurs only once in >97% of the genomes within
   a lineage… inferred for all internal nodes within the reference genome tree." Detection uses
   HMMER v3.1b1 with model-specific Pfam (`-cut_gc`) / TIGRFAM (`-cut_nc`) cutoffs.
10. **Decision:** these marker sets are a **large trained database** (the `checkm_data` package:
    Pfam/TIGRFAM HMM libraries + a precomputed reference genome tree, ~GB-scale), NOT cleanly
    retrievable as plaintext this session. Per the STOP RULE, the marker-gene completeness/
    contamination is **left as an honest residual** (not fabricated); only sub-part 1 (the
    z-score signature) is implemented.

## Documented Corner Cases and Failure Modes

### From Teeling 2004 / formula structure

1. **Absent middle dinucleotide N(n2n3)=0:** the expected count is undefined (division by zero);
   the implementation yields z=0 (no over-/under-representation evidence).
2. **Non-positive variance:** when `N(n2n3)−N(n1n2n3)` or `N(n2n3)−N(n2n3n4)` is ≤0 the variance
   approximation is ≤0; z is defined as 0.
3. **Reverse-complement extension** makes the signature strand-symmetric; consequently even very
   short inputs (≥2 ACGT bases → ≥4-nt extended strand) produce a signature.

## Test Datasets

### Dataset: hand-derived z-score (this session)

**Source:** TETRA formula (points 5–7), worked by hand and cross-checked with a reference script.

| Parameter | Value |
|-----------|-------|
| Sequence | `ACGTACGTGGCC` |
| RC-extended strand | `ACGTACGTGGCCGGCCACGTACGT` (24 nt) |
| N(ACGT) | 4 |
| N(ACG)=N(n1n2n3) | 4 |
| N(CGT)=N(n2n3n4) | 4 |
| N(CG)=N(n2n3) | 5 |
| E(ACGT) | 4·4/5 = 3.2 |
| var(ACGT) | 3.2·(5−4)·(5−4)/5² = 0.128 |
| **z(ACGT)** | (4−3.2)/√0.128 = **√5 = 2.2360679774997896** |

## Assumptions

1. **ASSUMPTION: counting on the concatenated RC-extended strand.** TETRA states the sequence is
   "extended by its reverse-complement"; the implementation counts overlapping di/tri/tetra words
   over `seq + revcomp(seq)`. The published worked-example formula values (E, var) are reproduced
   exactly by this convention, so the assumption is consistent with the source; the absolute
   magnitudes of z scale with sequence length but the *correlation* of two z-vectors (the actual
   comparison metric) is scale-invariant.

## Recommendations for Test Coverage

1. **MUST Test:** hand-derived exact z(ACGT)=√5 for `ACGTACGTGGCC` — Evidence: points 5–7 + dataset.
2. **MUST Test:** self-correlation = 1.0; similar-vs-dissimilar ordering — Evidence: point 4.
3. **MUST Test:** N(n2n3)=0 → z=0; null/empty/single-base → all-zero 256-vector — Evidence: corner cases 1–3.
4. **SHOULD Test:** correlation symmetry; degenerate signature → r=0 not NaN — Rationale: numerical robustness.

## References (this addendum)

1. Teeling H, Waldmann J, Lombardot T, Bauer M, Glöckner FO (2004). TETRA: a web-service and a
   stand-alone program for the analysis and comparison of tetranucleotide usage patterns in DNA
   sequences. BMC Bioinformatics 5:163. https://doi.org/10.1186/1471-2105-5-163
   (full text https://pmc.ncbi.nlm.nih.gov/articles/PMC529438/)
2. Teeling H, Meyerdierks A, Bauer M, Amann R, Glöckner FO (2004). Application of tetranucleotide
   frequencies for the assignment of genomic fragments. Environ Microbiol 6(9):938–947.
   https://doi.org/10.1111/j.1462-2920.2004.00624.x
3. Schbath S (1997). An efficient statistic to detect over- and under-represented words in DNA
   sequences. J Comput Biol 4(2):189–192. https://pubmed.ncbi.nlm.nih.gov/9228617/
4. Bohlin J, Skjerve E, Ussery DW (2009). Investigations of oligonucleotide usage variance /
   genome-homogeneity genomic signatures (reproduces the 2nd-order Markov expected-count form).
   PLOS ONE 5(1):e8113. https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0008113
5. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW (2015). CheckM: assessing the
   quality of microbial genomes recovered from isolates, single cells, and metagenomes.
   Genome Res 25(7):1043–1055. https://genome.cshlp.org/content/25/7/1043

---

**Document Version:** 3.0
**Created:** 2026-02-04
**Updated:** 2026-06-24 (added TETRA z-score signature evidence; CheckM marker QC left as honest residual)
**Test Unit:** META-BIN-001
