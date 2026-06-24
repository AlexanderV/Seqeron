# CHROM-CENT-001 Evidence: Centromere Analysis

## Test Unit
- **ID:** CHROM-CENT-001
- **Area:** Chromosome
- **Canonical Method:** `ChromosomeAnalyzer.AnalyzeCentromere(...)`

## Sources Consulted

### Primary Sources

1. **Wikipedia - Centromere**
   - URL: https://en.wikipedia.org/wiki/Centromere
   - Retrieved: 2026-01-31
   - Key information:
     - Centromere links sister chromatids during cell division
     - Creates short arm (p) and long arm (q)
     - Position classifications: Metacentric, Submetacentric, Subtelocentric, Acrocentric, Telocentric
     - Regional centromeres contain large arrays of repetitive DNA (alpha-satellite)
     - Human centromeric repeat unit is called α-satellite (alphoid), ~171 bp monomer
     - "Telocentric chromosomes... are not present in humans but can form through cellular chromosomal errors"

2. **Wikipedia - Karyotype**
   - URL: https://en.wikipedia.org/wiki/Karyotype
   - Retrieved: 2026-01-31
   - Key information:
     - Human chromosome groups based on centromere position
     - Group A (1-3): Large, metacentric or submetacentric
     - Group D (13-15): Medium-sized, acrocentric
     - Group G (21-22, Y): Very small, acrocentric

3. **Wikipedia - Chromosome**
   - URL: https://en.wikipedia.org/wiki/Chromosome
   - Retrieved: 2026-01-31
   - Key information:
     - Centromere region has constitutive heterochromatin with repetitive sequences
     - p arm (short) named from French "petit"
     - q arm (long) follows p in Latin alphabet

### Levan (1964) Classification — Arm Ratio Thresholds

Source: Levan A, Fredga K, Sandberg AA (1964). "Nomenclature for centromeric position on chromosomes". Hereditas. 52(2):201-220.

As cited on Wikipedia (Centromere article), the classification uses the arm ratio q/p (long arm / short arm):

| Centromere Position    | Arms length ratio (q/p) | Sign | Classification   |
|------------------------|-------------------------|------|------------------|
| Medial sensu stricto   | 1.0 – 1.6              | M    | Metacentric      |
| Medial region          | 1.7                     | m    | Metacentric      |
| Submedial              | 3.0                     | sm   | Submetacentric   |
| Subterminal            | 3.1 – 6.9              | st   | Subtelocentric   |
| Terminal region        | 7.0                     | t    | Acrocentric      |
| Terminal sensu stricto | ∞                       | T    | Telocentric      |

Implementation boundaries (derived from Levan table):
- ratio ≤ 1.7 → Metacentric (M + m)
- ratio (1.7, 3.0] → Submetacentric (sm)
- ratio (3.0, 7.0) → Subtelocentric (st)
- ratio ≥ 7.0 → Acrocentric (t)
- p = 0 → Telocentric (T)

**Boundary values:** 1.7 → Metacentric, 3.0 → Submetacentric, 7.0 → Acrocentric (per Levan table entries).

### Human Chromosome Centromere Positions (from Wikipedia)

| Chromosome | Centromere Position (Mb) | Type |
|------------|-------------------------|------|
| 1          | 125.0                   | Metacentric |
| 2          | 93.3                    | Submetacentric |
| 3          | 91.0                    | Metacentric |
| 13         | 17.9                    | Acrocentric |
| 21         | 13.2                    | Acrocentric |
| Y          | 12.5                    | Acrocentric |

Note: Practical human karyotype classifications are based on cytogenetic (microscopic) observation, not genomic coordinate ratios. The Levan thresholds applied to sequence-derived positions may yield slightly different classifications for borderline chromosomes.

### Alpha-Satellite DNA (from Wikipedia)

- Primary centromeric repeat in humans
- ~171 bp monomer repeat
- Forms higher-order repeat arrays
- Associated with heterochromatin

## Testing Methodology

Based on the literature:

1. **Arm-ratio-based classification test**: Verify classification thresholds match Levan (1964) nomenclature using arm ratio (q/p)
2. **Alpha-satellite detection**: Verify recognition of repetitive patterns characteristic of centromeric regions
3. **Boundary conditions**: Empty, short, and edge-case sequences
4. **Invariants**:
   - Start ≤ End when centromere found
   - Length = End - Start
   - Type is one of: Metacentric, Submetacentric, Subtelocentric, Acrocentric, Telocentric, Unknown

## Implementation Notes

The implementation uses:
- Sliding window approach with k-mer frequency analysis
- GC content variability as discriminating feature (centromeres have low GC variability)
- Repeat content estimation using k-mer counting (k=15)
- Arm ratio classification matching Levan (1964) nomenclature exactly

## Alpha-Satellite-Specific Detection (added 2026-06-24)

`AnalyzeCentromere` is a **generic tandem-repeat-density** heuristic — its `AlphaSatelliteContent`
is a repeat score, not an alpha-satellite-specific measurement. The methods
`DetectAlphaSatellite` and `FindCenpBBoxes` add **alpha-satellite-specific** detection based on
the two defining molecular signatures of human alphoid DNA, sourced below. Defaults of
`AnalyzeCentromere` and the meaning of `AlphaSatelliteContent` are unchanged (additive).

### Sources actually retrieved this session (2026-06-24)

1. **Hartley G, O'Neill RJ (2019). "Alpha satellite DNA biology: Finding function in the recesses of
   the genome." Genes (Basel).** — retrieved via PMC.
   - URL fetched: https://pmc.ncbi.nlm.nih.gov/articles/PMC6121732/
   - Verbatim: *"Alpha satellite DNA is composed of fundamental 171bp monomeric repeat units."*
   - Verbatim: the CENP-B box is *"a 17-bp sequence motif (5'-T/CTCGTTGGAAA/GCGGGA-3')"*; *"the CENP-B
     box is present in only a subset of alpha satellite monomers"* (e.g. D7Z1 has an alternating /
     every-other-monomer pattern; pentameric arrays irregular).
   - Verbatim: *"The individual monomers within a HOR unit have 50–70% identity …"*; monomers *"differ
     in sequence by 10–40%."*

2. **Masumoto H, Masukata H, Muro Y, Nozaki N, Okazaki T (1989). "A human centromere antigen
   (CENP-B) interacts with a short specific sequence in alphoid DNA, a human centromeric satellite."
   J Cell Biol 109(4):1963-1973.** — the primary CENP-B box source; canonical 17-bp consensus
   confirmed via the secondary retrieval below.
   - DOI/stable: https://doi.org/10.1083/jcb.109.4.1963 ; PubMed record retrieved at
     https://pubmed.ncbi.nlm.nih.gov/1730770/ (Masumoto 1992 follow-up; confirms "the 17-bp sequence,
     designated previously as CENP-B box").

3. **Centromere-formation / CENP-B-box review (PMC4843215, "CENP-B box … occurs in a New World
   monkey").** — retrieved via PMC; used to confirm the exact canonical 17-bp consensus string.
   - URL fetched: https://pmc.ncbi.nlm.nih.gov/articles/PMC4843215/
   - Verbatim canonical consensus (from Masumoto et al. 1989): **`YTTCGTTGGAARCGGGA`** (Y = C/T, R = A/G);
     broader/looser definition `NTTCGNNNNANNCGGGN` retains the TTCG and CGGG core recognition elements.

4. **Willard HF (1985); Waye JS, Willard HF (1987)** — original definition of the 171-bp alpha-satellite
   monomer and higher-order repeat (HOR) organization; referenced and quoted via PMC6121732 (source 1)
   and a literature search (Nature JHG; Genome Research 26:1301). Monomers are 50–70% identical and
   arranged into HOR units that are tandemly amplified.

### Derived detection parameters (all source-backed)

| Parameter | Value | Source |
|-----------|-------|--------|
| Alpha-satellite monomer length / tandem period | 171 bp | PMC6121732 (Willard 1985; Waye & Willard 1987) |
| Period search tolerance (±) | 5 bp | Indel/divergence allowance around the 171-bp period (monomers diverge 10–40%) |
| Min periodicity (base-level self-similarity) to call a tandem array | 0.50 | Lower bound of the 50–70% intra-array monomer identity (PMC6121732) |
| Min AT content (AT-richness signature) | > 0.50 | Alpha satellite described as "AT-rich 171-bp alphoid monomer" (PMC6121732); 0.5 = balance point |
| CENP-B box length | 17 bp | Masumoto et al. 1989 |
| CENP-B box consensus | `YTTCGTTGGAARCGGGA` (Y=C/T, R=A/G) | Masumoto et al. 1989 (PMC4843215, PMC6121732) |

**No consensus monomer sequence is embedded in the implementation.** Detection uses tandem-period
self-similarity + AT-richness + CENP-B box (IUPAC) matching, so no alphoid monomer string is invented.
The 62-bp `AlphaSatelliteConsensus` constant predates this work and remains unused by the new methods.

### Residual (out of scope)
Higher-order repeat (HOR) structure — the chromosome-specific organization of monomers into HOR units
(2–34 monomers) and the suprachromosomal family classification — is not modelled. `DetectAlphaSatellite`
detects the monomer-level tandem + AT + CENP-B signal, not the HOR hierarchy.
