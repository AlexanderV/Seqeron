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

### Residual (out of scope, before the 2026-06-24 HOR fix)
Higher-order repeat (HOR) structure — the organization of monomers into HOR units and the
suprachromosomal family classification — was not modelled. `DetectAlphaSatellite` detects the
monomer-level tandem + AT + CENP-B signal, not the HOR hierarchy.

## Higher-Order Repeat (HOR) Structure Detection (added 2026-06-24)

`DetectHigherOrderRepeat` adds **opt-in, additive** detection of the HOR organisation of an
alpha-satellite array. `AnalyzeCentromere`, `DetectAlphaSatellite`, `FindCenpBBoxes`, and the Levan
classification are unchanged.

### Sources actually retrieved this session (2026-06-24)

1. **McNulty SM, Sullivan BA (2018). "Alpha satellite DNA biology: finding function in the recesses
   of the genome." Chromosome Res 26:115-138.** — retrieved via PMC.
   - URL fetched: https://pmc.ncbi.nlm.nih.gov/articles/PMC6121732/
   - Verbatim: *"A defined number of individual monomers (black arrows) that are 50–70% identical in
     sequence are arranged tandemly to form a HOR unit."*
   - Verbatim: *"The individual monomers within a HOR unit have 50–70% identity and can be
     distinguished such that HOR unit length is determined by where the next monomer shows nearly
     total sequence identity to the first monomer in the HOR."*  (← the HOR-period definition)
   - Verbatim: *"… the HORs are repeated hundreds to thousands of times to create homogenous arrays in
     which HOR within a given array are 97–100% identical."*
   - Verbatim: *"HORs within a chromosome-specific array differ in sequence by only a few percent,
     however, HORs between non-homologous chromosomes are only 50–70% identical."*
   - Verbatim: *"Alpha satellite monomers differ in sequence by 10–40% …"*

2. **Rosandić M, Paar V, et al. (2024). "Novel Concept of Alpha Satellite Cascading Higher-Order
   Repeats (HORs) … in … T2T-CHM13 … Chromosome 15." (Int J Mol Sci / MDPI).** — retrieved via PMC.
   - URL fetched: https://pmc.ncbi.nlm.nih.gov/articles/PMC11050224/
   - Verbatim: *"sequences of n monomers, collectively known as nmer HORs."*
   - Verbatim: *"HOR copies are further organized in tandem, with minimal divergence between HOR
     copies, typically less than 5%."*
   - Verbatim: *"The degree of divergence between any two monomers within each HOR copy is significant,
     ranging from approximately 20% to 40%. … Monomers exhibiting less than 5% divergence are
     classified as the same monomer type."*

3. **Alkan C, et al. (2007). "Genome-wide characterization of centromeric satellites …" / ColorHOR
   (Paar et al. 2005, Bioinformatics 21(7):846).** — retrieved via Oxford Academic.
   - URL fetched: https://academic.oup.com/bioinformatics/article/21/7/846/268781
   - Verbatim: *"HORs exhibit mutual sequence divergence of <5%"*; *"Independent, ∼171 bp monomers of
     alpha satellite DNA generally exhibit substantial intermonomeric sequence divergence (20–40%)."*
   - Worked example (chr1 1866-bp 11-mer HOR): segmented into 11 monomers of lengths
     168,171,171,171,171,171,170,167,169,167,170 bp with **1.8%** divergence between HOR copies.

### Derived HOR-detection parameters (all source-backed)

| Parameter | Value | Source |
|-----------|-------|--------|
| Monomer length used to split the array | 171 bp | PMC6121732 (Willard 1985); Alkan 2007 worked example (167–171 bp monomers) |
| Inter-HOR (copy-to-copy) min identity to accept a HOR period | ≥ 95% (i.e. < 5% divergence) | PMC11050224 ("< 5%"); Alkan 2007 ("<5%"); PMC6121732 ("97–100% identical") |
| Intra-HOR monomer identity (distinct monomers within a unit) | 50–70% (10–40% divergence) | PMC6121732; PMC11050224 (20–40% divergence) |
| HOR period definition | smallest k with the next k-shifted monomer near-identical to the first | PMC6121732 ("HOR unit length is determined by where the next monomer shows nearly total sequence identity to the first monomer") |
| HOR-period consistency fraction across the array | 0.90 | Arrays are "homogenous" / "repeated hundreds to thousands of times" (PMC6121732) — periodicity must hold across (essentially) the whole array, not at one pair |

### Method (reuses the library aligner — no external trained data)
Split the array into consecutive 171-bp monomers; compute monomer-vs-monomer percent identity with
`SequenceAligner.GlobalAlign` + `SequenceAligner.CalculateStatistics` (Needleman-Wunsch + EMBOSS-style
identity); the HOR period is the **smallest** block size k (≥1) for which monomers k apart are ≥95%
identical across ≥90% of the array. Period 1 = homogeneous single-monomer (1-mer) array, **not** a
multi-monomer HOR. Report period (monomers/unit), unit length (k×171 bp), copy number
(⌊monomers/k⌋), and mean inter-HOR vs intra-HOR identity. No consensus monomer or family library is
embedded.

### Residual after this fix (honest, data-blocked)
**Suprachromosomal-family / specific α-satellite family (J1/J2/W/…) assignment** is still out of scope:
naming the family a HOR belongs to requires curated reference HOR libraries (chromosome-specific
consensus HORs), which are external trained/curated data the library does not embed. The HOR
*structure* (period, copy number, inter-/intra-HOR identity) is now detected; the *family label* is not.

## Suprachromosomal-Family (SF) Assignment (added 2026-06-25 — bundled CC0 reference)

`AssignSuprachromosomalFamily` adds **opt-in, additive** suprachromosomal-family assignment for a
detected alpha-satellite monomer/HOR, backed by a **bundled CC0 reference**. `AnalyzeCentromere`,
`DetectAlphaSatellite`, `FindCenpBBoxes`, `DetectHigherOrderRepeat`, and the Levan classification
are byte-unchanged.

### Sources actually retrieved this session (2026-06-25)

1. **McNulty SM, Sullivan BA (2018). "Alpha satellite DNA biology …" Chromosome Res; PMC6121732.** —
   page fetched 2026-06-25 (search query "alpha satellite suprachromosomal family SF1 SF2 monomer
   classes J1 J2"; page read).
   - SF taxonomy / periodicity: **SF1** dimeric (J1·J2), **SF2** dimeric (D1·D2), **SF3** pentameric
     (W1–W5), **SF4** monomeric (M1), **SF5** irregular (R1·R2).
   - Verbatim A/B-type assignment: *"A-type monomers include J1, D2, W4, W5, M1, and R2 monomers,
     while B-type consist of J2, D1, W1–W3, and R1 monomers. B-type monomers contain CENP-B boxes;
     A-type contain pJα binding sites."*
   - Chromosome distribution: SF1 → HSA1,3,5,6,7,10,12,16,19; SF2 → HSA2,4,8,9,13,14,15,18,20,21,22;
     SF3 → HSA1,11,17,X; SF4 → HSA13,14,15,21,22,Y; SF5 → HSA5,7,15,19.
   - *"W1–W5 were initially described as A–E monomers."*

2. **Shepelev VA, Uralsky LI, Alexandrov AA, Yurov YB, Rogaev EI, Alexandrov IA (2009). PLOS
   Genetics 5(9):e1000641.** — page fetched 2026-06-25.
   - Foundational SF classification; SF4 = single M1 class with no higher-order periodicity; SF5 =
     R1·R2 alternating irregularly. *"R1 represents the first appearance of novel class B monomers,
     which bind CENP-B protein."*  Licence: **CC BY 4.0**.

3. **Dfam alpha-satellite consensus monomers (bundled, CC0).** — REST API JSON fetched 2026-06-25:
   `https://www.dfam.org/api/families/DF000000029` (ALR), `…/DF000000014` (ALRa), `…/DF000000015`
   (ALRb). Consensus strings copied verbatim into `Resources/AlphaSatelliteReference.fasta`.
   - **ALR** (DF000000029, 171 bp), **ALRa** (DF000000014, 172 bp): no CENP-B box → **A-type**.
   - **ALRb** (DF000000015, 169 bp): contains the 17-bp CENP-B box `YTTCGTTGGAARCGGGA` (substring
     `CTTCGTTGGAAACGGGA`) at consensus position 126 → **B-type** (verified this session by matching
     the IUPAC CENP-B consensus against each Dfam string: ALRb hit at 126; ALR/ALRa no hit).
   - All three classified `root; Tandem_Repeat; Satellite; Centromeric`. Licence: **CC0** (Storer et
     al. 2021, *Mobile DNA* 12:2; dfam.org).

4. **T2T Consortium / marbl CHM13** — page fetched 2026-06-25; verbatim: *"All data is released to
   the public domain (CC0) and we encourage its reuse."* (corroborates the SF facts; CC0).

### Bundled reference (CC0) and derived parameters

| Item | Value | Source |
|------|-------|--------|
| Bundled reference monomers | Dfam ALR (A), ALRa (A), ALRb (B) | Dfam DF000000029/14/15 (CC0) |
| Min identity to call alpha-satellite | ≥ 60% to closest reference | **ASSUMPTION** — conservative gate vs the ~16% (consensus) / 50–70% (inter-monomer) divergence (Waye & Willard 1987; PMC6121732) |
| A/B-box typing | best-match reference box type; ALRb=B (CENP-B box), ALR/ALRa=A | PMC6121732 A/B rule + Masumoto 1989 box; Dfam sequences |
| SF3 ⇔ HOR period a positive multiple of 5 | period % 5 == 0 (≥5) | **ASSUMPTION** — pentameric-ancestry proxy (Waye & Willard 1986; Shepelev 2009) |
| SF4 ⇔ monomeric (period 1), A-type only | M1 is A-type, monomeric | PMC6121732 |
| {SF1,SF2} ⇔ dimeric (period 2), A+B per unit | both are A→B dimers | PMC6121732 |
| SF5 ⇔ irregular A/B mix, no regular period | R1·R2 irregular | PMC6121732; Shepelev 2009 |

### Method (reuses the library aligner + the validated HOR detector — no external trained model)
Split the array into 171-bp monomers; best-match each monomer to the bundled reference with
`SequenceAligner.GlobalAlign` + `CalculateStatistics().Identity` to (a) confirm alpha-satellite
identity (≥ 60%) and (b) type it A vs B; take the HOR period from `DetectHigherOrderRepeat`; assign
the SF from the period + the A/B composition of one HOR unit per the rules above.

### Residual after this fix (sharpened — honest, data-blocked)
The bundled CC0 reference resolves **SF3** (pentameric), **SF4** (monomeric A-type) and **SF5**
(irregular A/B), and narrows dimeric arrays to **{SF1, SF2}**. It does **NOT** separate **SF1 from
SF2** (both dimeric with the identical A→B box pattern), nor does it tag SF3 arrays whose detected
HOR period is not a multiple of 5 (e.g. the diverged-pentamer dodecameric DXZ1, period 12). Doing
either needs the **SF-resolved consensus monomer library** (J1/J2/D1/D2/W1–W5/M1/R1/R2). Those
SF-resolved consensus sequences are published only in supplementary matrices that are not
machine-retrievable as FASTA and in third-party HMM repositories (enigene/HumAS-HMMER,
logsdon-lab) that ship **no LICENSE file** — therefore not redistributable here (cf. the
TIGRFAM/LM22 non-redistribution rule). Callers who hold an SF-resolved consensus set can pass it to
`AssignSuprachromosomalFamily(sequence, reference)`.
