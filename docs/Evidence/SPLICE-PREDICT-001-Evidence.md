# Evidence Artifact: SPLICE-PREDICT-001

**Test Unit ID:** SPLICE-PREDICT-001
**Algorithm:** Gene Structure Prediction (Intron/Exon)
**Date Collected:** 2026-02-12

---

## Online Sources

### 1. Wikipedia — Intron

**URL:** https://en.wikipedia.org/wiki/Intron
**Accessed:** 2026-02-12
**Authority rank:** 4 (cites primary sources)

**Key Extracted Points:**

1. **GT-AG rule:** Spliceosomal introns are characterized by specific sequences at the
   boundaries between introns and exons. The 5' splice site (donor) contains GT (GU in RNA)
   and the 3' splice site (acceptor) contains AG — citing Padgett et al. (1986).
2. **Intron length range:** Shortest known metazoan intron is 30 bp (human MST1L gene,
   Piovesan et al. 2015). Longest exceeds 3.6 Mb (Drosophila DhDhc7, Reugels et al. 2000).
3. **Intron types:** Four main types: spliceosomal, tRNA, group I, group II. Major (U2)
   spliceosomal uses GT-AG; minor (U12) uses AT-AC.
4. **Branch point:** Located near 3' end of intron, becomes covalently linked to 5' end
   during splicing, generating a lariat structure.
5. **Splicing accuracy:** Under ideal conditions ~99.999% accurate per intron
   (Hsu & Hertel, 2009).

### 2. Wikipedia — Exon

**URL:** https://en.wikipedia.org/wiki/Exon
**Accessed:** 2026-02-12
**Authority rank:** 4

**Key Extracted Points:**

1. **Definition:** An exon is any part of a gene that will form a part of the final mature
   RNA after introns have been removed — Gilbert (1978).
2. **Average exon count:** Across eukaryotic genes, average 5.48 exons per protein coding
   gene — Sakharkar et al. (2002).
3. **Average exon size:** Encodes 30–36 amino acids (~90–108 nt) — Sakharkar et al. (2002).
4. **Exon types:** Initial (first, contains 5' UTR/start codon), Internal (middle),
   Terminal (last, contains stop codon/3' UTR), Single (single-exon genes have no introns).

### 3. Wikipedia — Gene Structure

**URL:** https://en.wikipedia.org/wiki/Gene_structure
**Accessed:** 2026-02-12
**Authority rank:** 4

**Key Extracted Points:**

1. **Eukaryotic gene structure:** Transcripts subdivided into exon and intron regions.
   Introns are spliced out during post-transcriptional processing — Matera & Wang (2014).
2. **Spliced sequence:** Once spliced together, exons form a single continuous
   protein-coding region; splice boundaries are not detectable in mature mRNA.
3. **Single-exon genes:** Some genes lack introns entirely (e.g., most prokaryotic genes,
   some eukaryotic genes like histone genes).

### 4. Breathnach & Chambon (1981)

**Citation:** Breathnach R, Chambon P. Organization and expression of eucaryotic split
genes coding for proteins. Annu Rev Biochem. 1981;50:349-383.
**Authority rank:** 1

**Key Extracted Points:**

1. **GT-AG consensus:** >99% of spliceosomal introns begin with GT and end with AG.
2. **Donor consensus:** MAG|GURAGU (M=A/C, R=A/G), where | marks the exon-intron boundary.
3. **Acceptor consensus:** (Y)nNCAG|G, where Y=pyrimidine, N=any nucleotide.

### 5. Shapiro & Senapathy (1987)

**Citation:** Shapiro MB, Senapathy P. RNA splice junctions of different classes of
eukaryotes: sequence statistics and functional implications in gene expression.
Nucleic Acids Res. 1987;15(17):7155-7174.
**Authority rank:** 1

**Key Extracted Points:**

1. **Splice site scoring:** Established position weight matrix approach for scoring
   splice sites based on nucleotide frequencies at each position.
2. **Statistical basis:** Analyzed >2,000 splice junctions from multiple species.

### 6. Burge, Tuschl & Sharp (1999)

**Citation:** Burge CB, Tuschl T, Sharp PA. Splicing of precursors to mRNAs by the
spliceosomes. In: Gesteland RF, Cech TR, Atkins JF, eds. The RNA World. 2nd ed.
Cold Spring Harbor Laboratory Press; 1999:525-560.
**Authority rank:** 1

**Key Extracted Points:**

1. **U2-type introns:** Use GT-AG termini, constitute ~99% of spliceosomal introns.
2. **U12-type introns:** Use AT-AC termini, constitute ~0.5% of introns.
3. **GC-AG variant:** ~0.5% of U2-type introns use GC instead of GT at 5' end.
4. **Intron pairing:** Each donor must pair with exactly one acceptor to define an intron.

---

## Documented Corner Cases and Failure Modes

### From Wikipedia (Intron)

1. **Extremely short introns:** Shortest metazoan intron is 30 bp. Implementation uses
   default minIntronLength=60, which is biologically reasonable for most organisms.
2. **Overlapping introns:** In a computational prediction, multiple donor-acceptor pairs
   may overlap; a greedy non-overlapping selection is needed.
3. **No splice sites found:** Sequences without GT-AG or AT-AC pairs produce no introns
   and should default to single-exon structure.

### From Breathnach & Chambon (1981)

1. **Non-canonical splice sites:** Rare non-GT-AG introns exist (~1%).
   Implementation handles GC-AG and AT-AC variants.

---

## Test Datasets

### Dataset 1: Constructed Two-Exon Gene

**Source:** Synthetic, based on consensus sequences from Breathnach & Chambon (1981)

| Parameter | Value |
|-----------|-------|
| Exon 1 | 35 nt: AUGCCCAAAGGGCCCUUUAAAGGGCCCUUUAAAGC |
| Donor | GUAAGU (GT consensus) |
| Intron body | 60 × A (filler) |
| PPT | UUUUUUUUUUUUUU (14 nt polypyrimidine tract) |
| Acceptor | CAG |
| Exon 2 | 35 nt: GCCUUUAAAGGGCCCUUUAAAGGGCCCUUUAAAGC |
| Total length | 35 + 6 + 60 + 14 + 3 + 35 = 153 nt |
| Expected introns | 1 (from donor to acceptor position) |
| Expected exons | 2 (Initial + Terminal) |
| Spliced seq length | 35 + 35 = 70 nt |

### Dataset 2: Single-Exon Gene (No Introns)

**Source:** Synthetic, based on histone gene structure (single-exon)

| Parameter | Value |
|-----------|-------|
| Sequence | 50 nt: AUGAAAGGGCCCUUUAAAGGGCCCUUUAAAGGGCCCUUUAAAGGGCCCUAA |
| Expected introns | 0 |
| Expected exons | 1 (Single type) |
| Exon span | 0 to length-1 |

### Dataset 3: Empty/Null Input

**Source:** Trivial edge case

| Parameter | Value |
|-----------|-------|
| Input | "" or null |
| Expected introns | 0 |
| Expected exons | 0 |
| Spliced sequence | empty |
| Overall score | 0 |

---

## Assumptions

1. **ASSUMPTION: Score threshold for greedy selection** — The greedy non-overlapping intron
   selection by descending score is an implementation-specific heuristic, not formally
   specified in literature. Most ab initio gene finders (e.g., GenScan, Augustus) use
   dynamic programming instead. The greedy approach is acceptable for a simplified predictor.

2. **ASSUMPTION: Exon phase calculation** — Phase = (sum of previous exon lengths) mod 3.
   This follows the standard convention for reading frame tracking across introns
   (Alberts et al., 2002, Molecular Biology of the Cell).

3. **ASSUMPTION: Overall score as mean of intron scores** — Using the arithmetic mean of
   individual intron scores as the overall gene structure quality metric is
   implementation-specific; no standard is defined in the literature.

4. **ASSUMPTION: Minimum intron length default (60)** — The default of 60 bp is reasonable
   for most eukaryotes but excludes the shortest known metazoan introns (30 bp).
   This is a documented scope limitation.

---

## Recommendations for Test Coverage

1. **MUST Test:** Empty/null input → empty GeneStructure — Evidence: trivial edge case
2. **MUST Test:** Single-exon gene (high threshold, no introns) → single exon type — Evidence: Gilbert (1978)
3. **MUST Test:** Two-exon gene with clear GT-AG intron → 1 intron, 2 exons — Evidence: Breathnach & Chambon (1981)
4. **MUST Test:** Intron has GT donor and AG acceptor — Evidence: GT-AG rule
5. **MUST Test:** PredictIntrons respects minIntronLength — Evidence: algorithm parameter
6. **MUST Test:** PredictIntrons respects maxIntronLength — Evidence: algorithm parameter
7. **MUST Test:** Spliced sequence excludes intron — Evidence: splicing definition
8. **MUST Test:** Exon types assigned correctly (Initial/Internal/Terminal/Single) — Evidence: Gilbert (1978)
9. **MUST Test:** Exon phase tracks reading frame — Evidence: Alberts et al. (2002)
10. **MUST Test:** Score in valid range [0, 1] — Evidence: implementation normalization
11. **SHOULD Test:** Non-overlapping intron selection — Evidence: algorithm invariant
12. **SHOULD Test:** DNA T-equivalence (T treated as U) — Evidence: implementation converts T→U
13. **SHOULD Test:** Intron type classification (U2 for GT-AG) — Evidence: Burge et al. (1999)
14. **SHOULD Test:** Short sequence (< minimum viable) handling — Evidence: guard clause
15. **COULD Test:** Multi-intron gene structure — Evidence: Sakharkar et al. (2002)
16. **COULD Test:** Case insensitivity — Evidence: implementation uses ToUpperInvariant

---

## References

1. Gilbert W (1978). Why genes in pieces? Nature 271:501.
2. Breathnach R, Chambon P (1981). Organization and expression of eucaryotic split genes
   coding for proteins. Annu Rev Biochem 50:349-383.
3. Shapiro MB, Senapathy P (1987). RNA splice junctions of different classes of eukaryotes.
   Nucleic Acids Res 15:7155-7174.
4. Burge CB, Tuschl T, Sharp PA (1999). Splicing of precursors to mRNAs by the spliceosomes.
   In: The RNA World, 2nd ed. CSHL Press, pp. 525-560.
5. Sakharkar MK et al. (2002). ExInt: an Exon Intron Database. Nucleic Acids Res 30:191-194.
6. Piovesan A et al. (2015). Identification of minimal eukaryotic introns through GeneBase.
   DNA Res 22:495-503.
7. Hsu SN, Hertel KJ (2009). Spliceosomes walk the line. RNA Biol 6:526-530.
8. Alberts B et al. (2002). Molecular Biology of the Cell, 4th ed. Garland Science.
9. Padgett RA et al. (1986). Splicing of messenger RNA precursors. Annu Rev Biochem 55:1119-1150.

---

## Change History

- **2026-02-12**: Initial documentation.
