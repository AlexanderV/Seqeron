# Evidence Artifact: SPLICE-PREDICT-001

**Test Unit ID:** SPLICE-PREDICT-001
**Algorithm:** Gene Structure Prediction (Intron/Exon)
**Date Collected:** 2026-02-12

---

## Online Sources

### 1. Wikipedia — Intron

**URL:** https://en.wikipedia.org/wiki/Intron
**Accessed:** 2026-02-12 · re-fetched 2026-10-14 (SRC citation audit: supports the claims above)
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
**Accessed:** 2026-02-12 · re-fetched 2026-10-14 (SRC citation audit: supports the claims above)
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
**Accessed:** 2026-02-12 · re-fetched 2026-10-14 (SRC citation audit: see corrections below)
**Authority rank:** 4

**Key Extracted Points:**

1. **Eukaryotic gene structure:** Transcripts subdivided into exon and intron regions.
   Introns are spliced out during post-transcriptional processing — Matera & Wang (2014).
2. **Spliced sequence:** Once spliced together, exons form a single continuous
   protein-coding region; splice boundaries are not detectable in mature mRNA.
3. **Single-exon genes:** Some genes lack introns entirely (e.g., most prokaryotic genes,
   some eukaryotic genes; the Gene structure page itself only states "Introns are extremely rare in prokaryotes").

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

No residual assumptions. All former assumptions have been resolved:

1. **~RESOLVED: Score threshold for greedy selection~** — Reclassified as Design Decision.
   External sources prescribe the biological invariant (introns don't overlap) but not
   the computational algorithm. Greedy-by-score is a valid simplified approach.

2. **~RESOLVED: Exon phase calculation~** — Mathematically trivial. Phase = (sum of
   previous exon lengths) mod 3; first exon has 0 preceding lengths = phase 0 by definition.
   Standard convention per Alberts et al. (2002).

3. **~RESOLVED: Overall score as mean of intron scores~** — Reclassified as Definition.
   No biological standard exists for a gene structure quality metric. This is a defined
   formula, not a biological assumption.

4. **~RESOLVED: Default branch point score of 0.3~** — Fixed in code. When no branch
   point is found, combined score = (donor.Score + acceptor.Score) / 2 instead of
   (donor.Score + acceptor.Score + 0.3) / 3. Eliminated arbitrary magic constant.

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


## B11 Review (2026-11) — sources opened and oracle numbers

Sources actually consulted: Ensembl REST (rest.ensembl.org, GRCh38; canonical transcripts ENST00000335295 HBB,
ENST00000269305 TP53, ENST00000003084 CFTR exons 2-5, neutral User-Agent), the SUPPA/rMATS event definitions
(Wang 2008; Shen 2014; Alamancos 2015 - from the published definitions, not re-fetched), Braunschweig 2014 (retained
intron features), Gao 2008 / Levine & Durbin 2001 / Dietrich 1997 (already in the donor/acceptor evidence).

- Independent Python oracle (pairing, types, Gao branch point, scores, ordering, greedy selection with exon tiling,
  phase, SE/A5SS/A3SS/MXE events, IsWithinCodingRegion for every position x frame): 0 mismatches on 1,400 random
  sequences (120-1,200 nt, GC 0.4-0.6, planted introns; ~540,000 candidate introns, ~3,300 selected, ~150,000 events,
  ~110,000 coding positions). Donor/acceptor site lists are the canonical (separately validated) methods' output.
- Ensembl recall (minIntron 30, minExon 30): HBB 2/2, TP53 10/10, CFTR exons 2-5 3/3 annotated introns are among the
  PredictIntrons candidates (exact coordinates, all typed U2; true-intron score rank 2-4,500 of 3,011 / 602,245 /
  1,293,492 candidates); the greedy PredictGeneStructure selection recovers 0/15 (it picks high-scoring false pairs):
  splice-motif only, see the wiki gotcha.

### A3 (B11 F45, 2026-11): splice + open-reading-frame dynamic program `PredictGeneStructureDp`

Sources actually opened: Ensembl REST (rest.ensembl.org, GRCh38; lookup/id + sequence/region, neutral User-Agent) for ENST00000335295 (HBB, 11:5225464-5227071:-1), ENST00000269305 (TP53, 17:7668421-7687490:-1, CDS from the first CDS base; 5' UTR intron 1 excluded as the model keeps UTRs intron-free), ENST00000003084 (CFTR exons 2-5, window 7:117504203..exon-5 end + 3 kb cut 100 nt after the first in-frame stop downstream of exon 5; the first complete codon of exon 2 is the start), plus 12 training genes (ACTB, GAPDH, HRAS, INS, GH1, TPI1, LDHA, SOD1, APOE, TTR, CYCS, GCK, ALB; MYC has no ATG at the Ensembl CDS start, HPRT1/PGK1 > 60 kb skipped). A web search snippet list for GENSCAN (Burge & Karlin 1997: separate exon / intron length models, phase states; training 238 genes) was the only access to that paper; GeneSplicer (Pertea 2001) and the papers themselves were NOT opened - the model description is from the audit brief and general knowledge, and no GENSCAN / GeneSplicer parameter is used. No published minimum exon length, intron prior or per-intron cost was found.

- **Model:** exact Bellman recursion (see the algorithm doc 5.5). Independent oracle: exhaustive enumeration of all non-overlapping intron subsets of random candidate sets (spliced-sequence based ORF test, independent code), 1,200 random sequences of 40-120 nt (random donors / acceptors with random bit scores, random min exon 1-8 / min intron 1-19 / max intron, a third with a pseudo-random exon coding scorer): 0 mismatches in optimum score (<= 1e-9) and no ORF / ORF status; 939 cases with a unique optimum also match structure-for-structure; 394 optima contain introns, 261 cases have no open frame (both return none).
- **Ensembl recall (exact intron coordinates, 14 coding introns: HBB 2, TP53 9, CFTR exons 2-5 3; the TP53 5' UTR intron is outside the model, 15 in the audit counting):** greedy `PredictGeneStructure` 0/14 (F32-F38). DP, open-frame constraint only (defaults: min exon 3, min intron 60, intron penalty 12): HBB 2/2 (2 predicted), TP53 6/9 (13 predicted), CFTR 0/3 (20 predicted; donors 3/3 and acceptors 3/3 each hit, paired wrongly): total 8/14 recall 57%, precision 8/35 = 23%. Boundary accuracy (donor / acceptor individually): HBB 2/2 + 2/2, TP53 9/9 + 8/9, CFTR 3/3 + 3/3.
- With a user exon coding score (sum of in-frame CPAT hexamer log-ratios in bits, `GenomeAnnotator.CalculateCodingPotential` x hexamers, table from `BuildHexamerTable` on the 12 training genes' CDS vs introns, embedded as `hexamer_train12.tsv`): HBB 2/2 (2 predicted), TP53 8/9 (10), CFTR 0/3 (17; donors 1/3, acceptors 0/3): 10/14 recall 71%, precision 10/29 = 34% (weight 0.5: 11/14). CFTR (22-kb introns) stays unresolved. Time: TP53 19 kb 0.09 s (candidate generation included), CFTR window 30 kb ~0.01 s; memory O(acceptors x 31) doubles.
- **Intron penalty / defaults (no published values):** the 12-gene training set (65 coding introns) with min exon 3 / site floor 0: penalty 0: 15/65 recovered of 772 predicted; 4: 16/676; 8: 19/435; 10: 26/305; 12: 32/208; 14: 34/119; 16: 25/68. Site floors 4 and min exon 20 gave 38/152 at 12 but lose true weak sites, so the defaults stay 0 / 3 (documented tunable; the penalty is the only tuned value and it was chosen without the evaluation genes). HBB and TP53 are recovered at penalties 8-16; CFTR needs the coding score.

---

## Change History

- **2026-02-12**: Initial documentation.
- **2026-11**: B11 review (F32-F38): typed U12 requires both U12 sites, canonical branch point, fixed-term score, exon tiling, SUPPA/rMATS events, ORF-based IsWithinCodingRegion.
- **2026-11 (audit A10, B11 F43)**: `FindRetainedIntronCandidates` rewritten (see B11.md F43). Opened: Braunschweig et al. 2014 Genome Res 24:1774, PMC4216919 full text (retained introns "on average shorter, more C/G rich, and associated with weaker splice sites", all P < 0.001; IR code AUC 0.79; PIR = 100 x avg(E1I, IE2)/(E1E2 + avg(E1I, IE2)), intron-body window; no numeric length/GC/MaxEnt cut-offs printed, hence no defaults), Jacob & Smith 2017 Nat Rev Genet 18:437 (not opened; used for the ORF-open / NMD statement via the search summary only). Oracle: independent Python frame/GC/score/length filter over the C# `PredictIntrons` output vs C# on 300 random sequences (50,343 predicted introns, 269 with candidates, random criteria combinations): 0 mismatches. Ensembl REST canonical transcripts (HBB 2/2, TP53 9/10, CFTR 15/26 introns <= 4 kb, +-99 nt flanks): C# == oracle in all windows; 0 known (constitutive) introns flagged (their frames contain stops or length not a multiple of 3), as expected.
- **2026-11 (audit A3, B11 F45)**: `PredictGeneStructureDp` (splice + open-reading-frame dynamic program), MCP `model=splice_dp`; oracle numbers above.
