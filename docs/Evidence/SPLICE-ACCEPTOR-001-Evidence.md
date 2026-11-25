# Evidence Artifact: SPLICE-ACCEPTOR-001

**Test Unit ID:** SPLICE-ACCEPTOR-001
**Algorithm:** Acceptor Site Detection
**Date Collected:** 2026-02-12

---

## Online Sources

### 1. Wikipedia — RNA splicing

**URL:** https://en.wikipedia.org/wiki/RNA_splicing
**Accessed:** 2026-02-12 · re-fetched 2026-10-14 (SRC citation audit: supports the claims above)
**Authority rank:** 4 (Wikipedia citing primary sources)

**Key Extracted Points:**

1. **3' splice site consensus:** Y-rich-N-C-A-G|G — the acceptor site terminates the intron with an almost invariant AG sequence (citing Molecular Biology of the Cell).
2. **Polypyrimidine tract (PPT):** Upstream of the AG, a region high in pyrimidines (C and U) is located. Further upstream is the branch point.
3. **GU-AG rule:** The major spliceosome splices introns containing GU at the 5' splice site and AG at the 3' splice site. This accounts for >99% of intron splicing (canonical splicing / lariat pathway).
4. **Non-canonical splicing:** When intronic flanking sequences do not follow the GU-AG rule, noncanonical splicing occurs via the minor spliceosome.

---

### 2. Wikipedia — Polypyrimidine tract

**URL:** https://en.wikipedia.org/wiki/Polypyrimidine_tract
**Accessed:** 2026-02-12 · re-fetched 2026-10-14 (SRC citation audit: see corrections below)
**Authority rank:** 4 (Wikipedia citing Lodish et al. 2004, Wagner & Garcia-Blanco 2001)

**Key Extracted Points:**

1. **PPT definition:** Region of pre-mRNA rich in pyrimidine nucleotides (especially uracil), usually 15–20 base pairs long.
2. **PPT location:** Located about 5–40 base pairs before the 3' end of the intron.
3. **U2AF binding:** the page states only that protein factors bind the PPT, "including the spliceosome component U2AF" (and PTB). The U2AF65 (PPT) / U2AF35 (3' AG) subunit split is NOT on this page (re-checked 2026-10-14); see Burge et al. (1999), source 5, for U2AF65 recruitment.
4. **PPT quality:** Strong PPT (continuous pyrimidines) produces robust splicing; weak PPT (interrupted by purines) can be suppressed by PTB binding.

---

### 3. Wikipedia — Spliceosome

**URL:** https://en.wikipedia.org/wiki/Spliceosome
**Accessed:** 2026-02-12 · re-fetched 2026-10-14 (SRC citation audit: see corrections below)
**Authority rank:** 4 (Wikipedia citing Will & Lührmann 2011, Burge et al. 1999)

**Key Extracted Points:**

1. **3' splice site elements:** 3' splice site defined by AG dinucleotide, preceded by polypyrimidine tract (PPT) and branch point sequence (BPS).
2. **Branch point distance:** not stated on this page; the distance used is from Mercer et al. (2015) ("most branch sites are between 19 and 35 nucleotides away from the 3′ss") and Gao et al. (2008) (−34..−21), see below.
3. **Minor spliceosome:** the page mentions the minor (U12-type) spliceosome (citing Patel & Steitz 2003) but does not state the AT-AC dinucleotides; the U12-type AC acceptor is per Burge et al. (1999), source 5 (corrected 2026-10-14, SRC citation audit).

---

### 4. Shapiro & Senapathy (1987)

**Citation:** Shapiro MB, Senapathy P. "RNA splice junctions of different classes of eukaryotes: sequence statistics and functional implications in gene expression." Nucleic Acids Research 15(17):7155–7174, 1987.
**Authority rank:** 1 (peer-reviewed paper)

**Opened (B11 review, 2026-11-25):** the scanned paper as PDF via WebFetch of
`europepmc.org/api/getPdf?pmcid=PMC306199`, pages rendered with `pdftoppm` (150/400 dpi) and read visually:
Methods p. 7157–7158, Table 1 p. 7159 (3' ss subtables), Table 5 p. 7168 (published 3' ss scores).

**Key Extracted Points:**

1. **Table 1 (3' ss, positions −14..+1, A C G T %, 1432 total / 542 primate / 378 rodent / 67 other mammal /
   202 other vertebrate / 86 invertebrate / 110 plant / 47 viral sites).** Primate: −14 9 31 15 45; −13 9 33 13 45;
   −12 7 31 11 51; −11 7 35 7 51; −10 10 35 7 47; −9 10 35 11 44; −8 7 43 7 42; −7 9 41 8 42; −6 6 39 6 48;
   −5 6 40 8 46; −4 23 29 23 24; −3 3 74 1 22; −2 100 0 0 0; −1 0 0 100 0; +1 28 13 49 10. All eight subtables are
   transcribed in `SpliceSitePredictor.ShapiroSenapathyAcceptorTables`.
2. **Methods (3' ss score, verbatim):** "to the left of the NCAGG consensus, the best 8 of 10 locations were used,
   since the total pyrimidine content of these 10 locations is about 80% and most 3' splice sites have at least two
   purines in these 10 locations. Thus a pyrimidine sequence containing two purines can score as well as one
   containing no purines. Based on the fact that less than 2% of the 3' splice sites contained an AG within 10 bases
   upstream of the true 3' ss (Figure 4), an added restriction was made that any sequence containing an AG in those
   10 locations be eliminated from consideration." `score = 100((t1−l1)/(h1−l1)+(t2−l2)/(h2−l2))/2`; t1 = sum of best
   8 of 10 percentages at pyrimidine positions; t2 = sum of 4 percentages at CAGG positions; l1 = sum of lowest 8 of
   lowest 10; l2 = sum of lowest 4 at CAGG; h1 = sum of highest 8 of highest 10; h2 = sum of highest 4 at CAGG.
3. **Position −2 (A), −1 (G):** 100 %; −3 C 72–79 % (plant 62 %); +1 G 49–66 %.

---

### 5. Burge, Tuschl & Sharp (1999)

**Citation:** Burge CB, Tuschl T, Sharp PA. "Splicing precursors to mRNAs by the spliceosomes." In: Gesteland RF, Cech TR, Atkins JF (eds.), The RNA World, 2nd ed. Cold Spring Harbor Laboratory Press, pp. 525–560, 1999.
**Authority rank:** 1 (peer-reviewed book chapter)

**Key Extracted Points:**

1. **3' splice site structure:** (Yn)NYAG|G — polypyrimidine tract followed by NYAG, where N is any nucleotide, Y is pyrimidine.
2. **PPT functional role:** Critical for U2AF65 recruitment and branch point recognition.
3. **GC-AG introns:** ~1% of mammalian introns use GC at the 5' splice site but retain AG at the 3' acceptor site.
4. **U12-type acceptor:** Uses AC instead of AG at the 3' splice site (~0.4% of human introns).

---

### 6. Yeo & Burge (2004)

**Citation:** Yeo G, Burge CB. "Maximum entropy modeling of short sequence motifs with applications to RNA splicing signals." Journal of Computational Biology 11(2–3):377–394, 2004. DOI 10.1089/1066527041410418.
**Authority rank:** 1 (peer-reviewed paper) for the model; rank 3 (reference implementation) for the factorisation + tables.
**Retrieved this session (2026-06-24):** the reference `score3` factorisation and the precomputed probability tables were fetched VERBATIM from the MIT-licensed maxentpy port (kepbod/maxentpy):
- factorisation: `https://raw.githubusercontent.com/kepbod/maxentpy/master/maxentpy/maxent.py` (WebFetch)
- tables: `https://raw.githubusercontent.com/kepbod/maxentpy/master/maxentpy/data/score3_matrix.txt` (curl, 82 560 records)
- licence: `https://raw.githubusercontent.com/kepbod/maxentpy/master/LICENSE` (MIT, fetched verbatim)

**Key Extracted Points:**

1. **MaxEntScan score3ss acceptor model:** Scores a 23-nt window (20 intronic + 3 exonic positions; conserved AG at 0-based positions 18-19) for the 3' splice site. The maximum-entropy model captures position dependencies, outperforming simple PWM models. The score is `log2( P_maxent(seq) / P_background(seq) )`.

2. **Factorisation (verbatim from maxentpy `score3`):** the conserved AG dinucleotide is scored separately and removed; the remaining 21-nt "rest" sequence is factorised over 9 overlapping sub-sequences. With `rest = window[0:18] + window[20:23]`:
   - AG term: `cons1_3[A] * cons2_3[G] / (bgd_3[A] * bgd_3[G])`
   - numerator (multiply): matrix0·rest[0:7], matrix1·rest[7:14], matrix2·rest[14:21], matrix3·rest[4:11], matrix4·rest[11:18]
   - denominator (divide): matrix5·rest[4:7], matrix6·rest[7:11], matrix7·rest[11:14], matrix8·rest[14:18]
   - final score = `log2(AG_term * rest_score)`; `hashseq` = base-4 encoding A=0,C=1,G=2,T=3 (T==U).
   - Probabilities: `bgd_3 = {A:0.27,C:0.23,G:0.23,T:0.27}`; `cons1_3 = {A:0.9903,C:0.0032,G:0.0034,T:0.0030}`; `cons2_3 = {A:0.0027,C:0.0037,G:0.9905,T:0.0030}`.

3. **Tables (`score3_matrix.txt`, 82 560 records):** flat `index hash probability` rows. The nine matrices' sizes match the factorisation: indices 0–4 are 7-nt (4^7 = 16384 entries each), index 5 = 3-nt (64), index 6 = 4-nt (256), index 7 = 3-nt (64), index 8 = 4-nt (256). Total 82 560 — confirms the embedded file is complete and matches the published factorisation.

4. **Worked examples (from the maxentpy `score3` docstring, the documented reference values):**
   - `score3('ttccaaacgaacttttgtAGgga')` → **2.89** (reproduced 2.886773 this session)
   - `score3('tgtctttttctgtgtggcAGtgg')` → **8.19** (reproduced 8.190965)
   - `score3('ttctctcttcagacttatAGcaa')` → **-0.08** (reproduced -0.080278)
   The canonical 2.89 value is the example given on the MaxEntScan website / maxentpy README.

5. **PROVENANCE + LICENCE (flagged, not buried):** the embedded table (`Data/maxent_score3.txt`) and the factorisation come from **kepbod/maxentpy, which is MIT-licensed** (redistribution permitted — full MIT text recorded in `Data/maxent_score3.LICENSE.md`). maxentpy's README notes the *original* Burge-lab Perl scripts/models carry **academic terms** (`http://genes.mit.edu/burgelab/maxent/download/READTHIS`); the artifact bundled here is the MIT-licensed port table, not the original. A maintainer seeking belt-and-suspenders clearance for commercial redistribution should review the upstream Burge-lab terms directly. This licence flag is repeated prominently in `Data/maxent_score3.LICENSE.md` and the validation report.

---

### 7. Gao, Masuda, Matsuura & Ohno (2008) — human branch-point consensus (PRIMARY)

**Citation:** Gao K, Masuda A, Matsuura T, Ohno K. "Human branch point consensus sequence is yUnAy." Nucleic Acids Research 36(7):2257–2267, 2008. DOI 10.1093/nar/gkn073.
**Authority rank:** 1 (peer-reviewed paper)
**Retrieved this session (2026-06-24):** WebSearch + WebFetch of PMC full text https://pmc.ncbi.nlm.nih.gov/articles/PMC2367711/ and the publisher page https://academic.oup.com/nar/article/36/7/2257/2410214

**Key Extracted Points (verbatim / closely paraphrased from the fetched text):**

1. **Consensus:** "the human consensus BPS is yUnAy" — "The consensus BPS thus becomes … yUnAy at positions −3 to +1", where the underlined A is "the branch point at position zero". The motif is 5 nt spanning positions −3,−2,−1,0,+1; the branch-point adenosine is the 4th character.
2. **Per-position conservation (lariat RT-PCR, n = 181 confirmed branch sites):**
   - Position −3 (y): C = 47.0%, U = 32.0% → pyrimidine = 79.0%
   - Position −2 (U): U = 74.6% ("We observed U at position −2 in 74.6% branch sites")
   - Position 0 (A, branch point): A = 92.3% ("A at position 0 in 92.3% branch sites")
   - Position +1 (y): C = 33.1%, U = 42.0% → pyrimidine = 75.1%
   - Position −1 = n (any; uninformative)
3. **Location upstream of 3' splice site:** "Among the 181 sites, 150 (83%) were at positions −34 to −21"; **median −26, mean −27.7 ± 7.6** nt upstream of the 3' end of the intron.
4. **Polypyrimidine tract:** "the polypyrimidine tract spans 4–24 nucleotides downstream of the branch point" (i.e. between the branch point and the AG).

---

### 8. Mercer et al. (2015) — genome-wide branch-point distribution (corroborating)

**Citation:** Mercer TR, Clark MB, Andersen SB, et al. "Genome-wide discovery of human splicing branchpoints." Genome Research 25(2):290–303, 2015. DOI 10.1101/gr.182899.114.
**Authority rank:** 1 (peer-reviewed paper)
**Retrieved this session (2026-06-24):** WebSearch returning the Genome Research / PMC4315302 record.

**Key Extracted Points:**

1. **59,359 high-confidence human branchpoints** in >10,000 genes.
2. **Branchpoints are predominantly adenosine, highly conserved, and closely distributed to the 3′ splice site**: "most branch sites are between 19 and 35 nucleotides away from the 3′ss" — corroborates the Gao (2008) −34..−21 core and supports the 18–40 nt search envelope used by the implementation.

---

## Documented Corner Cases and Failure Modes

### From Shapiro & Senapathy (1987) / Burge et al. (1999)

1. **Non-AG acceptors (U12):** Minor spliceosome introns use AC dinucleotide at 3' splice site instead of canonical AG.
2. **Weak PPT:** Acceptor sites with interrupted or short polypyrimidine tracts score lower and may be skipped during alternative splicing.
3. **Cryptic acceptor sites:** AG dinucleotides in intronic sequence can create decoy acceptor sites; context (PPT quality, AG position) distinguishes real from cryptic.
4. **Short sequences:** Acceptor site detection requires sufficient upstream context for PPT assessment — minimum ~20 nt.

### From implementation

5. **Sequence length < 20:** Returns empty — guard condition prevents scanning without sufficient PPT context.
6. **No AG dinucleotide in scannable range:** Returns empty.

---

## Test Datasets

### Dataset 1: Canonical CAG acceptor with strong PPT

**Source:** Derived from consensus (Shapiro & Senapathy 1987)

| Parameter | Value |
|-----------|-------|
| Sequence | `UUUUUUUUUUUUUUUUCAGGG` (21 nt) |
| PPT region | positions 0–15: continuous U (strong PPT) |
| AG position | positions 16–17 |
| Expected | 1 Acceptor at position 18, S&S primate score 1.0 (best 8 of 10 U + CAG|G all row maxima) |

### Dataset 2: No AG → empty result

**Source:** Negative control (trivially correct)

| Parameter | Value |
|-----------|-------|
| Sequence | `UUUUUUUUUUUUUUUUUUUUU` (21 nt) |
| Expected | Empty result (no AG dinucleotide) |

### Dataset 3: Strong vs weak PPT comparison

**Source:** PPT quality principle (Burge et al. 1999)

| Parameter | Value |
|-----------|-------|
| Strong PPT | `UUUUUUUUUUUUUUUUCAGGG` — continuous pyrimidines |
| Weak context | `AAAAAAAAAAAAAAAACAGGG` — all-purine region (0.5171875); `AAGAAG…CAGGG` is eliminated (AG in −14..−5, S&S 1987) |
| Expected | Strong PPT score > weak context score |

### Dataset 4: Branch-point detection (yUnAy) — Gao et al. (2008)

**Source:** hand-derived from the Gao (2008) consensus + per-position conservation.

The conservation-weighted branch-point score is `matched / maxScore` where the
informative positions contribute their conservation fractions:
`max = 0.790 (y@−3) + 0.746 (U@−2) + 0.923 (A@0) + 0.751 (y@+1) = 3.210`.

| Parameter | Sequence (5'→3') | AG G-index | Branch A index / distance | Motif | Expected score | PPT fraction |
|-----------|------------------|------------|---------------------------|-------|----------------|--------------|
| Perfect yUnAy | `UU CUUAC U×22 AG GGG` | 30 | 5 / 25 nt | `CUUAC` | `3.210/3.210 = 1.0` | 1.0 (all U) |
| No branch A | `UU CUUUC U×22 AG GGG` | 30 | — (best lacks A@0) | — | `(0.790+0.746+0+0.751)/3.210 = 0.712461` → < 0.8 ⇒ not found | — |
| Purine at −3 | `GG AUUAC G×22 AG GGG` | 30 | 5 / 25 nt | `AUUAC` | `(0+0.746+0.923+0.751)/3.210 = 0.753894` | 0.0 (all G) |
| Near edge 18 nt | `GG CUUAC G×15 AG GGG` | 23 | 5 / 18 nt | `CUUAC` | 1.0 ⇒ found | — |
| Just inside 17 nt | `GG CUUAC G×14 AG GGG` | 22 | 5 / 17 nt | — | < 18 nt ⇒ not found | — |
| Far edge 40 nt | `GG CUUAC G×37 AG GGG` | 45 | 5 / 40 nt | `CUUAC` | 1.0 ⇒ found | — |
| Just outside 41 nt | `GG CUUAC G×38 AG GGG` | 46 | 5 / 41 nt | — | > 40 nt ⇒ not found | — |

### Dataset 5: MaxEntScan score3ss worked examples — Yeo & Burge (2004) / maxentpy

**Source:** the maxentpy `score3` docstring (documented reference values); reproduced this session.

| 23-nt window | Expected (2 dp) | Reproduced (full precision) |
|--------------|-----------------|------------------------------|
| `ttccaaacgaacttttgtAGgga` | **2.89** | 2.886773 |
| `tgtctttttctgtgtggcAGtgg` | **8.19** | 8.190965 |
| `ttctctcttcagacttatAGcaa` | **-0.08** | -0.080278 |

The 2.89 value is the canonical example from the MaxEntScan website / maxentpy README; it is the
primary cross-check for the factorisation + embedded tables.

**Exhaustive oracle cross-check (B11 review, 2026-11-25; SPLICE-MAXENT3-001).** Opened this session:
`raw.githubusercontent.com/kepbod/maxentpy/master/maxentpy/maxent.py` + `.../data/score3_matrix.txt`
(byte-identical to the embedded `Data/maxent_score3.txt`, sha256 9e8a74dc…b3ae) + `LICENSE` (MIT), and the
CPAN port of the original score3.pl `Bio-RNA-SpliceSites-Scoring-MaxEntScan-0.05` (`MaxEntScan.pm::score3`,
`me2x3acc.pm` = the original me2x3acc1..9 tables), run under Perl 5.38. Table check: all 82 560 entries of
the 9 sub-matrices present in both, numerically and textually identical. 200 025 23-mers scored (seeded
random: 100 000 uniform ACGT + 100 000 with AG forced at 19-20, half in mixed case; all 16 dinucleotides on
the `ttccaaacgaacttttgt..GGA` backbone; the docstring examples; homopolymers):

| reference | 2-dp (as score3.pl prints) | full precision |
|-----------|----------------------------|----------------|
| maxentpy `score3` | 200 025 / 200 025 identical | 200 025 / 200 025 bit-identical (after F22; before: 62 349 differed by 1 ulp, max abs Δ 7.1e-15, `Math.Log2` vs ln/ln2 — incl. `tgtctttttctgtgtggcAGtgg` 8.19096548716202 vs 8.190965487162021) |
| score3.pl (CPAN port, me2x3acc1-9) | 200 025 / 200 025 identical | max abs Δ 1.4e-14 (Perl multiplies the 5 numerators, then divides by the product of the 4 denominators) |

Worked examples (score3.pl printed / maxentpy): `ttccaaacgaacttttgtAGgga` 2.89 / 2.8867730651152104;
`tgtctttttctgtgtggcAGtgg` 8.19 / 8.190965487162021; `ttctctcttcagacttatAGcaa` −0.08 / −0.08027774948780869;
range over the sample −61.11 (`TCAAGCAGAGGAGAAATGCTTGA`) … 13.97 (`CTTTCCCGGTCTCCCCGCAGGGC`). Non-AG
acceptors are scored, not rejected (`…gtACGGA` −5.18, `…gtGGGGA` −5.07). Edge cases vs reference: case folded
(maxentpy `.upper()`, score3.pl `uc`); wrong length → maxentpy `sys.exit('Wrong length of fa!')`, CPAN port
`invalid_length`, C# `ArgumentException`; N / R / `-` → maxentpy `ValueError`/`KeyError`, CPAN port
`invalid_alphabet`, C# `ArgumentException`; U accepted as T (C# extension; both references reject U).
`CalculateMaxEntScore(motif, Acceptor)` (MCP `maxent_score`) previously returned a PWM log-weight sum
(`TTCCAAACGAACTTTTGTAGGGA` → −26.81 instead of 2.89; any length/alphabet accepted), now delegates to
`ScoreAcceptorMaxEnt` (F21).

---

### Dataset 6: Shapiro & Senapathy 3' ss score, U12 / branch-point models — B11 review oracle (2026-11-25)

**Opened this session:** S&S 1987 (above); Gao et al. 2008 PDF (europepmc getPdf PMC2367711, `pdftotext`: Abstract
"Yeast carries a strictly conserved branch point sequence (BPS) of UACUAAC"; Table 4 nucleotide frequencies at the 181
branch sites, positions −5..+3 — A .254 .232 .083 .066 .166 .923 .182 .302 .201; C .210 .227 .470 .160 .448 .033 .331
.274 .391; G .254 .193 .127 .028 .177 .017 .066 .112 .112; U .282 .348 .320 .746 .210 .028 .420 .313 .296 — confirms the
0.790/0.746/0.923/0.751 weights; branch points −50..−5, 83 % at −34..−21, PPT +4..+24); Turunen et al. 2013 WIREs RNA
(europepmc getPdf PMC3584512: U12 introns lack a PPT, BPS-to-3'ss optimum 11–13 nt, AT-AC and GT-AG subtypes);
Levine thesis ch. 4 (ftp.sanger.ac.uk/pub/resources/theses/levine/chapter4.pdf = Levine & Durbin 2001: branch site
TTCCTTAA, branch-site-to-acceptor distance 8–21 bp); WebSearch snippets of Brock, Dietrich & Padgett 2008 (PMC2578861:
"fully consensus branch site sequence (UUCCUUAAC)", branch = second A, NOL1 intron F A90) and of Mount's RNAinfo page
(U12 5'ss RTATCCTY, branch TCCTRAY; no U12 3'ss consensus). No source supporting the former `YCCAC` U12 3'ss
consensus was found (WebSearch; Dietrich et al. 2001: penultimate A essential, terminal C most active — no extended
3'ss motif).

| check | result |
|---|---|
| S&S Table 5 published 3' ss scores (41 sites: primate 6, rodent 6, other mammal 6, other vertebrate 6, invertebrate 8, plant 9), formula above, AG rule | 40/41 identical at 1 dp; rabbit MHC intron 3 `AATCTTTCCCTCAG/A` 86.79 vs published 86.9 (Table 1 percentages are rounded) |
| Independent Python oracle (S&S 8 subtables + AG elimination, maxentpy score3 for MaxEnt mode, U12 UUCCUUAAC 8–21 nt, yUnAy / UACUAAC scanners, FindAcceptorBranchPoint incl. PPT fraction; confidences, motifs, windows, N/IUPAC/`-`, case, U, both overloads) vs C# | 8,104 inputs (8,000 random 0–70 nt with planted motifs + HBB/HBA1/TP53/CFTR genomic (Ensembl REST, GRCh38) × {primate, total, plant, maxent, yUnAy, yeast} + the 40 real acceptors × 2 thresholds): 1,174,365 records (76,918 Acceptor, 161,875 U12Acceptor, 934,141 Branch, 1,431 branch-point results), 0 mismatching inputs; Score / Confidence / PPT fraction bit-identical |
| Real human acceptors (Ensembl canonical, 40) | primate S&S 0.717–0.992 (mean 0.880; 1 eliminated by an AG in −14..−5); MaxEnt 3.21–13.37 bits (mean 9.19); `yUnAy` ≥ 0.5 found 18–40 nt upstream for 40/40; fraction of all AG with S&S ≥ 0.5: 0.36–0.55 (42–62 % of AG eliminated / unscorable) |

## Verified Design Decisions

All previously documented assumptions have been audited and resolved:

1. **Default acceptor score — Shapiro & Senapathy (1987) 3' ss score (B11 F27, 2026-11-25):** the former 8-offset
   "PWM" (weights not in any table of the paper; −15/−10 rows invented) + separate 12-nt PPT count + ad-hoc
   `(score/(count+1)+2)/4` normalisation was replaced by the published Table 1 / Methods score (Dataset 6).

2. **Normalization:** none beyond the published formula (score/100).

3. **U12 acceptor scoring — U12 branch site (B11 F29):** the former `YCCAC` weights (1, 1, 0.5) + PPT fraction /3.5
   had no source (no `YCCAC` consensus found; U12 introns have no PPT — Turunen 2013). Replaced by the best
   `UUCCUUAAC` match fraction with the branch A 8–21 nt upstream (Brock 2008; Levine & Durbin 2001), AC or AG terminus.

4. **Branch-point detection (`FindAcceptorBranchPoint`) — sourced from Gao et al. (2008):** the opt-in branch-point detector scans the 18–40 nt window upstream of the acceptor AG for the human `yUnAy` consensus (positions −3..+1, branch A at 0). The per-position conservation weights (y@−3 = 0.790, U@−2 = 0.746, A@0 = 0.923, y@+1 = 0.751) are the Gao (2008) lariat-RT-PCR frequencies; the 18–40 nt envelope brackets the Gao −34..−21 core and is corroborated by Mercer et al. (2015) "19–35 nt". The PPT-fraction window (4–24 nt downstream of the branch point) is also from Gao (2008). **No assumption — all constants are source-traceable.** The acceptor score itself (`FindAcceptorSites`, S&S 1987) has no branch-point term. `FindBranchPoints` (formerly an unsourced YNYURAC PWM) now uses the same yUnAy scorer, plus the yeast UACUAAC model (B11 F30).

5. **MaxEntScan score3ss maximum-entropy acceptor model — IMPLEMENTED (opt-in `ScoreAcceptorMaxEnt`):** the Yeo & Burge (2004) 23-nt maximum-entropy 3' model is now bundled. The precomputed probability tables were retrieved this session from the **MIT-licensed maxentpy port** (`score3_matrix.txt`, 82 560 records) and embedded as `Data/maxent_score3.txt` (provenance + MIT licence in `Data/maxent_score3.LICENSE.md`). The factorisation is the maxentpy `score3` (see source #6 above). The implementation reproduces all three documented worked examples EXACTLY (`2.89`/`8.19`/`-0.08`, full precision `2.886773`/`8.190965`/`-0.080278`), so a wrong table or factorisation would fail the 2.89 cross-check. The original Burge-lab `me2x3acc*` Perl files (academic terms) are NOT used — the bundled artifact is the MIT-licensed equivalent. `CalculateMaxEntScore(motif, Acceptor)` (MCP `maxent_score`) delegates to `ScoreAcceptorMaxEnt` (B11 F21, 2026-11-25; formerly a PWM log-weight sum), and the default `FindAcceptorSites` PWM+PPT scoring is unchanged; `ScoreAcceptorMaxEnt` is a new, additive, opt-in method.

---

## Recommendations for Test Coverage

1. **MUST Test:** Canonical AG acceptor detected — Evidence: Shapiro & Senapathy (1987)
2. **MUST Test:** No AG returns empty — Evidence: trivially correct
3. **MUST Test:** Empty/null input returns empty — Evidence: guard behavior
4. **MUST Test:** No full −14..+1 window (e.g. 14 nt) returns empty — Evidence: S&S scores whole subsequences
5. **MUST Test:** Strong PPT scores higher than weak PPT context — Evidence: Burge et al. (1999)
6. **MUST Test:** Score in [0, 1] range — Evidence: normalization formula
7. **MUST Test:** Confidence in [0, 1] range — Evidence: CalculateConfidence contract
8. **MUST Test:** DNA T↔U equivalence — Evidence: implementation T→U conversion
9. **MUST Test:** Case insensitivity — Evidence: implementation ToUpperInvariant
10. **MUST Test:** Multiple AG sites discovered — Evidence: scanning algorithm
11. **SHOULD Test:** U12 acceptor (UUCCUUAAC branch site 8–21 nt upstream) detected with includeNonCanonical=true — Evidence: Brock et al. (2008), Levine & Durbin (2001)
12. **SHOULD Test:** U12 acceptors excluded with includeNonCanonical=false — Evidence: default parameter
13. **SHOULD Test:** Position reported after AG (i+1) — Evidence: implementation spec
14. **SHOULD Test:** Motif non-empty and contains AG context — Evidence: implementation
15. **COULD Test:** Threshold filtering — higher minScore produces fewer results
16. **COULD Test:** AG with < 12 upstream nt not detected — Evidence: S&S window −14..+1

---

## References

1. Shapiro MB, Senapathy P (1987). RNA splice junctions of different classes of eukaryotes: sequence statistics and functional implications in gene expression. Nucleic Acids Research 15(17):7155–7174.
2. Burge CB, Tuschl T, Sharp PA (1999). Splicing precursors to mRNAs by the spliceosomes. In: The RNA World, 2nd ed. CSHL Press, pp. 525–560.
3. Yeo G, Burge CB (2004). Maximum entropy modeling of short sequence motifs with applications to RNA splicing signals. J Comput Biol 11(2–3):377–394.
4. Patel AA, Steitz JA (2003). Splicing double: insights from the second spliceosome. Nat Rev Mol Cell Biol 4(12):960–970.
5. Lodish H et al. (2004). Molecular Cell Biology, 5th ed. W.H. Freeman.
6. Hall SL, Padgett RA (1994). Conserved sequences in a class of rare eukaryotic nuclear introns with non-consensus splice sites. J Mol Biol 239(3):357–365.
7. Jackson IJ (1991). A reappraisal of non-consensus mRNA splice sites. Nucleic Acids Res 19(14):3795–3798.
8. Dietrich RC, Incorvaia R, Padgett RA (1997). Terminal intron dinucleotide sequences do not distinguish between U2- and U12-dependent introns. Molecular Cell 1(1):151–160.
9. Gao K, Masuda A, Matsuura T, Ohno K (2008). Human branch point consensus sequence is yUnAy. Nucleic Acids Research 36(7):2257–2267. DOI 10.1093/nar/gkn073.
10. Mercer TR, Clark MB, Andersen SB, et al. (2015). Genome-wide discovery of human splicing branchpoints. Genome Research 25(2):290–303. DOI 10.1101/gr.182899.114.
11. Langford CJ, Gallwitz D (1983). Evidence for an intron-contained sequence required for the splicing of yeast RNA polymerase II transcripts. Cell 33:519–527.
12. Brock JE, Dietrich RC, Padgett RA (2008). Mutational analysis of the U12-dependent branch site consensus sequence. RNA 14:2430–2439.
13. Levine A, Durbin R (2001). A computational scan for U12-dependent introns in the human genome sequence. Nucleic Acids Res 29:4006–4013.
14. Turunen JJ, Niemelä EH, Verma B, Frilander MJ (2013). The significant other: splicing by the minor spliceosome. WIREs RNA 4:61–76.

---

## Change History

- **2026-02-12**: Initial documentation.
- **2026-03-16**: Resolved all 3 assumptions. Replaced U12 fixed 0.6 with YCCAC consensus scoring (Hall & Padgett 1994). Verified PWM weights against Shapiro & Senapathy (1987). Documented normalization as design decision. Added references 6-8.
- **2026-06-24**: Added explicit branch-point detection (`FindAcceptorBranchPoint`) sourced from Gao et al. (2008) `yUnAy` consensus + Mercer et al. (2015) distribution (refs 9–10, dataset 4). Default PWM+PPT acceptor scoring unchanged.
- **2026-06-24 (MaxEntScan)**: IMPLEMENTED the Yeo & Burge (2004) MaxEntScan score3ss maximum-entropy 3' acceptor model as the opt-in `ScoreAcceptorMaxEnt`. Embedded the precomputed probability tables (`Data/maxent_score3.txt`, 82 560 records) retrieved verbatim this session from the MIT-licensed maxentpy port; recorded factorisation + provenance + MIT licence (source #6, design decision #5, dataset 5, `Data/maxent_score3.LICENSE.md`). Cross-checked the canonical `score3('...AGgga') == 2.89` plus two further reference values exactly. Default scorers unchanged.
- **2026-11-25 (B11 SPLICE-MAXENT3-001)**: `CalculateMaxEntScore(Acceptor)` / MCP `maxent_score` → `ScoreAcceptorMaxEnt` (F21); score3 log2 as ln/ln2 (F22). Oracle: 200 025 windows bit-identical to maxentpy, 2-dp identical to the score3.pl CPAN port; table bytes = maxentpy, values = original me2x3acc1-9.
- **2026-11-25 (B11 SPLICE-ACCEPTOR-001, F27–F31)**: default acceptor score = Shapiro & Senapathy (1987) 3' ss score
  (Table 1, all eight subtables; primate default; AG elimination); full-window-only scoring; `FindAcceptorSites`
  MaxEntScan mode; U12 acceptors by branch site; `FindBranchPoints` yUnAy / UACUAAC; Dataset 6 oracle.
