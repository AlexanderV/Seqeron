# Evidence Artifact: SPLICE-ACCEPTOR-001

**Test Unit ID:** SPLICE-ACCEPTOR-001
**Algorithm:** Acceptor Site Detection
**Date Collected:** 2026-02-12

---

## Online Sources

### 1. Wikipedia — RNA splicing

**URL:** https://en.wikipedia.org/wiki/RNA_splicing
**Accessed:** 2026-02-12
**Authority rank:** 4 (Wikipedia citing primary sources)

**Key Extracted Points:**

1. **3' splice site consensus:** Y-rich-N-C-A-G|G — the acceptor site terminates the intron with an almost invariant AG sequence (citing Molecular Biology of the Cell).
2. **Polypyrimidine tract (PPT):** Upstream of the AG, a region high in pyrimidines (C and U) is located. Further upstream is the branch point.
3. **GU-AG rule:** The major spliceosome splices introns containing GU at the 5' splice site and AG at the 3' splice site. This accounts for >99% of intron splicing (canonical splicing / lariat pathway).
4. **Non-canonical splicing:** When intronic flanking sequences do not follow the GU-AG rule, noncanonical splicing occurs via the minor spliceosome.

---

### 2. Wikipedia — Polypyrimidine tract

**URL:** https://en.wikipedia.org/wiki/Polypyrimidine_tract
**Accessed:** 2026-02-12
**Authority rank:** 4 (Wikipedia citing Lodish et al. 2004, Wagner & Garcia-Blanco 2001)

**Key Extracted Points:**

1. **PPT definition:** Region of pre-mRNA rich in pyrimidine nucleotides (especially uracil), usually 15–20 base pairs long.
2. **PPT location:** Located about 5–40 base pairs before the 3' end of the intron.
3. **U2AF binding:** U2AF65 binds to the PPT; U2AF35 binds at the 3' splice site AG. This promotes spliceosome assembly.
4. **PPT quality:** Strong PPT (continuous pyrimidines) produces robust splicing; weak PPT (interrupted by purines) can be suppressed by PTB binding.

---

### 3. Wikipedia — Spliceosome

**URL:** https://en.wikipedia.org/wiki/Spliceosome
**Accessed:** 2026-02-12
**Authority rank:** 4 (Wikipedia citing Will & Lührmann 2011, Burge et al. 1999)

**Key Extracted Points:**

1. **3' splice site elements:** 3' splice site defined by AG dinucleotide, preceded by polypyrimidine tract (PPT) and branch point sequence (BPS).
2. **Branch point distance:** BPS typically 18–50 nt upstream of the 3' splice site.
3. **Minor spliceosome:** U12-type introns use AT-AC dinucleotides instead of canonical GT-AG (Patel & Steitz 2003).

---

### 4. Shapiro & Senapathy (1987)

**Citation:** Shapiro MB, Senapathy P. "RNA splice junctions of different classes of eukaryotes: sequence statistics and functional implications in gene expression." Nucleic Acids Research 15(17):7155–7174, 1987.
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Nucleotide frequencies at 3' splice site:** Compiled from 3,700+ intron sequences across multiple eukaryotic species.
2. **Position -2 (A):** 100% conserved — the A of the CAG acceptor is nearly invariant.
3. **Position -1 (G):** 100% conserved — the G of the CAG acceptor is nearly invariant.
4. **Position -3:** Strong C preference (~65–70%), followed by U (~20%).
5. **Upstream positions (-5 to -15):** Enriched in pyrimidines (C+U ≈ 70–80%), reflecting the polypyrimidine tract.
6. **Position +1 (first exonic nucleotide):** G is most frequent (~50%).

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
| Expected | At least 1 acceptor site with Acceptor type |

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
| Weak context | `AAGAAGAAGAAGAAGAACAGGG` — purines interrupt PPT |
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

---

## Verified Design Decisions

All previously documented assumptions have been audited and resolved:

1. **PWM weights — verified against Shapiro & Senapathy (1987):** The AcceptorPwm values at key positions match published nucleotide frequencies: position -3 C=0.70 (source: C ~65–70%), position -2 A=1.00 (source: 100%), position -1 G=1.00 (source: 100%), position 0 G=0.50 (source: ~50%). Upstream positions reflect documented pyrimidine enrichment (C+U = 0.80, source: 70–80%). **No assumption — values verified.**

2. **Normalization formula — design decision:** The score normalization `(score/(count+1) + 2) / 4` is a heuristic linear mapping from composite PWM + PPT raw scores to [0, 1]. It preserves monotonic ordering (higher biological signal → higher normalized score) and clamps output to [0, 1]. Behavioral properties verified by tests M5, M6, S5. **Not a biological assumption — implementation design choice.**

3. **U12 acceptor scoring — resolved via YCCAC consensus:** Replaced the previous fixed 0.6 score with proper YCCAC consensus-based scoring per Hall & Padgett (1994) and Jackson (1991). The implementation checks the 3' splice site pattern (Y-C-C-A-C), assessing match quality at positions -3 (Y), -2 (C), -1 (C), plus PPT quality upstream. Normalized to [0, 1]. **No assumption — scoring based on published U12 consensus.**

4. **Branch-point detection (`FindAcceptorBranchPoint`) — sourced from Gao et al. (2008):** the opt-in branch-point detector scans the 18–40 nt window upstream of the acceptor AG for the human `yUnAy` consensus (positions −3..+1, branch A at 0). The per-position conservation weights (y@−3 = 0.790, U@−2 = 0.746, A@0 = 0.923, y@+1 = 0.751) are the Gao (2008) lariat-RT-PCR frequencies; the 18–40 nt envelope brackets the Gao −34..−21 core and is corroborated by Mercer et al. (2015) "19–35 nt". The PPT-fraction window (4–24 nt downstream of the branch point) is also from Gao (2008). **No assumption — all constants are source-traceable.** The default PWM+PPT acceptor scoring (`FindAcceptorSites`) is unchanged and still does not include a branch-point term in the acceptor score itself.

5. **MaxEntScan score3ss maximum-entropy acceptor model — IMPLEMENTED (opt-in `ScoreAcceptorMaxEnt`):** the Yeo & Burge (2004) 23-nt maximum-entropy 3' model is now bundled. The precomputed probability tables were retrieved this session from the **MIT-licensed maxentpy port** (`score3_matrix.txt`, 82 560 records) and embedded as `Data/maxent_score3.txt` (provenance + MIT licence in `Data/maxent_score3.LICENSE.md`). The factorisation is the maxentpy `score3` (see source #6 above). The implementation reproduces all three documented worked examples EXACTLY (`2.89`/`8.19`/`-0.08`, full precision `2.886773`/`8.190965`/`-0.080278`), so a wrong table or factorisation would fail the 2.89 cross-check. The original Burge-lab `me2x3acc*` Perl files (academic terms) are NOT used — the bundled artifact is the MIT-licensed equivalent. The existing `CalculateMaxEntScore` helper remains an explicitly named PWM-based approximation, and the default `FindAcceptorSites` PWM+PPT scoring is unchanged; `ScoreAcceptorMaxEnt` is a new, additive, opt-in method.

---

## Recommendations for Test Coverage

1. **MUST Test:** Canonical AG acceptor detected — Evidence: Shapiro & Senapathy (1987)
2. **MUST Test:** No AG returns empty — Evidence: trivially correct
3. **MUST Test:** Empty/null input returns empty — Evidence: guard behavior
4. **MUST Test:** Short sequence (< 20 nt) returns empty — Evidence: implementation guard
5. **MUST Test:** Strong PPT scores higher than weak PPT context — Evidence: Burge et al. (1999)
6. **MUST Test:** Score in [0, 1] range — Evidence: normalization formula
7. **MUST Test:** Confidence in [0, 1] range — Evidence: CalculateConfidence contract
8. **MUST Test:** DNA T↔U equivalence — Evidence: implementation T→U conversion
9. **MUST Test:** Case insensitivity — Evidence: implementation ToUpperInvariant
10. **MUST Test:** Multiple AG sites discovered — Evidence: scanning algorithm
11. **SHOULD Test:** U12 YCCAC acceptor detected with includeNonCanonical=true — Evidence: Hall & Padgett (1994), Patel & Steitz (2003)
12. **SHOULD Test:** U12 YCCAC excluded with includeNonCanonical=false — Evidence: default parameter
13. **SHOULD Test:** Position reported after AG (i+1) — Evidence: implementation spec
14. **SHOULD Test:** Motif non-empty and contains AG context — Evidence: implementation
15. **COULD Test:** Threshold filtering — higher minScore produces fewer results
16. **COULD Test:** AG at position < 15 not detected — Evidence: implementation scan start

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

---

## Change History

- **2026-02-12**: Initial documentation.
- **2026-03-16**: Resolved all 3 assumptions. Replaced U12 fixed 0.6 with YCCAC consensus scoring (Hall & Padgett 1994). Verified PWM weights against Shapiro & Senapathy (1987). Documented normalization as design decision. Added references 6-8.
- **2026-06-24**: Added explicit branch-point detection (`FindAcceptorBranchPoint`) sourced from Gao et al. (2008) `yUnAy` consensus + Mercer et al. (2015) distribution (refs 9–10, dataset 4). Default PWM+PPT acceptor scoring unchanged.
- **2026-06-24 (MaxEntScan)**: IMPLEMENTED the Yeo & Burge (2004) MaxEntScan score3ss maximum-entropy 3' acceptor model as the opt-in `ScoreAcceptorMaxEnt`. Embedded the precomputed probability tables (`Data/maxent_score3.txt`, 82 560 records) retrieved verbatim this session from the MIT-licensed maxentpy port; recorded factorisation + provenance + MIT licence (source #6, design decision #5, dataset 5, `Data/maxent_score3.LICENSE.md`). Cross-checked the canonical `score3('...AGgga') == 2.89` plus two further reference values exactly. Default scorers unchanged.
