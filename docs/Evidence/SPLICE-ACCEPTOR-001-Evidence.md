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

**Citation:** Yeo G, Burge CB. "Maximum entropy modeling of short sequence motifs with applications to RNA splicing signals." Journal of Computational Biology 11(2–3):377–394, 2004.
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **MaxEntScan acceptor model:** Uses 23-mer context (20 intronic + 3 exonic positions) for scoring 3' splice sites.
2. **Scoring approach:** Maximum entropy model captures dependencies between positions, outperforming simple PWM models.
3. **PPT contribution:** PPT quality (pyrimidine fraction) is a major determinant of acceptor site score.

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

---

## Verified Design Decisions

All previously documented assumptions have been audited and resolved:

1. **PWM weights — verified against Shapiro & Senapathy (1987):** The AcceptorPwm values at key positions match published nucleotide frequencies: position -3 C=0.70 (source: C ~65–70%), position -2 A=1.00 (source: 100%), position -1 G=1.00 (source: 100%), position 0 G=0.50 (source: ~50%). Upstream positions reflect documented pyrimidine enrichment (C+U = 0.80, source: 70–80%). **No assumption — values verified.**

2. **Normalization formula — design decision:** The score normalization `(score/(count+1) + 2) / 4` is a heuristic linear mapping from composite PWM + PPT raw scores to [0, 1]. It preserves monotonic ordering (higher biological signal → higher normalized score) and clamps output to [0, 1]. Behavioral properties verified by tests M5, M6, S5. **Not a biological assumption — implementation design choice.**

3. **U12 acceptor scoring — resolved via YCCAC consensus:** Replaced the previous fixed 0.6 score with proper YCCAC consensus-based scoring per Hall & Padgett (1994) and Jackson (1991). The implementation checks the 3' splice site pattern (Y-C-C-A-C), assessing match quality at positions -3 (Y), -2 (C), -1 (C), plus PPT quality upstream. Normalized to [0, 1]. **No assumption — scoring based on published U12 consensus.**

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

---

## Change History

- **2026-02-12**: Initial documentation.
- **2026-03-16**: Resolved all 3 assumptions. Replaced U12 fixed 0.6 with YCCAC consensus scoring (Hall & Padgett 1994). Verified PWM weights against Shapiro & Senapathy (1987). Documented normalization as design decision. Added references 6-8.
