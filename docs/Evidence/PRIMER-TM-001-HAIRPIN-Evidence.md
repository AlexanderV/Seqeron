# Evidence Artifact: PRIMER-TM-001-HAIRPIN

**Test Unit ID:** PRIMER-TM-001 (hairpin / secondary-structure Tm extension)
**Algorithm:** DNA self-folding hairpin (stem + loop) MFE detection + unimolecular hairpin Tm
**Date Collected:** 2026-06-25

---

## Online Sources

### SantaLucia & Hicks (2004) — "The Thermodynamics of DNA Structural Motifs", Annu Rev Biophys Biomol Struct 33:415–440

**URL:** https://users.cs.duke.edu/~reif/courses/molcomplectures/DNA.Thermodynamics&Kinetics/Annu._Rev._Biophys._Biomol._Struct._2004_SantaLucia_Jr.pdf
**Accessed:** 2026-06-25 (PDF fetched, saved locally, text extracted with `pdftotext`)
**Authority rank:** 1 (peer-reviewed review article)

**Key Extracted Points:**

1. **DNA Watson-Crick NN parameters (Table 1):** the unified set — AA/TT = −7.6 kcal/mol ΔH°, −21.3 e.u. ΔS°,
   −1.00 kcal/mol ΔG°37; AT/TA = −7.2/−20.4/−0.88; TA/AT = −7.2/−21.3/−0.58; CA/GT = −8.5/−22.7;
   GT/AC = −8.4/−22.4; CT/AG = −7.8/−21.0; GA/TC = −8.2/−22.2; CG/CG = −10.6/−27.2; GC/GC = −9.8/−24.4;
   GG/CC = −8.0/−19.9. (Identical to the repo's `NnUnifiedParams` reused for the stem stacks.)
2. **Hairpin loop ΔG°37 increments (Table 4 "Hairpin loops" column, 1 M NaCl), verbatim by loop size:**
   3 → 3.5; 4 → 3.5; 5 → 3.3; 6 → 4.0; 7 → 4.2; 8 → 4.3; 9 → 4.5; 10 → 4.6; 12 → 5.0; 14 → 5.1;
   16 → 5.3; 18 → 5.5; 20 → 5.7; 25 → 6.1; 30 → 6.3 (all kcal/mol).
3. **Loop ΔH°/ΔS° rule (Table 4 footnote a, verbatim):** "All loop ΔH° parameters are assumed to equal zero.
   The loop ΔS° increment may be calculated from: ΔS° = ΔG°37 × 1000/310.15." (Because the loop is
   destabilising, ΔG°37 > 0, so ΔS° = −ΔG°37·1000/310.15 once the sign of ΔG° = ΔH° − T·ΔS° is respected.)
4. **Steric minimum (verbatim):** "Hairpin loops with lengths shorter than 3 are sterically prohibited."
5. **Complete hairpin = stem stacks + loop (verbatim, p.428, Eq.10 discussion):** "To compute the stability of a
   complete hairpin + stem, one simply adds the salt-corrected base pair NN contributions (Table 1; Equation 3)
   to the loop energy from Equations 8–10." For loops of length ≥ 5, Eq.10: ΔG°37(total) =
   ΔG°37(Hairpin of N) + ΔG°37(terminal mismatch).
6. **Unimolecular hairpin Tm (Eq.11, verbatim):** "TM = ΔH° × 1000/ΔS° − 273.15, where hairpin loop ΔH° and
   ΔS° are computed with equations analogous to Equations 8–10." This is the two-state hairpin TM; it carries
   NO R·ln(C_T/x) strand-concentration term (contrast the bimolecular duplex Eq.3 which does).
7. **Length-3 / length-4 special handling (Eqs. 8–9):** length-3 hairpins add a triloop bonus + a closing-AT
   penalty (+0.5 kcal/mol, applied only when the loop is closed by an A·T pair); length-4 hairpins add a
   tetraloop bonus + a terminal-mismatch increment. The triloop/tetraloop bonus values and the
   terminal-mismatch increment are in the Annual-Reviews supplementary material (not in the article body).
8. **Jacobson-Stockmayer large-loop extrapolation (Eq.7, verbatim):** ΔG°37(loop-n) =
   ΔG°37(loop-x) + 2.44 × R × 310.15 × ln(n/x), where x is the longest tabulated loop ≤ n; the coefficient
   2.44 is from DNA kinetics measurements (preferred over the older 1.75).

### SantaLucia (1998) — PNAS 95:1460 unified NN set (reuse + cross-check)

**URL:** (NN table also cross-checked against Biopython `DNA_NN4` in the prior PRIMER-TM-001-NN evidence)
**Accessed:** 2026-06-25 (reused; the Table 1 values above are identical in SantaLucia 1998 and SantaLucia & Hicks 2004)
**Authority rank:** 1

**Key Extracted Points:**

1. The repo already embeds this exact unified NN ΔH°/ΔS° table (`PrimerDesigner.NnUnifiedParams`, validated
   under PRIMER-TM-001-NN). The hairpin folder reuses it verbatim for the stem stacks — no new NN values added.

### Search corroboration of concentration-independence of unimolecular hairpin transitions

**URL:** https://pubmed.ncbi.nlm.nih.gov/10423551/ (Vallone & Benight, melting studies of short DNA hairpins)
**Accessed:** 2026-06-25 (web search snippet)
**Authority rank:** 1

**Key Extracted Points:**

1. Hairpin melting temperatures are concentration-independent over strand concentrations 0.5–260 µM; the
   transitions are unimolecular. This independently confirms Eq.11's lack of a concentration term — an
   intramolecular hairpin Tm does not include R·ln(C_T/x).

---

## Documented Corner Cases and Failure Modes

### From SantaLucia & Hicks (2004)

1. **Loop < 3 nt:** sterically prohibited — no hairpin may close a loop of 0–2 nt.
2. **Length-3 / length-4 loops:** require the supplementary triloop/tetraloop bonus + (for length-4) terminal
   mismatch; these supplementary tables are NOT in the article body and are NOT bundled here. Without them the
   computed hairpin is the stem-stack + loop-initiation core (exact and sourced); the bonus is exposed as an
   opt-in caller-supplied increment.
3. **No complementary stem:** a homopolymer (e.g. poly-A) or any oligo with no Watson-Crick stem of ≥ 2 bp
   closing a ≥ 3-nt loop forms no hairpin → no structure / not computable.
4. **Bimolecular vs unimolecular:** the bimolecular duplex-initiation term (+0.2/−5.7) physically nucleates two
   separate strands and is excluded from a unimolecular hairpin; the loop-initiation term is the nucleation cost.

---

## Test Datasets

### Dataset: Hand-derived canonical hairpin GGGCTTTTGCCC (4-bp stem GGGC/GCCC, 4-nt TTTT loop)

**Source:** Derived from SantaLucia & Hicks (2004) Table 1 stem stacks + Table 4 hairpin-loop-of-4 increment.

| Parameter | Value |
|-----------|-------|
| Stem (5' arm) | GGGC; NN steps GG, GG, GC |
| Stem ΔH° | (−8.0) + (−8.0) + (−9.8) = −25.8 kcal/mol |
| Stem ΔS° | (−19.9) + (−19.9) + (−24.4) = −64.2 cal/(K·mol) |
| Loop size | 4 → ΔG°37 = 3.5 kcal/mol (Table 4) |
| Loop ΔH° | 0.0 (Table 4 footnote) |
| Loop ΔS° | −3.5 × 1000 / 310.15 = −11.28486216346929 cal/(K·mol) |
| Total ΔH° | −25.8 kcal/mol |
| Total ΔS° | −64.2 + (−11.28486216346929) = −75.48486216346927 cal/(K·mol) |
| Total ΔG°37 | −25.8 − 310.15·(−75.48486216346927)/1000 = −2.3883700000000054 kcal/mol |
| Hairpin Tm (Eq.11) | (−25.8 × 1000)/(−75.48486216346927) − 273.15 = 68.6403836682880 °C |

### Dataset: 5-nt loop hairpin GGGCAAAAAGCCC (4-bp stem, 5-nt AAAAA loop)

**Source:** Table 1 stem stacks + Table 4 hairpin-loop-of-5 increment (ΔG°37 = 3.3).

| Parameter | Value |
|-----------|-------|
| Stem ΔH° / ΔS° | −25.8 / −64.2 (same GGGC stem) |
| Loop ΔG°37 | 3.3 (size 5) |
| Loop ΔS° | −3.3 × 1000/310.15 = −10.6397... |
| Total ΔS° | −74.83968... |
| Total ΔG°37 | −2.5883700000000... |
| Hairpin Tm | 71.585... °C |

### Dataset: Non-hairpin poly-A (AAAAAAAAAAAA)

**Source:** SantaLucia & Hicks (2004) — no Watson-Crick stem possible.

| Parameter | Value |
|-----------|-------|
| Hairpin found | none (null) |
| Hairpin Tm | NaN |

---

## Assumptions

1. **ASSUMPTION: Bimolecular duplex-initiation term excluded from the hairpin.** The paper states the hairpin =
   "base pair NN contributions" + loop energy and does not add the +0.2/−5.7 bimolecular initiation; this term
   is physically a two-strand nucleation cost (the loop initiation is the unimolecular nucleation cost).
   The literature/UNAFold convention is consistent: a hairpin stem contributes only its NN stacks. This is the
   standard interpretation, not an invented value.
2. **ASSUMPTION: terminal-AT penalty not applied at the open stem end of a hairpin.** The article body's
   hairpin equations (8–10) add only stem NN stacks + loop; the terminal-AT penalty is a duplex-end term in
   Eq.3. It is therefore omitted for the hairpin core (kept exact and sourced); applying it is a refinement
   that would need the supplementary terminal handling.

---

## Recommendations for Test Coverage

1. **MUST Test:** the hand-derived GGGCTTTTGCCC hairpin returns ΔH° = −25.8, ΔS° = −75.48486216346927,
   ΔG°37 = −2.3883700000000054, and Tm = 68.6403836682880 °C exactly — Evidence: Table 1 + Table 4 derivation.
2. **MUST Test:** the folder FINDS the hairpin (stem length 4, loop size 4, stem span covering the whole oligo),
   not a worse partial structure — Evidence: MFE definition + Table 4.
3. **MUST Test:** poly-A returns no hairpin (null) and NaN Tm — Evidence: no Watson-Crick stem possible.
4. **MUST Test:** loop ΔS° follows ΔS° = −ΔG°37·1000/310.15 (a sign error or T error must fail) — Evidence:
   Table 4 footnote a.
5. **MUST Test:** the hairpin Tm has NO concentration term (Eq.11 Tm equals ΔH°·1000/ΔS° − 273.15 regardless of
   any strand concentration) — Evidence: Eq.11 + Vallone & Benight concentration-independence.
6. **SHOULD Test:** a 5-nt loop hairpin (GGGCAAAAAGCCC) uses the Table 4 size-5 increment (3.3) — Rationale:
   exercises a second loop size.
7. **SHOULD Test:** null / empty / non-ACGT input → null / NaN — Rationale: documented invalid input.
8. **COULD Test:** the Jacobson-Stockmayer extrapolation returns Table 4 values exactly at tabulated sizes and a
   monotonically increasing ΔG°37 between them — Rationale: Eq.7 boundary behaviour.

---

## References

1. SantaLucia J Jr, Hicks D (2004). The Thermodynamics of DNA Structural Motifs. Annu Rev Biophys Biomol Struct
   33:415–440. https://doi.org/10.1146/annurev.biophys.32.110601.141800
   (PDF: https://users.cs.duke.edu/~reif/courses/molcomplectures/DNA.Thermodynamics&Kinetics/Annu._Rev._Biophys._Biomol._Struct._2004_SantaLucia_Jr.pdf)
2. SantaLucia J Jr (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor
   thermodynamics. Proc Natl Acad Sci USA 95:1460–1465. https://doi.org/10.1073/pnas.95.4.1460
3. Vallone PM, Benight AS (1999). Melting studies of short DNA hairpins: influence of loop sequence and
   adjoining base pair identity on hairpin thermodynamic stability. https://pubmed.ncbi.nlm.nih.gov/10423551/

---

## Change History

- **2026-06-25**: Initial documentation (hairpin folding + unimolecular hairpin Tm extension of PRIMER-TM-001).
