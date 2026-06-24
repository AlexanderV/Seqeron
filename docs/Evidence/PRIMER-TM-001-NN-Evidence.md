# Evidence Artifact: PRIMER-TM-001

**Test Unit ID:** PRIMER-TM-001
**Algorithm:** Primer melting temperature — nearest-neighbour salt-corrected Tm (opt-in)
**Date Collected:** 2026-06-24

> This artifact covers the **nearest-neighbour (NN) salt-corrected design Tm** added under
> PRIMER-TM-001. The legacy Wallace / Marmur-Doty Tm and the Primer3 weighted-penalty objective are
> documented in `PRIMER-TM-001-Evidence.md` / `tests/TestSpecs/PRIMER-TM-001-Penalty.md`; they are
> unchanged. The NN unified ΔG°37 table (already used for 3'-end stability) is validated under
> SEQ-THERMO-001; this artifact adds the ΔH°/ΔS° → Tm path and the salt corrections.

---

## Online Sources

### SantaLucia & Hicks (2004) — "The Thermodynamics of DNA Structural Motifs"

**URL:** https://users.cs.duke.edu/~reif/courses/molcomplectures/DNA.Thermodynamics&Kinetics/Annu._Rev._Biophys._Biomol._Struct._2004_SantaLucia_Jr.pdf (Duke mirror of Annu Rev Biophys Biomol Struct 33:415-440)
**Accessed:** 2026-06-24 (fetched the PDF and read pages 419–423 directly)
**Authority rank:** 1 (peer-reviewed review reproducing the SantaLucia 1998 unified parameters)

**Key Extracted Points:**

1. **Table 1 (unified NN parameters at 1 M NaCl), verbatim ΔH° (kcal/mol) / ΔS° (e.u.):**
   AA/TT −7.6 / −21.3 ; AT/TA −7.2 / −20.4 ; TA/AT −7.2 / −21.3 ; CA/GT −8.5 / −22.7 ;
   GT/CA −8.4 / −22.4 ; CT/GA −7.8 / −21.0 ; GA/CT −8.2 / −22.2 ; CG/GC −10.6 / −27.2 ;
   GC/CG −9.8 / −24.4 ; GG/CC −8.0 / −19.9. **Initiation** +0.2 / −5.7 ; **Terminal AT penalty**
   +2.2 / +6.9 (per end closed by an A·T pair) ; **Symmetry correction** 0.0 / −1.4 (self-comp only).
2. **Tm equation (Eq. 3):** `Tm = ΔH° × 1000 / (ΔS° + R × ln(C_T/x)) − 273.15`, R = 1.9872 cal/(K·mol),
   x = 4 for non-self-complementary and x = 1 for self-complementary duplexes; C_T = total molar strand
   concentration.
3. **Worked example (p.419):** non-self-complementary duplex, ΔH° = −43.5 kcal/mol, ΔS° = −122.5 e.u.,
   0.2 mM each strand → `Tm = −43.5×1000/(−122.5 + 1.9872·ln(0.0004/4)) − 273.15 = 35.8 °C`.
4. **Sodium dependence (Eq. 5):** `ΔS°[Na⁺] = ΔS°[1 M NaCl] + 0.368 × N/2 × ln[Na⁺]`, N = total number
   of phosphates in the duplex; ΔH° independent of [Na⁺]; valid 0.05 M ≤ [Na⁺] ≤ 1.1 M, duplexes ≤ 16 bp.
   Worked example: a 6-bp duplex has N = 10 phosphates; ΔG°37[0.115 M] = −5.35 − 0.114×10/2×ln(0.115)
   = −4.12 kcal/mol (reproduced exactly during evidence collection).

### Owczarzy et al. (2004) — "Effects of sodium ions on DNA duplex oligomers"

**URL:** https://pubmed.ncbi.nlm.nih.gov/15035624/ (abstract); equation coefficients cross-confirmed via the Biopython source (method 6) and web search of Biochemistry 43:3537-54.
**Accessed:** 2026-06-24
**Authority rank:** 1 (peer-reviewed; full text paywalled — coefficients via the reference implementation below)

**Key Extracted Points:**

1. **Monovalent (Na⁺) quadratic 1/Tm correction:**
   `1/Tm[Na] = 1/Tm[1 M] + (4.29·f(GC) − 3.95)·1e-5·ln[Na⁺] + 9.40e-6·(ln[Na⁺])²`, f(GC) = GC fraction.
   The Tm vs ln[Na⁺] relationship is non-linear (quadratic) over 69 mM – 1.02 M.

### Owczarzy et al. (2008) — divalent Mg²⁺ correction (Biochemistry 47:5336-53)

**URL:** coefficients taken verbatim from the Biopython reference implementation (salt_correction method 7).
**Accessed:** 2026-06-24
**Authority rank:** 1 (peer-reviewed; coefficients via reference implementation)

**Key Extracted Points:**

1. **Divalent (Mg²⁺/dNTP) 1/Tm correction:** `1/Tm = 1/Tm[1 M] + corr`, with
   `corr = a + b·ln[Mg] + f(GC)·(c + d·ln[Mg]) + (1/(2(N−1)))·(e + f·ln[Mg] + g·ln[Mg]²)`,
   base coefficients (×1e-5) a=3.92, b=−0.911, c=6.26, d=1.42, e=−48.2, f=52.5, g=8.31; in the mixed
   regime (0.22 ≤ R < 6, R=√[Mg]/[Mon]) a, d, g are reparameterised by the published √[Mon]/ln[Mon]
   expressions; free Mg²⁺ is reduced by dNTP chelation via Ka = 3×10⁴.

### Biopython `Bio.SeqUtils.MeltingTemp` — reference implementation

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py
**Accessed:** 2026-06-24 (fetched the raw source)
**Authority rank:** 3 (well-maintained reference library)

**Key Extracted Points:**

1. **`DNA_NN4` (SantaLucia 1998)** verbatim, matching Table 1 above:
   `init (0.2, -5.7)`, `init_A/T (2.2, 6.9)`, `sym (0, -1.4)`, AA/TT (−7.6,−21.3), AT/TA (−7.2,−20.4),
   TA/AT (−7.2,−21.3), CA/GT (−8.5,−22.7), GT/CA (−8.4,−22.4), CT/GA (−7.8,−21.0), GA/CT (−8.2,−22.2),
   CG/GC (−10.6,−27.2), GC/CG (−9.8,−24.4), GG/CC (−8.0,−19.9). **Independent cross-check of Table 1.**
2. **`Tm_NN`:** R = 1.987; `Tm = (1000·ΔH)/(ΔS + R·ln(k)) − 273.15`; selfcomp ⇒ k = C_T and the `sym`
   term applied; else k uses the (dnac1 − dnac2/2) form (≡ x = 4 for equal strands).
3. **`salt_correction` method 5** (SantaLucia Eq. 5 entropy form): `corr = 0.368·(len−1)·ln(mon)`, added
   to ΔS. **method 6** (Owczarzy 2004): `corr = (4.29·gc − 3.95)·1e-5·ln(mon) + 9.40e-6·ln(mon)²`, applied
   as `1/Tm_new = 1/Tm_old + corr`. **method 7** (Owczarzy 2008 divalent) as above.

### SantaLucia & Hicks (2004) Table 2 / Table 3 — mismatch + dangling-end ΔH° (primary cross-check)

**URL:** https://users.cs.duke.edu/~reif/courses/molcomplectures/DNA.Thermodynamics&Kinetics/Annu._Rev._Biophys._Biomol._Struct._2004_SantaLucia_Jr.pdf
**Accessed:** 2026-06-24 (fetched the PDF, ran `pdftotext -layout`, read Tables 2 and 3 directly)
**Authority rank:** 1 (peer-reviewed review reproducing the original mismatch/dangling-end measurements)

**Key Extracted Points:**

1. **Internal single mismatches (Table 2).** The review tabulates only the ΔG°37 increments and states
   verbatim: *"Error bars and ΔH° and ΔS° parameters are provided in the original references."* It gives a
   **worked example** (Eq. 6): the duplex `5'-GGACTGACG-3'` / `3'-CCTGGCTGC-5'` (two internal mismatches)
   has predicted **ΔG°37 = −8.32 kcal/mol** summed as `init + sym + GG + GA + AC + CT + CG + GA + AC + CG`.
   - **Cross-check performed this session:** summing the Biopython `DNA_NN` + `DNA_IMM` ΔH°/ΔS° for this exact
     duplex and converting `ΔG°37 = ΔH° − 310.15·ΔS°/1000` gives **−8.349 kcal/mol** (ΔH° = −58.5, ΔS° = −161.7),
     within rounding of the paper's −8.32 (the paper sums independently-rounded ΔG°37 increments). This confirms
     the Biopython `DNA_IMM` table is a faithful transcription of the Allawi/SantaLucia/Peyret mismatch parameters.
2. **Dangling ends (Table 3).** The review gives **both ΔH° and ΔG°37** for all 32 single dangling ends, in the
   antiparallel `XB/C` notation (`XB/C` = `5'-XB-3'` over `3'-C-5'`, X unpaired). Example ΔH° rows:
   5'-dangling `XA/T`: A=0.2, C=0.6, G=−1.1, T=−6.9; `XC/G`: A=−6.3, C=−4.4, G=−5.1, T=−4.0; …
   3'-dangling `AX/T`: A=−0.5, C=4.7, G=−4.1, T=−3.8; `TX/A`: A=−0.7, C=4.4, G=−1.6, T=2.9; …
   - **Cross-check performed this session:** mapping every Biopython `DNA_DE` key to the SantaLucia Table 3 cell
     (16 `XY/.Z` keys = 5'-dangling; 16 `.X/YZ` keys = 3'-dangling) reproduced **all 32 ΔH° values exactly**
     (`scratchpad/de_check3.py`). The Biopython `DNA_DE` table is a verbatim transcription of Bommarito (2000) /
     SantaLucia & Hicks (2004) Table 3.

### Bommarito, Peyret & SantaLucia (2000) — dangling-end primary reference

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC103285/ (open-access PMC mirror of NAR 28(9):1929-1934)
**Accessed:** 2026-06-24
**Authority rank:** 1 (the primary dangling-end measurement)

**Key Extracted Points:**

1. All 32 single-nucleotide dangling ends measured by UV melting in **1 M NaCl**; ΔG°37 contributions range
   **+0.48 to −0.96 kcal/mol** (Table 2 of the paper; the numeric grid is published as figure images on PMC, so
   the ΔH°/ΔS° values were cross-confirmed against SantaLucia & Hicks (2004) Table 3 ΔH° and the Biopython
   `DNA_DE` table above).

### Biopython `DNA_IMM` / `DNA_DE` + `Tm_NN` mismatch/dangling path — reference implementation

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py
**Accessed:** 2026-06-24 (fetched the raw source; read the `Tm_NN` body verbatim)
**Authority rank:** 3 (well-maintained reference library transcribing the papers above)

**Key Extracted Points:**

1. **`DNA_IMM` references** (verbatim header): *"Allawi & SantaLucia (1997), Biochemistry 36: 10581-10594;
   (1998) 37: 9435-9444, 2170-2179, Nucl Acids Res 26: 2694-2701; Peyret et al. (1999), Biochemistry 38: 3468-3477."*
   The A/C/G/T single-mismatch entries (52 keys) were copied verbatim into `NnInternalMismatch` (the
   Watson-Crick `WC` placeholders and the inosine `I` entries are excluded — only true A/C/G/T mismatches kept).
2. **`DNA_DE` reference** (verbatim header): *"Bommarito et al. (2000), Nucl Acids Res 28: 1929-1934."* All 32
   keys copied verbatim into `NnDanglingEnd`.
3. **`Tm_NN` convention (read verbatim):** `c_seq = Seq(seq).complement()` (NOT reverse complement); the NN/
   mismatch key is `tmp_seq[i:i+2] + "/" + tmp_cseq[i:i+2]`, looked up forward then character-reversed
   (`neighbors[::-1]`); a left dangling end uses `tmp_seq[:2] + "/" + tmp_cseq[:2]` then strips the first column;
   a right dangling end uses `tmp_cseq[-2:][::-1] + "/" + tmp_seq[-2:][::-1]` then strips the last column;
   `ends = seq[0] + seq[-1]` (terminal-AT count uses the **un-dotted** top strand); `sym` applied only for a
   self-complementary duplex. The C# `CalculateNearestNeighborThermodynamicsMismatch` mirrors this exactly.

### OligoPool SantaLucia Tm tutorial (cross-check of equation form)

**URL:** https://oligopool.com/resources/tutorials/calculating-tm
**Accessed:** 2026-06-24
**Authority rank:** 4 (tutorial citing the primaries; used only to corroborate the equation form)

**Key Extracted Points:**

1. Confirms `Tm = ΔH° / (ΔS° + R × ln(Ct/4)) − 273.15 + Salt Correction`, R = 1.987 cal/mol·K, divisor
   4 for non-self-complementary sequences.

---

## Documented Corner Cases and Failure Modes

### From SantaLucia & Hicks (2004)

1. **Self-complementary duplex:** x = 1 (not 4); add symmetry ΔS° = −1.4; no terminal-AT penalty when
   ends are G·C.
2. **Terminal A·T penalty per end:** a duplex with both ends A·T incurs the +2.2/+6.9 term twice.
3. **Salt-correction validity:** Eq. 5 holds for 0.05–1.1 M [Na⁺] and duplexes ≤ 16 bp.

### From Biopython

1. **Non-ACGT base:** no NN parameter → thermodynamics not computable (here returns null / NaN).

---

## Test Datasets

### Dataset: Published worked example (SantaLucia & Hicks 2004, p.419)

**Source:** SantaLucia & Hicks (2004) Annu Rev Biophys 33:415-440, Eq. 3 example.

| Parameter | Value |
|-----------|-------|
| ΔH° | −43.5 kcal/mol |
| ΔS° | −122.5 cal/(K·mol) |
| Strand conc. (each) | 0.2 mM → C_T = 4×10⁻⁴, x = 4 |
| R | 1.9872 cal/(K·mol) |
| **Tm** | **35.8 °C** |

### Dataset: Hand-derived NN sums (from Table 1)

**Source:** Table 1 (above), summed per Eq. 1.

| Oligo | Self-comp? | ΔH° (kcal/mol) | ΔS° (e.u.) | Tm @ 1 M, C_T=0.5 µM (°C) |
|-------|-----------|----------------|------------|---------------------------|
| GCGCGC | yes (x=1) | −50.4 | −134.7 | 35.0473059911 |
| ATGCATGC | no (x=4) | −57.1 | −156.5 | 30.4338060665 |
| CGCGAATTCGCG | yes (x=1) | −100.6 | −272.1 | 61.1452300219 |

| Oligo | Owczarzy 2004 @ 50 mM Na (°C) | SantaLucia Eq.5 @ 50 mM Na (°C) |
|-------|-------------------------------|----------------------------------|
| GCGCGC | 28.1593085080 | 24.9976652723 |
| ATGCATGC | 18.1899960529 | — |

### Dataset: Internal-mismatch + dangling-end hand-derived sums (DNA_IMM / DNA_DE)

**Source:** Allawi/SantaLucia/Peyret `DNA_IMM` + Bommarito (2000) `DNA_DE`, summed by the Tm_NN convention.
Bottom strand written 3'→5' (complement direction). Full arithmetic in `tests/TestSpecs/PRIMER-TM-001-NN.md` §4.

| Duplex (top 5'→3' / bottom 3'→5') | Feature | ΔH° (kcal/mol) | ΔS° (e.u.) | Tm @ no-salt, C_T=0.5 µM (°C) |
|-----------------------------------|---------|----------------|------------|-------------------------------|
| CGTGAC / GCGCTG | one internal T·G mismatch | −35.5 | −101.5 | −6.4060879279 |
| AGCGCGC / .CGCGCG | one 5'-dangling A (on a C·G pair) | −51.9 | −136.4 | 35.8034921829 |
| GCGCGC / CGCGCG | perfect (equivalence check) | −50.4 | −134.7 | 35.0473059911 (= perfect-match path) |

Stack-by-stack (CGTGAC/GCGCTG): `CG/GC`(WC) `GT/CG`(IMM,G·T) `TG/GC`→`CG/GT`rev(IMM) `GA/CT`(WC) `AC/TG`→`GT/CA`rev(WC).
Stack-by-stack (AGCGCGC/.CGCGCG): left-DE `AG/.C`(−3.7,−10.0), then WC `GC/CG`×3, `CG/GC`×2; terminal-AT once (top ends A,C).

---

## Assumptions

1. **ASSUMPTION: Default C_T = 0.5 µM** — a common PCR primer working concentration, exposed as a
   parameter so the caller can override. Affects only the default operating point, not formula correctness.
2. **ASSUMPTION: Eq. 5 phosphate count N = 2·(length − 1)** — confirmed against the paper's 6-bp →
   N = 10 worked example (no 5'-terminal phosphates); the −4.12 kcal/mol example reproduces exactly.

---

## Recommendations for Test Coverage

1. **MUST Test:** ΔH°/ΔS° sums for a non-self-comp and a self-comp oligo match Table 1 term-by-term —
   Evidence: SantaLucia & Hicks 2004 Table 1.
2. **MUST Test:** the published 35.8 °C worked example reproduces via the Tm equation — Evidence: p.419.
3. **MUST Test:** Owczarzy 2004 correction lowers Tm and matches the hand-derived 1/Tm value —
   Evidence: Owczarzy 2004 + Biopython method 6.
4. **SHOULD Test:** SantaLucia Eq. 5 entropy correction; self-comp x=1 vs non-self-comp x=4; Mg²⁺ raises
   Tm; divalent-with-no-Mg falls back to monovalent.
5. **COULD Test:** monotonic Tm vs [Na⁺]; default-Tm-method-unchanged guard; invalid input → NaN.

---

## References

1. SantaLucia J (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. PNAS 95(4):1460-1465. https://doi.org/10.1073/pnas.95.4.1460
2. SantaLucia J, Hicks D (2004). The thermodynamics of DNA structural motifs. Annu Rev Biophys Biomol Struct 33:415-440. https://doi.org/10.1146/annurev.biophys.32.110601.141800
3. Owczarzy R, You Y, Moreira BG, et al. (2004). Effects of sodium ions on DNA duplex oligomers: improved predictions of melting temperatures. Biochemistry 43(12):3537-3554. https://doi.org/10.1021/bi034621r
4. Owczarzy R, Moreira BG, You Y, et al. (2008). Predicting stability of DNA duplexes in solutions containing magnesium and monovalent cations. Biochemistry 47(19):5336-5353. https://doi.org/10.1021/bi702363u
5. Biopython `Bio.SeqUtils.MeltingTemp` (DNA_NN4, DNA_IMM, DNA_DE, Tm_NN, salt_correction). https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/MeltingTemp.py (accessed 2026-06-24)
6. Allawi HT, SantaLucia J (1997). Thermodynamics and NMR of internal G·T mismatches in DNA. Biochemistry 36(34):10581-10594. https://doi.org/10.1021/bi962590c
7. Allawi HT, SantaLucia J (1998). Nearest-neighbor thermodynamic parameters for internal G·A mismatches in DNA. Biochemistry 37(8):2170-2179. https://doi.org/10.1021/bi9724873 ; Thermodynamics of internal C·T mismatches in DNA. Nucleic Acids Res 26(11):2694-2701. https://doi.org/10.1093/nar/26.11.2694 ; Nearest-neighbor parameters for internal A·C mismatches. Biochemistry 37(26):9435-9444. https://doi.org/10.1021/bi9803729
8. Peyret N, Seneviratne PA, Allawi HT, SantaLucia J (1999). Nearest-neighbor thermodynamics and NMR of DNA sequences with internal A·A, C·C, G·G and T·T mismatches. Biochemistry 38(12):3468-3477. https://doi.org/10.1021/bi9825091
9. Bommarito S, Peyret N, SantaLucia J (2000). Thermodynamic parameters for DNA sequences with dangling ends. Nucleic Acids Res 28(9):1929-1934. https://doi.org/10.1093/nar/28.9.1929

---

## Change History

- **2026-06-24**: Initial documentation — NN salt-corrected Tm added (opt-in) under PRIMER-TM-001.
- **2026-06-24**: Added the internal single-mismatch (Allawi/SantaLucia/Peyret) and single-dangling-end
  (Bommarito 2000) NN ΔH°/ΔS° tables and the `*Mismatch` Tm path (opt-in extension). The mismatch table was
  cross-checked against the SantaLucia & Hicks (2004) Table 2 worked example; all 32 dangling-end ΔH° values
  were cross-checked term-by-term against SantaLucia & Hicks (2004) Table 3.
