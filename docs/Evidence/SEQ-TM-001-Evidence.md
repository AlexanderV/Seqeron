# Evidence Artifact: SEQ-TM-001

**Test Unit ID:** SEQ-TM-001
**Algorithm:** Melting Temperature (Wallace rule / Marmur-Doty GC formula, and nearest-neighbor Tm)
**Date Collected:** 2026-06-14

---

## Online Sources

### Biopython — `Bio.SeqUtils.MeltingTemp` source

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py
**Accessed:** 2026-06-14 (fetched the raw source file in this session)
**Authority rank:** 3 (reference implementation in an established bioinformatics library)

**Key Extracted Points:**

1. **Wallace rule (`Tm_Wallace`):** "Tm = 4 degC * (G + C) + 2 degC * (A+T)". Docstring worked example: `>>> mt.Tm_Wallace('ACGTTGCAATGCCGTA')` returns `48.0`. Cited reference in the docstring: "The Wallace rule (Thein & Wallace 1986 …) is often used as rule of thumb for approximate Tm calculations for primers of 14 to 20 nt length."
2. **GC formula (`Tm_GC`):** general form "Tm = A + B(%GC) - C/N + salt correction - D(%mismatch)". Valueset 1 comment: "Tm = 69.3 + 0.41(%GC) - 650/N (Marmur & Doty 1962, J Mol Biol 5: 109-118; Chester & Marshak 1993, Anal Biochem 209: 284-290)". Valueset 2 ("QuikChange"): A,B,C,D = (81.5, 0.41, 675, 1); docstring example `>>> Tm_GC('CTGCTGATXGCACGAGGTTATGG', valueset=2)` → `69.20`.
3. **Nearest-neighbor (`Tm_NN`):** final return line `melting_temp = (1000 * delta_h) / (delta_s + (R * (math.log(k)))) - 273.15`, with `R = 1.987` Cal/(°C·mol). Concentration term `k = (dnac1 - (dnac2 / 2.0)) * 1e-9`; `if selfcomp: k = dnac1 * 1e-9`. Defaults: `dnac1=25, dnac2=25, Na=50, saltcorr=5, nn_table=DNA_NN3`. With `dnac1 = dnac2 = 25 nM`, `k = (25 − 12.5)·1e-9 = 12.5e-9` = (C_T = 50 nM)/4.
4. **NN worked example:** `>>> myseq = Seq('CGTTCCAAAGATGTGGGCATGAGCTTAC')` then `>>> print('%0.2f' % mt.Tm_NN(myseq))` → `60.32`.
5. **Salt correction method 5 (`salt_correction`):** "Correction for deltaS: 0.368 x (N-1) x ln[Na+] (SantaLucia (1998), Proc Natl Acad Sci USA 95: 1460-1465)" — the default used by `Tm_NN`.

### Biopython 1.76 API docs — `Bio.SeqUtils.MeltingTemp`

**URL:** https://biopython.org/docs/1.76/api/Bio.SeqUtils.MeltingTemp.html
**Accessed:** 2026-06-14 (fetched in this session)
**Authority rank:** 3 (reference-implementation documentation)

**Key Extracted Points:**

1. **Wallace docstring + example:** confirms "Tm = 4 degC * (G + C) + 2 degC * (A+T)" and `mt.Tm_Wallace('ACGTTGCAATGCCGTA')` → `48.0`; described as suitable for primers 14–20 nt.
2. **GC docstring + example:** confirms general form "Tm = A + B(%GC) - C/N + salt correction - D(%mismatch)" and `Tm_GC('CTGCTGATXGCACGAGGTTATGG', valueset=2)` → `69.20`.
3. **NN docstring + example:** confirms `Tm_NN(Seq('CGTTCCAAAGATGTGGGCATGAGCTTAC'))` → `60.32`, salt-correction method 5 default.

### SantaLucia (1998) PNAS — unified nearest-neighbor thermodynamics

**URL:** https://www.pnas.org/doi/abs/10.1073/pnas.95.4.1460 (DOI https://doi.org/10.1073/pnas.95.4.1460); equation form additionally confirmed via the search-result abstract excerpt for the same paper.
**Accessed:** 2026-06-14 (publisher PDF returned HTTP 403; the abstract/equation text was retrieved via web search of the paper in this session)
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **Oligonucleotide Tm equation:** Tm = ΔH° / (ΔS° + R·ln(C_T / x)) − 273.15, with R the universal gas constant 1.987 Cal/(°C·mol).
2. **Concentration factor x:** for non-self-complementary duplexes the denominator uses C_T/4 (x = 4); for self-complementary duplexes C_T/1 (x = 1).
3. **Salt correction:** ΔS°(salt) = ΔS°(1 M NaCl) + 0.368·(N−1)·ln[Na+], [Na+] in mol/L (this is Biopython's method 5, above).

### Allawi & SantaLucia (1997) Biochemistry — unified NN parameters (Table 1)

**URL:** parameters confirmed via the Biopython `DNA_NN3` table in `MeltingTemp.py` (fetched above), which the source comments attribute to "Allawi HT, SantaLucia J (1997), Biochemistry 36(34):10581-10594, Table 1". DOI of the paper: https://doi.org/10.1021/bi962590c
**Accessed:** 2026-06-14 (the numeric parameters were read from the retrieved Biopython source; the publisher article was not separately fetched)
**Authority rank:** 1 (peer-reviewed primary paper; parameters via rank-3 reference implementation)

**Key Extracted Points:**

1. **Unified NN ΔH (kcal/mol) / ΔS (cal/(mol·K)) at 1 M NaCl** as encoded in the repository implementation (AA/TT −7.9/−22.2; AT −7.2/−20.4; TA −7.2/−21.3; CA/TG −8.5/−22.7; GT/AC −8.4/−22.4; CT/AG −7.8/−21.0; GA/TC −8.2/−22.2; CG −10.6/−27.2; GC −9.8/−24.4; GG/CC −8.0/−19.9).
2. **Helix-initiation:** term-G·C ΔH/ΔS = 0.1 / −2.8; term-A·T = 2.3 / 4.1, applied once at each terminus.

---

## Documented Corner Cases and Failure Modes

### From Biopython `Bio.SeqUtils.MeltingTemp`

1. **Short oligos (Wallace):** Wallace rule is documented as a rule of thumb only for 14–20 nt primers; for longer sequences the GC or NN formula is preferred.
2. **Nearest-neighbor undefined for < 2 nt:** the NN model sums over overlapping dinucleotides, so a length-0 or length-1 input has no NN step.
3. **Concentration convention:** the NN Tm depends on total strand concentration C_T and whether the duplex is self-complementary (x = 1 vs 4); changing C_T or x changes Tm.

---

## Test Datasets

### Dataset: Published / reference-implementation worked examples

**Source:** Biopython `Bio.SeqUtils.MeltingTemp` (cited above).

| Input | Method | Parameters | Expected Tm (°C) |
|-------|--------|-----------|------------------|
| `ACGTTGCAATGCCGTA` (16 nt, 8 GC) | Wallace 2(A+T)+4(G+C) | — | 2·8 + 4·8 = 48.0 |
| `CGTTCCAAAGATGTGGGCATGAGCTTAC` (28 nt) | Nearest-neighbor (DNA_NN3) | Na = 50 mM, C_T = 50 nM (k = 12.5 nM), method 5 | 60.32 |

### Dataset: Marmur-Doty GC-formula hand derivation

**Source:** Marmur & Doty (1962) classic GC form Tm = 64.9 + 41·(GC − 16.4)/N (as encoded in `ThermoConstants`).

| Input | GC | N | Expected Tm (°C) |
|-------|----|---|------------------|
| `GCGCGCGCGCATATATATAT` | 10 | 20 | 64.9 + 41·(10 − 16.4)/20 = 51.78 |

---

## Assumptions

1. **ASSUMPTION: repository default total strand concentration C_T = 250 nM (vs Biopython's 50 nM).** The repository `CalculateThermodynamics` default `primerConcentration = 2.5e-7` (250 nM) with divisor F = 4 gives C_T/4 = 62.5 nM; Biopython's default `dnac1=dnac2=25 nM` gives k = 12.5 nM. The *formula* is identical (ΔH°/(ΔS° + R·ln(C_T/4)) − 273.15); only the default concentration differs and it is an explicit, documented parameter, not an invented constant. Passing `primerConcentration: 5e-8` reproduces Biopython's 60.32 exactly (verified by independent derivation this session). Not correctness-affecting for the formula; it is a documented default parameter.

---

## Recommendations for Test Coverage

1. **MUST Test:** Wallace rule on a published oligo (`ACGTTGCAATGCCGTA` → 48.0) — Evidence: Biopython `Tm_Wallace`.
2. **MUST Test:** Marmur-Doty GC formula exact value (51.78 for a 20-mer with 10 GC) — Evidence: Marmur & Doty 1962 form in Biopython valueset 1 / `ThermoConstants`.
3. **MUST Test:** Nearest-neighbor Tm reproduces the Biopython reference (`CGTTCCAAAGATGTGGGCATGAGCTTAC` → 60.32/60.3 at Na = 50 mM, C_T = 50 nM) — Evidence: Biopython `Tm_NN`.
4. **MUST Test:** NN ΔH/ΔS/ΔG exact tuple from DNA_NN3 + initiation + salt — Evidence: Allawi & SantaLucia 1997 Table 1; SantaLucia 1998 salt/Tm equations.
5. **MUST Test:** Empty / length-1 → all-zero (NN undefined) — Evidence: NN model defined over dinucleotides.
6. **SHOULD Test:** Higher [Na+] raises Tm (monotonic salt term) — Rationale: SantaLucia 1998 salt correction.
7. **SHOULD Test:** Case-insensitivity and Watson-Crick symmetry (AA == TT) — Rationale: DNA_NN3 symmetry.

> All of these cases are already implemented and passing in the canonical fixture
> `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateThermodynamics_Tests.cs`
> (delivered under SEQ-THERMO-001). This unit is consolidated, not re-implemented — see TestSpec §7.

---

## References

1. Allawi, H. T., & SantaLucia, J. (1997). Thermodynamics and NMR of internal G·T mismatches in DNA. *Biochemistry* 36(34):10581–10594. https://doi.org/10.1021/bi962590c
2. SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS* 95(4):1460–1465. https://doi.org/10.1073/pnas.95.4.1460
3. Marmur, J., & Doty, P. (1962). Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature. *J Mol Biol* 5:109–118. (GC formula, as cited in Biopython `Tm_GC` valueset 1)
4. Thein, S. L., & Wallace, R. B. (1986). The use of synthetic oligonucleotides as specific hybridization probes in the diagnosis of genetic disorders. In *Human Genetic Diseases: A Practical Approach*, 33–50. (Wallace rule, as cited in Biopython `Tm_Wallace`)
5. Cock, P. J. A. et al. Biopython, `Bio.SeqUtils.MeltingTemp` (`Tm_Wallace`, `Tm_GC`, `Tm_NN`, `salt_correction`, `DNA_NN3`). https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py and https://biopython.org/docs/1.76/api/Bio.SeqUtils.MeltingTemp.html (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation. Records that SEQ-TM-001 is a duplicate Registry entry for the two melting-temperature methods already delivered under SEQ-THERMO-001 (`CalculateMeltingTemperature`, `CalculateThermodynamics`); evidence independently re-retrieved this session, the implementation verified conformant, and the unit consolidated rather than re-implemented (see TestSpec §7).
