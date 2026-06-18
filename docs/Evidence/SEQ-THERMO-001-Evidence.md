# Evidence Artifact: SEQ-THERMO-001

**Test Unit ID:** SEQ-THERMO-001
**Algorithm:** DNA Duplex Thermodynamics (Nearest-Neighbor ΔH°/ΔS°/ΔG°/Tm)
**Date Collected:** 2026-06-13

---

## Online Sources

### Biopython `Bio.SeqUtils.MeltingTemp` (reference implementation)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py
**Accessed:** 2026-06-13 (fetched raw source file in this session)
**Authority rank:** 3 (reference implementation in an established bioinformatics library)

**Key Extracted Points:**

1. **NN table (DNA_NN3, Allawi & SantaLucia 1997):** retrieved verbatim —
   `init`: (0, 0); `init_A/T`: (2.3, 4.1); `init_G/C`: (0.1, -2.8);
   `AA/TT`: (-7.9, -22.2); `AT/TA`: (-7.2, -20.4); `TA/AT`: (-7.2, -21.3);
   `CA/GT`: (-8.5, -22.7); `GT/CA`: (-8.4, -22.4); `CT/GA`: (-7.8, -21.0);
   `GA/CT`: (-8.2, -22.2); `CG/GC`: (-10.6, -27.2); `GC/CG`: (-9.8, -24.4);
   `GG/CC`: (-8.0, -19.9). (dH in kcal/mol, dS in cal/(mol·K).)
2. **Two-terminus initiation:** the `Tm_NN` source computes
   `ends = seq[0] + seq[-1]; AT = ends.count("A")+ends.count("T"); GC = ends.count("G")+ends.count("C")`
   then adds `init_A/T` × AT and `init_G/C` × GC. Initiation is applied to BOTH the
   first and last base pair, not only the first.
3. **Tm equation:** `Tm = (1000 * delta_h) / (delta_s + R * ln(k)) - 273.15`, with
   gas constant `R = 1.987`, and `k = (dnac1 - (dnac2/2.0)) * 1e-9`. ΔH is multiplied
   by 1000 (kcal→cal) so units match ΔS.
4. **Strand concentration:** default `dnac1 = dnac2 = 25` nM ⇒ `k = 12.5` nM = C_T/4
   when C_T = 50 nM (i.e. division by F = 4 for two equimolar non-self-complementary strands).
5. **Salt correction (default method 5):** `corr = 0.368 * (len(seq) - 1) * math.log(mon)`
   with `mon = Na * 1e-3` ([Na+] in mol/L), added to ΔS.
6. **Worked example (docstring):** `Tm_NN(Seq('CGTTCCAAAGATGTGGGCATGAGCTTAC'))` prints `60.32`
   (defaults: dnac1 = dnac2 = 25 nM, Na = 50 mM).

### MELTING 5 User Guide (EBI; Dumousseau et al. 2012)

**URL:** https://www.ebi.ac.uk/biomodels/tools/melting/melting5-UserGuide.pdf
**Accessed:** 2026-06-13 (fetched PDF, converted to text in this session)
**Authority rank:** 3 (reference implementation / curated tool documentation)

**Key Extracted Points:**

1. **Tm equation (§4.2):** "Tm = ΔH / (ΔS + R·ln(C_T/F)) − 273.15", Tm in K for [Na+] = 1 M.
2. **Concentration factor F (§4.3, option `-F`):** "F is 1 in the case of self-complementary
   oligonucleotides… F is 4 if both strands are present in equivalent amount and 1 if one
   strand is in excess. The default factor value is 4." When C_max ≈ C_min, "C_max − C_min/2
   is equivalent to C_T/4, which is the default correction."
3. **Default DNA NN model:** "all97 (from Allawi and Santalucia 1997) (by default)".
4. **Helix-coil model (§4.1):** total ΔH/ΔS are the sum of each Crick's-pair (nearest-neighbor)
   contribution plus a model-specific initiation; for all97 the initiation parameters are added
   for each end.

### Wikipedia — Nucleic acid thermodynamics (cites SantaLucia 1998 as primary)

**URL:** https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics
**Accessed:** 2026-06-13 (fetched in this session)
**Authority rank:** 4 (Wikipedia citing primaries; used only to corroborate, not as sole authority)

**Key Extracted Points:**

1. **Per-terminus initiation:** "Terminal A/T base pair" and "Terminal G/C base pair" are listed
   as separate initiation contributions in Table 1 — both duplex ends receive an initiation term.
   Terminal A/T ΔH° = 9.6 kJ/mol (= 2.29 kcal/mol ≈ 2.3); Terminal G/C ΔH° = 0.4 kJ/mol
   (= 0.096 kcal/mol ≈ 0.1) — corroborates the DNA_NN3 `init_A/T`/`init_G/C` values.
2. **Tm formula:** Tm = ΔH° / (ΔS° + R·ln(C_T/x)); x = 4 for non-self-complementary, x = 1
   for self-complementary; primary reference cited is SantaLucia (1998), reference 13.

---

## Documented Corner Cases and Failure Modes

### From Biopython MeltingTemp

1. **Length < 2:** the nearest-neighbor model requires at least one dinucleotide step; a
   single base or empty input has no NN contribution and no defined duplex.
2. **Lowercase / mixed case:** sequences are processed case-insensitively (Biopython upper-cases
   the sequence); the repository implementation upper-cases via `ToUpperInvariant`.

### From MELTING 5 User Guide

1. **Self-complementary vs non-self-complementary:** the factor F changes the C_T term (1 vs 4);
   this unit uses the default non-self-complementary case F = 4.

---

## Test Datasets

### Dataset: Biopython Tm_NN docstring example

**Source:** Biopython `Bio.SeqUtils.MeltingTemp.Tm_NN` docstring (DNA_NN3 / Allawi & SantaLucia 1997).

| Parameter | Value |
|-----------|-------|
| Sequence | CGTTCCAAAGATGTGGGCATGAGCTTAC |
| dnac1 = dnac2 | 25 nM (⇒ k = 12.5 nM = C_T/4 with C_T = 50 nM) |
| Na+ | 50 mM (0.05 M) |
| Expected Tm | 60.32 °C (rounds to 60.3 at one decimal) |

### Dataset: Hand-derived short oligonucleotide (GCGC), repository defaults

**Source:** Derivation from DNA_NN3 NN/init parameters and SantaLucia (1998) salt + Tm equations.

| Parameter | Value |
|-----------|-------|
| Sequence | GCGC (Na+ = 0.05 M, C_T = 250 nM) |
| Init (both ends G/C) | ΔH 2×0.1 = 0.2; ΔS 2×(−2.8) = −5.6 |
| NN steps GC, CG, GC | ΔH (−9.8−10.6−9.8) = −30.2; ΔS (−24.4−27.2−24.4) = −76.0 |
| ΔH° | −30.0 kcal/mol |
| Salt ΔS | 0.368×3×ln(0.05) = −3.307 |
| ΔS° | −5.6 − 76.0 − 3.307 = −84.91 cal/(mol·K) |
| ΔG°₃₇ | −30.0 − 310.15×(−84.91)/1000 = −3.67 kcal/mol |
| Tm | (−30000)/(−84.91 + 1.987·ln(2.5e-7/4)) − 273.15 = −18.6 °C |

---

## Assumptions

<!-- Zero correctness-affecting assumptions: every constant and formula is source-backed. -->

1. **ASSUMPTION: Empty / length-1 input returns all-zero result** — Neither Biopython nor
   MELTING define an NN result for length < 2 (no dinucleotide exists). The repository contract
   returns `(0,0,0,0)` for such inputs. This is an API/edge-case convention, not a thermodynamic
   value; it does not alter any computed quantity for valid (length ≥ 2) input.

---

## Recommendations for Test Coverage

1. **MUST Test:** Biopython worked example Tm = 60.3 °C for CGTTCCAAAGATGTGGGCATGAGCTTAC under
   matching parameters (Na = 0.05 M, C_T = 50 nM) — Evidence: Biopython Tm_NN docstring.
2. **MUST Test:** Hand-derived GCGC → ΔH −30.0, ΔS −84.91, ΔG −3.67, Tm −18.6 (repo defaults) —
   Evidence: DNA_NN3 + SantaLucia (1998) equations.
3. **MUST Test:** Both-terminus initiation — a sequence with A/T at one end and G/C at the other
   (e.g. ATCG) yields ΔH −23.6, ΔS −71.81 — Evidence: Biopython two-end init logic.
4. **MUST Test:** Length-1 and empty input return `(0,0,0,0)` — Evidence: NN model undefined for length < 2.
5. **SHOULD Test:** Case-insensitivity (lowercase equals uppercase) — Rationale: implementation upper-cases input.
6. **SHOULD Test:** Self-complement symmetry of the NN table (AA/TT, CA/TG identical) reflected in equal results.
7. **COULD Test:** Salt-concentration monotonicity (higher [Na+] ⇒ higher Tm) — Rationale: ΔS salt term.

---

## References

1. Allawi HT, SantaLucia J Jr. (1997). Thermodynamics and NMR of internal G·T mismatches in DNA. Biochemistry 36(34):10581–10594. https://doi.org/10.1021/bi962590c
2. SantaLucia J Jr. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. PNAS 95(4):1460–1465. https://doi.org/10.1073/pnas.95.4.1460
3. Biopython `Bio.SeqUtils.MeltingTemp` (master). DNA_NN3 table and `Tm_NN`. https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py
4. Dumousseau M, Rodriguez N, Juty N, Le Novère N. (2012). MELTING, a flexible platform to predict the melting temperatures of nucleic acids. BMC Bioinformatics 13:101. User guide: https://www.ebi.ac.uk/biomodels/tools/melting/melting5-UserGuide.pdf
5. Wikipedia. Nucleic acid thermodynamics (corroborating, cites SantaLucia 1998). https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics

---

## Change History

- **2026-06-13**: Initial documentation.
