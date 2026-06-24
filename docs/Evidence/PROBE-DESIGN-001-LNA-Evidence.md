# Evidence Artifact: PROBE-DESIGN-001 (LNA-adjusted nearest-neighbour Tm)

**Test Unit ID:** PROBE-DESIGN-001
**Algorithm:** LNA (locked nucleic acid)-adjusted nearest-neighbour melting temperature; citable MGB design rules
**Date Collected:** 2026-06-24

---

## Online Sources

### McTigue, Peterson & Kahn (2004) — LNA-DNA nearest-neighbour thermodynamics (primary)

**URL:** https://pubs.acs.org/doi/abs/10.1021/bi035976d (DOI 10.1021/bi035976d); PubMed https://pubmed.ncbi.nlm.nih.gov/15122905/
**Accessed:** 2026-06-24
**Authority rank:** 1 (peer-reviewed paper, *Biochemistry* 43(18):5388–5405, 2004)

**Key Extracted Points:**

1. **Scope:** Hybridization ΔH°, ΔS° and Tm were measured from absorbance melting curves for 100 duplex oligonucleotides carrying **single internal LNA nucleotides on one strand**. Results are reported as ΔΔH°, ΔΔS° (and ΔΔG°, ΔTm) for **all 32** possible nearest neighbours for LNA+DNA:DNA hybridization (5'-MX_L and 5'-X_L N, where M, N, X ∈ {A,C,G,T} and X_L is the LNA monomer). (Abstract / PubMed record, retrieved 2026-06-24.)
2. **Direction (stabilization):** "Incorporation of single LNA nucleotides provides substantial increases in duplex stability"; "LNA provides the largest known increase in thermal stability of any modified DNA duplex" — i.e. an internal LNA substitution **raises Tm** of the duplex (effect is sequence-dependent). (Cross-confirmed via Levin et al. 2006, NAR 34(20):e142, ref. 22, retrieved 2026-06-24.)
3. **Application model:** The 32 LNA NN values are **increments (ΔΔH°, ΔΔS°) added to the underlying DNA nearest-neighbour stack** for the step containing the LNA base; the duplex ΔH°/ΔS° are first computed with the standard DNA NN model and the LNA increment is then added per LNA-containing NN step (the MELTING reference implementation realizes exactly this; see below).
4. **Terminal-LNA exclusion:** The parameters are for **internal** LNA substitutions; they are not established for a terminal (5'- or 3'-end) LNA position. (MELTING `McTigue04LockedAcid.isApplicable` rejects terminal-L positions, citing the paper.)

### MELTING 5 reference implementation — verbatim McTigue 2004 NN parameter table

**URL (data file):** GitHub `mohakjain/TmCalculator` → `MELTING5.2.0/Data/McTigue2004lockedmn.xml` (mirrored verbatim in `aravind-j/rmelting`, `bioc/rmelting`, `hrbrmstr/melting5jars` → `inst/extdata/Data/McTigue2004lockedmn.xml`)
**URL (model):** `MELTING5.2.0/src/melting/patternModels/specificAcids/McTigue04LockedAcid.java`
**URL (docs):** https://pmc.ncbi.nlm.nih.gov/articles/PMC3733425/ (Dumousseau et al. 2012, *MELTING, a flexible platform…*, BMC Bioinformatics 13:101)
**Accessed:** 2026-06-24 (fetched via `gh api repos/mohakjain/TmCalculator/contents/...`)
**Authority rank:** 3 (established reference implementation; the XML header states "contains the parameters for DNA/DNA locked acid nucleic of McTigue et al. (2004). Biochemistry 43 : 5388-5405")

**Key Extracted Points:**

1. **Units:** "the enthalpy and entropy values of each parameter are given in cal/mol" (MELTING paper). So the XML enthalpy `992.0` means ΔΔH° = +0.992 kcal/mol; entropy is in cal/(mol·K).
2. **Combination (verbatim from `McTigue04LockedAcid.computeThermodynamics`):**
   `result = computeThermodynamicsWithoutLockedNucleicAcid(...)` (the plain DNA NN sum), then for each NN step `enthalpy += lockedAcidValue.getEnthalpy(); entropy += lockedAcidValue.getEntropy();` — confirming the LNA value is an **additive increment** to the DNA NN ΔH°/ΔS°.
3. **Terminal-LNA rejection (verbatim from `isApplicable`):** if the LNA is at duplex position 0 or the last position, "The thermodynamics parameters for locked nucleic acids of McTigue (2004) are not established for terminal locked nucleic acids" → not applicable.
4. **Notation:** sequence keys are `<NN-with-L>/<complement>` where `L` immediately follows the locked base, e.g. `TTL/AA` = step `TT`, the 3' T is locked; `TLG/AC` = step `TG`, the 5' T is locked.

### MELTING 5 worked example (rmelting tutorial) — exact reproduction target

**URL:** https://aravind-j.github.io/rmelting/articles/Tutorial.html
**Accessed:** 2026-06-24
**Authority rank:** 3 (reference-implementation worked example)

**Key Extracted Points:**

1. Input `melting(sequence = "CCATTLGCTACC", nucleic.acid.conc = 0.0001, hybridisation.type = "dnadna", Na.conc = 1, method.locked = "mct04")` →
   **Tm = 63.61426 °C**, ΔH = −81100 cal/mol, ΔS = −222.5 cal/(mol·K).
2. **Notation confirmation:** "Locked nucleic acid bases are denoted with an 'L' suffix following the nucleotide letter" — so `CCATTLGCTACC` is the 11-nt DNA duplex `CCATTGCTACC` with the 5th base (the second `T`, 0-based index 4) locked.

---

## Documented Corner Cases and Failure Modes

### From McTigue 2004 / MELTING

1. **Terminal LNA not parameterised:** an LNA at the first or last duplex position has no McTigue increment — must be rejected (return not-computable), not silently treated as internal.
2. **Sequence-dependence:** the per-LNA ΔΔ varies widely by NN context (some ΔΔH increments are positive, some negative), but the net effect on Tm is stabilizing for internal LNA; tests must use the actual table values, not a single average.
3. **Non-ACGT base:** the underlying DNA NN lookup fails on a non-ACGT base → not computable.

---

## Test Datasets

### Dataset: MELTING worked example (single internal LNA)

**Source:** rmelting tutorial (MELTING 5 with `mct04`), retrieved 2026-06-24.

| Parameter | Value |
|-----------|-------|
| DNA duplex (top strand 5'→3') | CCATTGCTACC (11 nt) |
| LNA position (0-based) | 4 (the second T) |
| Strand conc. | 1e-4 M |
| [Na⁺] | 1 M (reference state) |
| MELTING Tm | 63.61426 °C |
| MELTING ΔH | −81100 cal/mol |
| MELTING ΔS | −222.5 cal/(mol·K) |

### Dataset: McTigue 2004 LNA NN increments (verbatim, cal/mol and cal/(mol·K))

**Source:** `McTigue2004lockedmn.xml` (MELTING 5), retrieved 2026-06-24. Key = `<NN-with-L>/<complement>`.

| Key | ΔΔH (cal/mol) | ΔΔS (cal/mol·K) | | Key | ΔΔH | ΔΔS |
|-----|------|------|---|-----|------|------|
| ALA/TT | 707 | 2.5 | | AAL/TT | 992 | 4.1 |
| ALT/TA | 2282 | 7.5 | | ATL/TA | 1816 | 6.9 |
| ALG/TC | 264 | 2.6 | | AGL/TC | -1200 | -1.8 |
| ALC/TG | 1131 | 4.1 | | ACL/TG | 2890 | 10.6 |
| TLA/AT | -46 | 1.6 | | TAL/AT | 1591 | 5.3 |
| TLT/AA | 1528 | 5.3 | | TTL/AA | 2326 | 8.1 |
| TLG/AC | -1540 | -3.0 | | TGL/AC | 2165 | 7.2 |
| TLC/AG | 1893 | 6.7 | | TCL/AG | 609 | 3.2 |
| GLA/CT | 3162 | 10.5 | | GAL/CT | 444 | 2.9 |
| GLT/CA | -212 | 0.1 | | GTL/CA | -635 | -0.3 |
| GLG/CC | -2844 | -6.7 | | GGL/CC | -943 | -0.9 |
| GLC/CG | -360 | -0.3 | | GCL/CG | -925 | -1.1 |
| CLA/GT | 1049 | 4.3 | | CAL/GT | 1358 | 4.4 |
| CLT/GA | 708 | 4.2 | | CTL/GA | -1671 | -4.1 |
| CLG/GC | 785 | 3.7 | | CGL/GC | -276 | -0.7 |
| CLC/GG | 2096 | 8.0 | | CCL/GG | 2063 | 7.6 |

---

## Worked example (hand-derived, independent of the implementation)

Duplex `CCATTGCTACC`, LNA at index 4 (second T), at C_T = 1e-4 M, x = 4 (non-self-comp), [Na⁺] = 1 M reference state.

**Base DNA NN (SantaLucia 1998 unified, the library's `NnUnifiedParams`):**
steps CC,CA,AT,TT,TG,GC,CT,TA,AC,CC → ΔH° = init(0.2) + Σstacks + terminal-AT (none: ends C,C) = **−80.8 kcal/mol**; ΔS° = init(−5.7) + Σstacks = **−221.7 cal/(mol·K)**.

**LNA increments for the locked index-4 (T):**
- step index 3 = `TT`, 3' base (index 4) locked → key `TTL/AA` = (+2326 cal, +8.1).
- step index 4 = `TG`, 5' base (index 4) locked → key `TLG/AC` = (−1540 cal, −3.0).

ΔH°_LNA = −80.8 + 2.326 − 1.540 = **−80.014 kcal/mol**; ΔS°_LNA = −221.7 + 8.1 − 3.0 = **−216.6 cal/(mol·K)**.

Tm = ΔH°·1000 / (ΔS° + R·ln(C_T/4)) − 273.15, R = 1.9872:
- **All-DNA (no LNA):** Tm = −80800 / (−221.7 + 1.9872·ln(1e-4/4)) − 273.15 = **59.692 °C**.
- **LNA-adjusted:** Tm = −80014 / (−216.6 + 1.9872·ln(1e-4/4)) − 273.15 = **63.528 °C**.
- **ΔTm = +3.84 °C** → confirms an internal LNA **raises** Tm (McTigue 2004).
- vs MELTING `mct04` (63.614 °C): agrees to **0.086 °C**; the small residual is MELTING pairing the McTigue increments with a *different* base DNA NN set than SantaLucia 1998 unified — the LNA increment application itself is exact.

---

## Assumptions

1. **ASSUMPTION: base DNA NN model = SantaLucia 1998 unified.** The library's LNA Tm adds the McTigue increments to the library's existing SantaLucia-1998-unified DNA NN ΔH°/ΔS° (`CalculateNearestNeighborThermodynamics`), not to McTigue's own reference DNA set. This keeps the LNA Tm consistent with the rest of the library's NN Tm (`CalculateMeltingTemperatureNN`); the ~0.09 °C offset vs MELTING is attributable to this base-model choice (documented), and the **increment values are applied exactly as published**. Non-correctness-affecting for the increment contribution itself; the base model is independently sourced (SantaLucia 1998/2004).

---

## Recommendations for Test Coverage

1. **MUST Test:** ΔH°/ΔS° of `CCATTGCTACC` with LNA at index 4 equals the hand-derived −80.014 kcal/mol / −216.6 cal/(mol·K) — Evidence: McTigue 2004 XML (`TTL/AA`, `TLG/AC`) + SantaLucia base NN.
2. **MUST Test:** LNA-adjusted Tm of that duplex (C=1e-4, Na=1, no extra salt correction) = 63.52759 °C (hand value), and is within ~0.1 °C of the MELTING `mct04` 63.61426 °C — Evidence: rmelting worked example.
3. **MUST Test:** adding the LNA monomer **raises** Tm vs the all-DNA duplex (63.53 > 59.69) — Evidence: McTigue 2004 stabilization.
4. **MUST Test:** a terminal LNA (index 0 or last) is rejected (not-computable / NaN) — Evidence: McTigue/MELTING terminal exclusion.
5. **SHOULD Test:** every one of the 32 increment keys present; spot-check a negative-ΔΔH key (`GLG/CC` = −2844) is applied with correct sign.
6. **SHOULD Test:** null/empty/short/non-ACGT and out-of-range LNA index → not-computable.
7. **COULD Test:** MGB design-rule check — shorter probe length window (13–20 nt) flagged; 3'-MGB placement guidance.

---

## References

1. McTigue PM, Peterson RJ, Kahn JD (2004). Sequence-dependent thermodynamic parameters for locked nucleic acid (LNA)-DNA duplex formation. *Biochemistry* 43(18):5388–5405. https://doi.org/10.1021/bi035976d ; https://pubmed.ncbi.nlm.nih.gov/15122905/
2. Dumousseau M, Rodriguez N, Juty N, Le Novère N (2012). MELTING, a flexible platform to predict the melting temperatures of nucleic acids. *BMC Bioinformatics* 13:101. https://pmc.ncbi.nlm.nih.gov/articles/PMC3733425/
3. MELTING 5 data file `McTigue2004lockedmn.xml` and `McTigue04LockedAcid.java` — GitHub `mohakjain/TmCalculator/MELTING5.2.0/`, mirrored in `aravind-j/rmelting/inst/extdata/Data/`. Retrieved 2026-06-24.
4. rmelting tutorial worked example. https://aravind-j.github.io/rmelting/articles/Tutorial.html (retrieved 2026-06-24).
5. Levin JD, Fiala D, Samala MF, Kahn JD, Peterson RJ (2006). Position-dependent effects of locked nucleic acid (LNA) on DNA sequencing and PCR primers. *Nucleic Acids Res* 34(20):e142 (cites McTigue 2004 ref. 22 for LNA stabilization). https://academic.oup.com/nar/article/34/20/e142/3100503
6. SantaLucia J (1998). PNAS 95(4):1460–65 (base DNA NN model used by the library). https://www.pnas.org/doi/10.1073/pnas.95.4.1460

---

## Change History

- **2026-06-24**: Initial documentation (LNA-adjusted NN Tm; McTigue 2004 increment table retrieved verbatim from MELTING 5; MGB design-rule citability assessed).
