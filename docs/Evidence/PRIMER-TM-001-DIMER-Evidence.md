# Evidence Artifact: PRIMER-TM-001-DIMER

**Test Unit ID:** PRIMER-TM-001 (self-/hetero-dimer Tm extension)
**Algorithm:** Self-dimer / hetero-dimer (intermolecular) Tm via thermodynamic alignment
**Date Collected:** 2026-06-25

---

## Online Sources

### SantaLucia & Hicks (2004) — A unified view of polymer, dumbbell, and oligonucleotide DNA NN thermodynamics

**URL:** https://doi.org/10.1146/annurev.biophys.32.110601.141800 (Annu Rev Biophys 33:415-440)
**Accessed:** 2026-06-25
**Authority rank:** 1 (peer-reviewed review; primary unified NN parameter set)

**Key Extracted Points:**

1. **Unified NN parameters (Table 1):** the 10 distinct Watson-Crick nearest-neighbour ΔH°/ΔS° at 1 M NaCl
   (already in the repo as `NnUnifiedParams`), plus duplex initiation ΔH°=+0.2 kcal/mol, ΔS°=−5.7 cal/(K·mol);
   terminal A·T penalty ΔH°=+2.2, ΔS°=+6.9 per A·T-closed end; symmetry ΔS°=−1.4 for a self-complementary duplex.
2. **Bimolecular Tm (Eq. 3):** `Tm = ΔH°·1000 / (ΔS° + R·ln(C_T/x)) − 273.15`, R = 1.9872 cal/(K·mol),
   x = 4 for non-self-complementary and x = 1 for self-complementary duplexes.
3. **Entropy salt correction (Eq. 5):** `ΔS°[Na⁺] = ΔS°[1 M] + 0.368·(N/2)·ln[Na⁺]`, N = total phosphates,
   so N/2 = number of NN stacks.

### Untergasser et al. (2012) — Primer3, new capabilities and interfaces

**URL:** https://doi.org/10.1093/nar/gks596 (Nucleic Acids Res 40(15):e115)
**Accessed:** 2026-06-25
**Authority rank:** 1 (peer-reviewed; describes the `ntthal` thermodynamic-alignment engine)

**Key Extracted Points:**

1. **Thermodynamic alignment:** Primer3 2.0+ scores oligo–oligo and oligo–template secondary structures
   (dimers, hairpins) by a thermodynamic alignment over the SantaLucia unified NN model, finding the most
   stable structure (the `ntthal` engine, exposed as `calcHomodimer` / `calcHeterodimer`).

### Primer3 `thal.c` (ntthal reference implementation, primer3-py vendored libprimer3)

**URL:** https://raw.githubusercontent.com/libnano/primer3-py/master/primer3/src/libprimer3/thal.c
**Accessed:** 2026-06-25 (fetched and read this session)
**Authority rank:** 3 (reference implementation source)

**Key Extracted Points:**

1. **Initiation (lines 588-589):** for a bimolecular hybridization `dplx_init_H = 200` (cal), `dplx_init_S = −5.7`.
2. **Terminal A·T penalty (lines 128-129):** `AT_H = 2200.0` (cal), `AT_S = 6.9`, applied at each A·T-closed
   helix end (`atPenaltyH`/`atPenaltyS` in `LSH`/`RSH`).
3. **Entropy salt correction (lines 623-624, 1042):**
   `saltCorrection = 0.368·ln((mv + 120·√max(0, dv − dntp))/1000)`, multiplied by N = number of stacks
   (the comment states "8bp dplx N=7"). With dv=dntp=0 this is `0.368·ln(mv/1000)`.
4. **Bimolecular Tm divisor (lines 590-593):** `RC = R·log(dna_conc/1e9)` when **both** strands are symmetric,
   else `R·log(dna_conc/4e9)`. `dna_conc` is in nM, so /1e9 → C_T (mol/L) with x = 1; /4e9 → x = 4.
   Tm comes from `T = (H + dplx_init_H)/(S + dplx_init_S + RC) − 273.15` form.
5. **Symmetry test (line 2771, `symmetry_thermo`):** a sequence is "symmetric" iff it is an even-length
   reverse-complement palindrome.
6. **Default conditions (lines 829, 844):** `dna_conc = 50` (nM), `temp = 310.15` K, used by `calc_homodimer` /
   `calc_heterodimer`.

### primer3-py 2.3.0 — reference numbers (the `ntthal` engine)

**URL:** https://pypi.org/project/primer3-py/ (installed `pip3 install --user primer3-py` → 2.3.0)
**Accessed:** 2026-06-25
**Authority rank:** 3 (reference-implementation output)

**Key Extracted Points (mv=50, dv=0, dntp=0, dna_conc=50 nM):**

1. `calc_homodimer("GCGCGCGC")` → ΔH=−70800 cal, ΔS=−192.6170, ΔG=−11059.84 cal, Tm=40.0906 °C.
2. `calc_homodimer("ACGTACGTACGT")` → ΔH=−92000, ΔS=−262.6267, Tm=37.6251 °C.
3. `calc_heterodimer("ATCGATCGATCG","CGATCGATCGAT")` → ΔH=−92000, ΔS=−264.7267, Tm=32.6107 °C.
4. `calc_homodimer("CGATCGATCG")` → ΔH=−78800, ΔS=−226.8219, Tm=29.6600 °C (CGATCGATCG is a palindrome → x=1).
5. `calc_homodimer("GCATGC")` → ΔH=−43600, ΔS=−125.8121, Tm=0.6859 °C.
6. `calc_heterodimer("GGGGCCCC","GGGGCCCC")` → ΔH=−57600, ΔS=−157.2170, Tm=29.0150 °C.
7. `calc_heterodimer("TGCATGCATG","CATGCATGCA")` → ΔH=−74100, ΔS=−211.8219, Tm=25.6596 °C (non-palindromic → x=4).
8. `calc_homodimer("AAAAAAAA")` → `structure_found = False`, Tm=0 (no stable self-dimer).

---

## Documented Corner Cases and Failure Modes

### From thal.c / primer3-py

1. **No stable dimer:** a homopolymer such as poly-A (`AAAAAAAA`) self-dimer returns `structure_found = False`
   (no Watson-Crick duplex can form). The corresponding method returns `null` / `NaN`.
2. **Symmetric vs non-symmetric concentration divisor:** x = 1 only when *both* aligned oligos are
   reverse-complement palindromes; otherwise x = 4. A self-dimer of a non-palindromic oligo uses x = 4.
3. **ntthal terminal-stack / overhang extension:** for some sequences (e.g. poly-A overhangs, `ATCGTTAC`/
   `GTAACGAT`) ntthal extends the duplex into terminal mismatches/overhangs with extra `tstack2` terms not
   captured by the plain contiguous-WC NN model; those deviate by a small terminal term (documented limit).

---

## Test Datasets

### Dataset: primer3-py 2.3.0 dimer reference (mv=50 mM, dv=0, dntp=0, dna_conc=50 nM)

**Source:** primer3-py 2.3.0 `calc_homodimer` / `calc_heterodimer` (ntthal); SantaLucia & Hicks (2004) Table 1.

| Pair (strand1 / strand2) | x | ΔH° (cal) | ΔS° (cal/K·mol) | Tm (°C) |
|--------------------------|---|-----------|------------------|---------|
| GCGCGCGC / GCGCGCGC (self) | 1 | −70800 | −192.6170 | 40.0906 |
| ACGTACGTACGT / ACGTACGTACGT (self) | 1 | −92000 | −262.6267 | 37.6251 |
| ATCGATCGATCG / CGATCGATCGAT | 1 | −92000 | −264.7267 | 32.6107 |
| CGATCGATCG / CGATCGATCG (self, palindrome) | 1 | −78800 | −226.8219 | 29.6600 |
| GCATGC / GCATGC (self) | 1 | −43600 | −125.8121 | 0.6859 |
| GGGGCCCC / GGGGCCCC | 1 | −57600 | −157.2170 | 29.0150 |
| TGCATGCATG / CATGCATGCA (non-palindromic) | 4 | −74100 | −211.8219 | 25.6596 |
| AAAAAAAA / AAAAAAAA (self) | — | (no structure) | — | NaN |

### Hand-derived exact values (from SantaLucia & Hicks 2004 Table 1, independent of primer3)

- **GCGCGCGC self-dimer (x=1):** stacks 4·GC(−9.8,−24.4) + 3·CG(−10.6,−27.2), init (+0.2,−5.7), no A·T end.
  ΔH° = −70.8 kcal/mol; ΔS°(no salt) = −184.9; salt = 7·0.368·ln(0.05) = −7.71700633667508;
  ΔS° = −192.61700633667505; Tm(C_T=50 nM, x=1) = 40.09064476882935 °C; ΔG°37 = −11.059835484680235 kcal/mol.
- **TGCATGCATG / CATGCATGCA (x=4):** stacks TG,GC,CA,AT,TG,GC,CA,AT,TG; init (+0.2,−5.7); one A·T end (T).
  ΔH° = −74.1 kcal/mol; ΔS° = −211.8218652900108; Tm(C_T=50 nM, x=4) = 25.659587124835923 °C;
  ΔG°37 = −8.403448480303155 kcal/mol.

---

## Assumptions

1. **ASSUMPTION: gapless alignment only.** The implemented thermodynamic alignment is gapless (no internal
   loops/bulges between WC runs) and scores only maximal contiguous Watson-Crick runs. ntthal additionally
   models internal loops and terminal overhang extension; those extra terms are not reproduced here. This is a
   documented model boundary, not a parameter assumption — every NN/init/penalty/salt constant is source-backed.

---

## Recommendations for Test Coverage

1. **MUST Test:** GCGCGCGC self-dimer ΔH°/ΔS°/Tm exactly (hand-derived and primer3-parity). — Evidence: SantaLucia & Hicks (2004) Table 1; primer3-py 2.3.0.
2. **MUST Test:** TGCATGCATG/CATGCATGCA non-palindromic hetero-dimer (x=4) ΔH°/ΔS°/Tm. — Evidence: hand-derived; primer3-py.
3. **MUST Test:** several primer3-parity cases (ACGTACGTACGT, ATCGATCGATCG/CGATCGATCGAT, CGATCGATCG, GCATGC, GGGGCCCC). — Evidence: primer3-py 2.3.0.
4. **MUST Test:** poly-A self-dimer / fully non-complementary pair → no dimer (null / NaN). — Evidence: primer3-py `structure_found=False`.
5. **MUST Test:** invalid input (null, < 2 bases, non-ACGT) → null / NaN. — Evidence: input-validation contract.
6. **SHOULD Test:** self-dimer convenience method equals the two-argument dimer with the same sequence twice. — Rationale: delegation.
7. **SHOULD Test:** monotonic salt dependence (lower [Na⁺] → lower Tm) for a fixed duplex. — Rationale: invariant of Eq. 5.
8. **COULD Test:** the most-stable selection prefers the higher-Tm of two possible alignments. — Rationale: ntthal max-Tm objective.

---

## References

1. SantaLucia J, Hicks D (2004). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Annu Rev Biophys Biomol Struct 33:415-440. https://doi.org/10.1146/annurev.biophys.32.110601.141800
2. Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG (2012). Primer3 — new capabilities and interfaces. Nucleic Acids Res 40(15):e115. https://doi.org/10.1093/nar/gks596
3. Primer3 `thal.c` (ntthal), primer3-py vendored libprimer3. https://raw.githubusercontent.com/libnano/primer3-py/master/primer3/src/libprimer3/thal.c
4. primer3-py 2.3.0 (`calc_homodimer`, `calc_heterodimer`). https://pypi.org/project/primer3-py/

---

## Change History

- **2026-06-25**: Initial documentation of the self-/hetero-dimer Tm extension under PRIMER-TM-001.
