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

## Full ntthal dimer DP (non-contiguous optima) — implemented 2026-06-25

The dimer alignment now ships the **complete** Primer3 `ntthal` oligo–oligo DP (mode ANY, `type==1`),
a verbatim port of `thal.c` (`fillMatrix`, `LSH`, `RSH`, `maxTM`, `calc_bulge_internal`, `traceback`,
`calcDimer`), so the most stable dimer can be **non-contiguous** — internal single mismatch, internal
loop, single/multi-base bulge, or terminal overhang/dangling end. Method:
`PrimerDesigner.CalculateDimerThermodynamicsNtthal`; `CalculateDimerMeltingTemperature` /
`CalculateSelfDimerMeltingTemperature` delegate to it. (The contiguous `FindMostStableDimer` scorer is
retained unchanged for its existing contiguous-WC results.)

### Recurrence (thal.c, retrieved this session — URL in References [3])

`EnthalpyDPT(i,j)` / `EntropyDPT(i,j)` = best ΔH°/ΔS° of a duplex closing at the WC pair `(i,j)`
(strand1 5′→3′ index `i`, reversed-strand2 index `j`). Built by `fillMatrix` (thal.c L1581-1629):
- **init** (`initMatrix` L1547): WC pair → (0, MinEntropy=−3224); non-pair → (∞, −1).
- **terminal/dangling open** `LSH(i,j)` (L1741): A·T-penalty + `tstack2[s2j][s2j-1][s1i][s1i-1]`, vs
  dangling-end alternatives (`dangle3`/`dangle5`), pick the higher-Tm.
- **stack extension** `maxTM(i,j)` (L1662): compare current vs `(i-1,j-1)` + `stack[..]`.
- **loops/bulges** `calc_bulge_internal(ii,jj,i,j)` (L2149), inner pair `(ii,jj)`:
  - bulge (one side 0 unpaired): size-1 → `bulge[ls] + stack[s1i][s1ii][s2j][s2jj]`; larger → `bulge[ls] + 2·A·T-pen`.
  - 1×1 internal mismatch → `stackmm[..] + stackmm[..]` (both terminal mismatches).
  - general internal loop → `interior[ls] + tstack[..] + tstack[..] + ILAS·|Δn|` (ILAS=−300/310.15 cal/K/mol; ILAH=0).
- **best terminal pair** over all `(i,j)` minimising ΔG (L710-723), then `RSH` closes the 3′ end.
- **N from traceback** (L2957) → `Tm = ΔH/(ΔS + N·saltCorrection + RC) − 273.15`,
  `saltCorrection = 0.368·ln(mv/1000)` (L1042), `RC = R·ln(dna_conc/x)`, x∈{1e9,4e9}.

### Parameter tables (verbatim from primer3 `primer3_config/*.dh,*.ds`)

`stack`, `stackmm` (internal-mismatch), **`tstack2`** (terminal-stacking, the previously-missing table),
`tstack` (`tstack_tm_inf.ds`/`tstack.dh`, internal-loop terminal), `dangle` (5′/3′), and the
interior/bulge loop-length parameters (`loops.dh`/`loops.ds`, lengths 1..30) are embedded in
`NtthalDimer.cs` exactly as `getStack`/`getStackint2`/`getTstack2`/`getTstack`/`getDangle`/`getLoop` read them.

### primer3-py 2.3.0 non-contiguous reference numbers (mv=50, dv=0, dntp=0, dna_conc=50 nM)

| Case | structure | ΔH (kcal/mol) | ΔS (cal/K/mol) | ΔG°37 (kcal/mol) | Tm (°C) |
|------|-----------|---------------|----------------|------------------|---------|
| GCGCATGCGC self | 2×2 internal loop | −84.40000 | −233.42187 | −12.00421 | 43.1572 |
| GCGCAAAGCGC / GCGCTTTGCGC | 3×3 internal loop | −92.30000 | −256.82429 | −12.64594 | 41.8816 |
| GCGCGCGC / GCGCAGCGC | 1-base bulge | −70.80000 | −205.50701 | −7.06200 | 19.8125 |
| GCGCACGCGC / GCGCTAGCGC | 2×2 internal loop | −68.40000 | −198.31701 | −6.89198 | 18.5604 |
| GCGCGCAAAA / AAAAGCGCGC | terminal overhang | −60.00000 | −165.31215 | −8.72844 | 24.6547 |

The C# port reproduces every value above to machine precision (ΔTm = 0, ΔΔG ≈ 1e-12 cal/mol in the
Python reference port that mirrors the C# DP), and still reproduces all contiguous-WC cases
(regression). Verified against `primer3.thermoanalysis.ThermoAnalysis.calc_homodimer` /
`calc_heterodimer`, primer3-py 2.3.0.

---

## Assumptions

1. **(RESOLVED 2026-06-25) gapless alignment only** — the dimer alignment now implements the full
   `ntthal` DP (internal mismatches, internal loops, bulges, `tstack2` terminal overhangs); see the
   section above. No correctness-affecting assumptions remain for the dimer Tm. The only capability
   not modelled is the optional caller-supplied tri/tetraloop & terminal-mismatch hairpin **bonus
   tables** (a hairpin/monomer feature, not a dimer one).

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
