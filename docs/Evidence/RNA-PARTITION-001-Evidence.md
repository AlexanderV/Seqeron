# Evidence Artifact: RNA-PARTITION-001

**Test Unit ID:** RNA-PARTITION-001
**Algorithm:** RNA Partition Function (McCaskill) and Boltzmann Structure Probability
**Date Collected:** 2026-06-14

---

## Online Sources

### McCaskill JS (1990) — PubMed record

**URL:** https://pubmed.ncbi.nlm.nih.gov/1695107/
**Accessed:** 2026-06-14 (fetched with WebFetch)
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **Citation:** J. S. McCaskill, "The equilibrium partition function and base pair binding probabilities for RNA secondary structure," *Biopolymers* 1990, **29**(6-7):1105-19. DOI 10.1002/bip.360290621, PMID 1695107.
2. **Partition function definition:** The method computes "the full equilibrium partition function for secondary structure and the probabilities of various substructures."
3. **Complexity:** Both the partition function and the probabilities of all base pairs are computed by "a recursive scheme of polynomial order N³ in the sequence length N."
4. **Base-pair probabilities:** the derived pair binding probabilities are presented in a "box matrix" display over the full ensemble of probable alternative equilibrium structures.

### S. Will, MIT 18.417 (Fall 2011) — McCaskill lecture slides (inside recursion)

**URL:** https://math.mit.edu/classes/18.417/Slides/mccaskill.pdf
**Accessed:** 2026-06-14 (fetched with WebFetch; PDF extracted with `pdftotext`)
**Authority rank:** 1 (graduate course built directly on the primary literature)

**Key Extracted Points:**

1. **Matrices:** `Q_ij := Z(P_ij)` (partition function of sub-ensemble), `Q^b_ij := Z(P^b_ij)` where `P^b_ij = { P ∈ P_ij | (i,j) ∈ P }`.
2. **Q recursion (verbatim):** `Q_ij = Q_{i,j-1} + Σ_{i≤k<j-m} Q_{i,k-1} · Q^b_{kj}`, with `Q_ij = 1` for `i ≥ j - m`.
3. **Total partition function:** "Partition function of the ensemble of S is `Z = Q_{1n}`."
4. **Complexity (verbatim):** "O(n²) space, O(n³) time (after bounding size of interior loops)."
5. **Probabilities (verbatim):** structure `Pr[P|S] = Z⁻¹ exp(−βE(P))`; base pair `Pr[(i,j)|S] = Z⁻¹ Σ_{P∋(i,j)} exp(−βE(P))`, with `β = 1/RT`.
6. **Unambiguous decomposition:** correctness "due to disjoint (=unambiguous) and independent decomposition."

### S. Will, MIT 18.417 (Fall 2011) — McCaskill base-pair probabilities (outside recursion)

**URL:** https://math.mit.edu/classes/18.417/Slides/mccaskill2.pdf
**Accessed:** 2026-06-14 (fetched with WebFetch; PDF extracted with `pdftotext`)
**Authority rank:** 1

**Key Extracted Points:**

1. **External base-pair probability (verbatim):** `p^E_kl = Z(P^E_kl)/Z = (Q_{1,k-1} · Q^b_{kl} · Q_{l+1,n}) / Q_{1n}`, where `(k,l)` is an external base pair (no enclosing pair).
2. **General case (outside recursion):** the full probability sums the external term PLUS, for every pair `(i,j)` that encloses `(k,l)`, the contribution `p_ij · Q^b_kl · (region between (i,j) and (k,l)) / (inside of (i,j))`. This enclosing term does **not** vanish in the flat fixed-per-pair model: a pair `(k,l)` that can be nested inside another pair gets probability strictly larger than its external term. Concretely the correct McCaskill probability is `p_kl = Q^b_kl · O_kl / Z`, where the **outside partition function** obeys
   `O_kl = Q_{1,k-1}·Q_{l+1,n} + Σ_{i<k, j>l, CanPair(i,j), j−i>m} w·Q_{i+1,k-1}·Q_{l+1,j-1}·O_ij`
   with `w = exp(−βE_bp)`. **Verified independently** this session against exhaustive Boltzmann-weighted enumeration on `GGGGCCCC` and `GGGAAACCC` (E_bp = 0 and −1) and on 300 random sequences: the outside recursion matches brute force to machine precision (max error 3.3e-16), whereas the external-only term is wrong for nestable pairs (e.g. `GGGAAACCC` P[2,6] = 6/20 = 0.30, not the external 1/20 = 0.05; `GGGGCCCC` P[1,5] = 3/16 = 0.1875, not 1/16). **CORRECTION (2026-06-16):** an earlier version of this Evidence doc claimed "only the external decomposition contributes, so p_kl = p^E_kl" in the flat model — that is FALSE and matched a bug in the implementation; the outside recursion above is required.

### Freiburg RNA Tools — McCaskill teaching tool (simplified energy model)

**URL:** https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=McCaskill
**Accessed:** 2026-06-14 (fetched with WebFetch)
**Authority rank:** 3 (educational reference implementation, built on McCaskill 1990)

**Key Extracted Points:**

1. **Partition function:** `Z = Σ_P exp(−E(P)/RT)` over all nested (pseudoknot-free) structures P.
2. **Simplified energy model (verbatim idea):** "each base pair of a structure contributes a fixed energy term `E_bp`."
3. **Parameters of the tool:** RNA sequence `S`, minimal loop length `l` (range 0–5), energy weight `E_bp`, normalized temperature `RT`.
4. **Pairing alphabet:** Watson-Crick and GU pairs are treated as complementary.

### ViennaRNA Package — Partition Function and Equilibrium Properties reference

**URL:** https://www.tbi.univie.ac.at/RNA/ViennaRNA/refman/pf_fold.html
**Accessed:** 2026-06-14 (fetched with WebFetch)
**Authority rank:** 3 (canonical reference implementation documentation)

**Key Extracted Points:**

1. **Boltzmann distribution (verbatim):** `p(s) ∝ e^(−βE(s))` with `β = 1/kT`, `k ≈ 1.987 × 10⁻³ kcal/(mol·K)`.
2. **Structure probability (verbatim):** `p(s) = e^(−βE(s)) / Z`.
3. **Default folding temperature** is 37 °C (= 310.15 K), the standard for the Turner/NNDB model.

---

## Documented Corner Cases and Failure Modes

### From MIT 18.417 slides

1. **Minimum loop length:** the recursion is only over `k` with `j − k > m`; pairs with too few enclosed bases are forbidden (`Q^b = 0`), preventing sterically impossible hairpins.
2. **Empty / base interval:** `Q_ij = 1` for `i ≥ j − m` — a sub-sequence too short to contain any pair has exactly one (empty) structure, contributing weight 1.

### From Freiburg tool

1. **Non-pairing alphabet:** only Watson-Crick (A-U, G-C) and GU pairs contribute; all other character pairs give `Q^b = 0`.

---

## Re-validation 2026-11-25 (B12, Stage A + Stage B)

### 8. Stage A verdict: PASS-WITH-NOTES; Stage B verdict: FAIL (fixed)

**Stage A.** McCaskill (1990), the MIT slides and the ViennaRNA reference are consistent with the recursions. The Freiburg tool is a *teaching*
model (fixed energy per pair) and cannot be the authority for a library function that is compared with ViennaRNA; the evidence points above that
call the fixed-`E_bp` model "documented ASSUMPTION" (old §Assumptions 1, Test Datasets "E_bp = 0 ⇒ Z counts structures") are superseded below.
The k constant "≈ 1.987e-3" of the pf_fold page is rounded; the implementation constant is ViennaRNA's `GASCONST` 1.98717 (`params/constants.h`).

**Stage B findings** (reference: ViennaRNA 2.7.2 Python `RNA`, `md.dangles = 2`; NNDB-precision parameter file = `RNA.params_save()` with the NNDB stacks,
TerminalAU 0.45 and the 27 NNDB-vs-default table differences, RNA-STRUCT-001 Evidence §12.2):

| ID | Defect | Repro before → after (oracle) | Source |
|----|--------|-------------------------------|--------|
| F-PF1 | `CalculatePartitionFunction` summed a **fixed energy per base pair** (`basePairEnergy`, default −1) instead of the Turner 2004 loop model; Z, `P(i,j)` and `G` were not the ensemble of the MFE model. | `GGGAAAACCC`: Z 407.70, G −3.704 → Z 32.77, G **−2.1535** (`fc.pf()` −2.1529); MFE structure/energy not reported → `(((....)))` −2.02, P(MFE) 0.805. | ViennaRNA `pf()`; McCaskill 1990 (energy of a structure is the sum of loop energies) |
| F-PF2 | Temperature only changed RT in the fixed-energy model; no Turner parameter was rescaled and the MFE functions had no `temperature` at all. | `GGGGAAAACCCC` MFE at 20 °C: −5.28 (T-independent) → **−7.22** (ViennaRNA −7.20); at 60 °C: −5.28 → **−2.66** (Vienna −2.64); differences = ViennaRNA's int truncation of rescaled terms. | ViennaRNA `params/params.c` `RESCALE_dG`, `default.c` `*_dH`, `intl*dH.h` |
| F-PF3 | No numerical scaling: `Z` is `+∞` and the probability map empty from n ≈ 760 (`Z ≈ 10^(0.40·n)` at `E_bp = −1`); for the Turner model (`ln Z ≈ 0.53 n`) an unscaled recursion would overflow from n ≈ 1300. | n = 1000 random (pre-fix build, measured): Z = ∞, 0 probabilities → ln Z = 535.2, G = −329.9, 186 612 probabilities. | ViennaRNA `pf_scale` (`params.c`: `exp(−sfact·MFE/(RT·n))`) |
| F-PF4 | Base-pair probabilities with the Turner model (`CalculateUnpairedProbabilities`, F7) used one constrained re-fold per pair: O(n⁵). | Release, random GC 0.5: n = 50 2.1 s, n = 100 24.7 s, n = 150 **186 s** → n = 150 0.3 s, n = 1000 6.5 s. | McCaskill outside recursion |
| F-PF5 | `CalculateStructureProbability` formed `exp(−E/RT)/exp(−G/RT)` → ∞/∞ = **NaN** for \|G\| ≳ 440 kcal/mol. | `(−2000, −2000.5)`: NaN → 0.4442 (= exp(−0.5/RT)). | ViennaRNA `pr_structure`: `exp((dG − e)/kT)` |
| F-PF6 | Result lacked ensemble free energy, unpaired probabilities, MFE structure and probability, ensemble diversity, centroid. | added (`PartitionFunctionResult`). | ViennaRNA `pf/bpp/mfe/pr_structure/mean_bp_distance/centroid` |
| F-PF7 | Gas constant 1.987 (rounded) vs ViennaRNA 1.98717 (relative 8.6e-5 in RT): the 0.0028 kcal/mol residual recorded in B12 F7. The NNDB loop extrapolations keep their published 1.75·1.987·T/1000 prefactor. | mean \|ΔG\| vs `fc.pf()` over 320 sequences 0.0004; the remaining max 0.0165 is Vienna's `SMOOTH()` quirk (§9). | `params/constants.h` |

No MCP wrapper existed for the partition function; `rna_partition_function` was added (additive) and `temperature` was added to `minimum_free_energy`,
`rna_mfe_structure`, `rna_structure_energy` (docs `docs/mcp/tools/analysis/`).

### 9. Oracle cross-check (ViennaRNA 2.7.2, d2)

Script `pfcmp3.py` (scratchpad): random ACGU sequences of length 10–150, our `CalculatePartitionFunction` vs `fc.pf()` / `fc.bpp()` / `fc.mfe()` /
`fc.pr_structure(our MFE structure)` / `fc.mean_bp_distance()` / `fc.centroid()`; `params_load(file)` followed by an explicit `md` passed to `fold_compound` (the
`params_load` gotcha).

| Set | N | max / mean \|ΔG_ens\| | max / mean \|ΔP_unpaired\| | max \|ΔP(i,j)\| | \|ΔMFE\| | MFE structure equal | centroid equal |
|-----|---|------|------|------|------|------|------|
| NNDB file, 37 °C, no repeated neighbours | 300 | 0.0165 / 0.0004 | 0.0065 / 0.0002 | 0.0067 | 0 | 298 / 300 | 300 / 300 |
| NNDB file, 37 °C, with repeated neighbours (NNDB-only terms Vienna lacks: special GU closure, all-C, bulge degeneracy / special C) | 320 | 3.04 / 0.30 | — | — | 2.2 / 0.19 | 268 / 320 | 182 / 320 |
| default Vienna file, 37 °C, no repeated neighbours | 320 | 0.54 / 0.107 | 0.15 / 0.035 | 0.15 | 0.58 / 0.12 | 279 / 320 | 242 / 320 |
| NNDB file, 5 / 20 / 55 / 80 / 90 °C | 100 each | 0.0068 / 0.0075 / 0.0058 / 0.0081 / 0.0040 | 0.0057 / 0.0029 / 0.0023 / 0.0045 / 0.0035 | 0.0057 / 0.0027 / 0.0023 / 0.0040 / 0.0029 | 0.25 / 0.15 / 0.23 / 0.17 / 0.13 (max; mean 0.03–0.10) | | |
| NNDB file, 15 / 45 / 70 °C | 150 each | 0.0051 / 0.0080 / 0.0015 | 0.0021 / 0.0049 / 0.0013 | 0.0021 / 0.0048 / 0.0014 | | | |

The default-file residual is parameter rounding (ViennaRNA rounds the NNDB stacks to 0.1 kcal/mol, 26 interior mismatches and one hairpin mismatch differ from NNDB,
TerminalAU 0.50 vs 0.45); the repeated-neighbour residual are NNDB terms with sourced values that ViennaRNA does not implement (RNA-STRUCT-001 Evidence §12).
**Exactness** is established independently of ViennaRNA by exhaustive enumeration (28 sequences ≤ 25 nt × 3 temperatures, every structure scored by
`CalculateFreeEnergy`): Z, G, MFE and all `P(i,j)` agree to 1e-9 (unit test `CalculatePartitionFunction_EqualsExhaustiveEnumeration`).

**Residual of ViennaRNA's `pf()` vs its own `eval_structure`.** Single-structure partition functions (hard constraints) differ from `eval_structure` by
−0.017 kcal/mol for exterior stems with a dangle at the sequence end whose dangle energy is 0 (`GGAAAAC` pair (1,6) with a 5' G/U: pf 4.4830, eval 4.5000): `pf` weights dangles with
`SMOOTH()` (`params/params.c`: "dangling ends should never be destabilizing … specific heat needs smooth function"). The exact Boltzmann sum over `eval_structure` energies — what McCaskill's Z
is — is implemented; the ≤ 0.017 difference is the largest `ΔG_ens` in the table above.

**Probability normalisation (item 3).** 40 random sequences 10–120 nt: every row of the base-pair-probability matrix sums to ≤ 1 (1 + 2e-16), `P(unpaired) + Σ P(pairs) = 1` to 1e-9,
`P(MFE structure) = exp((G − E_mfe)/RT)` to 1e-12 and equals `fc.pr_structure(mfe_struct)` within 2.7 % (max, mean 0.06 %; the Vienna value uses its own `eval`
energy and `SMOOTH`ed Z).

### 9.1 Temperature (item 1)

`ΔG(T) = ΔH − (ΔH − ΔG₃₇)·T/310.15` for every tabulated parameter with the ViennaRNA 2.7.2 Turner 2004 enthalpies (generated into
`Turner2004EnthalpyTables.cs` from `default.c`, `intl11dH.h`, `intl21dH.h`, `intl22dH.h`); `lxc`-extrapolations scaled by T; hairpin first-mismatch with its
terminal-AU enthalpy folded as ViennaRNA does (`mismatchH` includes the AU penalty); dangles / exterior / multiloop terminal mismatches clamped ≤ 0
(`params.c`). Loop-level check (single hairpins 3–40 nt × 6 closing pairs, 19 bulge/interior sizes × 36 pair combinations, exterior stems with all
neighbours) against single-structure ViennaRNA partition functions at 20 and 90 °C: hairpins 0.0000 (> 30 nt: 0.0015 = NNDB 1.75RT extrapolation vs
ViennaRNA's 30-nt anchor), interior loops 0.0000, bulges/GU-closure: only the NNDB-only terms.
Before the dangle/mismatch clamp and the hairpin-AU enthalpy fix the 100-sequence ensemble error at 20 °C was 0.46 (max, mean 0.15) and at 90 °C 0.11 (mean 0.044); after: ≤ 0.0081 at every temperature 5–90 °C (table above).
**Kept at 37 °C (no enthalpy source):** the NNDB-only terms special G·U closure (−2.2), all-C hairpin penalties (+1.5 / 0.3·n + 1.6), special-C bulge (−0.9) and the GGUC/CUGG three-stack (−4.12) have no ΔH in ViennaRNA 2.7.2 (which omits them) nor in the opened NNDB hairpin / bulge / GU pages (ΔG°37 only); the n = 1 bulge degeneracy −RT·ln(states) follows RT(T).

ViennaRNA's `fc.mfe()` at T ≠ 37 °C uses **int-truncated** rescaled parameters (`params.c` `(int)` assignment) — |ΔMFE| mean 0.05–0.09 (max 0.26) vs our unrounded MFE;
its `pf()` does not (`pf_smooth`), which our `G` matches.

### 9.2 Speed and stability (item 2), Release, one thread

| n | time (MFE + inside + outside) | ln Z | Z | pairs listed | min P(unpaired) |
|---|------|------|---|------|------|
| 150 | 0.2 s | 84.4 | 4.3e36 | 4 330 | 0.0011 |
| 200 | 0.23 s | 108.5 | 1.4e47 | 7 452 | 0.0018 |
| 500 | 1.3 s | 296.3 | 4.6e128 | 48 001 | 0.0001 |
| 1000 | 6.5 s | 535.2 | 2.8e232 | 186 612 | 0.0001 |
| 2000 | 52 s | 1030.7 | +∞ (beyond e^709; ln Z finite) | 732 580 | 0.0001 |

ViennaRNA 2.7.2 `mfe()+pf()+bpp()` (C, one thread, other random sequences): n = 500 0.9 s, 1000 4.6 s, 2000 29.3 s. No overflow, NaN or negative probability for any n;
Σ_pairs P + Σ P_unpaired = n to 1e-9 (n = 500 unit test, n = 1000 / 2000 heavy-tier tests). Memory: nine n×n double matrices during the outside pass (≈ 290 MB at n = 2000).

### 9.3 Parameter file notes

NNDB values are used where ViennaRNA's default file departs from NNDB (RNA-STRUCT-001 Evidence §12.3): the default-file ensemble free energy therefore differs by ≈ 0.1 kcal/mol (mean).

### 9.4 Locked ViennaRNA values (NNDB-precision file)

| Sequence | T (°C) | MFE structure | MFE | G (pf) | P(MFE) | ⟨d⟩ |
|----------|--------|---------------|-----|--------|--------|-----|
| GGGAAAACCC | 37 | `(((....)))` | −2.02 | −2.1529 | 0.8060 | 1.3195 |
| GGGCGAAAGCCC | 37 | `(((....)))`-type hairpin | −6.64 | −6.6562 | 0.9740 | 0.0520 |
| GGGGAAAACCCC | 37 | `((((....))))` | −5.28 | −5.3941 | 0.8310 | 1.7016 |
| GCGGAAACGCCGAAAGGCGC | 37 | `.......((((....))))` | −6.84 | −6.8736 | 0.9470 | 0.1876 |
| 27-nt interior/bulge stem | 37 | `((((((((((((...))))))))))))` | −23.57 | −23.5964 | 0.9580 | 0.0842 |
| 40-nt three-hairpin multiloop | 37 | `((((....((((..((((....))))..))))....))))` | −25.02 | −25.0805 | 0.9065 | 2.1116 |
| 76-nt random (no repeated neighbours) | 37 | multi-domain | −21.71 | −23.1022 | 0.1045 | 21.488 |
| 30-nt AU-rich | 37 | `((((((((.(((....))).))))))))..` | −3.95 | −4.2444 | 0.6203 | 3.1347 |
| further rows at 20 and 60 °C | | | | | | see the unit test table |

### 10. Sources opened (2026-11-25)

ViennaRNA 2.7.2 sdist (`pip download viennarna==2.7.2 --no-binary :all:`): `src/ViennaRNA/params/params.c` (`RESCALE_dG`, `RESCALE_BF`, `SMOOTH`, dangle / mismatch clamps,
`pf_scale`: lines 25–90, 400–420, 560–600, 740–830), `params/default.c` (all `*_dH`), `params/intl11dH.h`, `intl21dH.h`, `intl22dH.h`, `params/constants.h` (`GASCONST 1.98717`),
`probabilities/probs_structures.c` (`pr_structure`, `mean_bp_distance`), `structures/centroid.c`, `model.c` (`sfact 1.07`).

### 11. Generators (scratchpad, not in the repository)

`gen_enthalpy.py` — parses `PUBLIC int name[..] = {..};` arrays (INF = 10⁷, comments stripped) of `default.c` / `intl*dH.h` and emits `Turner2004EnthalpyTables.cs`
(Vienna layout, `int` spans, special tri/tetra/hexaloop strings with ΔH and ΔG₃₇); `mkpar.py` (NNDB-precision parameter file); `pfcmp3.py` (cross-check above);
`brute2.py` / `loops.py` / `ends.py` / `evcmp2.py` (single-structure and single-loop ViennaRNA partition functions via hard constraints `hc_add_from_db(db.replace('.', 'x'), CONSTRAINT_DB_DEFAULT | CONSTRAINT_DB_ENFORCE_BP)`).

---

## References

1. McCaskill, J. S. (1990). The equilibrium partition function and base pair binding probabilities for RNA secondary structure. *Biopolymers* 29(6-7):1105-1119. https://doi.org/10.1002/bip.360290621 (PMID 1695107)
2. Will, S. (2011). McCaskill algorithm lecture slides, MIT 18.417 Fall 2011. https://math.mit.edu/classes/18.417/Slides/mccaskill.pdf and https://math.mit.edu/classes/18.417/Slides/mccaskill2.pdf
3. Freiburg RNA Tools. McCaskill teaching tool. https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=McCaskill
4. ViennaRNA Package. Partition Function and Equilibrium Properties (pf_fold reference). https://www.tbi.univie.ac.at/RNA/ViennaRNA/refman/pf_fold.html
5. Nussinov, R., & Jacobson, A. B. (1980). Fast algorithm for predicting the secondary structure of single-stranded RNA. *PNAS* 77(11):6309-6313. https://doi.org/10.1073/pnas.77.11.6309 (minimum-loop and nesting conventions)
6. Mathews, D. H., Disney, M. D., Childs, J. L., Schroeder, S. J., Zuker, M., Turner, D. H. (2004). Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure. *PNAS* 101:7287-7292 (Turner 2004); NNDB https://rna.urmc.rochester.edu/NNDB/turner04/
7. ViennaRNA 2.7.2 source distribution (PyPI `viennarna-2.7.2.tar.gz`): `src/ViennaRNA/params/params.c`, `default.c`, `intl11dH.h`, `intl21dH.h`, `intl22dH.h`, `constants.h`, `probabilities/probs_structures.c`, `structures/centroid.c`; Lorenz R et al. (2011) *Algorithms Mol Biol* 6:26.

---

## Change History

- **2026-06-14**: Initial documentation (RNA-PARTITION-001).
- **2026-06-16**: Validation session. Corrected the base-pair-probability evidence: the prior claim that the external decomposition suffices in the flat model is FALSE for nestable pairs; documented the McCaskill outside recursion `P[i,j] = Q^b·O/Z`, added exact probability datasets (E_bp=0 and −1) for all pairs incl. nested ones, and the per-base ≤1 invariant. Implementation was fixed to use the outside recursion (validation report RNA-PARTITION-001).
- **2026-11-25**: B12 re-validation (Stage A PASS-WITH-NOTES, Stage B FAIL → fixed). The fixed-energy-per-pair model was replaced by the Turner 2004 McCaskill partition function (inside + outside recursion, `pf_scale`, ViennaRNA enthalpy-based temperature rescaling), the result type extended (ensemble free energy, unpaired probabilities, MFE structure / probability, diversity, centroid), `CalculateStructureProbability` evaluated in the exponent, GASCONST 1.98717; §8–§11 added, old fixed-model datasets and assumptions superseded.
