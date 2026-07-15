---
type: source
title: "Evidence: SEQ-THERMO-001 (DNA duplex thermodynamics — nearest-neighbor ΔH°/ΔS°/ΔG°/Tm via Biopython Tm_NN + DNA_NN3 / Allawi & SantaLucia 1997)"
tags: [validation, thermodynamics, sequence-statistics]
doc_path: docs/Evidence/SEQ-THERMO-001-Evidence.md
sources:
  - docs/Evidence/SEQ-THERMO-001-Evidence.md
source_commit: 41aaf8a4899a8795559ba409d32abc24d36553a1
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-THERMO-001

The validation-evidence artifact for test unit **SEQ-THERMO-001** — **DNA duplex
thermodynamics**: the full nearest-neighbor (NN) **ΔH° / ΔS° / ΔG° / Tm** 4-tuple of a
Watson-Crick DNA duplex, a verbatim port of Biopython `Bio.SeqUtils.MeltingTemp.Tm_NN` over the
**DNA_NN3** parameter set (Allawi & SantaLucia 1997). The synthesis of the algorithm lives on
[[dna-duplex-nearest-neighbor-thermodynamics]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern. This page records only what the
artifact adds (sources, oracles, corner cases).

## Distinct from the two other Tm surfaces in the wiki

This is neither of the previously-ingested Tm engines — it is the SEQ-family **full-tuple** engine:

- **NOT** the legacy Wallace / Marmur-Doty scalar Tm (`CalculateMeltingTemperature`, SEQ-TM-001,
  the length<14 rule-of-thumb bundled by [[seq-summary-001-evidence|SummarizeNucleotideSequence]]) —
  that is a `%GC`-only rule of thumb, not a thermodynamic model.
- **NOT** the PCR-primer NN engine [[primer-dimer-thermodynamics-tm]], which uses the **2004
  unified** SantaLucia & Hicks parameters (`NnUnifiedParams`) with a fixed duplex-initiation
  (+0.2 / −5.7) plus a terminal-A·T penalty (+2.2 / +6.9). This unit uses the **1997 DNA_NN3**
  vintage with a different initiation bookkeeping (per-terminus `init_A/T` / `init_G/C` counted off
  the two terminal bases) and — unlike the primer path — returns the whole **ΔH°/ΔS°/ΔG°** tuple,
  not just a Tm. See the concept page for the side-by-side.

## Online sources

- **Biopython `Bio.SeqUtils.MeltingTemp`** (rank 3, reference implementation) — the port target.
  - **DNA_NN3 table (Allawi & SantaLucia 1997), retrieved verbatim** (ΔH kcal/mol, ΔS cal/(mol·K)):
    `init` (0,0); `init_A/T` (2.3, 4.1); `init_G/C` (0.1, −2.8); `AA/TT` (−7.9,−22.2);
    `AT/TA` (−7.2,−20.4); `TA/AT` (−7.2,−21.3); `CA/GT` (−8.5,−22.7); `GT/CA` (−8.4,−22.4);
    `CT/GA` (−7.8,−21.0); `GA/CT` (−8.2,−22.2); `CG/GC` (−10.6,−27.2); `GC/CG` (−9.8,−24.4);
    `GG/CC` (−8.0,−19.9).
  - **Two-terminus initiation:** `Tm_NN` counts terminal bases —
    `ends = seq[0]+seq[-1]; AT = ends.count("A")+ends.count("T"); GC = ends.count("G")+ends.count("C")`
    — then adds `init_A/T × AT + init_G/C × GC`. Initiation is applied to **both** the first and
    last base pair, not only the first.
  - **Tm equation:** `Tm = (1000·ΔH)/(ΔS + R·ln k) − 273.15`, `R = 1.987`,
    `k = (dnac1 − dnac2/2)·1e−9`. ΔH×1000 converts kcal→cal so units match ΔS.
  - **Strand concentration:** default `dnac1 = dnac2 = 25` nM ⇒ `k = 12.5` nM `= C_T/4` (C_T = 50 nM),
    i.e. the **F = 4** factor for two equimolar non-self-complementary strands.
  - **Salt correction (default method 5, SantaLucia 1998):** `corr = 0.368·(len−1)·ln(mon)` with
    `mon = [Na⁺]·1e−3` (mol/L), added to ΔS.
  - **Docstring worked example:** `Tm_NN(Seq('CGTTCCAAAGATGTGGGCATGAGCTTAC'))` → **60.32 °C**
    (defaults dnac1 = dnac2 = 25 nM, Na = 50 mM).
- **MELTING 5 User Guide** (EBI; Dumousseau et al. 2012; rank 3) — corroborates the Tm equation
  `Tm = ΔH/(ΔS + R·ln(C_T/F)) − 273.15` (§4.2); F factor (§4.3) = 1 self-complementary / 4 both
  strands equimolar (default 4; `C_max − C_min/2 ≡ C_T/4`); default DNA NN model "all97 (Allawi and
  SantaLucia 1997)"; helix-coil model = Σ per-Crick-pair contributions + per-end initiation (§4.1).
- **Wikipedia — Nucleic acid thermodynamics** (rank 4, cites SantaLucia 1998) — corroborates the
  per-terminus initiation (Terminal A·T ΔH° ≈ 2.3, Terminal G·C ΔH° ≈ 0.1 kcal/mol) and the Tm
  formula with x = 4 (non-self-complementary) / x = 1 (self-complementary).

## Datasets (oracles)

- **Biopython Tm_NN docstring:** `CGTTCCAAAGATGTGGGCATGAGCTTAC`, dnac1 = dnac2 = 25 nM (k = 12.5 nM),
  Na = 50 mM → **Tm = 60.32 °C** (rounds to 60.3 at one decimal).
- **Hand-derived `GCGC`** (Na = 0.05 M, C_T = 250 nM): init both ends G/C ⇒ ΔH 2×0.1 = 0.2,
  ΔS 2×(−2.8) = −5.6; NN steps GC+CG+GC ⇒ ΔH −30.2, ΔS −76.0 ⇒ **ΔH° = −30.0**; salt
  ΔS = 0.368×3×ln(0.05) = −3.307 ⇒ **ΔS° = −84.91**; **ΔG°₃₇ = −30.0 − 310.15×(−84.91)/1000 = −3.67**;
  **Tm = (−30000)/(−84.91 + 1.987·ln(2.5e−7/4)) − 273.15 = −18.6 °C** (a valid negative-Tm result for
  a very short duplex at low concentration).
- **Both-terminus init `ATCG`** (A/T at one end, G/C at the other): **ΔH° = −23.6, ΔS° = −71.81** —
  exercises the per-terminus init logic (one A/T end, one G/C end).

## Corner cases / failure modes

- **Length < 2 ⇒ `(0,0,0,0)`** — the NN model needs at least one dinucleotide step; an empty or
  single-base input has no NN contribution and no defined duplex. (An API/edge-case convention, not a
  thermodynamic value — the one documented assumption; it changes no computed quantity for length ≥ 2.)
- **Case-insensitive** — the sequence is `ToUpperInvariant`-folded before lookup (Biopython
  upper-cases too); lowercase equals uppercase.
- **Self-complementary vs non-self-complementary** — this unit fixes the **non-self-complementary
  F = 4** case (the Biopython/MELTING default); a self-complementary duplex would use F = 1.

## Assumptions & deviations

The artifact declares **zero correctness-affecting assumptions** — every constant and formula is
source-backed. The single API/edge-case convention is the **length < 2 ⇒ `(0,0,0,0)`** return above
(neither Biopython nor MELTING define an NN result for length < 2). **No source contradictions** —
Biopython, MELTING 5, and Wikipedia (citing SantaLucia 1998) agree on the DNA_NN3 table, the
two-terminus initiation, and the Tm/salt equations.

Recommended coverage (from the artifact): MUST — Biopython `CGTTCCAAAGATGTGGGCATGAGCTTAC` → 60.3 °C;
hand-derived `GCGC` → ΔH −30.0 / ΔS −84.91 / ΔG −3.67 / Tm −18.6; both-terminus init `ATCG` →
ΔH −23.6 / ΔS −71.81; length-1 and empty → `(0,0,0,0)`. SHOULD — case-insensitivity; NN-table
self-complement symmetry (AA/TT, CA/GT equal). COULD — salt monotonicity (higher [Na⁺] ⇒ higher Tm).

## References

Allawi & SantaLucia (1997) *Biochemistry* 36(34):10581 (DNA_NN3 / all97); SantaLucia (1998) *PNAS*
95(4):1460 (unified NN + salt Eq.); Biopython `Bio.SeqUtils.MeltingTemp` (`Tm_NN`, DNA_NN3);
Dumousseau et al. (2012) *BMC Bioinformatics* 13:101 (MELTING 5). Full citations in the source doc.
</content>
</invoke>
