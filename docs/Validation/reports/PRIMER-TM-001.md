# Validation Report: PRIMER-TM-001 — Primer Tm + Primer3 weighted penalty objective

- **Validated:** 2026-06-24   **Area:** MolTools
- **Canonical method(s):**
  - `PrimerDesigner.CalculateMeltingTemperature(string)` / `…WithSalt(string,double)` (Wallace / Marmur-Doty + salt; SantaLucia NN validated separately under SEQ-THERMO-001)
  - `PrimerDesigner.CalculatePrimer3Penalty(Primer3PenaltyInputs, Primer3PenaltyWeights?, Primer3Optima?)` + `DefaultPrimer3Weights` / `DefaultPrimer3Optima` (added by fix commit e55e658c, B2)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm

1. **Primer3 reference source `libprimer3.cc`** (branch main, fetched 2026-06-23/24). Read the verbatim
   `p_obj_fn` accumulation (lines 3902–3976) and `pr_set_default_global_args_1` defaults, plus the
   `gc_content` definition (line 3856).
   - `p_obj_fn` left/right-primer terms (verbatim):
     - `if (weights.temp_gt && h->temp > opt_tm) sum += weights.temp_gt * (h->temp - opt_tm);`
     - `if (weights.temp_lt && h->temp < opt_tm) sum += weights.temp_lt * (opt_tm - h->temp);`
     - `gc_content_gt`: `sum += weights.gc_content_gt * (h->gc_content - opt_gc_content);` (and `_lt` symmetric)
     - `if (weights.length_lt && h->length < opt_size) sum += weights.length_lt * (opt_size - h->length);`
     - `if (weights.length_gt && h->length > opt_size) sum += weights.length_gt * (h->length - opt_size);`
     - `compl_any`: `sum += weights.compl_any * h->self_any;`  `compl_end`: `sum += weights.compl_end * h->self_end;`
     - `num_ns`: `sum += weights.num_ns * h->num_ns;`
   - Default weights (`pr_set_default_global_args_1`): `temp_gt=temp_lt=length_gt=length_lt=1`;
     `gc_content_gt=gc_content_lt=compl_any=compl_end=num_ns=0`. Optima: `opt_size=20`, `opt_tm=60.0`.
   - `h->gc_content = 100.0 * ((double)num_gc)/num_gcat;` (line 3856) → **GC is a percentage 0–100**.
2. **Primer3 manual §19 "HOW PRIMER3 CALCULATES THE PENALTY VALUE"** (primer3.org/manual.html). Confirms the
   per-primer `PRIMER_LEFT_4_PENALTY` formula structure (one-sided `WT_TM_GT*(TM-OPT_TM)` … `WT_GC_PERCENT_LT*(OPT_GC-GC%)` …
   always-added self/num-Ns), "right primers identical to left", "lower is better", and
   **`PRIMER_OPT_GC_PERCENT` default = 50.0** (the source `#define DEFAULT_OPT_GC_PERCENT` is UNDEFINED, but the
   GC weight is 0 by default so it is inert; the user-facing documented default is 50.0).
3. **Untergasser et al. 2012 NAR 40(15):e115** — primers/pairs ranked by minimising a penalty function (defers
   term details to the manual). Peer-reviewed backing for the objective-minimisation framing.
4. **SantaLucia 1998 PNAS 95:1460 / OpenWetWare Tm-methods** — the unified NN ΔG°37 table used by
   `Calculate3PrimeStability` (validated in detail under SEQ-THERMO-001; cross-checked here, see below).
5. **Wallace rule / Marmur-Doty** — `Tm = 2(A+T)+4(G+C)` (short oligos); `Tm = 64.9 + 41(GC-16.4)/N` (Marmur & Doty 1962).

### Formula check
- **Primer3 penalty:** the TestSpec/Evidence description matches `p_obj_fn` term-for-term, including the
  weight-gated (`if (weight && …)`) short-circuit, strict sign gates (`>`/`<`), the percentage GC unit, and
  one-sided `_gt`/`_lt` weights. Default weights/optima (TM/SIZE=1, GC/SELF/NUM_NS=0; 60/20/50) are exactly sourced.
- **Tm:** Wallace and Marmur-Doty formulas and the `16.6·log10(Na/1000)` salt correction match the cited sources.

### Edge-case semantics
- Penalty: term at optimum → strict gate excludes it (contributes 0); 0-weight term inert; total ≥ 0 since every
  term is weight·(non-negative deviation). All defined & sourced.
- Tm: empty/null → 0; non-ACGT/U/N ignored; <14 valid bases → Wallace, ≥14 → Marmur-Doty. Defined behaviours.

### Independent cross-check (numbers)
Hand-computed from the published `p_obj_fn` formula (default weights unless noted):
- A: Tm=60,len=20,GC=50 → 0.0 ; B: Tm=63 → 1·3 = 3.0 ; C: Tm=57,len=18 → 3+2 = 5.0 ; D: Tm=62.5,len=22 → 2.5+2 = 4.5
- E: GcGt=0.5,GC=60 → 0.5·10 = 5.0 ; J: Tm=62,len=22,GC=55(GcGt=0.5),selfAny=2(0.25),N=1(1) → 2+2+2.5+0.5+1 = 8.0
- S2 asymmetric: TmGt=2,Tm=62 → 2·2 = 4.0 (TmLt not applied)
Tm: ATATATAT → 2·8 = 16 ; GCGCGCGC → 4·8 = 32 ; Marmur-Doty 20bp/50%GC → 64.9+41·(−6.4)/20 = 51.78 ;
salt 50mM → 16.6·log10(0.05) = −21.597.
SantaLucia ΔG°37 in code matches the canonical unified table: AA/TT=−1.00, AT=−0.88, TA=−0.58, CA/TG=−1.45,
GT/AC=−1.44, CT/AG=−1.28, GA/TC=−1.30, CG=−2.17, GC=−2.24, GG/CC=−1.84; init term. G·C=+0.98, A·T=+1.03.

### Findings / divergences
- **Note (framing, not a defect):** the orchestration prompt states the *Tm formula* uses SantaLucia NN; in this
  code `CalculateMeltingTemperature` uses Wallace/Marmur-Doty only. SantaLucia NN appears solely in
  `Calculate3PrimeStability` (3'-end ΔG), which the Tm tests do not exercise. The NN values there are correct.
- **Note (Wallace −7):** TestSpec documents the deliberate omission of the Sigma-Aldrich −7 correction; Wallace as
  published has no −7. Acceptable documented divergence (Stage A PASS-WITH-NOTES on this older Tm portion only).

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs:553-592` (`CalculatePrimer3Penalty`),
  `:516-537` (`DefaultPrimer3Weights`/`DefaultPrimer3Optima`), `:690-723` (input/weight/optima record structs).
- `:197-235` (`CalculateMeltingTemperature` / `…WithSalt`); `ThermoConstants.cs:87-122` (Wallace/Marmur-Doty/salt);
  `:427-475` (`Calculate3PrimeStability` SantaLucia NN).

### Formula realised correctly?
- `CalculatePrimer3Penalty` reproduces every `p_obj_fn` term in the same order with the same gating:
  `w.X != 0 && deviation-sign` then `sum += w.X * deviation`. GC operates on `GcPercent` (0–100) → matches
  the percentage convention. One-sided `_gt`/`_lt` preserved. Result is the bare non-negative `sum`.
- Default structs hold exactly the sourced values (1/1/1/1, 0/0/0/0/0; 60/20/50).
- Tm: Wallace branch `< 14` valid ACGT bases, else Marmur-Doty (clamped ≥0); salt = `16.6·log10(Na/1000)` added.

### Cross-verification table recomputed vs code (tests run)
`PrimerDesigner_Primer3Penalty_Tests` (17) + `PrimerDesigner_MeltingTemperature_Tests` (34) = **51 passed, 0 failed**.
Each asserts the exact hand-computed value to 1e-10 (penalty) / exact decimals (Tm). Spot values M1–M10, S1–S4,
C1–C2 and Tm M1–M19 all match my independent hand-computations above.

### Variant/delegate consistency
- `DefaultPrimer3Weights`/`DefaultPrimer3Optima` are the only defaults; M11 asserts each field. `with`-expressions
  in tests exercise non-default weights. No `*Fast` variant for the penalty. Tm shares `ThermoConstants` with
  ProbeDesigner / SequenceStatistics (consistent constants).

### Numerical robustness
- Penalty: pure additions of finite products; no division, no overflow on stated ranges; `int` length deviations
  exact. Tm: Marmur-Doty guards `length==0`; salt uses `log10` (positive Na assumed).

### Test quality audit
- Assertions check exact sourced values with `Within(1e-10)`, deterministic, no tautologies. M11 locks the sourced
  constants (drift guard). C2 guards the GC-percent-vs-fraction failure mode. S2 guards the one-sided-weight pitfall.
  Covers every Stage-A edge case.

### Findings / defects
- None. No code changed.

## 2026-06-24 update — nearest-neighbour salt-corrected Tm added (opt-in)

The earlier PASS-WITH-NOTES "full salt/Mg²⁺-corrected Tm not offered" limitation is now resolved by an
**opt-in** method, leaving only the standard Watson-Crick-NN-model residual (no mismatch / dangling-end /
secondary-structure terms). The default Wallace/Marmur-Doty Tm and the Primer3 penalty objective are unchanged.

- **New API (`PrimerDesigner`):** `CalculateMeltingTemperatureNN(string, strandConcentrationMolar=0.5e-6,
  sodiumMolar=0.05, magnesiumMolar=0, dntpMolar=0, SaltCorrectionMode=Owczarzy2004Monovalent)` and
  `CalculateNearestNeighborThermodynamics(string) → (ΔH°, ΔS°, IsSelfComplementary)?`; `SaltCorrectionMode`
  = None / SantaLuciaEntropy / Owczarzy2004Monovalent / Owczarzy2008Divalent.
- **Sources retrieved this session:** SantaLucia 1998 PNAS 95:1460 + SantaLucia & Hicks 2004 Table 1/Eq.3/Eq.5
  (Duke PDF, pages read directly); Owczarzy 2004 (Biochemistry 43:3537) monovalent + Owczarzy 2008 (47:5336)
  divalent; Biopython `DNA_NN4` / `salt_correction` (cross-check). Unified NN table cross-checked verbatim
  against Biopython `DNA_NN4`.
- **Hand-derived test values (1e-8):** ΔH°/ΔS° sums for GCGCGC (−50.4 / −134.7, self-comp) and ATGCATGC
  (−57.1 / −156.5); Tm(no salt) GCGCGC=35.0473059911, ATGCATGC=30.4338060665, CGCGAATTCGCG=61.1452300219;
  Owczarzy-2004 @50 mM GCGCGC=28.1593085080, ATGCATGC=18.1899960529; SantaLucia-Eq.5 GCGCGC=24.9976652723.
  The paper's published worked example (ΔH°=−43.5, ΔS°=−122.5, 0.2 mM → 35.8 °C) reproduces via the same Tm
  equation. Tests: `PrimerDesigner_NearestNeighborTm_Tests.cs` (17 tests, all green).
- **Evidence:** `docs/Evidence/PRIMER-TM-001-NN-Evidence.md`; **TestSpec:** `tests/TestSpecs/PRIMER-TM-001-NN.md`;
  **Algorithm doc:** `docs/algorithms/MolTools/NearestNeighbor_Salt_Corrected_Tm.md`.
- **Checklist:** root-registry Status reset `☑`→`☐` for re-validation; Quick-Reference counts adjusted.

## 2026-06-24 update — NN internal-mismatch + dangling-end Tm terms added (opt-in extension)

The earlier residual "no mismatch / dangling-end terms" is now resolved by an **opt-in extension**;
only the hairpin / secondary-structure (folding-based) Tm residual remains. The perfect-match
`CalculateMeltingTemperatureNN` and the default Wallace/Marmur-Doty Tm are unchanged.

- **New API (`PrimerDesigner`):** `CalculateMeltingTemperatureNNMismatch(string top, string bottom3to5,
  …, SaltCorrectionMode)` and `CalculateNearestNeighborThermodynamicsMismatch(string top, string bottom3to5)
  → (ΔH°, ΔS°, IsSelfComplementary)?`. The bottom strand is supplied **3'→5'** (complement direction,
  aligned base-for-base under the top), with `.` marking a single dangling base — the Biopython
  `Tm_NN(imm_table=DNA_IMM, de_table=DNA_DE)` convention.
- **Sources retrieved this session:** Allawi & SantaLucia (1997 G·T; 1998 G·A/C·T/A·C) and Peyret et al.
  (1999 A·A/C·C/G·G/T·T) internal-mismatch NN ΔH°/ΔS°; Bommarito, Peyret & SantaLucia (2000, NAR 28:1929)
  single dangling-end NN ΔH°/ΔS°; transcribed verbatim from Biopython `DNA_IMM`/`DNA_DE`. **Primary-source
  cross-checks:** the mismatch table reproduces the SantaLucia & Hicks (2004) Table 2 worked example
  `5'-GGACTGACG-3'/3'-CCTGGCTGC-5'` (ΔG°37 = −8.35 vs the paper's −8.32 kcal/mol); all 32 dangling-end
  ΔH° values reproduce SantaLucia & Hicks (2004) Table 3 ΔH° exactly (term-by-term, `pdftotext` of the
  Duke PDF).
- **Hand-derived test values (1e-8):** internal T·G mismatch `CGTGAC/GCGCTG` → ΔH°=−35.5, ΔS°=−101.5,
  Tm(no salt)=−6.4060879279 °C; 5'-dangling A `AGCGCGC/.CGCGCG` → ΔH°=−51.9, ΔS°=−136.4,
  Tm(no salt)=35.8034921829 °C. **Perfect-match equivalence:** a fully paired `GCGCGC/CGCGCG` through the
  extension equals `CalculateMeltingTemperatureNN("GCGCGC")` (35.0473059911 °C) and the ΔH°/ΔS°/self-comp
  tuple. A tandem (adjacent) mismatch / unequal length / null → null/NaN. Tests:
  `PrimerDesigner_NearestNeighborTm_Tests.cs` (now 25 tests, all green).
- **Residual (genuinely open):** hairpin / secondary-structure (folding-based) Tm only — a `[scope]` item;
  use UNAFold / ViennaRNA / MELTING 5.
- **Checklist:** root-registry Status remains `☐` (re-validation); Quick-Reference counts unchanged
  (already counted as Not-Started).

## 2026-06-25 update — DNA hairpin folding + secondary-structure (hairpin) Tm added (opt-in)

The last open residual — hairpin / secondary-structure (folding-based) Tm — is now resolved by an
**opt-in** DNA self-folder + unimolecular hairpin Tm. The perfect-match `CalculateMeltingTemperatureNN`,
the NN mismatch/dangling-end extension, and the default Wallace/Marmur-Doty Tm are all UNCHANGED.

- **New API (`PrimerDesigner`):**
  `FindMostStableHairpin(string sequence, int minStemLength = 2, double loopBonusDeltaG37 = 0)
  → HairpinResult?` (record struct: `StemStart, StemEnd, StemLength, LoopSize, DeltaH, DeltaS, DeltaG37`)
  and `CalculateHairpinMeltingTemperature(string sequence, int minStemLength = 2,
  double loopBonusDeltaG37 = 0) → double`. The folder scans every Watson-Crick closing pair, extends each
  stem maximally, closes a hairpin loop (≥ 3 nt), and keeps the minimum-ΔG°37 (MFE) hairpin.
- **Model (sources retrieved & extracted this session):** SantaLucia & Hicks (2004) Annu Rev Biophys 33:415,
  "Hairpin Loops" — the article PDF (Duke mirror) was fetched and `pdftotext`-extracted verbatim:
  ΔG°37(hairpin) = Σ stem NN stacks (Table 1, the repo's existing `NnUnifiedParams`) + ΔG°37(hairpin loop of N)
  (Table 4 verbatim: 3→3.5, 4→3.5, 5→3.3, 6→4.0, 7→4.2, 8→4.3, 9→4.5, 10→4.6, 12→5.0, 14→5.1, 16→5.3,
  18→5.5, 20→5.7, 25→6.1, 30→6.3); loop ΔH° = 0; loop ΔS° = −ΔG°37·1000/310.15 (Table 4 footnote a);
  Jacobson-Stockmayer (Eq. 7, coeff 2.44) for non-tabulated sizes. The bimolecular duplex-initiation term is
  intentionally excluded (unimolecular structure — loop init is the nucleation cost). **Unimolecular Tm
  (Eq. 11, verbatim):** `Tm = ΔH°·1000/ΔS° − 273.15` — NO R·ln(C_T/x) strand-concentration term (the
  intramolecular transition is concentration-independent; cross-confirmed by Vallone & Benight 1999).
- **Hand-derived test values (1e-9):** canonical hairpin `GGGCTTTTGCCC` (4-bp stem GGGC/GCCC + 4-nt loop) →
  ΔH° = −25.8, ΔS° = −75.48486216346927, ΔG°37 = −2.3883700000000054, Tm = 68.6403836682880 °C;
  5-nt-loop `GGGCAAAAAGCCC` → ΔG°37 = −2.5883700000000054 (loop-of-5 ΔG°37 = 3.3); poly-A → no hairpin
  (null / NaN); `GCGC` → null (loop < 3 prohibited). Tests: `PrimerDesigner_HairpinTm_Tests.cs` (16, all green).
- **Evidence:** `docs/Evidence/PRIMER-TM-001-HAIRPIN-Evidence.md`; **TestSpec:** `tests/TestSpecs/PRIMER-TM-001-HAIRPIN.md`;
  **Algorithm doc:** `docs/algorithms/MolTools/DNA_Hairpin_Folding_Tm.md`.
- **Residual (genuinely open):** self-dimer / cross-dimer (intermolecular) Tm, and the not-bundled
  triloop/tetraloop + terminal-mismatch special-loop bonus tables (the length-3/4 supplementary tables) —
  exposed as a caller-supplied `loopBonusDeltaG37` increment; use UNAFold / ViennaRNA / MELTING 5 for those.
- **Checklist:** root-registry Status remains `☐` (re-validation); Quick-Reference counts unchanged
  (already counted as Not-Started).

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** (Tm portion: documented Wallace −7 omission + Marmur-Doty simplification; prompt's
  "Tm uses SantaLucia NN" is a framing inaccuracy — NN is in 3'-stability only). **Stage B: PASS.**
- **State: CLEAN.** The new Primer3 penalty objective (commit e55e658c) faithfully reproduces Primer3's `p_obj_fn`
  left/right-primer branch with exact sourced defaults; hand-computed penalties match the C# output to 1e-10. The
  legacy Tm calculation is unchanged and correct per its documented (simplified) model. The new opt-in NN
  salt-corrected Tm (2026-06-24) is added above; only the Watson-Crick-NN-model residual remains.
- No defects logged.
