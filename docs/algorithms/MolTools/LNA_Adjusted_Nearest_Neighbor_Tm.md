# LNA-Adjusted Nearest-Neighbour Melting Temperature

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PROBE-LNATM-001 (formerly filed under PROBE-DESIGN-001) |
| Related Projects | Seqeron.Genomics.MolTools, Seqeron.Genomics.Infrastructure |
| Implementation Status | Production |
| Last Reviewed | 2026-10-01 (review campaign 2026-09, B07 F26–F28) |

## 1. Overview

Nearest-neighbour (NN) thermodynamics and melting temperature of a DNA oligonucleotide carrying
**internal LNA (locked nucleic acid) monomers**, hybridised to an unmodified DNA strand (its perfect
complement or a supplied target with mismatches). Two published LNA·DNA NN models are implemented,
both on the SantaLucia (1998) unified DNA parameters, exactly as the MELTING 5 reference
implementation realises them [2,3]:

- **Owczarzy et al. (2011)** [6] — *default* (MELTING's default): complete NN parameters for steps
  with one LNA (32), two consecutive LNAs (16) and an LNA·DNA mismatch inside an LNA triplet.
- **McTigue, Peterson & Kahn (2004)** [1]: ΔΔH°/ΔΔS° increments for an isolated internal LNA added
  to the DNA stacks (32 increments).

It also provides a qualitative 3'-MGB (minor-groove binder) probe design-rule check from Kutyavin
et al. (2000) [4].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

LNA monomers (2'-O,4'-C-methylene ribose) lock the sugar in the C3'-endo conformation and stabilise
the duplex by a sequence-dependent amount (average ΔTm per single LNA: +4.4 °C C, +3.2 °C T,
+2.8 °C G, +2.1 °C A [1]). LNA-modified primers/probes (and LNA probes for allele discrimination)
need an LNA-aware Tm.

### 2.2 Core Model

Duplex ΔH°, ΔS° (1 M Na⁺) = initiation + Σ NN steps (+ symmetry, LNA-free self-complementary only):

| Term | Source / value |
|------|----------------|
| Initiation | SantaLucia (1998) unified = Allawi & SantaLucia (1997): per terminal A·T +2.3 kcal/mol / +4.1 eu, per terminal G·C +0.1 / −2.8 [5] (Biopython `DNA_NN3`, MELTING `all97`) |
| Step without LNA, Watson–Crick | unified stack (Allawi & SantaLucia 1997) |
| Step without LNA, internal mismatch | Allawi / SantaLucia / Peyret internal-mismatch NN (Biopython `DNA_IMM1`, MELTING `allsanpey`) |
| Step touching an **isolated** LNA, McTigue model | DNA stack + McTigue ΔΔH°/ΔΔS° (key = step + locked base, e.g. `TTL/AA` = TT with the 3' T locked) [1] |
| Step with one LNA, Owczarzy model (and LNAs of a run in the McTigue model) | complete Owczarzy single-LNA parameter, no DNA base added [6] |
| Step with two LNAs, Watson–Crick | Owczarzy consecutive-LNA parameter (e.g. `CLCL/GG` −15.399 / −36.375) |
| Step with two LNAs, mismatch opposite one | Owczarzy LNA-mismatch parameter (e.g. `GLTL/TA`) |

Tm = 1000·ΔH° / (ΔS° + R·ln(C_T/x)) − 273.15, x = 4 (x = 1 for an LNA-free self-complementary
duplex), followed by the selected salt correction (default Owczarzy et al. 2004 Eq. 22, MELTING's
default for Na⁺; Owczarzy 2008 for Mg²⁺/dNTP), computed by the canonical
`ThermoConstants.CalculateNearestNeighborTmFromThermodynamics`.

Note: Owczarzy's single-LNA ΔH° equal the unified stack + McTigue increment exactly (e.g.
`TTL/AA` −5.574 = −7.9 + 2.326); the ΔS° were refitted (−14.149 vs −14.1), so the two models agree
closely for isolated LNAs (`CCATT(L)GCTACC`, 1e-4 M: 63.614 vs 63.483 °C).

### 2.3 Modeling Assumptions

1. LNAs on one strand only (LNA·DNA duplex); the opposite strand is unmodified DNA.
2. NN additivity, two-state melting, 1 M Na⁺ reference state.
3. An LNA strand with its DNA complement is a heteroduplex — never treated as self-complementary
   (MELTING rejects `-self` with LNAs).

### 2.4 Properties and Invariants

| ID | Invariant |
|----|-----------|
| INV-01 | No LNA ⇒ equals the unified DNA NN model (Biopython `Tm_NN(nn_table=DNA_NN3)` with the same R). |
| INV-02 | ΔH°/ΔS° equal MELTING 5.2.0 bit-exactly for every computable duplex; Tm equal with R = 1.99 whenever MELTING's f(GC) convention coincides (see 5.4). |
| INV-03 | Not computable (null/NaN) exactly where MELTING reports missing parameters, and additionally for a terminal LNA (MELTING 5.2.0 evaluates it by extrapolation — see 5.4). |
| INV-04 | Position order/duplicates irrelevant; case-insensitive. |

## 3. Contract

### 3.1 Inputs and Parameters

| Parameter | Meaning |
|-----------|---------|
| `sequence` | LNA-modified oligo as DNA letters 5'→3' (≥ 2 ACGT). |
| `lnaPositions` | zero-based LNA positions (internal only). |
| `model` | `LnaNearestNeighborModel.Owczarzy2011` (default overloads) or `McTigue2004`. |
| `target` | opposite DNA strand 3'→5', same length (`null` = perfect complement). |
| `strandConcentrationMolar`, `sodiumMolar`, `magnesiumMolar`, `dntpMolar`, `saltMode` | as `CalculateMeltingTemperatureNN`. |
| `gasConstant` | R (default 1.9872; MELTING uses 1.99). |

### 3.2 Output / Return Value

`CalculateNearestNeighborThermodynamicsLna` → `(ΔH° kcal/mol, ΔS° cal/(K·mol), IsSelfComplementary)?`;
`CalculateMeltingTemperatureNNLna` → Tm (°C) or `NaN`.

### 3.3 Preconditions and Validation

- null `lnaPositions` → `ArgumentNullException`; unknown model / C_T ≤ 0 / [Na⁺] ≤ 0 / [Mg²⁺] < 0 /
  [dNTP] < 0 / R ≤ 0 → `ArgumentOutOfRangeException`.
- Not computable (null / NaN): empty, < 2 nt, non-ACGT strand; target of another length; terminal
  or out-of-range LNA; terminal mismatch; a mismatch opposite an isolated LNA or at the end of an
  LNA run; any LNA step without a published parameter.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate; build the bottom strand; mark LNA positions; reject terminal LNA / terminal mismatch.
2. Add the unified initiation terms of the two terminal pairs.
3. For each step: no LNA → WC/mismatch DNA value; isolated LNA under McTigue → stack + increment;
   otherwise Owczarzy single / consecutive / mismatch table.
4. Tm by the canonical bimolecular equation + salt correction.

### 4.2 Decision Rules / Reference Tables

Tables transcribed from MELTING 5.2.0 `McTigue2004lockedmn.xml`, `Owczarzy2011lockedmn.xml`,
`Owczarzy2011lockedTandemmn.xml`, `Owczarzy2011lockedmmn.xml` (cal/mol ÷ 1000).

### 4.3 Complexity

O(n) time, O(n) memory.

## 5. Implementation Notes

### 5.1 Location and Entry Points

`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs`:
`CalculateNearestNeighborThermodynamicsLna` (2 overloads), `CalculateMeltingTemperatureNNLna`
(2 overloads), `LnaNearestNeighborModel`. DNA step lookup:
`ThermoConstants.TryGetNearestNeighborDuplexStep`. MGB rules: `ProbeDesigner.EvaluateMgbProbeDesign`.
No MCP tool exposes these methods.

### 5.2 Current Behavior

Default overloads use Owczarzy (2011) with the perfect complement.

### 5.3 Conformance to Theory / Spec

Differential check vs melting5.jar 5.2.0 (`Main.getMeltingResults`, full precision; 4 200 random
duplexes: 6–30 nt, 1–8 LNAs incl. runs, LNA-triplet and DNA mismatches, C_T 0.25–100 µM,
Na⁺ 0.05–1 M, Mg²⁺ 0/3 mM, both models): ΔH°/ΔS° identical in every computable case, identical
not-computable set for internal LNAs, Tm identical (≤ 1e-8 °C, R = 1.99) except the f(GC) convention cases below.
Terminal LNAs are the one deliberate difference (re-run 2026-10-08, 600 random duplexes with ~15 % terminal
LNAs: all 67 differences are terminal-LNA duplexes that MELTING evaluates and this library returns NaN; 106
jointly not computable; see 5.4).
rmelting `test-method.locked.R`: 63.61426 (mct04), 63.48299 (owc11), 12.94323 (GALCLC) reproduced.

Agreement with measured Tm (MELTING-shipped data sets): McTigue (2004) 100 single-LNA duplexes
(5 µM, 1 M Na⁺) mean |error| 1.37 °C (Owczarzy) / 1.38 °C (McTigue); Owczarzy (2011) LNA-triplet
duplexes (2 µM) 1.01 °C perfect match, 2.87 °C central mismatch.

### 5.4 Deviations and Assumptions

- **f(GC) of the Owczarzy 2004/2008 salt correction** = G+C fraction of the DNA letters (Owczarzy
  2004 definition: fraction of G·C pairs; Biopython convention for mismatches). MELTING computes it
  with `isBasePairEqualTo("G","C")`, which does not recognise an LNA G·C pair ("GL") nor a G/C
  opposite a mismatch — a MELTING artefact; Tm differs from MELTING only when Na⁺ ≠ 1 M and such
  pairs are present.
- **R**: 1.9872 by default (SantaLucia & Hicks 2004); MELTING hard-codes 1.99 (`gasConstant: 1.99`
  for parity).
- **Data anomalies kept as MELTING** (paper tables behind the ACS paywall): mismatch table has no
  `GLAL/CA` but a double-mismatch key `CLAL/CA`; `GLTL/TA` listed twice (later −14.213/−40.041
  kept, as MELTING's HashMap); consecutive `TLCL/AG` ΔS −21.735 (an independent transcription,
  iCarrin/Bio_dpt `lna_tm.py`, has −21.535).
- **McTigue model with consecutive LNAs**: McTigue (2004) measured single internal LNAs only;
  runs use the Owczarzy (2011) tables, as MELTING does.
- **Terminal LNA (position 0 or n − 1) → not computable; MELTING 5.2.0 does compute it.** Neither
  paper parameterises an LNA at a duplex end (McTigue 2004: single internal LNAs; MELTING's own
  `isApplicable` of `McTigue04LockedAcid` / `Owczarzy11LockedAcid` / `LockedAcidNNMethod` is written to
  warn "The thermodynamics parameters for locked nucleic acids … are not established for terminal
  locked nucleic acids." and return false). That guard never fires in 5.2.0: it compares the terminal
  base pair with the literal pattern `"L"`/`"-"`, which an LNA base (`"CL"`, `"AL"`, …) never equals,
  so neither the warning nor the rejection is emitted and MELTING applies the internal-LNA doublet
  parameter at the end (e.g. 1e-4 M, 1 M Na⁺: `CLCATTGCTACC` owc11 66.65650883512683 °C (step
  `CLC/G G` −5.904/−11.904), mct04 66.66234775338887; `CCATTGCTACCL` owc11 66.17003475528679 °C).
  That value is an extrapolation of internal-LNA parameters outside their validated range, contrary
  to MELTING's documented intent, so it is deliberately not reproduced (B07 F66).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

See 3.3; mismatch opposite an isolated LNA has no published parameter (also missing in MELTING).

### 6.2 Limitations

- LNA on both strands, terminal LNAs, dangling ends and terminal mismatches are not parameterised
  by either paper (not computable). MELTING 5.2.0 also rejects terminal mismatches ("No method for
  terminal mismatches") but evaluates terminal LNAs through a non-firing guard (see 5.4).
- MGB: only the qualitative Kutyavin (2000) rules (3'-attachment, 12–20mer window). The quantitative
  MGB ΔTm is not computed — Kutyavin (2000) reports measured Tm only; the vendor model (Primer Express
  `MGB_dds` entropy term, Epoch/ELITech patent US 7,715,989) has no obtainable parameter values
  (LimitationPolicy `PROBE-DESIGN-001`).

## 7. Examples and Related Material

### 7.1 Worked Example

`CCATT(L)GCTACC` (LNA at index 4), C_T = 1e-4 M, 1 M Na⁺:

- McTigue: ΔH° = −81.1 (unified) + 2.326 (`TTL/AA`) − 1.540 (`TLG/AC`) = −80.314 kcal/mol;
  ΔS° = −222.5 + 8.1 − 3.0 = −217.4 eu; Tm(R = 1.99) = 63.614259 °C (MELTING `mct04` 63.61426).
- Owczarzy: `TTL/AA` (−5.574/−14.149), `TLG/AC` (−10.040/−25.744) → −80.314 / −217.493,
  Tm(R = 1.99) = 63.482987 °C (MELTING default 63.48299).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ProbeDesigner_LnaTm_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/ProbeDesigner_LnaTm_Tests.cs)
- TestSpec: [PROBE-LNATM-001.md](../../../tests/TestSpecs/PROBE-LNATM-001.md)
- Evidence: [PROBE-DESIGN-001-LNA-Evidence.md](../../../docs/Evidence/PROBE-DESIGN-001-LNA-Evidence.md)
- Review: `docs/Validation/review-2026-09/B07.md` F26–F28
- Related algorithms: [NearestNeighbor_Salt_Corrected_Tm](NearestNeighbor_Salt_Corrected_Tm.md), [Hybridization_Probe_Design](Hybridization_Probe_Design.md)

## 8. References

1. McTigue PM, Peterson RJ, Kahn JD. 2004. Sequence-dependent thermodynamic parameters for locked nucleic acid (LNA)-DNA duplex formation. Biochemistry 43(18):5388–5405. https://doi.org/10.1021/bi035976d
2. Dumousseau M, Rodriguez N, Juty N, Le Novère N. 2012. MELTING, a flexible platform to predict the melting temperatures of nucleic acids. BMC Bioinformatics 13:101.
3. MELTING 5.2.0 (`melting5.jar`, Java sources + data files) as shipped in Bioconductor rmelting (github.com/aravind-j/rmelting `inst/java`, `inst/extdata/Data`; tests `tests/testthat/test-method.locked.R`).
4. Kutyavin IV, Afonina IA, Mills A, et al. 2000. 3'-Minor groove binder-DNA probes increase sequence specificity at PCR extension temperatures. Nucleic Acids Res 28(2):655–661. https://doi.org/10.1093/nar/28.2.655
5. SantaLucia J. 1998. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. PNAS 95(4):1460–1465; Allawi HT, SantaLucia J. 1997. Biochemistry 36:10581.
6. Owczarzy R, You Y, Groth CL, Tataurov AV. 2011. Stability and mismatch discrimination of locked nucleic acid–DNA duplexes. Biochemistry 50(43):9352–9367. https://doi.org/10.1021/bi200904e
7. Owczarzy R, et al. 2004. Biochemistry 43:3537 (Na⁺ correction); 2008. Biochemistry 47:5336 (Mg²⁺).
