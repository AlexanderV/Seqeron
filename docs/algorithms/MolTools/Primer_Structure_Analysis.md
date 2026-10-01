# Primer Structure Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Molecular Tools |
| Test Unit ID | PRIMER-STRUCT-001 |
| Related Projects | N/A |
| Implementation Status | Complete (Primer3 semantics) |
| Last Reviewed | 2026-10-01 |

## 1. Overview

Primer structure analysis screens PCR primers for secondary structures and for features that
reduce amplification efficiency: hairpins, self- and cross-dimers (in particular 3′-end dimers
that the polymerase can extend), an over-stable 3′ end, and homopolymer / dinucleotide runs.
The implementation reproduces Primer3 (Untergasser et al. 2012; `libprimer3.cc`, `oligotm.c`,
`dpal.c`, `thal.c`) for every quantity that Primer3 defines:

| Quantity | Primer3 output tag | Method |
|----------|--------------------|--------|
| Hairpin / self-dimer / 3′ self-dimer Tm (thermodynamic, Primer3 default) | `PRIMER_*_HAIRPIN_TH`, `_SELF_ANY_TH`, `_SELF_END_TH` | `CalculatePrimer3OligoStructure` |
| Pair hetero-dimer / 3′ hetero-dimer Tm (thermodynamic, Primer3 default) | `PRIMER_PAIR_COMPL_ANY_TH`, `_COMPL_END_TH` | `CalculatePrimer3PairComplementarity` |
| Pair 3′ complementarity, alignment mode | `PRIMER_PAIR_COMPL_END` (`PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0`) | `CalculatePrimerDimerEndComplementarity`, `HasPrimerDimer` |
| Self 3′ complementarity, alignment mode | `PRIMER_*_SELF_END` | `CalculatePrimerSelfEndComplementarity` |
| 3′-end stability | `PRIMER_*_END_STABILITY` (opposite sign) | `Calculate3PrimeStability` |
| Longest mononucleotide run | `PRIMER_MAX_POLY_X` check | `FindLongestHomopolymer` |
| ntthal ANY / END1 / END2 dimer, hairpin with Mg²⁺/dNTP | `ntthal -a ANY/END1/END2/HAIRPIN` | `CalculateDimerThermodynamicsNtthal(…, mode, …)`, `CalculateHairpinThermodynamicsNtthal(…, dv, dntp)` |

Two sequence-only screens without a Primer3 counterpart remain: `HasHairpinPotential` (exact
Watson–Crick stem of ≥ `minStemLength` closing a loop of ≥ `minLoopLength`) and
`FindLongestDinucleotideRepeat`.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Hairpins form when a primer contains an inverted repeat (stem) closing a loop; loops shorter
than 3 nt are sterically excluded (SantaLucia & Hicks 2004). Primer-dimers arise when two primers
(or two copies of one) anneal; a dimer whose duplex contains a primer's 3′-terminal base can be
extended by the polymerase (Primer3 "END" alignments). Primer3's default since 2.x
(`PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1`) scores all of these as the melting temperature of the
most stable structure found by the `ntthal` nearest-neighbour dynamic programme (SantaLucia &
Hicks 2004 parameters) and rejects a primer/pair when a Tm exceeds 47 °C. The older alignment
mode scores 3′ complementarity with a `dpal` end-anchored alignment. The 3′-end stability is the
nearest-neighbour ΔG°37 of the last five bases (SantaLucia 1998).

### 2.2 Core Model

**Thermodynamic screen (Primer3 default).** For a primer `p` (5′→3′) Primer3's
`oligo_compl_thermod` / `oligo_hairpin` compute
`self_any = ntthal ANY(p, p)`, `self_end = ntthal END1(p, p)`, `hairpin = ntthal HAIRPIN(p)`; for a
pair (left `L`, right `R`, both 5′→3′) `characterize_pair` computes
`compl_any = ntthal ANY(L, R)` and
`compl_end = max(END1(L, R), END2(L, R), END1(rc R, rc L), END2(rc R, rc L))`.
Each value is the Tm (°C) of the optimal structure, reported as 0 when no structure forms or the
Tm is negative (`align_thermod`). Conditions are Primer3's primer conditions: 50 mM monovalent,
1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM oligo; Mg²⁺/dNTP enter ntthal only through
`saltCorrectS = 0.368·ln((mv + 120·√max(0, dv − dntp))/1000)`. END1 forces the 3′-terminal base
of the first strand into the terminal pair; END2 = END1 with the strands swapped. Default limits
`PRIMER_MAX_SELF_ANY_TH = PRIMER_MAX_SELF_END_TH = PRIMER_MAX_HAIRPIN_TH =
PRIMER_PAIR_MAX_COMPL_ANY_TH = PRIMER_PAIR_MAX_COMPL_END_TH = 47 °C`
(`PrimerDesigner.Primer3MaxStructureTm`); a value strictly greater fails.

**Alignment-mode 3′ complementarity.** `compl_end = max(align(L, rc R), align(R, rc L))`, where
`align` is `dpal` with flag `DPAL_GLOBAL_END` (the alignment must end at the last base of the
first sequence), match +1.00, mismatch −1.00, N −0.25, single-base gaps −2.00 (max gap 1), floored
at 0. Two primers whose 3′-terminal k bases are reverse complements score k; Primer3's default
limit `PRIMER_PAIR_MAX_COMPL_END = 3.00`, so `HasPrimerDimer(p1, p2, minComplementarity = 4)`
(score ≥ 4) is exactly the default Primer3 rejection for ACGT primers (integral scores).
`self_end = align(p, rc p)`.

**3′-end stability.** Primer3 `end_oligodg(seq, 5, santalucia)`: over the last five bases (the
whole primer if shorter) −ΔG = Σ SantaLucia (1998) NN −ΔG°37 − 1.96 − 0.05·(terminal A/T count)
− 0.43·(self-complementary); the library returns ΔG (negative = stable). For a 5-mer this equals
SantaLucia's "initiation with terminal G·C +0.98 / A·T +1.03" form. N uses Primer3's N row/column.
GCGCG = −6.86 (most stable), TATAT = −0.86.

**Poly-X.** Longest run of identical bases; N is a worst-case wildcard exactly as
`_pr_violates_poly_x` (forward scan assigning N to the preceding base, reverse scan to the next
base, longer run reported): ANA 3, GNGNG 5, ANGNG 4.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `HasHairpinPotential(...)` is `false` when the sequence is shorter than `2·minStemLength + minLoopLength` | Length guard |
| INV-02 | `Calculate3PrimeStability(...)` depends only on the last five bases; empty input → 0 | `end_oligodg` window |
| INV-03 | `FindLongestHomopolymer(...)` is 0 for empty input and ≥ 1 otherwise | Run scan |
| INV-04 | `FindLongestDinucleotideRepeat(...)` is 0 for inputs shorter than 4 nt | Short-circuit |
| INV-05 | `CalculatePrimerDimerEndComplementarity(a, b) = CalculatePrimerDimerEndComplementarity(b, a) ≥ 0` | Max over both orientations, floor 0 |
| INV-06 | Primer3 structure Tm values are ≥ 0 | `align_thermod` floors negative Tm at 0 |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `[HasHairpinPotential] minStemLength / minLoopLength` | `int` | 4 / 3 | Stem length and minimum loop of the sequence-only screen |
| `[HasPrimerDimer] minComplementarity` | `int` | 4 | Flag when the Primer3 alignment-mode compl_end score ≥ this value |
| `[CalculatePrimer3OligoStructure / PairComplementarity] monovalentMillimolar, divalentMillimolar, dntpMillimolar, dnaConcentrationNanomolar` | `double` | 50, 1.5, 0.6, 50 | Primer3 primer conditions |
| `[CalculateDimerThermodynamicsNtthal] mode` | `NtthalAlignmentMode` | — | `Any`, `End1`, `End2` |
| `[PrimerParameters] StructureScreen` | `PrimerStructureScreen` | `Primer3Thermodynamic` | Screen used by `EvaluatePrimer` / `DesignPrimers` (`Heuristic` = `HasHairpinPotential` + `HasPrimerDimer`) |
| `[PrimerParameters] MaxStructureTm` | `double` | 47 | Primer3 `*_TH` limit (0 → 47) |

### 3.2 Output / Return Value

| Method | Output |
|--------|--------|
| `CalculatePrimer3OligoStructure` | `Primer3OligoStructure(SelfAnyTh, SelfEndTh, HairpinTh)` in °C, or `null` for null/empty/non-ACGT |
| `CalculatePrimer3PairComplementarity` | `Primer3PairComplementarity(ComplAnyTh, ComplEndTh)` in °C, or `null` |
| `CalculatePrimerDimerEndComplementarity` / `CalculatePrimerSelfEndComplementarity` | Primer3 score (≥ 0; 0 for null/empty) |
| `HasPrimerDimer` | `bool` |
| `Calculate3PrimeStability` | ΔG°37 kcal/mol; 0 for null/empty; `NaN` when the 3′ window has a character other than ACGTN |
| `FindLongestHomopolymer` / `FindLongestDinucleotideRepeat` | `int` |
| `EvaluatePrimer` | `PrimerCandidate` with `SelfAnyTh`, `SelfEndTh`, `HairpinTh` (thermodynamic screen) and `HasHairpin` = `HairpinTh > MaxStructureTm` |

### 3.3 Preconditions and Validation

All methods are case-insensitive. The thermodynamic methods accept ACGT only. The alignment
methods score any non-ACGT character as N (Primer3 `p3_reverse_complement` turns it into N).

## 4. Algorithm

1. `EvaluatePrimer` computes the per-primer constraints; with the default screen it adds an issue
   for each of hairpin / self-dimer / 3′ self-dimer Tm above the limit.
2. `DesignPrimers` follows Primer3's `characterize_pair`: candidates are ranked by penalty; the
   per-primer structure screen runs lazily (once per primer, cached) only for primers that reach
   a pair passing the ΔTm check, followed by the pair `compl_any_th` / `compl_end_th` check.
3. `dpal` GLOBAL_END is a line-for-line port of `_dpal_long_nopath_maxgap1_global_end` (the routine
   Primer3 runs); for |X| ≤ 3 or |Y| = 1, where that C routine reads past the sequence end, the
   `_dpal_generic` GLOBAL_END recurrence is used.

### 4.3 Complexity

| Operation | Time |
|-----------|------|
| ntthal dimer / hairpin | O(n·m·L²), L = max loop 30 |
| dpal GLOBAL_END | O(n·m) |
| 3′ stability, poly-X | O(1) / O(n) |
| `HasHairpinPotential` | O(n²) below 100 nt, suffix tree at ≥ 100 nt |

## 5. Implementation Notes

**Implementation location:** [PrimerDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs),
[NtthalDimer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/NtthalDimer.cs) (END1/END2, divalent salt),
[NtthalHairpin.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/NtthalHairpin.cs) (divalent salt).

### 5.3 Conformance to Theory / Spec

- Cross-checked against Primer3 (primer3-py 2.3.1, and Primer3 C sources compiled locally):
  dpal compl_end 3000/3000 random pairs identical to `dpal.c`; 600/600 `END_STABILITY`, 600/600
  `SELF_END`, 300/300 `PAIR_COMPL_END` identical to `design_primers` (alignment mode);
  `end_oligodg` 2000/2000 identical (incl. N and primers < 5 nt); poly-X with N identical to
  `check_primers`.
- Thermodynamic values are as exact as the ntthal engines: the structure-screen formulas match
  `design_primers` whenever the engines match `calc_homodimer` / `calc_end_stability` /
  `calc_hairpin` (all END1/END2/dv code paths verified to 1e-9 on the engine-exact cases).
  Residual engine discrepancies (≈ 2 % of random dimers, ≈ 5 % of random hairpins, see B07
  report) are owned by PRIMER-DIMER-001 / PRIMER-HAIRPIN-001.
- `DesignPrimers` vs primer3-py `design_primers` (thermodynamic default, this library's per-primer
  limits): 574/600 random templates identical; every one of the 26 differences traced to an ntthal
  engine value differing from primer3-py.

## 6. Edge Cases and Limitations

| Case | Behaviour |
|------|-----------|
| Identical poly-A primers | Not a dimer (compl_end 0; ntthal no structure) |
| 3′ ends …GGCC / …GGCC | compl_end 4 → dimer (offset overlap) |
| Primer shorter than 5 nt | 3′ stability of the whole primer (Primer3 `end_oligodg`) |
| No ntthal structure / negative Tm | Structure Tm reported as 0 |
| Non-ACGT primer | No thermodynamic values (`null`); `EvaluatePrimer` already rejects it (Tm) |

`FindLongestDinucleotideRepeat` counts any repeated 2-mer (so AAAA counts as 2 "AA" units) and has
no Primer3 equivalent; `HasHairpinPotential` is a sequence-only screen kept for the `Heuristic`
mode and the `hairpin_potential` MCP tool.

## 8. References

1. Untergasser A et al. (2012) Primer3 — new capabilities and interfaces. NAR 40:e115.
2. Primer3 source (`libprimer3.cc` `characterize_pair`, `oligo_compl_thermod`, `oligo_hairpin`,
   `align`, `align_thermod`, `_pr_violates_poly_x`; `oligotm.c` `oligodg`, `end_oligodg`, `symmetry`;
   `dpal.c`; `thal.c`) — https://github.com/primer3-org/primer3 and the primer3-py vendored copy.
3. SantaLucia J (1998) PNAS 95:1460-65 (Table 1). SantaLucia J, Hicks D (2004) Annu Rev Biophys 33:415-40.
4. Rozen S, Skaletsky H (2000) Primer3 on the WWW for general users and for biologist programmers. Methods Mol Biol 132:365-86.
