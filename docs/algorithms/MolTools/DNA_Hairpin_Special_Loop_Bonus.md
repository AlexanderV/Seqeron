# DNA Hairpin Special Tri/Tetraloop Bonus (full ntthal hairpin DP)

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PRIMER-TM-001 |
| Related Projects | Seqeron.Genomics.MolTools |
| Implementation Status | Production |
| Last Reviewed | 2026-06-25 |

## 1. Overview

Computes the full Primer3 `ntthal` intramolecular-hairpin thermodynamics (ΔH°, ΔS°, ΔG°37, Tm) of a
DNA oligo, **automatically applying the bundled sequence-specific special triloop (3-nt) and
tetraloop (4-nt) stability bonus tables**. These tables are the extra stability that recognised loops
(e.g. GNRA tetraloops) gain over a generic loop of the same length; before this work the bonus had to
be supplied by hand. The method reproduces primer3-py `calc_hairpin` to machine precision and is
opt-in: the legacy SantaLucia & Hicks (2004) Table 4 hairpin model and the duplex/dimer Tm methods are
unchanged [1][3].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A DNA oligo can fold back on itself into a hairpin: a Watson-Crick stem closing an unpaired loop. The
loop free energy depends on its length and, for short loops (3 and 4 nt), on the *specific* loop
sequence including the closing base pair — certain loops are unusually stable and have measured bonus
energies [1].

### 2.2 Core Model

The hairpin free energy is the sum of (a) the stem nearest-neighbour stacks, (b) the size-keyed
hairpin-loop initiation, (c) for loops ≥ 4 nt a terminal-mismatch (`tstack2`) increment / for loops
of 3 nt a closing-A·T penalty, and (d) — the special-loop bonus — an additive ΔH°/ΔS° from the
triloop / tetraloop tables, keyed on the full loop string including the closing pair [3]:

- **Triloop bonus** (loop size 3): key = closing-5′ base + 3 loop nt + closing-3′ base (5 chars); ΔH in
  cal/mol (e.g. CGAAG = −2000), ΔS = 0 for all triloops [primer3 `triloop.dh`/`.ds`].
- **Tetraloop bonus** (loop size 4): key = closing-5′ base + 4 loop nt + closing-3′ base (6 chars);
  ΔH/ΔS in cal/mol and cal/(K·mol) (e.g. CGAAAG = −1100/0) [primer3 `tetraloop.dh`/`.ds`].

The unimolecular melting temperature is `Tm = ΔH°/(ΔS° + (N/2 − 1)·saltCorrection) − 273.15`, with
`saltCorrection = 0.368·ln([Na⁺])` and no strand-concentration term (the transition is intramolecular)
[3]. The bonus tables are the verbatim libprimer3 parameter files that `ntthal` loads; their values
trace to SantaLucia & Hicks (2004) "special hairpin loops" [1].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | A recognised tri/tetraloop hairpin's ΔH/ΔS/ΔG/Tm equals primer3 `calc_hairpin`. | The DP + bonus tables are a verbatim port of `thal.c` + the libprimer3 config files [3]. |
| INV-02 | A non-special-loop hairpin is identical with or without the bundled tables. | The bonus is added only on a bsearch match; non-special loops never match [3]. |
| INV-03 | No hairpin (homopolymer / too short) → `null`. | ntthal `no_structure` (`mh`/`ms` not finite) [3]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | DNA oligo (5′→3') | ACGT only; case-insensitive |
| sodiumMolar | double | 0.05 | Monovalent cation concentration (mol/L) | > 0; primer3 `mv=50` mM |

### 3.2 Output / Return Value

`HairpinThermodynamics?` — ΔH° (kcal/mol), ΔS° (cal/(K·mol), incl. salt term), ΔG°37 (kcal/mol),
Tm (°C), BasePairs (ntthal N/2). `null` when no hairpin forms or input is invalid.

### 3.3 Preconditions and Validation

null / empty / any non-ACGT character → `null`. The sequence is uppercased internally. A homopolymer
or an oligo too short to close a ≥ 3-nt loop returns `null` (matching primer3 `structure_found=False`).

## 4. Algorithm

### 4.1 High-Level Steps

1. 1-index the sequence (`numSeq1 == numSeq2`, the monomer convention).
2. `initMatrix2` / `fillMatrix2`: DP over closing pairs (i,j), stacking inward (`maxTM2`), internal
   loops/bulges (`CBI`/`calc_bulge_internal2`), and closing the loop (`calc_hairpin`).
3. `calc_hairpin` adds the size loop init + tstack2/AT-penalty + **the special triloop/tetraloop bonus**.
4. `calc_terminal_bp` / `END5_1..4`: best 5′ exterior assembly; `tracebacku` counts paired bases N.
5. `calcHairpin`: `Tm = ΔH°/(ΔS° + (N/2 − 1)·saltCorrection) − 273.15`.

### 4.2 Decision Rules, Scoring, Reference Tables

The special-loop bonus tables (16 triloops, 76 tetraloops) and the hairpin loop-length ΔS column are
embedded verbatim from the libprimer3 `primer3_config/{triloop,tetraloop,loops}.{dh,ds}` files with a
provenance header; the stem / terminal-mismatch / dangle / interior / bulge tables are reused from
`NtthalDimer` (same source set) [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| hairpin DP | O(n³) | O(n²) | n = oligo length; the standard ntthal monomer DP |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [NtthalHairpin.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/NtthalHairpin.cs),
[PrimerDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs)

- `PrimerDesigner.CalculateHairpinThermodynamicsNtthal(string, double)`: public entry; validates input,
  runs the ntthal hairpin DP, converts ntthal cal/mol → kcal/mol.
- `NtthalHairpin.Run(string, double)`: the full monomer DP + bundled special-loop bonus lookups.

### 5.2 Current Behavior

The bonus is applied only inside the new ntthal path. The legacy
`PrimerDesigner.FindMostStableHairpin` / `CalculateHairpinMeltingTemperature` (SantaLucia & Hicks 2004
Table 4 model) are unchanged and still accept a caller-supplied `loopBonusDeltaG37`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- The full ntthal hairpin DP and the triloop/tetraloop bonus tables (verbatim libprimer3 files),
  keyed on the full loop string including the closing pair, added to loop ΔH°/ΔS° [3].
- Unimolecular Tm with the `(N/2 − 1)·saltCorrection` entropy term and no concentration term [3].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Divalent (Mg²⁺) / dNTP salt terms: the public method exposes only monovalent `[Na⁺]` (dv=dntp=0),
  matching the captured primer3 conditions; **users should rely on:** primer3 directly for divalent salt.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| homopolymer (poly-A) | null | no Watson-Crick stem (ntthal no_structure) |
| too short (GCGC) | null | cannot close a ≥3-nt loop |
| null / empty / non-ACGT | null | input guards |
| non-special 4-nt loop | identical to no-bonus | bsearch no-match |

### 6.2 Limitations

Monovalent salt only (dv=dntp=0). Special bonus tables exist only for 3-nt and 4-nt loops (the
biology — there are no measured 5-nt+ special-loop tables).

## 7. Examples and Related Material

### 7.1 Worked Example

```csharp
// GNRA tetraloop CGAAAG → primer3 calc_hairpin parity
var h = PrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGGCGAAAGCCCC", 0.05);
// h.DeltaH*1000 == -40900 ; h.DeltaS == -114.1872884299936 ; h.TmCelsius == 85.03347700825856
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PrimerDesigner_HairpinSpecialLoop_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/PrimerDesigner_HairpinSpecialLoop_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`
- Evidence: [PRIMER-TM-001-SPECIAL-LOOP-Evidence.md](../../../docs/Evidence/PRIMER-TM-001-SPECIAL-LOOP-Evidence.md)
- Related algorithms: [DNA_Hairpin_Folding_Tm](DNA_Hairpin_Folding_Tm.md), [DNA_Dimer_Tm](DNA_Dimer_Tm.md)

## 8. References

1. SantaLucia J Jr, Hicks D. 2004. The Thermodynamics of DNA Structural Motifs. Annu Rev Biophys
   Biomol Struct 33:415–440. https://doi.org/10.1146/annurev.biophys.32.110601.141800
2. Untergasser A et al. 2012. Primer3 — new capabilities and interfaces. Nucleic Acids Res 40:e115.
   https://doi.org/10.1093/nar/gks596
3. libprimer3 `thal.c` and `primer3_config/{triloop,tetraloop,loops}.{dh,ds}`, vendored in primer3-py
   (GPL-2.0). https://github.com/libnano/primer3-py/tree/master/primer3/src/libprimer3
