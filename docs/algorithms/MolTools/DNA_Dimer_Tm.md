# DNA Self-Dimer / Hetero-Dimer Tm (Thermodynamic Alignment)

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools / Nucleic-acid thermodynamics |
| Test Unit ID | PRIMER-TM-001 (self-/hetero-dimer extension) |
| Related Projects | Seqeron.Genomics.MolTools |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-25 |

## 1. Overview

Computes the melting temperature of the most stable **intermolecular** duplex two DNA
oligonucleotides can form — a **self-dimer** (an oligo against a second copy of itself) or a
**hetero/cross-dimer** (two different oligos). It performs a Primer3 / `ntthal`-style
thermodynamic alignment over the SantaLucia & Hicks (2004) unified nearest-neighbour (NN) model:
it slides one strand against the other over every gapless antiparallel offset, sums the NN
stacking, terminal-A·T and duplex-initiation contributions of each contiguous Watson-Crick run,
and returns the duplex with the highest bimolecular Tm [1][2]. It is an **opt-in** design tool used
to flag primer self-/cross-dimers; the perfect-match duplex Tm, the hairpin Tm and the default
Wallace/Marmur-Doty Tm are unchanged.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

PCR primers and probes can hybridize to themselves or to each other instead of the template,
lowering amplification efficiency. The stability of such a dimer is governed by the same
nearest-neighbour duplex thermodynamics as a perfect duplex, applied to whatever contiguous
Watson-Crick helix the two strands can form when aligned antiparallel [1][2].

### 2.2 Core Model

For a contiguous Watson-Crick duplex of `L` base pairs (`L−1` NN stacks) [1]:

- ΔH° = ΔH°_init + Σ ΔH°(stack) + Σ_ends ΔH°_AT-penalty
- ΔS° = ΔS°_init + Σ ΔS°(stack) + Σ_ends ΔS°_AT-penalty + 0.368·(L−1)·ln[Na⁺]

with ΔH°_init = +0.2 kcal/mol, ΔS°_init = −5.7 cal/(K·mol); the terminal A·T penalty
(ΔH° = +2.2, ΔS° = +6.9) is applied once per duplex end that closes with an A·T pair; the
0.368 entropy salt correction (Eq. 5) is applied per NN stack [1]. The bimolecular Tm is [1]:

`Tm = ΔH°·1000 / (ΔS° + R·ln(C_T / x)) − 273.15`, R = 1.9872 cal/(K·mol),

where x = 1 when **both** oligos are reverse-complement palindromes and x = 4 otherwise [1][3].
The most stable duplex is the contiguous WC run (over all antiparallel offsets) with the highest Tm
— the `ntthal` selection rule [2][3].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Two-state, fixed-buffer NN thermodynamics at the SantaLucia 1 M reference, salt-corrected via Eq. 5 | Off in extreme salt / non-two-state transitions |
| ASM-02 | `FindMostStableDimer` scores a single gapless contiguous WC helix; the full `CalculateDimerThermodynamicsNtthal` DP additionally models internal mismatches/loops, bulges and terminal overhangs (`tstack2`) | The contiguous scorer underestimates stability of overhang/loop-stabilized dimers; the full DP does not |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | ΔH° is salt-independent and equals init + Σ stacks + A·T penalty | Salt enters only ΔS° (Eq. 5) [1] |
| INV-02 | x = 1 iff both oligos are reverse-complement palindromes, else x = 4 | ntthal `symmetry_thermo` + RC branch [3] |
| INV-03 | No WC duplex of ≥ 2 bp ⇒ no result (null / NaN) | A duplex needs ≥ 1 NN stack [2] |
| INV-04 | Lower [Na⁺] strictly lowers Tm for a fixed duplex | 0.368·(L−1)·ln[Na⁺] makes ΔS° more negative [1] |

### 2.5 Comparison with Related Methods

| Aspect | This method | `HasPrimerDimer` (legacy) |
|--------|-------------|---------------------------|
| Output | ΔH°/ΔS°/ΔG°37 + bimolecular Tm | boolean complementarity flag |
| Basis | SantaLucia NN thermodynamics | run-length complementarity heuristic |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| strand1 / strand2 | string | required | DNA oligos (5'→3'); same string twice ⇒ self-dimer | ≥ 2 ACGT bases, case-insensitive |
| sodiumMolar | double | 0.05 (50 mM) | [Na⁺] in mol/L; only the entropy salt term depends on it | > 0 |
| strandConcentrationMolar | double | 50e-9 (50 nM) | total strand concentration C_T (Primer3/ntthal convention) | > 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `DimerResult.Strand1Start` / `Strand2Start` | int | 0-based 5' indices of the aligned duplex on each strand |
| `DimerResult.BasePairs` | int | contiguous Watson-Crick base pairs in the duplex |
| `DimerResult.DeltaH` / `DeltaS` / `DeltaG37` | double | ΔH° (kcal/mol), ΔS° (cal/(K·mol), salt-corrected), ΔG°37 (kcal/mol) |
| Tm methods | double | bimolecular Tm in °C, or `NaN` if no dimer / invalid input |

### 3.3 Preconditions and Validation

0-based indexing; DNA alphabet A/C/G/T only (case-insensitive, upper-cased internally).
Null/empty, < 2-base, or non-ACGT input → `null` (`FindMostStableDimer`) / `NaN` (Tm methods).
A pair with no Watson-Crick duplex of ≥ 2 contiguous bp → `null` / `NaN`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate both strands (ACGT, ≥ 2 bases); decide x (1 if both palindromic else 4).
2. Read strand 2 in 3'→5' order so its index pairs base-for-base under strand 1 (5'→3').
3. For every gapless antiparallel offset, split the overlap into maximal contiguous WC runs.
4. Score each run of ≥ 2 bp: ΔH°/ΔS° = init + Σ stacks + A·T penalty + 0.368·(L−1)·ln[Na⁺]; Tm.
5. Keep the run with the highest Tm; return its `DimerResult` (and Tm from the Tm methods).

### 4.2 Decision Rules, Scoring, Reference Tables

Reuses the repository's existing source-cited tables/constants (no new parameters): the
SantaLucia (1998) unified NN ΔH°/ΔS° (`NnUnifiedParams`), duplex initiation (+0.2, −5.7),
terminal A·T penalty (+2.2, +6.9), the 0.368 entropy salt coefficient and R = 1.9872 — all
from SantaLucia & Hicks (2004) Table 1 [1], cross-checked against Primer3 `thal.c` [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindMostStableDimer(s1, s2)` | O(n·m) | O(1) extra | n=|s1|, m=|s2|; one linear scan per offset |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PrimerDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs)

- `PrimerDesigner.FindMostStableDimer(strand1, strand2, sodiumMolar, strandConcentrationMolar)`: the thermodynamic alignment; returns `DimerResult?`.
- `PrimerDesigner.CalculateDimerMeltingTemperature(strand1, strand2, sodiumMolar, strandConcentrationMolar)`: bimolecular Tm of the most stable dimer.
- `PrimerDesigner.CalculateSelfDimerMeltingTemperature(sequence, sodiumMolar, strandConcentrationMolar)`: self-dimer convenience wrapper.

### 5.2 Current Behavior

The dimer C_T defaults to **50 nM** (Primer3/ntthal `dna_conc`), which differs from the
0.5 µM convention of the single-strand primer Tm (`CalculateMeltingTemperatureNN`); pass an
explicit `strandConcentrationMolar` to change it. The alignment is gapless and scores only
maximal contiguous Watson-Crick runs (the most stable one by Tm). A search reuse evaluation was
made: the repository suffix tree was **not** used — dimer scoring is an O(n·m) thermodynamic
(score-based) alignment over offsets, not exact-substring matching, so the suffix tree does not fit.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- SantaLucia & Hicks (2004) unified NN ΔH°/ΔS°, duplex initiation, terminal-A·T penalty, the 0.368 entropy salt correction, and the bimolecular Tm Eq. 3 with x = 1 (both palindromic) / x = 4 [1].
- The Primer3/`ntthal` most-stable (highest-Tm) selection over antiparallel offsets [2][3]; reproduces primer3-py 2.3.0 ΔH°/ΔS°/Tm to machine precision for every pair whose optimum is a contiguous WC duplex (GCGCGCGC, ACGTACGTACGT, ATCGATCGATCG/CGATCGATCGAT, CGATCGATCG, GCATGC, GGGGCCCC, TGCATGCATG/CATGCATGCA).
- **The full `ntthal` dimer DP (mode ANY)** — `CalculateDimerThermodynamicsNtthal` (a verbatim port of `thal.c`: `fillMatrix`/`LSH`/`RSH`/`maxTM`/`calc_bulge_internal`/`traceback`/`calcDimer`) [2][3], with the `tstack2` terminal-stacking table, the `stackmm` internal-mismatch table, the `tstack` internal-loop terminal table, the 5′/3′ `dangle` tables and the interior/bulge loop-length parameters (all verbatim from `primer3_config/*.dh,*.ds`). It scores matched stacks, single internal mismatches, internal loops, single/multi-base bulges and terminal overhangs/dangling ends, and reproduces primer3-py 2.3.0 ΔH°/ΔS°/ΔG°/Tm to machine precision for **non-contiguous** optima (GCGCATGCGC internal loop; GCGCAAAGCGC/GCGCTTTGCGC 3×3 loop; GCGCGCGC/GCGCAGCGC bulge; GCGCGCAAAA/AAAAGCGCGC overhang) as well as all contiguous cases. `CalculateDimerMeltingTemperature`/`CalculateSelfDimerMeltingTemperature` delegate to it.

**Intentionally simplified:**

- Salt: only the Eq. 5 monovalent entropy correction (0.368) is applied in the dimer ΔS°; divalent (Mg²⁺) is not modelled here (use `CalculateMeltingTemperatureNN(..., Owczarzy2008Divalent)` for monomer divalent Tm). **Consequence:** dimer Tm is reported at the monovalent reference only.
- `FindMostStableDimer` (the public `DimerResult` record) keeps its original gapless contiguous-WC scorer for its `BasePairs`/spans/ΔH°/ΔS° fields; the loop/bulge/overhang model is exposed through `CalculateDimerThermodynamicsNtthal` and the Tm methods.

**Not implemented:**

- The optional caller-supplied tri/tetraloop & terminal-mismatch hairpin **special-loop bonus tables** (`triloop`/`tetraloop`) — a hairpin/monomer feature of `ntthal`, not part of the dimer model. **Users should rely on:** Primer3/`ntthal` directly if those hairpin bonuses are needed.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `FindMostStableDimer` reports the contiguous-WC optimum only | Simplification | Its `DimerResult` underestimates overhang/loop-stabilized dimers | accepted | The full DP (`CalculateDimerThermodynamicsNtthal`) and the Tm methods model loops/bulges/overhangs at full ntthal parity |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Poly-A self-dimer | `null` / `NaN` | no Watson-Crick duplex forms [2] |
| Fully non-complementary pair | `null` / `NaN` | no WC run of ≥ 2 bp (INV-03) |
| < 2 bases / non-ACGT / null | `null` / `NaN` | input contract |
| Non-palindromic self-dimer | x = 4 in the Tm | INV-02 |

### 6.2 Limitations

The full `CalculateDimerThermodynamicsNtthal` DP models internal mismatches/loops, bulges and
terminal overhangs at full ntthal parity; the legacy `FindMostStableDimer` record reports the
contiguous-WC optimum only. Both are monovalent-salt, two-state. For divalent buffers use the
monomer divalent Tm path; the only ntthal capability not ported is the optional tri/tetraloop &
terminal-mismatch hairpin bonus tables (a hairpin/monomer feature).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Self-dimer Tm of an EcoRI-palindrome oligo (defaults: 50 mM Na+, 50 nM strand).
double tm = PrimerDesigner.CalculateSelfDimerMeltingTemperature("GCGCGCGC"); // 40.0906 °C

// Hetero/cross-dimer between two primers.
var dimer = PrimerDesigner.FindMostStableDimer("TGCATGCATG", "CATGCATGCA");
// dimer.BasePairs == 10, dimer.DeltaH == -74.1, Tm == 25.6596 °C (x = 4)
```

**Numerical walk-through (GCGCGCGC self-dimer):** 7 stacks = 4·GC(−9.8,−24.4) + 3·CG(−10.6,−27.2);
ΔH° = 0.2 + (−71.0) = −70.8 kcal/mol; ΔS° = −5.7 + (−179.2) + 0.368·7·ln(0.05) = −192.617; no A·T end;
palindrome ⇒ x = 1; Tm = −70800/(−192.617 + 1.9872·ln(50e-9/1)) − 273.15 = 40.0906 °C [1].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PrimerDesigner_DimerTm_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/PrimerDesigner_DimerTm_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [PRIMER-TM-001-DIMER-Evidence.md](../../../docs/Evidence/PRIMER-TM-001-DIMER-Evidence.md)
- TestSpec: [PRIMER-TM-001-DIMER.md](../../../tests/TestSpecs/PRIMER-TM-001-DIMER.md)
- Related algorithms: [DNA_Hairpin_Folding_Tm](./DNA_Hairpin_Folding_Tm.md)

### 7.4 Change History

| Date | Version | Changes |
|------|---------|---------|
| 2026-06-25 | 1.0 | Initial self-/hetero-dimer Tm via thermodynamic alignment |

## 8. References

1. SantaLucia J, Hicks D. 2004. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Annu Rev Biophys Biomol Struct 33:415-440. https://doi.org/10.1146/annurev.biophys.32.110601.141800
2. Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG. 2012. Primer3 — new capabilities and interfaces. Nucleic Acids Res 40(15):e115. https://doi.org/10.1093/nar/gks596
3. Primer3 `thal.c` (ntthal), primer3-py vendored libprimer3. https://raw.githubusercontent.com/libnano/primer3-py/master/primer3/src/libprimer3/thal.c
