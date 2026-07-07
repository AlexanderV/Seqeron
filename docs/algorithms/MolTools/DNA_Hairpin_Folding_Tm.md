# DNA Hairpin Folding + Secondary-Structure (Hairpin) Tm

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PRIMER-TM-001 (hairpin / secondary-structure Tm extension) |
| Related Projects | Seqeron.Genomics.MolTools |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-25 |

## 1. Overview

Finds the most stable intramolecular **hairpin** (a single Watson-Crick stem closing one hairpin loop) that a
DNA primer/oligo can self-fold into, computes its ΔH°, ΔS°, and ΔG°37, and reports its **melting temperature**.
A hairpin is intramolecular, so its Tm is concentration-independent: `Tm = ΔH°·1000/ΔS° − 273.15` with no
strand-concentration term [1]. The model reuses the SantaLucia (1998) unified nearest-neighbour stem stacks [2]
and the SantaLucia & Hicks (2004) Table 4 hairpin-loop-initiation increments [1]. It is an **opt-in** addition:
the existing duplex Tm methods (`CalculateMeltingTemperatureNN`, the default Wallace/Marmur-Doty Tm) are
unchanged. The result is exact for the stem-stack + loop-initiation core; the supplementary triloop/tetraloop
bonus and terminal-mismatch increments are not bundled (see §5.3, §6.2).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A single-stranded DNA oligo can fold back on itself when an internal subsequence is complementary to a
downstream subsequence, forming a stem (a short duplex) that closes a loop of unpaired bases — a hairpin.
Hairpins in primers/probes inhibit hybridization and shift the effective melting behaviour, so primer design
tools screen for them [1]. Unlike a bimolecular duplex, a hairpin forms within one molecule.

### 2.2 Core Model

For a hairpin with a stem of paired bases closing a loop of `N` unpaired nucleotides (SantaLucia & Hicks
2004, "Hairpin Loops", Eqs 7–11 [1]):

- **Stem:** sum of nearest-neighbour stacking ΔH°/ΔS° over consecutive stem base-pair steps, using the
  SantaLucia (1998) unified Table 1 NN parameters [1][2]. The bimolecular duplex-initiation term (+0.2/−5.7) is
  **excluded** — it nucleates two separate strands; for a unimolecular hairpin the loop-initiation term is the
  nucleation cost (ASM-01).
- **Loop:** ΔG°37(loop of N) from Table 4 [1]; `ΔH°(loop) = 0`; `ΔS°(loop) = −ΔG°37·1000/310.15` (the loop is
  destabilising, ΔG°37 > 0 → ΔS° < 0) [1, Table 4 footnote a].
- **Total:** `ΔG°37 = Σ stem stacks + ΔG°37(loop)`; `ΔH° = Σ stem ΔH°`; `ΔS° = Σ stem ΔS° + ΔS°(loop)` [1, p.428].
- **Hairpin Tm (Eq.11):** `Tm = ΔH°·1000/ΔS° − 273.15` — **unimolecular, concentration-independent** [1][3].

Hairpin-loop ΔG°37 by size N (kcal/mol, 1 M NaCl, Table 4 [1]): 3→3.5, 4→3.5, 5→3.3, 6→4.0, 7→4.2, 8→4.3,
9→4.5, 10→4.6, 12→5.0, 14→5.1, 16→5.3, 18→5.5, 20→5.7, 25→6.1, 30→6.3. Non-tabulated sizes use the
Jacobson-Stockmayer extrapolation `ΔG°37(n) = ΔG°37(x) + 2.44·R·310.15·ln(n/x)` from the largest tabulated x ≤ n
(Eq.7 [1]).

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The bimolecular duplex-initiation term does not apply to a unimolecular hairpin (loop init is the nucleation cost). | ΔS° would carry an extra −5.7 e.u. and ΔH° +0.2 kcal/mol, biasing Tm. |
| ASM-02 | The terminal-AT penalty is not applied at the open stem end of the hairpin core (Eqs 8–10 add only stem stacks + loop). | A small ΔG°37 shift at A·T-closed open ends. |
| ASM-03 | Loop ΔH° = 0; all loop temperature dependence is entropic. | ΔH°/Tm error if a future loop ΔH° model is required. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The returned hairpin minimises ΔG°37 over all stem/loop placements (MFE). | The folder scans every closing pair, extends each stem maximally, and keeps the minimum ΔG°37. |
| INV-02 | An oligo with no Watson-Crick stem of ≥2 bp closing a ≥3-nt loop returns null. | No candidate satisfies the stem/loop constraints (e.g. homopolymer). |
| INV-03 | Hairpin Tm = ΔH°·1000/ΔS° − 273.15 with no concentration term. | Unimolecular transition (Eq.11 [1]; concentration-independence [3]). |
| INV-04 | Loop ΔH° = 0; loop ΔS° = −ΔG°37·1000/310.15. | Table 4 footnote a [1]. |
| INV-05 | Returned loop size is always ≥ 3. | Loops < 3 nt are sterically prohibited [1]. |
| INV-06 | Stem ΔH°/ΔS° = SantaLucia Table 1 NN stacks summed over the stem (no bimolecular init). | Eq.10 model [1]; ASM-01. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | DNA oligo, 5'→3' | A/C/G/T only (case-insensitive); else null |
| minStemLength | int | 2 | minimum stem base pairs (≥1 NN stack) | must be ≥ 2 |
| loopBonusDeltaG37 | double | 0.0 | opt-in caller-supplied terminal-mismatch / special-loop ΔG°37 increment (kcal/mol) | not bundled |

### 3.2 Output / Return Value

`FindMostStableHairpin` → `HairpinResult?`:

| Field | Type | Description |
|-------|------|-------------|
| StemStart / StemEnd | int | 0-based indices of the outermost stem pair on the input strand |
| StemLength | int | stem length in base pairs |
| LoopSize | int | number of unpaired loop nucleotides |
| DeltaH / DeltaS / DeltaG37 | double | hairpin ΔH° (kcal/mol), ΔS° (cal/(K·mol)), ΔG°37 (kcal/mol) |

`CalculateHairpinMeltingTemperature` → `double`: hairpin Tm in °C, or `NaN` if no hairpin / invalid input.

### 3.3 Preconditions and Validation

Null / empty / non-ACGT / `minStemLength < 2` → `null` (and `NaN` Tm). Sequences are upper-cased; indexing is
0-based; the alphabet is strict A/C/G/T (degenerate bases are rejected, consistent with the duplex NN methods).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate input (non-empty, ACGT-only, `minStemLength ≥ 2`).
2. For every candidate outermost closing pair (i, j) that is Watson-Crick: extend the stem inward as far as
   pairing allows, summing the NN stem stacks (ΔH°/ΔS°).
3. The remaining inner bases form the hairpin loop; reject loops < 3 nt and stems < `minStemLength`.
4. Add the loop initiation (ΔG°37 by size + optional bonus; ΔH° = 0; ΔS° = −ΔG°37·1000/310.15).
5. Keep the hairpin with the minimum total ΔG°37.
6. Tm (Eq.11, unimolecular): `ΔH°·1000/ΔS° − 273.15`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **Stem stacks:** `NnUnifiedParams` (SantaLucia 1998 Table 1 [2], shared with the duplex NN Tm).
- **Hairpin loop ΔG°37 by size:** `HairpinLoopInitiationDeltaG` (SantaLucia & Hicks 2004 Table 4 [1]).
- **Large loops:** Jacobson-Stockmayer (Eq.7 [1], coefficient 2.44).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindMostStableHairpin | O(n²) closing-pair scan, each extended O(n) → O(n³) worst case | O(1) | n = oligo length; primers are short (≤ ~40 nt). |

A suffix tree was **not** used: this is a scoring-based self-fold (thermodynamic minimisation), not exact
substring matching; the suffix tree fits exact-occurrence enumeration, not energy-weighted folding.

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PrimerDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs)

- `PrimerDesigner.FindMostStableHairpin(string, int, double)`: MFE hairpin (stem span/length, loop size, ΔH°/ΔS°/ΔG°37).
- `PrimerDesigner.CalculateHairpinMeltingTemperature(string, int, double)`: unimolecular hairpin Tm (Eq.11).
- `PrimerDesigner.HairpinResult` (record struct): the returned structure.

### 5.2 Current Behavior

The folder restricts structures to a single stem closing one hairpin loop (no bulges, internal loops, or
multibranch). It reuses the duplex NN ΔH°/ΔS° table verbatim for the stem and adds the Table 4 loop term. The
Tm uses Eq.11 with no concentration argument (an intramolecular transition). The optional `loopBonusDeltaG37`
lets a caller add the supplementary terminal-mismatch / triloop-tetraloop increment if they have it.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Stem NN stacks from SantaLucia (1998)/(2004) Table 1 [1][2].
- Hairpin-loop ΔG°37 by size from SantaLucia & Hicks (2004) Table 4 [1], with loop ΔH° = 0 and
  ΔS° = −ΔG°37·1000/310.15.
- Jacobson-Stockmayer large-loop extrapolation (Eq.7, coefficient 2.44) [1].
- Unimolecular hairpin Tm = ΔH°·1000/ΔS° − 273.15 (Eq.11), concentration-independent [1][3].

**Intentionally simplified:**

- Stem nucleation: the bimolecular duplex-initiation and terminal-AT terms are omitted (ASM-01, ASM-02);
  **consequence:** a small ΔG°37/Tm offset relative to a full duplex-end treatment.

**Not implemented:**

- The supplementary triloop (length-3) / tetraloop (length-4) bonus tables and the terminal-mismatch
  increment; **users should rely on:** the opt-in `loopBonusDeltaG37` (caller-supplied) or a full folding tool
  (UNAFold, ViennaRNA, MELTING 5) for those special-loop corrections.
- Bulges, internal loops, multibranch loops, and self-/cross-dimer (intermolecular) structures; **users should
  rely on:** UNAFold / ViennaRNA.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Bimolecular init excluded | Assumption | small ΔG°37/Tm offset | accepted | ASM-01 (unimolecular model) |
| 2 | Terminal-AT penalty omitted | Assumption | small ΔG°37 offset | accepted | ASM-02 |
| 3 | Triloop/tetraloop + terminal-mismatch tables not bundled | Deviation | length-3/4 special loops not auto-corrected | accepted | opt-in `loopBonusDeltaG37` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Homopolymer (poly-A) | null / NaN Tm | no Watson-Crick stem [1] |
| Loop would be < 3 nt | not returned (null if it is the only option) | loops < 3 sterically prohibited [1] |
| Non-ACGT base | null / NaN | strict alphabet, as duplex NN methods |
| null / empty | null / NaN | invalid input |
| minStemLength < 2 | null | a stem needs ≥ 1 NN stack |

### 6.2 Limitations

Single hairpin only — no bulges, internal loops, multibranch, or pseudoknots; no self-dimer/cross-dimer
(intermolecular) Tm; the length-3/4 special-loop bonuses and terminal mismatch are caller-supplied. For a full
secondary-structure energy minimisation use UNAFold, ViennaRNA, or MELTING 5.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var hp = PrimerDesigner.FindMostStableHairpin("GGGCTTTTGCCC");
// hp.StemLength == 4, hp.LoopSize == 4, hp.DeltaH == -25.8,
// hp.DeltaS == -75.48486216346927, hp.DeltaG37 == -2.3883700000000054
double tm = PrimerDesigner.CalculateHairpinMeltingTemperature("GGGCTTTTGCCC"); // 68.6403836682880 °C
```

**Numerical walk-through:** stem GGGC steps GG, GG, GC → ΔH° = −8.0−8.0−9.8 = −25.8; ΔS° = −19.9−19.9−24.4 =
−64.2. Loop of 4: ΔG°37 = 3.5 → ΔS°(loop) = −3.5·1000/310.15 = −11.28486. Total ΔS° = −75.48486;
ΔG°37 = −25.8 − 310.15·(−75.48486)/1000 = −2.38837; Tm = −25.8·1000/−75.48486 − 273.15 = 68.6404 °C.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PrimerDesigner_HairpinTm_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/PrimerDesigner_HairpinTm_Tests.cs) — covers `INV-01`–`INV-06`
- Evidence: [PRIMER-TM-001-HAIRPIN-Evidence.md](../../../docs/Evidence/PRIMER-TM-001-HAIRPIN-Evidence.md)
- Related algorithms: [NearestNeighbor_Salt_Corrected_Tm](./NearestNeighbor_Salt_Corrected_Tm.md)

## 8. References

1. SantaLucia J Jr, Hicks D. 2004. The Thermodynamics of DNA Structural Motifs. Annu Rev Biophys Biomol Struct 33:415–440. https://doi.org/10.1146/annurev.biophys.32.110601.141800
2. SantaLucia J Jr. 1998. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proc Natl Acad Sci USA 95:1460–1465. https://doi.org/10.1073/pnas.95.4.1460
3. Vallone PM, Benight AS. 1999. Melting studies of short DNA hairpins: influence of loop sequence and adjoining base pair identity on hairpin thermodynamic stability. Biochemistry. https://pubmed.ncbi.nlm.nih.gov/10423551/
