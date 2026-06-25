# LNA-Adjusted Nearest-Neighbour Melting Temperature

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PROBE-DESIGN-001 |
| Related Projects | Seqeron.Genomics.MolTools |
| Implementation Status | Production |
| Last Reviewed | 2026-06-24 |

## 1. Overview

This is an **opt-in** extension to the library's nearest-neighbour (NN) melting-temperature model
that supports DNA oligonucleotide probes carrying one or more **internal LNA (locked nucleic acid)
substitutions** on one strand. It adds the McTigue, Peterson & Kahn (2004) sequence-dependent
LNA-DNA NN increments (ΔΔH°, ΔΔS°) to the underlying SantaLucia (1998) DNA NN stack, so the design
Tm reflects the substantial duplex stabilization an internal LNA monomer confers [1]. It also
provides a qualitative 3'-MGB (minor-groove binder) probe **design-rule** check from Kutyavin et al.
(2000) [4]. The generic probe designer, the TaqMan rules, and all melting-temperature defaults are
unchanged; this is a separate additive entry point.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

An LNA monomer is a ribonucleotide whose 2'-O and 4'-C are bridged by a methylene group, locking the
sugar in a C3'-endo conformation. Incorporating LNA into a DNA probe greatly increases duplex
thermal stability — "LNA provides the largest known increase in thermal stability of any modified
DNA duplex" — with a magnitude that is strongly sequence-dependent [1]. A minor-groove binder (MGB)
conjugated to a probe's 3' end likewise stabilizes the duplex, allowing shorter, more specific
probes [4].

### 2.2 Core Model

The duplex thermodynamics are the standard NN model (initiation + per-step stacks + terminal-A·T
penalty + self-complementary symmetry) [5]. For each nearest-neighbour step that contains an LNA
base, an **additive increment** is applied [1]:

> ΔH° = ΔH°(DNA NN) + Σ ΔΔH°(LNA step) ; ΔS° = ΔS°(DNA NN) + Σ ΔΔS°(LNA step)

The Tm uses the bimolecular equation Tm = ΔH°·1000 / (ΔS° + R·ln(C_T/x)) − 273.15, R = 1.9872
cal/(K·mol), x = 4 (non-self-complementary) or 1 (self-complementary) [5], identical to the
perfect-match NN Tm.

The 32 LNA-DNA NN increments are keyed by the DNA dinucleotide step and which base of the step is
locked (the 5' base — McTigue X_L N notation — or the 3' base — MX_L) [1]. Values are in cal/mol and
cal/(mol·K) in the source data; the library stores them in kcal/mol (÷1000) [3].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Internal LNA only; McTigue (2004) did not parameterise a terminal LNA. | A terminal (end) LNA has no increment → the result is reported not-computable. |
| ASM-02 | The base DNA model is SantaLucia (1998) unified, not McTigue's own reference DNA NN set. | A ~0.09 °C offset vs the MELTING `mct04` Tm; the increment contribution itself is exact. |
| ASM-03 | Fixed-temperature two-state duplex melting; salt corrections per the chosen mode. | Same caveats as the perfect-match NN model. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | An empty LNA-position set reproduces the perfect-match NN Tm exactly. | No increment is added; the base computation is unchanged [1][3]. |
| INV-02 | Adding an internal LNA raises Tm vs the all-DNA duplex (stabilization). | McTigue (2004): internal LNA substantially stabilizes the duplex [1] (verified for the worked example, +3.84 °C). |
| INV-03 | A terminal or out-of-range LNA position is not computable. | McTigue (2004) has no terminal-LNA parameter [1][3]. |
| INV-04 | The applied increment equals the verbatim McTigue value for the (step, locked-position) key. | Table transcribed verbatim from MELTING `McTigue2004lockedmn.xml` [3]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | DNA oligo (one strand, 5'→3'). | ≥ 2 ACGT bases. |
| lnaPositions | IReadOnlyCollection&lt;int&gt; | required | Zero-based indices of internal LNA monomers. | Each in (0, length−1); order/duplicates tolerated. |
| strandConcentrationMolar | double | 0.5 µM | Total strand concentration C_T. | > 0 |
| sodiumMolar | double | 0.05 M | Monovalent cation concentration. | ≥ 0 |
| magnesiumMolar, dntpMolar | double | 0 | Divalent-mode inputs. | ≥ 0 |
| saltMode | SaltCorrectionMode | Owczarzy2004Monovalent | Salt correction. | enum |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (DeltaH, DeltaS, IsSelfComplementary)? | nullable tuple | LNA-adjusted duplex ΔH° (kcal/mol), ΔS° (cal/(K·mol)), self-complementarity; `null` if not computable. |
| Tm | double | LNA-adjusted Tm in °C; `double.NaN` if not computable. |

### 3.3 Preconditions and Validation

0-based indexing; ACGT only (case-insensitive, upper-cased internally). Null `sequence` → null/NaN;
null `lnaPositions` → `ArgumentNullException`. Empty / < 2 nt / non-ACGT → null/NaN. Any LNA position
≤ 0 or ≥ length−1 (terminal/out-of-range) → null/NaN. The MGB check throws `ArgumentNullException`
on a null probe.

## 4. Algorithm

### 4.1 High-Level Steps

1. Compute the base DNA NN (ΔH°, ΔS°, self-comp) via the existing perfect-match NN model [5].
2. Validate LNA positions: reject any terminal (0 or length−1) or out-of-range index.
3. For each step (i, i+1), if base i is locked add the 5'-locked increment for step seq[i..i+1]; if
   base i+1 is locked add the 3'-locked increment [1].
4. Apply the bimolecular Tm equation with the chosen salt correction [5].

### 4.2 Decision Rules / Reference Tables

The 32 McTigue (2004) LNA-DNA NN increments (ΔΔH° kcal/mol, ΔΔS° cal/(K·mol)) are transcribed
verbatim from MELTING 5 `McTigue2004lockedmn.xml` [3] (16 with the 5' base locked, 16 with the 3'
base locked). MGB design rule: length window 12–20 nt; MGB attached at the 3' end [4].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| LNA NN ΔH°/ΔS° | O(n) | O(1) | n = sequence length; one dictionary lookup per step. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PrimerDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs),
[ProbeDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs)

- `PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(string, IReadOnlyCollection<int>)`: LNA-adjusted ΔH°/ΔS°.
- `PrimerDesigner.CalculateMeltingTemperatureNNLna(string, IReadOnlyCollection<int>, …)`: LNA-adjusted Tm.
- `ProbeDesigner.EvaluateMgbProbeDesign(string)`: qualitative 3'-MGB design-rule check.

### 5.2 Current Behavior

The LNA methods reuse `CalculateNearestNeighborThermodynamics` (SantaLucia 1998 unified) for the base
duplex, then add increments per LNA-containing step, exactly as the MELTING reference implementation
combines them (DNA NN sum, then `enthalpy += lockedAcidValue`) [3]. The MGB method reports only
length-window conformance and 3'-attachment guidance; it does **not** compute a quantitative MGB ΔTm.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- All 32 McTigue (2004) LNA-DNA NN ΔΔH°/ΔΔS° increments, applied additively per LNA step [1][3].
- Internal-only LNA (terminal positions rejected) [1][3].
- The bimolecular NN Tm equation and salt corrections (reused, unchanged) [5].
- MGB 3'-attachment and 12–20mer design rules [4].

**Intentionally simplified:**

- Base DNA NN model is SantaLucia (1998), not McTigue's reference set; **consequence:** ~0.09 °C
  offset vs MELTING `mct04` (the LNA increments themselves are exact).

**Not implemented:**

- Quantitative MGB ΔTm; **users should rely on:** the empirical/proprietary MGB-Eclipse model or a
  chemistry-specific tool — Kutyavin (2000) gives no closed-form MGB ΔTm [4] (honest residual).
- LNA mismatch-discrimination parameters and consecutive/terminal-LNA models (e.g. IDT 2012);
  **users should rely on:** the dedicated chemistry tool.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Base DNA NN = SantaLucia 1998 | Assumption | ~0.09 °C vs MELTING mct04 | accepted | ASM-02 |
| 2 | Quantitative MGB ΔTm not modelled | Deviation | no MGB Tm value | accepted | empirical only [4]; see LIMITATIONS.md |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty LNA-position set | equals perfect-match NN Tm | INV-01 |
| Terminal LNA (index 0 / last) | null / NaN | INV-03; no McTigue parameter |
| Out-of-range / duplicate / unsorted index | out-of-range → null; dup/order tolerated | set semantics |
| non-ACGT base | null / NaN | base NN lookup fails |
| MGB length outside 12–20 | flagged in guidance | Kutyavin (2000) window |

### 6.2 Limitations

Quantitative MGB ΔTm is not computed (empirical, no published model) [4]. Only single/independent
internal LNA increments are modelled (additive); cooperative consecutive-LNA effects beyond the NN
sum are out of scope. Two-state melting assumption inherited from the NN model.

## 7. Examples and Related Material

### 7.1 Worked Example

`CCATTGCTACC` with an LNA at index 4 (the second T), C_T = 1e-4 M, [Na⁺] = 1 M, no salt correction:
base DNA ΔH° = −80.8, ΔS° = −221.7; add `TTL/AA`(+2.326, +8.1) and `TLG/AC`(−1.540, −3.0) →
ΔH° = −80.014 kcal/mol, ΔS° = −216.6 cal/(K·mol) → **Tm = 63.528 °C** (MELTING `mct04`: 63.614 °C).
The all-DNA Tm is 59.692 °C, so the single internal LNA raises Tm by **+3.84 °C** [1][3].

```csharp
double tm = PrimerDesigner.CalculateMeltingTemperatureNNLna(
    "CCATTGCTACC", new[] { 4 },
    strandConcentrationMolar: 1e-4, sodiumMolar: 1.0,
    saltMode: PrimerDesigner.SaltCorrectionMode.None); // 63.528 °C
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ProbeDesigner_LnaTm_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_LnaTm_Tests.cs) — covers INV-01..INV-04
- Evidence: [PROBE-DESIGN-001-LNA-Evidence.md](../../../docs/Evidence/PROBE-DESIGN-001-LNA-Evidence.md)
- Related algorithms: [NearestNeighbor_Salt_Corrected_Tm](NearestNeighbor_Salt_Corrected_Tm.md), [Hybridization_Probe_Design](Hybridization_Probe_Design.md)

## 8. References

1. McTigue PM, Peterson RJ, Kahn JD. 2004. Sequence-dependent thermodynamic parameters for locked nucleic acid (LNA)-DNA duplex formation. Biochemistry 43(18):5388–5405. https://doi.org/10.1021/bi035976d
2. Dumousseau M, Rodriguez N, Juty N, Le Novère N. 2012. MELTING, a flexible platform to predict the melting temperatures of nucleic acids. BMC Bioinformatics 13:101. https://pmc.ncbi.nlm.nih.gov/articles/PMC3733425/
3. MELTING 5 data file `McTigue2004lockedmn.xml` / `McTigue04LockedAcid.java` (reference implementation; verbatim McTigue 2004 increments). https://github.com/aravind-j/rmelting (inst/extdata/Data). Retrieved 2026-06-24.
4. Kutyavin IV, Afonina IA, Mills A, et al. 2000. 3'-Minor groove binder-DNA probes increase sequence specificity at PCR extension temperatures. Nucleic Acids Res 28(2):655–661. https://doi.org/10.1093/nar/28.2.655
5. SantaLucia J. 1998. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. PNAS 95(4):1460–1465. https://www.pnas.org/doi/10.1073/pnas.95.4.1460
