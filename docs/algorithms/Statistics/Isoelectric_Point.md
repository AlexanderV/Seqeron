# Isoelectric Point (pI) Calculation

| Field | Value |
|-------|-------|
| Algorithm Group | Statistics |
| Test Unit ID | SEQ-PI-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

The isoelectric point (pI) of a protein is the pH at which its net charge is zero. This algorithm computes pI from the amino-acid composition using the Henderson–Hasselbalch net-charge model with the EMBOSS Epk.dat pKa scale [1][2]. It is a deterministic, specification-driven calculation: the net charge as a function of pH is a smooth, monotonically non-increasing curve, and the pH at which it crosses zero is located by bisection over the standard window [0, 14] [1]. The model assumes each ionizable group titrates independently (no electrostatic coupling), so pI depends only on residue composition, not order [1].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Each protein carries ionizable groups: the N-terminal amino group and basic side chains (Arg, Lys, His) are positively charged when protonated; the C-terminal carboxyl group and acidic side chains (Asp, Glu, Cys, Tyr) are negatively charged when deprotonated [1][2]. The fraction of each group that is charged at a given pH follows the Henderson–Hasselbalch relation, governed by the group's pKa. The pI is the pH at which positive and negative contributions exactly cancel.

### 2.2 Core Model

Net charge at pH, summed over all ionizable groups [2]:

- Basic group (N-terminus, R, K, H): contributes `+1 / (1 + 10^(pH − pKa))`
- Acidic group (C-terminus, D, E, C, Y): contributes `−1 / (1 + 10^(pKa − pH))`

The isoelectric point is the pH where the total net charge equals zero [1]. Because the net-charge function is strictly decreasing in pH, the root is unique in [0, 14] and is found by bisection.

EMBOSS Epk.dat pKa values [1]: N-terminus 8.6, C-terminus 3.6, Cys 8.5, Asp 3.9, Glu 4.1, His 6.5, Lys 10.8, Arg 12.5, Tyr 10.1.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | No electrostatic interactions between groups; each group titrates independently [1] | Predicted pI deviates from experimental pI for proteins with strong charge coupling |
| ASM-02 | A single pKa per residue type (EMBOSS scale), independent of sequence context | Differs from position-dependent models (e.g., Bjellqvist), giving a slightly different pI [4] |
| ASM-03 | Both termini always present and ionizable | pI is undefined for a zero-length sequence (handled as an input guard) |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ pI ≤ 14 | Bisection is confined to [0, 14] [1] |
| INV-02 | pI is composition-only: permuting residues leaves pI unchanged | Charge is summed over residue counts, not positions [1] |
| INV-03 | Net charge is monotonically non-increasing in pH | Each Henderson–Hasselbalch term is non-increasing in pH [2] |
| INV-04 | Termini-only sequence → pI = (8.6 + 3.6) / 2 = 6.10 | Only the two terminal terms contribute; they cancel at the pKa midpoint [1] |

### 2.5 Comparison with Related Methods

| Aspect | EMBOSS scale (this) | Bjellqvist / ExPASy |
|--------|---------------------|---------------------|
| pKa parameterization | One pKa per residue type | Position-dependent pKa for terminal residues (17 parameters) [4] |
| pI of `ACDEFGHIKLMNPQRSTVWY` | 7.36 | 6.78454 [3] |
| Source | EMBOSS Epk.dat [1] | Bjellqvist et al. 1993 [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| proteinSequence | string | required | Single-letter amino-acid sequence | Case-insensitive; non-ionizable residues ignored; null/empty → sentinel |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | double | Isoelectric point in [0, 14], rounded to 2 decimal places |

### 3.3 Preconditions and Validation

Input is uppercased (case-insensitive). Null or empty returns the neutral sentinel 7.0 (pI is undefined for a zero-length protein; see ASM-03). Only the nine ionizable groups (7 side chains + 2 termini) contribute; any other character is ignored, so non-standard residues, gaps, or whitespace do not throw. No exceptions are raised for valid string input.

## 4. Algorithm

### 4.1 High-Level Steps

1. Guard: null/empty → return 7.0.
2. Count occurrences of each ionizable side-chain residue (D, E, C, Y, H, K, R).
3. Bisect pH over [0, 14]: at each midpoint compute net charge (termini + side chains); if positive, search higher pH, else lower.
4. Stop when the interval is below the precision threshold (0.01); round the pH to 2 decimals.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

EMBOSS Epk.dat pKa table (origin: EMBOSS iep [1]):

| Group | pKa | Sign |
|-------|-----|------|
| N-terminus | 8.6 | + |
| C-terminus | 3.6 | − |
| Cys (C) | 8.5 | − |
| Asp (D) | 3.9 | − |
| Glu (E) | 4.1 | − |
| His (H) | 6.5 | + |
| Lys (K) | 10.8 | + |
| Arg (R) | 12.5 | + |
| Tyr (Y) | 10.1 | − |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateIsoelectricPoint | O(n) | O(1) | One O(n) pass to count residues; bisection is a fixed number of iterations (≈ log2(14/0.01) ≈ 11), each over a constant-size alphabet |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceStatistics.CalculateIsoelectricPoint(string)`: public entry point; counts residues and bisects for the zero-charge pH.
- `SequenceStatistics.NetCharge(...)` (private): Henderson–Hasselbalch net charge at a given pH.

### 5.2 Current Behavior

pKa values are stored as named constants and a dictionary keyed by residue, each annotated with its EMBOSS source. The bisection runs to a fixed pH precision (0.01) and the result is rounded to 2 decimals, matching the resolution at which the pKa scale is meaningful. This is not a search/matching unit, so the repository suffix tree is not applicable (N/A).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- EMBOSS Epk.dat pKa values for the seven ionizable side chains and both termini [1].
- Henderson–Hasselbalch net-charge formula: basic `+1/(1+10^(pH−pKa))`, acidic `−1/(1+10^(pKa−pH))` [2].
- pI = pH where net charge crosses zero, located over [0, 14] [1].

**Intentionally simplified:**

- Single pKa per residue (EMBOSS scale) rather than the position-dependent Bjellqvist parameterization; **consequence:** pI differs slightly from ExPASy/Bjellqvist for the same sequence (e.g., 7.36 vs 6.78454 for the all-20 sequence) [3][4].

**Not implemented:**

- Phosphorylation / PTM charges and non-standard residue pKa; **users should rely on:** specialized tools (e.g., IPC, pIChemiSt) for modified peptides.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Empty/null → 7.0 | Assumption | pI undefined for zero-length protein | accepted | Input-guard sentinel; ASM-03 |
| 2 | EMBOSS pKa scale | Assumption | pI differs from Bjellqvist | accepted | ASM-02; chosen to match single-pKa model |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty | 7.0 | Input guard; pI undefined for zero-length (ASM-03) |
| Termini-only (`A`, `AG`) | 6.10 | Midpoint of N/C-term pKa (INV-04) |
| Acidic-only (`DDDD`) | low pI (3.23) | Acidic side chains pull pI down |
| Basic-only (`KKKK`) | high pI (11.27) | Basic side chains push pI up |
| Lowercase input | same as uppercase | Case-insensitive normalization |

### 6.2 Limitations

Predicted pI is a composition-based estimate that ignores 3-D structure and electrostatic coupling; it is less reliable for very small or highly basic proteins [5]. PTMs and non-standard residues are not modeled. Scale choice (EMBOSS vs Bjellqvist) shifts the result; this implementation uses the EMBOSS scale.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
double pi = SequenceStatistics.CalculateIsoelectricPoint("FLPVLAGLTPSIVPKLVCLLTKKC");
// pi ≈ 9.67 — a basic peptide (net charge remains +0.72 at pH 9 on the EMBOSS scale)
```

**Numerical walk-through:** For `A` (no ionizable side chains), net charge = `+1/(1+10^(pH−8.6)) − 1/(1+10^(3.6−pH))`. This is zero when `pH − 8.6 = −(3.6 − pH)`, i.e. pH = (8.6 + 3.6)/2 = 6.10.

### 7.2 Applications and Use Cases

- **2-D gel electrophoresis / IEF:** predicting the focusing position of a protein along an immobilized pH gradient [4].
- **Protein purification:** choosing buffer pH for ion-exchange chromatography based on predicted net charge.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_CalculateIsoelectricPoint_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateIsoelectricPoint_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Evidence: [SEQ-PI-001-Evidence.md](../../../docs/Evidence/SEQ-PI-001-Evidence.md)

### 7.4 Change History

| Date | Version | Changes |
|------|---------|---------|
| 2026-06-13 | 1.0 | Initial doc; pKa values corrected to EMBOSS Epk.dat scale (SEQ-PI-001) |

## 8. References

1. EMBOSS. iep — Calculate the isoelectric point of proteins. EMBOSS application documentation. https://emboss.sourceforge.net/emboss/apps/iep.html
2. Osorio D, Rondón-Villarreal P, Torres R. 2015. Peptides: A Package for Data Mining of Antimicrobial Peptides. The R Journal 7(1):4–14. Source `src/charge_pI.cpp`. https://github.com/cran/Peptides/blob/master/src/charge_pI.cpp
3. Charif D, Lobry JR. seqinr — computePI. CRAN documentation. https://rdrr.io/cran/seqinr/man/computePI.html
4. Bjellqvist B, Hughes GJ, Pasquali C, et al. 1993. The focusing positions of polypeptides in immobilized pH gradients can be predicted from their amino acid sequences. Electrophoresis 14:1023–1031. https://doi.org/10.1002/elps.11501401163
5. ExPASy. Compute pI/Mw documentation. https://web.expasy.org/compute_pi/pi_tool-doc.html
