# Isoelectric Point (pI) Calculation

| Field | Value |
|-------|-------|
| Algorithm Group | Statistics |
| Test Unit ID | SEQ-PI-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-09-28 |

## 1. Overview

The isoelectric point (pI) of a protein is the pH at which its net charge is zero. This algorithm computes pI with the Henderson–Hasselbalch net-charge model on one of two published pK sets: the **EMBOSS `iep`** set (default; EMBOSS 6.6.0 `Epk.dat`) [1][6] or the **Bjellqvist** set used by ExPASy Compute pI/Mw and Biopython (with N-/C-terminal-residue-specific pKs) [4][5][7]. The net charge as a function of pH is a smooth, monotonically non-increasing curve; the pH at which it crosses zero is located by bisection over [0, 14] to 1e-9 pH and rounded to 2 decimals. Each ionizable group titrates independently (no electrostatic coupling); on the EMBOSS scale pI therefore depends only on composition, on the Bjellqvist scale also on the identity of the two terminal residues.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Each protein carries ionizable groups: the N-terminal amino group and basic side chains (Arg, Lys, His) are positively charged when protonated; the C-terminal carboxyl group and acidic side chains (Asp, Glu, Cys, Tyr) are negatively charged when deprotonated [1][2]. The fraction of each group that is charged at a given pH follows the Henderson–Hasselbalch relation, governed by the group's pKa. The pI is the pH at which positive and negative contributions exactly cancel.

### 2.2 Core Model

Net charge at pH, summed over all ionizable groups [2]:

- Basic group (N-terminus, R, K, H): contributes `+1 / (1 + 10^(pH − pKa))`
- Acidic group (C-terminus, D, E, C, Y): contributes `−1 / (1 + 10^(pKa − pH))`

The isoelectric point is the pH where the total net charge equals zero [1]. Because the net-charge function is strictly decreasing in pH, the root is unique in [0, 14] and is found by bisection.

EMBOSS 6.6.0 `Epk.dat` pKa values [6]: N-terminus **7.5**, C-terminus 3.6, Cys 8.5, Asp 3.9, Glu 4.1, His 6.5, Lys 10.8, Arg 12.5, Tyr 10.1. (The listing printed on the EMBOSS iep web page shows `Amino 8.6`, as does the Peptides R package "EMBOSS" scale; that listing is stale — the same page's worked outputs, LACI_ECOLI 6.8385 and IFNA2_HUMAN 5.7240, are only reproduced with Amino 7.5, verified with the EMBOSS 6.6.0 `iep` binary.) EMBOSS also splits ambiguity codes B → D/N (5.5:4.3) and Z → E/Q (6.0:3.9) by Dayhoff frequency, `(int)(0.5 + n·f)` (`embIepCompC`).

Bjellqvist pK values [4][7]: side chains Cys 9.0, Asp 4.05, Glu 4.45, His 5.98, Lys 10.0, Arg 12.0, Tyr 10.0; N-terminus 7.5 unless the N-terminal residue is A 7.59, M 7.0, S 6.93, P 8.36, T 6.82, V 7.44, E 7.7; C-terminus 3.55 unless the C-terminal residue is D 4.55, E 4.75.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | No electrostatic interactions between groups; each group titrates independently [1] | Predicted pI deviates from experimental pI for proteins with strong charge coupling |
| ASM-02 | EMBOSS scale: a single pKa per residue type; Bjellqvist scale: terminal-residue-specific N/C-terminus pKs | Scale choice shifts pI (all-20: EMBOSS 6.97 vs Bjellqvist 6.78) |
| ASM-03 | Both termini always present and ionizable | pI is undefined for a zero-length sequence (handled as an input guard) |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ pI ≤ 14 | Bisection is confined to [0, 14] [1] |
| INV-02 | EMBOSS scale: pI is composition-only (permutations equal). Bjellqvist: invariant only for permutations that keep both terminal residues | Charge summed over counts; Bjellqvist terminal pKs depend on the end residues |
| INV-03 | Net charge is monotonically non-increasing in pH | Each Henderson–Hasselbalch term is non-increasing in pH [2] |
| INV-04 | Termini-only sequence → pI = (7.5 + 3.6) / 2 = 5.55 (EMBOSS iep "A" = 5.5500) | Only the two terminal terms contribute; they cancel at the pKa midpoint [6] |

### 2.5 Comparison with Related Methods

| Aspect | EMBOSS scale (this) | Bjellqvist / ExPASy |
|--------|---------------------|---------------------|
| pKa parameterization | One pKa per residue type | Terminal-residue-specific N/C-terminus pKa [4][7] |
| pI of `ACDEFGHIKLMNPQRSTVWY` | 6.97 (iep 6.9681) | 6.78 (Biopython 6.784552; seqinr 6.78454 [3]) |
| Source | EMBOSS 6.6.0 Epk.dat / embiep.c [6] | Bjellqvist et al. 1993/1994 [4]; Biopython [7] |
| API | `CalculateIsoelectricPoint(seq)` or `(seq, PkaScale.Emboss)` | `CalculateIsoelectricPoint(seq, PkaScale.Bjellqvist)` |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| proteinSequence | string | required | Single-letter amino-acid sequence | Case-insensitive; non-ionizable residues ignored; null/empty → sentinel |
| scale | `SequenceStatistics.PkaScale` | `Emboss` | pK set | `Emboss` or `Bjellqvist` |

`CalculateNetCharge(seq, pH, scale = Emboss)` returns the Henderson–Hasselbalch net charge at a pH (EMBOSS `embIepGetCharge` / Biopython `charge_at_pH`); 0 for null/empty.

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | double | Isoelectric point in [0, 14], rounded to 2 decimal places |

### 3.3 Preconditions and Validation

Input is uppercased (case-insensitive). Null or empty returns the neutral sentinel 7.0 (pI is undefined for a zero-length protein; see ASM-03). Only the nine ionizable groups (7 side chains + 2 termini) contribute (plus B/Z via the Dayhoff split on the EMBOSS scale); any other character is ignored, so non-standard residues, gaps, or whitespace do not throw. No exceptions are raised for valid string input.

## 4. Algorithm

### 4.1 High-Level Steps

1. Guard: null/empty → return 7.0.
2. Count occurrences of each ionizable side-chain residue (D, E, C, Y, H, K, R); EMBOSS: add Dayhoff-split B→D, Z→E; Bjellqvist: pick N/C-terminus pK from the first/last residue.
3. Bisect pH over [0, 14]: at each midpoint compute net charge (termini + side chains); if positive, search higher pH, else lower.
4. Stop when the interval is below 1e-9 pH; round the midpoint to 2 decimals (correctly rounded pI).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

EMBOSS 6.6.0 Epk.dat pKa table (origin: EMBOSS source [6]):

| Group | pKa | Sign |
|-------|-----|------|
| N-terminus | 7.5 | + |
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
| CalculateIsoelectricPoint | O(n) | O(1) | One O(n) pass to count residues; bisection is a fixed number of iterations (≈ log2(14/1e-9) ≈ 34), each over a constant-size alphabet |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceStatistics.CalculateIsoelectricPoint(string)`: EMBOSS-scale entry point.
- `SequenceStatistics.CalculateIsoelectricPoint(string, PkaScale)`: scale-selectable pI.
- `SequenceStatistics.CalculateNetCharge(string, double, PkaScale)`: Henderson–Hasselbalch net charge at a given pH.
- MCP `isoelectric_point` (Seqeron.Mcp.Sequence) delegates to the EMBOSS-scale overload.

### 5.2 Current Behavior

pKa values are stored as named constants and a dictionary keyed by residue, each annotated with its EMBOSS source. The bisection runs to 1e-9 pH and the result is rounded to 2 decimals (until 2026-09 the bracket stopped at 0.01 and the midpoint was rounded, which could misround the last digit). This is not a search/matching unit, so the repository suffix tree is not applicable (N/A).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- EMBOSS 6.6.0 Epk.dat pKa values for the seven ionizable side chains and both termini, and the B/Z Dayhoff split [6]; reproduces the `iep` binary.
- Bjellqvist pK set with terminal-residue-specific pKs [4][7]; reproduces Biopython (and the seqinr/ExPASy all-20 value).
- Henderson–Hasselbalch net-charge formula: basic `+1/(1+10^(pH−pKa))`, acidic `−1/(1+10^(pKa−pH))` [2].
- pI = pH where net charge crosses zero, located over [0, 14] [1].

**Intentionally different from the references:**

- Search window [0, 14]: EMBOSS searches [1, 14] (and reports none without a sign change), Biopython [4.05, 12] (returns the window edge for roots outside, e.g. Bjellqvist `DDDD` 3.52, `RRRRRRRR` 12.85). Inside those windows results are identical.
- EMBOSS `iep` options `-disulphides`, `-lysinemodified`, `-notermini`, custom `-pkdata` are not exposed.

**Not implemented:**

- Phosphorylation / PTM charges and non-standard residue pKa; **users should rely on:** specialized tools (e.g., IPC, pIChemiSt) for modified peptides.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Empty/null → 7.0 | Assumption | pI undefined for zero-length protein | accepted | Input-guard sentinel; ASM-03 |
| 2 | Default pK scale = EMBOSS | Assumption | pI differs from Bjellqvist/ExPASy | accepted | ASM-02; Bjellqvist available via `PkaScale.Bjellqvist` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty | 7.0 | Input guard; pI undefined for zero-length (ASM-03) |
| Termini-only (`A`, `AG`) | 5.55 (EMBOSS); `A` 5.57 (Bjellqvist, N-term A 7.59) | Midpoint of N/C-term pKa (INV-04) |
| Acidic-only (`DDDD`) | 3.23 (EMBOSS iep 3.2279); 3.52 (Bjellqvist) | Acidic side chains pull pI down |
| Basic-only (`KKKK`) | 11.28 (EMBOSS iep 11.2772); 10.48 (Bjellqvist) | Basic side chains push pI up |
| B / Z (EMBOSS) | `B` 3.75, `Z` 3.85 | Dayhoff split to D / E (iep) |
| Lowercase input | same as uppercase | Case-insensitive normalization |

### 6.2 Limitations

Predicted pI is a composition-based estimate that ignores 3-D structure and electrostatic coupling; it is less reliable for very small or highly basic proteins [5]. PTMs and non-standard residues are not modeled. Scale choice (EMBOSS vs Bjellqvist) shifts the result; the default is EMBOSS, Bjellqvist is selectable.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
double pi = SequenceStatistics.CalculateIsoelectricPoint("FLPVLAGLTPSIVPKLVCLLTKKC");
// 9.57 — EMBOSS iep 9.5678 (net charge +0.46 at pH 9)
double piB = SequenceStatistics.CalculateIsoelectricPoint("FLPVLAGLTPSIVPKLVCLLTKKC", SequenceStatistics.PkaScale.Bjellqvist);
// 9.39 — Biopython 9.3902
```

**Numerical walk-through:** For `A` (no ionizable side chains), EMBOSS net charge = `+1/(1+10^(pH−7.5)) − 1/(1+10^(3.6−pH))`. This is zero when `pH − 7.5 = −(3.6 − pH)`, i.e. pH = (7.5 + 3.6)/2 = 5.55.

### 7.2 Applications and Use Cases

- **2-D gel electrophoresis / IEF:** predicting the focusing position of a protein along an immobilized pH gradient [4].
- **Protein purification:** choosing buffer pH for ion-exchange chromatography based on predicted net charge.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_CalculateIsoelectricPoint_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceStatistics_CalculateIsoelectricPoint_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Evidence: [SEQ-PI-001-Evidence.md](../../../docs/Evidence/SEQ-PI-001-Evidence.md)

### 7.4 Change History

| Date | Version | Changes |
|------|---------|---------|
| 2026-06-13 | 1.0 | Initial doc; pKa values corrected to EMBOSS Epk.dat scale (SEQ-PI-001) |
| 2026-09-28 | 2.0 | Review 2026-09: EMBOSS N-terminus pKa 8.6 → 7.5 (actual EMBOSS 6.6.0 Epk.dat; confirmed with the `iep` binary); B/Z Dayhoff split; bisection to 1e-9 (correct rounding); Bjellqvist scale with terminal-residue pKs added (`PkaScale`); public `CalculateNetCharge` |

## 8. References

1. EMBOSS. iep — Calculate the isoelectric point of proteins. EMBOSS application documentation. https://emboss.sourceforge.net/emboss/apps/iep.html
2. Osorio D, Rondón-Villarreal P, Torres R. 2015. Peptides: A Package for Data Mining of Antimicrobial Peptides. The R Journal 7(1):4–14. Source `src/charge_pI.cpp`. https://github.com/cran/Peptides/blob/master/src/charge_pI.cpp
3. Charif D, Lobry JR. seqinr — computePI. CRAN documentation. https://rdrr.io/cran/seqinr/man/computePI.html
4. Bjellqvist B, Hughes GJ, Pasquali C, et al. 1993. The focusing positions of polypeptides in immobilized pH gradients can be predicted from their amino acid sequences. Electrophoresis 14:1023–1031. https://doi.org/10.1002/elps.11501401163
5. ExPASy. Compute pI/Mw documentation. https://web.expasy.org/compute_pi/pi_tool-doc.html
6. EMBOSS 6.6.0 source: `emboss/data/Epk.dat`, `nucleus/embiep.c` (https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/data/Epk.dat, …/nucleus/embiep.c); `iep` binary (Debian `emboss` 6.6.0).
7. Biopython 1.88 `Bio/SeqUtils/IsoelectricPoint.py` (Bjellqvist pK set, terminal-residue pKs).
