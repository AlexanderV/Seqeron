# RNA Free Energy Calculation

| Field | Value |
|-------|-------|
| Algorithm Group | RNA Structure |
| Test Unit ID | RNA-ENERGY-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

RNA free-energy calculation estimates the thermodynamic stability of RNA secondary-structure motifs and whole sequences under the nearest-neighbor model. The repository exposes local energy evaluators for stem stacking and hairpin loops, along with a dynamic-programming routine for minimum free energy (MFE). These routines are grounded in Turner 2004 parameter tables from NNDB and report physical energies in kcal/mol, but internal-loop coverage is partial rather than table-complete. The implementation is simplified relative to full production RNA folding packages because it omits pseudoknot thermodynamics, uses a reduced multibranch-loop treatment, and falls back to generic models for several larger internal-loop classes.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The nearest-neighbor model represents an RNA secondary structure as a composition of stacked base pairs and loop motifs whose thermodynamic contributions are approximately additive. Turner-style parameterizations remain a standard basis for RNA folding predictions at 37°C and standard buffer assumptions (Xia et al., 1998; Mathews et al., 2004; NNDB Turner 2004).

### 2.2 Core Model

At a high level, the free energy is decomposed as:

$$
\Delta G^\circ_{total} = \sum \Delta G^\circ_{stacking} + \sum \Delta G^\circ_{loops} + \sum \Delta G^\circ_{special}
$$

The current implementation includes:

- Stacking energies between adjacent Watson-Crick and wobble pairs.
- Hairpin, internal, bulge, and multibranch loop terms.
- Special cases such as named hairpin loops, terminal mismatch contributions, AU/GU penalties, and selected coaxial effects.

For the global MFE routine, the repository uses a Zuker-style dynamic program over subsequences, documented in the source as `W`, `V`, and `WM` recurrences with Turner 2004 parameters.

### 2.3 Modeling Assumptions (Optional)

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Thermodynamic parameters are interpreted at 37°C under the Turner 2004 / NNDB parameter regime | Energies can shift under different temperatures or ionic conditions |
| ASM-02 | The folding model is restricted to nested secondary structures for the MFE routine | Crossing interactions such as pseudoknots are outside the computed thermodynamic model |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Hairpin loops shorter than three nucleotides are disallowed in the MFE path | The implementation clamps `minLoopSize` to `3` and treats shorter generic hairpins as prohibitive |
| INV-02 | `CalculateStemEnergy` returns `0` when no base pairs are supplied | The stem-energy helper exits early on an empty base-pair list |
| INV-03 | Energies are rounded to two decimal places in the public helpers | The energy routines call `Math.Round(..., 2)` before returning |
| INV-04 | More negative energy corresponds to greater predicted stability | This is the nearest-neighbor thermodynamic interpretation used throughout the cited RNA folding literature |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[CalculateStemEnergy] sequence` | `string` | required | RNA sequence that provides the paired bases | Used only to interpret the supplied base-pair list |
| `[CalculateStemEnergy] basePairs` | `IReadOnlyList<BasePair>` | required | Ordered base-pair list defining the stem | Empty list returns `0` |
| `[CalculateHairpinLoopEnergy] loopSequence` | `string` | required | Unpaired loop nucleotides | Generic loops shorter than `3` receive prohibitive energy unless handled as a special loop |
| `[CalculateHairpinLoopEnergy] closingBase5` | `char` | required | 5' closing-pair base |
| `[CalculateHairpinLoopEnergy] closingBase3` | `char` | required | 3' closing-pair base |
| `[CalculateHairpinLoopEnergy] specialGUClosure` | `bool` | `false` | Enables the special `G-U` closure bonus when the documented motif applies | Used by the helper's sequence-dependent scoring branch |
| `[CalculateMinimumFreeEnergy] rnaSequence` | `string` | required | RNA sequence to fold | Empty or too-short sequences return `0` |
| `[CalculateMinimumFreeEnergy] minLoopSize` | `int` | `3` | Minimum allowed hairpin loop size | Values below `3` are clamped to `3` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `[CalculateStemEnergy] return value` | `double` | Stem stacking energy including terminal AU/GU penalties |
| `[CalculateHairpinLoopEnergy] return value` | `double` | Hairpin-loop free energy including sequence-dependent bonuses and penalties |
| `[CalculateMinimumFreeEnergy] return value` | `double` | Minimum free energy of the full RNA sequence under the implemented nested-structure model |

### 3.3 Preconditions and Validation

The public energy helpers uppercase sequence characters internally where needed. `CalculateMinimumFreeEnergy` clamps `minLoopSize` to `3`, returns `0` for empty or too-short sequences, and uses a constant `MAXLOOP = 30` for internal and bulge-loop enumeration. `CalculateHairpinLoopEnergy` returns a high prohibitive value for generic loops shorter than three nucleotides, reflecting the steric impossibility of such hairpins in the cited Turner/NNDB model.

## 4. Algorithm

### 4.1 High-Level Steps

1. For stem energy, walk adjacent base pairs and add the corresponding Turner 2004 stacking energies, including special context and terminal AU/GU penalties.
2. For hairpin energy, check special named loops first, then add initiation, mismatch, bonus, and penalty terms from the Turner 2004 tables.
3. For MFE, fill dynamic-programming states over increasing subsequence spans.
4. Evaluate hairpin, stack, internal/bulge, and multibranch alternatives for paired states.
5. Combine local motif energies into the global minimum free energy for the full sequence.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The implementation embeds Turner 2004 / NNDB tables for stacking, hairpin initiation, terminal mismatches, dangling ends, selected internal-loop classes, bulges, and multibranch parameters. It includes explicit handling for special hairpin loops, all-`C` penalties, terminal AU/GU penalties, and the documented `GGUC/CUGG` three-stack context. Internal-loop support is exact for 1x1 loops and then falls back to generic mismatch-based models for larger classes such as 2x3 and other wider internal loops; the source and test specifications explicitly note that dedicated int21 and int22 tables are not implemented. The MFE routine uses pooled `W`, `V`, and `WM` arrays, and its multibranch term currently sets the free-base cost to `0.0` as a simplification.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateStemEnergy` | `O(b)` | `O(1)` | `b` = number of base pairs in the supplied stem |
| `CalculateHairpinLoopEnergy` | `O(l)` | `O(1)` | `l` = loop length |
| `CalculateMinimumFreeEnergy` | `O(n^3)` | `O(n^2)` | `n` = RNA sequence length; internal/bulge enumeration is bounded by `MAXLOOP = 30` |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.CalculateStemEnergy(string, IReadOnlyList<BasePair>)`: Computes stem stacking energy.
- `RnaSecondaryStructure.CalculateHairpinLoopEnergy(string, char, char, bool)`: Computes hairpin-loop energy from Turner 2004 tables and sequence-dependent modifiers.
- `RnaSecondaryStructure.CalculateMinimumFreeEnergy(string, int)`: Computes sequence-level MFE by dynamic programming.

### 5.2 Current Behavior

The repository currently uses Turner 2004 tables sourced from NNDB and rounds returned energies to two decimals. The MFE routine is more detailed than the older document text implied: the source documents a Zuker-style dynamic program with hairpin, stacking, internal/bulge, and multibranch alternatives, plus dangling-end and AU/GU penalties. A classic baseline method is still retained internally as `CalculateMinimumFreeEnergyClassic`, but the public API uses the more detailed routine. Pseudoknot thermodynamics remain out of scope for the public MFE path.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Turner 2004 nearest-neighbor stacking and hairpin-parameter usage, plus selected internal-loop tables.
- Hairpin special-loop handling and sequence-dependent mismatch bonuses/penalties.
- Dynamic-programming MFE search over nested secondary structures.

**Intentionally simplified:**

- The multibranch-loop free-base term is fixed to `0.0` in the current MFE routine; **consequence:** multibranch energetics are less detailed than a full Turner implementation.
- Internal-loop coverage is partial: 1x1 loops use dedicated tables, while larger classes fall back to generic mismatch-based models; **consequence:** the implementation does not fully reproduce Turner/NNDB int21 and int22 parameterization.

**Not implemented:**

- Pseudoknot thermodynamics in the MFE routine; **users should rely on:** no current alternative in this class.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty base-pair list for `CalculateStemEnergy` | Returns `0` | Early-return branch in the helper |
| Empty or too-short sequence for `CalculateMinimumFreeEnergy` | Returns `0` | No valid loop or helix can be formed |
| `minLoopSize < 3` | Clamped to `3` | Turner rules prohibit shorter hairpin loops |
| Special tetraloops such as GNRA or UNCG | Return table-driven special energies rather than generic hairpin energy | The code checks the special-loop table first |

### 6.2 Limitations

The current implementation assumes Turner 2004 conditions at 37°C and standard buffer assumptions. It does not model pseudoknots, modified bases, or the full range of refinements present in mature folding packages such as ViennaRNA. Internal-loop coverage is also incomplete for larger classes because some dedicated Turner tables are replaced by generic mismatch models. As a result, the reported energies are useful and source-backed, but still represent a simplified thermodynamic surface.

## 7. Examples and Related Material

- [RNA-ENERGY-001](../../../tests/TestSpecs/RNA-ENERGY-001.md) documents the repository's RNA free-energy test specification.
- [RNA_Secondary_Structure.md](./RNA_Secondary_Structure.md) documents the higher-level structure-prediction API that consumes these energy models.

## 8. References

1. Turner, D. H. 2004. Thermodynamics of RNA. In The RNA World, 3rd Edition.
2. Xia, T., et al. 1998. Thermodynamic parameters for an expanded nearest-neighbor model for formation of RNA duplexes with Watson-Crick base pairs. Biochemistry 37(42):14719-14735.
3. Mathews, D. H., et al. 2004. Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure. Proceedings of the National Academy of Sciences 101(19):7287-7292.
4. NNDB Turner 2004. https://rna.urmc.rochester.edu/NNDB/
5. Wikipedia contributors. Nucleic acid thermodynamics. Wikipedia.

