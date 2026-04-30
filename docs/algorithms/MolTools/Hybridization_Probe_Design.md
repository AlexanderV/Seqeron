# Hybridization Probe Design

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PROBE-DESIGN-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Hybridization probe design generates oligonucleotide probes that detect specific nucleic acid sequences through complementary binding. In this repository, probe design supports multiple applications including FISH, DNA microarrays, Northern blots, qPCR, and Southern blots. The documented implementation is a heuristic candidate-generation and ranking workflow driven by GC content, melting temperature, self-complementarity, secondary structure, simple repeats, and application-specific parameter ranges.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Hybridization probes are used in sequence-detection assays such as FISH, DNA microarrays, Northern blots, and qPCR. Probe quality depends on a balance among target complementarity, melting temperature, GC content, and avoidance of self-structure or repetitive motifs that can reduce specificity. The original document ties these considerations to nucleic-acid thermodynamics and to application-specific probe lengths and Tm targets. Sources: Wikipedia (Hybridization probe, Fluorescence in situ hybridization, DNA microarray, Nucleic acid thermodynamics), SantaLucia (1998), Breslauer et al. (1986).

### 2.2 Core Model

The candidate score starts at `1.0` and is reduced by a fixed set of penalties:

| Factor | Penalty |
|--------|---------|
| GC content outside the configured range | `-0.3` |
| Tm outside the configured range | `-0.3` |
| Homopolymer run above the configured maximum | `-0.2` |
| Self-complementarity above the configured maximum | `-0.2` |
| Secondary-structure potential | `-0.15` |
| Simple repeats | `-0.1` |
| Starts with `G/C` | `-0.02` |
| Ends with `G/C` | `-0.02` |

The original document describes the thermodynamic background with a Wallace rule for short oligos and a salt-adjusted longer-probe formula. The current source implements the Wallace rule for lengths below 14 and a salt-adjusted formula through `ThermoConstants.CalculateSaltAdjustedTm(...)` for longer probes.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The base ranking pass retains only raw-score-positive candidates before optional specificity rescaling | `DesignProbesOptimized(...)` rejects candidates when `score <= 0`, but the suffix-tree overload can later rescale shortlisted scores |
| INV-02 | Probe GC content is computed as a fraction of length | The source uses `gcCount / length` |
| INV-03 | `CheckSpecificity(...)` returns `0` for no hits, `1` for a unique hit, and `1 / hits` otherwise | That mapping is explicit in source |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `targetSequence` | `string` | required | Sequence from which probe candidates are generated | Uppercased before processing |
| `parameters` | `ProbeParameters?` | `Defaults.Microarray` | Application-specific probe design limits | Includes length, Tm, GC, homopolymer, and self-complementarity thresholds |
| `maxProbes` | `int` | `10` | Maximum number of returned probes | Applied after ranking |
| `genomeIndex` | `ISuffixTree` | required for specificity overload | Pre-built suffix tree for genome-wide uniqueness filtering | Used only by the overload with specificity checking |
| `requireUnique` | `bool` | `true` | Whether non-unique probes are excluded when `genomeIndex` is provided | Filters candidates with specificity `< 1.0` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Sequence` | `string` | Probe sequence |
| `Start` / `End` | `int` | Probe coordinates in the source sequence |
| `Tm` | `double` | Melting-temperature estimate |
| `GcContent` | `double` | GC fraction |
| `Score` | `double` | Heuristic quality score |
| `Type` | `ProbeType` | Probe category such as `Standard`, `Tiling`, `Antisense`, `LNA`, or `MolecularBeacon` |
| `Warnings` | `IReadOnlyList<string>` | Quality warnings recorded during evaluation |

### 3.3 Preconditions and Validation

`DesignProbes(...)` returns no probes when the target sequence is null, empty, or shorter than the configured minimum length. All sequences are converted to uppercase before processing. The suffix-tree overload first builds a larger raw-score shortlist than the final requested count, then either filters that shortlist for uniqueness or scales shortlisted scores by the specificity value returned by `CheckSpecificity(...)`; when `requireUnique` is `false`, a shortlisted probe can therefore remain in the output with final score `0` if specificity is `0`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input sequence to uppercase.
2. Precompute GC prefix sums for the target sequence.
3. Enumerate all candidate windows within the configured length range.
4. Evaluate each candidate for GC content, Tm, homopolymers, self-complementarity, secondary structure, repeats, and terminal G/C penalties.
5. Keep only candidates with positive raw scores, sort them by score, and return the top results.
6. In the genome-index overload, apply suffix-tree specificity filtering or post-shortlist score adjustment before yielding the final probes.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Application-specific defaults from `ProbeDesigner.Defaults`:

| Application | Length (bp) | Tm (°C) | GC Content | Notes |
|-------------|-------------|---------|------------|-------|
| Microarray | 50-60 | 75-85 | 0.40-0.60 | Default parameter set |
| FISH | 200-500 | 70-90 | 0.35-0.65 | Higher self-complementarity tolerance |
| Northern Blot | 100-300 | 65-80 | 0.40-0.60 | Intermediate probe sizes |
| qPCR | 20-30 | 68-72 | 0.40-0.60 | Shortest default probes |
| Southern Blot | 150-500 | 65-75 | 0.35-0.65 | Long-probe setting |

The implementation uses prefix sums for O(1) GC lookup and a suffix-tree overload for O(m) uniqueness checking of each candidate probe.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DesignProbes` | `O(n × m)` | `O(k)` | `n` is sequence length, `m` is the scanned length range, and `k` is the number of retained candidates |
| `CheckSpecificity` | `O(m)` | `O(1)` | Uses suffix-tree lookups per probe |
| `DesignTilingProbes` | `O(n)` over fixed-length windows | `O(k)` | Produces overlapping probes for coverage |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProbeDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs), [ThermoConstants.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/ThermoConstants.cs)

- `ProbeDesigner.DesignProbes(string, ProbeParameters?, int)`: Main probe-generation and ranking routine.
- `ProbeDesigner.DesignProbes(string, ISuffixTree, ProbeParameters?, int, bool)`: Uniqueness-aware overload using a suffix tree.
- `ProbeDesigner.DesignTilingProbes(...)`: Generates overlapping tiling probes for coverage.
- `ProbeDesigner.CheckSpecificity(string, ISuffixTree)`: Maps suffix-tree hit counts to a specificity score.

### 5.2 Current Behavior

The implementation evaluates candidates with prefix-sum GC optimization and begins from a raw-score-positive shortlist. Probe sequences are uppercased before evaluation. The suffix-tree overload does not rescore the full candidate universe: it rechecks only a larger raw-score shortlist, currently `maxProbes * 5`, can enforce uniqueness on that shortlist, or scales shortlisted scores by `1 / hitCount`; if `requireUnique` is `false` and specificity is `0`, a shortlisted probe can remain in the final output with score `0`, so final ordering still inherits the initial raw-score pass. `DesignTilingProbes(...)` includes suboptimal probes when needed for coverage and reports coverage, mean Tm, and Tm range. The source also defines probe types `Standard`, `Tiling`, `Antisense`, `LNA`, and `MolecularBeacon`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Application-specific probe-length, GC, and Tm windows.
- Heuristic penalties for GC, Tm, self-complementarity, secondary structure, and repeats.
- Genome-index-based uniqueness checking through a suffix tree.

**Intentionally simplified:**

- Probe quality is represented by a fixed additive penalty score rather than a thermodynamic binding model; **consequence:** scores rank candidates heuristically but are not direct hybridization probabilities.
- Self-complementarity and secondary structure are detected with simple sequence rules; **consequence:** complex folding behavior is not explicitly modeled.
- The current Tm implementation uses the repository's fixed formula helpers rather than a full nearest-neighbor duplex calculation; **consequence:** fine-grained context effects are not captured.
- Genome-index specificity is applied only after an initial raw-score shortlist is formed; **consequence:** uniqueness-aware results are specificity-filtered or specificity-scaled subsets of the top raw-score candidates rather than a full-candidate rerank.

**Not implemented:**

- Database-style alignment or experimentally calibrated hybridization prediction; **users should rely on:** external probe-validation workflows when those are required.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Longer-probe Tm calculation is delegated to `ThermoConstants.CalculateSaltAdjustedTm(...)` | Assumption | Repository behavior follows the shared thermodynamic helper rather than the single explicit longer-probe formula shown in the original document | accepted | Source uses `81.5 + 16.6 * log10([Na+]) + 41 * gcFraction - 600 / length` through the shared helper |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null or empty target sequence | Returns no probes | Explicit early return in source |
| Sequence shorter than `MinLength` | Returns no probes | No valid candidate window exists |
| Candidate with score `<= 0` | Rejected | The evaluator returns `null` for non-positive scores |
| Specificity check with no genome hits | Returns `0` | The probe does not match the indexed genome |

### 6.2 Limitations

The current implementation uses heuristic penalties, simple self-structure detection, and shared Tm helpers instead of a full thermodynamic or database-backed specificity model. It is suitable for fast candidate generation and filtering, but not for high-confidence experimental validation by itself.

## 8. References

1. Wikipedia. "Nucleic acid thermodynamics." https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics
2. Wikipedia. "Hybridization probe." https://en.wikipedia.org/wiki/Hybridization_probe
3. Wikipedia. "Fluorescence in situ hybridization." https://en.wikipedia.org/wiki/Fluorescence_in_situ_hybridization
4. Wikipedia. "DNA microarray." https://en.wikipedia.org/wiki/DNA_microarray
5. SantaLucia, J. (1998). "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." PNAS 95(4):1460-5.
6. Breslauer, K.J. et al. (1986). "Predicting DNA Duplex Stability from the Base Sequence." PNAS 83:3746-3750.
