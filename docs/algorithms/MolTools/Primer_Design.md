# Primer Pair Design

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PRIMER-DESIGN-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Primer pair design selects forward and reverse oligonucleotides that can amplify a target DNA region by PCR. In this repository, primer design combines candidate enumeration in flanking regions with per-primer quality evaluation and a simple pair-compatibility check based on melting-temperature agreement and primer-dimer avoidance. The current implementation is a heuristic selector rather than a full Primer3-style combinatorial optimizer.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

PCR primer design balances primer length, GC content, melting temperature, repetitive sequence content, and 3' end properties so that both primers bind specifically and amplify the desired product. The original document cites Primer3 and Addgene for standard design ranges and also notes that 3' terminal stability affects extension efficiency. Sources: Primer3 Manual, Addgene primer-design guidance, Wikipedia (Primer (molecular biology)), SantaLucia (1998).

### 2.2 Core Model

The repository designs primers in four stages: search-region definition, candidate generation, greedy best-candidate selection, and pair compatibility checking. Candidate scoring uses:

$$
score = 100 - 2\lvert length - optimalLength \rvert - 2\lvert T_m - optimalT_m \rvert - 0.5\lvert GC\% - 50 \rvert - 5 \times homopolymerLength + bonus_{GC clamp}
$$

where the GC-clamp bonus is `+5` when the final base is `G` or `C` in the current source. Pair compatibility requires a Tm difference of at most `5°C` and no primer-dimer signal.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `DesignPrimers(...)` returns `IsValid = false` when either side has no valid candidates | The source returns an invalid `PrimerPairResult` when either best candidate is missing |
| INV-02 | Pair validity requires both `|Tm_f - Tm_r| <= 5` and `!HasPrimerDimer(...)` | That conjunction is explicit in source |
| INV-03 | `ProductSize = reverse.Position + reverse.Sequence.Length - forward.Position` | The source computes product size directly from the chosen candidates |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `template` | `DnaSequence` | required | Template DNA sequence |
| `targetStart` | `int` | required | Start of the target region | Must satisfy `targetStart >= 0` |
| `targetEnd` | `int` | required | End of the target region | Must satisfy `targetEnd < template.Length` and `targetStart < targetEnd` |
| `parameters` | `PrimerParameters?` | `PrimerDesigner.DefaultParameters` | Primer design thresholds | Defaults are `18-25` bp length, `40-60%` GC, `57-63°C` Tm, `OptimalLength = 20`, `OptimalTm = 60`, `MaxHomopolymer = 4`, `MaxDinucleotideRepeats = 4`, `Avoid3PrimeGC = false`, and `Check3PrimeStability = true` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Forward` | `PrimerCandidate?` | Selected forward primer or `null` |
| `Reverse` | `PrimerCandidate?` | Selected reverse primer or `null` |
| `IsValid` | `bool` | Pair validity flag |
| `Message` | `string` | Result explanation |
| `ProductSize` | `int` | Predicted amplicon size |

### 3.3 Preconditions and Validation

`DesignPrimers(...)` throws `ArgumentException` when the requested target region is invalid. Forward candidates are searched up to 200 bp upstream of `targetStart`; reverse candidates are searched up to 200 bp downstream of `targetEnd`. Reverse-primer candidates are reverse-complemented before evaluation so that they are scored in primer orientation.

## 4. Algorithm

### 4.1 High-Level Steps

1. Define a forward search region up to 200 bp upstream of the target start.
2. Define a reverse search region up to 200 bp downstream of the target end.
3. Enumerate all candidate primers within the configured length range.
4. Evaluate each candidate for GC content, Tm, homopolymers, dinucleotide repeats, hairpin potential, and 3' stability.
5. Select the highest-scoring forward and reverse candidates independently.
6. Accept the pair only if the Tm difference is at most `5°C` and no primer-dimer is detected.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Parameter ranges documented in the original file and current source:

| Parameter | Min | Optimal | Max | Notes |
|-----------|-----|---------|-----|-------|
| Length (bp) | 18 | 20 | 25 | Matches current source defaults |
| GC Content (%) | 40 | 50 | 60 | Matches current source defaults |
| Melting Temp (°C) | 55 | 60 | 65 | Original document summary; current source defaults narrow this to `57-63` |
| Homopolymer Run | N/A | N/A | 4 | Current source default |
| Dinucleotide Repeats | N/A | N/A | 4 | Current source default |
| `Avoid3PrimeGC` | N/A | N/A | `false` | When enabled, the current check requires at least one `G`/`C` in the last two bases |
| `Check3PrimeStability` | N/A | N/A | `true` | Gates whether `EvaluatePrimer(...)` records the `ΔG < -9` issue |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DesignPrimers` | `O(n²)` | `O(k)` | Candidate enumeration over positions and lengths |
| `EvaluatePrimer` | `O(n)` | `O(1)` | Per-primer scan and helper calculations |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PrimerDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs)

- `PrimerDesigner.DesignPrimers(DnaSequence, int, int, PrimerParameters?)`: Designs and validates a primer pair around a target region.
- `PrimerDesigner.EvaluatePrimer(string, int, bool, PrimerParameters?)`: Scores a single primer candidate.
- `PrimerDesigner.CalculatePrimerScore(...)`: Applies the heuristic candidate score.

### 5.2 Current Behavior

Forward primers are taken directly from the template in the forward orientation. Reverse primers are extracted from the downstream region and reverse-complemented before evaluation so they represent the sequence that binds the reverse strand. The source uses `57-63°C` as the default acceptable Tm range, `40-60%` GC, and `18-25` bp length. The evaluator also computes 3' stability and can emit a threshold-based issue for it when `Check3PrimeStability` is enabled. `Avoid3PrimeGC` is disabled by default; when callers enable it, the current implementation applies a GC-clamp-style check on the last two bases by flagging primers whose final dinucleotide contains no `G` or `C`. Pair selection itself is driven by the aggregate heuristic score plus the later compatibility checks. Pair selection is greedy: the highest-scoring forward and highest-scoring reverse candidates are chosen independently and then checked for compatibility.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Candidate screening by length, GC content, melting temperature, repetitive sequence content, and hairpin risk, with 3' stability calculated and exposed as an additional heuristic diagnostic.
- A pair-compatibility rule based on Tm agreement and primer-dimer avoidance.
- Reverse-primer evaluation in reverse-complement orientation.

**Intentionally simplified:**

- Pair selection is greedy rather than a full pairwise optimization; **consequence:** the highest-scoring individual primers are not guaranteed to be the globally best pair.
- The scoring model is a fixed additive heuristic; **consequence:** Primer3-style pair penalties and more detailed thermodynamic interactions are not represented.
- 3' stability is reported through a simple issue flag rather than a richer thermodynamic pair model; **consequence:** callers can inspect the value, but final ranking is still dominated by the additive heuristic score and compatibility checks.

**Not implemented:**

- Full Primer3-style combinatorial primer-pair optimization and richer laboratory constraints; **users should rely on:** external primer-design tools when those features are required.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Current source defaults use `57-63°C` rather than the broader `55-65°C` summary in the original document | Deviation | Default filtering is slightly narrower than the narrative summary | accepted | Confirmed from `PrimerDesigner.DefaultParameters` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Invalid target region | Throws `ArgumentException` | Explicit source guard |
| No valid forward or reverse candidates | Returns an invalid `PrimerPairResult` with null candidates | Explicit fallback in source |
| Large Tm difference | Returns `IsValid = false` | Pair compatibility requires `<= 5°C` |
| Primer-dimer detected | Returns `IsValid = false` | Pair compatibility requires no dimer signal |

### 6.2 Limitations

The current design surface is a fast heuristic filter. It does not implement exhaustive pair optimization, full Primer3 thermodynamics, or genome-wide specificity analysis, and it inherits the simplified melting-temperature and structure screens defined elsewhere in the repository.

## 7. Examples and Related Material

### 7.2 Applications and Use Cases (Optional)

Related material called out in the original document:

- `PRIMER-TM-001`: Melting temperature calculation (prerequisite).
- `PRIMER-STRUCT-001`: Hairpin and dimer detection (used in evaluation).

## 8. References

1. [Primer (molecular biology)](https://en.wikipedia.org/wiki/Primer_(molecular_biology)) - Standard primer design criteria.
2. [How to Design a Primer](https://www.addgene.org/protocols/primer-design/) - Addgene protocol guidance.
3. [primer3.org/manual.html](https://primer3.org/manual.html) - Primer3 manual.
4. SantaLucia JR (1998). "A unified view of polymer, dumbbell and oligonucleotide DNA nearest-neighbor thermodynamics", PNAS 95:1460-65.
