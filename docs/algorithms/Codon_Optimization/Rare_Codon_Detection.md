# Rare Codon Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Codon Optimization |
| Test Unit ID | CODON-RARE-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Rare codon detection identifies codons whose reference frequency in a target organism falls below a chosen threshold.[1][2][3] In this repository, `CodonOptimizer.FindRareCodons` streams the positions, codons, amino-acid translations, and table frequencies for codons whose frequency is strictly less than the threshold. The method is designed for sequence inspection and optimization workflows, especially when expressing genes in a heterologous host. It is deterministic, linear in the number of complete codons, and relies entirely on the supplied `CodonUsageTable`.[5]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Rare codons can slow translation because they are associated with lower-abundance tRNAs, and they are often discussed in the context of expression optimization and co-translational folding.[2][3][4] The original document recorded these common rare codons for E. coli K12 from Kazusa-backed frequencies:

| Codon | Amino Acid | Frequency | tRNA Gene |
|-------|------------|-----------|-----------|
| `AGA` | Arginine | `0.04` | `argU` |
| `AGG` | Arginine | `0.02` | `argW` |
| `CGA` | Arginine | `0.06` | `argW` |
| `CUA` | Leucine | `0.04` | `leuV` |

### 2.2 Core Model

For each codon in the normalized coding sequence, the repository looks up its frequency in the supplied codon-usage table and emits a result when:

```text
frequency(codon) < threshold
```

Each reported item is a tuple of the nucleotide position, codon, translated amino acid, and frequency.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The provided codon-usage table reflects the organism or host of interest. | Reported rare codons may not be biologically rare in the actual target system. |
| ASM-02 | The selected threshold matches the user's screening goal. | Results can be overly conservative or overly broad depending on the cutoff. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported position is a multiple of `3`. | Positions are emitted as `codonIndex * 3`. |
| INV-02 | Every reported frequency is strictly less than `threshold`. | The implementation uses `<`, not `<=`. |
| INV-03 | Every reported codon has length `3`. | Only complete triplets are split into codons. |
| INV-04 | Re-running the method with the same inputs yields the same output. | The method is a deterministic single pass over normalized codons. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `codingSequence` | `string` | required | DNA or RNA coding sequence to inspect. | Empty or `null` input yields no results. |
| `table` | `CodonUsageTable` | required | Reference codon-usage table used for frequency lookups. | Unknown codons default to frequency `0`. |
| `threshold` | `double` | `0.15` | Frequency cutoff for reporting a codon as rare. | Comparison is strict: only frequencies `< threshold` are reported. |

### 3.2 Output / Return Value

| Name | Type | Description |
|------|------|-------------|
| `Position` | `int` | Nucleotide index of the codon start. |
| `Codon` | `string` | Normalized RNA codon. |
| `AminoAcid` | `string` | Single-letter amino-acid translation, or `X` for an unknown codon. |
| `Frequency` | `double` | Frequency retrieved from the supplied codon-usage table. |

### 3.3 Preconditions and Validation

The sequence is uppercased and normalized from DNA to RNA. Codons are extracted in complete triplets only; trailing bases are ignored. Frequencies are looked up with `GetValueOrDefault`, so a codon absent from the table is assigned `0` and therefore reported whenever `threshold > 0`.[5]

## 4. Algorithm

### 4.1 High-Level Steps

1. Return an empty sequence for `null` or empty input.
2. Normalize the sequence to uppercase RNA notation.
3. Split the sequence into complete codons.
4. For each codon, read its table frequency.
5. If the frequency is strictly less than `threshold`, translate the codon and yield `(position, codon, aminoAcid, frequency)`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The original document recorded the following threshold-selection guidance:

| Threshold | Description | Use Case |
|-----------|-------------|----------|
| `0.10` | Very rare | Critical optimization |
| `0.15` | Default | Standard analysis |
| `0.20` | Moderately rare | Broad screening |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindRareCodons` | `O(n)` | `O(1)` auxiliary | Streams results while scanning complete codons. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [CodonOptimizer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs)

- `CodonOptimizer.FindRareCodons(string, CodonUsageTable, double)`

### 5.2 Current Behavior

The method converts `T` to `U`, ignores trailing incomplete codons, and reports nucleotide positions as `i * 3`. Frequency lookup uses `GetValueOrDefault`, so unknown codons are treated as frequency `0`. The threshold comparison is strict `<`, which means a codon exactly at the threshold is not reported.[5]

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Codons are flagged as rare when their reference frequency is below the chosen threshold.[2][5]
- The result includes position, codon, amino acid, and frequency, matching the repository test specification.[5]

**Intentionally simplified:**

- Detection is table-threshold-based and does not model ribosome dynamics or local codon context; **consequence:** the output is a screening list rather than a direct translation-efficiency prediction.
- Unknown codons default to frequency `0`; **consequence:** malformed codons are surfaced as rare rather than rejected.

**Not implemented:**

- Rare-codon cluster scoring or context-aware pause modeling; **users should rely on:** no current alternative.
- Automatic integration with codon optimization inside this method; **users should rely on:** [Sequence_Optimization.md](Sequence_Optimization.md).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns no results. | Explicit early exit. |
| Default-threshold call | Uses `0.15`. | Default parameter value in the public API. |
| Frequency exactly equal to threshold | Not reported. | The comparison is `<`, not `<=`. |
| Unknown codon | Reported with frequency `0` and amino acid `X` when `threshold > 0`. | Frequency lookup defaults to `0` and translation defaults to `X`. |
| Incomplete trailing codon | Ignored. | Codon splitting only emits complete triplets. |

### 6.2 Limitations

Rare-codon detection here is a simple table lookup over codon triplets. It does not consider rare-codon clustering, local sequence context, initiation effects, or organism-specific regulation beyond the supplied frequency table.

## 7. Examples and Related Material

### 7.1 Worked Example

```csharp
string sequence = "AUGAGAAGGCGA";
var rareList = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10);
// Returns: (3, "AGA", "R", 0.04), (6, "AGG", "R", 0.02), (9, "CGA", "R", 0.06)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [CodonOptimizer_FindRareCodons_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/CodonOptimizer_FindRareCodons_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Test specification: [CODON-RARE-001.md](../../../tests/TestSpecs/CODON-RARE-001.md)
- Related algorithms: [Sequence_Optimization.md](Sequence_Optimization.md), [CAI_Calculation.md](CAI_Calculation.md)

## 8. References

1. Sharp PM, Li WH. 1987. The codon adaptation index-a measure of directional synonymous codon usage bias. Nucleic Acids Research. N/A
2. Shu P, Dai H, Gao W, Goldman E. 2006. Inhibition of translation by consecutive rare leucine codons in E. coli. Gene Expression. N/A
3. Plotkin JB, Kudla G. 2011. Synonymous but not the same: causes and consequences of codon bias. Nature Reviews Genetics. N/A
4. Kane JF. 1995. Effects of rare codon clusters on high-level expression of heterologous proteins in E. coli. Current Opinion in Biotechnology. N/A
5. Test specification: [CODON-RARE-001.md](../../../tests/TestSpecs/CODON-RARE-001.md)
