# Codon Adaptation Index (CAI) Calculation

| Field | Value |
|-------|-------|
| Algorithm Group | Codon Optimization |
| Test Unit ID | CODON-CAI-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

The Codon Adaptation Index (CAI) measures how strongly a coding sequence favors codons that are preferred in a reference organism.[1][3] In this repository, `CodonOptimizer.CalculateCAI` computes CAI from organism-specific codon usage tables using the geometric mean of relative adaptiveness values and a logarithmic accumulation strategy for numeric stability. The implementation normalizes DNA to RNA notation, skips stop codons, and returns `0` when there are no evaluable codons. The result is organism-specific because the codon frequencies come from the supplied `CodonUsageTable`.[1][4][5][6]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Synonymous codon usage bias reflects organism-specific translation preferences and is widely used to study gene expression potential, codon optimization, and sequence adaptation.[1][3] The original repository document records these practical properties:

| Property | Meaning |
|----------|---------|
| Organism specificity | The same coding sequence can have different CAI values in different organisms. |
| Geometric-mean sensitivity | A single rare codon can substantially lower the overall CAI. |
| Range | CAI is bounded by `0` and `1`, with `1` representing exclusively optimal codons. |

The repository ships three predefined reference tables in `CodonOptimizer`.[4][5][6]

| Table | API Symbol | Example preference noted in the original document |
|-------|------------|-----------------------------------------------|
| E. coli K12 | `CodonOptimizer.EColiK12` | Leucine strongly favors `CUG` (`0.50`). |
| S. cerevisiae | `CodonOptimizer.Yeast` | Leucine favors `UUA` and `UUG`. |
| H. sapiens | `CodonOptimizer.Human` | Leucine still favors `CUG` (`0.40`), with weaker bias than E. coli. |

### 2.2 Core Model

For each codon `i` encoding amino acid `a`, the relative adaptiveness is:

```text
w_i = f_i / max(f_j)  for all synonymous codons j of amino acid a
```

where `f_i` is the frequency of codon `i` in the reference table. CAI is then the geometric mean over the `L` non-stop codons in the sequence:

```text
CAI = (product(w_i))^(1 / L)
CAI = exp((1 / L) * sum(ln(w_i)))
```

**Single-codon amino acids.** Sharp & Li (1987) — quoted verbatim by Jansen et al. (2003)[8] —
state that "codon families containing a single codon (e.g. AUG and UGG in the standard genetic
code) should be excluded in computing CAI", because their `w` is always `1` regardless of bias and
including them inflates CAI for Met/Trp-rich genes. The method exposes both conventions through the
optional `excludeSingleCodonAminoAcids` flag: the default (`false`) *includes* Met/Trp with
`w = 1.0` (historical behaviour); `true` *excludes* them so `L` and the product run only over the
multi-codon amino-acid positions, matching the canonical Sharp & Li / Jansen definition.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The provided `CodonUsageTable` is a meaningful reference for the target organism or condition. | CAI becomes a score against the wrong codon-preference landscape. |
| ASM-02 | The input sequence is interpreted in coding-frame triplets. | Trailing bases outside a complete codon are ignored and the score reflects only the retained codons. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 <= CAI <= 1`. | `w_i` values are normalized by the maximum synonymous frequency and the result is a geometric mean over evaluated codons. |
| INV-02 | In the default mode, single-codon amino acids such as Met (`AUG`) and Trp (`UGG`) contribute `w_i = 1`; when `excludeSingleCodonAminoAcids` is `true` they are dropped from `L` and the product entirely. | Their synonymous set contains only the codon being evaluated, so `w ≡ 1`; Sharp & Li (1987)/Jansen (2003) prescribe excluding them.[1][8] |
| INV-03 | Stop codons do not affect the result. | The implementation skips codons whose translated amino acid is `*`. |
| INV-04 | The same sequence can yield different CAI values for different organisms. | Frequencies are taken from the caller-supplied codon usage table. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `codingSequence` | `string` | required | Coding sequence in DNA or RNA notation. | `null` or empty input returns `0`. |
| `table` | `CodonUsageTable` | required | Reference codon usage table used to compute relative adaptiveness. | Must provide codon-frequency data for the desired organism or reference set. |
| `excludeSingleCodonAminoAcids` | `bool` | `false` | When `true`, codons of single-codon amino acids (Met/`AUG`, Trp/`UGG`) are excluded from the geometric mean per Sharp & Li (1987)/Jansen (2003).[1][8] | Optional; default `false` preserves the historical inclusive behaviour. |

### 3.2 Output / Return Value

| Name | Type | Description |
|------|------|-------------|
| `CAI` | `double` | Geometric-mean codon adaptation score for the evaluated non-stop codons. |

### 3.3 Preconditions and Validation

The implementation uppercases the sequence and converts `T` to `U` before splitting into triplets. Incomplete trailing bases are ignored because codons are extracted only when three bases are available. If the sequence is empty, or if no non-stop codons remain after filtering, the method returns `0`. Codons that translate to a non-standard amino acid and therefore lack synonymous-frequency data are skipped.

## 4. Algorithm

### 4.1 High-Level Steps

1. Return `0` for `null` or empty input.
2. Normalize the sequence to uppercase RNA notation.
3. Split the sequence into complete codons.
4. For each codon, translate it to an amino acid and skip stop codons. When `excludeSingleCodonAminoAcids` is `true`, also skip codons of single-codon amino acids (Met/`AUG`, Trp/`UGG`).
5. Compute relative adaptiveness `w_i` against the maximum synonymous frequency in the supplied table.
6. Accumulate `ln(w_i)` and count the evaluated codons.
7. Return `exp(logSum / count)` when at least one codon was evaluated; otherwise return `0`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Relative adaptiveness is calculated against the maximum synonymous frequency for the codon's amino acid. When the amino acid has no frequency data in the table, the codon is skipped. When the amino acid is represented in the table but the specific codon has frequency `0`, the current implementation clamps `w_i` to `1e-6` rather than allowing `0`, as documented in the test specification.[7]

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateCAI` | `O(n)` | `O(1)` | `n` is the sequence length; reference tables are constant-size. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [CodonOptimizer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs)

- `CodonOptimizer.CalculateCAI(string, CodonUsageTable, bool excludeSingleCodonAminoAcids = false)`
- `CodonOptimizer.CalculateRelativeAdaptiveness(...)` (private helper)
- `CodonOptimizer.SingleCodonAminoAcids` (private set derived from the standard genetic code: amino acids with exactly one codon, i.e. Met/`M` and Trp/`W`)

### 5.2 Current Behavior

The method returns `0` for `null`, empty, or no-evaluable-codon input. It converts DNA input to RNA notation, ignores incomplete trailing bases, and excludes stop codons from the codon count `L`. Unknown or non-standard codons translate to `X`; because no synonymous-frequency set exists for `X`, those codons are skipped rather than assigned a frequency. Zero-frequency codons inside an otherwise populated synonymous set are clamped to `1e-6` before taking the logarithm.[7]

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Relative adaptiveness uses `w_i = f_i / max(f_j)` across synonymous codons.[1]
- CAI is computed as a geometric mean using the equivalent logarithmic form `exp((1/L) * sum(ln(w_i)))`.[1]
- The single-codon-amino-acid exclusion of Sharp & Li (1987)/Jansen (2003) is available verbatim via `excludeSingleCodonAminoAcids: true` (Met/`AUG`, Trp/`UGG` dropped from `L` and the product).[1][8]

**Intentionally simplified:**

- The implementation relies entirely on the provided codon-usage table and does not infer organism context; **consequence:** interpretation is only as good as the caller-supplied reference table.
- Codons without recognizable synonymous-frequency data are skipped instead of causing an error; **consequence:** malformed or non-standard input can reduce the effective codon count without an exception.
- The Met/Trp exclusion is **opt-in**, not the default; **consequence:** by default CAI includes Met/Trp with `w = 1.0`, which can bias CAI upward for Met/Trp-rich genes relative to strict Sharp & Li. Pass `excludeSingleCodonAminoAcids: true` for the canonical convention.

**Not implemented:**

- Confidence intervals, bootstrap uncertainty, or statistical significance for CAI scores; **users should rely on:** no current alternative.
- Automatic validation that the input sequence is a biologically valid CDS; **users should rely on:** caller-side validation.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Zero-frequency codons are clamped to `1e-6` when other synonymous codons exist in the table. | Deviation | CAI remains positive instead of becoming exactly `0` for those codons. | accepted | Documented in [CODON-CAI-001.md](../../../tests/TestSpecs/CODON-CAI-001.md). |
| 2 | Single-codon amino acids (Met/Trp) are **included** by default; the canonical Sharp & Li exclusion is opt-in via `excludeSingleCodonAminoAcids: true`. | Convention (opt-in) | Default biases CAI slightly upward for Met/Trp-rich genes vs strict Sharp & Li; the opt-in mode is exact. | accepted | Sharp & Li (1987)/Jansen (2003).[1][8] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns `0`. | Explicit early return. |
| Sequence containing only stop codons | Returns `0`. | Stop codons are excluded, leaving no evaluated codons. |
| DNA input | Treated as RNA after `T -> U` normalization. | The method normalizes notation before splitting. |
| Lowercase input | Handled identically to uppercase input. | The sequence is uppercased first. |
| Incomplete trailing codon | Ignored. | Only complete triplets are split into codons. |
| Sequence of only Met/Trp, `excludeSingleCodonAminoAcids: true` | Returns `0`. | All codons are excluded, leaving no evaluated codons (`L = 0`). |
| Sequence with no Met/Trp, either mode | Identical result in both modes. | The exclusion flag only affects single-codon amino-acid positions. |

### 6.2 Limitations

CAI in this repository is a table-driven codon-bias score. It does not model tRNA abundance explicitly, does not validate full biological correctness of the coding sequence, and does not account for context effects such as codon pairs or mRNA structure. Interpretation remains specific to the chosen codon-usage table. The single-codon amino-acid exclusion of Sharp & Li (1987) is supported but opt-in (`excludeSingleCodonAminoAcids: true`); the default keeps Met/Trp for backward compatibility.

## 7. Examples and Related Material

### 7.1 Worked Example

The original document included these hand-worked E. coli examples:

```text
Optimal sequence: AUGCUGCCGACC
Codons: AUG(1.0), CUG(1.0), CCG(1.0), ACC(1.0)
CAI = 1.0

Suboptimal sequence: AUGCUACCAACU
Codons: AUG(1.0), CUA(0.08), CCA(0.358), ACU(0.364)
CAI ≈ 0.31
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [CodonOptimizer_CAI_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_CAI_Tests.cs) — covers `INV-01`, `INV-02` (both default and exclusion modes), `INV-03`, `INV-04`
- Test specification: [CODON-CAI-001.md](../../../tests/TestSpecs/CODON-CAI-001.md)
- Related algorithms: [Codon_Usage_Analysis.md](Codon_Usage_Analysis.md), [Sequence_Optimization.md](Sequence_Optimization.md)

## 8. References

1. Sharp PM, Li WH. 1987. The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research. N/A
2. Wikipedia contributors. 2026. Codon Adaptation Index. Wikipedia. N/A
3. Plotkin JB, Kudla G. 2011. Synonymous but not the same. Nature Reviews Genetics. N/A
4. Kazusa Codon Usage Database. E. coli K-12 substr. W3110, species 316407. https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=316407
5. Kazusa Codon Usage Database. Saccharomyces cerevisiae, species 4932. https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932
6. Kazusa Codon Usage Database. Homo sapiens, species 9606. https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
7. Test specification: [CODON-CAI-001.md](../../../tests/TestSpecs/CODON-CAI-001.md)
8. Jansen R, Bauer P, Stadler PF. 2003. An Improved Implementation of the Codon Adaptation Index. Nucleic Acids Research. Retrieved 2026-06-24 from https://pmc.ncbi.nlm.nih.gov/articles/PMC2684136/
