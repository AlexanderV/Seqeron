# Codon Adaptation Index (CAI) Calculation

| Field | Value |
|-------|-------|
| Algorithm Group | Codon Optimization |
| Test Unit ID | CODON-CAI-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-09-28 (review campaign B02: F12–F14) |

## 1. Overview

The Codon Adaptation Index (CAI) measures how strongly a coding sequence favors codons that are preferred in a reference organism.[1][3] In this repository there is one CAI implementation, the canonical core `CodonUsageAnalyzer.CalculateCai` (reference RSCU/w table, genetic-code aware); `CodonOptimizer.CalculateCAI` (organism codon-usage frequency table) delegates to it. Both follow the CodonW reference implementation (`cai_out`, Peden 1999): geometric mean of relative adaptiveness over the scored codons, stop codons and single-codon families (Met/Trp) excluded, `w < 0.0001` scored as `0.01`, and `0` when no codon is scored. The result is organism-specific because the codon frequencies come from the supplied `CodonUsageTable`.[1][4][5][6]

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

**Single-codon amino acids and stops.** Sharp & Li (1987) — quoted by Xia (2007)[8] — state that
"codon families containing a single codon (e.g. AUG and UGG in the standard genetic code) should be
excluded in computing CAI", because their `w` is always `1` regardless of bias. CodonW `cai_out`
("Non-synonymous codons and termination codons (genetic code dependent) are excluded"), seqinr `cai`
and Biopython `CodonAdaptationIndex.calculate` (ATG/TGG skipped) do the same. This is the default
(`excludeSingleCodonAminoAcids: true`); `false` scores Met/Trp with `w = 1`, as EMBOSS
`ajCodCalcCaiSeq` does (EMBOSS also counts stops; that part is not followed).

**Codons absent from the reference.** Sharp & Li assign a count of 0.5 to codons absent from the
reference set when *building* w (Biopython 1.88: "Following the description in the original paper, we
use a value of 0.5 for codons that do not appear in the reference sequences"; CodonW
`highest_x`/"adjust its frequency to 0.5"). When *scoring* a gene against a given w table, CodonW
replaces any `w < 0.0001` by `0.01` ("To prevent a codon having a relative adaptiveness value of zero,
which could result in a CAI of zero"; seqinr `zero.threshold = 0.0001, zero.to = 0.01`, Bulmer 1988).
The implementation applies the latter, since it receives a w/frequency table, not counts.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The provided `CodonUsageTable` is a meaningful reference for the target organism or condition. | CAI becomes a score against the wrong codon-preference landscape. |
| ASM-02 | The input sequence is interpreted in coding-frame triplets. | Trailing bases outside a complete codon are ignored and the score reflects only the retained codons. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 <= CAI <= 1`. | `w_i` values are normalized by the maximum synonymous frequency and the result is a geometric mean over evaluated codons. |
| INV-02 | By default single-codon amino acids (Met `AUG`, Trp `UGG` in table 1) are dropped from `L` and the product; with `excludeSingleCodonAminoAcids: false` they contribute `w = 1`. | Sharp & Li (1987)/Xia (2007); CodonW.[1][8][9] |
| INV-03 | Stop codons do not affect the result. | Stop codons (per genetic code) are never scored. |
| INV-05 | Every scored codon has `w ≥ 0.0001`, so CAI > 0 whenever a codon is scored. | CodonW `w < 0.0001 → 0.01`.[9] |
| INV-04 | The same sequence can yield different CAI values for different organisms. | Frequencies are taken from the caller-supplied codon usage table. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `codingSequence` | `string` | required | Coding sequence in DNA or RNA notation. | `null` or empty input returns `0`. |
| `table` | `CodonUsageTable` | required | Reference codon usage table used to compute relative adaptiveness. | Must provide codon-frequency data for the desired organism or reference set. |
| `excludeSingleCodonAminoAcids` | `bool` | `true` | `true`: Met/`AUG`, Trp/`UGG` excluded (Sharp & Li 1987; Xia 2007; CodonW).[1][8][9] `false`: scored with `w = 1` (EMBOSS-style). | Optional. Changed from `false` in review 2026-09 (F12). |

### 3.2 Output / Return Value

| Name | Type | Description |
|------|------|-------------|
| `CAI` | `double` | Geometric-mean codon adaptation score for the evaluated non-stop codons. |

### 3.3 Preconditions and Validation

Input is case-insensitive DNA or RNA (`U` read as `T`). Triplets containing any other symbol (IUPAC
codes, gaps) are skipped without shifting the frame (CodonW `ident_codon` gives them codon 0, not
counted); a trailing partial codon is ignored. `null`/empty input or no scored codon returns `0`.
Table keys may be RNA or DNA spelling. A negative or non-finite frequency/RSCU value throws
`ArgumentOutOfRangeException` (CodonW exits on w outside [0, 1]).

## 4. Algorithm

### 4.1 High-Level Steps

1. Return `0` for `null` or empty input.
2. Build `w` per synonymous family of the genetic code (Standard for `CodonOptimizer`): skip the stop family and (by default) single-codon families; `w = value / family max`; a family whose maximum is 0 has no reference data and is not scored; `w < 0.0001 → 0.01`.
3. Normalize the sequence (upper case, `U→T`) and read frame-0 triplets.
4. For each triplet that has a `w`, accumulate `ln w` and count it (`L`).
5. Return `exp(Σ ln w / L)` when `L > 0`, otherwise `0`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

`CodonOptimizer.CalculateCAI` converts the table keys to DNA spelling and calls the internal overload
of `CodonUsageAnalyzer.CalculateCai` with `GeneticCode.Standard`. Codons that are context-dependent
stops in NCBI tables 27/28/31 are scored in the family of the amino acid they encode (as for RSCU).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateCAI` | `O(n)` | `O(1)` | `n` is the sequence length; reference tables are constant-size. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

- `CodonOptimizer.CalculateCAI(string, CodonUsageTable, bool excludeSingleCodonAminoAcids = true)` — [CodonOptimizer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs)
- Canonical core: `CodonUsageAnalyzer.CalculateCai(string | DnaSequence, IReadOnlyDictionary<string,double>, GeneticCode)` and the Standard-code overloads `CalculateCai(string | DnaSequence, Dictionary<string,double>)` — [CodonUsageAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs)
- MCP: `cai_from_organism_table` → `CodonOptimizer.CalculateCAI`; `codon_adaptation_index` → `CodonUsageAnalyzer.CalculateCai`.

### 5.2 Current Behavior

See §3.3/§4. Cross-checked (review 2026-09) against a verbatim Python port of CodonW 1.4.4 `cai_out`
(itself matched to the compiled CodonW binary: 6 656 cases × 8 genetic codes × 4 w tables incl. zeros,
0 mismatches at CodonW's 3 dp) — C# vs port: 0 mismatches (max |Δ| 5.6e-16) over 832 inputs × 27 NCBI
tables and for `CodonOptimizer`; vs Biopython 1.88 `CodonAdaptationIndex` (400 genes, 40 indices built
with the 0.5 pseudo-count): max |Δ| 0.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):** `w = X/X_max` per synonymous family, geometric
mean via log-sum, stop and single-codon-family exclusion (genetic-code dependent), CodonW `0.01`
substitution.

**Conventions (documented):**
- A family with no reference data is not scored (no w is defined).
- Values are rescaled by the family maximum, so RSCU and w tables are both accepted.
- `excludeSingleCodonAminoAcids: false` (EMBOSS-style) is a non-Sharp & Li opt-in.
- Tables from `CreateCodonTableFromSequence` hold frequencies, not counts, so the Sharp & Li 0.5
  pseudo-count for absent codons (used when building w from counts) is not applied there; the scoring
  rule `w < 0.0001 → 0.01` applies instead.

**Not implemented:** confidence intervals / expected CAI (Xia 2007 eCAI); CDS validity checks.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Zero-frequency codons clamped to `1e-6` | Deviation | each absent codon lowered CAI by an extra factor (10⁻⁴)^(1/L) vs CodonW | **fixed 2026-09 (F13)** | now `w < 0.0001 → 0.01` (CodonW) |
| 2 | Met/Trp included by default | Deviation | CAI inflated for Met/Trp-rich genes | **fixed 2026-09 (F12)** | default now excludes; `false` = opt-in |
| 3 | Two CAI implementations disagreed | Duplication | MCP tools gave different CAI for the same data | **fixed 2026-09 (F14)** | CodonOptimizer delegates to CodonUsageAnalyzer |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns `0`. | Explicit early return. |
| Sequence containing only stop codons | Returns `0`. | Stop codons are excluded, leaving no evaluated codons. |
| DNA input | Treated as RNA after `T -> U` normalization. | The method normalizes notation before splitting. |
| Lowercase input | Handled identically to uppercase input. | The sequence is uppercased first. |
| Incomplete trailing codon | Ignored. | Only complete triplets are split into codons. |
| Sequence of only Met/Trp (default) | Returns `0`. | All codons are excluded (`L = 0`); CodonW `ATGTGG` → 0.000. |
| Codon absent from the table (synonym present) | Scored with `w = 0.01`. | CodonW `cai_out`. |
| Ambiguous triplet (`NNN`, `CUR`) | Skipped, frame kept. | CodonW `ident_codon`. |
| Sequence with no Met/Trp, either mode | Identical result in both modes. | The exclusion flag only affects single-codon amino-acid positions. |

### 6.2 Limitations

CAI in this repository is a table-driven codon-bias score. It does not model tRNA abundance explicitly, does not validate full biological correctness of the coding sequence, and does not account for context effects such as codon pairs or mRNA structure. Interpretation remains specific to the chosen codon-usage table. The Sharp & Li single-codon exclusion is the default.

## 7. Examples and Related Material

### 7.1 Worked Example

The original document included these hand-worked E. coli examples:

```text
Optimal sequence: AUGCUGCCGACC
Codons: AUG (not scored), CUG(1.0), CCG(1.0), ACC(1.0)
CAI = 1.0

Suboptimal sequence: AUGCUACCAACU
Codons: AUG (not scored), CUA(0.04/0.50), CCA(0.19/0.53), ACU(0.16/0.44)
CAI = 0.2184799938153881   (was ≈0.31 when AUG was scored with w = 1)
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
8. Xia X. 2007. An Improved Implementation of Codon Adaptation Index. Evolutionary Bioinformatics 3:53-58. PMC2684136 (formerly mis-cited here as Jansen et al. 2003). Statement seen via search snippet of journals.sagepub.com/doi/full/10.1177/117693430700300028 (2026-09-28).
9. Peden JF. 1999. CodonW 1.4.4, `codon_us.c` `cai_out` / `codonW.h` `cai[]` (original tarball, compiled and run 2026-09); seqinr `R/cai.R`, `man/cai.Rd` (raw.githubusercontent.com/cran/seqinr); Bulmer M. 1988. J. Evol. Biol. 1:15-26.
10. Biopython 1.88 `Bio/SeqUtils/__init__.py` `CodonAdaptationIndex`; Biopython 1.79 `Bio/SeqUtils/CodonUsage.py`, `CodonUsageIndices.py` (`SharpEcoliIndex`).
11. EMBOSS `ajcod.c` `ajCodCalcCaiSeq`, `codGetWstat`.
