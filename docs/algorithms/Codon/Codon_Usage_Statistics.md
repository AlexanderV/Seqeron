# Codon Usage Statistics

| Field | Value |
|-------|-------|
| Algorithm Group | Codon |
| Test Unit ID | CODON-STATS-001 |
| Related Projects | Seqeron.Genomics.MolTools |
| Implementation Status | Production |
| Last Reviewed | 2026-09-28 |

## 1. Overview

Codon Usage Statistics aggregates the standard descriptors of synonymous codon bias for a coding DNA sequence: the codon-count table, Relative Synonymous Codon Usage (RSCU), the Effective Number of Codons (ENC), the G+C content at each codon position (GC1/GC2/GC3), the G+C content at synonymous third positions (GC3s), and the total number of in-frame codons. It also provides the Codon Adaptation Index (CAI) of a sequence against a reference codon-usage set, and two reference tables (E. coli, human). The statistics are exact, specification-driven summaries used to study codon bias, infer expression level, and guide codon optimization [1][2].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The genetic code is degenerate: most amino acids are encoded by several synonymous codons. Genomes and individual genes use synonyms unequally; the degree and direction of this bias correlates with gene expression and base composition [1][8]. The descriptors below quantify that bias from a single coding sequence.

### 2.2 Core Model

- **Codon counts:** the number of occurrences of each in-frame ACGT codon.
- **RSCU** (Sharp, Tuohy & Mosurski 1986): for codon j of an amino acid with n synonymous codons, `RSCU_j = x_j / ((1/n) Σ_k x_k) = n·x_j / Σ_k x_k` [8].
- **ENC** (Wright 1990): the effective number of codons, 20 ≤ ENC ≤ 61 [9] (see the dedicated [Effective_Number_of_Codons](Effective_Number_of_Codons.md) doc).
- **GC1/GC2/GC3:** the fraction (here ×100) of in-frame codons with G or C at codon position 1, 2, 3 respectively ("1st/2nd/3rd letter GC"), over all valid codons including termination codons, as EMBOSS `cusp` (`ajcod.c`) [5]. CodonW's GC/GC1-3 exclude stop codons, so the two differ only for genes containing stops.
- **GC3s:** "the frequency of G or C nucleotides present at the third position of synonymous codons (i.e. excluding Met, Trp and termination codons)" — it counts only codons whose amino acid has more than one codon **in the selected genetic code** (CodonW 1.4.4 `gc_out`/`how_synon`, `-code`) [2].
- **CAI** (Sharp & Li 1987): relative adaptiveness `w_i = f_i / max(f_j)` over a codon's synonymous family in a reference set; `CAI = (∏_{i=1}^{L} w_i)^{1/L} = exp[(1/L) Σ ln w_i]`. Non-synonymous (single-codon Met, Trp) and termination codons are excluded [1][3][4].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The reference set passed to CAI represents the target codon preference (ideally highly expressed genes). | CAI no longer tracks expression adaptation; it just measures similarity to the supplied table. |
| ASM-02 | The sequence is a real coding sequence in frame 1 (multiple-of-three, ATG…stop). | Out-of-frame nucleotides shift codon assignment and distort all statistics. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 ≤ CAI ≤ 1` | geometric mean of `w_i ∈ [0,1]` [1] |
| INV-02 | CAI of an all-optimal sequence = 1 | every codon equals its family's w-max (w=1) [1] |
| INV-03 | GC3s excludes stop codons and single-codon amino acids of the genetic code (Met ATG, Trp TGG in table 1) | definition in Peden §1.8.2.1.3 [2]; CodonW `gc_out` |
| INV-04 | `0 ≤ GC1, GC2, GC3, GC3s ≤ 100` | each is a count/positions ratio ×100 |
| INV-05 | `TotalCodons` = number of valid ACGT codons in frame | codons with non-ACGT characters are skipped [5] |
| INV-06 | `OverallGc = (GC1+GC2+GC3)/3` | record-derived property |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | `string` / `DnaSequence` | required | coding DNA sequence | case-insensitive; non-ACGT codons skipped; read in frame 1, step 3 |
| referenceRscu | `Dictionary<string,double>` | required (CAI) | reference codon weights (RSCU or w) | keyed by upper-case DNA codon |
| code | `GeneticCode` | `GeneticCode.Standard` | genetic code for RSCU, ENC, GC3s (and CAI) | non-null |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `CodonCounts` | `IReadOnlyDictionary<string,int>` | per-codon occurrence counts |
| `Rscu` | `IReadOnlyDictionary<string,double>` | RSCU per codon |
| `Enc` | `double` | Effective Number of Codons (20–61) |
| `TotalCodons` | `int` | number of valid in-frame codons |
| `Gc1`, `Gc2`, `Gc3` | `double` | % G/C at codon positions 1/2/3 |
| `Gc3s` | `double` | % G/C at synonymous third positions |
| `OverallGc` | `double` | (Gc1+Gc2+Gc3)/3 |
| `CalculateCai` | `double` | CAI in [0,1] |

### 3.3 Preconditions and Validation

Input is read 0-based in steps of 3 (frame 1); a trailing partial codon (< 3 nt) is ignored. Sequences are upper-cased; codons containing any non-ACGT character are skipped (not errors) without shifting the frame. A `null` `DnaSequence`, `null` `GeneticCode` or a `null` reference table throws `ArgumentNullException`; a `null`/empty `string` returns a zeroed `CodonUsageStatistics` (CAI 0). RNA `U` is read as `T` (review 2026-09, CODON-RSCU-001 F10); no IUPAC degeneracy.

## 4. Algorithm

### 4.1 High-Level Steps

1. Count in-frame ACGT codons.
2. Compute RSCU per synonymous family and ENC (Wright 1990).
3. From the codon counts, accumulate G/C at positions 1/2/3 (all codons); for synonymous codons (sense codons whose amino acid has degeneracy > 1 in the genetic code) also accumulate the GC3s numerator/denominator.
4. Convert counts to percentages; `OverallGc` = mean of GC1/GC2/GC3.
5. CAI: build w = referenceRscu / family-max over the synonymous families of the genetic code (skipping single-codon families and stops; w < 0.0001 → 0.01 as CodonW `cai_out`), then geometric mean over scorable codons via log-sum. Canonical core shared with `CodonOptimizer.CalculateCAI` (see [CAI_Calculation.md](../Codon_Optimization/CAI_Calculation.md)).

### 4.2 Decision Rules, Scoring, Reference Tables

- **`EColiOptimalCodons`** — Sharp & Li 1987 E. coli relative-adaptiveness (w) values, transcribed from Biopython v1.79 `SharpEcoliIndex` (e.g. CTG=1, GCC=0.122, CGT=1, AGG=0.002, TTT=0.296); stop codons listed as 0.0 [1][6].
- **`HumanOptimalCodons`** — RSCU derived from the Kazusa H. sapiens [gbpri] per-thousand frequencies (93,487 CDS), `RSCU_j = n·x_j/Σx_k` (e.g. CTG=2.3713, GCC=1.5988, ATG=1.0) [7][8].
- Because CAI rescales each reference value by its family maximum, passing the w table (max 1.0) reproduces w; passing the RSCU table is equivalent (the family-max cancels).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| GetStatistics | O(n) | O(1) | n = sequence length; codon alphabet is constant (64) |
| CalculateCai | O(n) | O(1) | one pass after a constant-size w table |

No substring search / pattern matching is involved, so the repository suffix tree is **not applicable** (N/A) to this unit.

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [CodonUsageAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs)

- `CodonUsageAnalyzer.GetStatistics(string|DnaSequence[, GeneticCode])`: returns the `CodonUsageStatistics` record (default code: NCBI table 1).
- `CodonUsageAnalyzer.CalculateCai(string|DnaSequence, Dictionary<string,double>)`: returns CAI in [0,1].
- `CodonUsageAnalyzer.EColiOptimalCodons` / `HumanOptimalCodons`: reference tables.

### 5.2 Current Behavior

GC1/GC2/GC3 and GC3s are reported as percentages (0–100). GC3s uses only codons whose amino acid is degenerate (degeneracy > 1) in the chosen genetic code, excluding the stop codons (table 1: ATG, TGG and TAA/TAG/TGA excluded), per [2] and CodonW `gc_out`. Cross-check (review 2026-09): 5664 gene × code comparisons against CodonW 1.4.4 (GC3s, L_sym, L_aa, GC1-3) and Biopython `GC123` — 0 mismatches (see Evidence). CAI skips single-codon families and stop codons; a codon whose relative adaptiveness is below 0.0001 is scored as 0.01 (CodonW, Bulmer 1988); when no codon is scorable it returns 0. The suffix tree was not used (no search; single linear scan).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `w_i = f_i / max(f_j)` and `CAI = exp[(1/L) Σ ln w_i]` [1].
- CAI exclusion of non-synonymous (Met, Trp) and termination codons [3][4].
- GC3s over synonymous third positions excluding Met/Trp/stop [2].
- GC1/GC2/GC3 as per-position G/C content [5]; RSCU per Sharp et al. 1986 [8]; ENC per Wright 1990 [9].
- E. coli reference w values transcribed from Sharp & Li 1987 / Biopython [1][6].

**Intentionally simplified:**

- GC3s reported as a percentage; **consequence:** value is 100× the CodonW fraction; the synonymous subset is identical.

**Not implemented:**

- Per-amino-acid "optimal codon" lists, Fop, CBI, GC skew; **users should rely on:** dedicated CodonW/EMBOSS tools or future units.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | GC3s as percentage | Assumption | display units only | accepted | ASM in TestSpec §6 |
| 1a | GC1/GC2/GC3/OverallGc include stop codons | Convention | differs from CodonW GC/GC1-3 for genes with stops | accepted | EMBOSS cusp definition [5] |
| 1b | GC3s ignored `GeneticCode` (table 1 only) | Defect | wrong GC3s/RSCU/ENC for non-standard codes | fixed 2026-09 | `GetStatistics(…, GeneticCode)`, CODON-STATS-001 F18 |
| 2 | zero-w codon skipped | Deviation | edge-case CAI on absent codons | resolved 2026-09 | now w < 0.0001 → 0.01 (CodonW/seqinr/Bulmer 1988), CODON-CAI-001 F13 |
| 3 | reference tables replaced | Deviation (fix) | prior values untraceable | fixed | now Sharp&Li 1987 / Kazusa |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `null` DnaSequence / `null` reference | `ArgumentNullException` | input contract |
| `null` / empty string | zeroed statistics; CAI 0 | input contract |
| trailing partial codon (< 3 nt) | ignored | in-frame parsing |
| only Met/Trp/stop codons | CAI 0; GC3s 0 | no scorable / synonymous codon [1][2] |
| non-ACGT codon | skipped, frame kept | [5] |
| RNA / lower case | U read as T, case-insensitive | CodonW `ident_codon` |

### 6.2 Limitations

DNA or RNA alphabet (no IUPAC ambiguity resolution: ambiguous triplets are skipped). Frame 1 only. CAI quality depends on the supplied reference set (ASM-01); the bundled human table is whole-genome RSCU (Kazusa), not a curated highly-expressed set, so its absolute CAI values are descriptive rather than expression-predictive.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var stats = CodonUsageAnalyzer.GetStatistics("CTGGTTAAA");
// stats.TotalCodons == 3; stats.Gc1 ≈ 66.6667; stats.Gc2 == 0; stats.Gc3 ≈ 33.3333
double cai = CodonUsageAnalyzer.CalculateCai("GCTGCC", CodonUsageAnalyzer.EColiOptimalCodons);
// cai == sqrt(1 * 0.122) == 0.34928498393146 (Ala GCT w=1, GCC w=0.122)
```

**Numerical walk-through:** `ATGGCA` → Met (ATG) excluded from GC3s; Ala (GCA) third base = A (not G/C) ⇒ GC3s = 0/1 = 0%, whereas GC3 over all positions = (G of ATG, A of GCA) = 1/2 = 50%. This is exactly the Met/Trp/stop exclusion of Peden §1.8.2.1.3 [2].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [CodonUsageAnalyzer_GetStatistics_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonUsageAnalyzer_GetStatistics_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [CODON-STATS-001-Evidence.md](../../../docs/Evidence/CODON-STATS-001-Evidence.md)
- Related algorithms: [Relative_Synonymous_Codon_Usage](Relative_Synonymous_Codon_Usage.md), [Effective_Number_of_Codons](Effective_Number_of_Codons.md)

## 8. References

1. Sharp PM, Li W-H. 1987. The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research 15(3):1281–1295. https://doi.org/10.1093/nar/15.3.1281
2. Peden JF. 1999. Analysis of Codon Usage (CodonW reference, §1.8.2.1.3). https://codonw.sourceforge.net/JohnPedenThesisPressOpt_water.pdf
3. CodonW. Codon usage indices. https://codonw.sourceforge.net/Indices.html
4. Charif D, Lobry JR. seqinr `cai` documentation. https://search.r-project.org/CRAN/refmans/seqinr/html/cai.html
5. EMBOSS. `cusp` application documentation. https://www.bioinformatics.nl/cgi-bin/emboss/help/cusp
6. Biopython v1.79. Bio.SeqUtils.CodonUsageIndices (`SharpEcoliIndex`). https://raw.githubusercontent.com/biopython/biopython/biopython-179/Bio/SeqUtils/CodonUsageIndices.py
7. Nakamura Y, Gojobori T, Ikemura T. 2000. Codon usage tabulated from international DNA sequence databases. Nucleic Acids Research 28(1):292. Kazusa, Homo sapiens [gbpri]. https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
8. Sharp PM, Tuohy TMF, Mosurski KR. 1986. Codon usage in yeast. Nucleic Acids Research 14(13):5125–5143. https://doi.org/10.1093/nar/14.13.5125
9. Wright F. 1990. The 'effective number of codons' used in a gene. Gene 87(1):23–29. https://doi.org/10.1016/0378-1119(90)90491-9
