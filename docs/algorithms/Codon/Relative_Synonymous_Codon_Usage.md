# Relative Synonymous Codon Usage (RSCU)

| Field | Value |
|-------|-------|
| Algorithm Group | Codon |
| Test Unit ID | CODON-RSCU-001 |
| Related Projects | Seqeron.Genomics.MolTools, Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-09-28 |

## 1. Overview

Relative Synonymous Codon Usage (RSCU) quantifies, for each codon, how often it is used relative to the usage expected if all synonymous codons of the same amino acid were used equally [1]. It is an exact, specification-driven ratio (not heuristic): a value of 1 indicates no bias, greater than 1 over-representation, and less than 1 under-representation [2][4]. RSCU is the standard normalization for comparing codon preference between genes of different amino-acid composition and is the input to downstream measures such as the Codon Adaptation Index. This unit covers `CalculateRscu` and the supporting `CountCodons`.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The genetic code is degenerate: most amino acids are encoded by 2–6 synonymous codons (Met and Trp by one each; the three stop codons form one family). Synonymous codons are not used with equal frequency; the bias correlates with tRNA abundance and gene expression level [1]. Raw codon counts cannot be compared across genes because amino-acid composition differs; RSCU removes amino-acid composition by normalizing within each synonymous family [1][3].

### 2.2 Core Model

For amino acid *i* with degeneracy *n_i* (number of synonymous codons) and observed count *x_{i,j}* of its *j*-th codon, the RSCU is [3] (equivalent form in [2]):

```
RSCU_{i,j} = x_{i,j} / ((1/n_i) * Σ_{k=1}^{n_i} x_{i,k})
           = (n_i * x_{i,j}) / Σ_{k=1}^{n_i} x_{i,k}
```

The denominator `(1/n_i)·Σx` is the expected count under uniform synonymous usage; the numerator is the observed count [1][4]. RSCU = 1 means observed equals expected (no bias) [2][4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | RSCU = (n_i·x_{i,j})/Σx for a present family | Definition [3] |
| INV-02 | Equal usage within a family ⇒ RSCU = 1 | x_{i,j}=Σx/n_i ⇒ ratio = 1 [2][4] |
| INV-03 | 0 ≤ RSCU ≤ n_i | x_{i,j} ∈ [0, Σx]; max when one codon used [3] |
| INV-04 | Σ_j RSCU_{i,j} = n_i for a present family | Σ_j n_i·x_{i,j}/Σx = n_i·(Σx)/Σx [3] |
| INV-05 | Single-codon family (Met, Trp) ⇒ RSCU = 1 when present | n_i = 1 ⇒ ratio = x/x = 1 [3][5] |

### 2.5 Comparison with Related Methods

| Aspect | RSCU | Raw codon frequency (per-1000) |
|--------|------|--------------------------------|
| Normalization | within synonymous family | over all codons |
| Comparable across genes of different aa composition | yes | no |
| No-bias reference value | 1.0 | n/a |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | `DnaSequence` or `string` | required | Coding DNA or RNA sequence | Read as non-overlapping triplets from offset 0; `string` overload uppercases, reads U as T (CodonW `ident_codon`) and skips triplets with any other symbol without shifting the frame; `DnaSequence` rejects invalid bases at construction |
| code | `GeneticCode` | `GeneticCode.Standard` | Genetic code defining the synonymous families | any NCBI table (`GeneticCode.GetByTableNumber`) |
| codonCounts | `IReadOnlyDictionary<string,int>` | — | Alternative input: codon counts keyed by uppercase DNA codon | other keys ignored |
| pseudocount | `double` | 0 (overloads without it) | Added to every codon count before the per-family normalisation: RSCU = n·(x + c)/(Σx + n·c) | finite, ≥ 0 (`ArgumentOutOfRangeException`); 0 = CodonW, bit-identical to the overloads without it |
| level | `RscuFamilyLevel` | `AminoAcid` (overloads without it) | Synonymous grouping: `AminoAcid` (Sharp & Li / CodonW families) or `Subfamily` (cubar `level = "subfam"`: amino acid + first two codon bases) | defined enum value (`ArgumentOutOfRangeException`) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return (CalculateRscu) | `Dictionary<string,double>` | all 64 codons (DNA spelling) → RSCU value; codons of an absent family map to 0 (CodonW) |
| return (CountCodons) | `Dictionary<string,int>` | codon → occurrence count over counted triplets |

### 3.3 Preconditions and Validation

Null `DnaSequence`, `GeneticCode` or count table throws `ArgumentNullException`. Empty/null `string` returns an empty dictionary. Input is case-insensitive (string overload uppercases; `DnaSequence` normalizes at construction). Codons are 0-based non-overlapping triplets; trailing 1–2 bases are ignored. RNA U is read as T (CodonW 1.4.4 `codon_us.c` `ident_codon`). Any triplet containing a character outside {A,C,G,T,U} is excluded from counts.

## 4. Algorithm

### 4.1 High-Level Steps

1. Count non-overlapping triplets (`CountCodons`); U is read as T; exclude triplets containing other non-ACGT characters.
2. Group all 64 codons by the amino acid `GeneticCode.CodonTable` assigns them (default Standard, NCBI table 1: Met/Trp single-codon; stop codons '*' as one family). Context-dependent stops of NCBI tables 27/28/31 belong to the family of the amino acid they encode (Biopython `forward_table`).
3. For each family, compute the family total Σx and divide it equally over the n_i synonymous codons to get the expected count.
4. RSCU = observed / expected for each codon; if the family total is 0 (absent family) set every codon of that family to 0.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The genetic code (codon → amino acid) is the only reference table; it determines the synonymous families (for table 1: degeneracies 1,2,3,4,6 and the 3-fold stop family; for table 2: e.g. Met {ATA,ATG}, Trp {TGA,TGG}, stop {TAA,TAG,AGA,AGG}). Families are taken from the canonical `GeneticCode` (no private copy), as in CodonW ("RSCU values are genetic code dependent", `rscu_usage_out`). Family degeneracy n_i is derived from the table, not hard-coded per amino acid.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CountCodons | O(n) | O(k) | n = sequence length; k = distinct codons (≤64) |
| CalculateRscu | O(n) | O(64) | one pass to count + constant-size family aggregation |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [CodonUsageAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs)

- `CodonUsageAnalyzer.CalculateRscu(DnaSequence)` / `CalculateRscu(string)`: RSCU per codon under the Standard code.
- `CodonUsageAnalyzer.CalculateRscu(DnaSequence, GeneticCode)` / `CalculateRscu(string, GeneticCode)`: RSCU under any NCBI table.
- `CodonUsageAnalyzer.CalculateRscu(IReadOnlyDictionary<string,int>, GeneticCode)`: RSCU from a codon-count table (canonical core; e.g. pooled counts over many genes).
- `CodonUsageAnalyzer.CalculateRscu(DnaSequence | string | IReadOnlyDictionary<string,int>, GeneticCode, double pseudocount)`: RSCU with a pseudocount (codonbias / cubar smoothing).
- `CodonUsageAnalyzer.CalculateRscu(DnaSequence | string | IReadOnlyDictionary<string,int>, GeneticCode, double pseudocount, RscuFamilyLevel level)`: RSCU with a pseudocount at the amino-acid or cubar sub-family level.
- `CodonUsageAnalyzer.CountCodons(DnaSequence)` / `CountCodons(string)`: counts non-overlapping codon occurrences.

### 5.2 Current Behavior

Counting uses a simple linear scan over non-overlapping triplets (`CountCodonsCore`). This is not a substring-search problem (no pattern is being located inside the text; each fixed-width triplet is read once at a known offset), so the repository suffix tree is **not** used — it offers no benefit for a single linear pass and would add O(n) construction overhead with no query reuse. RSCU groups codons by amino acid and computes `observed / (total/n_i)`, which is algebraically `(n_i·observed)/total`. Absent families (total = 0) return 0 for every member.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- RSCU = (n_i · x_{i,j}) / Σ_{k} x_{i,k} for present families [1][2][3].
- No-bias value 1, range [0, n_i], family-sum n_i, single-codon ⇒ 1 (INV-01..INV-05) [2][3][4][5].
- CodonW conventions [6]: all 64 codons reported, stop codons as one family, absent family → 0, families per genetic code, U ≡ T. Cross-checked against the CodonW 1.4.4 binary (built from source) on 155 genes × 8 codes (79,360 values, max |Δ| = 5e-4 = its 3-decimal rounding) and against a Python port with Biopython tables for all 27 NCBI tables (472 inputs, max |Δ| 9e-16).

**Conventions and options:**

- Absent-family 0/0 case: the default pseudocount is c = 0, as CodonW (`rscu_usage_out` prints 0.000 for every codon of an absent amino acid), and such a family is reported as 0 (codonbias / cubar give NaN at c = 0); a pseudocount option is available (below), with which an absent family is uniform (1). Only absent families are affected, never the RSCU of an observed codon.

- Optional pseudocount `c` (overloads with `double pseudocount`): RSCU_j = n·(x_j + c)/Σ_k(x_k + c), applied per codon before the family normalisation, exactly as codonbias 0.5.0 `RelativeSynonymousCodonUsage(pseudocount=c)` (`counts + pseudocount`, normalised per amino acid, divided by the uniform 1/n) and cubar `est_rscu(pseudo_cnt = c, level = "amino_acid")` (`prop = (cts + c)/sum(cts + c)`, `rscu = prop/mean(prop)`) [5][7]. Families are whole amino acids plus the stop family (codonbias `ignore_stop=False`, cubar `incl_stop = TRUE`); with c > 0 an absent family is uniform (1). 0 mismatches > 1e-12 against codonbias over 200 random genes × all 27 NCBI tables × c = 0.25/1/3 (1 036 800 values, max |Δ| 4.4e-16). cubar's defaults differ: `level = "subfam"` splits 6-fold amino acids by the first two bases (Leu CT*/TT*, Ser TC*/AG*, Arg CG*/AG*), stops are excluded, `pseudo_cnt = 1`; codonbias defaults to `pseudocount = 1` and excludes stops. With c = 0 both tools give NaN (0/0) for absent families, where this library keeps CodonW's 0.

- Optional cubar sub-family level (`RscuFamilyLevel.Subfamily`): cubar `get_codon_table()` labels each codon `subfam = paste(amino_acid, substr(codon, 1, 2), sep = "_")` and `est_rscu(level = "subfam")` (cubar's default) normalises within each sub-family [5]. Every amino acid of the chosen genetic code is split by that rule — table 1: Leu_TT {TTA, TTG} / Leu_CT {CTN}, Ser_TC {TCN} / Ser_AG {AGT, AGC}, Arg_CG {CGN} / Arg_AG {AGA, AGG}, stops \*_TA {TAA, TAG} / \*_TG {TGA}; Ile_AT and Met_AT stay separate (different amino acids); table 2: stops \*_TA and \*_AG {AGA, AGG}. Stops are reported (cubar `incl_stop = TRUE`, all 64 codons). Verified against an independent Python port of cubar @350b295 (`R/codon_info.R`, `R/codon_weight.R`): 7 genes × tables 1/2 × c = 0/0.5/1 × 64 codons locked (`TestData/Codon/cubar_rscu_subfam.tsv`, 1e-12); cubar's NaN for an absent sub-family at c = 0 is 0 here.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Absent family returns 0 by default (pseudocount 0) | Assumption | Only affects amino acids absent from input; no canonical value exists there | accepted | Documented in Evidence; pass `pseudocount` > 0 for codonbias / cubar smoothing [5][7] |
| 2 | Stop codons grouped as one family (3-fold in table 1; at the sub-family level split by their first two bases) | Convention | Does not affect any amino-acid RSCU | accepted | CodonW does the same [6]; CodonU / codon-bias / cubar report the 61 sense codons only [5] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null `DnaSequence` | `ArgumentNullException` | input guard |
| Empty / null `string` | empty dictionary | input guard |
| Trailing 1–2 bases | ignored | non-overlapping triplets |
| RNA input (U) | read as T; same result as the DNA spelling | CodonW `ident_codon` [6] |
| Non-ACGT triplet (string overload) | excluded from counts, frame kept | `IsValidCodon` over {A,C,G,T} after U→T |
| Alternative genetic code | families from that table | CodonW `-code`, Biopython tables |
| Single-codon family present | RSCU = 1 | n_i = 1 [3][5] |
| Absent family | RSCU = 0 for all members (pseudocount c > 0: 1) | 0/0 convention (§5.3) |

### 6.2 Limitations

The MCP `rscu` tool takes `genetic_code`, `pseudocount` and `family_level` (`AminoAcid` default, `Subfamily` = cubar `subfam`). RSCU is a within-family bias measure and does not by itself indicate expression level or adaptation (use CAI for that).

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through (Leu, 6-fold):** input `CTGCTGCTGCTA` → CTG×3, CTA×1, family total = 4, n_i = 6.
RSCU(CTG) = 6·3/4 = 4.5; RSCU(CTA) = 6·1/4 = 1.5; RSCU(TTA)=RSCU(TTG)=RSCU(CTT)=RSCU(CTC) = 0; Σ over the family = 6.0 (= n_i, INV-04).

**API usage example:**

```csharp
var rscu = CodonUsageAnalyzer.CalculateRscu("CTGCTGCTGCTA");
// rscu["CTG"] == 4.5, rscu["CTA"] == 1.5
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [CodonUsageAnalyzer_CalculateRscu_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonUsageAnalyzer_CalculateRscu_Tests.cs) — covers INV-01..INV-05
- Tests: [CodonUsageAnalyzer_RscuPseudocount_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonUsageAnalyzer_RscuPseudocount_Tests.cs) — pseudocount, locked to codonbias 0.5.0 (`TestData/Codon/codonbias_rscu_pseudocount.tsv`)
- Tests: [CodonUsageAnalyzer_RscuSubfamily_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonUsageAnalyzer_RscuSubfamily_Tests.cs) — cubar sub-family level (`TestData/Codon/cubar_rscu_subfam.tsv`)
- Evidence: [CODON-RSCU-001-Evidence.md](../../../docs/Evidence/CODON-RSCU-001-Evidence.md)

## 8. References

1. Sharp P.M., Tuohy T.M.F., Mosurski K.R. 1986. Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed genes. Nucleic Acids Research 14(13):5125-5143. https://doi.org/10.1093/nar/14.13.5125
2. Suzuki H. et al. / LIRMM. RSCU RS — Measuring the bias in codon usage. https://www.lirmm.fr/~rivals/rscu/
3. GenomicSig (CRAN). RSCU: Relative Synonymous Codon Usage. https://rdrr.io/cran/GenomicSig/man/RSCU.html
4. Charif D., Lobry J.R. seqinr — `uco`: Codon usage indices. https://search.r-project.org/CRAN/refmans/seqinr/html/uco.html
5. cubar (CRAN). `est_rscu`: Estimate Relative Synonymous Codon Usage. https://rdrr.io/cran/cubar/man/est_rscu.html; source R/codon_info.R (`get_codon_table`: `subfam`), R/codon_weight.R (`est_rscu`), https://github.com/mt1022/cubar (commit 350b295).
6. Peden J.F. 1999. CodonW 1.4.4 source, `codon_us.c` (`rscu_usage_out`, `cutab_out`, `ident_codon`). http://archive.ubuntu.com/ubuntu/pool/universe/c/codonw/codonw_1.4.4.orig.tar.gz
7. codonbias 0.5.0 (PyPI `codon-bias`), `codonbias/scores.py` `RelativeSynonymousCodonUsage` (`_weights_and_counts`), `codonbias/stats.py` `CodonCounter.get_aa_table`; cubar source `R/codon_weight.R` `est_rscu` (github.com/mt1022/cubar @ 350b295).

- Related algorithm: [Relative_Synonymous_Codon_Usage.md](../Annotation/Relative_Synonymous_Codon_Usage.md) (ANNOT-CODONUSAGE-001 — the Annotation-side implementation of the same concept).
