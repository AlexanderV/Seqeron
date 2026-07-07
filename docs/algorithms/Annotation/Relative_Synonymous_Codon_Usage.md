# Relative Synonymous Codon Usage (RSCU)

| Field | Value |
|-------|-------|
| Algorithm Group | Annotation |
| Test Unit ID | ANNOT-CODONUSAGE-001 |
| Related Projects | Seqeron.Genomics.Annotation, Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Relative Synonymous Codon Usage (RSCU) measures the bias in usage of synonymous codons within a set of coding sequences, independent of amino-acid composition and overall codon frequency [1]. For each codon it reports the observed count relative to what would be expected if all synonymous codons for the same amino acid were used equally; a value of 1.0 means no bias, above 1.0 a preferred codon, below 1.0 an under-represented one [2][3]. It is an exact, deterministic ratio (not a heuristic), used to characterise translational selection, optimise heterologous expression, and compare codon strategies across genes or organisms.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The standard genetic code is degenerate: 61 sense codons encode 20 amino acids, so most amino acids have several synonymous codons (NCBI translation table 1) [5]. Organisms do not use synonymous codons uniformly; RSCU normalises away amino-acid composition so that only the relative preference within each synonymous family is captured [1].

### 2.2 Core Model

For an amino acid *i* with `n_i` synonymous codons, let `x_{i,j}` be the observed count of its *j*-th codon. The RSCU of that codon is [3] (originating with Sharp & Li 1986 [1]):

```
RSCU_{i,j} = (n_i · x_{i,j}) / Σ_j x_{i,j}
```

where the denominator sums the observed counts over all `n_i` synonymous codons of amino acid *i* [3]. Equivalently, RSCU is the ratio of the observed codon frequency to the expected frequency under equal synonymous usage [2]. Counts are pooled across the whole reference set of coding sequences before the ratio is computed [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | For every observed synonymous family, Σ RSCU over its codons = `n_i` | Σ_j (n_i·x_{i,j}/Σx) = n_i·(Σx/Σx) = n_i [3] |
| INV-02 | A single-codon amino acid (Met=ATG, Trp=TGG) always has RSCU = 1.0 | `n_i = 1` ⇒ 1·x/x = 1 [3][5] |
| INV-03 | Each RSCU value lies in [0, `n_i`] | bounded ratio per the formula [3] |
| INV-04 | Stop codons never appear in the output | RSCU is defined over sense codons only [4][5] |
| INV-05 | Uniform usage within a family ⇒ every member RSCU = 1.0 | equal counts ⇒ n_i·x/(n_i·x) = 1 [2][3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `codingSequences` | `IEnumerable<string>` | required | Coding (CDS) DNA sequences | Read in frame from index 0 in steps of 3; partial trailing codon ignored; only A/C/G/T codons counted; case-insensitive |
| `code` | `GeneticCode` | `GeneticCode.Standard` (table 1) | Genetic code defining synonymous families | non-null |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `IReadOnlyDictionary<string, double>` | Map from each sense codon (uppercase DNA) to its RSCU value |

### 3.3 Preconditions and Validation

Null `codingSequences` or null `code` throw `ArgumentNullException`. Null/empty individual sequences are skipped. Input is uppercased (case-insensitive). Codons are read 0-based in non-overlapping triplets; a partial trailing codon (length not a multiple of 3) is ignored. The genetic-code table is RNA-keyed (U); DNA codons (T) are reconciled internally (T↔U). Stop codons are excluded. An entirely unobserved synonymous family yields RSCU 0.0 for each of its codons (no division by zero; the CAI 0.5 pseudocount is not applied — see 5.4).

## 4. Algorithm

### 4.1 High-Level Steps

1. Pool codon counts: for each coding sequence, step through in-frame triplets and tally A/C/G/T codons across all sequences [4].
2. Build synonymous families: group the genetic code's sense codons by the amino acid they encode (stop codons excluded) [5].
3. For each family compute the family total Σ_j x_{i,j} and set RSCU_{i,j} = n_i·x_{i,j}/Σ for each member (0.0 if the family total is 0) [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| RSCU over all sequences | O(N) | O(1) | N = total nucleotides; the 61 sense codons / 20 families are a fixed constant |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomeAnnotator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs)

- `GenomeAnnotator.GetCodonUsage(IEnumerable<string>)`: RSCU over the Standard genetic code (NCBI table 1).
- `GenomeAnnotator.GetCodonUsage(IEnumerable<string>, GeneticCode)`: RSCU over a caller-supplied genetic code.
- `GenomeAnnotator.GetCodonUsage(string)`: pre-existing raw codon-count method (not RSCU); retained unchanged for backward compatibility.

### 5.2 Current Behavior

Counts are pooled across the input collection before computing RSCU (transcriptome-wide RSCU), matching the reference implementation [4]. The genetic-code table from `Seqeron.Genomics.Core.GeneticCode` is RNA-keyed; each codon key is converted U→T so DNA inputs map correctly. No substring search is performed (this is a single linear in-frame scan), so the repository suffix tree is **not** applicable here.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- The RSCU formula RSCU = n_i·x/Σx exactly as stated [3] (Sharp & Li 1986 [1]).
- Sense-codon-only families with stop codons excluded [4][5].
- Pooling of counts over the whole reference set [4].

**Intentionally simplified:**

- (none)

**Not implemented:**

- CAI 0.5 pseudocount for zero-count codons; **users should rely on:** a dedicated Codon Adaptation Index routine — the 0.5 substitution is a CAI convention (Sharp & Li 1987), not part of plain RSCU [4].

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Empty synonymous family → RSCU 0.0 for each member | Assumption | Only affects amino acids with no observed codon; avoids division by zero | accepted | Base RSCU is undefined when Σ=0; reporting 0.0 (no preferred codon) keeps output total. CAI 0.5 pseudocount deliberately not used [4]. |
| 2 | Default genetic code = Standard (table 1) | Assumption | API default; non-standard codes available via overload | accepted | Overload accepts a `GeneticCode` [5] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null `codingSequences` | `ArgumentNullException` | input validation |
| Empty collection / all-empty strings | empty dictionary | nothing observed |
| Partial trailing codon | ignored | in-frame triplet stepping [4] |
| Lower-case input | counted (uppercased first) | case-insensitive contract |
| Stop codon present (e.g. TAA) | absent from output | sense codons only [4][5] |
| Met/Trp codon | RSCU = 1.0 | n_i = 1 [3][5] |

### 6.2 Limitations

RSCU is computed only over A/C/G/T codons; codons containing ambiguity codes (N, IUPAC) are not counted. The metric describes within-family bias only and says nothing about absolute expression or amino-acid composition. For codon-adaptation scoring with pseudocounts, use a dedicated CAI implementation rather than this RSCU output.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var rscu = GenomeAnnotator.GetCodonUsage(new[] { "CTTCTTCTGTTA" });
// Leucine family (n_i = 6): CTT=3.0, CTG=1.5, TTA=1.5, TTG=0, CTC=0, CTA=0
```

**Numerical walk-through:** CDS `CTTCTTCTGTTA` → codons CTT, CTT, CTG, TTA (all Leu). Family counts CTT=2, CTG=1, TTA=1; Σ=4; n_i=6. RSCU(CTT)=6·2/4=3.0, RSCU(CTG)=6·1/4=1.5, RSCU(TTA)=6·1/4=1.5, others 0.0. Σ over family = 6.0 = n_i (INV-01) [3].

### 7.2 Applications and Use Cases

- **Heterologous expression / codon optimisation:** preferred codons (RSCU > 1) guide synthetic-gene design.
- **Comparative genomics:** RSCU profiles distinguish translational-selection regimes between genes or genomes [1].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomeAnnotator_GetCodonUsage_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/GenomeAnnotator_GetCodonUsage_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [ANNOT-CODONUSAGE-001-Evidence.md](../../../docs/Evidence/ANNOT-CODONUSAGE-001-Evidence.md)

## 8. References

1. Sharp P.M., Li W.-H. 1986. Codon usage in regulatory genes in Escherichia coli does not reflect selection for 'rare' codons. Nucleic Acids Research 14(19):7737–7749. https://doi.org/10.1093/nar/14.19.7737
2. Analysis of synonymous codon usage and evolution of begomoviruses. PMC2528880. https://pmc.ncbi.nlm.nih.gov/articles/PMC2528880/
3. Rivals E. et al. RSCU RS: Measuring the bias in codon usage. LIRMM, Université de Montpellier. https://www.lirmm.fr/~rivals/rscu/
4. SouradiptoC. CodonU, `CodonU/analyzer/internal_comp.py`, function `rscu`. https://github.com/SouradiptoC/CodonU/blob/master/CodonU/analyzer/internal_comp.py
5. NCBI. The Genetic Codes — Standard Code (transl_table=1). https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

- Related algorithm: [Relative_Synonymous_Codon_Usage.md](../Codon/Relative_Synonymous_Codon_Usage.md) (CODON-RSCU-001 — canonical RSCU doc; this page covers the Annotation-side `GenomeAnnotator` computation of the same concept).
