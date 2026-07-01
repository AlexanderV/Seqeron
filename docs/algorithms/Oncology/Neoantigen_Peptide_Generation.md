# Neoantigen Candidate Peptide Window Generation

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / Neoantigen Prediction |
| Test Unit ID | ONCO-NEO-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Given a somatic missense mutation (a single amino-acid substitution) in a protein, this algorithm enumerates
the candidate MHC class I neoantigen peptides: every fixed-length window (8–14-mer by default) of the mutant
protein that *spans* the substituted residue, each paired with the wild-type peptide occupying the same
coordinates (the agretope). It is the deterministic, specification-driven *windowing* step of neoantigen
prediction [1][2]; it does not score MHC binding (IC50 / presentation rank), which requires a trained model
(NetMHCpan) and is caller-supplied or out of scope [3][4]. The implementation status is **Framework**: it
produces the well-defined peptide candidates that a downstream binding predictor consumes.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Tumour-specific somatic mutations can create neoantigens — mutant peptides presented on MHC molecules and
recognised as non-self by T cells. For MHC class I, presented ligands are short peptides, predominantly 9-mers
but ranging 8–11 residues, with length preference varying by HLA allele [3][4]. A missense mutation alters one
residue of the protein; only peptides that contain that residue can differ from self, so neoantigen prediction
considers exactly the peptides spanning the mutation [2].

### 2.2 Core Model

Let the wild-type protein be `P` of length `L` with 1-based positions, and a missense mutation at position `p`
(1 ≤ p ≤ L) changing residue `P[p]` to a different amino acid `a`. The mutant protein `P'` equals `P` except
`P'[p] = a`.

For a peptide length `k`, a window starting at 1-based position `s` covers positions `s … s+k-1`. It **spans**
the mutation iff `s ≤ p ≤ s+k-1`, i.e.

> s ∈ [ max(1, p − k + 1), min(p, L − k + 1) ]

The candidate mutant peptides of length `k` are `P'[s … s+k-1]` for each such `s`; the paired wild-type peptide
(agretope) is `P[s … s+k-1]` at the same coordinates [2]. The default length range is k ∈ {8, 9, 10, 11, 12,
13, 14} for MHC class I — the NetMHCpan-4.1 class I peptide window (Reynisson et al. 2020) [1]. The shorter
8–11-mers are exactly the set contained in the 21-mer ±10-flank window that ProGeo-neo constructs around the
substitution [2]; longer 12–14-mers extend symmetrically as the protein bounds permit.

When `p` is at least `k − 1` residues from both termini and `L ≥ p + k − 1`, the bounds give exactly `k`
windows of length `k`. Near a terminus the count is truncated to the windows that fit [2].

The number of mutant/wild-type pairs that differ from one another equals the number of windows (each window
contains the substituted residue), and within each pair the peptides differ at exactly one position — the
mutation offset `p − s`. The agretopicity / differential agretopic index downstream compares the binding of
these two peptides [5].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every peptide length is in [minLength, maxLength] (8–14 by default) | enumeration loop ranges over the requested lengths only [1] |
| INV-02 | Every peptide spans the mutation: `StartPosition + MutationOffset == p` and `0 ≤ MutationOffset < Length` | start bounds `s ∈ [max(1,p−k+1), min(p,L−k+1)]` [2] |
| INV-03 | Mutant and wild-type peptides have equal length and differ at exactly the mutation offset | `P'` differs from `P` only at position `p`, which lies in the window |
| INV-04 | `MutantPeptide[offset] == a` (mutant residue); `WildTypePeptide[offset] == P[p]` (original) | substitution applied only at `p` [5] |
| INV-05 | Interior mutation (≥ k−1 from both ends) yields exactly k windows of length k | bound width `min(p,L−k+1) − max(1,p−k+1) + 1 = k` |
| INV-06 | Output ordered by length ascending then start ascending | nested loops k outer, start inner |

### 2.5 Comparison with Related Methods

| Aspect | Window generation (this unit) | MHC binding prediction (ONCO-MHC-001) |
|--------|-------------------------------|----------------------------------------|
| Nature | Deterministic, combinatorial | Learned model (neural network) |
| Output | Candidate peptide + agretope pairs | IC50 / %rank binding score |
| Requires trained weights | No | Yes (NetMHCpan) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| wildTypeProtein | string | required | Wild-type protein, one-letter AA codes | non-null, non-empty |
| mutantResidue | char | required | Substituted (mutant) amino acid | must differ from WT residue at the position |
| mutationPosition | int | required | 1-based position of the substitution | 1 ≤ position ≤ length |
| minLength | int | 8 | Minimum peptide length | ≥ 1 |
| maxLength | int | 11 | Maximum peptide length | ≥ minLength |

### 3.2 Output / Return Value

`IReadOnlyList<NeoantigenPeptide>` ordered by length ascending then start position ascending.

| Field | Type | Description |
|-------|------|-------------|
| Length | int | Peptide length k |
| StartPosition | int | 1-based position of the window's first residue |
| MutantPeptide | string | k-mer from the mutant protein |
| WildTypePeptide | string | k-mer from the wild-type protein at the same coordinates (agretope) |
| MutationOffset | int | 0-based offset of the substituted residue within the window |

### 3.3 Preconditions and Validation

1-based protein coordinates. Null protein → `ArgumentNullException`. Empty protein, `mutantResidue` equal to
the wild-type residue (not a substitution), `minLength < 1`, or `maxLength < minLength` → `ArgumentException`.
`mutationPosition` outside [1, length] → `ArgumentOutOfRangeException`. A requested length `k > L` is skipped
(no window fits); if no length yields a window the result is empty. Sequences are treated as opaque
one-letter-code strings (no alphabet validation; case preserved).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; confirm a genuine substitution at `mutationPosition`.
2. Build the mutant protein `P'` by copying `P` and replacing the residue at the mutation position.
3. For each length `k` in [minLength, maxLength] (skip if `k > L`): compute the spanning start range
   `s ∈ [max(0, p0−k+1), min(p0, L−k)]` (0-based, `p0 = p−1`).
4. For each start, emit the mutant and wild-type k-mers and the mutation offset.
5. Return all peptides ordered by length then start.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Parameter | Value | Source |
|-----------|-------|--------|
| Default MHC class I min length | 8 | Hundal et al. (2020) [1]; NetMHCpan-4.1 [4] |
| Default MHC class I max length | 11 | Hundal et al. (2020) [1] |
| Window-spanning rule | start ∈ [max(1,p−k+1), min(p,L−k+1)] | ProGeo-neo 21-mer ±10-flank [2] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| GenerateNeoantigenPeptides | O(Σ_k k²) ≈ O((max−min+1)·k_max²) | O(total peptide chars) | per length, ≤ k windows each of k chars; constant for fixed 8–14 range |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.GenerateNeoantigenPeptides(string, char, int, int, int)`: enumerates the candidate peptides.
- `OncologyAnalyzer.NeoantigenPeptide`: record struct holding one mutant/wild-type peptide pair.
- `OncologyAnalyzer.MhcClassIMinPeptideLength` / `MhcClassIMaxPeptideLength`: 8 / 11.

### 5.2 Current Behavior

The mutant protein is materialised once per call; each length's windows are extracted by `Substring`. Output
order is length-ascending then start-ascending. No suffix tree is used (see 5.4). Sequences are not
alphabet-validated — any one-letter codes are accepted and case is preserved.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Enumeration of every length-k window spanning the substituted residue, for k = 8–14 by default [1][2].
- Wild-type/mutant agretope pairing at identical coordinates [2][5].

**Intentionally simplified:**

- (none)

**Not implemented:**

- MHC binding-affinity / IC50 / %rank scoring; **users should rely on:** a trained predictor (NetMHCpan,
  ONCO-MHC-001) supplied with these peptides [3][4].
- Frameshift / indel / fusion neopeptide generation; **users should rely on:** dedicated indel/fusion
  translation (out of scope for this missense unit) [1].

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Suffix tree not used | Deviation | none | accepted | The "find windows spanning a position" task is a bounded arithmetic range over one short protein, not multi-query exact-match search; the repository suffix tree does not fit. Correctness unaffected. |
| 2 | Class named OncologyAnalyzer, not NeoantigenPredictor | Assumption | none (API naming) | accepted | Follows existing Oncology project layout and the task requirement; checklist placeholder name superseded. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Mutation at N-/C-terminus | Only the windows that fit while spanning it (truncated count) | ProGeo-neo builds the flanked window "if possible" [2] |
| Requested length > protein length | That length skipped; shorter lengths still returned | no window of that length fits |
| mutantResidue == wild-type residue | ArgumentException | not a missense substitution [1] |
| Single length range (min==max) | Only that length returned | range respected |

### 6.2 Limitations

Single-residue missense substitutions only. No binding prediction, no proteasomal cleavage / TAP transport
modelling, no expression/VAF filtering, no alphabet validation. Output is the candidate set, not ranked or
filtered neoantigens.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Wild-type MKTAYIAKQRSTVWLNDEFGH, missense Y5C; default 8–14-mers.
var peptides = OncologyAnalyzer.GenerateNeoantigenPeptides("MKTAYIAKQRSTVWLNDEFGH", 'C', 5);
// 35 peptides (5 per length 8..14). First 8-mer: MutantPeptide "MKTACIAK", WildTypePeptide "MKTAYIAK",
// StartPosition 1, MutationOffset 4.
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_GenerateNeoantigenPeptides_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_GenerateNeoantigenPeptides_Tests.cs) — covers `INV-01`–`INV-06`
- Evidence: [ONCO-NEO-001-Evidence.md](../../../docs/Evidence/ONCO-NEO-001-Evidence.md)

## 8. References

1. Hundal J, Kiwala S, McMichael J, et al. 2020. pVACtools: A Computational Toolkit to Identify and Visualize Cancer Neoantigens. Cancer Immunology Research 8(3):409–420. https://doi.org/10.1158/2326-6066.CIR-19-0401
2. Li Y, Wang G, Tan X, et al. 2020. ProGeo-neo: a customized proteogenomic workflow for neoantigen prediction and selection. BMC Medical Genomics 13:52. https://doi.org/10.1186/s12920-020-0683-4
3. Jurtz V, Paul S, Andreatta M, et al. 2017. NetMHCpan-4.0: Improved Peptide-MHC Class I Interaction Predictions. Journal of Immunology 199(9):3360–3368. https://doi.org/10.4049/jimmunol.1700893
4. NetMHCpan-4.1 web service. DTU Health Tech. https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
5. Wells DK, van Buuren MM, Dang KK, et al. 2020. Key Parameters of Tumor Epitope Immunogenicity Revealed Through a Consortium Approach Improve Neoantigen Prediction (TESLA). Cell 183(3):818–834. https://doi.org/10.1016/j.cell.2020.09.015
