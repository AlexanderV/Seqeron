# Signal Peptide Cleavage-Site Prediction (von Heijne Weight Matrix)

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinMotif |
| Test Unit ID | PROTMOTIF-SP-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Predicts the cleavage site between an N-terminal signal peptide and the mature secreted/membrane protein, using von Heijne's (1986) position-specific weight-matrix method [1]. Each candidate position is scored by summing log-odds residue weights over a 15-residue window (positions −13..+2 relative to the cleavage site), and the highest-scoring position is reported as the predicted mature-protein start. The scoring is a faithful re-implementation of the EMBOSS `sigcleave` reference program [2][3], using the eukaryotic (161-sequence) matrix by default and the prokaryotic (36-sequence) matrix on request. The method is a statistical heuristic: it always returns a best site for an in-window sequence, and a configurable weight threshold (default 3.5) flags whether that site is a likely signal peptide [2].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Secretory and membrane proteins carry a transient N-terminal signal peptide that is removed by signal peptidase after translocation. The cleavage site shows a conserved pattern of small residues at positions −1 and −3 (the "(−3,−1) rule") and a hydrophobic core upstream [1]. von Heijne (1986) captured this pattern as a residue-frequency matrix over the 15 positions flanking the cleavage site [1].

### 2.2 Core Model

For an aligned set of signal peptides, let `C(a, p)` be the count of residue `a` at window position `p ∈ {−13,…,−1,+1,+2}`, and let `E(a)` be the residue's background ("Expect") count. The position-specific log-odds weight is [1][3]

```
W(a, p) = ln( C(a, p) / E(a) )
```

with a zero count replaced by a pseudocount before the log: `1.0 × 10⁻¹⁰` at the conserved columns −3 and −1 (forcing a strong penalty when the conserved small residue is absent), and `1.0` at all other columns [3]. The score of a candidate cleavage site whose +1 residue is sequence index `i` is the sum over the window [3]:

```
S(i) = Σ_{p = −13}^{+2} W( seq(i+p'), p )      (p' maps window position p to sequence offset)
```

The predicted cleavage site is `argmax_i S(i)`; cleavage occurs between the −1 and +1 residues, so the mature protein begins at index `i` [2][3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `S(i)` is the sum of log-odds weights over 15 window columns | matches `sigcleave.c` scoring loop [3] |
| INV-02 | Returned site is the global argmax of `S(i)` | `maxweight`/`maxsite` update [3] |
| INV-03 | Cleavage is between mature-start−1 and mature-start (1-based) | EMBOSS `+1`/`−1` convention [2][3] |
| INV-04 | Output is case-independent | input is upper-cased before scoring |
| INV-05 | `IsLikelySignalPeptide ⇔ S ≥ minWeight` (default 3.5) | EMBOSS `-minweight` semantics [2] |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | von Heijne weight matrix (this) | Hydrophobicity scan (Kyte-Doolittle) |
|--------|--------------------------------|--------------------------------------|
| Output | single cleavage position + score | hydrophobic regions |
| Basis | position-specific residue frequencies [1] | average hydropathy [—] |
| Use | signal-peptide cleavage prediction | transmembrane / hydrophobic core |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `proteinSequence` | `string` | required | Amino-acid sequence, one-letter codes | case-insensitive; ≥ 15 residues for a result |
| `prokaryote` | `bool` | `false` | Use the prokaryotic matrix instead of eukaryotic | — |
| `minWeight` | `double` | `3.5` | Acceptance threshold for the likelihood flag | EMBOSS recommends ≥ 3.5 [2] |

### 3.2 Output / Return Value

`SignalPeptide?` — `null` when no result, otherwise:

| Field | Type | Description |
|-------|------|-------------|
| `CleavagePosition` | `int` | 1-based first residue of the mature protein (cleavage between `CleavagePosition−1` and `CleavagePosition`) |
| `Score` | `double` | von Heijne weight-matrix score at the best site |
| `SignalSequence` | `string` | residues `1..CleavagePosition−1` (the predicted signal peptide) |
| `WindowSequence` | `string` | the scoring window (positions −13..+2), up to 15 residues |
| `IsLikelySignalPeptide` | `bool` | `Score ≥ minWeight` |

### 3.3 Preconditions and Validation

Null, empty, or sequences shorter than one full 15-residue window return `null`. Input is upper-cased (case-insensitive). Indexing is 0-based internally; `CleavagePosition` is reported 1-based. Residues outside the 20 standard one-letter codes are not in the matrix and contribute 0 to the score (no exception), matching EMBOSS's handling of unmapped residues [3]. The accepted alphabet is protein (no DNA/RNA translation).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate input; reject sequences shorter than 15 residues.
2. Select the eukaryotic (default) or prokaryotic log-odds weight matrix (built once from the count matrices).
3. For each candidate index `i`, sum the 15 window weights `W(residue, column)`, skipping out-of-range columns and unmapped residues.
4. Track the maximum score and its index (argmax).
5. Report the best site: 1-based mature start, score, signal sequence, window, and the `Score ≥ minWeight` flag.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

The eukaryotic count matrix (161 signal peptides) and prokaryotic count matrix (36 signal peptides) are reproduced verbatim from EMBOSS 6.6.0 `data/Esig.euk` and `data/Esig.pro` [4], each attributed to von Heijne (1986) [1]. The log-odds transform (with the −3/−1 pseudocount special case) is applied once at static initialization. The default acceptance threshold is 3.5 [2].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `PredictSignalPeptide` | O(n × 15) = O(n) | O(1) | n = sequence length; window is fixed-width |
| Matrix build (once) | O(20 × 15) | O(20 × 15) | static, amortized |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProteinMotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)

- `ProteinMotifFinder.PredictSignalPeptide(string, bool, double)`: scores all positions and returns the argmax site.
- `ProteinMotifFinder.BuildWeightMatrix(int[][], double[])`: builds the log-odds matrix from counts with EMBOSS pseudocounts (private).

### 5.2 Current Behavior

The two count matrices and their background ("Expect") vectors are stored as static literals copied from the EMBOSS data files; the log-odds matrices are computed once. Scoring uses a fixed 15-column window (`pval = −13`, 15 columns) exactly as in `sigcleave.c`. Search/matching reuse: **not applicable** — this is a fixed-width position-specific weight-matrix scan, not substring/occurrence matching, so the repository suffix tree does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- The eukaryotic (161) and prokaryotic (36) von Heijne (1986) count matrices [1][4].
- The `ln(count/expect)` log-odds transform with the `1.0e-10` penalty at columns −3/−1 and `1.0` elsewhere [3].
- Window positions −13..+2 and the argmax single-site selection [3].
- The `-minweight` default of 3.5 as the likelihood threshold [2].

**Intentionally simplified:**

- Minimum input length is one full 15-residue window; EMBOSS scores any length by skipping off-window columns. **Consequence:** sequences shorter than 15 residues return `null` rather than a truncated-window score (a partial-window score is not meaningful for cleavage prediction).

**Not implemented:**

- The hidden-Markov / neural refinements of later predictors (SignalP, Phobius); **users should rely on:** those external tools for higher accuracy. The von Heijne matrix is a classical baseline with ≈75–80% cleavage accuracy [2].

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Min length = 15 (full window) | Assumption | very short inputs return null | accepted | ASM/INV-05; in-scope signal peptides are ≥ 15 aa |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null / empty | `null` | precondition |
| < 15 residues | `null` | one full window required |
| Exactly 15 residues | non-null (one window scored) | window fits once |
| Lower-case input | identical to upper-case | input upper-cased |
| Non-standard residue (X, B, Z, *) | contributes 0; no exception | not in matrix [3] |
| `minWeight` above the best score | best site returned, `IsLikelySignalPeptide = false` | flag is independent of selection [2] |

### 6.2 Limitations

The method is a statistical heuristic (≈75–80% cleavage-site accuracy [2]); it always returns a best site even for proteins without a signal peptide (the `IsLikelySignalPeptide` flag, not the presence of a result, indicates a probable signal peptide). It does not model the n/h/c region lengths explicitly, signal anchors, or lipoprotein cleavage. The matrices are organism-class generic (eukaryote vs prokaryote), not species-specific.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
var sp = ProteinMotifFinder.PredictSignalPeptide(achDromeSequence);
// sp.Value.CleavagePosition == 42   (mature protein starts at residue 42)
// sp.Value.Score          ≈ 13.739  (von Heijne weight)
// sp.Value.IsLikelySignalPeptide == true
```

**Numerical walk-through:** For UniProt P17644 (ACH2_DROME), the maximum weight-matrix score over all positions is 13.739 at mature-protein start residue 42; the 15-residue window ending at position −1 is `…LLVLLLLCETVQA`, matching the EMBOSS `sigcleave` reference output [2]. The runner-up site (mature start 39) scores 12.135, confirming the argmax.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ProteinMotifFinder_PredictSignalPeptide_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_PredictSignalPeptide_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [PROTMOTIF-SP-001-Evidence.md](../../../docs/Evidence/PROTMOTIF-SP-001-Evidence.md)
- Related algorithms: [Domain_Prediction](../ProteinMotif/Domain_Prediction.md)

## 8. References

1. von Heijne, G. 1986. A new method for predicting signal sequence cleavage sites. Nucleic Acids Research 14(11):4683–4690. https://doi.org/10.1093/nar/14.11.4683
2. Rice, P., Longden, I., Bleasby, A. 2000. EMBOSS `sigcleave` application documentation (6.6.0). https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/sigcleave.html
3. EMBOSS 6.6.0 source, `emboss/sigcleave.c` (scoring loop and `sigcleave_readSig` matrix transform). https://raw.githubusercontent.com/lauringlab/CodonShuffle/master/lib/EMBOSS-6.6.0/emboss/sigcleave.c
4. EMBOSS 6.6.0 data files `data/Esig.euk` (161 eukaryotic) and `data/Esig.pro` (36 prokaryotic). https://raw.githubusercontent.com/lauringlab/CodonShuffle/master/lib/EMBOSS-6.6.0/emboss/data/Esig.euk
5. UniProt Consortium. ACH2_DROME (P17644). https://rest.uniprot.org/uniprotkb/P17644.fasta
