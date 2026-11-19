# Sequence Similarity (k-mer Jaccard Index)

| Field | Value |
|-------|-------|
| Algorithm Group | Analysis |
| Test Unit ID | GENOMIC-SIMILARITY-001 |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-11-19 (audit round 1: line-by-line pass WP26) |

## 1. Overview

`GenomicAnalyzer.CalculateSimilarity` measures the similarity of two DNA sequences using the
**Jaccard index of their k-mer sets**. Each sequence is reduced to the set of its distinct
length-k substrings (k-mers); the Jaccard index is the fraction of k-mers shared between the
two sets out of all distinct k-mers in either set [1][2]. By default the result is exact (the distinct k-mer
sets, not heuristic or probabilistic) and is reported as a percentage in [0, 100]. Optionally the index is
estimated from a MinHash sketch (Mash bottom-s, `sketchSize`; sourmash FracMinHash, `scaled`) and
`CompareSimilarity` returns the full set of derived quantities (Jaccard distance, Mash distance, Jaccard ANI,
containment indices and containment ANI; audit A2-31). It is an alignment-free comparison suitable for quick
resemblance estimates between sequences.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Alignment-free sequence comparison summarizes a sequence by its k-mer composition and compares
those compositions instead of aligning residues. The Jaccard index is the canonical set-resemblance
measure used for this purpose; Mash popularized k-mer Jaccard (via MinHash sketches) for genome and
metagenome distance estimation [2].

### 2.2 Core Model

For two finite sets A and B, the **Jaccard index** is [1]:

```
J(A,B) = |A ∩ B| / |A ∪ B| = |A ∩ B| / (|A| + |B| − |A ∩ B|)
```

with `0 ≤ J(A,B) ≤ 1` [1]. Here A and B are the **sets of distinct k-mers** of the two sequences.
Applied to k-mer sets, J is "the fraction of shared k-mers out of all distinct k-mers in A and B" [2].
This method reports `J × 100` as a percentage.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 ≤ result ≤ 100` | `0 ≤ J ≤ 1` [1], scaled ×100 |
| INV-02 | Identical non-empty sequences → 100.0 | A = B ⇒ A∩B = A∪B ⇒ J = 1 [1] |
| INV-03 | Disjoint k-mer sets → 0.0 | A∩B = ∅, A∪B ≠ ∅ ⇒ J = 0 [1] |
| INV-04 | Symmetric: result(a,b,k) = result(b,a,k) | ∩ and ∪ are commutative [1] |
| INV-05 | k-mers compared as sets (within-sequence repeats counted once) | Jaccard over distinct k-mers [2] |

### 2.5 Comparison with Related Methods

| Aspect | k-mer Jaccard (this) | Alignment identity |
|--------|----------------------|--------------------|
| Basis | set resemblance of k-mer composition | residue-by-residue over an alignment |
| Cost | O(n+m) | O(n·m) (dynamic programming) |
| Sensitivity to order | low (composition only) | high (positional) |
| Sketching | none by default (exact sets); optional Mash bottom-s / sourmash FracMinHash estimators [2] | n/a |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence1 | `DnaSequence` or `string` | required | first sequence | non-null; `DnaSequence`: alphabet A/C/G/T (enforced); `string`: a k-mer window with any other symbol (N, IUPAC, gap) is skipped (Mash `-n` / Jellyfish rule) |
| sequence2 | `DnaSequence` or `string` | required | second sequence | as `sequence1` |
| kmerSize | int | 5 | k-mer length | ≥ 1 (≤ 32 with a sketch) |
| canonical | bool | false | compare canonical k-mers min(w, revcomp w) (Mash/sourmash) instead of literal k-mers; not in the 3-argument `DnaSequence` overload | — |
| sketchSize | int | none | Mash bottom-s sketch size s (`CalculateSimilarity(…, canonical, sketchSize)`, `CompareSimilarity`) | ≥ 1; not with `scaled` |
| scaled | long | none | sourmash FracMinHash scaled factor S (`CompareSimilarity` only) | 1 … 4294967295; not with `sketchSize` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return of `CalculateSimilarity`) | double | k-mer Jaccard index ×100 (exact or sketch estimate), a percentage in [0, 100] |
| (return of `CalculateJaccardDistance`) | double | `1 − J` in [0, 1] |
| (return of `CompareSimilarity`) | `SimilarityEstimate` | `Estimator` (Exact / MashBottomSketch / FracMinHash), `KmerSize`, `Canonical`, `Jaccard`, `JaccardDistance`, `MashDistance`, `JaccardAni`, `ContainmentAInB`, `ContainmentBInA`, `MaxContainment`, `ContainmentAniAInB`, `ContainmentAniBInA`, `SharedHashes`, `DenominatorHashes`, `MashPValue` (all fractions in [0, 1]; containment is null for bottom-s sketches) |

### 3.3 Preconditions and Validation

- `sequence1` / `sequence2` null → `ArgumentNullException`.
- `kmerSize < 1` → `ArgumentOutOfRangeException`; an invalid `sketchSize` / `scaled` (or both given) → `ArgumentOutOfRangeException` / `ArgumentException`.
- `DnaSequence` inputs are upper-case A/C/G/T (normalized and validated); `string` inputs are case-insensitive and may contain other symbols (windows containing one are skipped).
- Empty sequences (or sequences shorter than k) produce empty k-mer sets; if both are empty the union is empty and the method returns `0.0` (see §5.4).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate arguments (non-null sequences; `kmerSize ≥ 1`; sketch parameters).
2. Build the distinct-k-mer set of each sequence (sliding window, `HashSet<string>`; canonical: min(w, revcomp w)).
3. Compute `|A ∩ B|` (count of A's k-mers present in B) and `|A ∪ B| = |A| + |B| − |A ∩ B|`.
4. Return `0` if the union is empty; otherwise `|A∩B| / |A∪B| × 100`.

With `sketchSize` / `scaled` steps 2-3 are replaced by the canonical sketch builders (§5.3): the Mash bottom-s estimate is x/s over the s smallest hashes of the union, the FracMinHash estimate is the shared fraction of the hashes below 2^64 / S.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateSimilarity (exact) | O(n + m) | O(n + m) | n,m = sequence lengths; k-mer extraction + hash-set intersection, each k-mer O(k) treated as O(1) for fixed k |
| CalculateSimilarity / CompareSimilarity (sketch) | O(n + m) hashing + the sketch comparison | O(s) bottom-s hashes, or the hashes kept by the scaled filter | MurmurHash3-x64-128 of every k-mer (canonical: of the smaller strand), then the sketch comparison |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

- `GenomicAnalyzer.CalculateSimilarity(DnaSequence, DnaSequence, int)`: literal k-mer Jaccard index ×100.
- `GenomicAnalyzer.CalculateSimilarity(DnaSequence, DnaSequence, int, bool canonical)`: same, optionally over canonical k-mers.
- `GenomicAnalyzer.CalculateSimilarity(…, int, bool, int sketchSize)`, `CalculateJaccardDistance(…)`, `CompareSimilarity(…, int?, long?)` (`DnaSequence` and `string` overloads; [GenomicAnalyzer.Similarity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.Similarity.cs), audit A2-31 / A2-35): the sketch estimators and the derived distances / ANI.
- The exact overloads are facades: validation, then `KmerAnalyzer.JaccardSimilarity(a.Sequence, b.Sequence, k, new KmerCountingOptions(Canonical: canonical)) * 100` (the canonical k-mer Jaccard, DUP_MAP §16). The former private `GetKmers` duplicate was removed (B09 F32).
- MCP `calculate_similarity` delegates to these overloads (`canonical`, `sketchSize`, `scaled`, estimates; IUPAC strings accepted).

### 5.2 Current Behavior

- k-mers are held in a distinct-k-mer set (`KmerAnalyzer.DistinctKmers`), so duplicate k-mers within a sequence collapse to one element — exactly the distinct-set semantics required by the Jaccard definition [1][2].
- Union is `|A|+|B|−|A∩B|`. Default mode = literal (strand-specific) k-mers = Python `set` / `scipy.spatial.distance.jaccard` on presence vectors = `mash -n`; `canonical: true` = Mash default / sourmash (exact, unsketched; a sequence and its reverse complement score 100). Both modes are oracle-checked on 2086 rows (scipy 1.18.1, sourmash 4.9.4).
- Case: `DnaSequence` folds to upper case, so the comparison is case-insensitive; no non-ACGT symbol can occur, so the Mash `-n`/Jellyfish ACGT-only rule holds by construction.
- **Suffix tree evaluated but not used.** This unit computes a set intersection/union of distinct k-mers, not an occurrence-search ("find all positions of X in Y") problem. A `HashSet<string>` gives O(n+m) construction and O(1) membership, which is optimal for set resemblance; the repository `SuffixTree` would add O(n) construction cost without providing the set-overlap query this metric needs. The suffix tree is therefore correctly not applied here.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `J(A,B) = |A∩B| / |A∪B|` over the distinct-k-mer sets [1].
- k-mer Jaccard as the fraction of shared k-mers out of all distinct k-mers [2].
- Range and extremes: J=1 for identical sets, J=0 for disjoint sets [1].

**Intentionally simplified:**

- Percentage reporting: result is `J × 100` rather than the raw [0,1] coefficient; **consequence:** values are 0–100, not 0–1. Relative ordering is unchanged.

**Added by audit A2-31 / A2-35 (formerly "not implemented"):**

- **MinHash / FracMinHash estimators.** `CalculateSimilarity(a, b, k, canonical, sketchSize)` and `CompareSimilarity(a, b, k, canonical, sketchSize | scaled)` estimate J from a Mash bottom-s sketch (MurmurHash3-x64-128, seed 42, `mash dist`'s x/s over the smallest s hashes of the union; k ≤ 32) or a sourmash FracMinHash sketch (`scaled = S`); the sketches and comparisons are the canonical `KmerAnalyzer.CreateMinHashSketch` / `CompareMinHashSketches` / `CreateFracMinHashSketch` / `CompareFracMinHashSketches` (B06). A sketch at least as large as the union reproduces the exact index.
- **Distances and ANI** (`SimilarityEstimate`): Jaccard distance `1 − J` (`CalculateJaccardDistance`), Mash distance `−ln(2J/(1+J))/k` (`KmerAnalyzer.MashDistance` / `MashDistanceFromJaccard`), the sourmash Jaccard ANI `(2J/(1+J))^(1/k)`, the containment indices `C(A,B) = |A∩B|/|A|`, `C(B,A)`, max-containment, and the containment ANI `C^(1/k)` (sourmash `containment_to_distance` point estimate; Hera et al. 2023).
- **String input**: a k-mer window with a symbol other than A/C/G/T is skipped (Mash `-n` / Jellyfish rule), so N runs are allowed; an A/C/G/T string equals the `DnaSequence` overload.
- Verification: 4290 oracle rows — exact canonical sets (Python) == sourmash `MinHash(scaled=1)`, sourmash 4.9.4 FracMinHash (scaled 2 … 10: jaccard, contained_by both ways, max_containment) and bottom-s MinHash (`num` 5 … 500, k = 17 … 31, where Mash and sourmash hash identically; Mash truncates the hash to 32 bits for k ≤ 16), `containment_to_distance` / `jaccard_to_distance` point estimates (`scripts/analysis/generate_similarity_sketch_oracle.py`, `GenomicAnalyzer_SimilaritySketch_Oracle_Tests`).
- Not provided: confidence intervals of the ANI estimates (sourmash `estimate_ci`), `mash dist` p-value beyond `MashPValue`.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Empty-union returns 0.0 | Assumption | Affects both-empty / both-shorter-than-k inputs | accepted | Jaccard undefined for empty union [1]; 0 chosen as "no similarity" (ASM-1 in TestSpec) |
| 2 | ×100 scaling | Assumption | Output in [0,100] not [0,1] | accepted | presentation convention only |
| 3 | Empty-union value 0 | Convention | scipy / Mash imply J = 1, sourmash 0.0 | accepted | no source fixes it; 0 = "no shared content", same as `KmerAnalyzer.JaccardSimilarity` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Identical non-empty sequences | 100.0 | J = 1 [1] |
| Disjoint k-mer sets | 0.0 | J = 0 [1] |
| Both sequences empty | 0.0 | empty union; undefined → impl returns 0 [1] |
| Both sequences shorter than k | 0.0 | both k-mer sets empty → empty union |
| One sequence empty, other non-empty | 0.0 | empty intersection over non-empty union |
| Within-sequence repeated k-mers | counted once | distinct-set semantics [2] |
| Null sequence | ArgumentNullException | input validation |
| kmerSize < 1 | ArgumentOutOfRangeException | input validation |

### 6.2 Limitations

- Composition-only: insensitive to k-mer order/position, so two sequences with the same k-mer composition but different arrangement score identically. For positional/edit comparison use alignment.
- Exact-match k-mers only — no mismatch tolerance; a single substitution alters up to k k-mers. The sketch estimators are estimates (the error shrinks with a larger s / a smaller S); the Mash sketch needs k ≤ 32 and there is no confidence interval or p-value beyond Mash's `MashPValue`.
- Default is strand-specific: a sequence and its reverse complement are generally dissimilar; pass `canonical: true` for strand-agnostic comparison.
- Choice of k (default 5) controls resolution: larger k is more specific but more sensitive to errors; the formula is unchanged.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var a = new DnaSequence("ACGTACGT");
var b = new DnaSequence("ACGTACGA");
double pct = GenomicAnalyzer.CalculateSimilarity(a, b, kmerSize: 3); // 80.0
```

**Numerical walk-through (k=3, ACGTACGT vs ACGTACGA):**

- A = distinct 3-mers of ACGTACGT = {ACG, CGT, GTA, TAC} (|A| = 4; ACG and CGT repeat but collapse).
- B = distinct 3-mers of ACGTACGA = {ACG, CGT, GTA, TAC, CGA} (|B| = 5).
- A ∩ B = {ACG, CGT, GTA, TAC} → |A∩B| = 4.
- A ∪ B = {ACG, CGT, GTA, TAC, CGA} → |A∪B| = 5.
- J = 4/5 = 0.8 → ×100 = **80.0**.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomicAnalyzer_CalculateSimilarity_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_CalculateSimilarity_Tests.cs) — covers INV-01..INV-05; the sketch estimators, distances and ANI: `GenomicAnalyzer_SimilaritySketch_Oracle_Tests` (4290 sourmash / Mash rows, `scripts/analysis/generate_similarity_sketch_oracle.py`)
- Evidence: [GENOMIC-SIMILARITY-001-Evidence.md](../../../docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md)

## 8. References

1. Jaccard, Paul. 1901. Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société vaudoise des sciences naturelles 37(142):547–579. https://en.wikipedia.org/wiki/Jaccard_index
2. Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM. 2016. Mash: fast genome and metagenome distance estimation using MinHash. Genome Biology 17:132. https://doi.org/10.1186/s13059-016-0997-x
3. Mash documentation — Distance Estimation. marbl/Mash. https://github.com/marbl/Mash/blob/master/doc/sphinx/distances.rst
