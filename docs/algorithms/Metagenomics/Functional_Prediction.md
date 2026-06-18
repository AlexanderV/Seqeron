# Functional Prediction (Homology Transfer + Pathway Enrichment)

| Field | Value |
|-------|-------|
| Algorithm Group | Metagenomics |
| Test Unit ID | META-FUNC-001 |
| Related Projects | Seqeron.Genomics.Metagenomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Functional prediction assigns biological function to predicted genes/proteins in a metagenome and then asks which biological pathways are statistically over-represented among them. This unit implements two specification-driven, deterministic pieces: (1) `PredictFunctions`, which transfers an annotation (function, pathway, KO) from a signature database to a query protein and reports the homology significance as a BLAST bit score and E-value computed from Karlin-Altschul statistics [1][2][3]; and (2) `FindPathwayEnrichment`, which scores each pathway with the hypergeometric over-representation test [4]. The numerical core (bit score, E-value, hypergeometric p-value) is exact with respect to the cited formulas; the matching step is a simplified exact-signature model (see §5.3).

## 2. Scientific / Formal Basis

> A = Homology-based annotation transfer (BLAST statistics), B = Pathway over-representation (hypergeometric)

### 2.1 Domain Context

Homology-based annotation transfer is the standard way to annotate uncharacterized genes: a query is compared to characterized sequences, and the function of the best significant hit is transferred [1]. The significance of a local alignment of score S is quantified by the E-value — the expected number of equally good or better alignments arising by chance in a database of the given size [1]. Over-representation analysis (ORA) then asks whether the annotated query set contains more members of a pathway than expected from random sampling of the background gene universe [4].

### 2.A Homology-based annotation transfer

#### Core Model

For an ungapped alignment of raw score S the expected number of high-scoring segment pairs (HSPs) with score at least S is `E = K·m·n·e^(−λ·S)`, where m and n are the lengths of the two sequences and K, λ are statistical parameters of the scoring system [1]. The normalized **bit score** is `S' = (λ·S − ln K) / ln 2`, and the E-value corresponding to a bit score is `E = m·n·2^(−S')` [1]. For the BLOSUM62 scoring system, the ungapped parameters are `λ = 0.3176`, `K = 0.134` (and `H = 0.4012`), from NCBI's `blast_stat.c` `blosum62_values` table [2]. The raw self-alignment score of a matched segment is the sum of the BLOSUM62 diagonal scores over its residues (e.g. W = 11, C = 9, A = 4) [3].

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Bit score S' is strictly increasing in raw score S | S' is linear in S with slope λ/ln2 > 0 [1] |
| INV-02 | `K·m·n·e^(−λ·S) = m·n·2^(−S')` | Algebraic identity of the two cited forms [1] |
| INV-03 | E-value is strictly decreasing in S | E ∝ e^(−λ·S), λ > 0 [1] |
| INV-06 | The annotation transferred is the matching DB entry with the lowest E-value | Best-hit rule; E-value ranks significance [1] |

### 2.B Pathway over-representation (hypergeometric)

#### Core Model

Given a background universe of N genes containing M members of a pathway, and a query sample of n genes of which x fall in the pathway, the right-tail over-representation p-value is `P(X ≥ x) = 1 − Σ_{i=0}^{x−1} C(M,i)·C(N−M, n−i) / C(N, n)` [4]. X follows the hypergeometric distribution (sampling n without replacement from N) [4].

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-04 | p-value ∈ [0, 1] | It is a probability; result is clamped to [0,1] [4] |
| INV-05 | p-value = 1 when x = 0, or M = 0, or n = 0 | Empty sum / degenerate margins ⇒ no over-representation [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| proteins | `IEnumerable<(string GeneId, string ProteinSequence)>` | required | Query genes (single-letter amino acids) | non-null; empty/whitespace sequences are skipped |
| functionDatabase | `IReadOnlyDictionary<string,(string Function,string Pathway,string Ko)>` | required | signature → annotation | non-null |
| queryGenes | `IEnumerable<string>` | required | gene set of interest | non-null |
| pathwayDatabase | `IReadOnlyDictionary<string, IReadOnlyCollection<string>>` | required | pathway id → member genes | non-null |
| backgroundGenes | `IEnumerable<string>?` | null | background universe; null/empty ⇒ union of all pathway members | — |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `FunctionalAnnotation` | record struct | Best-hit annotation per matched gene: Function, Pathway, KoNumber, CogCategory, EValue, BitScore |
| `PathwayEnrichment` | record struct | Per-pathway Overlap, PathwaySize, QuerySize, BackgroundSize, PValue (right-tail) |

### 3.3 Preconditions and Validation

Null `proteins`/`functionDatabase`/`queryGenes`/`pathwayDatabase` throw `ArgumentNullException`. Empty or whitespace protein sequences and empty signatures are skipped (no annotation). Signature matching is case-sensitive ordinal; BLOSUM62 lookup upper-cases residues and treats unknown residues as score 0. Genes that match no signature yield no annotation. The query is always included in the background universe.

## 4. Algorithm

> A = Homology transfer, B = Pathway enrichment

### 4.A Homology transfer

#### High-Level Steps

1. For each gene with a non-empty protein sequence, scan every database signature.
2. If the signature occurs exactly in the protein, compute the raw BLOSUM62 self-score S of the signature, then S' = (λS − ln K)/ln2 and E = K·m·n·e^(−λS) (m = protein length, n = signature length).
3. Keep the candidate with the lowest E-value (best hit) and emit one `FunctionalAnnotation` per matched gene.

### 4.B Pathway enrichment

#### High-Level Steps

1. Build the query set and the background universe (explicit, else union of all pathway members; the query is always added).
2. For each pathway, intersect its members with the background to get M, and count the overlap x with the query (n = query size, N = background size).
3. Compute P(X ≥ x) in log-space (log-Gamma); sort pathways ascending by p-value.

#### Decision Rules / Reference Tables

- Ungapped BLOSUM62 Karlin-Altschul parameters: λ = 0.3176, K = 0.134 [2].
- BLOSUM62 diagonal scores: A4 R5 N6 D6 C9 Q5 E5 G6 H8 I4 L4 K5 M5 F6 P7 S4 T5 W11 Y7 V4 [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `PredictFunctions` | O(g × d × L) | O(1) per hit | g genes × d database signatures × substring scan cost L; matches Registry O(n×g) [META-FUNC-001] |
| `FindPathwayEnrichment` | O(P·M̄ + Σ_p x_p) | O(N) | P pathways; hypergeometric tail sums up to x_p terms; log-Gamma O(1) per term |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MetagenomicsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs)

- `MetagenomicsAnalyzer.PredictFunctions(...)`: best-hit homology annotation transfer with BLAST bit score / E-value.
- `MetagenomicsAnalyzer.FindPathwayEnrichment(...)`: hypergeometric pathway over-representation.
- `MetagenomicsAnalyzer.FunctionalBitScore(rawScore)`, `ExpectedValue(rawScore, m, n)`, `Blosum62SelfScore(segment)`, `HypergeometricUpperTail(x, N, M, n)`: source-backed numeric helpers.

### 5.2 Current Behavior

Matching uses an exact, ordinal `string.Contains` of the database signature within the query protein. The repository suffix tree was **evaluated and not used**: each gene is scanned against the database once (a single short search per gene×signature, not many queries against one fixed text), and the operative cost is the per-residue BLOSUM62 scoring and the E-value computation, not occurrence enumeration — so the suffix tree's O(n) build + O(m) query advantage does not apply here. The hypergeometric tail is summed in log-space via a Lanczos log-Gamma to avoid factorial overflow for large backgrounds (validated to N = 8000).

### 5.3 Conformance to Theory / Spec

#### 5.3.A Homology transfer

**Implemented (verbatim from the cited theory/spec):**

- Bit score `S' = (λS − ln K)/ln 2` and E-value `E = K·m·n·e^(−λS) = m·n·2^(−S')` with ungapped BLOSUM62 λ = 0.3176, K = 0.134 [1][2].
- BLOSUM62 diagonal self-match scoring of the matched segment [3].
- Best-hit (lowest E-value) annotation transfer [1].

**Intentionally simplified:**

- Matching: exact signature occurrence rather than gapped local alignment; **consequence:** only verbatim sub-sequence hits are found; divergent homologs that a full Smith-Waterman/BLAST search would detect are missed (ASM-01).
- CogCategory: inferred by keyword from the function string (a labeling convenience, not a scoring input); **consequence:** the COG letter is heuristic and does not affect bit score, E-value, or which hit is selected.

**Not implemented:**

- Gapped alignment, gapped λ/K, composition-based statistics, multiple-testing correction (e.g. Benjamini-Hochberg) of pathway p-values; **users should rely on:** external tools (NCBI BLAST+, DIAMOND, clusterProfiler) for production-scale annotation and FDR control.

#### 5.3.B Pathway enrichment

**Implemented (verbatim from the cited theory/spec):**

- Hypergeometric right-tail p-value `P(X ≥ x) = 1 − Σ_{i=0}^{x−1} C(M,i)·C(N−M,n−i)/C(N,n)` [4].
- Degenerate-margin / empty-sum handling (p = 1) [4].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Multiple-testing correction across pathways; **users should rely on:** applying BH/Bonferroni to the returned p-values.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Exact-signature matching (ASM-01) | Assumption | Affects which hits are found, not the bit-score/E-value formulas | accepted | See §5.3.A; full alignment is out of scope for this unit |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty/whitespace protein sequence | No annotation for that gene | Nothing to match |
| Gene matches no signature | No annotation emitted | Homology transfer needs a hit |
| Empty or null function database | No annotations | No signatures to match |
| Null `proteins`/`functionDatabase`/`queryGenes`/`pathwayDatabase` | `ArgumentNullException` | Input validation |
| Pathway overlap x = 0, or M = 0, or query n = 0 | p-value = 1.0 | Empty sum / degenerate margins [4] |

### 6.2 Limitations

Exact-substring matching does not detect divergent homologs; there is no gapped alignment, no composition-adjusted statistics, no FDR control, and no protein-domain (HMM/Pfam) modeling. The bit-score/E-value model assumes the ungapped BLOSUM62 scoring system; using other matrices would require their own λ/K.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**Numerical walk-through (homology transfer):** A signature "WWW" matched in a 3-residue protein. Raw score S = 3 × 11 = 33 (BLOSUM62 W diagonal [3]). Bit score S' = (0.3176·33 − ln 0.134)/ln 2 = 18.0202932787533. E-value (m = n = 3) = 0.134·3·3·e^(−0.3176·33) = 3.3852730346546 × 10⁻⁵, equal to 3·3·2^(−18.0202932787533) [1].

**Numerical walk-through (enrichment):** N = 8000, M = 400, n = 100, x = 20 ⇒ P(X ≥ 20) = 7.88 × 10⁻⁸ [4].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MetagenomicsAnalyzer_PredictFunctions_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_PredictFunctions_Tests.cs) — covers INV-01..INV-06.
- Evidence: [META-FUNC-001-Evidence.md](../../../docs/Evidence/META-FUNC-001-Evidence.md)

## 8. References

1. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. 1990. Basic local alignment search tool. J Mol Biol 215(3):403–410. NCBI BLAST tutorial "The Statistics of Sequence Similarity Scores". https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
2. NCBI C BLAST Toolkit. blast_stat.c (`blosum62_values` ungapped Karlin-Altschul parameters). https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/algo/blast/core/blast_stat.c
3. Henikoff S, Henikoff JG. 1992. Amino acid substitution matrices from protein blocks (BLOSUM62). NCBI matrix file. https://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62
4. PNNL Computational Mass Spectrometry. Proteomics Data Analysis in R/Bioconductor, §8.2 Over-Representation Analysis. https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html
