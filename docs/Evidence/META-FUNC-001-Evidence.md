# Evidence Artifact: META-FUNC-001

**Test Unit ID:** META-FUNC-001
**Algorithm:** Functional Prediction (homology-based annotation transfer + pathway over-representation)
**Date Collected:** 2026-06-13

---

## Online Sources

### Altschul et al. BLAST Tutorial — "The Statistics of Sequence Similarity Scores" (NCBI)

**URL:** https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
**Accessed:** 2026-06-13
**Authority rank:** 2 (official NCBI specification / canonical BLAST statistics tutorial)
**Retrieved by:** WebFetch of the URL above; facts extracted from the returned page text.

**Key Extracted Points:**

1. **E-value formula:** "the expected number of HSPs with score at least S is given by the formula" `E = K m n e^(-λS)`, where E is the expected number of high-scoring segment pairs (HSPs) with score ≥ S; K and λ are statistical parameters characterizing the scoring system; m and n are the lengths of the two sequences being compared; S is the raw alignment score.
2. **Bit score (normalized score):** "By normalizing a raw score using the formula" `S' = (λS − ln K) / ln 2`, where S' is the bit score in standard units, λ and K are the same parameters, S is the raw score.
3. **E-value from bit score:** "The E-value corresponding to a given bit score is simply" `E = m n 2^(−S')`. This shows that bit scores relate directly to E-values without needing K and λ individually.

### NCBI `blast_stat.c` — ungapped Karlin-Altschul parameters for BLOSUM62

**URL:** https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/algo/blast/core/blast_stat.c
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation — NCBI C BLAST toolkit source)
**Retrieved by:** WebFetch of the URL above; values extracted from the `blosum62_values` array text.

**Key Extracted Points:**

1. **Ungapped BLOSUM62 parameters:** the first row of `blosum62_values` (gap open = gap extend = INT2_MAX, i.e. the ungapped case) holds `λ = 0.3176`, `K = 0.134`, `H = 0.4012` (the 4th, 5th and 6th array elements).

### NCBI BLAST matrices — BLOSUM62 substitution matrix

**URL:** https://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62
**Accessed:** 2026-06-13
**Authority rank:** 2 (official NCBI data file)
**Retrieved by:** WebFetch of the URL above; diagonal (self-match) scores extracted from the matrix text.

**Key Extracted Points:**

1. **BLOSUM62 diagonal (self-match) scores:** A 4, R 5, N 6, D 6, C 9, Q 5, E 5, G 6, H 8, I 4, L 4, K 5, M 5, F 6, P 7, S 4, T 5, W 11, Y 7, V 4. A self-alignment of a protein segment scores the sum of these diagonal entries over its residues.

### Over-Representation Analysis (ORA) — hypergeometric test (PNNL Proteomics Data Analysis tutorial)

**URL:** https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html
**Accessed:** 2026-06-13
**Authority rank:** 4 (course material citing the hypergeometric test; the underlying distribution is standard, see Wikipedia/Fisher references)
**Retrieved by:** WebFetch of the URL above; formula and worked example extracted from the page text.

**Key Extracted Points:**

1. **Enrichment p-value (right tail):** `P(X ≥ x) = 1 − P(X ≤ x−1) = 1 − Σ_{i=0}^{x−1} [ C(M,i) · C(N−M, n−i) ] / C(N, n)`.
2. **Parameter definitions:** N = total number of background genes; n = number of "interesting" genes (the query set); M = number of genes annotated to the gene set/pathway S; x = number of interesting genes annotated to S.
3. **Worked example:** N = 8000 background genes, M = 400 annotated to S, n = 100 interesting genes, x = 20 interesting genes annotated to S → `P(X ≥ 20) = 7.88 × 10⁻⁸`. R equivalent: `phyper(q = 20 − 1, m = 400, n = 8000 − 400, k = 100, lower.tail = FALSE)`.

---

## Documented Corner Cases and Failure Modes

### From Altschul BLAST tutorial / Karlin-Altschul theory

1. **Score must be positive for the EVD to apply:** Karlin-Altschul statistics require a scoring system with positive expected maximal score (i.e. an HSP with S > 0); an empty or non-matching query yields no HSP. For an exact self-match of a non-empty segment over BLOSUM62 the raw score is strictly positive (all diagonal entries ≥ 4), so S > 0 always holds for the matched region.

### From ORA hypergeometric test

1. **No interesting genes annotated (x = 0):** `P(X ≥ 0) = 1` (the empty sum is 0, so the p-value is 1 − 0 = 1). A pathway with no overlap is never enriched.
2. **Degenerate margins:** if the pathway has M = 0 annotated genes, or the query set is empty (n = 0), no over-representation is possible and the p-value is 1.
3. **Right-tail only:** ORA tests over-representation, so only the upper tail P(X ≥ x) is reported; under-representation is a separate (lower-tail) test.

---

## Test Datasets

### Dataset: ORA hypergeometric worked example

**Source:** PNNL Proteomics Data Analysis tutorial §8.2 (Over-Representation Analysis), URL above.

| Parameter | Value |
|-----------|-------|
| N (background genes) | 8000 |
| M (genes annotated to set S) | 400 |
| n (interesting genes) | 100 |
| x (interesting genes in S) | 20 |
| P(X ≥ 20) | 7.88 × 10⁻⁸ |

### Dataset: BLAST bit-score / E-value self-consistency (derived from the two cited formulas)

**Source:** Altschul BLAST tutorial formulas + NCBI BLOSUM62 ungapped parameters (λ = 0.3176, K = 0.134) and BLOSUM62 diagonals (W = 11).

| Parameter | Value |
|-----------|-------|
| Matched segment | "WWW" (3 × W, diagonal score 11) |
| Raw score S | 33 |
| Bit score S' = (λS − ln K)/ln 2 | 18.0202932787533… |
| E = K·m·n·e^(−λS), m = n = 3 | 3.3852730346546 × 10⁻⁵ |
| E = m·n·2^(−S'), m = n = 3 | 3.3852730346546 × 10⁻⁵ |

Both E-value expressions agree to machine precision, confirming the algebraic identity `K·m·n·e^(−λS) = m·n·2^(−S')` for the cited parameters.

---

## Assumptions

1. **ASSUMPTION: Ungapped exact-match scoring model.** `PredictFunctions` transfers function from a database reference by detecting an exact occurrence of the reference signature in the query protein and scoring that ungapped self-match with BLOSUM62. The cited sources define the bit-score/E-value formulas for general HSPs; restricting to ungapped exact matches is a deliberate simplification so that λ, K (ungapped BLOSUM62, source-backed) apply directly and no gap model is needed. Justification: the unit's canonical complexity is O(n × g) (per-gene × per-database-entry scan), consistent with signature matching rather than full alignment. This affects which hits are found, not the bit-score/E-value formulas, which are used exactly as cited.

---

## Recommendations for Test Coverage

1. **MUST Test:** Bit score `S' = (λS − ln K)/ln 2` for an exact BLOSUM62 self-match (e.g. "WWW" → S = 33 → S' = 18.0202932787533) — Evidence: Altschul tutorial bit-score formula + NCBI λ, K, W diagonal.
2. **MUST Test:** E-value `E = m·n·2^(−S') = K·m·n·e^(−λS)` equals 3.3852730346546 × 10⁻⁵ for the "WWW" hit (m = n = 3) — Evidence: Altschul tutorial E-value formulas.
3. **MUST Test:** Hypergeometric ORA p-value for N=8000, M=400, n=100, x=20 equals 7.88 × 10⁻⁸ — Evidence: PNNL ORA worked example.
4. **MUST Test:** Best-hit selection — when multiple database entries match one gene, the annotation transferred is the one with the lowest E-value (highest bit score) — Evidence: BLAST best-hit / homology transfer convention (E-value ranks significance).
5. **SHOULD Test:** ORA p-value = 1 when x = 0 or M = 0 or n = 0 (no possible over-representation) — Rationale: empty-sum / degenerate-margin corner cases from ORA.
6. **SHOULD Test:** Null/empty inputs (null gene list, empty database, empty protein sequence) produce no annotations / no enrichment without throwing — Rationale: robustness, mirrors sibling analyzer methods.
7. **COULD Test:** Monotonicity invariant — larger raw score ⇒ larger bit score ⇒ smaller E-value — Rationale: follows directly from the strictly monotone formulas.

---

## References

1. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990). Basic local alignment search tool. J Mol Biol 215(3):403–410. NCBI BLAST tutorial "The Statistics of Sequence Similarity Scores": https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
2. NCBI C++ BLAST Toolkit. `blast_stat.c` (ungapped Karlin-Altschul parameter tables, `blosum62_values`). https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/algo/blast/core/blast_stat.c
3. Henikoff S, Henikoff JG (1992). Amino acid substitution matrices from protein blocks (BLOSUM62). NCBI matrix file: https://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62
4. PNNL Computational Mass Spectrometry. Proteomics Data Analysis in R/Bioconductor, §8.2 Over-Representation Analysis (hypergeometric test formula and worked example). https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html

---

## Change History

- **2026-06-13**: Initial documentation.
