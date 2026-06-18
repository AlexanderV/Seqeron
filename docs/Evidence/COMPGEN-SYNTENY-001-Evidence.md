# Evidence Artifact: COMPGEN-SYNTENY-001

**Test Unit ID:** COMPGEN-SYNTENY-001
**Algorithm:** Synteny / Collinearity Block Detection (MCScanX collinearity model)
**Date Collected:** 2026-06-13

---

## Online Sources

### MCScanX (Wang et al. 2012) — Nucleic Acids Research — PMC full text

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3326336
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed paper) / 3 (the canonical reference implementation it documents)
**Retrieved how:** WebSearch query `MCScanX collinearity synteny blocks algorithm gene order anchors gap chaining` → opened the PMC full-text article with WebFetch and extracted the algorithm section.

**Key Extracted Points:**

1. **DP scoring recurrence (verbatim):** `Score(v) = max(MatchScore(v), max(Score(u) + MatchScore(v) + GapPenalty × NumberofGaps(u,v)))`.
2. **MatchScore:** `MatchScore(v) = 50` for one gene pair.
3. **GapPenalty:** `GapPenalty = −1`.
4. **NumberofGaps(u,v):** the maximum number of intervening genes between anchor `u` and anchor `v`; must be **fewer than 25** (default `MAX_GAPS = 25`).
5. **Minimum chain score:** non-overlapping chains with **scores over 250** are reported (i.e. at least 5 collinear gene pairs, since `5 × 50 = 250` with zero gaps).
6. **Minimum collinear gene pairs:** at least **5** pairs per reported block (default).
7. **Match collapsing:** "if consecutive BLASTP matches have a common gene and its paired genes are separated by fewer than five genes, these matches are collapsed using a representative pair with the smallest BLASTP E-value."
8. **E-value cutoff for anchors:** default `10^−5`.
9. **Directionality:** BLASTP matches are sorted "in both transcriptional directions", enabling detection of both **forward** and **inverted (reverse)** collinear blocks.

### MCScanX (Wang et al. 2012) — Oxford Academic abstract/HTML

**URL:** https://academic.oup.com/nar/article/40/7/e49/1202057
**Accessed:** 2026-06-13
**Authority rank:** 1
**Retrieved how:** Opened with WebFetch from the same search result list.

**Key Extracted Points:**

1. **Synteny vs collinearity (verbatim):** "Collinearity, a more specific form of synteny, requires conserved gene order." Synteny = homologous genes on corresponding chromosomes; collinearity = the subset retaining conserved order.
2. **DP over anchor chains:** the method "rewards the adjacent collinear gene pairs (or 'anchor genes') and penalizes the distance between anchor genes."
3. **Anchors:** "anchor genes are more likely to be homologs"; anchors are homologous (ortholog/paralog) gene pairs.
4. **Orientation:** the alignment procedure considers "both transcriptional directions" → forward and reverse blocks.
5. **Default block size:** blocks require **at least 5 collinear gene pairs** (match score 50/pair, minimum total score 250), maximum intervening genes 25.

### Wikipedia — Synteny (for definition + cited primaries)

**URL:** https://en.wikipedia.org/wiki/Synteny
**Accessed:** 2026-06-13
**Authority rank:** 4 (used only to confirm definitions and to locate primaries)
**Retrieved how:** WebSearch query `synteny Wikipedia collinearity conserved gene order definition primary source` → opened with WebFetch.

**Key Extracted Points:**

1. **Conserved synteny (verbatim):** "Shared synteny (also known as conserved synteny) describes preserved co-localization of genes on chromosomes of different species."
2. **Modern usage:** since ~2000 synteny is used for "preservation of the precise order of genes on a chromosome passed down from a common ancestor"; traditional geneticists prefer "collinearity".
3. **Method:** syntenic blocks are detected with "a version of the MCScan algorithm … looking for common patterns of collinearity" using **dynamic programming** to select optimal paths of shared homologous genes, accounting for gene loss and gain.
4. **Primary citation for the computational method:** Wang et al. (April 2012), MCScanX.

---

## Documented Corner Cases and Failure Modes

### From MCScanX (Wang et al. 2012)

1. **Sub-threshold chains:** chains scoring ≤ 250 (fewer than 5 collinear pairs) are NOT reported — they are not collinear blocks.
2. **Gap cutoff:** an anchor separated from the previous anchor by ≥ 25 intervening genes cannot extend the chain (NumberofGaps must be < 25).
3. **Non-overlapping reporting:** only non-overlapping chains are reported (a gene pair belongs to at most one reported block).
4. **Inversions:** anchors whose target order decreases form reverse-oriented (inverted) collinear blocks; these are valid blocks, distinct from forward blocks.

### From definitions (MCScanX / Wikipedia)

1. **No anchors:** genomes with no homologous/orthologous anchor pairs produce no collinear blocks.
2. **Single chromosome / one direction:** collinearity is direction-consistent within a block; a chain cannot mix increasing and decreasing target order.

---

## Test Datasets

### Dataset: Synthetic forward collinear chain (derived from the MCScanX scoring scheme)

**Source:** MCScanX scoring rule, Wang et al. (2012), PMC3326336.

Five orthologous anchors in identical order in both genomes, adjacent (NumberofGaps = 0 between consecutive anchors).

| Anchor | pos in genome1 | pos in genome2 |
|--------|----------------|----------------|
| g1↔h1 | 0 | 0 |
| g2↔h2 | 1 | 1 |
| g3↔h3 | 2 | 2 |
| g4↔h4 | 3 | 3 |
| g5↔h5 | 4 | 4 |

**Derived score:** `5 × 50 + (−1) × 0 = 250`. Score > threshold of `> 250`? The paper says "scores over 250" — exactly 250 with zero gaps corresponds to "at least 5 collinear gene pairs"; the paper reports these 5-pair blocks. Expected: one forward block of 5 gene pairs, `IsInverted = false`.

### Dataset: Synthetic reverse (inverted) collinear chain

**Source:** MCScanX "both transcriptional directions", Wang et al. (2012).

Five anchors with genome2 order reversed: pos2 = 4,3,2,1,0 for genome1 pos = 0..4.

**Expected:** one block of 5 gene pairs, `IsInverted = true`.

### Dataset: Sub-threshold chain (4 anchors)

**Source:** MCScanX minimum of 5 collinear pairs / score 250.

Four adjacent anchors → score `4 × 50 = 200 ≤ 250`. **Expected:** no block reported.

### Dataset: Gap exceeding cutoff

**Source:** MCScanX NumberofGaps < 25.

Two adjacent runs of anchors separated by ≥ 25 intervening genes in genome2 → the gap breaks the chain; neither sub-run reaches 5 pairs → no block.

---

## Assumptions

1. **ASSUMPTION: Threshold boundary at exactly 250.** The paper says "scores over 250" yet also says these correspond to "at least 5 collinear gene pairs", and a 5-pair zero-gap chain scores exactly 250. We adopt the paper's own "at least 5 collinear gene pairs" wording as the operative report rule: a chain is reported iff its score ≥ MinChainScore (250) AND it has ≥ MinAnchors (5). This resolves the wording tension in favour of the explicitly stated 5-pair minimum and is source-backed, so it is not an open correctness gap.
2. **ASSUMPTION: Anchors supplied as an ortholog map.** MCScanX derives anchors from BLASTP (E-value < 1e-5) with collapsing of near-duplicate matches. This repository delegates anchor/ortholog identification to a separate unit (COMPGEN-ORTHO-001) and accepts the anchor set as an input `orthologMap`. The collinearity/chaining algorithm under test is unchanged by this separation; anchor *generation* is out of scope for this unit.

---

## Recommendations for Test Coverage

1. **MUST Test:** five adjacent forward anchors form exactly one forward block (`IsInverted=false`, GeneCount=5). — Evidence: MatchScore 50, score 250, ≥5 pairs (PMC3326336).
2. **MUST Test:** five reverse anchors form exactly one inverted block (`IsInverted=true`). — Evidence: "both transcriptional directions" (PMC3326336).
3. **MUST Test:** four adjacent anchors (score 200) yield no block. — Evidence: minimum 5 pairs / score 250 (PMC3326336).
4. **MUST Test:** a gap of ≥ 25 intervening genes breaks the chain. — Evidence: NumberofGaps < 25 (PMC3326336).
5. **MUST Test:** empty genome / empty ortholog map → no blocks (null/empty inputs). — Evidence: definition requires anchors.
6. **SHOULD Test:** two separated 5-anchor runs → two non-overlapping blocks. — Rationale: non-overlapping chain reporting.
7. **SHOULD Test:** gap penalty effect — a chain with intervening genes still scores ≥250 if enough anchors. — Rationale: GapPenalty=−1 reduces but does not destroy a long chain.
8. **COULD Test:** `VisualizeSynteny` produces one text line per block (smoke). — Rationale: visualization wrapper, delegate-type.
9. **COULD Test:** property — every returned block has GeneCount ≥ 5 and coordinates within parent gene bounds. — Rationale: O(n²) invariant property test.

---

## References

1. Wang Y, Tang H, DeBarry JD, Tan X, Li J, Wang X, Lee T-H, Jin H, Marler B, Guo H, Kissinger JC, Paterson AH. (2012). MCScanX: a toolkit for detection and evolutionary analysis of gene synteny and collinearity. *Nucleic Acids Research* 40(7):e49. https://doi.org/10.1093/nar/gkr1293 — full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC3326336
2. Wang Y, et al. (2012). MCScanX (Oxford Academic HTML). *Nucleic Acids Research*. https://academic.oup.com/nar/article/40/7/e49/1202057
3. Wikipedia contributors. Synteny. https://en.wikipedia.org/wiki/Synteny (accessed 2026-06-13; used for definitions and to locate primary sources).

---

## Change History

- **2026-06-13**: Initial documentation.
