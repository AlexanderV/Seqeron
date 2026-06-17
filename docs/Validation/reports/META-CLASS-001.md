# Validation Report: META-CLASS-001 — Kraken k-mer / LCA taxonomic classifier

- **Validated:** 2026-06-17   **Area:** Metagenomics
- **Canonical method(s):** `TaxonomyTree` (`Lca`, `GetPathToRoot`, `IsAncestorOf`, `GetDepth`),
  `MetagenomicsAnalyzer.BuildKmerDatabase`, `MetagenomicsAnalyzer.ClassifyReads`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN (no defect; tests independently re-derived & mutation-checked; no code/test change required)

This is an **independent re-validation** of the Kraken-style rewrite committed by a different
(implementer) session at `9839268`. Per the governing rule, every expected classification, LCA, RTL
path score, and C/Q confidence below was **hand-derived from the Kraken algorithm definition applied
to the test's taxonomy** (a from-scratch Python reference written only from the paper rules), then the
code was checked against that derivation — never the reverse.

## Stage A — Description (primary source)

**Primary source retrieved this session (open-access PMC mirror, WebFetched):**
Wood DE, Salzberg SL. "Kraken: ultrafast metagenomic sequence classification using exact alignments."
*Genome Biology* 2014, **15:R46**. doi:10.1186/gb-2014-15-3-r46.
- PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC4053813/ (Springer/BMC copy is paywalled behind IDP).

**C/Q confidence** (Kraken 2 manual, WebFetched):
https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown — confidence scoring section.

### Verbatim rules confirmed

1. **DB build (LCA collapse).** "For each sequence, the taxon associated with it is used to set the
   stored LCA values of all k-mers in the sequence. As sequences are processed, **if a k-mer from a
   sequence has had its LCA value previously set, then the LCA of the stored value and the current
   sequence's taxon is calculated** and that LCA is stored for the k-mer."
2. **Classification (RTL = sum of weights).** "Each node in the classification tree is weighted with
   the number of k-mers in K(S) that mapped to the taxon associated with that node. Then, **each
   root-to-leaf (RTL) path in the classification tree is scored by calculating the sum of all node
   weights along the path.**"
3. **Max path + tie-break.** "**The maximum scoring RTL path in the classification tree is the
   classification path, and S is assigned the label corresponding to its leaf** (if there are multiple
   maximally scoring paths, **the LCA of all those paths' leaves is selected**)."
4. **Unclassified.** "Sequences for which none of the k-mers in K(S) are found in any genome are left
   unclassified by this algorithm." (Reported here as the root taxon `TaxonomyTree.RootId`.)
5. **Confidence C/Q** (Kraken 2). "A sequence label's score is a fraction **C/Q, where C is the number
   of k-mers mapped to LCA values in the clade rooted at the label, and Q is the number of k-mers in
   the sequence that lack an ambiguous nucleotide** (i.e., they were queried against the database)."

The implementer's in-code algorithm comment (`MetagenomicsAnalyzer.cs` lines 151-170) quotes these
rules; each quote was checked against the fetched source text and matches. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/TaxonomyTree.cs` — `Lca(a,b)` (lines 163-192)
  collects a's ancestor set then walks b upward to the first hit (root if disjoint); `Lca(set)`
  folds pairwise; `GetPathToRoot`/`IsAncestorOf`/`GetDepth` as expected.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs`
  - `ClassifyRead` (lines 211-292): per-taxon hit counts + Q (non-ambiguous canonical k-mers);
    leaves = hit taxa that are **not** a proper ancestor of another hit taxon; each leaf scored by
    summing hit-weights over its **root path** (so ancestor weights are added to every descendant
    path — the true RTL sum, not deepest-hit); single max wins, tie → `Lca(bestLeaves)`; C = sum of
    hit-weights in the clade rooted at the call; confidence = C/Q.
  - `BuildKmerDatabase` (lines 363-397): first reference sets a canonical k-mer's taxon, each later
    reference replaces it with `Lca(existing, taxonId)`.
  - Canonical k-mer (lines 407-411): `min(kmer, revcomp)` by ordinal; ambiguous (non-ACGT) skipped.

### Independent hand/Python derivation vs the rewritten tests (all match)

Reference Python implementation written **only from the paper rules** on the test taxonomy
(`root(1) → {Archaea(5); Bacteria(2)→Proteobacteria(3)→Gammaproteobacteria(4)→Enterobacteriaceae(10)
→ {Escherichia(20)→{E.coli(100), E.fergusonii(101)}, Salmonella(21)→S.enterica(200)}}`):

**`TaxonomyTree.Lca` primitive** — `Lca(100,101)=20` (siblings→genus), `Lca(100,200)=10`
(cross-genus→family), `Lca(20,100)=20` (ancestor/descendant→ancestor), `Lca(100,100)=100` (self),
`Lca(5,100)=1` (disjoint→root), `Lca({100,101,200})=10` (fold), `GetPathToRoot(100)=[100,20,10,4,3,2,1]`.
All match the code and the test assertions.

**ClassifyReads** (k=4, every key self-canonical unless noted):

| Test | hit weights | leaves | RTL path scores | assigned | C | Q | conf |
|------|-------------|--------|-----------------|----------|---|---|------|
| SingleSpecies | 100×4 | {100} | RTL(100)=4 | **100** | 4 | 4 | 1.0 |
| SplitWithinGenus | 100×2,101×2 | {100,101} | both 2 → tie | **20** = Lca(100,101) | 4 | 4 | 1.0 |
| SplitAcrossGenera | 100×2,200×2 | {100,200} | both 2 → tie | **10** = Lca(100,200) | 4 | 4 | 1.0 |
| RtlAncestorWeight | 100×1,101×2,20×1 | {100,101} (20 internal) | RTL(100)=1+1=2, RTL(101)=2+1=**3** | **101** (unique max) | 2 | 4 | 0.5 |
| NoHits | — | — | — | **1** (root) | 0 | 4 | 0.0 |
| Empty / Short | — | — | — | **1** | 0 | 0 | 0.0 |
| Ambiguous (NNNN…) | — | — | — | **1** | 0 | 0 | 0.0 |
| CanonicalLookup `AGGTT` vs DB `AACC` | window GGTT→canon AACC→100 | {100} | 1 | **100** | 1 | 2 | 0.5 |

The **RtlAncestorWeight** case is the load-bearing one for "sum along the path, not deepest hit": the
genus weight w(20)=1 is added to **both** species paths, so 101 wins 3-to-2. Confirmed by hand and by
code.

**BuildKmerDatabase**
- SharedKmer: refs `(100,"AGCTAAAA")`, `(101,"AGCTCCCC")` → `AGCT` owned by both → `Lca(100,101)=20`;
  `GCTA`→100 (E.coli-only); `GCTC` canonicalizes to `GAGC`→101 (E.fergusonii-only). All match.
- SingleReference `(100,"AAAACAA")` → 4 distinct canonical k-mers all →100. Match.

### Beyond-the-tests adversarial probes (code vs my Python, both agree)
Two cases **not** in the suite were checked by a throwaway in-tree test that I derived independently and
then deleted:
- **Genus-heavy + species:** weights {20×5, 100×1}; only leaf is 100, RTL(100)=w(100)+w(20)=1+5=**6**
  → assigned **100** (species), C=clade(100)=1, Q=6. Code matched (taxon 100, RTL 6, C 1, Q 6).
- **Family tie with ancestor weight:** weights {100×2, 200×2, 10×1}; RTL(100)=RTL(200)=2+1=3 tie →
  `Lca(100,200)=10`; C=clade(10)=5, Q=5. Code matched (taxon 10, RTL 3, C 5, Q 5).
These confirm the RTL/C/Q logic is correct on inputs the rewritten tests don't exercise.

### Test-quality audit (HARD gate)
- **Not code-echoes.** Every asserted value was reproduced by an independent from-the-paper Python
  derivation *before* running the C# suite; all agreed. No value was read off the code.
- **Real assertions.** Exact taxon ids, ranks, RTL scores, C, Q, and confidences (not "no-throw").
- **Mutation-checked this session.** Reverting the tie-break to "first leaf" fails **2** tests
  (SplitWithinGenus, SplitAcrossGenera); reverting DB-build LCA-collapse to "first-wins" fails **1**
  test (SharedKmer). Both restored; tree clean.
- Coverage includes the requested LCA primitive cases (siblings→parent, ancestor/descendant→ancestor,
  self→self, disjoint→root), DB-build LCA collapse, canonical (revcomp) lookup, and tie→LCA.

### Findings / defects
**None.** The code faithfully realises all five Kraken rules; the rewritten tests encode source-correct,
mutation-discriminating, hand-derived values. No code or test change was required.

## Verdict & follow-ups
- **Stage A: PASS** — primary source retrieved; all five rules verbatim-confirmed.
- **Stage B: PASS** — code matches an independent paper-derived reference on every test read plus two
  out-of-suite adversarial cases; tests are real and mutation-checked.
- **End-state: ✅ CLEAN.** No defect logged. (Pre-existing: the net9.0 MCP metagenomics project cannot
  build on the in-env net8 SDK — unrelated to this unit; the net8 core library and the genomics test
  project build clean and the full genomics suite is green: **6773 passed, 0 failed**.)
