# Validation Report: META-CLASS-001 — Kraken k-mer / LCA taxonomic classifier

- **Validated:** 2026-06-24   **Area:** Metagenomics
- **Canonical method(s):** `TaxonomyTree` (`Lca`, `GetPathToRoot`, `IsAncestorOf`, `GetDepth`, ctor validation),
  `MetagenomicsAnalyzer.BuildKmerDatabase`, `MetagenomicsAnalyzer.ClassifyReads`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN (no defect; tests independently re-derived & cross-checked; no code/test change)

Independent re-validation of the Kraken-style redesign committed at `9839268` (approved breaking change:
`TaxonomyTree`+`Lca`, `BuildKmerDatabase`, `ClassifyReads`). Every expected LCA, RTL path score,
assigned taxon, and C/Q below was **hand-derived from the five Kraken rules** via a from-scratch Python
reference written only from the paper text (`/tmp/kraken_check.py`), then the C# code and tests were
checked against that derivation — never the reverse.

## Stage A — Description

**Primary source retrieved this session (WebFetch, open-access PMC mirror):**
Wood DE, Salzberg SL. "Kraken: ultrafast metagenomic sequence classification using exact alignments."
*Genome Biology* 2014, **15:R46**. doi:10.1186/gb-2014-15-3-r46 — https://pmc.ncbi.nlm.nih.gov/articles/PMC4053813/

### The five rules, verbatim-confirmed against the fetched source

1. **DB build = LCA-collapse.** "At the core of Kraken is a database that contains records consisting
   of a k-mer and **the LCA of all organisms whose genomes contain that k-mer**." (Incremental fold:
   when a k-mer's LCA value is already set, store the LCA of the stored value and the current taxon.)
2. **Node-weighted classification tree.** "**Each node in the classification tree is weighted with the
   number of k-mers in K(S) that mapped to the taxon associated with that node.**"
3. **RTL = sum of node weights along root→leaf.** "**Each root-to-leaf (RTL) path in the classification
   tree is scored by calculating the sum of all node weights along the path.**"
4. **Max path → leaf; tie → LCA of leaves.** "**The maximum scoring RTL path in the classification tree
   is the classification path, and S is assigned the label corresponding to its leaf (if there are
   multiple maximally scoring paths, the LCA of all those paths' leaves is selected).**"
5. **No hits → unclassified.** "**Sequences for which none of the k-mers in K(S) are found in any genome
   are left unclassified by this algorithm.**" (Reported here as `TaxonomyTree.RootId`.)

**Confidence C/Q** (Kraken 2 manual, confirmed in the Evidence doc and prior fetch): "A sequence label's
score is a fraction **C/Q, where C is the number of k-mers mapped to LCA values in the clade rooted at
the label, and Q is the number of k-mers in the sequence that lack an ambiguous nucleotide**." Canonical
(min of k-mer / reverse-complement) k-mers in both build and query; ambiguous (non-ACGT) k-mers skipped.

The in-code algorithm comment (`MetagenomicsAnalyzer.cs` lines 151-170) quotes these rules; each matches
the fetched text. TestSpec, Evidence doc, and the checklist invariants all agree. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
- `TaxonomyTree.cs` — `Lca(a,b)` (163-192): collect a's ancestor set, walk b upward to first hit, root
  if disjoint; `Lca(set)` (201-212) folds pairwise; `GetPathToRoot`/`IsAncestorOf`/`GetDepth` straightforward;
  ctor (49-84) enforces exactly-one self-parented root, no duplicate ids, present parents.
- `MetagenomicsAnalyzer.cs` — `ClassifyRead` (211-292): per-taxon hit counts + Q (non-ambiguous canonical
  k-mers); leaves = hit taxa not a proper ancestor of another hit; each leaf scored by summing hit-weights
  over its **root path** (true RTL sum — ancestor weight added to every descendant path, not deepest-hit);
  single max wins, tie → `Lca(bestLeaves)`; C = Σ hit-weights with `IsAncestorOf(assigned, taxon)` (clade
  rooted at the call); confidence = C/Q. `BuildKmerDatabase` (363-397): first owner sets a canonical
  k-mer's taxon, each later owner replaces with `Lca(existing, taxonId)`. Canonical = `min(kmer, revcomp)`
  ordinal; ambiguous skipped and excluded from Q (407-411, 399-405).

### Independent from-the-paper derivation vs the code/tests (all match)

Reference written only from the five rules on the test taxonomy
`root(1) → {Archaea(5); Bacteria(2)→Proteobacteria(3)→Gammaproteobacteria(4)→Enterobacteriaceae(10)
→ {Escherichia(20)→{E.coli(100), E.fergusonii(101)}, Salmonella(21)→S.enterica(200)}}`:

**LCA primitive** — `Lca(100,101)=20`, `Lca(100,200)=10`, `Lca(20,100)=20`, `Lca(100,100)=100`,
`Lca(5,100)=1`, `Lca({100,101,200})=10`, `GetPathToRoot(100)=[100,20,10,4,3,2,1]`. All match code + tests.

**ClassifyReads** (k=4, keys self-canonical unless noted):

| Test | hit weights | leaves | RTL scores | assigned | C | Q | conf |
|------|-------------|--------|-----------|----------|---|---|------|
| SingleSpecies | 100×4 | {100} | RTL(100)=4 | **100** | 4 | 4 | 1.0 |
| SplitWithinGenus | 100×2,101×2 | {100,101} | 2,2 → tie | **20**=Lca(100,101) | 4 | 4 | 1.0 |
| SplitAcrossGenera | 100×2,200×2 | {100,200} | 2,2 → tie | **10**=Lca(100,200) | 4 | 4 | 1.0 |
| RtlAncestorWeight | 100×1,101×2,20×1 | {100,101} (20 internal) | RTL(100)=1+1=2, RTL(101)=2+1=**3** | **101** (unique max) | 2 | 4 | 0.5 |
| NoHits | — | — | — | **1** (root) | 0 | 4 | 0.0 |
| Empty/Short | — | — | — | **1** | 0 | 0 | 0.0 |
| Ambiguous (NNNN…) | — | — | — | **1** | 0 | 0 | 0.0 |
| CanonicalLookup `AGGTT` vs DB `AACC` | GGTT→canon AACC→100 | {100} | 1 | **100** | 1 | 2 | 0.5 |

**RtlAncestorWeight is the load-bearing case for "sum along the path, not deepest hit":** the genus weight
w(20)=1 is added to **both** species paths, so 101 wins 3-to-2 — confirmed by hand and by code.

**BuildKmerDatabase** — SharedKmer: `(100,"AGCTAAAA")`,`(101,"AGCTCCCC")` → palindromic `AGCT` owned by
both → `Lca(100,101)=20`; `GCTA`→100; `GCTC` canonicalizes to `GAGC`→101 (verified by computing RC).
SingleReference `(100,"AAAACAA")` → 4 distinct canonical k-mers all →100. Ambiguous/mixed-case and
unknown-taxon (KeyNotFoundException) paths correct.

### Beyond-the-tests adversarial probes (from-paper reference)
- **Genus-heavy + species** {20×5,100×1}: only leaf 100, RTL(100)=1+5=**6** → assign **100**, C=clade(100)=1, Q=6.
- **Family tie with ancestor weight** {100×2,200×2,10×1}: RTL(100)=RTL(200)=2+1=3 tie → `Lca(100,200)=10`,
  C=clade(10)=5, Q=5.

Both reproduced by my Python reference and consistent with the code's logic (path-sum RTL, clade-C).

### Test-quality audit
- **Not code-echoes** — every value reproduced by an independent from-paper Python derivation that agreed.
- **Real assertions** — exact taxon ids, ranks, RTL scores, C, Q, confidences (not no-throw / tautology).
- **Mutation note honoured (TestSpec)** — first-wins DB build would fail B4 (SharedKmer); first-leaf
  tie-break would fail C2+C3; both encoded by tests with discriminating expected values.
- Coverage includes the requested LCA primitives, DB-build LCA collapse, canonical revcomp lookup, RTL
  ancestor-weight single-winner, both tie→LCA cases, unclassified/empty/short/ambiguous, ctor validation,
  null/k≤0 argument rejection, and the [0,1]/C≤Q invariants.

### Findings / defects
**None.** The code faithfully realises all five Kraken rules plus the C/Q confidence; the tests encode
source-correct, hand-derived, mutation-discriminating values.

## Verdict & follow-ups
- **Stage A: PASS** — primary source (PMC4053813) retrieved; all five rules + C/Q verbatim-confirmed.
- **Stage B: PASS** — code matches an independent paper-derived reference on every test read plus two
  out-of-suite adversarial cases; canonical-k-mer claims verified by RC computation; tests are real and
  mutation-discriminating.
- **End-state: ✅ CLEAN.** No code or test change. Test project builds (0 warnings/errors); the unit's
  fixture is green (**26 passed, 0 failed**).
