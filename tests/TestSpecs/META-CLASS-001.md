# META-CLASS-001: Taxonomic Classification Test Specification

## Test Unit Information

| Field | Value |
|-------|-------|
| **ID** | META-CLASS-001 |
| **Area** | Metagenomics |
| **Canonical Methods** | `MetagenomicsAnalyzer.ClassifyReads`, `MetagenomicsAnalyzer.BuildKmerDatabase`, `TaxonomyTree` (`Lca`, parent chains) |
| **Complexity** | DB build `O(Σ m·h)`; classify `O(r·(l + leaves·h))` |
| **Invariant** | \|output\| = \|input reads\|; 0 ≤ Confidence ≤ 1; Confidence = C/Q; assigned taxon ∈ tree |

## Methods Under Test

| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `ClassifyReads(reads, kmerDb, taxonomy, k)` | MetagenomicsAnalyzer | Canonical (Kraken RTL/LCA) | Deep |
| `BuildKmerDatabase(refs, taxonomy, k)` | MetagenomicsAnalyzer | DB construction (k-mer LCA) | Deep |
| `TaxonomyTree.Lca`, `GetPathToRoot`, `IsAncestorOf`, ctor validation | TaxonomyTree | Data model | Deep |

> **Classifier scope:** This is the **faithful Kraken algorithm** (Wood & Salzberg 2014). The DB maps
> each canonical k-mer to the **LCA of its owning taxa**; per-read classification builds the
> classification tree (hit taxa + ancestors, weighted by k-mer count), assigns the leaf of the
> **maximum-scoring root-to-leaf (RTL) path**, and breaks ties by the **LCA of the maximally-scoring
> leaves**; no-hit reads are unclassified (root). Confidence = C/Q.

## Evidence Sources

1. **Wood & Salzberg (2014)** — Kraken paper; DB-build LCA, classification tree, RTL max-weight path,
   LCA-of-leaves tie-break, unclassified rule (verbatim quotes recorded in the Evidence doc, fetched
   via PMC4053813).
2. **Kraken 1 / 2 Manuals** — k=31 default, ambiguous-k-mer filtering, canonical k-mers, C/Q score.

## Hand-Built Taxonomy (for hand-derived expectations)

```
root(1) ── Archaea(5, domain)
        └─ Bacteria(2) ─ Proteobacteria(3) ─ Gammaproteobacteria(4) ─ Enterobacteriaceae(10)
              ├─ Escherichia(20, genus) ─ {E.coli(100), E.fergusonii(101)}  [species]
              └─ Salmonella(21, genus)  ─ S.enterica(200)                    [species]
```

Hand-derived LCAs: Lca(100,101)=20; Lca(100,200)=10; Lca(20,100)=20; Lca(100,100)=100;
Lca(5,100)=1; Lca({100,101,200})=10.

## Test Categories

### TaxonomyTree / LCA (hand-derived)

| ID | Test | Expected |
|----|------|----------|
| L1 | siblings → parent | Lca(100,101)=20 (order-independent) |
| L2 | ancestor/descendant → ancestor | Lca(20,100)=20 |
| L3 | same node → itself | Lca(100,100)=100 |
| L4 | disjoint branches → root | Lca(5,100)=1 |
| L5 | different genera → family | Lca(100,200)=10 |
| L6 | LCA of a set folds pairwise | Lca({100,101,200})=10 |
| L7 | path/depth/ancestry helpers | GetPathToRoot(100)=[100,20,10,4,3,2,1]; IsAncestorOf |
| L8 | ctor rejects malformed trees | no-root / two-roots / duplicate-id → ArgumentException |

### BuildKmerDatabase

| ID | Test | Evidence |
|----|------|----------|
| B1 | empty input → empty DB | robustness |
| B2 | reference shorter than k → empty | cannot extract k-mers |
| B3 | single reference → each canonical k-mer maps to that taxon | core |
| B4 | **shared k-mer collapses to LCA** (AGCT owned by 100,101 → 20) | Kraken DB-build LCA |
| B5 | ambiguous k-mers skipped; mixed case uppercased | Kraken filtering |
| B6 | reference taxon not in tree → KeyNotFoundException | validation |

### ClassifyReads (RTL / LCA)

| ID | Test | Expected (hand-derived) |
|----|------|-------------------------|
| C1 | read hitting one species | → 100, RtlScore 4, C=Q=4, conf 1.0 |
| C2 | split equally within a genus (100,101) | → genus 20 (LCA of tied leaves), C=4,Q=4 |
| C3 | split equally across genera (100,200) | → family 10 (LCA of tied leaves) |
| C4 | RTL ancestor-weight single winner (100×1,101×2,genus20×1) | → 101, RtlScore 3, C=2,Q=4,conf 0.5 |
| C5 | no k-mer hits | → root(1), C=0, conf 0, Q recorded |
| C6 | empty / shorter-than-k read | → root(1), Q=0 |
| C7 | all-ambiguous read | Q=0 → root(1) |
| C8 | canonical (reverse-complement) lookup | RC window GGTT canon AACC → 100 |
| C9 | output count & order preserved | one result per read, in order |

### Invariants & Validation

| ID | Test |
|----|------|
| I1 | Confidence ∈ [0,1]; C ≤ Q for every read |
| I2 | unclassified → C=0 and TaxonId=root |
| I3 | null reads/db/tree → ArgumentNullException; k≤0 → ArgumentOutOfRangeException |

## Mutation Note

Reverting the DB-build LCA to "first-wins" fails B4; reverting the RTL tie-break to "first leaf"
fails C2 and C3. Both mutations were applied and confirmed to fail the suite.

## Test File

**Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_TaxonomicClassification_Tests.cs`
(26 tests; LCA unit tests + DB-build + RTL classification + invariants/validation).

## Open Questions / Decisions

None. Algorithm rules are taken verbatim from Kraken (Wood & Salzberg 2014).
