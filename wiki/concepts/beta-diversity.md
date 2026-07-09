---
type: concept
title: "Beta diversity (between-sample dissimilarity: Bray-Curtis + Jaccard)"
tags: [metagenomics, algorithm]
sources:
  - docs/Evidence/META-BETA-001-Evidence.md
  - docs/algorithms/Metagenomics/Beta_Diversity.md
source_commit: 1eb99ff7a5a3f4819f81562f353676b522f5769f
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: meta-beta-001-evidence
      evidence: "Test Unit ID: META-BETA-001; Algorithm: CalculateBetaDiversity; Area: Metagenomics"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:alpha-diversity
      source: meta-beta-001-evidence
      evidence: "Whittaker 1960 α/β/γ framework: beta diversity = between-sample turnover, the complement of alpha = within-sample diversity; both cited by META-BETA-001-Evidence and Beta_Diversity.md §1"
      confidence: high
      status: current
---

# Beta diversity (between-sample dissimilarity: Bray-Curtis + Jaccard)

**Beta diversity** measures the dissimilarity *between two community samples* — how much the taxon
composition **turns over** from one site to another — as opposed to the diversity *within* a single
sample. It is the between-sample term of Whittaker's α/β/γ framework (Whittaker 1960): alpha =
within-sample, **beta = between-sample turnover**, gamma = regional. This is the **second ingested unit
of the Metagenomics family** and the direct sibling of [[alpha-diversity]]. Validated under test unit
**META-BETA-001**; the record is [[meta-beta-001-evidence]], [[test-unit-registry]] tracks the unit,
and [[algorithm-validation-evidence]] describes the artifact pattern.

The single entry point
`MetagenomicsAnalyzer.CalculateBetaDiversity(string sample1Name, IReadOnlyDictionary<string,double> sample1, string sample2Name, IReadOnlyDictionary<string,double> sample2)`
takes **two taxon → abundance** maps and returns one `BetaDiversity` record carrying **two dissimilarity
metrics** — `BrayCurtis` and `JaccardDistance` — plus set-composition counts (`SharedSpecies`,
`UniqueToSample1`, `UniqueToSample2`), the two sample labels, and a placeholder `UniFracDistance`. The
computation is one pass over the **union of taxon keys**; missing keys default to abundance `0`, and a
taxon counts as **present** iff its abundance is strictly `> 0`. `O(u)` time/space in the size of the
key union; deterministic.

## The two measures

For samples `j` and `k` with per-taxon abundances `N_ij`, `N_ik`:

```
Bray-Curtis      BC_jk = 1 − 2·C_jk / (S_j + S_k)                  (Bray & Curtis 1957)
                       = 1 − 2·Σ min(N_ij, N_ik) / Σ(N_ij + N_ik)
                 where C_jk = Σ min(N_ij, N_ik), S_j = Σ N_ij, S_k = Σ N_ik

Jaccard distance J = 1 − |A ∩ B| / |A ∪ B|                        (Jaccard 1901)
                   = (u_A + u_B) / (s + u_A + u_B)   (presence/absence form)
                 where A,B = the sets of taxa PRESENT (abundance > 0) in each sample,
                       s = shared, u_A / u_B = unique to A / B
```

- **Bray-Curtis** retains **abundance** information: it is the fraction of abundance mass *not* shared
  between the two samples. Range `[0,1]` (0 = identical composition, 1 = no shared species). It is
  **NOT a true distance metric** — it violates the triangle inequality — so it is a *dissimilarity*,
  not a distance, and should not be fed to algorithms that assume metric axioms without care.
- **Jaccard distance** collapses each sample to **presence/absence** and measures turnover in species
  *membership* only (abundance magnitude discarded). Range `[0,1]` (0 = identical taxon sets, 1 = no
  shared taxa). Unlike Bray-Curtis, Jaccard distance **is a true distance metric** (satisfies the
  triangle inequality).

The abundance-vs-presence contrast is the whole point of returning both: Bray-Curtis is sensitive to
abundance shifts even when the same taxa are present, whereas Jaccard changes only when the *set* of
present taxa changes.

## Invariants and edge cases

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 ≤ BrayCurtis ≤ 1`, `0 ≤ JaccardDistance ≤ 1` | both are normalized ratios of shared to total |
| INV-02 | Identical profiles → `BC = 0`, `J = 0` | shared equals the whole; `|A∩B| = |A∪B|` |
| INV-03 | No shared positive-abundance taxon (disjoint, non-empty) → `BC = 1`, `J = 1` | `C_jk = 0` / `|A∩B| = 0` over a positive denominator |
| INV-04 | Symmetric: `D(A,B) = D(B,A)` for both | `min`, `∩`, `∪` are commutative |
| INV-05 | Abundance sensitivity: Bray-Curtis uses abundances; Jaccard uses presence/absence only | by construction |

- **Empty / all-zero denominators are guarded.** When both samples are empty (or all abundances are
  `0`) both metrics return **`0`** (identical inputs ⇒ zero dissimilarity). Per the algorithm doc's
  edge-case table, when *one* sample is empty and the other has present taxa, the distance is **`1`**
  (maximum), consistent with the `scipy.spatial.distance` convention.
- **Presence rule `> 0`:** a taxon present in one sample with value `0` in the other is counted as
  unique to the positive-abundance sample; zero/negative abundance = absent for the Jaccard branch.
- **No internal normalization** of Bray-Curtis inputs: the formula is applied directly to the supplied
  values, so the caller controls whether they are raw counts or relative abundances.

## Worked oracles (from [[meta-beta-001-evidence]])

Wikipedia fish-tank example — Tank 1 `{Goldfish:6, Guppy:7, Rainbow:4}` (S=17) vs Tank 2
`{Goldfish:10, Guppy:0, Rainbow:6}` (S=16): shared min-sum `C = 6+0+4 = 10` →
**BC = 1 − 20/33 = 13/33 ≈ 0.3939**; Jaccard shared `{Goldfish, Rainbow} = 2` / union `= 3` →
**J = 1 − 2/3 = 1/3 ≈ 0.3333** (Guppy absent in Tank 2). The symmetry oracle gives exact
**BC = 3/5, Jaccard = 2/3**. Identical samples → `0/0`; disjoint non-empty → `1/1`.

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for the two classical
beta-diversity measures over a *supplied* pair of abundance maps. Both formulas match their primary
sources exactly (Bray & Curtis 1957, Jaccard 1901; Whittaker 1960 for the framing). The single accepted
**deviation** (algorithm doc §5.4) is the **`UniFracDistance` placeholder**: it is hard-coded to `0`
because the phylogenetic UniFrac measure needs a tree that is outside this method's contract — the field
exists on the record for future extension only. No phylogenetic weighting, no sampling-effort model, no
rarefaction; the Jaccard branch can be dominated by rare taxa because abundance is discarded.

## Relation to siblings

- **Sibling of [[alpha-diversity]]** — the same Whittaker 1960 α/β/γ decomposition: alpha summarizes
  diversity *within* one abundance map (Shannon/Simpson/Chao1/…); beta summarizes turnover *between*
  two maps. They are the within- vs between-sample halves of the metagenomics diversity topic and share
  the taxon→abundance input shape.
- **Ecological Jaccard vs k-mer Jaccard.** The Jaccard-distance branch uses the *same index math*
  (`1 − |A∩B|/|A∪B|`) as the alignment-free [[kmer-jaccard-similarity]] unit, but over a different
  domain: here A and B are sets of **taxa present in two community samples** (ecological turnover),
  whereas there they are sets of **k-mers of two sequences** (sequence resemblance). Same formula,
  ecologically vs compositionally distinct — not interchangeable.
</content>
