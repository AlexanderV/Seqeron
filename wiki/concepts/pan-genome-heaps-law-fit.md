---
type: concept
title: "Pan-genome Heaps'-law fit (power-law new-gene curve → open/closed)"
tags: [comparative-genomics, pan-genome, algorithm]
mcp_tools:
  - fit_heaps_law
sources:
  - docs/algorithms/PanGenome/Pan_Genome_Growth_Model.md
  - docs/Evidence/PANGEN-HEAP-001-Evidence.md
source_commit: cc39f4dd59b853125c9eb4985774c4f4df018ad8
created: 2026-07-10
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: pangen-heap-001-evidence
      evidence: "Test Unit ID: PANGEN-HEAP-001 ... Algorithm: Pan-Genome Growth Model (Heaps' law fit + gene presence/absence matrix)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:pan-genome-core-accessory-partition
      source: pangen-heap-001-evidence
      evidence: "micropan open/closed rule 'if alpha<1.0 the pan-genome is open, if alpha>1.0 it is closed' — the decay exponent this unit fits is exactly the open/closed classifier the ConstructPanGenome partition (PANGEN-CORE-001) reports"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:pan-genome-gene-clustering
      source: pangen-heap-001-evidence
      evidence: "Coverage rec 7: 'dictionary overload clusters then fits and agrees with the matrix overload' — the convenience wrapper clusters genes into families, then binarizes the presence/absence matrix this unit fits"
      confidence: medium
      status: current
---

# Pan-genome Heaps'-law fit (power-law new-gene curve → open/closed)

The **Heaps'-law fit** is the **pan-genome family**'s (`PANGEN-*`) dedicated **growth-model** unit:
given a gene **presence/absence matrix** over N genomes, it counts how many **new** gene clusters
each added genome contributes (the rarefaction / new-gene curve), fits the **power law**
`y = K·x^(-alpha)` to that curve by least squares, and classifies the pan-genome **open vs closed**
by the fitted **decay exponent alpha**. It is the fitting *engine* underneath the open/closed facet
of the occupancy partition [[pan-genome-core-accessory-partition]] (PANGEN-CORE-001, which reports the
same alpha classification as one output among core/accessory/unique + fluidity) and consumes the gene
families produced upstream by [[pan-genome-gene-clustering]] (PANGEN-CLUSTER-001). Validated under
test unit **PANGEN-HEAP-001**; the validation record is [[pangen-heap-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the artifact
pattern.

## The new-gene curve (micropan `heaps()` first-appearance rule)

The input is a **presence/absence matrix** (genomes × gene clusters). It is first **binarized** — any
count > 0 becomes 1, so copy-number and duplicate gene ids collapse to mere presence (presence/absence,
not abundance). For a given genome ordering, a cluster is counted as **new at genome i** (i ≥ 2) when
its cumulative presence is **1 at row i and 0 at row i−1** — i.e. it first appears at genome i:

```
cm    = cumsum over the ordered rows          # micropan: apply(pan.matrix, 2, cumsum)
new_i = rowSums( (cm[i] == 1) & (cm[i-1] == 0) )   # first-appearance of a cluster at genome i
```

The **genome index starts at 2** (the first genome has no predecessor, contributes no "new" count),
so the fit uses N = 2..G. Because the ordering is arbitrary, micropan **pools the curve over many
random permutations** (`n.perm`, default 100 — "certainly a minimum") and fits all pooled points at
once. This implementation uses a fixed seed and the natural input order for single-permutation
deterministic tests; the averaging principle matches the source.

## The power-law fit (bounded least squares)

Heaps' law here is `y = K·x^(-alpha)` — new clusters `y` as a decreasing function of genome index `x`.
The two parameters are fitted by minimizing the least-squares objective over box constraints:

```
J(K, alpha) = sqrt( Σ (y − K·x^(-alpha))² ) / |x|      # micropan objectFun
K     ∈ [0, 10000]      alpha ∈ [0, 2]                 # L-BFGS-B lower/upper bounds
start:  K0 = mean(y at x==2),   alpha0 = 1
```

The return is the named pair `(Intercept = K, alpha)`. micropan calls R's `optim(method="L-BFGS-B")`;
this implementation minimizes the **identical** objective over the **identical** bounds from the
**identical** start point using a deterministic bounded coordinate descent. For data lying exactly on a
power curve within the bounds the global minimum is unique, so both optimizers reach it; recovered
`(K, alpha)` matches the analytic solution to < 1e-9 (see oracles).

## Open vs closed (the alpha classifier)

micropan's verbatim rule: *"if alpha<1.0 the pan-genome is open, if alpha>1.0 it is closed."*

- **Open** ⟺ **alpha < 1** — new genes keep accumulating, no asymptote (e.g. *E. coli*).
- **Closed** ⟺ **alpha > 1** — few new genes per added genome, pan-genome size approaches a limit.
- The boundary **alpha = 1** is treated by the strict inequality as **not open** (closed).

This is the same decay-exponent openness that the Tettelin et al. 2008 review frames as a power law
(rather than only Tettelin 2005's exponential new-gene decay `F_s = κs·exp[−n/τs] + tg(θ)`).

## Worked oracles (exact power-curve fits)

Two points on a single power curve determine `(K, alpha)` uniquely with objective J = 0:

| x (genome) | y (new clusters) | fit | classification |
|---|---|---|---|
| 2, 3 | 8, 4 | alpha = ln 2 / ln(3/2) ≈ **1.70951**, K = 8·2^alpha ≈ **26.16400** | alpha > 1 → **closed** |
| 2, 3 | 1, 1 | constant curve → **alpha = 0**, **K = 1** | alpha < 1 → **open** |

The closed oracle solves `8/4 = (3/2)^alpha`; the open (boundary) oracle is the best power fit to a
constant new-gene curve. Both are realized by fixed-order presence/absence matrices where genome 2 (and
genome 3) introduce exactly the stated number of previously-absent clusters.

Qualitative anchor (Tettelin et al. 2005, *S. agalactiae*): the 2nd genome added **161** new genes,
the 5th added **54**, the core is ≈ 80% of any single genome, and the asymptotic new genes per genome
tg(θ) = **33 ± 3.5** is nonzero (p < 6×10⁻⁴ that it equals zero) → an **open** pan-genome.

## Edge cases

- **Fewer than 2 genomes** → the curve `x = 2:ng` is empty; no fit is defined. The contract returns a
  degenerate fit (Intercept = 0, predictor → 0), **not** an exception; empty/null input likewise.
- **Genome index starts at 2** — genome 1 has no predecessor, so it contributes no new-gene count.
- **Binary presence only** — copy-number > 1 or duplicate gene ids are collapsed to 1 before counting.
- **alpha boundary** — alpha is clamped to [0, 2]; alpha = 1 is the open/closed boundary and, by strict
  inequality, classifies as not-open (closed).

## Assumptions (source-backed)

1. **Optimizer method (assumption).** The optimization *method* (deterministic bounded coordinate
   descent vs micropan's L-BFGS-B) is non-correctness-affecting: the objective, box constraints, and
   start point — which determine the result — are copied verbatim from `powerlaw.R`, and for exact
   power-curve data within the bounds the minimum is unique and both reach it (< 1e-9 agreement).
2. **Permutation RNG (assumption).** micropan uses R's `sample()`; the pooled curve is
   permutation-dependent except for permutation-invariant matrices. This implementation fixes the seed
   and uses natural input order for the first permutation, so single-permutation fixed-order fits are
   exactly reproducible; the pool-over-orderings averaging principle matches the source.

## Contract and implementation (`PanGenomeAnalyzer`)

The primary spec `docs/algorithms/PanGenome/Pan_Genome_Growth_Model.md` marks the unit
**Production**. Entry points live in `PanGenomeAnalyzer.cs` (`Seqeron.Genomics.Metagenomics`):

- **`FitHeapsLaw(IEnumerable<GenePresenceRow>, int permutations = 100)`** — the canonical fit over a
  presence/absence matrix. Presence is read from `GenePresenceRow.GenePresence` (true = present) and
  collapsed to binary; cluster columns are stabilized by first appearance so the matrix is
  deterministic. `permutations` is clamped to ≥ 1.
- **`FitHeapsLaw(IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>,
  double identityThreshold = 0.9, int permutations = 100)`** — convenience overload that clusters the
  genomes (CD-HIT-style at `identityThreshold`), builds the matrix, then delegates to the canonical fit.
- **`CreatePresenceAbsenceMatrix(genomes, clusters)`** — builds the binary matrix the fit consumes.

The return is `(Intercept = K, Alpha = α, IsOpen = α < 1, PredictNewGenes)`. **`PredictNewGenes`** is a
`Func<int,double>` realizing the fitted predictor `N ↦ K·N^(−α)` — the expected number of new gene
clusters at the N-th genome; it is **non-increasing in N** for α ≥ 0 (INV-06). Null/empty input or
fewer than 2 genomes yields the degenerate fit `(0, 0, false, predictor→0)` rather than an exception.

Invariants **INV-01..INV-06** (open ⇔ α<1; first-appearance new-gene rule; binary presence; fitted
α∈[0,2] and K∈[0,10000]; exact-power-curve recovery; predictor monotonicity) are covered by
`PanGenomeAnalyzer_FitHeapsLaw_Tests.cs`. The repo **suffix tree is deliberately not used** — this unit
counts set first-appearances and fits a power curve, doing no substring/occurrence search.

Complexity: new-gene curve over all permutations O(P·G·C), objective evaluation O(P·G) per point pool,
bounded minimization O(I·P·G) with I deterministic coordinate-descent iterations (P permutations,
G genomes, C clusters).

## Reference tools

Definitions trace to **micropan `heaps()`** (Snipen & Liland, `R/powerlaw.R`) for the model
`y = K·x^(-alpha)`, the first-appearance new-gene counting rule, the least-squares objective + bounds,
and the verbatim open/closed criterion; **Tettelin, Riley, Cattuto, Medini 2008** (Curr Opin Microbiol
11(5):472–477) for the power-law framing of pan-genome growth and the open/closed dichotomy; and
**Tettelin et al. 2005** (PNAS 102(39):13950–13955, the paper that coined "pan-genome") for the
qualitative *S. agalactiae* new-gene anchor (161/54 new genes, tg(θ) = 33, ≈ 80% core, open). No source
contradictions — the power-law model, counting rule, and openness criterion are mutually consistent.
