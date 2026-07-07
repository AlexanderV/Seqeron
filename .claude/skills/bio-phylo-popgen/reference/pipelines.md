# bio-phylo-popgen — fuller recipes & parameter guidance

Progressive-disclosure detail for `SKILL.md`. Rigor rules (parse-with-a-tool, envelope, provenance,
cross-check, 0-based coords, alpha caveat) come from **`bio-rigor`** — not repeated here. Full I/O
schemas live in `docs/mcp/tools/{phylogenetics,population}/*.md`; formulas in
`docs/algorithms/{Phylogenetics,Population_Genetics}/*.md`.

## Distance models (phylogenetics) — pick deliberately

`pairwise_distance`, `distance_matrix`, and `build_phylogenetic_tree` all take a `method`/`distanceMethod`
in `{PDistance, JukesCantor, Kimura2Parameter, Hamming}` (default **`JukesCantor`**):

| Model | Meaning | Use when |
|---|---|---|
| `Hamming` | raw count of differing sites | already-aligned, closely related; want a raw count |
| `PDistance` | proportion of differing sites = differences / comparable sites | quick uncorrected divergence |
| `JukesCantor` (JC69) | `d = −3/4·ln(1 − 4p/3)`; corrects for multiple hits | standard nucleotide correction |
| `Kimura2Parameter` (K80) | separates transitions (S) from transversions (V) | ts/tv bias matters |

- **Saturation:** `JukesCantor` and `Kimura2Parameter` return **`+Infinity`** when divergence is too high
  (`p ≥ 0.75` for JC69). If a matrix contains `+Inf`, tree construction is unreliable — fall back to
  `PDistance`/`Hamming` or report the saturated pairs.
- Gaps (`-`) and non-standard/ambiguous bases are **skipped per column**; case is ignored.
- All sequences must be the **same length (aligned)** for `distance_matrix` / `build_phylogenetic_tree`
  (error 1002 otherwise). Align first with **`bio-alignment`** (`multiple_align`).

## Tree method — UPGMA vs Neighbor-Joining

`treeMethod` ∈ `{UPGMA, NeighborJoining}` (default **`UPGMA`**):

- **UPGMA** — rooted, ultrametric, strictly bifurcating. Assumes a **molecular clock** (equal rates).
  Fine for a quick dendrogram or clock-like data; misleading when lineages evolve at different rates.
- **Neighbor-Joining** (Saitou & Nei 1987) — **unrooted**; the central node is a trifurcation of the last
  three OTUs. No clock assumption — the safer default for real divergence data. Root it externally
  (e.g. with an outgroup) if you need a rooted tree.

`build_phylogenetic_tree` does distance + tree in one call. Use `build_tree_from_matrix` when you already
have a matrix (e.g. from an external tool or a Wikipedia worked example) — supply `taxa` in matrix order
and a **square symmetric** matrix; sizes must match (errors 1001/1002).

## Tree stats & comparison

- `tree_length` = Σ branch lengths (total tree length). `tree_depth` = max **edges** root→leaf (a count,
  not a distance). `tree_leaves` = taxon nodes with their terminal branch lengths.
- `patristic_distance`(taxonA, taxonB) = Σ branch lengths along the unique path (a real distance);
  `mrca`(taxonA, taxonB) = their most recent common ancestor node. Both need a **rooted** tree.
- `robinson_foulds_distance`(tree1, tree2) = symmetric clade-difference between two **rooted** trees;
  0 ⇒ identical topology. Use it to compare an NJ tree vs a UPGMA tree, or a bootstrap replicate vs the
  reference. `bootstrap_support` gives per-clade resampling support (report the replicate count).
- Newick round-trips: `to_newick`(tree) ⇄ `parse_newick`(newick). Newick strings end with `;`.

## Population genetics — inputs & conventions

- **Genotype encoding is `0/1/2`** everywhere (0 = hom-ref, 1 = het, 2 = hom-alt) for
  `minor_allele_frequency` and `linkage_disequilibrium` (`{geno1,geno2}` pairs). `allele_frequencies`
  and `hardy_weinberg_test` take **genotype counts** (n(AA), n(Aa), n(aa)) instead.
- **Allele frequencies:** p = (2·nAA + nAa)/(2N), q = (2·naa + nAa)/(2N), p + q = 1; all-zero input → (0,0).
- **MAF** is folded: `min(altFreq, 1−altFreq)` ∈ [0, 0.5]. `filter_variants_by_maf` filters on the stored
  `alleleFrequency` field and **does not recompute** frequencies from genotypes — make sure that field is
  populated. Both bounds are inclusive; input order preserved.
- **Hardy-Weinberg:** χ² goodness-of-fit with **1 df** vs expected `p²n, 2pqn, q²n`; `inEquilibrium` iff
  `pValue ≥ significanceLevel` (default 0.05). Zero sample → χ²=0, p=1, in equilibrium (degenerate).

## Diversity & neutrality — the k̂ trap

- `diversity_statistics` computes π, Watterson's θ, Tajima's D, S, n, and observed/expected
  heterozygosity in **one pass** from aligned sequences — prefer it when you want the whole panel.
- If you compute Tajima's D from parts, **`tajimas_d` expects k̂ = average pairwise differences (NOT
  per-site)**. `nucleotide_diversity` returns per-site π ∈ [0,1]; convert with **k̂ = π · L** (L = alignment
  length) before feeding `tajimas_d`. This is the single most common mistake in this domain.
- `D = (k̂ − S/a₁) / √(var)`, a₁ = Σ_{i=1}^{n−1} 1/i. Returns **0** when S = 0, n < 3, or variance ≤ 0.
  **Sign:** D<0 ⇒ excess rare variants (population expansion / purifying selection); D>0 ⇒ deficit of rare
  variants (balancing selection / population contraction).
- Watterson's θ = S / (a₁·L). π and θ estimate the same quantity under neutrality — their divergence is
  what Tajima's D measures; report both.

## Differentiation (Fst)

- `fst`(pop1, pop2) — Wright's variance-based `Fst = σ²_S / (p̄(1−p̄))` summed over loci. Each population is
  a list of per-variant `{alleleFreq, sampleSize}`; the two lists **must have the same per-locus count**
  (error 1001). Empty input → 0. Fst = 0 panmixia, 1 fixed difference.
- `f_statistics` adds **Fis** and **Fit** to Fst (needs genotype-level input per the tool schema).
  For >2 populations use `pairwise_fst` → symmetric matrix; the diagonal is 0.

## Linkage & haplotypes

- `linkage_disequilibrium` — from per-individual `{geno1,geno2}` (0/1/2): **r²** = squared Pearson
  correlation of genotype values (clamped [0,1]); **D'** from genotype covariance normalized by D_max
  (clamped [0,1]). Identical vectors ⇒ r²=D'=1; independent design ⇒ 0. Monomorphic/empty ⇒ 0 (no NaN).
  `distance` (bp between variants) is echoed unchanged.
- `haplotype_blocks` — chains adjacent-variant LD into blocks. `runs_of_homozygosity` → ROH segments →
  `inbreeding_from_roh` → genomic inbreeding coefficient F_ROH.

## Advanced (precomputed inputs required)

- `integrated_haplotype_score`(ehh0, ehh1, positions) — iHS from **precomputed EHH curves**, not raw reads.
- `estimate_ancestry` — proportions vs **fixed reference populations** you supply.
- `scan_selection_signals`(regions, thresholds) — combines Tajima's D / Fst / iHS **already computed
  per region**; supply those statistics, this does not derive them.

## Cross-checking playbook (satisfy `bio-rigor` rule 4)

| Result | Independent corroboration |
|---|---|
| `build_phylogenetic_tree` matrix | equals `distance_matrix`(same method) on the same aligned seqs |
| NJ vs UPGMA tree | `robinson_foulds_distance` = 0 ⇒ same topology; >0 ⇒ methods disagree (report) |
| `diversity_statistics` Tajima's D | recompute via `nucleotide_diversity`→k̂=π·L + `wattersons_theta` + `tajimas_d` |
| `allele_frequencies` p,q | p+q must equal 1; consistent with HWE expected `2pqn` |
| `fst` between two pops | pooled p̄ heterozygosity sanity: Fst ∈ [0,1]; identical pops ⇒ 0 |
| MAF | `min(altFreq, 1−altFreq)` ≤ 0.5; consistent with `allele_frequencies.minorFreq` |

## Units & coordinates (report explicitly)

- Distances are **model-specific** (Hamming = a count; PDistance/JC69/K2P = per-site, may be `+Inf`).
  `tree_length`/`patristic_distance` are in the same distance units; `tree_depth` is an **edge count**.
- π, θ (per site), MAF, Fst, D', r², heterozygosity are **fractions**; Hamming distance and S are **counts**;
  Tajima's D is a dimensionless statistic. Do not mix scales.
- Genotype counts vs genotype vectors are different input shapes — see "inputs & conventions" above.

## Envelope note

The Phylogenetics/Population tools used here are **not** among the 9 `LimitationPolicy`-guarded units
([`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md) — guarded units are
PARSE-FASTQ, CHROM-CENT, DISORDER-REGION, MIRNA-*, ONCO-MHC, ONCO-IMMUNE, META-BIN, PROBE-DESIGN,
RNA-STRUCT), so no `MinimumMode` gate applies to a normal phylo/popgen call. If a task chains into a
guarded unit elsewhere, `bio-rigor`'s envelope rule takes over — stop on a `SeqeronLimitationException`
rather than forcing output.
