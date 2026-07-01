# bio-metagenomics — fuller recipes & parameter guidance

Progressive-disclosure detail for `SKILL.md`. Rigor rules (parse-with-a-tool, envelope, provenance,
cross-check, units, alpha caveat) come from **`bio-rigor`** — not repeated here.

## The classification chain (database → classify → profile)

Three tools compose into a metagenome workflow; **k must be identical across the chain** (default 31):

1. **`build_kmer_database`**(referenceGenomes=`[{taxonId,sequence}]`, taxonomy=`[{id,name,rank,parentId}]`, k)
   → `entries[{kmer,taxonId}]`, `count`. Each length-`k` window is **canonicalised** (lexicographic
   min of the window and its reverse complement). A canonical k-mer contributed by several taxa is
   collapsed to the **LCA** of those taxa. Non-ACGT windows are skipped; input is upper-cased. A
   reference whose `taxonId` is absent from the taxonomy tree raises error 4001 (`KeyNotFoundException`).
2. **`classify_reads`**(reads=`[{id,sequence}]`, kmerDatabase=entries, taxonomy, k)
   → `items[]` with `taxonId`, `taxonName`, `rank`, `rtlScore`, `confidence`, `matchedKmers` (C),
   `totalKmers` (Q), and a `kingdom…species` lineage. Per read: collect canonical k-mer hits, build
   the Kraken classification tree weighted by k-mer count, score every root-to-leaf path, assign the
   leaf of the max-scoring path; ties → LCA of the tied leaves; no hits → root (unclassified).
   **Confidence** is Kraken 2's C/Q = clade k-mers / non-ambiguous k-mers queried. `k` must be
   positive (error 1002); null args → 1003.
3. **`taxonomic_profile`**(classifications=items) → `kingdom/phylum/genus/speciesAbundance` (`{name,fraction}`),
   `shannonDiversity`, `simpsonDiversity`, `totalReads`, `classifiedReads`. Reads with empty or
   `"Unclassified"` kingdom are **excluded from the abundance denominator and from `classifiedReads`**,
   but `totalReads` counts every input read. Empty input → empty profile (no error).

**Gotchas:**
- A `classify_reads`→`taxonomic_profile` species `name` is the species epithet as stored in taxonomy;
  match names consistently when you later feed abundances into diversity tools.
- Because unclassified reads are dropped from the profile denominator, `Σ fraction = 1` over
  *classified* reads only — state `classifiedReads/totalReads` in provenance.

## Alpha diversity (`alpha_diversity`)

Input `abundances=[{name,fraction}]`; **raw counts or fractions both accepted** (normalised
internally); zero-abundance entries dropped. Outputs:

| Field | Formula | Note |
|---|---|---|
| `shannonIndex` | `H = -Σ pᵢ ln pᵢ` | natural log (nats) |
| `simpsonIndex` | `λ = Σ pᵢ²` | concentration (higher = less diverse) |
| `inverseSimpson` | `1/λ` | Hill number of order 2 |
| `chao1Estimate` | Chao1 richness | **needs integer counts**; with fractions falls back to observed species |
| `observedSpecies` | count of `pᵢ>0` | richness |
| `pielouEvenness` | `J = H / ln(S)` | 0 when `S ≤ 1` |

**Cross-check:** for a profile from the classification chain, `alpha_diversity(speciesAbundance)`'s
`shannonIndex`/`simpsonIndex` must equal `taxonomic_profile`'s `shannonDiversity`/`simpsonDiversity`
(different code paths, same species vector) — a mismatch flags a name/normalisation bug.

## Beta diversity (`beta_diversity`) & differential abundance

`beta_diversity`(sample1Name, sample1, sample2Name, sample2) with two `[{name,fraction}]` vectors:
- **Bray-Curtis** `1 − 2·Σ min(aᵢ,bᵢ) / Σ(aᵢ+bᵢ)` — uses abundances.
- **Jaccard distance** `1 − |shared|/|union|` — presence/absence, "present" = abundance `>0`.
- `sharedSpecies`, `uniqueToSample1`, `uniqueToSample2`; `uniFracDistance` is **always 0** (no tree supplied).
- Empty vectors → zero distances (no error).

For a **per-taxon test across two condition groups** (not a single pairwise distance), use
`differential_abundance` instead of `beta_diversity`.

## Binning (`bin_contigs`) — ⚠ META-BIN-001, min mode `Moderate`

Input `contigs=[{contigId,sequence,coverage}]`; optional `numBins` (=k, default 10),
`minBinSize` (bp, default 500000), `expectedGenomeSize` (bp, default 4000000). Bins by **k-means over
three features**: GC content, normalised coverage, and tetranucleotide (TETRA) frequency (Teeling
2004, Pearson distance). Centroids are **deterministically seeded** across GC-sorted contigs →
reproducible. Per bin with total length `≥ minBinSize`:

- `completeness = min(totalLength / expectedGenomeSize × 100, 100)`
- `contamination = min(stddev(GC)/0.5 × 100, 100)` (within-bin GC variance, CheckM-style, Parks 2014)
- `gcContent`, `coverage`, member `contigIds`. `predictedTaxonomy` is reserved/empty.
- Bins below `minBinSize` are dropped; empty/all-small input → empty list.

**Envelope / STOP rule.** `MetagenomicsAnalyzer.BinContigs` is guarded by **META-BIN-001** (domain-level
CheckM completeness/contamination). Minimum access mode **`Moderate`** (the default): allowed under
`Moderate`/`Permissive`, **throws `SeqeronLimitationException` under `Strict`**. If a binning call
throws, **stop and report the limitation** — do not force output. In C#, only bootstrap
`Moderate`/`Permissive` when the caller explicitly accepts the domain-level approximation. Always
label completeness/contamination as *domain-level estimates*, not marker-based CheckM values.
Source of truth: [`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md);
algorithm: [`docs/algorithms/Metagenomics/Genome_Binning.md`](../../../../docs/algorithms/Metagenomics/Genome_Binning.md).

## Function, resistance, functional diversity

- **`predict_functions`** — predict functional annotations for a protein set; feed its output to
  `functional_diversity`.
- **`functional_diversity`** — functional richness + Shannon over functions (an alpha-style index in
  function space rather than taxon space).
- **`find_resistance_genes`** — screen genes for antibiotic-resistance markers
  ([`Antibiotic_Resistance_Detection.md`](../../../../docs/algorithms/Metagenomics/Antibiotic_Resistance_Detection.md)).

## Pan-genome (PanGenomeAnalyzer) — comparative multi-genome

Typical order: `construct_pangenome` → then any of `cluster_genes`, `core_gene_clusters`,
`accessory_genes`, `gene_presence_absence_matrix`, `find_genome_specific_genes`,
`core_genome_alignment`, `select_phylogenetic_markers`, `fit_heaps_law`. Use
`select_phylogenetic_markers` (single-copy core) as input to a downstream tree — hand markers to
**`bio-alignment`** (MSA) then **`bio-phylo-popgen`** (tree). `fit_heaps_law` classifies the
pan-genome as **open** (new genes keep appearing) or **closed**. Background:
[`PanGenome_Core_Accessory.md`](../../../../docs/algorithms/Metagenomics/PanGenome_Core_Accessory.md).

## Cross-checking playbook (satisfy `bio-rigor` rule 4)

| Result | Independent corroboration |
|---|---|
| profile Shannon/Simpson | `alpha_diversity(speciesAbundance)` must reproduce them |
| classification confidence | recompute C/Q = `matchedKmers`/`totalKmers` from the same item |
| beta Jaccard | equals `1 − sharedSpecies/(sharedSpecies+uniqueToSample1+uniqueToSample2)` |
| bin completeness | `totalLength / expectedGenomeSize` — sanity vs the reported percentage |
| core vs accessory | core ∪ accessory clusters ≈ full pan-genome cluster set |

## Units & conventions (report explicitly)

- Abundances are **fractions in [0,1]** (or raw counts, normalised internally). Diversity indices are
  unitless; Shannon is in **nats** (natural log). Simpson higher = less diverse; inverse Simpson higher
  = more diverse.
- Lengths (`totalLength`, `minBinSize`, `expectedGenomeSize`) are in **bp**; `completeness`/
  `contamination` are **percentages** (0–100).
- `k` (k-mer length) must be **identical** for `build_kmer_database` and `classify_reads`.
- State `classifiedReads/totalReads` whenever you report a profile — the abundance denominator excludes
  unclassified reads.
