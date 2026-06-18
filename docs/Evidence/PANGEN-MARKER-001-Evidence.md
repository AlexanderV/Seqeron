# Evidence Artifact: PANGEN-MARKER-001

**Test Unit ID:** PANGEN-MARKER-001
**Algorithm:** Phylogenetic Marker Selection (single-copy core genes ranked by parsimony-informative sites)
**Date Collected:** 2026-06-13

---

## Online Sources

### Ding, Baumdicker & Neher (2018), *Nucleic Acids Research* — panX

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5758898/
**Accessed:** 2026-06-13 (WebSearch query "core genome phylogenetic marker selection single-copy core genes pan-genome alignment Roary"; then WebFetch of the PMC full-text page after following the 301 redirect from ncbi.nlm.nih.gov to pmc.ncbi.nlm.nih.gov)
**Authority rank:** 1 (peer-reviewed primary paper, Nucleic Acids Research)

**Key Extracted Points:**

1. **Single-copy core gene definition (verbatim):** "those gene clusters **in which all strains are represented exactly once**" — i.e. a cluster is a phylogenetic marker only when every genome contributes exactly one gene to it (present in all strains, no paralogs / no missing strains).
2. **Marker → tree pipeline (verbatim):** "PanX **extracts all variable positions from the nucleotide alignments of all single copy core genes** (those gene clusters in which all strains are represented exactly once) to construct a core genome phylogenetic tree using FastTree."
3. **Why variable positions:** invariant (fully conserved) columns carry no phylogenetic signal; the SNP/variable-position matrix is what distinguishes strains. (Conserved-column removal is the panX selection step before tree building.)
4. **Caveat:** the core-genome tree "may not reflect true evolutionary history due to homologous recombination affecting even core genes" — selection picks markers, it does not guarantee a correct tree.

### Page et al. (2015), *Bioinformatics* — Roary

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4817141/
**Accessed:** 2026-06-13 (same WebSearch; WebFetch of the PMC full-text after following the 301 redirect)
**Authority rank:** 1 (peer-reviewed primary paper; reference pan-genome pipeline)

**Key Extracted Points:**

1. **Core gene definition (verbatim):** "Core is defined as a gene being in **at least 99% of samples**, which allows for some assembly errors in very large datasets." (Default core fraction = 0.99.)
2. **Paralog resolution (verbatim):** "Using conserved gene neighborhood information, homologous groups containing paralogs are split into groups of true orthologs." — phylogenetic markers must be single-copy orthologs, not paralog mixtures.

### Roary documentation (Sanger Pathogens) — core gene alignment behaviour

**URL:** https://sanger-pathogens.github.io/Roary/
**Accessed:** 2026-06-13 (WebFetch of the documentation page)
**Authority rank:** 3 (canonical project documentation for the Roary reference implementation)

**Key Extracted Points:**

1. **Paralog filtering (verbatim):** "By default **any cluster containing paralogous genes gets filtered out of the final core gene alignment**." — markers exclude clusters with more than one gene from any single genome.
2. **Alignment basis:** the core gene alignment is built from the per-cluster member sequences (MAFFT / PRANK); phylogenetic signal is read off the aligned columns.

### "Informative site" (Wikipedia) → Zvelebil & Baum, *Understanding Bioinformatics* (2008)

**URL:** https://en.wikipedia.org/wiki/Informative_site
**Accessed:** 2026-06-13 (WebSearch query "parsimony informative sites phylogenetic marker selection variable sites core gene alignment"; then WebFetch of the article, which cites Zvelebil 2008 as the primary source for the definition)
**Authority rank:** 4 (Wikipedia article; the numeric definition is attributed to the Zvelebil 2008 textbook — primary, authority rank 1)

**Key Extracted Points:**

1. **Parsimony-informative site definition (verbatim):** "a position in the relevant set of aligned sequences at which there are **at least two different character states and each of those states occurs in at least two of the sequences**." (Zvelebil & Baum 2008.)
2. **Exclusions (verbatim sense):** fully conserved columns and singleton columns (a variant appearing in only one sequence) are **not** parsimony-informative, "since these show the same number of evolutionary changes regardless of tree topology."
3. **Numeric criterion:** ≥ 2 distinct states in the column AND each of at least two of those states appears in ≥ 2 sequences.

---

## Documented Corner Cases and Failure Modes

### From panX / Roary

1. **Not single-copy:** a core cluster where some genome contributes 0 or ≥ 2 genes is excluded from the marker set (panX "exactly once"; Roary paralog filtering).
2. **Below core threshold:** a cluster absent from ≥ 1% of genomes is not core (Roary 99% default), so it is not a phylogenetic-marker candidate.

### From Zvelebil 2008 (parsimony-informative sites)

1. **Monomorphic column:** all sequences share one state → 0 informative sites (no signal).
2. **Singleton column:** exactly one sequence differs → not informative (the variant occurs in only one sequence).
3. **Two states, each ≥ 2 sequences:** informative (the minimal informative pattern).
4. **Gaps / unequal lengths:** parsimony-informative counting is defined over aligned columns; sequences must share a common length (a true alignment) for column-wise comparison.

---

## Test Datasets

### Dataset: Parsimony-informative-site worked columns (Zvelebil 2008 definition)

**Source:** Definition from https://en.wikipedia.org/wiki/Informative_site (Zvelebil & Baum, *Understanding Bioinformatics*, 2008).

Four aligned sequences, evaluated column by column:

| Column | s1 | s2 | s3 | s4 | States (count) | Parsimony-informative? | Reason |
|--------|----|----|----|----|----------------|------------------------|--------|
| 1 | A | A | A | A | A:4 | No | monomorphic (1 state) |
| 2 | A | A | A | C | A:3, C:1 | No | singleton (C occurs once) |
| 3 | A | A | C | C | A:2, C:2 | **Yes** | 2 states, each ≥ 2 |
| 4 | A | C | G | T | A:1,C:1,G:1,T:1 | No | 4 singletons (no state ≥ 2) |
| 5 | A | A | G | G | A:2, G:2 | **Yes** | 2 states, each ≥ 2 |

Sequences: s1=`AAAAA`, s2=`AAACA`, s3=`AACCG`, s4=`ACCTG`.
**Expected parsimony-informative sites = 2** (columns 3 and 5).

### Dataset: Single-copy core marker selection (panX/Roary definitions)

**Source:** Ding et al. 2018 (panX) "all strains represented exactly once"; Page et al. 2015 (Roary) 99% core + paralog filtering.

Three genomes g1,g2,g3; total genomes = 3; core threshold = 0.99 → a marker cluster must be present in all 3 genomes with exactly one gene per genome and ≥ 1 parsimony-informative site.

| Cluster | Genomes (count) | Genes per genome | Single-copy core? | Aligned members | PIS |
|---------|-----------------|------------------|-------------------|-----------------|-----|
| informative | g1,g2,g3 (3) | 1 each | Yes | AAC, AAC, GGC → cols A/A/G(singleton? no: A:2,G:1 → singleton)… see below | — |
| paralog | g1(2),g2,g3 (3) | g1 has 2 genes | No (paralog in g1) | — | excluded |
| not-core | g1,g2 (2) | 1 each | No (absent in g3) | — | excluded |
| conserved | g1,g2,g3 (3) | 1 each | Yes but 0 PIS | AAA,AAA,AAA | 0 → excluded |

Marker member alignment for `informative` (one sequence per genome): `ACGT`, `ACGT`, `TCGA`.
Columns: col1 A,A,T (A:2,T:1 → singleton, not PI); col2 C,C,C (mono); col3 G,G,G (mono); col4 T,T,A (singleton). PIS = 0 → would be excluded. To make a positive marker, use members `ACGT`,`ACGT`,`TCGG`,`TCGG` across **four** single-copy genomes: col1 A,A,T,T (A:2,T:2 → PI); col4 T,T,G,G (PI) → PIS = 2.

(The test fixture uses explicit alignments whose PIS is computed directly from the Zvelebil definition; selection keeps only single-copy core clusters with PIS ≥ 1 and ranks by descending PIS.)

---

## Assumptions

1. **ASSUMPTION: Per-cluster member alignment = equal-length member sequences (ungapped columns).** panX/Roary align each single-copy core cluster (MAFFT/PRANK) before reading columns. This unit has no in-repo multiple aligner, so parsimony-informative sites are counted directly over the cluster's member sequences when they share a common length (treating them as already aligned, position-by-position) — the same ungapped, position-wise convention CD-HIT global identity uses elsewhere in `PanGenomeAnalyzer`. When member sequences differ in length, no common alignment exists and PIS is defined as 0 (the cluster carries no usable column-wise signal here). This affects only how the alignment is obtained, not the parsimony-informative criterion itself (which is copied verbatim from Zvelebil 2008); the criterion and the single-copy-core selection rule are fully source-backed.

---

## Recommendations for Test Coverage

1. **MUST Test:** `CountParsimonyInformativeSites` on the 5-column worked alignment (s1..s4) returns exactly 2 (cols 3,5) — Evidence: Zvelebil 2008 definition + worked dataset.
2. **MUST Test:** monomorphic column → 0; singleton column → 0; two-state-each-≥2 column → 1 — Evidence: Zvelebil 2008 exclusions.
3. **MUST Test:** `SelectPhylogeneticMarkers` excludes non-single-copy clusters (paralog: a genome with ≥ 2 genes) — Evidence: panX "exactly once"; Roary paralog filtering.
4. **MUST Test:** excludes non-core clusters (absent from ≥ 1 genome below threshold) — Evidence: Roary 99% core; panX "all strains".
5. **MUST Test:** excludes single-copy core clusters with 0 parsimony-informative sites (fully conserved) — Evidence: panX "extracts all variable positions"; conserved columns carry no signal.
6. **MUST Test:** ranks retained markers by descending parsimony-informative-site count and caps at `maxMarkers` — Evidence: panX uses the most informative variable positions; informativeness ↔ PIS (Wikipedia/Zvelebil).
7. **MUST Test:** null/empty inputs (no clusters, no genomes) → empty marker set, no exception — Evidence: corner cases (degenerate input).
8. **SHOULD Test:** member sequences of unequal length → PIS 0 (no common alignment) → cluster not selected — Evidence: Assumption 1 (alignment precondition).
9. **COULD Test:** PIS is symmetric to sequence/row order and to relabeling states (A↔C) — Rationale: column-content property, independent of row order.

---

## References

1. Ding W, Baumdicker F, Neher RA (2018). panX: pan-genome analysis and exploration. *Nucleic Acids Research* 46(1):e5. https://doi.org/10.1093/nar/gkx977 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC5758898/)
2. Page AJ, Cummins CA, Hunt M, Wong VK, Reuter S, Holden MTG, Fookes M, Falush D, Keane JA, Parkhill J (2015). Roary: rapid large-scale prokaryote pan genome analysis. *Bioinformatics* 31(22):3691–3693. https://doi.org/10.1093/bioinformatics/btv421 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC4817141/)
3. Sanger Pathogens. Roary: the pan genome pipeline (documentation). https://sanger-pathogens.github.io/Roary/ (accessed 2026-06-13).
4. Zvelebil M, Baum JO (2008). *Understanding Bioinformatics*. Garland Science, New York — definition of parsimony-informative site (via https://en.wikipedia.org/wiki/Informative_site, accessed 2026-06-13).

---

## Change History

- **2026-06-13**: Initial documentation.
