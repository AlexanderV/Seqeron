---
type: synthesis
title: "Silent traps in an OLC → ANI → RBH → variant → phylogeny pipeline (and how to work around each)"
tags: [assembly, comparative-genomics, variants, phylogenetics, gotcha]
sources:
  - docs/algorithms/Assembly/Overlap_Layout_Consensus.md
  - docs/algorithms/Comparative_Genomics/Average_Nucleotide_Identity.md
  - docs/algorithms/Comparative_Genomics/Reciprocal_Best_Hits.md
  - docs/algorithms/Variants/Variant_Detection.md
  - docs/algorithms/Variants/Indel_Detection.md
  - docs/algorithms/Phylogenetics/Distance_Matrix.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# Silent traps in an OLC → ANI → RBH → variant → phylogeny pipeline

A common comparative-genomics pipeline — **assemble reads with OLC → confirm species by ANI to a
reference → find orthologs by RBH → call SNPs/indels against the reference → build a phylogeny on
p-distance** — runs end-to-end without raising a single error, yet each stage hides a semantic trap
that yields a *formally plausible but scientifically wrong* result. This page collects the five
per-stage traps, the concrete workaround for each, and one cross-cutting cascade effect that no
single gotcha page captures. Each trap has its own atomic gotcha page; this is the pipeline-level
view.

## The five per-stage traps and their workarounds

| Stage | Trap | Workaround |
|---|---|---|
| **1. OLC assembly** | `assemble_olc` is **greedy-SCS**, order-dependent, and assumes **error-free, exact-match (identity 1.0)** reads; a different read order can change the contig set, and error-containing reads don't overlap under the exact-match path. | Error-correct reads first ([[kmer-spectrum-error-correction]]); don't treat one greedy run as globally optimal; inspect the layout via `find_overlap` / `find_all_overlaps`. See [[olc-assembly-is-greedy-not-optimal]]. |
| **2. ANI to reference** | `calculate_ani` fragments **only the query**, so ANI(A,B) ≠ ANI(B,A); the ~95 % species cutoff depends on direction and two directions can straddle the boundary. Non-conserved fragments are dropped (not scored 0); a reference shorter than one 1020-nt fragment returns ANI = 0 ("not measurable", not "0 %"). | Use **`CalculateReciprocalAni`** (mean of both directions) for species delimitation / clustering / distance matrices; if single-direction, record which genome was query. See [[ani-is-directional-use-reciprocal]]. |
| **3. RBH orthologs** | `find_orthologs` / `find_reciprocal_best_hits` rank best hits by **5-mer Jaccard** (alignment-free, Mash-style), not BLAST bit-score; among near-identical candidates the "winning" hit can differ from a BLAST-RBH pipeline. Reciprocity rule, coverage gate, and min-similarity threshold are source-backed — only the ranking metric is substituted. | Use as a fast alignment-free screen; confirm ambiguous best-hit pairs with an alignment-based reciprocal search. Identical sequences still score 1.0. See [[rbh-orthologs-use-kmer-jaccard-not-blast]]. |
| **4. SNP/indel calling** | `find_snps` / `find_indels` call from a **single query-vs-reference global alignment**, not a read pileup — no depth, no base quality, no diploid genotypes (no GT/DP/GQ, no allele fractions). Indels are **not left-aligned**; in tandem-repeat / low-complexity regions the reported coordinate can differ from a normalized VCF (same count/type, different position). | Left-align + parsimony-normalize **both sides** before comparing coordinates across callers; for read-based genotyping use a real pileup caller. See [[variant-calling-is-alignment-not-pileup]]. |
| **5. p-distance phylogeny** | p-distance is uncorrected Hamming proportion — systematically **underestimates** divergence, so branch lengths/topology from `build_tree_from_matrix` are wrong for divergent input. Switching to correction gives the opposite failure: **JC69 → +∞ at p ≥ 3/4**, **K2P → +∞** when the log argument goes non-positive, breaking NJ/UPGMA if unhandled. | Use JC69/K2P (or richer) for divergent sequences and handle `+∞` as an explicit saturation sentinel; reserve p-distance for closely related sequences. Model choice is a scientific decision, not a formatting toggle. See [[p-distance-is-uncorrected-and-jc-k2p-saturate]]. |

## The cross-cutting cascade (not on any single gotcha page)

Stages 3 and 4 both operate **against the assembly from stage 1**, which is itself a
non-deterministic greedy approximation built on an error-free assumption. If the reads carried
sequencing errors, the assembly is already shifted, and every downstream coordinate — indel
positions, ANI fragment placements, ortholog boundaries — inherits that shift **even after each
individual trap above is fixed**. The correct order of trust is the reverse of the pipeline order:
validate the assembly first, then everything that stands on it.

## Umbrella limitation

All of the above sits under [[research-grade-limitations]]: Seqeron is **beta, not for clinical or
diagnostic use**, with no external audit, and many algorithms are faithful but **simplified / subset**
realisations of fuller published methods (scope in `docs/Validation/LIMITATIONS.md`, guarded at
runtime by `LimitationPolicy`). Before relying on any output with real data, independently verify it
against established tools for the specific use case.

## See also

Full models: [[overlap-layout-consensus-assembly]], [[average-nucleotide-identity]],
[[ortholog-detection-reciprocal-best-hits]], [[germline-variant-calling-snp-indel]],
[[evolutionary-distance-matrix]].
