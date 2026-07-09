---
type: concept
title: "Taxonomic profile (community abundance profiling)"
tags: [metagenomics, algorithm]
sources:
  - docs/Evidence/META-PROF-001-Evidence.md
  - docs/algorithms/Metagenomics/Taxonomic_Profile.md
source_commit: 02f28f4a5999dc16f47afa3db37eca0cb1eda2ee
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: meta-prof-001-evidence
      evidence: "Test Unit ID: META-PROF-001, Area: Metagenomics, Methods GenerateTaxonomicProfile"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:taxonomic-classification
      source: meta-prof-001-evidence
      evidence: "Input: IEnumerable<TaxonomicClassification> — the profile aggregates the per-read classification output into per-taxon counts and relative abundances; 'taxonomic profiling aggregates classification results into abundance distributions' (Wikipedia Metagenomics)."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:alpha-diversity
      source: meta-prof-001-evidence
      evidence: "The profile computes Shannon H=−Σpᵢln(pᵢ) and Simpson λ=Σpᵢ² inline from its species-level abundances; alpha-diversity computes the same two (plus four more) indices over a supplied taxon→abundance map."
      confidence: high
      status: current
---

# Taxonomic profile (community abundance profiling)

A **taxonomic profile** is the community-composition summary of a metagenome: it **aggregates
per-read taxonomic classifications into relative-abundance estimates per taxon**, at several
taxonomic ranks. This is the aggregation/estimation step *downstream* of per-read assignment —
`MetagenomicsAnalyzer.GenerateTaxonomicProfile(IEnumerable<TaxonomicClassification>)` takes the
output of [[taxonomic-classification]] (Kraken k-mer/LCA per-read labels) and returns a
`TaxonomicProfile` record. It is a sibling of the diversity pair [[alpha-diversity]] /
[[beta-diversity]] (which summarize the abundance profile it produces), [[metagenomic-binning]],
and [[functional-prediction]] / [[pathway-enrichment-ora]]. Validated under test unit
**META-PROF-001**; the record is [[meta-prof-001-evidence]], [[test-unit-registry]] tracks the
unit, and [[algorithm-validation-evidence]] describes the artifact pattern.

This concept is the **profiling** unit that [[taxonomic-classification]] explicitly left out of
scope: classification assigns one label *per read*; profiling turns the *bag* of those labels
into a normalized community-abundance distribution. Canonical profilers: MetaPhlAn (Segata et al.
2012, clade-specific markers) and AMPHORA — relative abundances at kingdom → phylum → … → species.

## What the profile contains

The `TaxonomicProfile` record carries three things (impl `MetagenomicsAnalyzer.cs` lines ~229–280):

1. **Abundance maps at four ranks** — kingdom, phylum, genus, species. Each is a
   `taxon → relative-abundance` map where **relative abundance = count(taxon) / Σcount(all
   classified taxa at that rank)**, so the values sum to ≈ 1.0 per rank (for classified reads).
2. **Diversity metrics** — `Shannon` `H = −Σ pᵢ ln(pᵢ)` and `Simpson` `λ = Σ pᵢ²`, computed from
   the **species-level** abundance distribution (natural log ⇒ nats; λ = concentration index, not
   Gini-Simpson `1−λ` nor inverse `1/λ`). The fuller six-index within-sample summary is the
   separate [[alpha-diversity]] unit.
3. **Read counts** — `TotalReads` (all input classifications) and `ClassifiedReads` (those with a
   real taxon, i.e. excluding `Unclassified`).

## Counting rules

- **"Unclassified" reads are excluded** from every abundance denominator (only classified reads
  contribute to the community composition).
- **Empty rank strings are filtered** from the rank-specific maps — a read that lacks (say) a
  genus label does not create an empty-key abundance entry.
- **Relative-abundance normalization** divides by the total *classified* count, so each rank's
  map is a proper probability distribution before Shannon/Simpson are taken.

## Invariants and edge cases

- **Sum invariant:** Σ(abundance values) ≈ 1.0 at each rank (over classified reads).
- **Count invariant:** `ClassifiedReads ≤ TotalReads`.
- **Consistency:** `ClassifiedReads = Σ(counts at any rank)` (all classified reads carry every
  populated rank).
- **Diversity bounds:** `Shannon ≥ 0`; `0 ≤ Simpson ≤ 1.0`.
- **Single taxon:** abundance `1.0`, `Shannon = 0` (no uncertainty), `Simpson = 1.0` (certainty of
  drawing the same species).
- **Uniform distribution over `n`:** `Shannon = ln(n)`, `Simpson = 1/n`; skewed → low Shannon,
  high Simpson.
- **Empty input:** `TotalReads = 0`, `ClassifiedReads = 0`, empty abundance maps.
- **All Unclassified:** `ClassifiedReads = 0`, no abundance entries.
- **Empty-sum convention:** with no terms, `Shannon = −Σ(∅) = 0` and `Simpson = Σ(∅) = 0` — a
  mathematical fact (note: `Simpson = 0` for the *empty* profile, vs `1.0` for a *single* taxon).

Worked oracles (from [[meta-prof-001-evidence]]): 3 uniform species → `Shannon = ln 3`; the
non-uniform count vector `[2,1,1]` → `Simpson = 0.375`; `TotalReads = 3 / ClassifiedReads = 2`
when one of three reads is Unclassified.

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for read-tally community
profiling. The Evidence file lists **no source contradictions** and only design decisions that are
mathematical facts (Shannon in nats via `Math.Log`; Simpson as the concentration index `λ = Σpᵢ²`;
empty input → both diversity metrics 0). It is a **count-based** profiler: it tallies the reads a
classifier already assigned and normalizes them — it does **not** implement MetaPhlAn's
marker-gene coverage estimation, does no genome-size / copy-number correction, and inherits the
accuracy of whatever [[taxonomic-classification]] produced upstream. The richer within-sample index
set (inverse-Simpson, Pielou evenness, Chao1) lives in [[alpha-diversity]]; between-sample turnover
in [[beta-diversity]]; and the between-**group** differential-abundance test that consumes these
per-sample abundance vectors is [[significant-taxa-detection]] (Mann–Whitney U / Wilcoxon rank-sum).
