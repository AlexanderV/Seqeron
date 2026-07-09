---
type: concept
title: "Metagenomic binning (contigs → MAGs: k-means on composition + coverage, TETRA z-score signature)"
tags: [metagenomics, algorithm]
sources:
  - docs/Evidence/META-BIN-001-Evidence.md
  - docs/Evidence/META-BIN-001-MarkerQC-Evidence.md
  - docs/algorithms/Metagenomics/Genome_Binning.md
source_commit: c21e0be5032ffdc39a0e53405a4c7fbd09482958
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: meta-bin-001-evidence
      evidence: "Test Unit ID: META-BIN-001; Algorithm: BinContigs (+ CalculateTetranucleotideZScores); Area: Metagenomics"
      confidence: high
      status: current
---

# Metagenomic binning (contigs → MAGs)

**Metagenomic binning** groups the contigs of a mixed-community assembly back to their **genome of
origin**, producing **Metagenome-Assembled Genomes (MAGs)**. A mixed sample assembles into thousands
of fragments from many organisms; binning re-clusters them using signals that are roughly
genome-constant across a genome but differ between genomes — **compositional** signals (GC content,
tetranucleotide frequencies) and **coverage** signals (read depth, differing across samples). This is
the **third ingested unit of the Metagenomics family**, alongside the diversity siblings
[[alpha-diversity]] and [[beta-diversity]] (which quantify *community* diversity rather than
reconstructing genomes). Validated under test unit **META-BIN-001**; the record is
[[meta-bin-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

The entry point
`MetagenomicsAnalyzer.BinContigs(IEnumerable<(string ContigId, string Sequence, double Coverage)> contigs, int numBins = 10, double minBinSize = 500000, double expectedGenomeSize = 4_000_000)`
returns a set of `GenomeBin` records (`BinId, ContigIds, TotalLength, GcContent, Coverage,
Completeness, Contamination, PredictedTaxonomy`). An opt-in companion pair,
`CalculateTetranucleotideZScores` / `TetranucleotideZScoreCorrelation`, exposes the published **TETRA
z-score signature** used to compare tetranucleotide usage between fragments.

## The binning algorithm (k-means on a composite feature distance)

Each contig is reduced to three features, and contigs are clustered by **k-means** over a composite
distance:

```
composite distance(a,b) = |ΔGC| + |Δcoverage| + TNF_PearsonDistance(a,b)
```

1. **Features:** GC content, tetranucleotide frequency (TNF) vector, coverage (normalized).
2. **Composite distance:** sum of the absolute GC difference, the absolute normalized-coverage
   difference, and the **Pearson distance between the two TNF vectors**. All three contribute — TNF
   is no longer computed-but-ignored (see *Deviations resolved*).
3. **Centroid initialization by GC-sorted spread** — picks diverse starting centroids deterministically
   (no random seed), so binning is reproducible.
4. **Iterate** assignment/update to convergence, capped at **50 iterations**.
5. **Filter** bins whose `TotalLength` is below `minBinSize` (default 500 kb).

Composition and coverage are complementary: two organisms with **similar GC** can still be separated
if their coverage differs, and vice-versa — which is exactly why real binners (MetaBAT2, MaxBin2,
CONCOCT) combine both. When both signals coincide (similar GC *and* similar coverage), contigs from
different genomes can be **mis-binned together** — a documented failure mode.

## Quality metrics: completeness & contamination

Each bin carries two quality proxies, both clamped to `[0,100]`:

```
Completeness  = min(totalLength / expectedGenomeSize × 100, 100)     (default expectedGenomeSize = 4 Mbp)
Contamination = min(gcStdDev / 0.5 × 100, 100)                        (0.5 = theoretical max GC std-dev)
```

- **Completeness** is the bin's length as a fraction of the expected genome size (parameterized, not a
  hardcoded constant). A 2 MB bin against a 4 MB expectation → **50.0** (oracle M5).
- **Contamination** is the GC-content dispersion within the bin, normalized by the theoretical maximum
  GC std-dev of 0.5: uniform GC → **0**, maximal GC variance → **100** (oracle M6).

> **GOTCHA — the `GenomeBin.Completeness`/`Contamination` fields are still proxies; CheckM marker-gene
> QC is now a SEPARATE opt-in path.** The default `BinContigs` output above remains the length-ratio /
> GC-variance **proxies** — useful ordering signals, not marker calls. The published CheckM
> (Parks 2015) single-copy marker-gene metrics are now **built and validated** (test unit
> [[meta-bin-001-markerqc-evidence|META-BIN-001 MarkerQC addendum]]) but exposed through a
> **distinct API** (`EstimateBinQualityFromMarkerCounts` / `EstimateBinQualityFromMarkers` /
> `DetectMarkers`, see *Marker-gene QC* below), NOT wired into `BinContigs`. So a `GenomeBin` from
> `BinContigs` still carries proxy values; to get CheckM-style completeness/contamination you run the
> marker-QC path over a bin's contigs with a marker HMM set.

## Marker-gene QC — CheckM completeness & contamination (opt-in)

The MarkerQC addendum ([[meta-bin-001-markerqc-evidence]]) implements the published **CheckM**
(Parks et al. 2015) genome-quality metrics over **collocated single-copy marker sets** `M`, verbatim
against the CheckM reference `MarkerSet.genomeCheck`:

```
Completeness  = 100 · ( Σ_{s∈M} |s ∩ G_M| / |s| ) / |M|      (CheckM Eq. 1)
Contamination = 100 · ( Σ_{s∈M} Σ_{g∈s} C_g / |s| ) / |M|    (CheckM Eq. 2), C_g = N−1 for a gene seen N≥1×, else 0
```

A **multi-copy marker** counts **once** toward `present` (completeness uses set intersection, not the
copy count) and contributes `N−1` to contamination. Completeness ≈ fraction of unique single-copy
genes present; contamination ≈ how many are present in multiple copies. Grouping into marker *sets*
down-weights correlated, consistently-collocated genes.

- **API:** `EstimateBinQualityFromMarkerCounts` (feed a marker→copy-count map directly),
  `EstimateBinQualityFromMarkers` (over an end-to-end HMM-detected count), `DetectMarkers` (HMM search).
  Profiles load via `LoadMarkerHmms` (caller-supplied) or bundled loaders.
- **Bundled marker sets (licence-gated):** the 9 universal ribosomal Pfams (S2/S7/S8/S9/S10/S11/S19/
  L1/L3, Xu 2022) plus the **CC0 Pfam subsets** of GTDB's domain-level universal sets —
  `LoadBundledBacterialMarkerHmms`/`BundledBacterialMarkerSets` (**bac120**, 6 Pfams) and
  `LoadBundledArchaealMarkerHmms`/`BundledArchaealMarkerSets` (**ar122**, 35 Pfams). Each bundled
  family is treated as its own **singleton set** (|s|=1). **TIGRFAM-defined** bac120/ar122 members are
  **CC BY-SA 4.0 (share-alike) and deliberately NOT bundled** — caller supplies them via
  `LoadMarkerHmms`; only the CC0 Pfam subsets ship (39 distinct Pfam accessions across both bundles).
- **Detection gate:** a marker is present iff its **glocal Plan7 Viterbi** bit score ≥ the Pfam
  **GA1** per-sequence gathering threshold. The GA1 gate is sourced from the HMM file, but this
  engine's whole-sequence Viterbi differs from HMMER's local + null2-corrected `hmmsearch`, so
  absolute bit scores diverge (documented engine difference); the true-positive separation is
  decisive (E. coli uS8 vs PF00410 = +176 bits vs all other ribosomal families < 0).
- **Worked oracles:** synthetic 3-set bin → **Completeness = 250/3 ≈ 83.333%, Contamination =
  100/9 ≈ 11.111%** (pins Eqs. 1–2); one bundled family detected out of 9/6/35 singleton sets →
  **100/9 ≈ 11.111% / 100/6 ≈ 16.667% / ≈2.857%** respectively, contamination 0.
- **Guards:** empty marker set excluded (no div-by-zero); `|M|=0` → both metrics 0; missing marker
  contributes 0.
- **Residual:** CheckM's operon-based **collocation grouping** of markers into non-singleton sets
  needs the full lineage DB — not built (bundled families stay singleton sets); the QC path is also
  not auto-wired into `BinContigs` (call it explicitly over a bin's contigs).

## The TETRA z-score tetranucleotide signature (opt-in)

`CalculateTetranucleotideZScores` implements the **TETRA** method (Teeling et al. 2004) for a
genome-fragment compositional signature. The sequence is first **extended by its reverse-complement**
(`seq + revcomp(seq)`) so the signature is strand-symmetric, then over the extended strand:

```
E(n1n2n3n4)   = N(n1n2n3)·N(n2n3n4) / N(n2n3)              maximal-order (2nd-order) Markov expectation
var(n1n2n3n4) = E · [N(n2n3)−N(n1n2n3)] · [N(n2n3)−N(n2n3n4)] / N(n2n3)²      (Schbath 1997 approximation)
z(n1n2n3n4)   = (N(n1n2n3n4) − E) / √var
```

Two fragments are compared by the **Pearson correlation of their 256-dimensional z-vectors**
(`TetranucleotideZScoreCorrelation`); because correlation is **scale-invariant**, the absolute z
magnitudes (which scale with length) drop out and only the *pattern* of over-/under-representation
matters. Self-correlation = 1.0.

- **Worked oracle:** `ACGTACGTGGCC` RC-extends to `ACGTACGTGGCCGGCCACGTACGT` (24 nt); with N(ACGT)=4,
  N(ACG)=N(CGT)=4, N(CG)=5 → E=3.2, var=0.128, **z(ACGT)=√5=2.2360679774997896**.
- **Guards:** absent middle dinucleotide `N(n2n3)=0` (E undefined) → **z=0**; non-positive variance →
  **z=0**; degenerate signature → correlation **r=0**, not NaN.

## Invariants and edge cases

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `BinId`s are unique | one id per emitted bin |
| INV-02 | `Completeness ∈ [0,100]`, `Contamination ∈ [0,100]` | both formulas are `min(…, 100)` clamped, non-negative |
| INV-03 | `GcContent ∈ [0,1]` | `(G+C)/(A+T+G+C)` (oracle M7: pure GC→1.0, pure AT→0.0) |
| INV-04 | bin `TotalLength` = Σ contig lengths in the bin | additive aggregation |
| INV-05 | **partition** — no contig appears in two bins | k-means hard assignment (oracle M11) |
| INV-06 | bin `Coverage` = arithmetic mean of member coverages | oracle M12: mean(20,30,25)=25.0 |

- **Empty input → empty output.** A single contig may form its own bin or be size-filtered.
- **Bins below `minBinSize` are excluded** (default 500 kb), so short spurious clusters are dropped.
- **Determinism:** GC-sorted centroid initialization + fixed iteration cap ⇒ reproducible bins for a
  given input (no random seed).

## Scope and limitations

A [[research-grade-limitations|research-grade]] binner: k-means over composition+coverage with proxy
quality metrics, matching the compositional/coverage feature families the literature documents
(Wikipedia binning, TETRA). Known limits, several source-documented (Maguire 2020): chromosomes bin
well (82–94% recovered) but **plasmids (1–29%) and genomic islands (38–44%) are poorly recovered**
because their copy number and composition diverge from the host chromosome; organisms with **similar
GC and coverage** co-bin; high-contamination bins mix genomes. The **CheckM marker-gene QC is now
built** (opt-in path — see *Marker-gene QC* above; the default `BinContigs` metrics stay proxies), but
only the **CC0 Pfam** subsets of the universal / bac120 / ar122 marker sets are bundled (each as a
singleton set — no lineage-specific collocation grouping), and TIGRFAM-defined markers must be
caller-supplied. There is still no supervised taxonomic assignment, no differential-coverage model
across many samples, and no DAS-Tool-style bin refinement.

## Relation to siblings

- **Metagenomics family sibling of [[alpha-diversity]] and [[beta-diversity]].** Those quantify
  *community* structure — within-sample diversity indices and between-sample dissimilarity over a
  taxon→abundance profile — whereas binning *reconstructs the genomes themselves* from contigs. Same
  domain (`MetagenomicsAnalyzer`, mixed-community sequencing), different question: who is there / how
  different are two samples, vs. what are their genomes.
- **Tetranucleotide vs. k-mer composition.** The TETRA signature is a specialized 4-mer (256-word)
  compositional fingerprint with a Markov-expectation z-score normalization — related in spirit to the
  general k-mer composition captured by the K-mer family, but distinct in its expectation model and
  RC-extension.
</content>
