---
type: source
title: "Evidence: META-BIN-001 (metagenomic binning — k-means MAGs + TETRA z-score signature)"
tags: [validation, metagenomics]
doc_path: docs/Evidence/META-BIN-001-Evidence.md
sources:
  - docs/Evidence/META-BIN-001-Evidence.md
source_commit: 0627f8b8dd6ecb12eedb6143693e51bfe06dac31
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: META-BIN-001

The validation-evidence artifact (Document Version 3.0) for test unit **META-BIN-001** —
**metagenomic binning**, the grouping of assembled contigs into Metagenome-Assembled Genomes
(MAGs) computed by `MetagenomicsAnalyzer.BinContigs`, plus a 2026-06-24 **addendum** covering the
opt-in **TETRA z-score tetranucleotide signature** (`CalculateTetranucleotideZScores` /
`TetranucleotideZScoreCorrelation`). Third ingested unit of the Metagenomics family, alongside the
diversity siblings [[meta-alpha-001-evidence|META-ALPHA-001]] and
[[meta-beta-001-evidence|META-BETA-001]]. One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The method is synthesized in its own
concept, [[metagenomic-binning]]; [[test-unit-registry]] tracks the unit.

## What this file records

Binning groups contigs to their genome of origin using **compositional** and **coverage** features,
scored for quality with completeness/contamination metrics. See `docs/Evidence/META-BIN-001-Evidence.md`.

- **Online sources (mutually consistent):**
  - **Wikipedia — Binning (metagenomics):** binning groups assembled contigs to their genomes of
    origin (results = **MAGs**); uses compositional features (GC content, tetranucleotide
    frequencies) and cross-sample coverage; supervised or unsupervised. Names real tools
    TETRA / MEGAN / PhyloPythia / MetaBAT2 / MaxBin2 / CONCOCT / DAS Tool.
  - **TETRA (Teeling et al. 2004, BMC Bioinformatics 5:163):** tetranucleotide-usage classification;
    256 tetramers (4⁴); z-scores flag over-/under-representation vs. a maximal-order Markov
    expectation; sequences RC-extended for strand symmetry; sequences compared by **Pearson
    correlation of their z-score vectors**. Variance uses the **Schbath (1997)** word-count approximation.
  - **CheckM (Parks et al. 2014/2015, Genome Res 25(7):1043):** genome-quality assessment —
    **completeness** from lineage-specific single-copy marker genes present, **contamination** from
    duplicated marker genes, strain heterogeneity indicator. HMMER v3.1b1 with Pfam `-cut_gc` /
    TIGRFAM `-cut_nc` cutoffs.
  - **Maguire et al. (2020, Microb Genom):** MAG-binning limitations — 82–94% chromosomes recovered,
    only 38–44% of genomic islands, only 1–29% of plasmids; variable copy number/composition is
    problematic.

- **Documented feature families:** compositional (GC content, tetranucleotide frequency, codon
  usage); coverage (read depth across samples, differential-coverage binning); quality
  (completeness, contamination, strain heterogeneity).

## Implementation notes (from the Evidence file)

`BinContigs(IEnumerable<(string ContigId, string Sequence, double Coverage)> contigs, int numBins = 10,
double minBinSize = 500000, double expectedGenomeSize = 4_000_000)` runs **k-means on a composite
feature distance**:

1. Per-contig features: **GC content**, **tetranucleotide frequency (TNF)**, **normalized coverage**.
2. Composite distance `= |ΔGC| + |Δcoverage| + TNF Pearson distance` (all three features contribute —
   TNF is no longer computed-but-unused).
3. Centroid initialization by **GC-sorted spread** → deterministic, diverse starting points.
4. Iterative assignment/update to convergence (**max 50 iterations**).
5. Filter bins below **minBinSize**.
6. **Completeness** `= min(totalLength / expectedGenomeSize × 100, 100)` — parameterized genome size,
   no hardcoded constant.
7. **Contamination** `= min(gcStdDev / 0.5 × 100, 100)` — GC std-dev normalized by the theoretical
   maximum GC std-dev (0.5).

Output is a `GenomeBin` record: `BinId, ContigIds, TotalLength, GcContent, Coverage, Completeness,
Contamination, PredictedTaxonomy`.

### TETRA z-score signature addendum (2026-06-24)

The opt-in signature counts overlapping di/tri/tetra words over the **RC-extended strand**
`seq + revcomp(seq)` (strand symmetry). Formulas (Teeling 2004 / Schbath 1997, reproduced from
open-access sources):

```
E(n1n2n3n4)   = N(n1n2n3)·N(n2n3n4) / N(n2n3)                 (maximal-order / 2nd-order Markov)
var(n1n2n3n4) = E · [N(n2n3)−N(n1n2n3)] · [N(n2n3)−N(n2n3n4)] / N(n2n3)²
z(n1n2n3n4)   = (N(n1n2n3n4) − E) / √var
comparison    = Pearson correlation of the two 256-dim z-vectors     (scale-invariant)
```

**Hand-derived oracle:** for `ACGTACGTGGCC`, RC-extended to `ACGTACGTGGCCGGCCACGTACGT` (24 nt),
N(ACGT)=4, N(ACG)=N(CGT)=4, N(CG)=5 → E=3.2, var=0.128, **z(ACGT)=√5=2.2360679774997896**.

## Design decisions and deviations

- **Deviations from literature: None.** Three previously documented assumptions were all resolved:
  completeness was a length ratio to a hardcoded 2 Mbp genome → now parameterized via
  `expectedGenomeSize` (default 4 Mbp); contamination used an arbitrary GC-variance formula → now GC
  std-dev normalized by the theoretical maximum (0.5); TNF was computed but unused → now included in
  the k-means composite distance as a Pearson distance term.
- **CheckM marker-gene QC is an honest residual (NOT implemented).** CheckM's lineage-specific
  single-copy marker sets are a large trained database (`checkm_data`: Pfam/TIGRFAM HMM libraries +
  a precomputed reference-genome tree, ~GB-scale) not cleanly retrievable as plaintext. Per the STOP
  RULE the evidence file implements **only** the z-score signature; the marker-gene completeness/
  contamination is left explicitly unbuilt rather than fabricated. The `Completeness`/`Contamination`
  the code returns are the **length-ratio / GC-variance proxies**, not CheckM marker calls.
- **z-score guards:** absent middle dinucleotide `N(n2n3)=0` → E undefined → **z=0**; non-positive
  variance (`N(n2n3)−N(n1n2n3) ≤ 0` etc.) → **z=0**; degenerate signature → correlation `r=0` not NaN.
  RC-extension means even a ≥2-base input yields a signature.

## Source-verified invariants and oracles

**Invariants:** unique `BinId`s; `Completeness ∈ [0,100]`; `Contamination ∈ [0,100]`;
`GcContent ∈ [0,1]`; a bin's `TotalLength` = Σ its contig lengths; empty input → empty output; bins
below `minBinSize` excluded (partition — no contig in two bins).

**Worked oracles (evidence test table):**

| Test | Oracle | Basis |
|------|--------|-------|
| M5 | Completeness = 50.0 for a 2 MB bin / 4 MB genome | `min(L/E×100, 100)` |
| M6 | Contamination = 0 (uniform GC), = 100 (max variance) | `stddev/0.5×100` |
| M7 | GC = 1.0 (pure G/C), = 0.0 (pure A/T) | `(G+C)/(A+T+G+C)` |
| M10 | Extreme-GC populations separate into distinct bins | k-means on compositional distance |
| M11 | Bins are disjoint (no contig in two bins) | partition invariant |
| M12 | Coverage = 25.0 = mean(20,30,25) | arithmetic mean |
| TETRA | z(ACGT) = √5 for `ACGTACGTGGCC`; self-correlation = 1.0 | Markov E/var + Pearson |

## Edge cases from literature

Empty input → empty results; single contig → own bin or size-filtered; minimum-bin-size filtering;
similar GC content across organisms → incorrect co-binning; variable-coverage regions → assembly/
binning issues; high-contamination bins hold multiple genomes. No source contradictions — the
Wikipedia / Teeling 2004 / CheckM primary definitions are mutually consistent; the only honest gap
is the un-built CheckM marker-gene QC.
