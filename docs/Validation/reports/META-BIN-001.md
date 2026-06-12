# Validation Report: META-BIN-001 — Genome Binning (MAG assembly)

- **Validated:** 2026-06-12   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.BinContigs(IEnumerable<(string,string,double)>, int numBins, double minBinSize, double expectedGenomeSize)`
  (helpers: `KMeansCluster`, `CompositeDistance`, `TnfPearsonDistance`, `CalculateTetraNucleotideFrequency`, `CalculateGcContent`, `CalculateContamination`)
  in `src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **End state:** CLEAN (no defect; scope is honestly declared)

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — Binning (metagenomics)** (https://en.wikipedia.org/wiki/Binning_(metagenomics)):
  Compositional binning groups contigs into MAGs using GC-content, tetranucleotide/k-mer
  frequencies, and coverage across samples. Confirms tetranucleotide space is
  **4⁴ = 256** ("there can be 4⁴ = 256 different fragments of four consecutive nucleotides").
- **Teeling et al. (2004), TETRA** (BMC Bioinformatics, doi:10.1186/1471-2105-5-163;
  Env. Microbiol. doi:10.1111/j.1462-2920.2004.00624.x): counts all **256** tetranucleotides,
  derives expected frequencies via a maximal-order (di-/trinucleotide) Markov model, converts
  observed-vs-expected divergence to **z-scores**, then compares two sequences by the
  **Pearson correlation coefficient of their z-scores**. Confirms correlation > GC-content alone
  in discriminatory power.
- **Parks et al. (2014), CheckM** (Genome Research 25:1043): completeness = fraction of
  lineage-specific single-copy marker genes (SCGs) present; contamination = fraction of those
  SCGs found in multiple copies. Both are **marker-gene** based, not composition based.
- **Maguire et al. (2020)** (doi:10.1099/mgen.0.000436): 82–94% chromosomes recovered, but
  GIs/plasmids poorly recovered — establishes that real binning is imperfect and feature-driven.

### Method the spec defines
Bin contigs by compositional + coverage signature similarity: features = GC content,
tetranucleotide frequency (TNF) vector, coverage; similarity via composite distance
(|ΔGC| + |Δcoverage| + TNF Pearson distance); clustering via k-means. Quality reported as
completeness (length ratio) and contamination (GC std dev). This matches the standard binning
**feature set** in the literature.

### Honest-scope assessment (important)
The spec/Evidence declare this a **simple composition-distance k-means clusterer**, NOT a
MetaBAT2/CONCOCT/MaxBin-grade binner. Two simplifications relative to the cited sources are
genuine but **explicitly documented**, so the scope is honest:
1. **TNF distance uses raw 4-mer frequencies, not Teeling z-scores.** The code computes Pearson
   correlation directly on observed frequency vectors, omitting the Markov-model expected-frequency
   normalisation/z-score step of TETRA. (Comment says "per TETRA methodology"; this is a
   simplified Pearson-on-frequencies variant, not full TETRA.)
2. **CheckM metrics are proxies.** Completeness = `min(len/expectedGenomeSize×100, 100)` (length
   ratio), contamination = `min(gcStdDev/0.5×100, 100)` (within-bin GC dispersion). Neither uses
   single-copy marker genes. The Evidence doc states this directly.

These are surfaced in the source XML-doc comments and the Evidence "Deviations" section, so they
do not constitute a misadvertised capability.

### Feature dimensionality hand-check
TNF space is **256** (4⁴), matching Teeling/Wikipedia. The implementation stores TNF as a sparse
`Dictionary<string,double>` keyed by the literal 4-mer over the ACGT alphabet (non-ACGT 4-mers
skipped), so its support is a subset of the canonical 256 — dimensionally correct (no
reverse-complement collapsing to 136; full 256-space, as the spec/Evidence claim).

### Worked example (separation)
- Contig A `GCGC…` (100% GC): only 4-mers `GCGC`,`CGCG` each freq 0.5.
- Contig B `ATAT…` (0% GC): only `ATAT`,`TATA` each 0.5.
- Union keys = {GCGC,CGCG,ATAT,TATA}; va=[.5,.5,0,0], vb=[0,0,.5,.5]; mean=0.25.
  Pearson r = −1 ⇒ TNF distance = (1−(−1))/2 = **1.0**; plus |ΔGC| = 1.0 ⇒ large composite
  distance ⇒ they bin apart. Near-identical profiles (r≈1, ΔGC≈0) ⇒ distance≈0 ⇒ bin together.
  Matches the qualitative claim in the prompt and test M10.

### Findings
Stage A passes with the two declared simplifications noted (TNF z-score omission; proxy quality
metrics). The 256-dimension claim, feature set, Pearson-correlation similarity concept, and
GC-separation behaviour are all correct against authoritative sources.

## Stage B — Implementation

### Code path reviewed
`MetagenomicsAnalyzer.cs:464-678` — `BinContigs`, `KMeansCluster`, `CompositeDistance`,
`TnfPearsonDistance`, `CalculateTetraNucleotideFrequency`, `CalculateGcContent`,
`CalculateContamination`.

### Formula realised correctly?
- **GC content** = `CalculateGcFractionFast` = (G+C)/(A+C+G+T) ∈ [0,1]. ✓
- **TNF** = sliding 4-mer counts over ACGT, normalised to relative frequencies. ✓ (256-space subset)
- **TNF Pearson distance** = (1 − r)/2 over the union of keys (r=0 when undefined). ✓ matches worked example.
- **Composite distance** = |ΔGC| + |Δcoverage(normalised)| + TNF distance. ✓ (coverage normalised to [0,1] by max).
- **k-means**: deterministic GC-sorted centroid spread, assignment/update, max 50 iters, early stop. ✓ deterministic.
- **Completeness** = `min(totalLength/expectedGenomeSize×100, 100)`. ✓
- **Contamination** = `min(gcStdDev/0.5×100, 100)`, 0 for <2 contigs. ✓ (population std dev).

### Cross-verification table (hand-computed vs code/tests)
| Quantity | Input | Hand value | Test | Result |
|---|---|---|---|---|
| TNF Pearson dist (GC vs AT contig) | GCGC… vs ATAT… | 1.0 | M10 separation | ✓ |
| GC pure-high / pure-low | GCGC… / ATAT… | 1.0 / 0.0 | M7 | ✓ |
| Coverage mean | {20,30,25} | 25.0 | M12 | ✓ |
| Completeness | 2 Mb / 4 Mb | 50.0 | M5 | ✓ |
| Contamination uniform | all same GC | 0.0 | M6 | ✓ |
| Contamination max-variance | 10×GC=1, 10×GC=0 → σ=0.5 | 100.0 | M6 | ✓ |

### Edge cases
- Empty input → `yield break` (M1). ✓
- Single/below-min-size → filtered out, empty (M2, M9). ✓
- Identical contigs → cluster together / few bins (C1). ✓
- Dissimilar (extreme GC) → separate bins (M10). ✓
- Very short sequences (<4 bp) → empty TNF, no throw, size-filtered (C2). ✓
- Zero coverage → `maxCov` guarded to 1, no div-by-zero (ZeroCoverage test). ✓

### Variant/delegate consistency
Single public entry point; helpers are internal and consistent. No `*Fast` variant to diverge.

### Test quality audit
18 tests assert exact sourced values (GC 1.0/0.0, coverage 25.0, completeness 50.0,
contamination 0/100), disjointness, uniqueness, length accounting, and k-means k-bound — not
mere "no throw". Deterministic centroid init makes assertions reproducible.

### Findings / defects
No defect. The only divergences from full TETRA/CheckM are the documented simplifications
(z-score omission; proxy completeness/contamination) — these are declared in code comments,
the TestSpec, and the Evidence doc, so they are honest scope rather than defects.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES — features/dimensionality/similarity correct; TNF z-score step and
  marker-gene quality metrics simplified but honestly declared.
- **Stage B:** PASS-WITH-NOTES — code faithfully realises the (simplified) description; all
  hand-computed cross-checks match; edge cases handled.
- **State:** CLEAN. No code changes required. Tests: 18/18 GenomeBinning pass; full
  `Seqeron.Genomics.Tests` = 4484 passed, 0 failed.
- **Follow-up (optional, not a defect):** if MetaBAT-grade fidelity is ever required, add
  z-score normalisation (Markov expected frequencies) to TNF and replace the length/GC proxies
  with single-copy-marker-gene completeness/contamination.
