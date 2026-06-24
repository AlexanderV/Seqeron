# Validation Report: META-BIN-001 — Genome Binning (MAG assembly)

- **Validated:** 2026-06-24   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.BinContigs(IEnumerable<(string,string,double)>, int numBins, double minBinSize, double expectedGenomeSize)`
  (helpers: `KMeansCluster`, `CompositeDistance`, `TnfPearsonDistance`, `CalculateTetraNucleotideFrequency`, `CalculateGcContent`, `CalculateContamination`)
  in `src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs:647-861`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **End state:** CLEAN (no defect; the simplifications are honestly declared in spec, evidence, and code)

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — Binning (metagenomics)**: binning groups assembled contigs into MAGs using
  compositional features (GC content, tetranucleotide/k-mer frequencies) and/or read-mapping
  coverage across samples; supervised or unsupervised classification. Confirms the feature set
  the spec claims.
- **MetaBAT (Kang et al., PeerJ 2015)**: "Metagenome Binning based on Abundance and
  Tetranucleotide frequency" — computes per-contig-pair distances from TNF and mean-base
  coverage, then integrates them into one **composite distance**. This is exactly the
  composite-distance paradigm the implementation uses.
- **CONCOCT**: clusters contigs on combined sequence composition + coverage features (PCA + GMM).
  Confirms composition+coverage is the standard de-novo binning feature space.
- **TETRA (Teeling et al. 2004)**: tabulates all **256 (4⁴)** tetramer frequencies, derives
  expected frequencies via a maximal-order Markov model, converts observed-vs-expected to
  **z-scores**, and compares two sequences by the **Pearson correlation of their z-score
  vectors**.
- **CheckM (Parks et al. 2014, Genome Research 25:1043)**: completeness = fraction of
  lineage-specific single-copy marker genes (SCGs) present; contamination = fraction of those
  SCGs found in multiple copies. Both are **marker-gene** based, not length/composition based.

### Method the spec defines
Features per contig: GC content, TNF vector, coverage. Similarity via composite distance
(|ΔGC| + |Δcoverage_normalized| + TNF-Pearson distance); clustering via deterministic k-means.
Quality reported as completeness (length ratio) and contamination (within-bin GC std dev).
The feature set and composite-distance + clustering structure match the literature.

### Honest-scope assessment (two declared simplifications — not defects)
1. **TNF distance uses raw 4-mer relative frequencies, not Teeling z-scores.** Pearson is taken
   directly on observed frequency vectors, omitting the Markov-model expected-frequency / z-score
   normalisation of full TETRA. The XML-doc comment says "per TETRA methodology"; this is a
   simplified Pearson-on-frequencies variant. Declared in the Evidence "Deviations" note.
2. **CheckM metrics are proxies.** Completeness = `min(len/expectedGenomeSize×100, 100)` (length
   ratio); contamination = `min(gcStdDev/0.5×100, 100)` (within-bin GC dispersion). Neither uses
   single-copy marker genes. Stated explicitly in the TestSpec, Evidence doc, and code comments.

Because these are surfaced in code comments + TestSpec + Evidence, the scope is honest rather than
misadvertised.

### Feature dimensionality hand-check
TNF space = 256 (4⁴), matching Teeling/Wikipedia. Code stores TNF sparsely as
`Dictionary<string,double>` keyed by the literal ACGT 4-mer (non-ACGT 4-mers skipped); support is a
subset of the 256-space — dimensionally correct (full 256-space, no reverse-complement collapsing
to 136, as the spec claims).

### Findings
Stage A passes with the two declared simplifications noted. Feature set, 256-dimension TNF space,
Pearson-correlation similarity, composite-distance binning, and GC-separation behaviour are all
correct against authoritative sources.

## Stage B — Implementation

### Code path reviewed
`MetagenomicsAnalyzer.cs:647-861` — `BinContigs`, `KMeansCluster`, `CompositeDistance`,
`TnfPearsonDistance`, `CalculateTetraNucleotideFrequency`, `CalculateGcContent`,
`CalculateContamination`.

### Formula realised correctly?
- **GC content** = `CalculateGcFractionFast` = (G+C)/(A+C+G+T) ∈ [0,1]; empty → 0. ✓
- **TNF** = sliding 4-mer counts over ACGT, normalised to relative frequencies. ✓ (256-space subset)
- **TNF Pearson distance** = (1 − r)/2 over union of keys (r=0 when denom undefined; 1.0 when no
  keys). ✓
- **Composite distance** = |ΔGC| + |Δcoverage_normalized| + TNF distance (coverage normalised to
  [0,1] by max, guarded ≥1). ✓
- **k-means**: deterministic centroid init by GC-sorted spread (`sortedIndices[c*n/k]`),
  assignment/update steps, max 50 iters, early stop on no change; empty clusters skip update. ✓
  deterministic and order-stable.
- **Completeness** = `min(totalLength/expectedGenomeSize×100, 100)`. ✓
- **Contamination** = `min(stdDev/0.5×100, 100)`, population std dev (÷count), 0 for <2 contigs. ✓

### Cross-verification table (independently hand-computed this session vs code/tests)
| Quantity | Input | Hand value | Test | Result |
|---|---|---|---|---|
| TNF Pearson dist (GC vs AT contig) | keys{GCGC,CGCG,ATAT,TATA}, va=[.5,.5,0,0], vb=[0,0,.5,.5] → r=−1 → (1−(−1))/2 | 1.0 | M10 | ✓ |
| GC pure-high / pure-low | GCGC… / ATAT… | 1.0 / 0.0 | M7 | ✓ |
| Coverage mean | {20,30,25} | 25.0 | M12 | ✓ |
| Completeness | 2 Mb / 4 Mb × 100 | 50.0 | M5 | ✓ |
| Contamination uniform | all same GC → σ=0 | 0.0 | M6 | ✓ |
| Contamination max-variance | 10×GC=1, 10×GC=0 → mean .5, var (20·.25)/20=.25, σ=.5 → .5/.5·100 | 100.0 | M6 | ✓ |

### Edge cases
- Empty input → `yield break` (M1). ✓
- Single / below-min-size → size filtered, empty (M2, M9). ✓
- Identical contigs → cluster together / few bins (C1). ✓
- Extreme GC → separate bins, no mixed bin (M10). ✓
- Very short (<4 bp) sequences → empty TNF, no throw, size filtered (C2). ✓
- Zero coverage → `maxCov` guarded to 1, no div-by-zero (ZeroCoverage test). ✓

### Variant/delegate consistency
Single public entry point; helpers internal and mutually consistent. No `*Fast` variant to diverge.

### Test quality audit
18 tests assert exact sourced values (GC 1.0/0.0, coverage 25.0, completeness 50.0,
contamination 0/100), disjointness (M11), unique bin IDs (M4), length accounting (M8), and the
k-means k-bound (S2) — not mere "no throw". Deterministic centroid init makes assertions
reproducible.

### Findings / defects
No defect. The only divergences from full TETRA/CheckM are the documented simplifications
(z-score omission; proxy completeness/contamination), declared in code comments, TestSpec, and
Evidence — honest scope, not a defect.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES — features/dimensionality/similarity/composite-distance correct;
  TNF z-score step and marker-gene quality metrics simplified but honestly declared.
- **Stage B:** PASS-WITH-NOTES — code faithfully realises the (simplified) description; all
  hand-computed cross-checks match; edge cases handled.
- **State:** CLEAN. No code changes. Tests: GenomeBinning = 18/18 pass (this session);
  prior full `Seqeron.Genomics.Tests` baseline = 4484/0.
- **Follow-up (optional, not a defect):** for MetaBAT/CheckM-grade fidelity, add z-score
  normalisation (Markov expected frequencies) to TNF and replace the length/GC proxies with
  single-copy-marker-gene completeness/contamination.
