# Validation Report: META-BIN-001 — Genome Binning (MAG assembly)

- **Validated:** 2026-06-24   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.BinContigs(IEnumerable<(string,string,double)>, int numBins, double minBinSize, double expectedGenomeSize)`
  (helpers: `KMeansCluster`, `CompositeDistance`, `TnfPearsonDistance`, `CalculateTetraNucleotideFrequency`, `CalculateGcContent`, `CalculateContamination`)
  in `src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs:647-861`
- **Stage A verdict:** PASS-WITH-NOTES (re-confirmed 2026-06-25 — see fresh re-validation addendum below)
- **Stage B verdict:** PASS-WITH-NOTES (re-confirmed 2026-06-25 — see fresh re-validation addendum below)
- **End state:** CLEAN (no defect; the simplifications are honestly declared in spec, evidence, and code).
  Re-validated fresh on 2026-06-25 after the TETRA z-score signature (→ META-TETRA-001) and CheckM
  marker QC (→ META-CHECKM-001) were split out as separate units. This unit's OWN canonical surface
  is the TNF+coverage composite-distance k-means binner (`BinContigs` + its helpers). One test-coverage
  gap fixed (differential-coverage separation), no code change. CLEAN.

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

## Update 2026-06-24 — z-score-normalised TETRA signature added (opt-in)

The first follow-up above is now implemented as an **opt-in** API (the default raw-frequency
binning path is unchanged):

- **New methods** (`MetagenomicsAnalyzer.cs`):
  `CalculateTetranucleotideZScores(string)` → 256-component TETRA z-score signature;
  `TetranucleotideZScoreCorrelation(string, string)` → Pearson correlation of two z-score vectors.
- **Method realised verbatim from Teeling et al. (2004) / Schbath (1997):** reverse-complement
  strand extension; maximal-order (2nd-order) Markov expected count
  `E(n1n2n3n4) = N(n1n2n3)·N(n2n3n4)/N(n2n3)`; Schbath variance
  `var = E·[N(n2n3)−N(n1n2n3)][N(n2n3)−N(n2n3n4)]/N(n2n3)²`; `z = (N−E)/√var`; signatures compared
  by Pearson correlation of the 256-component z-score vectors.
- **Evidence retrieved this session:** TETRA BMC paper (PMC529438), the expected-count form from an
  open-access reproduction (PLOS ONE pone.0008113), the full z/variance form from independent
  search + a reference implementation, and the Schbath 1997 primary. URLs recorded in the Evidence
  addendum.
- **Hand-derived worked example asserted exactly:** for `ACGTACGTGGCC` (RC-extended to 24 nt),
  N(ACGT)=4, N(ACG)=4, N(CGT)=4, N(CG)=5 ⇒ E=3.2, var=0.128, **z(ACGT)=√5=2.2360679774997896**
  (asserted `Within(1e-10)`). Also asserted: self-correlation = 1.0; similar > dissimilar; z=0 when
  N(n2n3)=0; symmetry; degenerate→0 (not NaN). 9 tests in
  `MetagenomicsAnalyzer_TetranucleotideZScore_Tests.cs`.
- **CheckM single-copy-marker-gene QC — HONEST RESIDUAL (not implemented):** the CheckM
  lineage-specific Pfam/TIGRFAM HMM marker sets + reference genome tree are a large trained database
  (`checkm_data`, ~GB) with no clean plaintext source; per the no-fabrication rule the length-ratio
  completeness and GC-stddev contamination proxies are retained and the marker-gene QC is left as a
  declared residual. Use CheckM for marker-gene-based MAG quality.
- **Registry:** META-BIN-001 Status reset `☑`→`☐` in the root `ALGORITHMS_CHECKLIST_V2.md` for the
  re-validation pass (Completed 221→220, Not Started 13→14).

## Update 2026-06-25 — CheckM-style marker-gene completeness/contamination implemented (opt-in)

The prior residual is now implemented as an **opt-in** API (the default `BinContigs` TNF/proxy path
and its defaults are unchanged):

- **Formula retrieved verbatim** from Parks et al. (2015) *Genome Res* 25:1043 (open-access PMC4484387)
  and cross-confirmed against the CheckM reference implementation `MarkerSet.genomeCheck`
  (Ecogenomics/CheckM, `checkm/markerSets.py`). For collocated marker sets `M`:
  - Completeness = `100·(1/|M|)·Σ_{s∈M} |s∩G_M|/|s|` (Eq. 1).
  - Contamination = `100·(1/|M|)·Σ_{s∈M} Σ_{g∈s} C_g/|s|`, `C_g = N−1` for a marker found `N≥1`
    times, 0 if missing (Eq. 2). A multi-copy marker counts **once** toward completeness.
- **New methods** (`MetagenomicsAnalyzer.cs`): `EstimateBinQualityFromMarkerCounts(markerSets,
  markerCounts)` (pure formula), `DetectMarkers(proteins, markerHmms)` (Plan7 Viterbi bits ≥ GA1 ⇒
  copy count, reusing `Plan7ProfileHmm`), `EstimateBinQualityFromMarkers(...)` (detect + formula),
  `LoadBundledRibosomalMarkerHmms()` / `BundledRibosomalMarkerSets()` (9 bundled CC0 markers),
  `LoadMarkerHmms(readers, thresholds?)` (caller-supplied loader). New records `MarkerSet`,
  `MarkerHmm`, `BinMarkerQuality`.
- **Bundled marker set — CC0:** 9 universal single-copy ribosomal-protein Pfam HMMs (S2 PF00318,
  S7 PF00177, S8 PF00410, S9 PF00380, S10 PF00338, S11 PF00411, S19 PF00203, L1 PF00687, L3 PF00297),
  retrieved 2026-06-25 from the EMBL-EBI InterPro Pfam HMM API (verbatim HMMER3/f), licence CC0
  (public domain). Provenance + licence in `Resources/README.md` and the Evidence addendum.
- **Verified exactly:** hand-derived synthetic bin (`M={{A,B},{C,D,E},{F}}`, counts A1 B0 C2 D1 E1 F1)
  ⇒ Completeness = 250/3 % = 83.333…, Contamination = 100/9 % = 11.111… (asserted `Within(1e-10)`);
  triplicated marker ⇒ Contamination 200 %; HMM detection — bundled PF00410 detects E. coli uS8
  (UniProt P0A7W7) at ≈176 bits (GA1 24) and **none** of the 8 other ribosomal families do, so the
  end-to-end completeness over 9 singleton sets = 100/9 %. 14 tests in
  `MetagenomicsAnalyzer_MarkerGeneQuality_Tests.cs`.
- **HONEST RESIDUAL (narrowed):** the full CheckM `checkm_data` lineage-specific marker sets +
  reference-genome tree (tree-based lineage placement, operon-based marker-set collocation) are a
  large gated trained DB, not bundled — run CheckM itself or supply lineage-specific marker HMMs via
  `LoadMarkerHmms`.
- **Status:** META-BIN-001 remains `☐` (re-validation pass ongoing); Quick-Reference counts unchanged.
- **Tests:** marker-QC fixture 14/14 pass.

## Update 2026-06-25 — full domain-level universal marker set bundled (GTDB bac120/ar122 Pfam, CC0)

The prior "small 9-ribosomal set" is expanded to a **published domain-level universal single-copy
marker set**, additive and opt-in (the 9-ribosomal set, TNF/proxy `BinContigs`, and all defaults are
unchanged):

- **Marker set + source (verbatim):** the **Pfam-defined members of the GTDB `bac120` (Bacteria) and
  `ar122` (Archaea)** universal single-copy marker sets — Parks et al. (2018) *Nat Biotechnol*
  36:996, GTDB-Tk. Marker lists retrieved 2026-06-25 verbatim from the Ecogenomics/GTDBTk repo
  (`tests/.../gtdbtk.{bac120,ar122}.marker_info.tsv`) via the GitHub contents API. bac120 = 120
  markers (6 Pfam + 114 TIGRFAM); ar122 = 122 markers (35 Pfam + 87 TIGRFAM).
- **HMM provenance + LICENCE:** the **6 bac120 + 35 ar122 Pfam HMMs** (39 distinct, 3 already in the
  ribosomal set ⇒ 36 newly fetched) were retrieved 2026-06-25 from the EMBL-EBI InterPro Pfam HMM
  API (HMMER3/f ASCII, decompressed verbatim), licence **CC0 (public domain)** — embedded as
  resources. The **TIGRFAM-defined members (114 + 87)** are licensed **CC BY-SA 4.0** (TIGRFAMs at
  NCBI; reusabledata.org/tigrfams), confirmed **NOT public domain**, so per the no-non-redistributable
  rule they are **NOT bundled** — caller-supplied via `LoadMarkerHmms`. Hence the bundled domain sets
  are the **Pfam (CC0) subsets** of bac120/ar122. Provenance/licence in `Resources/README.md`.
- **New methods** (`MetagenomicsAnalyzer.cs`): `LoadBundledBacterialMarkerHmms()` (6),
  `BundledBacterialMarkerSets()` (6), `LoadBundledArchaealMarkerHmms()` (35),
  `BundledArchaealMarkerSets()` (35). The CheckM formula + Plan7 engine
  (`EstimateBinQualityFromMarkers`, `DetectMarkers`, `EstimateBinQualityFromMarkerCounts`) are reused
  unchanged.
- **Verified exactly (new tests):** bundled bac120 subset = the 6 expected Pfam ids, PF01025 GA1=25.8
  bits / LENG=165; ar122 subset = the 35 expected Pfam ids, PF00410 GA1=24. Real-marker detection:
  E. coli GrpE (UniProt P09372) hits **only** bac120 PF01025 (1 copy, 0 for the other 5) ⇒ bin
  completeness 100/6 % (≈16.667), 0 contamination; E. coli uS8 (P0A7W7) hits the universal PF00410 in
  the archaeal set ⇒ completeness 100/35 %, 0 contamination. Regression guard confirms the 9-ribosomal
  accessor is unchanged. 6 new tests (20 total in `MetagenomicsAnalyzer_MarkerGeneQuality_Tests.cs`).
- **HONEST RESIDUAL (narrowed):** the **per-lineage-specific** CheckM marker refinement + the
  **reference genome tree** for lineage placement (gated `checkm_data`; tree-based placement,
  operon-based collocation) are not bundled, plus the CC BY-SA TIGRFAM markers — all caller-supplied.
  Domain-level completeness/contamination now works out of the box.
- **Status:** META-BIN-001 remains `☐` (re-validation pass ongoing); Quick-Reference counts unchanged.
- **Tests:** marker-QC fixture 20/20 pass.

## Update 2026-06-25 — FRESH RE-VALIDATION of META-BIN-001's own binning surface (post-split)

This is the campaign re-validation requested after the unit was reset to ⬜ pending. The z-score
TETRA signature and the CheckM marker-gene QC are now their OWN units (META-TETRA-001,
META-CHECKM-001) and are **referenced but not re-litigated** here. META-BIN-001's canonical surface
is the **TNF + coverage composite-distance k-means binner**:
`BinContigs(IEnumerable<(ContigId, Sequence, Coverage)>, numBins, minBinSize, expectedGenomeSize)`
and its private helpers `KMeansCluster`, `CompositeDistance`, `TnfPearsonDistance`,
`CalculateTetraNucleotideFrequency`, `CalculateGcContent`, `CalculateContamination`
(`src/.../Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs:855-1069`).

### Stage A — description (re-confirmed against external sources retrieved THIS session)
- **MetaBAT (Kang et al., PeerJ 2015; peerj.com/articles/1165):** "Metagenome Binning based on
  Abundance and Tetranucleotide frequency" — per contig-pair distances from TNF and mean-base
  **coverage** are integrated into **one composite distance**, then clustered (modified k-medoid).
  This is exactly the composite-distance + clustering paradigm `BinContigs` implements.
- **CONCOCT (Alneberg et al., Nat Methods 2014; nature.com/articles/nmeth.3103):** clusters contigs
  on a feature vector of **tetra-mer composition + coverage across samples** (GMM). Confirms
  composition+coverage as the standard de-novo feature space and the **differential-coverage**
  signal (genomes with similar composition separated by distinct abundance).
- **Wikipedia — Binning (metagenomics):** "based on either compositional sequence features (such as
  **GC-content or tetranucleotide frequencies**) or sequence read mapping **coverage** … or both";
  TETRA compares tetramer-derived vectors **pair-wise**. Confirms `BinContigs`'s [GC, TNF, coverage]
  feature set and the pairwise correlation similarity.
- **Two declared simplifications (NOT defects, and not this unit's surface to fix):** the default
  TNF path uses raw 4-mer relative frequencies, not Teeling z-scores (full z-score signature =
  META-TETRA-001); completeness/contamination are length-ratio / GC-stddev proxies (marker-gene QC
  = META-CHECKM-001). Surfaced in code comments + TestSpec + Evidence → honest scope.
- **Stage A verdict: PASS-WITH-NOTES** — feature set, composite distance, pairwise-correlation TNF
  similarity, and the GC + differential-coverage separation behaviour are all correct against the
  primary binning literature.

### Stage B — implementation (re-read + independent cross-check THIS session)
- Code path re-read (`MetagenomicsAnalyzer.cs:855-1069`): `CompositeDistance = |ΔGC| + |Δcov_norm| +
  TnfPearsonDistance`; `TnfPearsonDistance = (1 − r)/2`; coverage normalised to [0,1] by max (guarded
  ≥1); deterministic GC-sorted-spread centroid init; max 50 iters, early stop; empty-cluster skip;
  min-bin-size filter; completeness `min(len/expected·100,100)`; contamination `min(σ/0.5·100,100)`,
  population σ, 0 for <2 contigs. Faithfully realises the validated description.
- **Independent cross-check (reflection probe against the built assembly + public API, exact numbers):**
  | Quantity | Input | Hand value | C# value | Match |
  |---|---|---|---|---|
  | `TnfPearsonDistance` idealised disjoint | va={GCGC:.5,CGCG:.5}, vb={ATAT:.5,TATA:.5} → r=−1 → (1−(−1))/2 | **1.0** | 1.0 (exact) | ✓ |
  | `TnfPearsonDistance` real GCGC×16 vs ATAT×16 | 29 tetramers, 15/29 & 14/29 split → r≈−1 | ≈1.0 | 0.99881376 | ✓ (finite-length boundary) |
  | `TnfPearsonDistance` self | identical vector → r=1 | **0.0** | 1.11e−16 | ✓ |
  | `CalculateGcContent` GCGC… / ATAT… | (G+C)/total | 1.0 / 0.0 | 1.0 / 0.0 | ✓ |
  | **Composition separation** | 10×high-GC(cov 50) + 10×low-GC(cov 5) | 2 unmixed bins | bin.1 GC=1 cov=50; bin.2 GC=0 cov=5; mixed=false | ✓ |
  | **Differential-coverage separation** | 20× identical GCTA repeat, cov 5 vs 500 | 2 unmixed bins (only coverage differs) | bin.1 cov=5; bin.2 cov=500; mixed=false | ✓ |
- The differential-coverage case is the decisive cross-check: with identical composition (|ΔGC|=0,
  TNF dist=0) the composite distance reduces to the coverage term, and the binner still separates the
  two abundance groups — confirming the coverage feature does real work (MetaBAT/CONCOCT signal).
- **Test-quality gate / defect found & fixed:** the existing 18-test fixture covered composition
  (GC) separation but had **no test isolating the differential-coverage signal** — a gap against the
  defining feature of the binning literature. Added **M13
  `BinContigs_DifferentialCoverage_SeparatesSameCompositionGenomes`**: identical GCTA composition,
  coverage 5× vs 500× → asserts ≥2 bins, no mixed bin, and exact per-bin mean coverage 5.0 / 500.0
  (values traced to the hand cross-check, not a code echo). Fixture now **19/19**.
- **Stage B verdict: PASS-WITH-NOTES** — code realises the validated (simplified) description; every
  hand-computed cross-check reproduced exactly; edge cases (empty, single/below-min, short <4 bp,
  zero coverage, identical contigs, extreme GC) handled; differential-coverage path now tested.

### Re-validation result
- **State: CLEAN.** No code change. One test-coverage defect (missing differential-coverage test)
  fully fixed. Per-lineage CheckM DB remains a documented out-of-scope data boundary (META-CHECKM-001).
- **Tests:** GenomeBinning fixture 19/19; full unfiltered `dotnet test Seqeron.sln -c Debug` =
  Failed 0 (`Seqeron.Genomics.Tests` 18805/0), 0 warnings.

## Runtime enforcement (LimitationPolicy)

This unit's guarded branch — domain-level CheckM completeness/contamination (no lineage-specific refinement) — has **minimum access mode `Moderate`** (`Seqeron.Genomics.Core.LimitationCatalog`). Under the default `LimitationPolicy.DefaultMode = Moderate` it is **allowed** (this guarded branch throws only under `Strict`); see [LIMITATIONS.md](../LIMITATIONS.md) › Runtime enforcement. Additive policy layer; the validated contract and `✅ CLEAN` verdict are unchanged.
