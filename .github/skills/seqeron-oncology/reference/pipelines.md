# seqeron-oncology — fuller C#-API pipelines & parameter guidance

C#-API only (no MCP oncology server). All calls are `Seqeron.Genomics.Oncology.*`; Method IDs +
partial files are in [`api-conventions.md`](api-conventions.md). Rigor (provenance, cross-check,
units, clinical caveat) is owned by [`bio-rigor`](../../bio-rigor/SKILL.md); C# construction / policy
mechanics by [`seqeron-dev`](../../seqeron-dev/SKILL.md). Paths from repo root `../../../../`.

**Clinical caveat (repeat on every reportable result):** these outputs are alpha / research-grade
and **not for clinical or diagnostic use** — independently validate before relying on any number.

## 1. Purity + ploidy (ASCAT) and purity QC (CNAqc)
`FitPurityPloidy` (`OncologyAnalyzer.CopyNumberPloidy.cs`) is a port of ASCAT `runASCAT` (grid search +
goodness of fit); defaults match ASCAT: `purityMin=0.1, purityMax=1.05, purityStep=0.01, ploidyMin=1.5,
ploidyMax=5.5, ploidyStep=0.05, gamma=AscatSequencingGamma` (1.0, sequencing data). Input is
`IReadOnlyList<AlleleSpecificSegmentSummary>` from `SegmentAlleleSpecificAspcf(loci, penalty = 70)`
(port of `ascat.aspcf`). `SegmentAlleleSpecific` is a compatibility alias that now runs the same ASPCF
(its threshold arguments are ignored).

- **Germline-aware path (preferred when germline genotypes are known):**
  `SegmentAlleleSpecificAspcf(loci, germlineHeterozygous)` → `AspcfSegmentation` →
  `FitPurityPloidyFromAspcf(segmentation)` (homozygous segments included; ploidy over all probes, as
  ASCAT). Multi-sample: `SegmentAlleleSpecificAsMultiPcf(samples, germlineHeterozygous?)`.
- **Male samples:** pass `AscatSexModel.MaleWithXNonPar(ReferenceGenome)` to the `FitPurityPloidy*`
  overloads (haploid X/Y, runASCAT `gender = "XY"`) and `AscatMaleXGenotyping` to `SegmentAlleleSpecificAspcf`.
  Default is `AscatSexModel.Female` (X/Y excluded from the fit).
- **Cross-check** the fitted purity against `EstimatePurity` / `EstimatePurityFromVAF`
  (`SomaticCalling.cs`). Divergence flags model strain (subclonal CN, WGD).
- **Purity QC (CNAqc):** `AnalyzePurityPeaks(mutations, purity, options?)` (simple karyotypes, KDE +
  peakPick, `purity_error` 0.05 → PASS/FAIL), `AnalyzeComplexKaryotypePeaks`, `AnalyzeSubclonalPurityPeaks`;
  BMix peaks via `PurityPeakOptions.FitMixturePeaks` / `FitBinomialMixture(successes, trials, seed)`;
  `PurityPeakOptions.BootstrapCount` = CNAqc `n_bootstrap`. CNAqc QCs a supplied purity — it is not a
  point estimator.
- **Ploidy / WGD:** `EstimatePloidy(segments)` (bp-weighted) or `EstimatePloidy(segments, probeCounts)`
  (ASCAT probe weighting); `DetectWholeGenomeDoubling(segments)` (facets-suite: denominator = sample span),
  `DetectWholeGenomeDoubling(segments, ReferenceGenome)` (reference autosome length),
  `ComputeAscatGenomeMetrics(segments)` (ASCAT `ascat.metrics`: mode of nMajor WGD, GI, LOH).
- **Sub-clonal copy number (Battenberg):** `FitSubclonalCopyNumber(summaries, purity, ploidy)`
  (`determine_copynumber`, maxdist only); `FitSubclonalCopyNumberWithSnpTest` (per-SNP t-test,
  siglevel 0.05); `FitSubclonalCopyNumberWithBootstrap(…, seed)` (solutions A–F, SDfrac, CIs). Phased
  path: `SegmentPhasedBaf` → `BuildBattenbergSegments`; full `callSubclones` driver:
  `CallBattenbergSubclones(segmentedSnps, logR, purity, ploidy, seed)` (uses `MergeBattenbergSegments`
  and `MaskHighCopyNumberSegments`). Haplotype imputation (IMPUTE2/Beagle5) is **not** built in — supply
  caller-phased BAFs.
- **Envelope:** ONCO-ASCAT-001 / ONCO-PURITY-001 have no runtime guard; caveats are in
  `Allele_Specific_Copy_Number_Derivation.md` / `Tumor_Purity_Estimation.md` §5–6.

## 2. Copy-number calling
Convert coverage log2-ratios: `Log2RatioToCopyNumber(log2Ratio, ploidy=2.0)` (CNVkit; add `purity` for
CNVkit's purity-rescaled `_log2_ratio_to_absolute`), then `ClassifyCopyNumbers` against
`DefaultCopyNumberThresholds` (−1.1, −0.25, 0.2, 0.7). All in `OncologyAnalyzer.CopyNumberPloidy.cs`.
- **Focal events:** per-segment `DetectFocalAmplifications` (GISTIC2 focal filter; marker-unit arm
  fraction when `MarkerCount`/`ArmMarkerCount` are set); GISTIC2 ziggurat deconstruction
  `DeconstructZiggurat(chromosomes, samples, options?)` → `DetectFocalAmplificationEvents(deconstruction)`.
  Arm boundaries: `GetChromosomeArmLengths(ReferenceGenome)` (UCSC cytoBand).
- **Deletions:** `DetectHomozygousDeletions` (CNVkit absolute CN 0; optional `purity`).
- **Genes:** `IdentifyAmplifiedOncogenes` / `IdentifyDeletedTumorSuppressors` — arm-label overloads or
  `GeneLocus` locus-overlap overloads (GISTIC2 `genes_at`), with default panels `DefaultOncogeneLoci` / `DefaultTumorSuppressorLoci` (GRCh38).
- For chromothripsis-like signals combine `CountCopyNumberStateOscillations` + `TestBreakpointClustering`
  (`ClinicalMisc.cs`).

## 3. Mutational signatures (SBS96) — `OncologyAnalyzer.Signatures.cs`
1. **Catalog.** For each SNV with its flanking bases, `ClassifySbsContext(fivePrime, ref, alt,
   threePrime)` → one of 96 channels; `EnumerateSbs96Channels` gives the canonical
   ordering to bin counts into the 96-vector. (Read `SBS96_Trinucleotide_Context_Catalog.md` for the
   pyrimidine-strand convention — do not re-derive it.)
2. **Fit vs extract.** Known catalog → `FitSignatures(catalog, referenceSignatures)` (NNLS
   refit). De-novo discovery → `ExtractSignatures(...)` (NMF). Label de-novo signatures
   with `MatchToReferenceSignatures`.
3. **Confidence.** `BootstrapExposures` → per-signature CIs. **Report exposures with their
   intervals**, never point estimates alone. Optionally `ClassifyMutationalProcess`.

## 4. Fusions — `OncologyAnalyzer.Fusions.cs`
`DetectFusions(candidates, thresholds?)` applies `FusionDetectionThresholds` (support/
read-through filters). Per call: `IsInFrame(fivePrimeCodingBases, threePrimeStartPhase)`,
`MatchKnownFusions` for curated hits, `AnalyzeBreakpoint` + `PredictFusionProtein`
for the junction protein. `ComputeTotalSupport` gives the evidence count used
by the threshold gate — a good cross-check.

## 5. TMB / MSI — `OncologyAnalyzer.TmbMsi.cs`
- **TMB.** `CalculateTMB(mutationCount, targetRegionMb)` or from `SomaticCall`s;
  result is **mutations/Mb** — report the unit and the panel size used. `ClassifyTMB` →
  `TmbStatus`. Panel-vs-WES calibration is out of scope; state the panel.
- **MSI.** `CalculateMSIScore(unstableLoci, totalLoci)` → fraction; `ClassifyMSIStatus`
  or the marker-count `ClassifyBethesdaPanel(unstableMarkers, totalMarkers)`.
  `DetectMSI(locusUnstableFlags)` wraps flags → `MsiResult`.

## 6. Neoantigens + MHC binding (guarded scorer) — `OncologyAnalyzer.Neoantigen.cs`
1. `GenerateNeoantigenPeptides(wildTypeProtein, mutantResidue, mutationPosition, minLength=8,
   maxLength=11)` → mutant peptide windows spanning the substitution (MHC-I lengths;
   validate lengths with `IsValidPeptideLength`).
2. **Binding — prefer the ANN.** `MhcflurryAffinityPredictor.LoadWeightPack(stream)`
   (`MhcflurryAffinityPredictor.cs`) then `PredictIc50(networks, peptide, allele)` or
   `PredictAndClassify`. The **caller supplies the weight pack** (not bundled).
3. `ClassifyMhcBinding(peptideLength, ic50Nm, mhcClass)` turns any IC50 into a
   `BindingStrength` — this classifier is **ungated** (it only classifies a supplied value).
4. **Matrix scorers are GUARDED (ONCO-MHC-001):** `PredictIc50Smm` and
   `PredictBindingHalfLifeBimas` each call `LimitationPolicy.Enforce("ONCO-MHC-001")` and
   require **caller-supplied `PmhcScoringMatrix` coefficients** (no redistributable matrix ships).
   Under `Strict` they throw `SeqeronLimitationException`. **STOP rule:** if no licensed matrix is
   available, do not fabricate one — report ONCO-MHC-001 and use the vendor NetMHCpan-4.1 model or
   the ANN path instead.

## 7. Immune infiltration / deconvolution (GUARDED ONCO-IMMUNE-001)
`ImmuneAnalyzer.EstimateInfiltration` (`ImmuneAnalyzer.cs`) is the signature-score entry.
`EstimateTumorPurity(estimateScore)` and `DeconvoluteImmuneCells[NuSvr]` each
call `LimitationPolicy.Enforce("ONCO-IMMUNE-001")` — `MinimumMode = Moderate`,
so they run by default but throw under `Strict`. The bundled matrix is **ABIS, not exact
CIBERSORT-LM22**, and the ESTIMATE purity transform is Affymetrix-calibrated (Yoshihara 2013). Report
both caveats; supply your own matrix via `LoadSignatureMatrix` for other platforms.

## 8. Clonal structure → subclones → phylogeny
1. `ClassifyClonality(clonalityVariants, purity)` (`Clonality.cs`) → clonal vs subclonal per variant.
2. `EstimateCcf(vaf, purity, tumorCopyNumber, multiplicity)`; multiplicity from
   `DeriveMultiplicity(vaf, purity, totalCN, majorCN)` (`CopyNumberPloidy.cs`, R half-to-even rounding).
3. `ClusterCcfValues(ccfValues, clusterCount)` (Ckmeans.1d.dp exact 1-D k-means) or
   `ClusterCcfValues(ccfValues, minClusters, maxClusters)` (k chosen by BIC, `select_levels`);
   `InferSubclones` / `AnalyzeHeterogeneity` (`PhylogenyHeterogeneity.cs`) for a summary.
   ITH: `CalculateITH(vafs)` (maftools MATH) or `CalculateITH(vafs, vafCutOff, minMutations = 5)`
   (maftools `math.score` filters, VAF cut-off 0.075); `CalculateCloneShannonDiversity(cloneFractions)`.
4. `ReconstructPhylogeny(clusters, tolerance = 0.1)` (LICHeE; default ε 0.1) or, with per-cluster
   member CCFs, `ReconstructPhylogenyFromClusterSummaries(summaries)` built with
   `CcfClusterSummary.FromMembers(id, memberCcfs)` (LICHeE 1.96·sd/√n margins + `fixNetwork`).
   `IdentifyTrunkMutations` / `IdentifyBranchMutations` split truncal vs branch events.

## 9. HRD / LOH / HLA-LOH / ctDNA
- HRD (`HrdLoh.cs`): `DetectHRD(segments, tai, lst)`, `CalculateHrdLohScore(segments)`.
- LOH (`HrdLoh.cs`): `DetectLOH(segments)`, `CalculateLOHFraction(segments, chromosome)`.
- HLA (`ClinicalMisc.cs`): `ParseHlaAllele`/`TryParseHlaAllele`, `DetectHlaLoh(alleleCopyNumber)`
  — HLA typing is clinical-flavoured; caveat applies strongly.
- ctDNA / MRD (`CtdnaMrdChip.cs`): `IntegratedMutantAlleleFractionV2(loci)` — INVAR-style IMAF.

## Provenance template (fill per pipeline)
```
Provenance
1) Seqeron.Genomics.Oncology.<Class>.<Method>(args) → <output>
2) <independent cross-check method>(args) → <corroborating output>
Envelope: <ONCO-* id(s)> — <guarded@Moderate / documented-only>; STOP under Strict for guarded units.
Coordinates/units: <e.g. TMB = mut/Mb; VAF ∈ [0,1]; CCF ∈ [0,1]>.
Caveat: alpha, research-grade — NOT for clinical/diagnostic use; validate independently.
```
