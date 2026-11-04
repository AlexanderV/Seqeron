---
name: seqeron-oncology
description: >-
  Cancer-genomics analysis with the Seqeron.Genomics C# API — tumor purity /
  ploidy (ASCAT-like), allele-specific copy number, mutational signatures
  (SBS96 + NMF/NNLS fitting), gene fusions, tumor mutational burden (TMB),
  microsatellite instability (MSI), neoantigen candidate peptides + MHC-binding
  classification, clonal / subclonal structure and tumor phylogeny, HRD / LOH,
  ctDNA / MRD. C#-API ONLY — there is NO MCP oncology server. Triggers:
  "estimate tumor purity / ploidy", "fit purity and ploidy (ASCAT)",
  "allele-specific copy number", "extract mutational signatures", "SBS96
  catalog", "fit signature exposures", "call gene fusions", "is this fusion
  in-frame", "compute TMB", "TMB status", "MSI status / score", "predict
  neoantigens", "MHC-peptide binding", "classify clonal vs subclonal", "cancer
  cell fraction / CCF", "reconstruct tumor phylogeny", "HRD score", "LOH",
  "ctDNA fraction / tumor fraction", "MRD / minimal residual disease".
  (Germline SNP / indel calling → bio-annotation; germline read-evidence SV / CNV →
  seqeron-structural-variants; chromosome-arm-scale amplification / deletion / aneuploidy →
  bio-chromosome; oncology owns the tumor allele-specific / clonal layer.) All
  results are RESEARCH-GRADE, ALPHA — NOT for clinical or diagnostic use.
allowed-tools: Read, Bash, Grep, Glob
---

# seqeron-oncology — cancer genomics via the Seqeron.Genomics C# API

Routing + orchestration skill for the **Oncology** algorithm family (~37 documented contracts,
`docs/algorithms/Oncology/*`). **C#-API ONLY: there is no MCP oncology server**, so every pipeline
below is a real `Seqeron.Genomics` call chain (verified `class.Method` + `file:line`), not an MCP
tool name. This is the oncology-specific extension of [`seqeron-dev`](../seqeron-dev/SKILL.md) — it
delegates *all* C# mechanics (namespaces, `TryCreate`, the 3-tier `LimitationPolicy`,
`SeqeronLimitationException`, the Permissive test bootstrap) there and does not restate them.

- **CLINICAL CAVEAT IS LOAD-BEARING HERE.** Purity, HLA/MHC typing, TMB/MSI, driver/pathogenicity,
  fusion, and neoantigen outputs read as clinical. They are **alpha / research-grade and NOT for
  clinical or diagnostic use** — surface this on every reportable result (see
  [`bio-rigor`](../bio-rigor/SKILL.md) rule 6).
- **Rigor is delegated** to [`bio-rigor`](../bio-rigor/SKILL.md): tool-only computation, provenance,
  envelope, cross-check, units/coords. Do not restate — it applies by default.
- **Point, don't duplicate.** Contracts / invariants live in `docs/algorithms/Oncology/*.md`; the
  operating envelope in `docs/Validation/LIMITATIONS.md`. Link, never copy.

## Namespace & shape (verified in source)

Everything below lives in **one namespace `Seqeron.Genomics.Oncology`**, across three static classes:

| Class | Covers | Source (from repo root `../../../`) |
|---|---|---|
| `OncologyAnalyzer` | purity/ploidy/ASCAT/Battenberg, copy number, signatures, fusions, TMB/MSI, neoantigen windows, matrix pMHC, clonal/phylogeny, HRD/LOH, ctDNA | `static partial` class split over `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.*.cs` (e.g. `.CopyNumberPloidy.cs`, `.Clonality.cs`, `.PhylogenyHeterogeneity.cs`, `.SomaticCalling.cs`, `.Signatures.cs`) |
| `ImmuneAnalyzer` | immune-infiltration / deconvolution + ESTIMATE-purity (**guarded ONCO-IMMUNE-001**) | `src/.../Seqeron.Genomics.Oncology/ImmuneAnalyzer.cs` |
| `MhcflurryAffinityPredictor` | ANN IC50 prediction (caller supplies the weight pack) | `src/.../Seqeron.Genomics.Oncology/MhcflurryAffinityPredictor.cs` |

Locate any method by name (line numbers drift, so this skill cites method + partial file, not lines):
`grep -n "public static .* <Name>(" src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/*.cs`.
For the full curated Method-ID index (with partial file + algorithm doc) see
[`reference/api-conventions.md`](reference/api-conventions.md).

## Guarded units — STOP rule (in scope)

Both throw `SeqeronLimitationException` when the effective mode is **below `Moderate`** (i.e. under
`Strict`). Default mode is `Moderate`, so they run by default; a `Strict` process must **STOP and
report** (name the unit + workaround), never widen to force output. Read
[`bio-rigor` envelope](../bio-rigor/reference/envelope.md) + `docs/Validation/LIMITATIONS.md`.

- **ONCO-MHC-001** — matrix (SMM/BIMAS) pMHC scorers `PredictIc50Smm` / `PredictBindingHalfLifeBimas`
  each call `LimitationPolicy.Enforce("ONCO-MHC-001")` (`OncologyAnalyzer.Neoantigen.cs`).
  No redistributable matrix ships — the **caller supplies the coefficients**. Workaround: use a
  vendor model (NetMHCpan-4.1) or `ClassifyMhcBinding` on a caller-supplied IC50.
- **ONCO-IMMUNE-001** — `ImmuneAnalyzer.EstimateTumorPurity` / `DeconvoluteImmuneCells[NuSvr]`
  call `LimitationPolicy.Enforce("ONCO-IMMUNE-001")` (`ImmuneAnalyzer.cs`).
  Ships a bundled ABIS matrix, **not** exact CIBERSORT-LM22; ESTIMATE purity is Affymetrix-calibrated.
- **ONCO-ASCAT-001 / ONCO-PURITY-001** have **no runtime guard** — they are ports of ASCAT
  (`runASCAT`, `ascat.aspcf`, `ascat.asmultipcf`, `ascat.metrics`), Battenberg and CNAqc, locked against
  the R reference; report the model caveat (algorithm doc §5–6), nothing throws.

## Decision guide — task → C# entry point

| Task | Method ID (`OncologyAnalyzer.*` unless noted) · partial file |
|---|---|
| Fit tumor purity + ploidy jointly (ASCAT `runASCAT`) | `FitPurityPloidy` (segment summaries) / `FitPurityPloidyFromAspcf` (germline-aware ASPCF, all-probe ploidy) · `CopyNumberPloidy` |
| Male (haploid X/Y) ASCAT model | `AscatSexModel.MaleWithXNonPar(genome)` passed to the `FitPurityPloidy*` / `SegmentAlleleSpecificAspcf` overloads · `CopyNumberPloidy` |
| Quick purity from VAF cluster | `EstimatePurity` / `EstimatePurityFromVAF` · `SomaticCalling` |
| QC a supplied purity (CNAqc peaks) | `AnalyzePurityPeaks` / `AnalyzeComplexKaryotypePeaks` / `AnalyzeSubclonalPurityPeaks`; BMix `FitBinomialMixture` · `CopyNumberPloidy` |
| Ploidy / WGD from allele-specific segments | `EstimatePloidy` (bp; probe-weighted overload) / `DetectWholeGenomeDoubling` / `ComputeAscatGenomeMetrics` · `CopyNumberPloidy` |
| Segment allele-specific CN (ASPCF) | `SegmentAlleleSpecificAspcf` / multi-sample `SegmentAlleleSpecificAsMultiPcf` · `CopyNumberPloidy` |
| Subclonal copy number (Battenberg) | `FitSubclonalCopyNumber` → `…WithSnpTest` / `…WithBootstrap`; driver `CallBattenbergSubclones` · `CopyNumberPloidy` |
| log2-ratio → copy number; classify amp/del | `Log2RatioToCopyNumber` → `ClassifyCopyNumbers` (optional `purity` overloads) · `CopyNumberPloidy` |
| GISTIC2 focal / broad events | `DeconstructZiggurat` → `DetectFocalAmplificationEvents`; per-segment `DetectFocalAmplifications` · `CopyNumberPloidy` |
| Homozygous deletions; affected genes | `DetectHomozygousDeletions`; `IdentifyAmplifiedOncogenes` / `IdentifyDeletedTumorSuppressors` (`GeneLocus` overlap overloads) · `CopyNumberPloidy` |
| SBS96 channel of a mutation / all 96 channels | `ClassifySbsContext` / `EnumerateSbs96Channels` · `Signatures` |
| Extract signatures de novo (NMF) | `ExtractSignatures` · `Signatures` |
| Fit / refit known signatures (NNLS) | `FitSignatures` → match: `MatchToReferenceSignatures` · `Signatures` |
| Bootstrap exposure confidence intervals | `BootstrapExposures` · `Signatures` |
| Call fusions; is it in-frame? | `DetectFusions` / `IsInFrame` · `Fusions` |
| Known-fusion DB lookup; breakpoint analysis | `MatchKnownFusions` / `AnalyzeBreakpoint` · `Fusions` |
| TMB (count or calls) + status | `CalculateTMB` → `ClassifyTMB` · `TmbMsi` |
| MSI score + status (Bethesda) | `CalculateMSIScore` → `ClassifyMSIStatus` / `ClassifyBethesdaPanel` · `TmbMsi` |
| Neoantigen candidate peptide windows | `GenerateNeoantigenPeptides` · `Neoantigen` |
| Classify a supplied IC50 → binding strength | `ClassifyMhcBinding` · `Neoantigen` |
| ANN IC50 (caller weight pack) | `MhcflurryAffinityPredictor.PredictIc50` · `MhcflurryAffinityPredictor.cs` |
| Clonal vs subclonal; CCF | `ClassifyClonality` / `EstimateCcf` · `Clonality` |
| CCF clustering (Ckmeans.1d.dp) | `ClusterCcfValues(ccf, k)` or BIC-selected `ClusterCcfValues(ccf, minK, maxK)` · `Clonality` |
| Reconstruct tumor phylogeny (LICHeE) | `ReconstructPhylogeny` / `ReconstructPhylogenyFromClusterSummaries` · `PhylogenyHeterogeneity` |
| Heterogeneity (MATH, Shannon) | `CalculateITH` / `CalculateCloneShannonDiversity` / `AnalyzeHeterogeneity` · `PhylogenyHeterogeneity` |
| HRD score / LOH | `CalculateHrdLohScore` / `DetectLOH` · `HrdLoh` |

## Canonical C#-API pipelines

### (a) Tumor purity + ploidy from allele-specific segments (ASCAT)
1. `SegmentAlleleSpecificAspcf(loci, germlineHeterozygous)` → `AspcfSegmentation` (ASCAT `ascat.aspcf`, penalty 70).
2. `FitPurityPloidyFromAspcf(segmentation)` → `PurityPloidyFit` (ASCAT `runASCAT` grid; defaults
   `purityMin=0.1, purityMax=1.05, purityStep=0.01, ploidyMin=1.5, ploidyMax=5.5, ploidyStep=0.05, gamma=1.0`).
   Without germline genotypes: `SegmentAlleleSpecificAspcf(loci)` → summaries → `FitPurityPloidy(summaries)` (same defaults).
3. Cross-check purity via an independent VAF path: `EstimatePurity(variants)` — magnitudes should agree;
   `AnalyzePurityPeaks(mutations, purity)` gives CNAqc's PASS/FAIL purity QC.
```
Provenance
1) OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci, germlineHet) → segmentation
2) OncologyAnalyzer.FitPurityPloidyFromAspcf(segmentation) → {purity, ploidy, goodness of fit}
3) OncologyAnalyzer.EstimatePurity(vafVariants) → purity (independent cross-check); AnalyzePurityPeaks → QC
Envelope: ONCO-ASCAT-001 / ONCO-PURITY-001 — no runtime guard; ASCAT/CNAqc ports (see algorithm docs §5–6).
Caveat: alpha, research-grade — NOT for clinical use; purity is decision-relevant, validate independently.
```

### (b) Mutational signatures — SBS96 catalog → fit → confidence
1. Per mutation: `ClassifySbsContext(fivePrime, ref, alt, threePrime)` → one of 96 channels; build the 96-vector catalog (`EnumerateSbs96Channels` gives the canonical order).
2. `FitSignatures(catalog, referenceSignatures)` → NNLS exposures; or `ExtractSignatures(...)` for de-novo NMF.
3. `MatchToReferenceSignatures(...)` to label; `BootstrapExposures(...)` for CIs — report exposures **with** their intervals.

### (c) Gene fusions — call → frame → annotate
1. `DetectFusions(candidates, thresholds?)` → `FusionCall[]`.
2. Per call: `IsInFrame(fivePrimeCodingBases, threePrimeStartPhase)`; `MatchKnownFusions(...)` for DB hits; `AnalyzeBreakpoint(...)` for junction detail.

### (d) TMB + MSI status
1. `CalculateTMB(mutationCount, targetRegionMb)` → `ClassifyTMB(tmb)` → `TmbStatus`.
2. `CalculateMSIScore(unstableLoci, totalLoci)` → `ClassifyMSIStatus(score)` (or `ClassifyBethesdaPanel`).
- Both are clinical-flavoured — attach the caveat; report units (TMB = mutations/Mb).

### (e) Neoantigen candidates + MHC binding (guarded scorer)
1. `GenerateNeoantigenPeptides(wildTypeProtein, mutantResidue, mutationPosition, minLength=8, maxLength=11)` → mutant peptide windows.
2. Score binding: **preferred** ANN `MhcflurryAffinityPredictor.PredictIc50(networks, peptide, allele)` with a caller weight pack; then `ClassifyMhcBinding(peptideLength, ic50Nm, mhcClass)`.
3. Matrix scorers `PredictIc50Smm` / `PredictBindingHalfLifeBimas` are **guarded (ONCO-MHC-001)** and need caller-supplied coefficients — if unavailable, STOP and report.

## End-to-end grounded example

**Task.** From allele-specific segments + somatic VAFs, estimate purity/ploidy, then place mutations on the clonal timeline.
1. `SegmentAlleleSpecificAspcf(loci, germlineHeterozygous)` → segmentation.
2. `FitPurityPloidyFromAspcf(segmentation)` → `{purity, ploidy}`.
3. `ClassifyClonality(clonalityVariants, purity)` → clonal vs subclonal; per variant `EstimateCcf(vaf, purity, tumorCN, multiplicity)`.
4. `ClusterCcfValues(ccfValues, 1, 9)` (BIC-selected k) → build `CcfClusterSummary.FromMembers(...)` per cluster →
   `ReconstructPhylogenyFromClusterSummaries(summaries)` (LICHeE, default ε 0.1) for the subclone tree.
```
Provenance
1) SegmentAlleleSpecificAspcf(loci, germlineHet) → segmentation
2) FitPurityPloidyFromAspcf(segmentation) → purity, ploidy
3) ClassifyClonality(variants, purity); EstimateCcf(vaf, purity, cn, mult) per variant
4) ClusterCcfValues(ccf, minK, maxK) → CcfClusterSummary.FromMembers → ReconstructPhylogenyFromClusterSummaries → tree
Envelope: ONCO-ASCAT-001 / ONCO-PURITY-001 — no runtime guard (ASCAT / LICHeE ports; caveats in algorithm docs).
Caveat: alpha, research-grade — NOT for clinical/diagnostic use; validate independently.
```

## Reference

- **Curated Method-ID map (all tasks; NOT in domain-map.json → NO generated slice):**
  [`reference/api-conventions.md`](reference/api-conventions.md)
- **Fuller recipes + parameter guidance:** [`reference/pipelines.md`](reference/pipelines.md)
- **Contracts / invariants (link, don't copy):** [`docs/algorithms/Oncology/`](../../../docs/algorithms/Oncology/)
  — e.g. [`Tumor_Purity_Estimation.md`](../../../docs/algorithms/Oncology/Tumor_Purity_Estimation.md),
  [`Allele_Specific_Copy_Number_Derivation.md`](../../../docs/algorithms/Oncology/Allele_Specific_Copy_Number_Derivation.md),
  [`SBS96_Trinucleotide_Context_Catalog.md`](../../../docs/algorithms/Oncology/SBS96_Trinucleotide_Context_Catalog.md),
  [`Fusion_Gene_Detection.md`](../../../docs/algorithms/Oncology/Fusion_Gene_Detection.md),
  [`Tumor_Mutational_Burden.md`](../../../docs/algorithms/Oncology/Tumor_Mutational_Burden.md),
  [`MHC_Peptide_Binding_Classification.md`](../../../docs/algorithms/Oncology/MHC_Peptide_Binding_Classification.md)
- **Operating envelope / guarded units:** [`LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) (ONCO-MHC-001, ONCO-IMMUNE-001); ONCO-ASCAT-001 / ONCO-PURITY-001 caveats live in their algorithm docs
- **Cross-cutting:** [`seqeron-dev`](../seqeron-dev/SKILL.md) (C# API mechanics — this skill's parent) ·
  [`bio-rigor`](../bio-rigor/SKILL.md) (rigor + clinical caveat) ·
  [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool/algorithm lookup)
