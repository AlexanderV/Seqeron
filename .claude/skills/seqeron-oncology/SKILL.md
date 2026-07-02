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
| `OncologyAnalyzer` | purity/ploidy/ASCAT, copy number, signatures, fusions, TMB/MSI, neoantigen windows, matrix pMHC, clonal/phylogeny, HRD/LOH, ctDNA | `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:30` |
| `ImmuneAnalyzer` | immune-infiltration / deconvolution + ESTIMATE-purity (**guarded ONCO-IMMUNE-001**) | `src/.../Seqeron.Genomics.Oncology/ImmuneAnalyzer.cs:23` |
| `MhcflurryAffinityPredictor` | ANN IC50 prediction (caller supplies the weight pack) | `src/.../Seqeron.Genomics.Oncology/MhcflurryAffinityPredictor.cs:45` |

`OncologyAnalyzer` is a large static class; find any method with
`grep -n "public static" src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`.
For the full curated Method-ID index see [`reference/api-conventions.md`](reference/api-conventions.md).

## Guarded units — STOP rule (in scope)

Both throw `SeqeronLimitationException` when the effective mode is **below `Moderate`** (i.e. under
`Strict`). Default mode is `Moderate`, so they run by default; a `Strict` process must **STOP and
report** (name the unit + workaround), never widen to force output. Read
[`bio-rigor` envelope](../bio-rigor/reference/envelope.md) + `docs/Validation/LIMITATIONS.md`.

- **ONCO-MHC-001** — matrix (SMM/BIMAS) pMHC scorers `PredictIc50Smm` / `PredictBindingHalfLifeBimas`
  each call `LimitationPolicy.Enforce("ONCO-MHC-001")` (`OncologyAnalyzer.cs:8523`, `:8490`).
  No redistributable matrix ships — the **caller supplies the coefficients**. Workaround: use a
  vendor model (NetMHCpan-4.1) or `ClassifyMhcBinding` on a caller-supplied IC50.
- **ONCO-IMMUNE-001** — `ImmuneAnalyzer.EstimateTumorPurity` / `DeconvoluteImmuneCells[NuSvr]`
  call `LimitationPolicy.Enforce("ONCO-IMMUNE-001")` (`ImmuneAnalyzer.cs:477`, `:525`, `:680`).
  Ships a bundled ABIS matrix, **not** exact CIBERSORT-LM22; ESTIMATE purity is Affymetrix-calibrated.
- **ONCO-ASCAT-001 / ONCO-PURITY-001** are **documented-only** (no runtime guard) — the purity/ploidy
  grid fit is exact for its stated ASCAT-style model; report the model caveat, nothing throws.

## Decision guide — task → C# entry point

| Task | Method ID (`OncologyAnalyzer.*` unless noted) · file:line |
|---|---|
| Fit tumor purity + ploidy jointly (ASCAT grid) | `FitPurityPloidy` · `:7411` |
| Quick purity from VAF cluster | `EstimatePurity` · `:515` / `EstimatePurityFromVAF` · `:451` |
| Ploidy from allele-specific segments | `EstimatePloidy` · `:6554` |
| Segment allele-specific CN (ASPCF) | `SegmentAlleleSpecificAspcf` · `:7641` / `SegmentAlleleSpecific` · `:7282` |
| log2-ratio → copy number; classify amp/del | `Log2RatioToCopyNumber` · `:5900` → `ClassifyCopyNumbers` · `:5995` |
| SBS96 channel of a mutation / all 96 channels | `ClassifySbsContext` · `:2659` / `EnumerateSbs96Channels` · `:2696` |
| Extract signatures de novo (NMF) | `ExtractSignatures` · `:3431` |
| Fit / refit known signatures (NNLS) | `FitSignatures` · `:2933` → match: `MatchToReferenceSignatures` · `:4609` |
| Bootstrap exposure confidence intervals | `BootstrapExposures` · `:4774` |
| Call fusions; is it in-frame? | `DetectFusions` · `:5421` / `IsInFrame` · `:5382` |
| Known-fusion DB lookup; breakpoint analysis | `MatchKnownFusions` · `:5558` / `AnalyzeBreakpoint` · `:5703` |
| TMB (count or calls) + status | `CalculateTMB` · `:1508`/`:1540` → `ClassifyTMB` · `:1557` |
| MSI score + status (Bethesda) | `CalculateMSIScore` · `:1634` → `ClassifyMSIStatus` · `:1664` / `ClassifyBethesdaPanel` · `:1688` |
| Neoantigen candidate peptide windows | `GenerateNeoantigenPeptides` · `:8026` |
| Classify a supplied IC50 → binding strength | `ClassifyMhcBinding` · `:8292` |
| ANN IC50 (caller weight pack) | `MhcflurryAffinityPredictor.PredictIc50` · `MhcflurryAffinityPredictor.cs:619` |
| Clonal vs subclonal; CCF | `ClassifyClonality` · `:6789` / `EstimateCcf` · `:7014` |
| CCF clustering → subclones; heterogeneity | `ClusterCcfValues` · `:7061` / `AnalyzeHeterogeneity` · `:11042` |
| Reconstruct tumor phylogeny | `ReconstructPhylogeny` · `:10617` |
| HRD score / LOH | `CalculateHrdLohScore` · `:2017` / `DetectLOH` · `:1979` |

## Canonical C#-API pipelines

### (a) Tumor purity + ploidy from allele-specific segments (ASCAT-like)
1. `SegmentAlleleSpecificAspcf(...)` → `AlleleSpecificSegmentSummary[]` (`OncologyAnalyzer.cs:7641`).
2. `FitPurityPloidy(segments, purityMin=0.05, purityMax=1.0, purityStep=0.01, ploidyMin=1.5, ploidyMax=5.0, ploidyStep=0.05, gamma=…)` → `PurityPloidyFit` (`:7411`).
3. Cross-check purity via an independent VAF path: `EstimatePurity(variants)` (`:515`) — magnitudes should agree.
```
Provenance
1) OncologyAnalyzer.SegmentAlleleSpecificAspcf(logR,baf,positions,…) → segments
2) OncologyAnalyzer.FitPurityPloidy(segments, grid defaults) → {purity, ploidy, gof}
3) OncologyAnalyzer.EstimatePurity(vafVariants) → purity (independent cross-check)
Envelope: ONCO-ASCAT-001 / ONCO-PURITY-001 documented-only — exact for the ASCAT-style grid model.
Caveat: alpha, research-grade — NOT for clinical use; purity is decision-relevant, validate independently.
```

### (b) Mutational signatures — SBS96 catalog → fit → confidence
1. Per mutation: `ClassifySbsContext(fivePrime, ref, alt, threePrime)` → one of 96 channels (`:2659`); build the 96-vector catalog (`EnumerateSbs96Channels` `:2696` gives the canonical order).
2. `FitSignatures(catalog, referenceSignatures)` → NNLS exposures (`:2933`); or `ExtractSignatures(...)` for de-novo NMF (`:3431`).
3. `MatchToReferenceSignatures(...)` (`:4609`) to label; `BootstrapExposures(...)` (`:4774`) for CIs — report exposures **with** their intervals.

### (c) Gene fusions — call → frame → annotate
1. `DetectFusions(candidates, thresholds?)` → `FusionCall[]` (`:5421`).
2. Per call: `IsInFrame(fivePrimeCodingBases, threePrimeStartPhase)` (`:5382`); `MatchKnownFusions(...)` (`:5558`) for DB hits; `AnalyzeBreakpoint(...)` (`:5703`) for junction detail.

### (d) TMB + MSI status
1. `CalculateTMB(mutationCount, targetRegionMb)` (`:1508`) → `ClassifyTMB(tmb)` → `TmbStatus` (`:1557`).
2. `CalculateMSIScore(unstableLoci, totalLoci)` (`:1634`) → `ClassifyMSIStatus(score)` (`:1664`) (or `ClassifyBethesdaPanel` `:1688`).
- Both are clinical-flavoured — attach the caveat; report units (TMB = mutations/Mb).

### (e) Neoantigen candidates + MHC binding (guarded scorer)
1. `GenerateNeoantigenPeptides(wildTypeProtein, mutantResidue, mutationPosition, minLength=8, maxLength=11)` → mutant peptide windows (`:8026`).
2. Score binding: **preferred** ANN `MhcflurryAffinityPredictor.PredictIc50(networks, peptide, allele)` (`MhcflurryAffinityPredictor.cs:619`) with a caller weight pack; then `ClassifyMhcBinding(peptideLength, ic50Nm, mhcClass)` (`:8292`).
3. Matrix scorers `PredictIc50Smm` (`:8516`) / `PredictBindingHalfLifeBimas` (`:8483`) are **guarded (ONCO-MHC-001)** and need caller-supplied coefficients — if unavailable, STOP and report.

## End-to-end grounded example

**Task.** From allele-specific segments + somatic VAFs, estimate purity/ploidy, then place mutations on the clonal timeline.
1. `SegmentAlleleSpecificAspcf(...)` → segments (`:7641`).
2. `FitPurityPloidy(segments, defaults)` → `{purity, ploidy}` (`:7411`).
3. `ClassifyClonality(clonalityVariants, purity)` → clonal vs subclonal (`:6789`); per variant `EstimateCcf(vaf, purity, tumorCN, multiplicity)` (`:7014`).
4. `ClusterCcfValues(ccfValues, clusterCount)` (`:7061`) → `ReconstructPhylogeny(...)` (`:10617`) for the subclone tree.
```
Provenance
1) SegmentAlleleSpecificAspcf(logR,baf,pos) → segments
2) FitPurityPloidy(segments) → purity,ploidy
3) ClassifyClonality(variants,purity); EstimateCcf(vaf,purity,cn,mult) per variant
4) ClusterCcfValues(ccf,k) → ReconstructPhylogeny(...) → subclone tree
Envelope: ONCO-ASCAT-001/ONCO-PURITY-001 documented-only (ASCAT-style model, exact for that model).
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
- **Operating envelope / guarded units:** [`LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) (ONCO-MHC-001, ONCO-IMMUNE-001, ONCO-ASCAT-001, ONCO-PURITY-001)
- **Cross-cutting:** [`seqeron-dev`](../seqeron-dev/SKILL.md) (C# API mechanics — this skill's parent) ·
  [`bio-rigor`](../bio-rigor/SKILL.md) (rigor + clinical caveat) ·
  [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool/algorithm lookup)
