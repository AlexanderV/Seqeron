# seqeron-oncology — fuller C#-API pipelines & parameter guidance

C#-API only (no MCP oncology server). All calls are `Seqeron.Genomics.Oncology.*`; Method IDs +
`file:line` are in [`api-conventions.md`](api-conventions.md). Rigor (provenance, cross-check,
units, clinical caveat) is owned by [`bio-rigor`](../../bio-rigor/SKILL.md); C# construction / policy
mechanics by [`seqeron-dev`](../../seqeron-dev/SKILL.md). Paths from repo root `../../../../`.

**Clinical caveat (repeat on every reportable result):** these outputs are alpha / research-grade
and **not for clinical or diagnostic use** — independently validate before relying on any number.

## 1. Purity + ploidy (ASCAT-style grid)
`FitPurityPloidy` (`OncologyAnalyzer.cs:7411`) does a joint grid search over the region fits; defaults:
`purityMin=0.05, purityMax=1.0, purityStep=0.01, ploidyMin=1.5, ploidyMax=5.0, ploidyStep=0.05,
gamma=AscatSequencingGamma`. Input is `IReadOnlyList<AlleleSpecificSegmentSummary>` — build it from
`SegmentAlleleSpecificAspcf` (`:7641`, ASPCF-style joint segmentation of logR + BAF) or the simpler
`SegmentAlleleSpecific` (`:7282`).

- Cross-check the fitted purity against an independent VAF-based estimate `EstimatePurity` (`:515`)
  or `EstimatePurityFromVAF` (`:451`). Divergence flags model strain (subclonal CN, WGD).
- WGD signal: `DetectWholeGenomeDoublingFromSuppliedLength` (`:6659`).
- **Envelope:** ONCO-ASCAT-001 / ONCO-PURITY-001 are documented-only (no runtime guard) — the grid
  fit is exact for the stated ASCAT-style NN model; report the model caveat, nothing throws.

## 2. Copy-number calling
Convert coverage log2-ratios: `Log2RatioToCopyNumber(log2Ratio, ploidy=2.0)` (`:5900`), then
`ClassifyCopyNumbers` (`:5995`) against `DefaultCopyNumberThresholds` (`:5835`). Focal events:
`DetectFocalAmplifications` (`:6157`), `DetectHomozygousDeletions` (`:6297`). For chromothripsis-like
signals combine `CountCopyNumberStateOscillations` (`:11823`) + `TestBreakpointClustering` (`:11850`).

## 3. Mutational signatures (SBS96)
1. **Catalog.** For each SNV with its flanking bases, `ClassifySbsContext(fivePrime, ref, alt,
   threePrime)` (`:2659`) → one of 96 channels; `EnumerateSbs96Channels` (`:2696`) gives the canonical
   ordering to bin counts into the 96-vector. (Read `SBS96_Trinucleotide_Context_Catalog.md` for the
   pyrimidine-strand convention — do not re-derive it.)
2. **Fit vs extract.** Known catalog → `FitSignatures(catalog, referenceSignatures)` (`:2933`, NNLS
   refit). De-novo discovery → `ExtractSignatures(...)` (`:3431`/`:3485`, NMF). Label de-novo signatures
   with `MatchToReferenceSignatures` (`:4609`).
3. **Confidence.** `BootstrapExposures` (`:4774`) → per-signature CIs. **Report exposures with their
   intervals**, never point estimates alone. Optionally `ClassifyMutationalProcess` (`:5189`).

## 4. Fusions
`DetectFusions(candidates, thresholds?)` (`:5421`) applies `FusionDetectionThresholds` (support/
read-through filters). Per call: `IsInFrame(fivePrimeCodingBases, threePrimeStartPhase)` (`:5382`),
`MatchKnownFusions` (`:5558`) for curated hits, `AnalyzeBreakpoint` (`:5703`) + `PredictFusionProtein`
(`:5750`) for the junction protein. `ComputeTotalSupport` (`:5365`) gives the evidence count used
by the threshold gate — a good cross-check.

## 5. TMB / MSI
- **TMB.** `CalculateTMB(mutationCount, targetRegionMb)` (`:1508`) or from `SomaticCall`s (`:1540`);
  result is **mutations/Mb** — report the unit and the panel size used. `ClassifyTMB` (`:1557`) →
  `TmbStatus`. Panel-vs-WES calibration is out of scope; state the panel.
- **MSI.** `CalculateMSIScore(unstableLoci, totalLoci)` (`:1634`) → fraction; `ClassifyMSIStatus`
  (`:1664`) or the marker-count `ClassifyBethesdaPanel(unstableMarkers, totalMarkers)` (`:1688`).
  `DetectMSI(locusUnstableFlags)` (`:1725`) wraps flags → `MsiResult`.

## 6. Neoantigens + MHC binding (guarded scorer)
1. `GenerateNeoantigenPeptides(wildTypeProtein, mutantResidue, mutationPosition, minLength=8,
   maxLength=11)` (`:8026`) → mutant peptide windows spanning the substitution (MHC-I lengths;
   validate lengths with `IsValidPeptideLength` `:8273`).
2. **Binding — prefer the ANN.** `MhcflurryAffinityPredictor.LoadWeightPack(stream)`
   (`MhcflurryAffinityPredictor.cs:516`) then `PredictIc50(networks, peptide, allele)` (`:619`) or
   `PredictAndClassify` (`:673`). The **caller supplies the weight pack** (not bundled).
3. `ClassifyMhcBinding(peptideLength, ic50Nm, mhcClass)` (`:8292`) turns any IC50 into a
   `BindingStrength` — this classifier is **ungated** (it only classifies a supplied value).
4. **Matrix scorers are GUARDED (ONCO-MHC-001):** `PredictIc50Smm` (`:8516`) and
   `PredictBindingHalfLifeBimas` (`:8483`) each call `LimitationPolicy.Enforce("ONCO-MHC-001")` and
   require **caller-supplied `PmhcScoringMatrix` coefficients** (no redistributable matrix ships).
   Under `Strict` they throw `SeqeronLimitationException`. **STOP rule:** if no licensed matrix is
   available, do not fabricate one — report ONCO-MHC-001 and use the vendor NetMHCpan-4.1 model or
   the ANN path instead.

## 7. Immune infiltration / deconvolution (GUARDED ONCO-IMMUNE-001)
`ImmuneAnalyzer.EstimateInfiltration` (`ImmuneAnalyzer.cs:402`) is the signature-score entry.
`EstimateTumorPurity(estimateScore)` (`:473`) and `DeconvoluteImmuneCells[NuSvr]` (`:516`/`:658`) each
call `LimitationPolicy.Enforce("ONCO-IMMUNE-001")` (`:477`/`:525`/`:680`) — `MinimumMode = Moderate`,
so they run by default but throw under `Strict`. The bundled matrix is **ABIS, not exact
CIBERSORT-LM22**, and the ESTIMATE purity transform is Affymetrix-calibrated (Yoshihara 2013). Report
both caveats; supply your own matrix via `LoadSignatureMatrix` (`:813`) for other platforms.

## 8. Clonal structure → subclones → phylogeny
1. `ClassifyClonality(clonalityVariants, purity)` (`:6789`) → clonal vs subclonal per variant.
2. `EstimateCcf(vaf, purity, tumorCopyNumber, multiplicity)` (`:7014`); multiplicity from
   `DeriveMultiplicity(vaf, purity, totalCN, majorCN)` (`:7579`).
3. `ClusterCcfValues(ccfValues, clusterCount)` (`:7061`) → clusters; `InferSubclones` (`:11009`);
   `AnalyzeHeterogeneity` (`:11042`) for a summary.
4. `ReconstructPhylogeny(...)` (`:10617`) → subclone tree; `IdentifyTrunkMutations` /
   `IdentifyBranchMutations` (`:10727`/`:10760`) split truncal vs branch events.

## 9. HRD / LOH / HLA-LOH / ctDNA
- HRD: `DetectHRD(segments, tai, lst)` (`:1896`), `CalculateHrdLohScore(segments)` (`:2017`).
- LOH: `DetectLOH(segments)` (`:1979`), `CalculateLOHFraction(segments, chromosome)` (`:2033`).
- HLA: `ParseHlaAllele`/`TryParseHlaAllele` (`:11281`/`:11350`), `DetectHlaLoh(alleleCopyNumber)`
  (`:11384`) — HLA typing is clinical-flavoured; caveat applies strongly.
- ctDNA / MRD: `IntegratedMutantAlleleFractionV2(loci)` (`:9087`) — INVAR-style IMAF.

## Provenance template (fill per pipeline)
```
Provenance
1) Seqeron.Genomics.Oncology.<Class>.<Method>(args) → <output>
2) <independent cross-check method>(args) → <corroborating output>
Envelope: <ONCO-* id(s)> — <guarded@Moderate / documented-only>; STOP under Strict for guarded units.
Coordinates/units: <e.g. TMB = mut/Mb; VAF ∈ [0,1]; CCF ∈ [0,1]>.
Caveat: alpha, research-grade — NOT for clinical/diagnostic use; validate independently.
```
