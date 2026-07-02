# seqeron-protein-features — pipelines & parameters

Fuller recipes for the protein-feature family. All tools are on the **analysis** server; input is a
**protein (amino-acid) sequence**. Rigor (tool-only computation, provenance, cross-check, the alpha
caveat) is delegated to [`bio-rigor`](../../bio-rigor/SKILL.md). Schemas: `docs/mcp/tools/analysis/<tool>.md`.
Algorithm contracts: [`ProteinPred/`](../../../../docs/algorithms/ProteinPred/) ·
[`ProteinMotif/`](../../../../docs/algorithms/ProteinMotif/). Envelope: [`LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md).

Coordinates: **0-based residue indices**; regions `[start,end]` inclusive. Report the base explicitly.

## 1. Disorder + IDRs + MoRFs

- `predict_disorder`(**sequence**, `windowSize`?=21, `disorderThreshold`?=0.542, `minRegionLength`?=5)
  → `residuePredictions[{position,residue,disorderScore,isDisordered}]`,
  `disorderedRegions[{start,end,meanScore,confidence,regionType}]`, `overallDisorderContent`, `meanDisorderScore`.
  A residue is disordered when its normalized TOP-IDP score exceeds `disorderThreshold`; a region needs
  ≥ `minRegionLength` residues. `regionType` ∈ {Proline-rich, Acidic, Basic, Ser/Thr-rich, Long IDR, Standard IDR} (heuristic labels).
- `predict_morfs`(**sequence**, `minLength`?=10, `maxLength`?=25) → MoRF intervals (short binding-prone segments inside disorder).
- Single residue: `disorder_propensity`(**aminoAcid**) → normalized TOP-IDP propensity; `is_disorder_promoting`(**aminoAcid**) → bool.
- **C#:** `DisorderPredictor.PredictDisorder` · `.PredictMoRFs` · `.GetDisorderPropensity` · `.IsDisorderPromoting`.
- **Envelope (DISORDER-REGION-001):** the per-residue / per-region `confidence` is uncalibrated and guarded.
  Report the **boundaries** + `overallDisorderContent` (both follow the validated TOP-IDP threshold). See §6.
- **Cross-check:** an IDR should overlap a hydrophobicity **trough** (§5) and usually a protein SEG region (§4).

## 2. Signal peptide + transmembrane topology

- `predict_signal_peptide`(**proteinSequence**, `prokaryote`?=false) → signal-peptide presence + cleavage-site position.
  Set `prokaryote=true` for bacterial/archaeal sequences (different scoring matrix).
- `predict_transmembrane_helices`(**proteinSequence**, `windowSize`?=19, `threshold`?=1.6) → TM segments + in/out topology.
- **C#:** `ProteinMotifFinder.PredictSignalPeptide` · `.PredictTransmembraneHelices`.
- **Cross-check:** each TM segment should sit on a hydrophobicity **peak** (§5). A called signal peptide occupies
  the extreme N-terminus — do not double-count it as the first TM helix.

## 3. Coiled-coils + domains + PROSITE / common motifs

- `predict_coiled_coils`(**proteinSequence**, `windowSize`?=28, `threshold`?=0.5) → coiled-coil intervals + score.
- `find_protein_domains`(**proteinSequence**) → domain calls.
- `find_protein_motifs`(**proteinSequence**) → common protein sequence motifs.
- Targeted PROSITE: `prosite_to_regex`(**prositePattern**) → inspect the compiled regex, then
  `find_motif_by_prosite`(**proteinSequence**, **prositePattern**, `motifName`?=Custom) → hits.
- **C#:** `ProteinMotifFinder.PredictCoiledCoils` · `.FindDomains` · `.FindCommonMotifs` · `.ConvertPrositeToRegex` · `.FindMotifByProsite`.
- **Scope reminder:** these are **protein** motifs. For **DNA** motifs (exact/degenerate/PWM, `MotifFinder.*`)
  use [`bio-annotation`](../../bio-annotation/SKILL.md); do not mix the two engines.

## 4. Protein low-complexity (SEG)

Two equivalent SEG entry points (protein, entropy-based):

- `predict_low_complexity_seg`(**sequence**, `triggerWindow`?=12, `triggerThreshold`?=2.2, `extensionThreshold`?=2.5, `minLength`?=1) — `DisorderPredictor.PredictLowComplexityRegions`.
- `find_protein_low_complexity_regions`(**proteinSequence**, `windowSize`?=12, `triggerComplexity`?=2.2, `extensionComplexity`?=2.5) — `ProteinMotifFinder.FindLowComplexityRegions`.
- Use as a **filter**: down-weight domain / PROSITE / motif hits that fall inside SEG regions (compositionally
  biased, prone to spurious matches). This is the protein analogue of DUST/SEG masking — do **not** reach for
  the nucleotide `SequenceComplexity.*` tools (those are DNA, in `bio-annotation`).

## 5. Hydrophobicity profile (supporting evidence)

- `hydrophobicity_profile`(**proteinSequence**, `windowSize`?=9) → per-window Kyte–Doolittle score — `SequenceStatistics.CalculateHydrophobicityProfile`.
- Not a standalone call in most workflows; it is the **cross-check backbone**: peaks corroborate TM helices (§2),
  troughs corroborate IDRs (§1). Report which windows support which feature.

## 6. Envelope — DISORDER-REGION-001 STOP rule

- **Unit:** DISORDER-REGION-001 — `predict_disorder` / `DisorderPredictor.PredictDisorder`, the per-residue /
  per-region **confidence** output.
- **MinimumMode:** **Permissive** (blocked in Strict & Moderate; library default is Moderate).
- **Why:** no disorder predictor publishes a calibrated confidence standard. The region **boundaries** are
  validated (TOP-IDP threshold), only the confidence value is non-ideal.
- **STOP rule:** if a task needs a **validated per-region confidence** and the call would run below its
  MinimumMode (throws `SeqeronLimitationException`), **STOP** and report the limitation — do **not** raise the
  LimitationPolicy mode to force output (bio-rigor). **Named alternative:** report the validated boundaries +
  `overallDisorderContent` via `DisorderPredictor.PredictDisorderRegions` and treat `confidence` as advisory
  only. Setting the policy mode in code (Permissive bootstrap) is a `seqeron-dev` concern.
- All other tools in this family (signal peptide, TM, coiled-coil, domains, motifs, SEG, hydrophobicity,
  propensity) are **not** guarded.

## Provenance template

```
Provenance
1) predict_disorder(sequence, minRegionLength=5) → disorderedRegions [0-based; boundaries validated, confidence advisory], overallDisorderContent
2) predict_signal_peptide(proteinSequence, prokaryote=?) → signal peptide + cleavage site
3) predict_transmembrane_helices(proteinSequence, threshold=1.6) → TM segments + topology
4) predict_coiled_coils(proteinSequence); find_protein_domains(proteinSequence); find_motif_by_prosite(...) → structure/motif calls
5) find_protein_low_complexity_regions(proteinSequence) → SEG intervals (down-weight overlapping motif hits)
6) hydrophobicity_profile(proteinSequence, windowSize=9) → peaks↔TM(2), troughs↔IDRs(1)
Coordinates: 0-based residue indices; regions [start,end] inclusive.
Envelope: DISORDER-REGION-001 — validated boundaries + overallDisorderContent only; STOP if a validated per-region confidence is required.
Caveat: alpha software; not for clinical use — independently validate.
```
