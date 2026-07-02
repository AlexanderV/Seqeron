---
name: seqeron-protein-features
description: >-
  Predict features of a PROTEIN sequence with Seqeron (MCP tools OR the C# API).
  Use for intrinsic disorder / IDRs (TOP-IDP) and MoRFs, signal-peptide
  prediction + cleavage site, transmembrane helices & membrane topology,
  coiled-coil regions, protein domains, PROSITE / common protein motifs,
  protein low-complexity (SEG) regions, and Kyte–Doolittle hydrophobicity
  profiles. Triggers: "predict disorder in this protein", "find the intrinsically
  disordered regions", "find the signal peptide / where is it cleaved",
  "transmembrane regions of this protein / membrane topology", "coiled-coil
  prediction", "what domains / PROSITE motifs are in this protein",
  "hydrophobicity plot", "is this residue disorder-promoting". Input is an
  AMINO-ACID sequence. Server: analysis.
allowed-tools: Read, Bash, Grep, Glob
---

# seqeron-protein-features — disorder/IDRs, signal peptide, TM, coiled-coil, domains/motifs, hydrophobicity

Routing + orchestration skill for the **protein-feature-prediction** family on the **Analysis** server
(backing classes `DisorderPredictor.*`, `ProteinMotifFinder.*`, plus `SequenceStatistics` for the
hydrophobicity profile). Input is always a **protein (amino-acid) sequence**. It picks the right tool and
gives a **dual-mode** recipe (MCP tool call **and** the equivalent `Seqeron.Genomics` C# `Method ID`).

- **Rigor is delegated.** Tool-only computation, envelope, provenance, cross-check, coordinate reporting,
  and the alpha / not-for-clinical caveat are owned by **[`bio-rigor`](../bio-rigor/SKILL.md)** — it applies
  by default; do not restate its rules.
- **Don't know a tool name?** Use **[`seqeron-discovery`](../seqeron-discovery/SKILL.md)**
  (`python3 scripts/skills/find-tool.py <kw> --server analysis`) — never guess.
- **Point, don't duplicate.** Full I/O schemas live in `docs/mcp/tools/analysis/<tool>.md`; algorithm
  invariants in `docs/algorithms/ProteinPred/` + `docs/algorithms/ProteinMotif/`. This skill links, it does not copy.

## Scope boundary — protein features vs DNA annotation

This skill OWNS the **protein-SEQUENCE** feature family. The overloaded **[`bio-annotation`](../bio-annotation/SKILL.md)**
name-drops it but its real subject is **DNA-shaped** annotation. Route by input:

- **PROTEIN input → here.** Disorder, MoRFs, signal peptide, TM helices, coiled-coils, protein domains,
  **protein** PROSITE / common motifs, protein SEG low-complexity, hydrophobicity.
- **DNA/RNA input → `bio-annotation`.** ORFs/genes/promoters, variants, **DNA** motifs (exact/degenerate/PWM,
  `MotifFinder.*`), **nucleotide** low-complexity/DUST (`SequenceComplexity.*`), repeats, k-mers, splicing,
  methylation, RNA structure. Do **not** confuse `ProteinMotifFinder.FindMotifByProsite` (protein, here) with
  the DNA motif tools there. If you only have DNA, translate to protein first via `bio-qc` (`translate`) before entering this skill.

## Decision guide — which tool for which question

| Question (protein input) | Tool ([MCP] / `Method ID`) |
|---|---|
| Per-residue disorder + IDRs + disorder content (TOP-IDP) | `predict_disorder` / `DisorderPredictor.PredictDisorder` ⚠ confidence guarded |
| Molecular-recognition features within disorder | `predict_morfs` / `DisorderPredictor.PredictMoRFs` |
| Single amino-acid disorder propensity / is-promoting | `disorder_propensity` · `is_disorder_promoting` / `DisorderPredictor.GetDisorderPropensity` · `.IsDisorderPromoting` |
| Protein low-complexity (SEG) — from the disorder engine | `predict_low_complexity_seg` / `DisorderPredictor.PredictLowComplexityRegions` |
| Protein low-complexity (SEG) — from the motif engine | `find_protein_low_complexity_regions` / `ProteinMotifFinder.FindLowComplexityRegions` |
| N-terminal signal peptide + cleavage site | `predict_signal_peptide` / `ProteinMotifFinder.PredictSignalPeptide` |
| Transmembrane helices + membrane topology | `predict_transmembrane_helices` / `ProteinMotifFinder.PredictTransmembraneHelices` |
| Coiled-coil regions (heptad / MTIDK-style) | `predict_coiled_coils` / `ProteinMotifFinder.PredictCoiledCoils` |
| Protein domains | `find_protein_domains` / `ProteinMotifFinder.FindDomains` |
| Common protein sequence motifs | `find_protein_motifs` / `ProteinMotifFinder.FindCommonMotifs` |
| Scan a PROSITE pattern in a protein | `find_motif_by_prosite` / `ProteinMotifFinder.FindMotifByProsite` |
| Compile a PROSITE pattern → regex | `prosite_to_regex` / `ProteinMotifFinder.ConvertPrositeToRegex` |
| Kyte–Doolittle hydrophobicity profile | `hydrophobicity_profile` / `SequenceStatistics.CalculateHydrophobicityProfile` |

**⚠ Envelope note.** `predict_disorder`'s per-residue / per-region **confidence** is **DISORDER-REGION-001**
(uncalibrated). See the STOP rule below. The region **boundaries** themselves follow the validated TOP-IDP
threshold and are fine to report.

## Canonical dual-mode pipelines

Coordinates: protein positions are **0-based residue indices**; regions are `[start,end]` inclusive per the
tool doc. Report the base explicitly.

### (a) Disorder + IDRs → MoRFs (intrinsic-disorder characterization)
1. **[MCP]** `predict_disorder`(sequence, windowSize?=21, disorderThreshold?=0.542, minRegionLength?=5) → `residuePredictions`, `disorderedRegions[{start,end,meanScore,confidence,regionType}]`, `overallDisorderContent`.
2. **[MCP]** `predict_morfs`(sequence, minLength?=10, maxLength?=25) → MoRF intervals within disordered context.
- **[C# API]** `DisorderPredictor.PredictDisorder(...)` → `.PredictMoRFs(...)`.
- **Envelope:** report region **boundaries** + `overallDisorderContent`; treat `confidence` as uncalibrated (DISORDER-REGION-001). If a validated per-region call is required, use the `DisorderPredictor.PredictDisorderRegions` branch and **STOP** rather than raising the mode (below).
- **Cross-check:** an IDR should overlap a hydrophobicity **minimum** from pipeline (e) and (usually) a protein SEG low-complexity region.

### (b) Membrane / secretion topology → signal peptide + TM helices
1. **[MCP]** `predict_signal_peptide`(proteinSequence, prokaryote?=false) → signal-peptide presence + cleavage site.
2. **[MCP]** `predict_transmembrane_helices`(proteinSequence, windowSize?=19, threshold?=1.6) → TM segments + in/out topology.
- **[C# API]** `ProteinMotifFinder.PredictSignalPeptide(...)` → `.PredictTransmembraneHelices(...)`.
- **Cross-check:** each TM segment should sit on a hydrophobicity **peak** from pipeline (e); a called signal peptide occupies the N-terminus and should not be double-counted as TM1.

### (c) Structural / functional motifs → coiled-coils + domains + PROSITE
1. **[MCP]** `predict_coiled_coils`(proteinSequence, windowSize?=28, threshold?=0.5) → coiled-coil intervals + score.
2. **[MCP]** `find_protein_domains`(proteinSequence) → domain calls; `find_protein_motifs`(proteinSequence) → common motifs.
3. **[MCP]** targeted PROSITE: `find_motif_by_prosite`(proteinSequence, prositePattern, motifName?) → hits; `prosite_to_regex`(prositePattern) inspects the compiled regex first.
- **[C# API]** `ProteinMotifFinder.PredictCoiledCoils(...)` · `.FindDomains(...)` · `.FindCommonMotifs(...)` · `.FindMotifByProsite(...)` · `.ConvertPrositeToRegex(...)`.
- **Cross-check:** confirm a PROSITE hit's residue span against `prosite_to_regex` before interpreting.

### (d) Protein low-complexity (SEG) masking before motif interpretation
1. **[MCP]** `predict_low_complexity_seg`(sequence, triggerWindow?=12, triggerThreshold?=2.2, extensionThreshold?=2.5, minLength?=1) **or** `find_protein_low_complexity_regions`(proteinSequence, windowSize?=12, triggerComplexity?=2.2, extensionComplexity?=2.5) → SEG intervals.
- **[C# API]** `DisorderPredictor.PredictLowComplexityRegions(...)` **or** `ProteinMotifFinder.FindLowComplexityRegions(...)`.
- Both are protein SEG. Use as a filter: down-weight domain/PROSITE hits that fall inside SEG regions (compositionally biased, prone to spurious matches).

### (e) Hydrophobicity profile (supporting evidence)
1. **[MCP]** `hydrophobicity_profile`(proteinSequence, windowSize?=9) → per-window Kyte–Doolittle score.
- **[C# API]** `SequenceStatistics.CalculateHydrophobicityProfile(...)`.
- Feeds cross-checks: peaks corroborate TM helices (b); troughs corroborate IDRs (a).

## Envelope — STOP rule (guarded unit in scope)

- **DISORDER-REGION-001** (`predict_disorder` / `DisorderPredictor.PredictDisorder`, output the per-residue /
  per-region **confidence**) is **Permissive-only** (blocked in Strict & Moderate; default is Moderate). No
  disorder predictor publishes a calibrated confidence standard. If a task requires a **validated per-region
  confidence** and the call would run below its MinimumMode (throws `SeqeronLimitationException`), **STOP** and
  report the envelope — do **not** raise the mode to force output (bio-rigor rule). **Named alternative:** the
  region **boundaries** are validated (TOP-IDP threshold), so report those + `overallDisorderContent` via
  `DisorderPredictor.PredictDisorderRegions`; treat `confidence` as advisory only. The other tools in this
  family (signal peptide, TM, coiled-coil, domains, motifs, SEG, hydrophobicity) are **not** guarded.

## End-to-end grounded example — "characterize this protein"

**Task.** Given one protein sequence, (1) map intrinsically disordered regions, (2) test for a signal peptide,
(3) find transmembrane helices, (4) call coiled-coils + domains, (5) mask protein low-complexity, (6) profile
hydrophobicity to corroborate (1) and (3).

Tool / `Method ID` chain (MCP names; C# path in parentheses):
1. `predict_disorder`(sequence, minRegionLength=5) → `disorderedRegions` (report boundaries + `overallDisorderContent`; confidence advisory). (`DisorderPredictor.PredictDisorder`)
2. `predict_signal_peptide`(proteinSequence, prokaryote=false) → N-terminal signal + cleavage site. (`ProteinMotifFinder.PredictSignalPeptide`)
3. `predict_transmembrane_helices`(proteinSequence, threshold=1.6) → TM segments + topology. (`ProteinMotifFinder.PredictTransmembraneHelices`)
4. `predict_coiled_coils`(proteinSequence) + `find_protein_domains`(proteinSequence) → coiled-coil + domain calls. (`ProteinMotifFinder.PredictCoiledCoils` / `.FindDomains`)
5. `find_protein_low_complexity_regions`(proteinSequence) → SEG intervals to down-weight motif hits. (`ProteinMotifFinder.FindLowComplexityRegions`)
6. `hydrophobicity_profile`(proteinSequence, windowSize=9) → peaks over TM (step 3), troughs over IDRs (step 1). (`SequenceStatistics.CalculateHydrophobicityProfile`)

```
Provenance
1) predict_disorder(sequence, minRegionLength=5) → disorderedRegions [0-based, boundaries validated; confidence advisory], overallDisorderContent
2) predict_signal_peptide(proteinSequence, prokaryote=false) → signal peptide + cleavage site
3) predict_transmembrane_helices(proteinSequence, threshold=1.6) → TM segments + topology
4) predict_coiled_coils(proteinSequence); find_protein_domains(proteinSequence) → coiled-coil + domains
5) find_protein_low_complexity_regions(proteinSequence) → SEG intervals (mask/down-weight)
6) hydrophobicity_profile(proteinSequence, windowSize=9) → peaks corroborate TM (3), troughs corroborate IDRs (1)
Coordinates: 0-based residue indices; regions [start,end] inclusive. Units: normalized scores / bits / kcal-scale KD.
Envelope: DISORDER-REGION-001 — report validated boundaries + overallDisorderContent; confidence uncalibrated (STOP if a validated per-region confidence is required).
Caveat: alpha software; not for clinical use — independently validate before any decision use.
```

## Reference

- **Tool map (this family's protein tools, one-liners + Method ID + doc) — the curated index (no generated slice):** [`reference/tool-map.md`](reference/tool-map.md)
- **Fuller recipes + parameter guidance + envelope STOP rule:** [`reference/pipelines.md`](reference/pipelines.md)
- **Algorithm background (invariants/formulas — link, don't copy):**
  [`Disorder_Prediction.md`](../../../docs/algorithms/ProteinPred/Disorder_Prediction.md) ·
  [`Disordered_Region_Detection.md`](../../../docs/algorithms/ProteinPred/Disordered_Region_Detection.md) ·
  [`MoRF_Prediction.md`](../../../docs/algorithms/ProteinPred/MoRF_Prediction.md) ·
  [`Signal_Peptide_Prediction.md`](../../../docs/algorithms/ProteinMotif/Signal_Peptide_Prediction.md) ·
  [`Transmembrane_Helix_Prediction.md`](../../../docs/algorithms/ProteinMotif/Transmembrane_Helix_Prediction.md) ·
  [`Coiled_Coil_Prediction.md`](../../../docs/algorithms/ProteinMotif/Coiled_Coil_Prediction.md) ·
  [`Domain_Prediction.md`](../../../docs/algorithms/ProteinMotif/Domain_Prediction.md) ·
  [`PROSITE_Pattern_Matching.md`](../../../docs/algorithms/ProteinMotif/PROSITE_Pattern_Matching.md) ·
  [`Low_Complexity_Region_Detection.md`](../../../docs/algorithms/ProteinMotif/Low_Complexity_Region_Detection.md)
- **Operating envelope / guarded units:** [`LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) (DISORDER-REGION-001).
- **Cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor guardrail) ·
  [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup) ·
  [`bio-annotation`](../bio-annotation/SKILL.md) (DNA-shaped annotation & DNA motifs — overlap boundary above).
