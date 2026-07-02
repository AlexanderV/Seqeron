# seqeron-protein-features — tool map (protein-feature family)

Server: **analysis**. Backing classes: `DisorderPredictor.*` (intrinsic disorder / MoRFs / propensity /
SEG), `ProteinMotifFinder.*` (signal peptide / TM / coiled-coil / domains / motifs / PROSITE / SEG), and
`SequenceStatistics` for the hydrophobicity profile. Input is always a **protein (amino-acid) sequence**.

> This skill is **not** in `domain-map.json`, so there is **no** generated `_generated/tools.md`. This
> curated table **is** the index. Verify every schema in `docs/mcp/tools/analysis/<tool>.md` before use.
>
> Scope: PROTEIN features only. DNA motifs (`MotifFinder.*`), nucleotide low-complexity/DUST
> (`SequenceComplexity.*`), and DNA/RNA annotation live in [`bio-annotation`](../../bio-annotation/SKILL.md) — not here.
>
> Coordinates: 0-based residue indices; regions `[start,end]` inclusive per the tool doc.

## Intrinsic disorder — `DisorderPredictor.*`

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `predict_disorder` ⚠ | `DisorderPredictor.PredictDisorder` | TOP-IDP per-residue disorder + IDRs (`start,end,meanScore,confidence,regionType`) + `overallDisorderContent`. **Confidence guarded (DISORDER-REGION-001).** | `predict_disorder.md` |
| `predict_morfs` | `DisorderPredictor.PredictMoRFs` | Molecular-recognition features within disordered context (`minLength`/`maxLength`). | `predict_morfs.md` |
| `disorder_propensity` | `DisorderPredictor.GetDisorderPropensity` | TOP-IDP normalized propensity of a single amino-acid letter. | `disorder_propensity.md` |
| `is_disorder_promoting` | `DisorderPredictor.IsDisorderPromoting` | Whether a single amino-acid letter is disorder-promoting. | `is_disorder_promoting.md` |
| `predict_low_complexity_seg` | `DisorderPredictor.PredictLowComplexityRegions` | Protein SEG low-complexity (trigger/extension entropy thresholds) — from the disorder engine. | `predict_low_complexity_seg.md` |

## Motif / domain / topology — `ProteinMotifFinder.*`

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `predict_signal_peptide` | `ProteinMotifFinder.PredictSignalPeptide` | N-terminal signal peptide + cleavage site (`prokaryote` matrix flag). | `predict_signal_peptide.md` |
| `predict_transmembrane_helices` | `ProteinMotifFinder.PredictTransmembraneHelices` | TM helices + membrane topology (`windowSize`, hydropathy `threshold`). | `predict_transmembrane_helices.md` |
| `predict_coiled_coils` | `ProteinMotifFinder.PredictCoiledCoils` | Coiled-coil regions + score (`windowSize`, `threshold`). | `predict_coiled_coils.md` |
| `find_protein_domains` | `ProteinMotifFinder.FindDomains` | Protein domain calls. | `find_protein_domains.md` |
| `find_protein_motifs` | `ProteinMotifFinder.FindCommonMotifs` | Common protein sequence motifs. | `find_protein_motifs.md` |
| `find_motif_by_prosite` | `ProteinMotifFinder.FindMotifByProsite` | Scan a PROSITE pattern in a protein (`prositePattern`, `motifName`). | `find_motif_by_prosite.md` |
| `prosite_to_regex` | `ProteinMotifFinder.ConvertPrositeToRegex` | Compile a PROSITE pattern to a regex (inspect before scanning). | `prosite_to_regex.md` |
| `find_protein_low_complexity_regions` | `ProteinMotifFinder.FindLowComplexityRegions` | Protein SEG low-complexity — from the motif engine (`triggerComplexity`/`extensionComplexity` bits). | `find_protein_low_complexity_regions.md` |

## Composition — `SequenceStatistics`

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `hydrophobicity_profile` | `SequenceStatistics.CalculateHydrophobicityProfile` | Kyte–Doolittle sliding-window hydrophobicity (`windowSize`). Peaks corroborate TM; troughs corroborate IDRs. | `hydrophobicity_profile.md` |
