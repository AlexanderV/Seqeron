---
name: bio-moldesign
description: >-
  Design and QC molecular-biology reagents with Seqeron's MolTools server (47
  tools) — PCR primer pairs and hybridization probes, CRISPR guide RNAs with
  off-target and specificity scoring, codon optimization of a CDS for a target
  organism (CAI/GC before-after), restriction-site finding and digest
  simulation, and melting-temperature (Tm) with salt correction. Use for
  "design primers for this template", "find/design CRISPR guides for this
  target and check off-targets", "codon-optimize this gene for E. coli / yeast
  / human", "where does EcoRI cut / simulate a digest", "design a qPCR/FISH
  probe", "compute the salt-corrected Tm", primer-dimer / hairpin / 3'-stability
  QC, and choosing compatible restriction enzymes. Dual-mode: MCP tool calls +
  equivalent Seqeron.Genomics C# Method IDs.
allowed-tools: Read, Grep, Glob, Bash
---

# bio-moldesign — molecular design & QC (MolTools)

Routing + orchestration layer over the **MolTools** MCP server (47 tools) and the
`Seqeron.Genomics.MolTools` C# namespace. It gives dual-mode recipes; it does **not**
restate schemas — read the per-tool doc for exact I/O, and delegate all rigor to `bio-rigor`.

- **Compute with tools, never by eye.** Parse/validate/compute via a tool or `Method ID`. See
  [`bio-rigor`](../bio-rigor/SKILL.md) for the six rules (tool-only, envelope, provenance, cross-check, units, disclaimer).
- **Don't guess tool names.** Confirm with `python3 scripts/skills/find-tool.py <kw> --server moltools`
  or via [`seqeron-discovery`](../seqeron-discovery/SKILL.md).
- Coordinates in MolTools tools are **0-based** (target/region start/end); Tm is **°C**, GC is **%** (or fraction in codon tools).

## Envelope / STOP rule

One MolTools unit is **guarded**: `PROBE-DESIGN-001` (probe design), MinimumMode **Moderate**
(default). It is *correct-but-incomplete*: it does **not** model quantitative **MGB ΔTm** or
dual-quencher labelling (empirical/proprietary, no closed form). If a task needs a real MGB /
dual-quencher probe Tm, **STOP and report the limitation** — do not force a number. See
[`docs/Validation/LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) (row PROBE-DESIGN-001)
and `bio-rigor` [`reference/envelope.md`](../bio-rigor/reference/envelope.md).

## Decision guide — task → tool family

| You need to… | Family | Entry tool(s) |
|---|---|---|
| Design a PCR primer **pair** for a target region | primer | `design_primers` → QC below |
| Enumerate / score candidate primers, or QC one primer | primer | `generate_primer_candidates`, `evaluate_primer` |
| Primer QC: dimer / hairpin / 3′ stability / repeats | primer | `primer_dimer`, `hairpin_potential`, `three_prime_stability`, `longest_homopolymer`, `longest_dinucleotide_repeat` |
| Melting temperature (basic / salt-corrected) | primer | `primer_melting_temperature`, `primer_melting_temperature_salt` |
| Design **CRISPR** guide RNAs for a region | crispr | `design_guide_rnas` → `find_off_targets` → `crispr_specificity_score` |
| Find PAM sites / evaluate one guide / system info | crispr | `find_pam_sites`, `evaluate_guide_rna`, `crispr_system_info` |
| **Codon-optimize** a CDS for an organism | codon | `optimize_codons` (reports CAI+GC before/after) |
| Codon usage analysis (CAI, RSCU, ENC, rare codons) | codon | `codon_adaptation_index`, `rscu`, `effective_number_of_codons`, `find_rare_codons` |
| **Restriction** sites / digest / enzyme choice | restriction | `find_restriction_sites`, `restriction_digest`, `digest_summary`, `compatible_enzymes` |
| Design / validate a hybridization **probe** | probe | `design_probes`, `validate_probe`, `analyze_oligo` (see STOP rule for MGB) |

Full per-tool map: [`reference/tool-map.md`](reference/tool-map.md). Fuller recipes + parameter
guidance: [`reference/pipelines.md`](reference/pipelines.md).

## Canonical dual-mode pipelines

Each pipeline is ordered; MCP tool = MolTools server, C# = `Seqeron.Genomics.MolTools`.

### (a) PCR primer pair design + QC
Goal: a valid, well-behaved primer pair flanking a target, with a Tm-difference report.
- **[MCP]** `design_primers`(template, target_start, target_end) → for each of fwd/rev:
  `primer_melting_temperature_salt`, `primer_dimer`(fwd,rev), `hairpin_potential`, `three_prime_stability`.
- **[C# API]** `PrimerDesigner.DesignPrimers` → `PrimerDesigner.CalculateMeltingTemperatureWithSalt`,
  `PrimerDesigner.HasPrimerDimer`, `PrimerDesigner.HasHairpinPotential`, `PrimerDesigner.Calculate3PrimeStability`.
- Key params: 0-based inclusive `target_start`/`target_end`; defaults 18–25 nt, 40–60% GC, 57–63 °C;
  pair "compatible" when ΔTm ≤ 5 °C and no 3′ dimer. Report `|Tm_fwd − Tm_rev|` (°C, 1 dp).
- **Provenance:** design_primers → primer_melting_temperature_salt(×2) → primer_dimer → hairpin_potential, with all params.

### (b) CRISPR guide design → off-targets → specificity → rank
Goal: rank guide RNAs targeting a region by on-target score and off-target specificity.
- **[MCP]** `design_guide_rnas`(sequence, region_start, region_end, system_type) → for each guide:
  `find_off_targets`(guide, genome, max_mismatches) → `crispr_specificity_score`(guide, genome) → rank by (score, specificity).
- **[C# API]** `CrisprDesigner.DesignGuideRnas` → `CrisprDesigner.FindOffTargets` → `CrisprDesigner.CalculateSpecificityScore`.
- Key params: 0-based inclusive region; default system `SpCas9`; `max_mismatches` 0..5 (default 3);
  specificity uses maxMismatches=4, seed = PAM-proximal 12 bp. Keep genome ≲ 1 Mb (O(genome×guide)).
- **Provenance:** design_guide_rnas → find_off_targets → crispr_specificity_score per guide, + system + params.

### (c) Codon-optimize a CDS for a target organism
Goal: raise CAI for the host while staying in a GC band; report CAI/GC before→after.
- **[MCP]** `optimize_codons`(coding_sequence, target_organism, strategy, gc_target_min/max) → returns
  original/optimized CAI + GC + changed codons; optionally `find_rare_codons` first to inspect.
- **[C# API]** `CodonOptimizer.OptimizeSequence`; standalone CAI via `CodonUsageAnalyzer.CalculateCai`
  or `CodonOptimizer.CalculateCAI` (organism table).
- Key params: presets `EColiK12`/`Yeast`/`Human` (or custom RNA table); strategy default `BalancedOptimization`
  (`MaximizeCAI` deterministic; `HarmonizeExpression` **non-deterministic** — set seed/repeat); GC band default 0.40–0.60.
- **Provenance:** optimize_codons(preset, strategy, GC band) → report ΔCAI, ΔGC, changedCodons.

### (d) Restriction analysis: find sites / simulate digest
Goal: locate cut sites for chosen enzymes and produce fragment list.
- **[MCP]** `find_restriction_sites`(sequence, enzyme_names) → `restriction_digest`(sequence, enzyme_names)
  → `digest_summary`; choose enzymes with `compatible_enzymes` / `get_enzyme`.
- **[C# API]** `RestrictionAnalyzer.FindSites` → `RestrictionAnalyzer.Digest` → `RestrictionAnalyzer.GetDigestSummary`.
- Key params: linear molecule; `k` distinct forward cuts → `k+1` fragments. Site positions are 0-based.
- **Provenance:** find_restriction_sites(enzymes) → restriction_digest(enzymes) → fragment table.

### (e) Probe design + salt-corrected Tm
Goal: rank hybridization probes for a target and report a defensible Tm.
- **[MCP]** `design_probes`(target_sequence, parameters=preset, max_probes) → `validate_probe` /
  `analyze_oligo` per candidate → `primer_melting_temperature_salt` for a salt-aware Tm.
- **[C# API]** `ProbeDesigner.DesignProbes` → `ProbeDesigner.ValidateProbe` / `ProbeDesigner.AnalyzeOligo`
  → `PrimerDesigner.CalculateMeltingTemperatureWithSalt`.
- Key params: presets `Microarray`(default)/`FISH`/`NorthernBlot`/`qPCR`/`SouthernBlot`; salt Tm adds
  `16.6·log10([Na+]/1000)`. **STOP** if the assay needs MGB ΔTm / dual-quencher (PROBE-DESIGN-001 envelope).
- **Provenance:** design_probes(preset) → validate_probe → primer_melting_temperature_salt([Na+]).

## End-to-end grounded example (extends the README primer-QC workflow)

The README shows validate-DNA → GC% → Tm → ΔTm for a *given* primer pair
([`docs/mcp/README.md`](../../../docs/mcp/README.md) §"PCR primer QC"). Extend it into a full
**design + QC** task: *"Design a primer pair flanking target [100,150) in this template and QC it."*

Tool / Method-ID chain (MCP → C#):
1. `design_primers`(template, 100, 150) → `PrimerDesigner.DesignPrimers` — yields best fwd + rev, `isValid`, `productSize`.
2. `primer_melting_temperature_salt`(fwd) and (rev) → `PrimerDesigner.CalculateMeltingTemperatureWithSalt` — Tm at chosen [Na+].
3. `primer_dimer`(fwd, rev) → `PrimerDesigner.HasPrimerDimer` — cross-dimer check.
4. `hairpin_potential`(fwd) and (rev) → `PrimerDesigner.HasHairpinPotential`.

Expected-shape output (Markdown table + a compatibility line):

```
| id | length | gc_percent | tm_c | dimer | hairpin |
|---|---:|---:|---:|:---:|:---:|
| FWD | 25 | 52.00 | 59.3 | no | no |
| REV | 25 | 52.00 | 59.3 | no | no |

product_size_bp = 180 ; tm_diff_c = 0.0 ; pair_compatible = true
```
Then a **Provenance** block listing tools/Method-IDs in call order with params, plus the
alpha / not-for-clinical-use caveat from `bio-rigor` if the assay is decision-relevant.

## Reference

- Generated tool index for this domain: [`_generated/tools.md`](_generated/tools.md) *(built by
  `scripts/skills/gen-catalog.py`; do not hand-edit)*.
- Grouped tool map + Method IDs: [`reference/tool-map.md`](reference/tool-map.md). Fuller recipes:
  [`reference/pipelines.md`](reference/pipelines.md).
- Algorithm background (link, don't copy): [`docs/algorithms/Codon_Optimization/`](../../../docs/algorithms/Codon_Optimization),
  [`docs/algorithms/MolTools/`](../../../docs/algorithms/MolTools) (Primer_Design, Guide_RNA_Design,
  Off_Target_Analysis, Restriction_*, Hybridization_Probe_Design),
  [`docs/algorithms/MolTools/`](../../../docs/algorithms/MolTools) (Melting_Temperature, PAM_Site_Detection).
- Rigor + discovery: [`bio-rigor`](../bio-rigor/SKILL.md), [`seqeron-discovery`](../seqeron-discovery/SKILL.md).
