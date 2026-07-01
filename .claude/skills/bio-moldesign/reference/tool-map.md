# MolTools tool map (47 tools)

Grouped by sub-domain. All on the **MolTools** MCP server; C# in `Seqeron.Genomics.MolTools`.
One line = tool · one-line purpose · `Method ID`. Open the per-tool doc under
[`docs/mcp/tools/moltools/`](../../../../docs/mcp/tools/moltools) for the exact I/O schema — do not
guess. Source of truth for names/status: `docs/mcp/MCP_STATUS.md`.

## Primer design & QC (`PrimerDesigner`)

| Tool | Purpose | Method ID |
|---|---|---|
| `design_primers` | Best fwd/rev primer pair flanking a target | `PrimerDesigner.DesignPrimers` |
| `generate_primer_candidates` | Enumerate candidate primers in a window | `PrimerDesigner.GeneratePrimerCandidates` |
| `evaluate_primer` | Score/QC a single primer vs parameters | `PrimerDesigner.EvaluatePrimer` |
| `primer_dimer` | Detect a 3′ primer-dimer between two primers | `PrimerDesigner.HasPrimerDimer` |
| `hairpin_potential` | Detect hairpin secondary structure | `PrimerDesigner.HasHairpinPotential` |
| `three_prime_stability` | 3′-end stability (ΔG-style) | `PrimerDesigner.Calculate3PrimeStability` |
| `longest_homopolymer` | Longest single-base run | `PrimerDesigner.FindLongestHomopolymer` |
| `longest_dinucleotide_repeat` | Longest dinucleotide repeat | `PrimerDesigner.FindLongestDinucleotideRepeat` |
| `primer_melting_temperature` | Basic Wallace/Marmur–Doty Tm (°C) | `PrimerDesigner.CalculateMeltingTemperature` |
| `primer_melting_temperature_salt` | Salt-corrected Tm (Schildkraut–Lifson) | `PrimerDesigner.CalculateMeltingTemperatureWithSalt` |

## CRISPR guide design (`CrisprDesigner`)

| Tool | Purpose | Method ID |
|---|---|---|
| `design_guide_rnas` | Design + score guides for a region | `CrisprDesigner.DesignGuideRnas` |
| `evaluate_guide_rna` | Score/QC one guide | `CrisprDesigner.EvaluateGuideRna` |
| `find_pam_sites` | Enumerate PAM sites for a system | `CrisprDesigner.FindPamSites` |
| `find_off_targets` | Naïve genome scan for off-targets | `CrisprDesigner.FindOffTargets` |
| `crispr_specificity_score` | Aggregate off-targets → 0..100 specificity | `CrisprDesigner.CalculateSpecificityScore` |
| `crispr_system_info` | System parameters (PAM, guide length…) | `CrisprDesigner.GetSystem` |

## Codon optimization & usage (`CodonOptimizer`, `CodonUsageAnalyzer`)

| Tool | Purpose | Method ID |
|---|---|---|
| `optimize_codons` | Optimize a CDS for a host (CAI+GC before/after) | `CodonOptimizer.OptimizeSequence` |
| `build_codon_table` | Codon usage table from a reference sequence | `CodonOptimizer.CreateCodonTableFromSequence` |
| `compare_codon_usage` | Compare two codon-usage profiles | `CodonOptimizer.CompareCodonUsage` |
| `find_rare_codons` | Locate rare codons below a threshold | `CodonOptimizer.FindRareCodons` |
| `reduce_secondary_structure` | Rewrite codons to reduce mRNA structure | `CodonOptimizer.ReduceSecondaryStructure` |
| `remove_restriction_sites` | Silently remove restriction sites from a CDS | `CodonOptimizer.RemoveRestrictionSites` |
| `cai_from_organism_table` | CAI vs an organism table | `CodonOptimizer.CalculateCAI` |
| `codon_adaptation_index` | CAI of a sequence (Sharp & Li) | `CodonUsageAnalyzer.CalculateCai` |
| `codon_usage_statistics` | Codon usage summary stats | `CodonUsageAnalyzer.GetStatistics` |
| `count_codons` | Codon counts | `CodonUsageAnalyzer.CountCodons` |
| `rscu` | Relative synonymous codon usage | `CodonUsageAnalyzer.CalculateRscu` |
| `effective_number_of_codons` | ENC (Nc) | `CodonUsageAnalyzer.CalculateEnc` |

## Restriction analysis (`RestrictionAnalyzer`)

| Tool | Purpose | Method ID |
|---|---|---|
| `find_restriction_sites` | Sites for named enzymes | `RestrictionAnalyzer.FindSites` |
| `find_all_restriction_sites` | Sites for all known enzymes | `RestrictionAnalyzer.FindAllSites` |
| `restriction_digest` | Simulate a linear digest → fragments | `RestrictionAnalyzer.Digest` |
| `restriction_map` | Build a restriction map | `RestrictionAnalyzer.CreateMap` |
| `digest_summary` | Summarize a digest (fragment sizes…) | `RestrictionAnalyzer.GetDigestSummary` |
| `get_enzyme` | Look up one enzyme's properties | `RestrictionAnalyzer.GetEnzyme` |
| `blunt_cutters` | Enzymes producing blunt ends | `RestrictionAnalyzer.GetBluntCutters` |
| `sticky_cutters` | Enzymes producing sticky ends | `RestrictionAnalyzer.GetStickyCutters` |
| `enzymes_by_cut_length` | Enzymes by recognition-site length | `RestrictionAnalyzer.GetEnzymesByCutLength` |
| `compatible_enzymes` | Enzymes compatible with a given one | `RestrictionAnalyzer.FindCompatibleEnzymes` |
| `enzymes_compatible` | Are two enzymes compatible? | `RestrictionAnalyzer.AreCompatible` |

## Probe / oligo design (`ProbeDesigner`) — see PROBE-DESIGN-001 envelope

| Tool | Purpose | Method ID |
|---|---|---|
| `design_probes` | Ranked hybridization probes for a target | `ProbeDesigner.DesignProbes` |
| `design_tiling_probes` | Tiling probe set across a target | `ProbeDesigner.DesignTilingProbes` |
| `design_antisense_probes` | Antisense probes | `ProbeDesigner.DesignAntisenseProbes` |
| `design_molecular_beacon` | Molecular-beacon probe design | `ProbeDesigner.DesignMolecularBeacon` |
| `validate_probe` | Validate a probe vs parameters | `ProbeDesigner.ValidateProbe` |
| `analyze_oligo` | Oligo properties (GC, Tm, structure) | `ProbeDesigner.AnalyzeOligo` |
| `oligo_extinction_coefficient` | Extinction coefficient | `ProbeDesigner.CalculateExtinctionCoefficient` |
| `oligo_concentration_from_absorbance` | Concentration from A260 | `ProbeDesigner.CalculateConcentration` |

**Envelope:** `PROBE-DESIGN-001` (MinimumMode **Moderate**) omits quantitative MGB ΔTm and
dual-quencher labelling — no closed form. STOP and report if the assay needs those; see
[`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md).
