# MolTools tool map (55 tools)

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
| `crispr_specificity_score` | Aggregate off-targets → 0..100 specificity (this library's own heuristic, not a published model) | `CrisprDesigner.CalculateHeuristicSpecificityScore` |
| `crispr_mit_specificity_score` | Published MIT / Hsu 2013 aggregate guide specificity (CRISPOR `calcMitGuideScore`) | `CrisprDesigner.CalculateMitSpecificityScore` |
| `crispr_cfd_specificity_score` | Published Doench 2016 CFD aggregate guide specificity (CRISPOR) | `CrisprDesigner.CalculateCfdSpecificityScore` |
| `calculate_on_target_doench2014` | Doench 2014 Rule Set 1 on-target score (0..100) from the 30-nt context | `CrisprDesigner.CalculateOnTargetDoench2014` |
| `calculate_on_target_rule_set2` | Doench 2016 Rule Set 2 / Azimuth on-target score (0..1) from the 30-nt context | `CrisprDesigner.CalculateOnTargetRuleSet2` |
| `crispr_system_info` | System parameters (PAM, guide length…) | `CrisprDesigner.GetSystem` |

## Codon optimization & usage (`CodonOptimizer`, `CodonUsageAnalyzer`)

| Tool | Purpose | Method ID |
|---|---|---|
| `optimize_codons` | Optimize a CDS for a host (CAI+GC before/after) | `CodonOptimizer.OptimizeSequence` |
| `build_codon_table` | Codon usage table from a reference sequence | `CodonOptimizer.CreateCodonTableFromSequence` |
| `compare_codon_usage` | Compare two codon-usage profiles | `CodonOptimizer.CompareCodonUsage` |
| `find_rare_codons` | Locate rare codons below a threshold | `CodonOptimizer.FindRareCodons` |
| `find_rare_codon_clusters` | Rare-codon clusters (Sherlocc: ≥4 rare in a 7-codon window) | `CodonOptimizer.FindRareCodonClusters` |
| `codon_minmax_profile` | %MinMax sliding-window codon-usage profile (Clarke & Clark) | `CodonOptimizer.CalculateMinMaxProfile` |
| `find_conserved_rare_codon_clusters` | Conserved rare-codon clusters of a codon alignment (Sherlocc 1.0) | `CodonOptimizer.FindConservedRareCodonClusters` |
| `reduce_secondary_structure` | Rewrite codons to reduce mRNA structure | `CodonOptimizer.ReduceSecondaryStructure` |
| `remove_restriction_sites` | Silently remove restriction sites from a CDS | `CodonOptimizer.RemoveRestrictionSites` |
| `cai_from_organism_table` | CAI vs an organism table | `CodonOptimizer.CalculateCAI` |
| `codon_adaptation_index` | CAI of a sequence (Sharp & Li) | `CodonUsageAnalyzer.CalculateCai` |
| `codon_usage_statistics` | Codon usage summary stats | `CodonUsageAnalyzer.GetStatistics` |
| `count_codons` | Codon counts | `CodonUsageAnalyzer.CountCodons` |
| `rscu` | Relative synonymous codon usage | `CodonUsageAnalyzer.CalculateRscu` |
| `effective_number_of_codons` | ENC (Nc) | `CodonUsageAnalyzer.CalculateEnc` |

## Restriction analysis (`RestrictionAnalyzer`)

Every sequence-taking tool here accepts **IUPAC ambiguity codes** anywhere in the sequence, so an
N-containing contig or a soft-masked assembly can be scanned, digested or mapped as it is (B08 H27).

| Tool | Purpose | Method ID |
|---|---|---|
| `find_restriction_sites` | Sites for named enzymes | `RestrictionAnalyzer.FindSites` |
| `find_all_restriction_sites` | Sites for every enzyme of a set (`enzyme_set`: 39 common by default, `"rebase"` = 754) | `RestrictionAnalyzer.FindAllSites` |
| `restriction_digest` | Simulate a digest (`topology` Linear \| Circular) → fragments | `RestrictionAnalyzer.Digest` |
| `partial_digest` | Partial digest: every fragment with any subset of sites uncut (+ missed-cut count; no yield) | `RestrictionAnalyzer.PartialDigest` |
| `restriction_map` | Build a restriction map | `RestrictionAnalyzer.CreateMap` |
| `digest_summary` | Summarize a digest (fragment sizes…) | `RestrictionAnalyzer.GetDigestSummary` |
| `get_enzyme` | Look up one enzyme's properties | `RestrictionAnalyzer.GetEnzyme` |
| `blunt_cutters` | Enzymes producing blunt ends | `RestrictionAnalyzer.GetBluntCutters` |
| `sticky_cutters` | Enzymes producing sticky ends | `RestrictionAnalyzer.GetStickyCutters` |
| `enzymes_by_cut_length` | Enzymes by recognition-site length; `max_length` makes it an inclusive range, so "all 4-6 bp cutters" is one call | `RestrictionAnalyzer.GetEnzymesByCutLength` |
| `find_pattern_occurrences` | The raw IUPAC pattern scan under `find_restriction_sites`: start indices on one strand, no enzyme and no cut arithmetic | `RestrictionAnalyzer.FindPatternOccurrences` |
| `find_sites_custom_enzyme` | Sites of a caller-supplied enzyme definition (site + the two cut offsets, optional Type IIB second pair) — for an enzyme not in the table, or one of the 334 REBASE enzymes with no defined cut | `RestrictionAnalyzer.FindSites` |
| `restriction_digest_custom_enzyme` | Digest with a caller-supplied enzyme definition | `RestrictionAnalyzer.Digest` |
| `compatible_enzymes` | Enzymes compatible with a given one | `RestrictionAnalyzer.FindCompatibleEnzymes` |
| `enzymes_compatible` | Are two enzymes compatible? | `RestrictionAnalyzer.AreCompatible` |

## Probe / oligo design (`ProbeDesigner`) — see PROBE-DESIGN-001 envelope

| Tool | Purpose | Method ID |
|---|---|---|
| `design_probes` | Ranked hybridization probes for a target | `ProbeDesigner.DesignProbes` |
| `design_probes_primer3` | Primer3 `pick_hyb_probe_only` probe picker (primer3-py parity) | `ProbeDesigner.DesignProbesPrimer3` |
| `design_tiling_probes` | Tiling probe set across a target | `ProbeDesigner.DesignTilingProbes` |
| `design_antisense_probes` | Antisense probes | `ProbeDesigner.DesignAntisenseProbes` |
| `design_molecular_beacon` | Molecular-beacon probe design | `ProbeDesigner.DesignMolecularBeacon` |
| `validate_probe` | Validate a probe: k-mismatch hits, Primer3 ntthal self-structure, Kane (2000) non-target cross-hybridization | `ProbeDesigner.ValidateProbe` |
| `analyze_oligo` | Oligo properties (GC, Tm, structure) | `ProbeDesigner.AnalyzeOligo` |
| `oligo_extinction_coefficient` | Extinction coefficient | `ProbeDesigner.CalculateExtinctionCoefficient` |
| `oligo_concentration_from_absorbance` | Concentration from A260 | `ProbeDesigner.CalculateConcentration` |

**Envelope:** `PROBE-DESIGN-001` (MinimumMode **Moderate**) omits quantitative MGB ΔTm and
dual-quencher labelling — no closed form. STOP and report if the assay needs those; see
[`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md).
