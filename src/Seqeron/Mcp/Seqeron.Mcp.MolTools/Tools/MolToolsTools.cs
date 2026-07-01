using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;

namespace Seqeron.Mcp.MolTools.Tools;

[McpServerToolType]
public class MolToolsTools
{
    #region PrimerDesigner

    [McpServerTool(Name = "design_primers", Title = "MolTools — Design PCR Primer Pair", ReadOnly = true), Description("Designs forward and reverse PCR primers flanking a target region in a DNA template; picks the highest-scoring valid candidates from a 200 bp flanking window on each side and reports product size and pair compatibility. target_start/target_end are 0-based; the region must satisfy 0 <= target_start < target_end < template.Length.")]
    public static PrimerPairResult design_primers(
        [Description("DNA template (A/C/G/T).")] string template,
        [Description("0-based inclusive start of target region.")] int target_start,
        [Description("0-based inclusive end of target region.")] int target_end,
        [Description("Optional primer design parameters (lengths, GC%, Tm, repeats, GC-clamp/3' stability checks). Defaults are used if null.")] PrimerParameters? parameters = null)
    {
        if (string.IsNullOrEmpty(template))
            throw new System.ArgumentException("Template cannot be null or empty.", nameof(template));
        if (target_start < 0)
            throw new System.ArgumentException("Target start must be non-negative.", nameof(target_start));
        if (target_end >= template.Length)
            throw new System.ArgumentException("Target end must be within the template.", nameof(target_end));
        if (target_start >= target_end)
            throw new System.ArgumentException("Target start must be strictly less than target end.", nameof(target_start));

        return PrimerDesigner.DesignPrimers(new DnaSequence(template), target_start, target_end, parameters);
    }

    [McpServerTool, Description("Evaluates a single primer sequence against quality criteria: length, GC%, Tm, homopolymer/dinucleotide-repeat runs, hairpin potential, 3' end stability, optional 3' GC clamp.")]
    public static PrimerCandidate evaluate_primer(
        [Description("Primer sequence to evaluate.")] string sequence,
        [Description("0-based location of the primer in the template (informational).")] int position,
        [Description("True if this is a forward primer; false for reverse.")] bool is_forward,
        [Description("Optional primer design parameters.")] PrimerParameters? parameters = null)
    {
        return PrimerDesigner.EvaluatePrimer(sequence, position, is_forward, parameters);
    }

    [McpServerTool, Description("Computes primer melting temperature using Wallace rule (<14 valid bases) or Marmur–Doty (≥14 bases). Non-ACGT characters are ignored.")]
    public static TmResult primer_melting_temperature(
        [Description("Primer sequence.")] string primer)
    {
        return new TmResult(PrimerDesigner.CalculateMeltingTemperature(primer));
    }

    [McpServerTool, Description("Primer Tm with Schildkraut–Lifson salt correction (16.6·log10([Na+]/1000)) added to Wallace/Marmur–Doty Tm.")]
    public static TmResult primer_melting_temperature_salt(
        [Description("Primer sequence.")] string primer,
        [Description("Na+ concentration in mM (default 50).")] double na_concentration = 50)
    {
        return new TmResult(PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, na_concentration));
    }

    [McpServerTool, Description("Length of the longest run of identical nucleotides in a sequence (case-insensitive).")]
    public static HomopolymerLengthResult longest_homopolymer(
        [Description("Nucleotide sequence.")] string sequence)
    {
        return new HomopolymerLengthResult(PrimerDesigner.FindLongestHomopolymer(sequence));
    }

    [McpServerTool, Description("Number of repeat units in the longest dinucleotide tandem repeat (e.g. ATAT…).")]
    public static DinucleotideRepeatResult longest_dinucleotide_repeat(
        [Description("Nucleotide sequence.")] string sequence)
    {
        return new DinucleotideRepeatResult(PrimerDesigner.FindLongestDinucleotideRepeat(sequence));
    }

    [McpServerTool, Description("Detects whether a sequence can form a hairpin with a stem of min_stem_length and a loop of at least min_loop_length. Uses O(n²) scan for short sequences and a suffix-tree based scan for sequences ≥100 bp.")]
    public static HairpinPotentialResult hairpin_potential(
        [Description("Nucleotide sequence.")] string sequence,
        [Description("Minimum stem length (default 4).")] int min_stem_length = 4,
        [Description("Minimum loop length (default 3).")] int min_loop_length = 3)
    {
        return new HairpinPotentialResult(PrimerDesigner.HasHairpinPotential(sequence, min_stem_length, min_loop_length));
    }

    [McpServerTool, Description("Heuristic 3' primer-dimer check between two primers using up-to-8-bp 3'-end complementarity. Returns the boolean dimer flag plus the count of complementary positions in the 3' window.")]
    public static PrimerDimerResult primer_dimer(
        [Description("First primer sequence.")] string primer1,
        [Description("Second primer sequence.")] string primer2,
        [Description("Minimum number of complementary 3'-end bases to flag a dimer (default 4).")] int min_complementarity = 4)
    {
        bool hasDimer = PrimerDesigner.HasPrimerDimer(primer1, primer2, min_complementarity);

        // Count complementary 3'-end positions (mirrors the inner loop in HasPrimerDimer).
        int complementary = 0;
        if (!string.IsNullOrEmpty(primer1) && !string.IsNullOrEmpty(primer2))
        {
            string seq1 = primer1.ToUpperInvariant();
            string seq2 = DnaSequence.GetReverseComplementString(primer2.ToUpperInvariant());
            int checkLength = System.Math.Min(8, System.Math.Min(seq1.Length, seq2.Length));
            string end1 = seq1.Substring(seq1.Length - checkLength);
            string end2 = seq2.Substring(0, checkLength);
            for (int i = 0; i < checkLength; i++)
            {
                if (IsComplementary(end1[i], end2[i]))
                    complementary++;
            }
        }
        return new PrimerDimerResult(hasDimer, complementary);

        static bool IsComplementary(char c1, char c2) =>
            (c1 == 'A' && c2 == 'T') || (c1 == 'T' && c2 == 'A') ||
            (c1 == 'G' && c2 == 'C') || (c1 == 'C' && c2 == 'G');
    }

    [McpServerTool, Description("SantaLucia (1998) ΔG°37 of the last 5 bases (NN parameters + initiation, 1 M NaCl), matching Primer3 PRIMER_MAX_END_STABILITY. More negative ΔG = more stable 3' end.")]
    public static ThreePrimeStabilityResult three_prime_stability(
        [Description("Primer sequence.")] string sequence)
    {
        return new ThreePrimeStabilityResult(PrimerDesigner.Calculate3PrimeStability(sequence));
    }

    [McpServerTool, Description("Enumerates all primer candidates of admissible lengths within a region of the template, evaluating each one. Useful when the caller wants to pick by custom criteria.")]
    public static PrimerCandidateListResult generate_primer_candidates(
        [Description("DNA template sequence.")] string template,
        [Description("0-based start of the search region (inclusive).")] int region_start,
        [Description("0-based end of the search region (exclusive).")] int region_end,
        [Description("True for forward-strand candidates, false for reverse-complement candidates.")] bool forward = true,
        [Description("Optional primer design parameters.")] PrimerParameters? parameters = null)
    {
        var list = PrimerDesigner
            .GeneratePrimerCandidates(new DnaSequence(template), region_start, region_end, forward, parameters)
            .ToList();
        return new PrimerCandidateListResult(list);
    }

    #endregion

    #region RestrictionAnalyzer

    [McpServerTool, Description("Looks up a built-in restriction enzyme by name (case-insensitive). Returns null when not found.")]
    public static EnzymeLookupResult get_enzyme(
        [Description("Enzyme name (e.g. EcoRI, BamHI).")] string name)
    {
        return new EnzymeLookupResult(RestrictionAnalyzer.GetEnzyme(name));
    }

    [McpServerTool, Description("Lists all built-in enzymes whose recognition sequence has the specified length (e.g. 4-cutters, 6-cutters, 8-cutters).")]
    public static EnzymeListResult enzymes_by_cut_length(
        [Description("Recognition-sequence length.")] int length)
    {
        return new EnzymeListResult(RestrictionAnalyzer.GetEnzymesByCutLength(length).ToList());
    }

    [McpServerTool, Description("Lists all built-in enzymes that produce blunt ends.")]
    public static EnzymeListResult blunt_cutters()
    {
        return new EnzymeListResult(RestrictionAnalyzer.GetBluntCutters().ToList());
    }

    [McpServerTool, Description("Lists all built-in enzymes that produce sticky (cohesive) ends.")]
    public static EnzymeListResult sticky_cutters()
    {
        return new EnzymeListResult(RestrictionAnalyzer.GetStickyCutters().ToList());
    }

    [McpServerTool, Description("Finds restriction sites for one or more named built-in enzymes on both strands of a DNA sequence. IUPAC degenerate codes in recognition sequences are matched against ACGT input.")]
    public static RestrictionSiteListResult find_restriction_sites(
        [Description("DNA sequence to scan.")] string sequence,
        [Description("Names of restriction enzymes to look for.")] string[] enzyme_names)
    {
        var sites = RestrictionAnalyzer.FindSites(new DnaSequence(sequence), enzyme_names).ToList();
        return new RestrictionSiteListResult(sites);
    }

    [McpServerTool, Description("Finds sites for every built-in restriction enzyme on both strands.")]
    public static RestrictionSiteListResult find_all_restriction_sites(
        [Description("DNA sequence to scan.")] string sequence)
    {
        var sites = RestrictionAnalyzer.FindAllSites(new DnaSequence(sequence)).ToList();
        return new RestrictionSiteListResult(sites);
    }

    [McpServerTool, Description("Simulates a restriction digest with one or more named enzymes and yields the resulting linear fragments in order.")]
    public static DigestResult restriction_digest(
        [Description("DNA sequence to digest.")] string sequence,
        [Description("Names of restriction enzymes to use (must be non-empty).")] string[] enzyme_names)
    {
        var fragments = RestrictionAnalyzer.Digest(new DnaSequence(sequence), enzyme_names).ToList();
        return new DigestResult(fragments);
    }

    [McpServerTool, Description("Aggregate statistics over a simulated digest: total fragment count, fragment sizes (descending), largest/smallest, average size, and enzymes used.")]
    public static DigestSummary digest_summary(
        [Description("DNA sequence to digest.")] string sequence,
        [Description("Names of restriction enzymes to use.")] string[] enzyme_names)
    {
        return RestrictionAnalyzer.GetDigestSummary(new DnaSequence(sequence), enzyme_names);
    }

    [McpServerTool, Description("Builds a restriction map: forward+reverse sites, sites grouped by enzyme, unique cutters (single forward-strand site), and non-cutters from the queried enzyme set. Empty enzyme_names considers all built-in enzymes.")]
    public static RestrictionMap restriction_map(
        [Description("DNA sequence to map.")] string sequence,
        [Description("Names of restriction enzymes to consider; empty to use every built-in enzyme.")] string[] enzyme_names)
    {
        return RestrictionAnalyzer.CreateMap(new DnaSequence(sequence), enzyme_names);
    }

    [McpServerTool, Description("Enumerates pairs of built-in enzymes whose ends can be ligated to each other (matching overhang type and sequence, or both blunt).")]
    public static CompatibleEnzymesResult compatible_enzymes()
    {
        var pairs = RestrictionAnalyzer.FindCompatibleEnzymes()
            .Select(t => new EnzymeCompatibilityPair(t.Enzyme1, t.Enzyme2, t.CompatibleEnd))
            .ToList();
        return new CompatibleEnzymesResult(pairs);
    }

    [McpServerTool, Description("Returns whether two named built-in enzymes have compatible ends. Unknown enzyme names yield false (no exception).")]
    public static EnzymeCompatibilityResult enzymes_compatible(
        [Description("First enzyme name.")] string enzyme1_name,
        [Description("Second enzyme name.")] string enzyme2_name)
    {
        return new EnzymeCompatibilityResult(RestrictionAnalyzer.AreCompatible(enzyme1_name, enzyme2_name));
    }

    #endregion

    #region CodonUsageAnalyzer

    [McpServerTool, Description("Counts occurrences of each ACGT codon in a coding DNA sequence (frame 0, multiples of 3). Trailing partial codons and non-ACGT codons are skipped silently.")]
    public static CodonCountsResult count_codons(
        [Description("Coding DNA sequence (frame 0).")] string sequence)
    {
        return new CodonCountsResult(CodonUsageAnalyzer.CountCodons(sequence));
    }

    [McpServerTool, Description("Relative Synonymous Codon Usage (RSCU) per codon: observed / expected-if-uniform among synonymous codons. RSCU=1 → no bias; >1 over-represented; <1 under-represented.")]
    public static RscuResult rscu(
        [Description("Coding DNA sequence (frame 0).")] string sequence)
    {
        return new RscuResult(CodonUsageAnalyzer.CalculateRscu(sequence));
    }

    [McpServerTool, Description("Codon Adaptation Index (Sharp & Li 1987) using a caller-supplied reference RSCU table (typically derived from highly expressed genes). Output range 0..1; codons absent from the reference contribute w=0 and are excluded from the geometric mean.")]
    public static CaiResult codon_adaptation_index(
        [Description("Coding DNA sequence (frame 0).")] string sequence,
        [Description("Reference RSCU table: codon (DNA alphabet) → RSCU value.")] Dictionary<string, double> reference_rscu)
    {
        return new CaiResult(CodonUsageAnalyzer.CalculateCai(sequence, reference_rscu));
    }

    [McpServerTool, Description("Effective Number of Codons (Wright Nc), measuring departure from uniform synonymous usage. Range 20..61 (clamped). Empty/very-short sequence → 0 because under-occupied amino-acid groups fall back to their degeneracy.")]
    public static EncResult effective_number_of_codons(
        [Description("Coding DNA sequence (frame 0).")] string sequence)
    {
        return new EncResult(CodonUsageAnalyzer.CalculateEnc(sequence));
    }

    [McpServerTool, Description("Aggregate codon-usage report: per-codon counts, RSCU, ENC, total codons, GC% at codon positions 1/2/3, GC3s, and overall GC.")]
    public static CodonUsageStatistics codon_usage_statistics(
        [Description("Coding DNA sequence (frame 0).")] string sequence)
    {
        return CodonUsageAnalyzer.GetStatistics(sequence);
    }

    #endregion

    #region CodonOptimizer

    [McpServerTool, Description("Optimizes a coding sequence for expression in a target organism using one of five strategies (MaximizeCAI, BalancedOptimization (default), HarmonizeExpression, MinimizeSecondary, AvoidRareCodeons). Internally trims to whole codons and converts T→U. Note: HarmonizeExpression is non-deterministic.")]
    public static OptimizationResultDto optimize_codons(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Target organism: a preset id (EColiK12 | Yeast | Human) or an inline custom table (organismName + codonFrequencies in RNA alphabet).")] CodonUsageTableInput target_organism,
        [Description("Optimization strategy.")] CodonOptimizer.OptimizationStrategy strategy = CodonOptimizer.OptimizationStrategy.BalancedOptimization,
        [Description("Lower bound for target GC fraction (default 0.40).")] double gc_target_min = 0.40,
        [Description("Upper bound for target GC fraction (default 0.60).")] double gc_target_max = 0.60,
        [Description("Frequency threshold below which a codon is considered rare (default 0.15).")] double rare_codon_threshold = 0.15)
    {
        var table = ResolveCodonUsageTable(target_organism);
        var result = CodonOptimizer.OptimizeSequence(coding_sequence, table, strategy, gc_target_min, gc_target_max, rare_codon_threshold);
        return new OptimizationResultDto(
            result.OriginalSequence,
            result.OptimizedSequence,
            result.ProteinSequence,
            result.OriginalCAI,
            result.OptimizedCAI,
            result.GcContentOriginal,
            result.GcContentOptimized,
            result.ChangedCodons,
            result.Changes.Select(c => new CodonChange(c.Position, c.Original, c.Optimized)).ToList());
    }

    [McpServerTool, Description("Computes CAI for a sequence against an organism CodonUsageTable (frequency table — distinct from codon_adaptation_index, which takes a reference RSCU dictionary). Codons without synonymous-group data are skipped; per-codon adaptiveness is clamped at 1e-6 to avoid ln(0) on incomplete custom tables.")]
    public static CaiResult cai_from_organism_table(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Target organism: preset id (EColiK12 | Yeast | Human) or inline custom table.")] CodonUsageTableInput target_organism)
    {
        var table = ResolveCodonUsageTable(target_organism);
        return new CaiResult(CodonOptimizer.CalculateCAI(coding_sequence, table));
    }

    [McpServerTool, Description("Synonymously rewrites codons to eliminate the listed restriction recognition sequences while preserving the protein. Site strings may be DNA or RNA; converted internally to RNA. Sites with no synonymous alternative are left in place.")]
    public static OptimizedSequenceResult remove_restriction_sites(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Restriction recognition sequences to eliminate.")] string[] restriction_sites,
        [Description("Target organism: preset id or inline custom table (used for synonymous-codon lookup).")] CodonUsageTableInput target_organism)
    {
        var table = ResolveCodonUsageTable(target_organism);
        return new OptimizedSequenceResult(CodonOptimizer.RemoveRestrictionSites(coding_sequence, restriction_sites, table));
    }

    [McpServerTool, Description("Greedy synonymous-codon swap that lowers a heuristic local self-complementarity score within a sliding window. Sequences shorter than the window are returned unchanged.")]
    public static OptimizedSequenceResult reduce_secondary_structure(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Target organism: preset id or inline custom table.")] CodonUsageTableInput target_organism,
        [Description("Sliding-window size in nucleotides (default 40).")] int window_size = 40)
    {
        var table = ResolveCodonUsageTable(target_organism);
        return new OptimizedSequenceResult(CodonOptimizer.ReduceSecondaryStructure(coding_sequence, table, window_size));
    }

    [McpServerTool, Description("Reports every codon whose frequency in the target table is below threshold.")]
    public static RareCodonsResult find_rare_codons(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Target organism: preset id or inline custom table.")] CodonUsageTableInput target_organism,
        [Description("Frequency threshold (default 0.15).")] double threshold = 0.15)
    {
        var table = ResolveCodonUsageTable(target_organism);
        var rare = CodonOptimizer.FindRareCodons(coding_sequence, table, threshold)
            .Select(t => new RareCodon(t.Position, t.Codon, t.AminoAcid, t.Frequency))
            .ToList();
        return new RareCodonsResult(rare);
    }

    [McpServerTool, Description("Codon-frequency similarity between two sequences: 1 − ½·Σ|f1−f2| ∈ [0,1] (1 = identical distribution, 0 = disjoint). Either input empty or with no codon hits → 0.")]
    public static SimilarityResult compare_codon_usage(
        [Description("First coding sequence.")] string sequence1,
        [Description("Second coding sequence.")] string sequence2)
    {
        return new SimilarityResult(CodonOptimizer.CompareCodonUsage(sequence1, sequence2));
    }

    [McpServerTool, Description("Derives a per-organism CodonUsageTable from a reference coding sequence by computing per-amino-acid relative codon frequencies (RNA alphabet).")]
    public static CodonUsageTableDto build_codon_table(
        [Description("Reference coding sequence (DNA or RNA).")] string reference_sequence,
        [Description("Organism name to attach to the resulting table.")] string organism_name)
    {
        var table = CodonOptimizer.CreateCodonTableFromSequence(reference_sequence, organism_name);
        return new CodonUsageTableDto(
            table.OrganismName,
            new Dictionary<string, double>(table.CodonFrequencies),
            new Dictionary<string, string>(table.CodonToAminoAcid));
    }

    /// <summary>
    /// Resolves a <see cref="CodonUsageTableInput"/> (preset id or inline custom
    /// table) to a <see cref="CodonOptimizer.CodonUsageTable"/>. For inline tables,
    /// supplies an empty <c>CodonToAminoAcid</c> dictionary because the underlying
    /// optimizer never reads that field (translation uses its private genetic-code map).
    /// </summary>
    private static CodonOptimizer.CodonUsageTable ResolveCodonUsageTable(CodonUsageTableInput input)
    {
        if (input is null)
            throw new System.ArgumentNullException(nameof(input));

        if (!string.IsNullOrEmpty(input.Preset))
        {
            return input.Preset switch
            {
                "EColiK12" => CodonOptimizer.EColiK12,
                "Yeast"    => CodonOptimizer.Yeast,
                "Human"    => CodonOptimizer.Human,
                _ => throw new System.ArgumentException(
                    $"Unknown codon-usage preset '{input.Preset}'. Expected: EColiK12 | Yeast | Human.",
                    nameof(input))
            };
        }

        if (input.CodonFrequencies is null)
            throw new System.ArgumentException(
                "CodonUsageTableInput must provide either a preset or a custom codonFrequencies table.",
                nameof(input));

        return new CodonOptimizer.CodonUsageTable(
            input.OrganismName ?? "Custom",
            new Dictionary<string, double>(input.CodonFrequencies),
            new Dictionary<string, string>());
    }

    #endregion

    #region CrisprDesigner

    [McpServerTool, Description("Returns the metadata record (PAM, guide length, PAM placement, description) for a known CRISPR system.")]
    public static CrisprSystem crispr_system_info(
        [Description("CRISPR system: SpCas9 | SpCas9_NAG | SaCas9 | Cas12a | AsCas12a | LbCas12a | CasX.")] CrisprSystemType system_type)
    {
        return CrisprDesigner.GetSystem(system_type);
    }

    [McpServerTool, Description("Finds all PAM matches (forward + reverse strand) for the chosen CRISPR system. PAM matching honours IUPAC codes (e.g. NGG, NNGRRT, TTTV). Sites whose target window falls outside the sequence are skipped.")]
    public static PamSitesResult find_pam_sites(
        [Description("DNA sequence to scan.")] string sequence,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9)
    {
        var sites = CrisprDesigner.FindPamSites(sequence, system_type).ToList();
        return new PamSitesResult(sites);
    }

    [McpServerTool, Description("Generates and scores guide-RNA candidates whose Cas9/Cas12a cut site falls inside the requested region. Candidates below parameters.minScore are filtered out. Region indices are 0-based; regionEnd is inclusive.")]
    public static GuideRnasResult design_guide_rnas(
        [Description("DNA sequence containing the target region.")] string sequence,
        [Description("0-based start of the target region.")] int region_start,
        [Description("0-based inclusive end of the target region.")] int region_end,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9,
        [Description("Optional guide-RNA design parameters (minGcContent, maxGcContent, minScore, avoidPolyT, checkSelfComplementarity). Defaults are used when null.")] GuideRnaParameters? parameters = null)
    {
        var guides = CrisprDesigner
            .DesignGuideRnas(new DnaSequence(sequence), region_start, region_end, system_type, parameters)
            .ToList();
        return new GuideRnasResult(guides);
    }

    [McpServerTool, Description("Scores a single guide RNA: GC%, seed-region GC%, polyT terminator presence, self-complementarity, common-restriction-site presence; returns a 0..100 score plus an issues list. Position is -1 for ad-hoc evaluation.")]
    public static GuideRnaCandidate evaluate_guide_rna(
        [Description("Guide RNA sequence.")] string guide_sequence,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9,
        [Description("Optional guide-RNA design parameters.")] GuideRnaParameters? parameters = null)
    {
        return CrisprDesigner.EvaluateGuideRna(guide_sequence, system_type, parameters);
    }

    [McpServerTool, Description("Naïve genome scan: enumerates all PAM sites in the genome and reports those whose target differs from the guide by ≤max_mismatches. Off-target score weights mismatches inside the seed region more heavily. O(genome × guide) — recommend genome ≤ ~1 Mb.")]
    public static OffTargetsResult find_off_targets(
        [Description("Guide RNA sequence (length must match the system's guide length).")] string guide_sequence,
        [Description("Genome / reference sequence to scan.")] string genome,
        [Description("Maximum allowed mismatches (range 0..5; default 3).")] int max_mismatches = 3,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9)
    {
        var hits = CrisprDesigner
            .FindOffTargets(guide_sequence, new DnaSequence(genome), max_mismatches, system_type)
            .ToList();
        return new OffTargetsResult(hits);
    }

    [McpServerTool, Description("Aggregates find_off_targets (≤4 mismatches) into a single specificity score in 0..100 (100 = no off-targets).")]
    public static SpecificityResult crispr_specificity_score(
        [Description("Guide RNA sequence (length must match the system's guide length).")] string guide_sequence,
        [Description("Genome / reference sequence to scan.")] string genome,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9)
    {
        return new SpecificityResult(
            CrisprDesigner.CalculateSpecificityScore(guide_sequence, new DnaSequence(genome), system_type));
    }

    #endregion

    #region ProbeDesigner

    [McpServerTool, Description("Designs hybridization probes by scanning the target for length-window candidates and ranking by GC%, Tm, homopolymers, self-complementarity, and structure heuristics. Use one of the ProbeParameters presets (Microarray | FISH | NorthernBlot | qPCR | SouthernBlot) or pass custom values; default = Microarray. Returns up to max_probes top-scoring probes.")]
    public static ProbesResult design_probes(
        [Description("Target DNA sequence.")] string target_sequence,
        [Description("Optional probe-design parameters (lengths, Tm range, GC range, max homopolymer, self-complementarity threshold). Defaults to Microarray when null.")] ProbeDesigner.ProbeParameters? parameters = null,
        [Description("Maximum probes to return (default 10).")] int max_probes = 10)
    {
        var probes = ProbeDesigner.DesignProbes(target_sequence, parameters, max_probes).ToList();
        return new ProbesResult(probes);
    }

    [McpServerTool, Description("Generates fixed-length probes covering the entire target with a configurable overlap. Sub-optimal candidates are still emitted (with a warning) so coverage is preserved.")]
    public static ProbeDesigner.TilingProbeSet design_tiling_probes(
        [Description("Target DNA sequence.")] string target_sequence,
        [Description("Tiling probe length in bp (default 60).")] int probe_length = 60,
        [Description("Overlap between adjacent probes in bp (default 20).")] int overlap = 20,
        [Description("Optional probe-design parameters (Tm/GC bounds etc.).")] ProbeDesigner.ProbeParameters? parameters = null)
    {
        return ProbeDesigner.DesignTilingProbes(target_sequence, probe_length, overlap, parameters);
    }

    [McpServerTool, Description("Reverse-complements the supplied mRNA-sense sequence and runs design_probes on it; tags returned probes with type=Antisense.")]
    public static ProbesResult design_antisense_probes(
        [Description("mRNA-sense sequence (will be reverse-complemented).")] string mrna_sequence,
        [Description("Optional probe-design parameters.")] ProbeDesigner.ProbeParameters? parameters = null,
        [Description("Maximum probes to return (default 5).")] int max_probes = 5)
    {
        var probes = ProbeDesigner.DesignAntisenseProbes(mrna_sequence, parameters, max_probes).ToList();
        return new ProbesResult(probes);
    }

    [McpServerTool, Description("Designs a hairpin molecular-beacon probe (GC-rich complementary stems flanking a target-specific loop) for real-time detection. Returns probe=null when target is shorter than probe_length or no acceptable loop is found.")]
    public static MolecularBeaconResult design_molecular_beacon(
        [Description("Target DNA sequence.")] string target_sequence,
        [Description("Loop (target-specific) length in bp (default 25).")] int probe_length = 25,
        [Description("Stem length in bp (default 5).")] int stem_length = 5)
    {
        return new MolecularBeaconResult(
            ProbeDesigner.DesignMolecularBeacon(target_sequence, probe_length, stem_length));
    }

    [McpServerTool, Description("Validates a probe against a set of reference sequences using k-mismatch approximate matching. Reports off-target hits, self-complementarity, secondary-structure flag, and a 0..1 specificity score (1.0 = unique match, 1/N for N hits, 0 if no hits).")]
    public static ProbeDesigner.ProbeValidation validate_probe(
        [Description("Probe sequence to validate.")] string probe_sequence,
        [Description("Reference sequences to scan for off-target hits.")] string[] reference_sequences,
        [Description("Maximum allowed mismatches (default 3).")] int max_mismatches = 3,
        [Description("Self-complementarity warning threshold (default 0.3).")] double self_complementarity_threshold = 0.3)
    {
        return ProbeDesigner.ValidateProbe(probe_sequence, reference_sequences, max_mismatches, self_complementarity_threshold);
    }

    [McpServerTool, Description("Returns Tm, GC fraction, molecular weight (Da), and 260 nm extinction coefficient (M⁻¹·cm⁻¹) for a short oligonucleotide.")]
    public static OligoAnalysisResult analyze_oligo(
        [Description("Oligonucleotide sequence.")] string sequence)
    {
        var (tm, gc, mw, eps) = ProbeDesigner.AnalyzeOligo(sequence);
        return new OligoAnalysisResult(tm, gc, mw, eps);
    }

    [McpServerTool, Description("Sums per-base 260 nm molar extinction contributions. Unknown bases contribute the fallback constant 10000.")]
    public static ExtinctionCoefficientResult oligo_extinction_coefficient(
        [Description("Oligonucleotide sequence.")] string sequence)
    {
        return new ExtinctionCoefficientResult(ProbeDesigner.CalculateExtinctionCoefficient(sequence));
    }

    [McpServerTool, Description("Beer–Lambert concentration (µM) from absorbance at 260 nm, extinction coefficient, and path length. No input validation — extinction_coefficient = 0 yields infinity.")]
    public static ConcentrationResult oligo_concentration_from_absorbance(
        [Description("Absorbance at 260 nm (A260).")] double absorbance260,
        [Description("Extinction coefficient ε in M⁻¹·cm⁻¹.")] double extinction_coefficient,
        [Description("Path length in cm (default 1.0).")] double path_length = 1.0)
    {
        return new ConcentrationResult(
            ProbeDesigner.CalculateConcentration(absorbance260, extinction_coefficient, path_length));
    }

    #endregion
}
