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

    [McpServerTool, Description("Designs forward and reverse PCR primers flanking a target region in a DNA template; picks the highest-scoring valid candidates from a 200 bp flanking window on each side and reports product size and pair compatibility.")]
    public static PrimerPairResult design_primers(
        [Description("DNA template (A/C/G/T).")] string template,
        [Description("0-based inclusive start of target region.")] int target_start,
        [Description("0-based inclusive end of target region.")] int target_end,
        [Description("Optional primer design parameters (lengths, GC%, Tm, repeats, GC-clamp/3' stability checks). Defaults are used if null.")] PrimerParameters? parameters = null)
    {
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

    [McpServerTool(Name = "primer_melting_temperature", Title = "MolTools — Primer Melting Temperature", ReadOnly = true), Description("Computes a primer's melting temperature (Tm, °C): Wallace rule Tm = 2·(A+T) + 4·(G+C) for < 14 valid bases, or Marmur–Doty Tm = 64.9 + 41·(GC−16.4)/N for ≥ 14 valid bases. Non-ACGT characters are ignored. Call for a quick Tm estimate of a short oligo/primer.")]
    public static TmResult primer_melting_temperature(
        [Description("Primer sequence.")] string primer)
    {
        if (string.IsNullOrEmpty(primer))
            throw new System.ArgumentException("Primer cannot be null or empty.", nameof(primer));

        return new TmResult(PrimerDesigner.CalculateMeltingTemperature(primer));
    }

    [McpServerTool(Name = "primer_melting_temperature_salt", Title = "MolTools — Salt-Corrected Primer Tm", ReadOnly = true), Description("Primer Tm with a Schildkraut–Lifson salt correction: adds 16.6·log10([Na+]/1000) to the Wallace/Marmur–Doty Tm, rounded to one decimal. Call when a monovalent-cation ([Na+]) adjusted primer Tm is needed. Na+ concentration is in mM (default 50).")]
    public static TmResult primer_melting_temperature_salt(
        [Description("Primer sequence.")] string primer,
        [Description("Na+ concentration in mM (default 50).")] double na_concentration = 50)
    {
        if (string.IsNullOrEmpty(primer))
            throw new System.ArgumentException("Primer cannot be null or empty.", nameof(primer));
        if (na_concentration <= 0)
            throw new System.ArgumentException("Na+ concentration must be positive.", nameof(na_concentration));

        return new TmResult(PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, na_concentration));
    }

    [McpServerTool(Name = "longest_homopolymer", Title = "MolTools — Longest Homopolymer Run", ReadOnly = true), Description("Returns the length of the longest run of identical consecutive nucleotides (e.g. AAAA = 4) in a sequence, case-insensitive. Call to flag homopolymer stretches that hurt primer/probe quality.")]
    public static HomopolymerLengthResult longest_homopolymer(
        [Description("Nucleotide sequence.")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        return new HomopolymerLengthResult(PrimerDesigner.FindLongestHomopolymer(sequence));
    }

    [McpServerTool(Name = "longest_dinucleotide_repeat", Title = "MolTools — Longest Dinucleotide Repeat", ReadOnly = true), Description("Returns the number of repeat units in the longest dinucleotide tandem repeat (e.g. ATATAT = 3 units of AT), case-insensitive. Sequences shorter than 4 nt return 0. Call to flag microsatellite-like dinucleotide repeats in a primer/probe.")]
    public static DinucleotideRepeatResult longest_dinucleotide_repeat(
        [Description("Nucleotide sequence.")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        return new DinucleotideRepeatResult(PrimerDesigner.FindLongestDinucleotideRepeat(sequence));
    }

    [McpServerTool(Name = "hairpin_potential", Title = "MolTools — Hairpin Potential", ReadOnly = true), Description("Detects whether a sequence can fold into a hairpin: a self-complementary stem of at least min_stem_length separated by a loop of at least min_loop_length. Uses an O(n²) scan for short sequences and a suffix-tree scan for sequences ≥ 100 bp. Call to screen a primer/probe for secondary structure.")]
    public static HairpinPotentialResult hairpin_potential(
        [Description("Nucleotide sequence.")] string sequence,
        [Description("Minimum stem length (default 4).")] int min_stem_length = 4,
        [Description("Minimum loop length (default 3).")] int min_loop_length = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (min_stem_length <= 0)
            throw new System.ArgumentException("Minimum stem length must be positive.", nameof(min_stem_length));
        if (min_loop_length < 0)
            throw new System.ArgumentException("Minimum loop length cannot be negative.", nameof(min_loop_length));

        return new HairpinPotentialResult(PrimerDesigner.HasHairpinPotential(sequence, min_stem_length, min_loop_length));
    }

    [McpServerTool(Name = "primer_dimer", Title = "MolTools — Primer-Dimer Check", ReadOnly = true), Description("Heuristic 3'-end primer-dimer check between two primers: reverse-complements primer2 and counts complementary positions in an up-to-8-bp 3'-end window. Flags a dimer when at least min_complementarity positions are complementary. Returns the boolean flag plus the complementary-base count. Call to screen a primer pair for 3'-dimer formation.")]
    public static PrimerDimerResult primer_dimer(
        [Description("First primer sequence.")] string primer1,
        [Description("Second primer sequence.")] string primer2,
        [Description("Minimum number of complementary 3'-end bases to flag a dimer (default 4).")] int min_complementarity = 4)
    {
        if (string.IsNullOrEmpty(primer1))
            throw new System.ArgumentException("First primer cannot be null or empty.", nameof(primer1));
        if (string.IsNullOrEmpty(primer2))
            throw new System.ArgumentException("Second primer cannot be null or empty.", nameof(primer2));

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

    [McpServerTool(Name = "three_prime_stability", Title = "MolTools — Primer 3' End Stability (ΔG°37)", ReadOnly = true), Description("SantaLucia (1998) nearest-neighbor ΔG°37 (kcal/mol) of a primer's last 5 bases, including initiation terms (1 M NaCl), matching Primer3 PRIMER_MAX_END_STABILITY. More negative ΔG = a more stable (more problematic) 3' end. Sequences shorter than 5 bases return 0. Call to assess primer 3'-end stability for mispriming risk.")]
    public static ThreePrimeStabilityResult three_prime_stability(
        [Description("Primer sequence.")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

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

    [McpServerTool(Name = "get_enzyme", Title = "MolTools — Look Up Restriction Enzyme", ReadOnly = true), Description("Looks up a built-in restriction enzyme by name (case-insensitive) and returns its recognition sequence, cut positions and organism. Returns enzyme=null when the name is not in the built-in database. Call to fetch a single enzyme's properties by name.")]
    public static EnzymeLookupResult get_enzyme(
        [Description("Enzyme name (e.g. EcoRI, BamHI).")] string name)
    {
        if (string.IsNullOrWhiteSpace(name))
            throw new System.ArgumentException("Enzyme name cannot be null or blank.", nameof(name));

        return new EnzymeLookupResult(RestrictionAnalyzer.GetEnzyme(name));
    }

    [McpServerTool(Name = "enzymes_by_cut_length", Title = "MolTools — Enzymes by Recognition Length", ReadOnly = true), Description("Lists all built-in restriction enzymes whose recognition sequence has exactly the specified length in base pairs (e.g. 4 for frequent cutters, 6 for typical cloning enzymes, 8 for rare cutters). Call to pick enzymes by how often they cut.")]
    public static EnzymeListResult enzymes_by_cut_length(
        [Description("Recognition-sequence length in bp (must be positive).")] int length)
    {
        if (length <= 0)
            throw new System.ArgumentException("Recognition-sequence length must be positive.", nameof(length));

        return new EnzymeListResult(RestrictionAnalyzer.GetEnzymesByCutLength(length).ToList());
    }

    [McpServerTool(Name = "blunt_cutters", Title = "MolTools — Blunt-End Restriction Enzymes", ReadOnly = true), Description("Lists all built-in restriction enzymes that produce blunt ends (both strands cut at the same position). Call when the user wants blunt-cutting enzymes for blunt-end cloning.")]
    public static EnzymeListResult blunt_cutters()
    {
        return new EnzymeListResult(RestrictionAnalyzer.GetBluntCutters().ToList());
    }

    [McpServerTool(Name = "sticky_cutters", Title = "MolTools — Sticky-End Restriction Enzymes", ReadOnly = true), Description("Lists all built-in restriction enzymes that produce sticky (cohesive) ends — a staggered cut leaving a 5' or 3' single-stranded overhang. Call when the user wants overhang-generating enzymes for directional / sticky-end cloning.")]
    public static EnzymeListResult sticky_cutters()
    {
        return new EnzymeListResult(RestrictionAnalyzer.GetStickyCutters().ToList());
    }

    [McpServerTool(Name = "find_restriction_sites", Title = "MolTools — Find Restriction Sites", ReadOnly = true), Description("Finds restriction sites for one or more named built-in enzymes on both strands of a DNA sequence. IUPAC degenerate codes in recognition sequences are matched against ACGT input. Palindromic enzymes report a forward and a reverse site at the same position. Call to locate where specific enzymes cut a sequence.")]
    public static RestrictionSiteListResult find_restriction_sites(
        [Description("DNA sequence to scan.")] string sequence,
        [Description("Names of restriction enzymes to look for.")] string[] enzyme_names)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (enzyme_names is null || enzyme_names.Length == 0)
            throw new System.ArgumentException("At least one enzyme name is required.", nameof(enzyme_names));

        var sites = RestrictionAnalyzer.FindSites(new DnaSequence(sequence), enzyme_names).ToList();
        return new RestrictionSiteListResult(sites);
    }

    [McpServerTool(Name = "find_all_restriction_sites", Title = "MolTools — Find All Restriction Sites", ReadOnly = true), Description("Finds sites for EVERY built-in restriction enzyme on both strands of a DNA sequence. Call for a comprehensive restriction-site scan when the enzyme set is not known in advance.")]
    public static RestrictionSiteListResult find_all_restriction_sites(
        [Description("DNA sequence to scan.")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        var sites = RestrictionAnalyzer.FindAllSites(new DnaSequence(sequence)).ToList();
        return new RestrictionSiteListResult(sites);
    }

    [McpServerTool(Name = "restriction_digest", Title = "MolTools — Restriction Digest", ReadOnly = true), Description("Simulates a restriction digest of a linear DNA molecule with one or more named enzymes and yields the resulting fragments in 5'→3' order (each with its sequence, start, length and flanking enzymes). With k distinct forward-strand cut positions a linear molecule yields k+1 fragments. Call to predict the fragment pattern of a digest.")]
    public static DigestResult restriction_digest(
        [Description("DNA sequence to digest.")] string sequence,
        [Description("Names of restriction enzymes to use (must be non-empty).")] string[] enzyme_names)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (enzyme_names is null || enzyme_names.Length == 0)
            throw new System.ArgumentException("At least one enzyme name is required.", nameof(enzyme_names));

        var fragments = RestrictionAnalyzer.Digest(new DnaSequence(sequence), enzyme_names).ToList();
        return new DigestResult(fragments);
    }

    [McpServerTool(Name = "digest_summary", Title = "MolTools — Restriction Digest Summary", ReadOnly = true), Description("Aggregate statistics over a simulated linear restriction digest: total fragment count, fragment sizes (descending), largest/smallest fragment, average fragment size, and enzymes used. Call for a gel-like summary of a digest without the full fragment sequences.")]
    public static DigestSummary digest_summary(
        [Description("DNA sequence to digest.")] string sequence,
        [Description("Names of restriction enzymes to use.")] string[] enzyme_names)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (enzyme_names is null || enzyme_names.Length == 0)
            throw new System.ArgumentException("At least one enzyme name is required.", nameof(enzyme_names));

        return RestrictionAnalyzer.GetDigestSummary(new DnaSequence(sequence), enzyme_names);
    }

    [McpServerTool(Name = "restriction_map", Title = "MolTools — Restriction Map", ReadOnly = true), Description("Builds a restriction map of a DNA sequence: all forward+reverse sites, sites grouped by enzyme, the total forward-strand site count, unique cutters (enzymes with exactly one forward-strand site), and non-cutters from the queried enzyme set. An empty enzyme_names list considers every built-in enzyme. Call to plan cloning around single-cutter enzymes.")]
    public static RestrictionMap restriction_map(
        [Description("DNA sequence to map.")] string sequence,
        [Description("Names of restriction enzymes to consider; empty to use every built-in enzyme.")] string[] enzyme_names)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        return RestrictionAnalyzer.CreateMap(new DnaSequence(sequence), enzyme_names ?? System.Array.Empty<string>());
    }

    [McpServerTool(Name = "compatible_enzymes", Title = "MolTools — Compatible Enzyme Pairs", ReadOnly = true), Description("Enumerates all pairs of built-in restriction enzymes whose ends can be ligated to each other — either both produce blunt ends, or both produce the same overhang type and overhang sequence. Call when the user needs enzyme pairs that yield compatible (ligatable) ends for cloning.")]
    public static CompatibleEnzymesResult compatible_enzymes()
    {
        var pairs = RestrictionAnalyzer.FindCompatibleEnzymes()
            .Select(t => new EnzymeCompatibilityPair(t.Enzyme1, t.Enzyme2, t.CompatibleEnd))
            .ToList();
        return new CompatibleEnzymesResult(pairs);
    }

    [McpServerTool(Name = "enzymes_compatible", Title = "MolTools — Enzyme End Compatibility", ReadOnly = true), Description("Returns whether two named built-in restriction enzymes produce ligatable (compatible) ends — both blunt, or the same overhang type and sequence. An unknown enzyme name yields false (no error). Call to check a specific pair for cloning compatibility.")]
    public static EnzymeCompatibilityResult enzymes_compatible(
        [Description("First enzyme name.")] string enzyme1_name,
        [Description("Second enzyme name.")] string enzyme2_name)
    {
        if (string.IsNullOrWhiteSpace(enzyme1_name))
            throw new System.ArgumentException("First enzyme name cannot be null or blank.", nameof(enzyme1_name));
        if (string.IsNullOrWhiteSpace(enzyme2_name))
            throw new System.ArgumentException("Second enzyme name cannot be null or blank.", nameof(enzyme2_name));

        return new EnzymeCompatibilityResult(RestrictionAnalyzer.AreCompatible(enzyme1_name, enzyme2_name));
    }

    #endregion

    #region CodonUsageAnalyzer

    [McpServerTool(Name = "count_codons", Title = "MolTools — Count Codons", ReadOnly = true), Description("Counts occurrences of each ACGT codon in a coding DNA sequence (frame 0, non-overlapping triplets). Trailing partial codons and codons containing non-ACGT characters are skipped silently. Call to get the raw codon-frequency table of a gene.")]
    public static CodonCountsResult count_codons(
        [Description("Coding DNA sequence (frame 0).")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        return new CodonCountsResult(CodonUsageAnalyzer.CountCodons(sequence));
    }

    [McpServerTool(Name = "rscu", Title = "MolTools — Relative Synonymous Codon Usage", ReadOnly = true), Description("Relative Synonymous Codon Usage (RSCU) per codon: observed count / count-expected-if-uniform among its synonymous codons. RSCU = 1 means no bias, > 1 over-represented, < 1 under-represented. Call to quantify codon bias per codon in a coding sequence.")]
    public static RscuResult rscu(
        [Description("Coding DNA sequence (frame 0).")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        return new RscuResult(CodonUsageAnalyzer.CalculateRscu(sequence));
    }

    [McpServerTool(Name = "codon_adaptation_index", Title = "MolTools — Codon Adaptation Index (RSCU)", ReadOnly = true), Description("Codon Adaptation Index (Sharp & Li 1987) using a caller-supplied reference RSCU table (typically derived from highly expressed genes), codons in the DNA alphabet. Output range 0..1; per-codon relative adaptiveness w = RSCU/max-synonymous-RSCU, and CAI is their geometric mean. Single-codon amino acids (Met/Trp), stop codons, and codons with w=0 are excluded. Call to score how well a gene matches an organism's preferred codons given an RSCU reference (distinct from cai_from_organism_table, which takes a frequency table).")]
    public static CaiResult codon_adaptation_index(
        [Description("Coding DNA sequence (frame 0).")] string sequence,
        [Description("Reference RSCU table: codon (DNA alphabet) → RSCU value.")] Dictionary<string, double> reference_rscu)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (reference_rscu is null)
            throw new System.ArgumentException("Reference RSCU table cannot be null.", nameof(reference_rscu));

        return new CaiResult(CodonUsageAnalyzer.CalculateCai(sequence, reference_rscu));
    }

    [McpServerTool(Name = "effective_number_of_codons", Title = "MolTools — Effective Number of Codons (ENC)", ReadOnly = true), Description("Effective Number of Codons (Wright's Nc), measuring how far a gene departs from uniform synonymous-codon usage. Result is clamped to 20..61 (20 = extreme bias, 61 = no bias). Call to summarise a gene's overall codon bias with a single number.")]
    public static EncResult effective_number_of_codons(
        [Description("Coding DNA sequence (frame 0).")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        return new EncResult(CodonUsageAnalyzer.CalculateEnc(sequence));
    }

    [McpServerTool(Name = "codon_usage_statistics", Title = "MolTools — Codon-Usage Statistics", ReadOnly = true), Description("Aggregate codon-usage report for a coding sequence: per-codon counts, RSCU, Effective Number of Codons (ENC), total codons, GC% at codon positions 1/2/3, GC3s (synonymous third-position GC), and overall GC. Call for a one-shot codon-usage summary of a gene.")]
    public static CodonUsageStatistics codon_usage_statistics(
        [Description("Coding DNA sequence (frame 0).")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

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

    [McpServerTool(Name = "cai_from_organism_table", Title = "MolTools — CAI from Codon-Usage Table", ReadOnly = true), Description("Computes the Codon Adaptation Index (Sharp & Li 1987) for a coding sequence against an organism codon-usage FREQUENCY table (distinct from codon_adaptation_index, which takes a reference RSCU dictionary). CAI is the geometric mean of per-codon relative adaptiveness w = f(codon)/max f(synonymous). Codons without synonymous-group data are skipped; w is clamped at 1e-6 to avoid ln(0) on incomplete custom tables. Call when scoring how well a gene matches an organism's preferred codons.")]
    public static CaiResult cai_from_organism_table(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Target organism: preset id (EColiK12 | Yeast | Human) or inline custom table.")] CodonUsageTableInput target_organism)
    {
        if (string.IsNullOrEmpty(coding_sequence))
            throw new System.ArgumentException("Coding sequence cannot be null or empty.", nameof(coding_sequence));

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

    [McpServerTool(Name = "compare_codon_usage", Title = "MolTools — Compare Codon Usage", ReadOnly = true), Description("Codon-frequency similarity between two coding sequences: 1 − ½·Σ|f1−f2| ∈ [0,1] (1 = identical codon distribution, 0 = disjoint). An input that is empty or contains no complete codons contributes 0 similarity. Call to compare the codon usage of two genes/organisms.")]
    public static SimilarityResult compare_codon_usage(
        [Description("First coding sequence.")] string sequence1,
        [Description("Second coding sequence.")] string sequence2)
    {
        if (sequence1 is null)
            throw new System.ArgumentException("First sequence cannot be null.", nameof(sequence1));
        if (sequence2 is null)
            throw new System.ArgumentException("Second sequence cannot be null.", nameof(sequence2));

        return new SimilarityResult(CodonOptimizer.CompareCodonUsage(sequence1, sequence2));
    }

    [McpServerTool(Name = "build_codon_table", Title = "MolTools — Build Codon-Usage Table", ReadOnly = true), Description("Derives a per-organism CodonUsageTable from a reference coding sequence by computing per-amino-acid relative codon frequencies (RNA alphabet, U not T). Call when the user wants a custom codon-usage table built from their own reference gene(s).")]
    public static CodonUsageTableDto build_codon_table(
        [Description("Reference coding sequence (DNA or RNA).")] string reference_sequence,
        [Description("Organism name to attach to the resulting table.")] string organism_name)
    {
        if (string.IsNullOrEmpty(reference_sequence))
            throw new System.ArgumentException("Reference sequence cannot be null or empty.", nameof(reference_sequence));
        if (string.IsNullOrWhiteSpace(organism_name))
            throw new System.ArgumentException("Organism name cannot be null or blank.", nameof(organism_name));

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

    [McpServerTool(Name = "crispr_system_info", Title = "MolTools — CRISPR System Metadata", ReadOnly = true), Description("Returns the metadata record (name, PAM sequence, guide length, PAM placement relative to target, description) for a known CRISPR nuclease system. Call when the user needs the PAM/guide parameters of SpCas9, SaCas9, Cas12a, CasX, etc.")]
    public static CrisprSystem crispr_system_info(
        [Description("CRISPR system: SpCas9 | SpCas9_NAG | SaCas9 | Cas12a | AsCas12a | LbCas12a | CasX.")] CrisprSystemType system_type)
    {
        return CrisprDesigner.GetSystem(system_type);
    }

    [McpServerTool(Name = "find_pam_sites", Title = "MolTools — Find CRISPR PAM Sites", ReadOnly = true), Description("Finds all PAM matches (forward + reverse strand) for the chosen CRISPR system. PAM matching honours IUPAC codes (e.g. NGG, NNGRRT, TTTV). Each site reports the PAM, the adjacent guide/target window, its position and strand; sites whose target window falls outside the sequence are skipped. Call to enumerate targetable protospacers in a sequence.")]
    public static PamSitesResult find_pam_sites(
        [Description("DNA sequence to scan.")] string sequence,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

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

    [McpServerTool(Name = "crispr_specificity_score", Title = "MolTools — CRISPR Guide Specificity Score", ReadOnly = true), Description("Aggregates off-target hits (≤4 mismatches) for a guide RNA against a genome into a single specificity score in 0..100 (100 = no off-targets; the score drops as more/seed-region off-targets are found). Call to judge how genome-specific a candidate guide is. Guide length must match the system's guide length.")]
    public static SpecificityResult crispr_specificity_score(
        [Description("Guide RNA sequence (length must match the system's guide length).")] string guide_sequence,
        [Description("Genome / reference sequence to scan.")] string genome,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9)
    {
        if (string.IsNullOrEmpty(guide_sequence))
            throw new System.ArgumentException("Guide sequence cannot be null or empty.", nameof(guide_sequence));
        if (string.IsNullOrEmpty(genome))
            throw new System.ArgumentException("Genome cannot be null or empty.", nameof(genome));

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

    [McpServerTool(Name = "oligo_extinction_coefficient", Title = "MolTools — Oligo Extinction Coefficient", ReadOnly = true), Description("Sums per-base 260 nm molar extinction contributions (A=15400, C=7400, G=11500, T=8700, U=9900 M⁻¹·cm⁻¹; any other base = 10000) for an oligonucleotide. Call to estimate an oligo's ε₂₆₀ for concentration calculations.")]
    public static ExtinctionCoefficientResult oligo_extinction_coefficient(
        [Description("Oligonucleotide sequence.")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        return new ExtinctionCoefficientResult(ProbeDesigner.CalculateExtinctionCoefficient(sequence));
    }

    [McpServerTool(Name = "oligo_concentration_from_absorbance", Title = "MolTools — Oligo Concentration (Beer–Lambert)", ReadOnly = true), Description("Computes oligonucleotide concentration in µM from the Beer–Lambert law: c = A₂₆₀ / (ε · path) · 1e6. Call to convert a spectrophotometer A260 reading into a molar concentration given the oligo's extinction coefficient.")]
    public static ConcentrationResult oligo_concentration_from_absorbance(
        [Description("Absorbance at 260 nm (A260).")] double absorbance260,
        [Description("Extinction coefficient ε in M⁻¹·cm⁻¹.")] double extinction_coefficient,
        [Description("Path length in cm (default 1.0).")] double path_length = 1.0)
    {
        if (extinction_coefficient <= 0)
            throw new System.ArgumentException("Extinction coefficient must be positive.", nameof(extinction_coefficient));
        if (path_length <= 0)
            throw new System.ArgumentException("Path length must be positive.", nameof(path_length));

        return new ConcentrationResult(
            ProbeDesigner.CalculateConcentration(absorbance260, extinction_coefficient, path_length));
    }

    #endregion
}
