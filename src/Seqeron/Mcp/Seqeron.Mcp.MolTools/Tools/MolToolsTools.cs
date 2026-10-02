using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;

namespace Seqeron.Mcp.MolTools.Tools;

[McpServerToolType]
public class MolToolsTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private MolToolsTools() { }

    #region PrimerDesigner

    [McpServerTool(Name = "design_primers", Title = "MolTools — Design PCR Primer Pair", ReadOnly = true), Description("Designs forward/reverse PCR primers flanking a target region with Primer3's pair search (verified against primer3-py design_primers): candidates on either side of the target (never overlapping it) are kept when they pass the per-primer limits (length, GC%, Primer3-default SantaLucia Tm at 50 mM Na+/1.5 mM Mg2+/0.6 mM dNTP/50 nM, poly-X, dinucleotide repeat) and, by default, Primer3's thermodynamic secondary-structure screen (ntthal self-dimer, 3' self-dimer and hairpin Tm <= 47 °C per primer); pairs must have a product size in product_size_range (Primer3 PRIMER_PRODUCT_SIZE_RANGE, default 100-300 bp, ranges tried in order), |Tm_f - Tm_r| <= max_tm_difference (default 5 °C) and pair hetero-dimer / 3' hetero-dimer ntthal Tm <= 47 °C; the pair with the lowest Primer3 pair penalty (sum of per-primer penalties) is returned, with product Tm (Primer3 long_seq_tm), pair complementarity Tm values, optionally an internal hybridization oligo (pick_internal_oligo, Primer3 PRIMER_PICK_INTERNAL_OLIGO) and up to num_return ranked pairs (PRIMER_NUM_RETURN). The target is the half-open 0-based interval [target_start, target_end) with 0 <= target_start < target_end < template.Length.")]
    public static DesignPrimersResult design_primers(
        [Description("DNA template (A/C/G/T).")] string template,
        [Description("0-based inclusive start of target region.")] int target_start,
        [Description("0-based exclusive end of target region (primers never overlap [target_start, target_end)).")] int target_end,
        [Description("Optional primer design parameters (lengths, GC%, Tm, repeats, GC-clamp/3' stability checks, structure screen). Defaults are used if null.")] PrimerParameters? parameters = null,
        [Description("PRIMER_PRODUCT_SIZE_RANGE in Primer3 syntax, e.g. \"100-300\" or \"150-250 100-400\" (ranges in order of preference; default 100-300).")] string? product_size_range = null,
        [Description("PRIMER_PAIR_MAX_DIFF_TM: maximum |Tm_forward - Tm_reverse| in °C (default 5; Primer3's own default is 100).")] double max_tm_difference = PrimerDesigner.MaxPairTmDifference,
        [Description("PRIMER_NUM_RETURN: number of ranked pairs listed in 'pairs' (default 1).")] int num_return = 1,
        [Description("PRIMER_PICK_INTERNAL_OLIGO: also pick an internal hybridization oligo between the primers (Primer3 PRIMER_INTERNAL_* defaults; default false).")] bool pick_internal_oligo = false)
    {
        if (string.IsNullOrEmpty(template))
            throw new System.ArgumentException("Template cannot be null or empty.", nameof(template));
        if (target_start < 0)
            throw new System.ArgumentException("Target start must be non-negative.", nameof(target_start));
        if (target_end >= template.Length)
            throw new System.ArgumentException("Target end must be within the template.", nameof(target_end));
        if (target_start >= target_end)
            throw new System.ArgumentException("Target start must be strictly less than target end.", nameof(target_start));
        if (num_return < 1)
            throw new System.ArgumentException("num_return must be at least 1.", nameof(num_return));
        if (!(max_tm_difference >= 0))
            throw new System.ArgumentException("max_tm_difference must be non-negative.", nameof(max_tm_difference));

        var options = PrimerPairOptions.Default with
        {
            MaxTmDifference = max_tm_difference,
            NumReturn = num_return,
            PickInternalOligo = pick_internal_oligo,
        };
        if (product_size_range is not null)
            options = options with { ProductSizeRanges = ParseProductSizeRanges(product_size_range) };

        var dna = new DnaSequence(template);
        var pairs = PrimerDesigner.DesignPrimerPairs(dna, target_start, target_end, parameters, options);
        var best = pairs.Count > 0 ? pairs[0] : PrimerDesigner.DesignPrimers(dna, target_start, target_end, parameters, options);
        return new DesignPrimersResult(
            best.Forward, best.Reverse, best.IsValid, best.Message, best.ProductSize,
            best.PairPenalty, best.ProductTm, best.ComplAnyTh, best.ComplEndTh, best.InternalOligo, pairs);
    }

    // Primer3 PRIMER_PRODUCT_SIZE_RANGE syntax: space-separated "min-max" ranges.
    private static List<ProductSizeRange> ParseProductSizeRanges(string text)
    {
        var ranges = new List<ProductSizeRange>();
        foreach (var token in text.Split(new[] { ' ', ',', ';' }, System.StringSplitOptions.RemoveEmptyEntries))
        {
            var parts = token.Split('-');
            if (parts.Length != 2
                || !int.TryParse(parts[0], System.Globalization.NumberStyles.None, System.Globalization.CultureInfo.InvariantCulture, out int min)
                || !int.TryParse(parts[1], System.Globalization.NumberStyles.None, System.Globalization.CultureInfo.InvariantCulture, out int max)
                || min < 1 || max < min)
                throw new System.ArgumentException($"Invalid product size range '{token}' (expected min-max with 1 <= min <= max).", nameof(text));
            ranges.Add(new ProductSizeRange(min, max));
        }
        if (ranges.Count == 0)
            throw new System.ArgumentException("product_size_range must contain at least one min-max range.", nameof(text));
        return ranges;
    }

    [McpServerTool(Name = "evaluate_primer", Title = "MolTools — Evaluate Primer", ReadOnly = true), Description("Evaluates a single primer sequence against quality criteria and returns a scored candidate: length, GC%, Tm (Primer3-default SantaLucia 1998 NN Tm), longest homopolymer, the Primer3 thermodynamic secondary-structure Tm values (hairpinTh / selfAnyTh / selfEndTh = primer3 calc_hairpin / calc_homodimer / calc_end_stability Tm at 50 mM Na+, 1.5 mM Mg2+, 0.6 mM dNTP, 50 nM; hasHairpin = hairpinTh > 47 °C, PRIMER_MAX_HAIRPIN_TH), 3'-end stability, an issues list, validity flag, an informational numeric score and the Primer3 per-primer penalty. Call to QC one primer (position/strand are informational).")]
    public static PrimerCandidate evaluate_primer(
        [Description("Primer sequence to evaluate.")] string sequence,
        [Description("0-based location of the primer in the template (informational).")] int position,
        [Description("True if this is a forward primer; false for reverse.")] bool is_forward,
        [Description("Optional primer design parameters.")] PrimerParameters? parameters = null)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        return PrimerDesigner.EvaluatePrimer(sequence, position, is_forward, parameters);
    }

    [McpServerTool(Name = "primer_melting_temperature", Title = "MolTools — Primer Melting Temperature", ReadOnly = true), Description("Computes a primer's melting temperature (Tm, °C): Wallace rule Tm = 2·(A+T) + 4·(G+C) for < 14 valid bases, or Marmur–Doty Tm = 64.9 + 41·(GC−16.4)/N for ≥ 14 valid bases. A/C/G/T/U are counted (U read as T, as Biopython); other characters are ignored. Call for a quick Tm estimate of a short oligo/primer.")]
    public static TmResult primer_melting_temperature(
        [Description("Primer sequence.")] string primer)
    {
        if (string.IsNullOrEmpty(primer))
            throw new System.ArgumentException("Primer cannot be null or empty.", nameof(primer));

        return new TmResult(PrimerDesigner.CalculateMeltingTemperature(primer));
    }

    [McpServerTool(Name = "primer_melting_temperature_salt", Title = "MolTools — Salt-Corrected Primer Tm", ReadOnly = true), Description("Salt-adjusted primer Tm (OligoCalc, Kibbe 2007): < 14 valid bases Tm = 2·(A+T) + 4·(G+C) + 16.6·log10([Na+]/0.050 M); ≥ 14 valid bases Tm = 100.5 + 41·(G+C)/N − 820/N + 16.6·log10([Na+] M); rounded to one decimal. Call when a monovalent-cation ([Na+]) adjusted primer Tm is needed. Na+ concentration is in mM (default 50).")]
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

    [McpServerTool(Name = "longest_homopolymer", Title = "MolTools — Longest Homopolymer Run", ReadOnly = true), Description("Returns the length of the longest run of identical consecutive nucleotides (e.g. AAAA = 4) in a sequence, case-insensitive — the quantity Primer3 limits with PRIMER_MAX_POLY_X; as in Primer3, N counts as the worst-case base (ANA = 3). Call to flag homopolymer stretches that hurt primer/probe quality.")]
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

    [McpServerTool(Name = "primer_dimer", Title = "MolTools — Primer-Dimer Check", ReadOnly = true), Description("Primer3 alignment-mode 3'-end primer-dimer check (PRIMER_PAIR_COMPL_END with PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0): the 3'-anchored dpal alignment score of each primer against the other's reverse complement (+1 per complementary pair, -1 mismatch, -2 per single-base gap; max of both orientations). Flags a dimer when the score is at least min_complementarity (default 4 = Primer3's default PRIMER_PAIR_MAX_COMPL_END 3.00 exceeded). Returns the flag, the integer score and the exact score. Call to screen a primer pair for 3'-dimer formation; for Primer3's default thermodynamic check use the C# API PrimerDesigner.CalculatePrimer3PairComplementarity.")]
    public static PrimerDimerResult primer_dimer(
        [Description("First primer sequence (5'->3').")] string primer1,
        [Description("Second primer sequence (5'->3').")] string primer2,
        [Description("Minimum 3'-end complementarity score to flag a dimer (default 4).")] int min_complementarity = 4)
    {
        if (string.IsNullOrEmpty(primer1))
            throw new System.ArgumentException("First primer cannot be null or empty.", nameof(primer1));
        if (string.IsNullOrEmpty(primer2))
            throw new System.ArgumentException("Second primer cannot be null or empty.", nameof(primer2));

        double score = PrimerDesigner.CalculatePrimerDimerEndComplementarity(primer1, primer2);
        return new PrimerDimerResult(
            PrimerDesigner.HasPrimerDimer(primer1, primer2, min_complementarity),
            (int)System.Math.Floor(score),
            score);
    }

    [McpServerTool(Name = "three_prime_stability", Title = "MolTools — Primer 3' End Stability (ΔG°37)", ReadOnly = true), Description("Primer3 3'-end stability (oligotm.c end_oligodg): SantaLucia (1998) nearest-neighbor ΔG°37 (kcal/mol, 1 M NaCl) of the primer's last 5 bases (the whole primer if shorter), with initiation (+1.96, +0.05 per terminal A·T, +0.43 if self-complementary). Primer3 reports the same magnitude with the opposite sign as PRIMER_*_END_STABILITY. More negative ΔG = a more stable (more problematic) 3' end. N is accepted (Primer3 N parameters); other characters in the 3' window are rejected. Call to assess primer 3'-end stability for mispriming risk.")]
    public static ThreePrimeStabilityResult three_prime_stability(
        [Description("Primer sequence.")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        double dg = PrimerDesigner.Calculate3PrimeStability(sequence);
        if (double.IsNaN(dg))
            throw new System.ArgumentException("The 3'-terminal 5 bases may contain only A, C, G, T or N.", nameof(sequence));
        return new ThreePrimeStabilityResult(dg);
    }

    [McpServerTool(Name = "generate_primer_candidates", Title = "MolTools — Generate Primer Candidates", ReadOnly = true), Description("Enumerates all primer candidates of admissible lengths (parameters.MinLength..MaxLength) at every start position within a region of the template and evaluates each one (candidates are emitted in generation order, NOT sorted by score). region_start is 0-based inclusive, region_end is exclusive; for a reverse request each candidate sequence is the reverse complement of the template substring. Useful when the caller wants the full candidate set to pick by custom criteria.")]
    public static PrimerCandidateListResult generate_primer_candidates(
        [Description("DNA template sequence.")] string template,
        [Description("0-based start of the search region (inclusive).")] int region_start,
        [Description("0-based end of the search region (exclusive).")] int region_end,
        [Description("True for forward-strand candidates, false for reverse-complement candidates.")] bool forward = true,
        [Description("Optional primer design parameters.")] PrimerParameters? parameters = null)
    {
        if (string.IsNullOrEmpty(template))
            throw new System.ArgumentException("Template cannot be null or empty.", nameof(template));
        if (region_start < 0)
            throw new System.ArgumentException("Region start must be non-negative.", nameof(region_start));
        if (region_end > template.Length)
            throw new System.ArgumentException("Region end must not exceed the template length.", nameof(region_end));
        if (region_start >= region_end)
            throw new System.ArgumentException("Region start must be strictly less than region end.", nameof(region_start));

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

    [McpServerTool(Name = "codon_adaptation_index", Title = "MolTools — Codon Adaptation Index (RSCU)", ReadOnly = true), Description("Codon Adaptation Index (Sharp & Li 1987) using a caller-supplied reference RSCU table (typically derived from highly expressed genes), codons in the DNA alphabet. Output range 0..1; per-codon relative adaptiveness w = RSCU/max-synonymous-RSCU, and CAI is their geometric mean. Single-codon amino acids (Met/Trp) and stop codons are excluded; a codon with w < 0.0001 (absent from the reference) is scored with w = 0.01 (CodonW convention); non-ACGT(U) triplets are skipped. Call to score how well a gene matches an organism's preferred codons given an RSCU reference (distinct from cai_from_organism_table, which takes a frequency table).")]
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

    [McpServerTool(Name = "effective_number_of_codons", Title = "MolTools — Effective Number of Codons (ENC)", ReadOnly = true), Description("Effective Number of Codons (Wright 1990 Nc; CodonW enc_out conventions), measuring how far a gene departs from uniform synonymous-codon usage under the standard genetic code. Range 20..61 (20 = extreme bias, 61 = no bias; overshoot re-adjusted to 61). Amino acids seen once, or with every observed codon used once, are not estimable; a missing isoleucine class is replaced by the mean of the 2- and 4-fold classes. Returns 0 when Nc cannot be calculated (some other synonymous class has no estimable amino acid - gene too short or amino-acid usage too skewed; CodonW prints *****). DNA or RNA, case-insensitive; non-ACGT(U) triplets are skipped. Call to summarise a gene's overall codon bias with a single number.")]
    public static EncResult effective_number_of_codons(
        [Description("Coding DNA sequence (frame 0).")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        return new EncResult(CodonUsageAnalyzer.CalculateEnc(sequence));
    }

    [McpServerTool(Name = "codon_usage_statistics", Title = "MolTools — Codon-Usage Statistics", ReadOnly = true), Description("Aggregate codon-usage report for a coding sequence: per-codon counts, RSCU, Effective Number of Codons (ENC; 0 when not calculable, see effective_number_of_codons), total codons, GC% at codon positions 1/2/3, GC3s (synonymous third-position GC), and overall GC. Call for a one-shot codon-usage summary of a gene.")]
    public static CodonUsageStatistics codon_usage_statistics(
        [Description("Coding DNA sequence (frame 0).")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        return CodonUsageAnalyzer.GetStatistics(sequence);
    }

    #endregion

    #region CodonOptimizer

    [McpServerTool(Name = "optimize_codons", Title = "MolTools — Optimize Codons for Expression", ReadOnly = true), Description("Optimizes a coding sequence for expression in a target organism using one of five strategies (MaximizeCAI, BalancedOptimization (default), HarmonizeExpression, MinimizeSecondary, AvoidRareCodeons). Internally trims to whole codons and converts T→U; stop codons and single-codon amino acids (Met/Trp) are left unchanged. Returns the original/optimized RNA, translated protein, original/optimized CAI, GC fractions, the number of changed codons, and each codon change. HarmonizeExpression matches the target codon-usage profile as closely as integer codon counts allow (DNA Chisel match_codon_usage) and is deterministic.")]
    public static OptimizationResultDto optimize_codons(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Target organism: a preset id (EColiK12 | Yeast | Human) or an inline custom table (organismName + codonFrequencies in RNA alphabet).")] CodonUsageTableInput target_organism,
        [Description("Optimization strategy.")] CodonOptimizer.OptimizationStrategy strategy = CodonOptimizer.OptimizationStrategy.BalancedOptimization,
        [Description("Lower bound for target GC fraction (default 0.40).")] double gc_target_min = 0.40,
        [Description("Upper bound for target GC fraction (default 0.60).")] double gc_target_max = 0.60,
        [Description("Frequency threshold below which a codon is considered rare (default 0.15).")] double rare_codon_threshold = 0.15)
    {
        if (string.IsNullOrEmpty(coding_sequence))
            throw new System.ArgumentException("Coding sequence cannot be null or empty.", nameof(coding_sequence));
        if (gc_target_min < 0 || gc_target_min > 1)
            throw new System.ArgumentException("GC target minimum must be a fraction in [0, 1].", nameof(gc_target_min));
        if (gc_target_max < 0 || gc_target_max > 1)
            throw new System.ArgumentException("GC target maximum must be a fraction in [0, 1].", nameof(gc_target_max));
        if (gc_target_min > gc_target_max)
            throw new System.ArgumentException("GC target minimum must not exceed the maximum.", nameof(gc_target_min));

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

    [McpServerTool(Name = "cai_from_organism_table", Title = "MolTools — CAI from Codon-Usage Table", ReadOnly = true), Description("Computes the Codon Adaptation Index (Sharp & Li 1987) for a coding sequence against an organism codon-usage FREQUENCY table (distinct from codon_adaptation_index, which takes a reference RSCU dictionary). CAI is the geometric mean of per-codon relative adaptiveness w = f(codon)/max f(synonymous). Stop codons and single-codon amino acids (Met/Trp) are excluded (Sharp & Li; CodonW); a codon with w < 0.0001 (absent from the table while a synonym is present) is scored with w = 0.01 (CodonW); amino acids without frequency data are skipped. Same core as codon_adaptation_index. Call when scoring how well a gene matches an organism's preferred codons.")]
    public static CaiResult cai_from_organism_table(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Target organism: preset id (EColiK12 | Yeast | Human) or inline custom table.")] CodonUsageTableInput target_organism)
    {
        if (string.IsNullOrEmpty(coding_sequence))
            throw new System.ArgumentException("Coding sequence cannot be null or empty.", nameof(coding_sequence));

        var table = ResolveCodonUsageTable(target_organism);
        return new CaiResult(CodonOptimizer.CalculateCAI(coding_sequence, table));
    }

    [McpServerTool(Name = "remove_restriction_sites", Title = "MolTools — Remove Restriction Sites", ReadOnly = true), Description("Synonymously rewrites codons to eliminate the listed restriction recognition sequences from a coding sequence while preserving the encoded protein (RNA-alphabet output). Site strings may be DNA or RNA and may contain IUPAC ambiguity codes; both strands are cleared (a non-palindromic site such as BsaI GGTCTC is also removed where its reverse complement occurs); each removal changes exactly one codon, choosing the highest-frequency synonymous codon in the supplied table; sites with no synonymous alternative are left in place. Call to make a gene compatible with a cloning strategy.")]
    public static OptimizedSequenceResult remove_restriction_sites(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Restriction recognition sequences to eliminate.")] string[] restriction_sites,
        [Description("Target organism: preset id or inline custom table (used for synonymous-codon lookup).")] CodonUsageTableInput target_organism)
    {
        if (string.IsNullOrEmpty(coding_sequence))
            throw new System.ArgumentException("Coding sequence cannot be null or empty.", nameof(coding_sequence));
        if (restriction_sites is null || restriction_sites.Length == 0)
            throw new System.ArgumentException("At least one restriction site is required.", nameof(restriction_sites));

        var table = ResolveCodonUsageTable(target_organism);
        return new OptimizedSequenceResult(CodonOptimizer.RemoveRestrictionSites(coding_sequence, restriction_sites, table));
    }

    [McpServerTool(Name = "reduce_secondary_structure", Title = "MolTools — Reduce mRNA Secondary Structure", ReadOnly = true), Description("Greedy synonymous-codon swap that lowers a heuristic local self-complementarity score (canonical pairs incl. G\u00b7U wobble) within a sliding window, reducing mRNA secondary structure while preserving the protein. Output is upper-case RNA trimmed to whole codons; sequences shorter than window_size are returned unchanged apart from that normalisation. Heuristic only \u2014 not a thermodynamic folding model (see rna_minimum_free_energy for that). Call to relax strong secondary structure in a coding sequence.")]
    public static OptimizedSequenceResult reduce_secondary_structure(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Target organism: preset id or inline custom table.")] CodonUsageTableInput target_organism,
        [Description("Sliding-window size in nucleotides (default 40).")] int window_size = 40)
    {
        if (string.IsNullOrEmpty(coding_sequence))
            throw new System.ArgumentException("Coding sequence cannot be null or empty.", nameof(coding_sequence));
        if (window_size <= 0)
            throw new System.ArgumentException("Window size must be positive.", nameof(window_size));

        var table = ResolveCodonUsageTable(target_organism);
        return new OptimizedSequenceResult(CodonOptimizer.ReduceSecondaryStructure(coding_sequence, table, window_size));
    }

    [McpServerTool(Name = "find_rare_codons", Title = "MolTools — Find Rare Codons", ReadOnly = true), Description("Reports every codon in a coding sequence whose frequency in the target organism's codon-usage table is below the threshold (default 0.15), with its 0-based position, codon (RNA), amino acid (Standard code, * for stop) and frequency. Frame-0 complete triplets only; ambiguous triplets (N, R, Y, …) are skipped without shifting the frame. Call to locate translation-slowing rare codons before optimization.")]
    public static RareCodonsResult find_rare_codons(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Target organism: preset id or inline custom table.")] CodonUsageTableInput target_organism,
        [Description("Frequency threshold (default 0.15).")] double threshold = 0.15)
    {
        if (string.IsNullOrEmpty(coding_sequence))
            throw new System.ArgumentException("Coding sequence cannot be null or empty.", nameof(coding_sequence));

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

    [McpServerTool(Name = "build_codon_table", Title = "MolTools — Build Codon-Usage Table", ReadOnly = true), Description("Derives a per-organism CodonUsageTable from a reference coding sequence by computing per-amino-acid relative codon frequencies (RNA alphabet, U not T). All 64 codons are returned: a codon absent from the reference set is counted as 0.5 (Sharp & Li 1987, as in Biopython CodonAdaptationIndex), so the relative adaptiveness derived from the table matches Biopython exactly. Call when the user wants a custom codon-usage table built from their own reference gene(s).")]
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
    /// table) to a <see cref="CodonOptimizer.CodonUsageTable"/>. Inline tables are built with
    /// <see cref="CodonOptimizer.CreateCodonUsageTable"/>, so they carry the canonical Standard
    /// genetic-code mapping in <c>CodonToAminoAcid</c> just like the presets.
    /// </summary>
    private static CodonOptimizer.CodonUsageTable ResolveCodonUsageTable(CodonUsageTableInput input)
    {
        ArgumentNullException.ThrowIfNull(input);

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

        return CodonOptimizer.CreateCodonUsageTable(
            input.OrganismName ?? "Custom",
            input.CodonFrequencies);
    }

    #endregion

    #region CrisprDesigner

    [McpServerTool(Name = "crispr_system_info", Title = "MolTools — CRISPR System Metadata", ReadOnly = true), Description("Returns the metadata record (name, PAM sequence, guide length, PAM placement relative to target, description) for a known CRISPR nuclease system. Call when the user needs the PAM/guide parameters of SpCas9, SaCas9, Cas12a, CasX, etc.")]
    public static CrisprSystem crispr_system_info(
        [Description("CRISPR system: SpCas9 | SpCas9_NAG | SaCas9 | Cas12a | AsCas12a | LbCas12a | CasX.")] CrisprSystemType system_type)
    {
        return CrisprDesigner.GetSystem(system_type);
    }

    [McpServerTool(Name = "find_pam_sites", Title = "MolTools — Find CRISPR PAM Sites", ReadOnly = true), Description("Finds all PAM matches (forward + reverse strand) for the chosen CRISPR system. PAM matching honours IUPAC codes (e.g. NGG, NNGRRT, TTTV). Each site reports the PAM, the adjacent guide/target window, its position and strand; sites whose target window falls outside the sequence are skipped. Coordinates (position, targetStart) are always 0-based forward-strand; the PAM and guide sequences are read on the protospacer strand (CRISPOR convention). Call to enumerate targetable protospacers in a sequence.")]
    public static PamSitesResult find_pam_sites(
        [Description("DNA sequence to scan.")] string sequence,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        var sites = CrisprDesigner.FindPamSites(sequence, system_type).ToList();
        return new PamSitesResult(sites);
    }

    [McpServerTool(Name = "design_guide_rnas", Title = "MolTools — Design CRISPR Guide RNAs", ReadOnly = true), Description("Generates and scores guide-RNA candidates whose Cas9/Cas12a cut site falls inside the requested region, ranked best-first. The cut site follows the CRISPOR convention: 3 bp 5\u0027 of the PAM on the PAM-bearing strand for Cas9, after the 18th protospacer base for Cas12a - on both strands. Candidates scoring below parameters.MinScore are filtered out. 20-nt NGG guides with 4 nt of 5\u0027 and 3 nt of 3\u0027 flanking context additionally report context30Mer and onTargetScore (the published Doench 2016 Rule Set 2 / Azimuth on-target efficacy score, 0..1), and parameters.ranking = OnTargetRuleSet2 ranks by it; grafMotif flags the Graf 2019 TT-/GCC- inefficiency motifs. Region indices are 0-based; region_end is inclusive and must satisfy 0 <= region_start <= region_end < sequence.Length. Call to enumerate high-quality guides targeting a locus.")]
    public static GuideRnasResult design_guide_rnas(
        [Description("DNA sequence containing the target region.")] string sequence,
        [Description("0-based start of the target region.")] int region_start,
        [Description("0-based inclusive end of the target region.")] int region_end,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9,
        [Description("Optional guide-RNA design parameters (minGcContent, maxGcContent, minScore, avoidPolyT, checkSelfComplementarity, ranking). Defaults are used when null.")] GuideRnaParameters? parameters = null)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (region_start < 0 || region_start >= sequence.Length)
            throw new System.ArgumentException("Region start must be within the sequence.", nameof(region_start));
        if (region_end < region_start || region_end >= sequence.Length)
            throw new System.ArgumentException("Region end must satisfy region_start <= region_end < sequence.Length.", nameof(region_end));

        var guides = CrisprDesigner
            .DesignGuideRnas(new DnaSequence(sequence), region_start, region_end, system_type, parameters)
            .ToList();
        return new GuideRnasResult(guides);
    }

    [McpServerTool(Name = "evaluate_guide_rna", Title = "MolTools — Evaluate Guide RNA", ReadOnly = true), Description("Scores a single guide RNA against on-target quality heuristics: overall GC%, seed-region GC%, polyT (Pol III terminator) presence, self-complementarity, and common-restriction-site presence; returns a 0..100 score plus an issues list. avoidPolyT / checkSelfComplementarity switch off the corresponding penalties. For NGG systems the Graf 2019 TT-/GCC- inefficiency motif is reported in grafMotif and in issues (a warning only, no deduction, as in CRISPOR). The published Doench 2016 Rule Set 2 score needs 30 nt of genomic context and is therefore only reported by design_guide_rnas or calculate_on_target_rule_set2. Position is -1 for ad-hoc evaluation. Call to QC one guide sequence.")]
    public static GuideRnaCandidate evaluate_guide_rna(
        [Description("Guide RNA sequence.")] string guide_sequence,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9,
        [Description("Optional guide-RNA design parameters.")] GuideRnaParameters? parameters = null)
    {
        if (string.IsNullOrEmpty(guide_sequence))
            throw new System.ArgumentException("Guide sequence cannot be null or empty.", nameof(guide_sequence));

        return CrisprDesigner.EvaluateGuideRna(guide_sequence, system_type, parameters);
    }

    [McpServerTool(Name = "find_off_targets", Title = "MolTools — CRISPR Off-Target Scan", ReadOnly = true), Description("Naïve genome scan: enumerates all PAM sites in the genome and reports those whose target differs from the guide by 1..max_mismatches (the 0-mismatch on-target is excluded). Off-target score weights mismatches inside the seed region more heavily. O(genome × guide) — recommend genome ≤ ~1 Mb. Guide length must match the system's guide length. Call to list a guide's candidate off-target sites.")]
    public static OffTargetsResult find_off_targets(
        [Description("Guide RNA sequence (length must match the system's guide length).")] string guide_sequence,
        [Description("Genome / reference sequence to scan.")] string genome,
        [Description("Maximum allowed mismatches (range 0..5; default 3).")] int max_mismatches = 3,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9)
    {
        if (string.IsNullOrEmpty(guide_sequence))
            throw new System.ArgumentException("Guide sequence cannot be null or empty.", nameof(guide_sequence));
        if (string.IsNullOrEmpty(genome))
            throw new System.ArgumentException("Genome cannot be null or empty.", nameof(genome));
        if (max_mismatches < 0 || max_mismatches > 5)
            throw new System.ArgumentException("Maximum mismatches must be in the range 0..5.", nameof(max_mismatches));

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

    [McpServerTool(Name = "calculate_on_target_doench2014", Title = "MolTools — On-Target Score (Doench 2014 Rule Set 1)", ReadOnly = true), Description("Doench et al. 2014 \"Rule Set 1\" on-target efficacy score for an SpCas9 guide, returned on a 0..100 scale (higher = predicted more active). This is the published logistic linear model (intercept + GC term over the protospacer + position-specific single/di-nucleotide weights, through a sigmoid). Input is the model's 30-nt context: 4 nt upstream + 20 nt protospacer + 3 nt PAM (must be NGG) + 3 nt downstream, A/C/G/T only. Call to rank guides by predicted cutting efficiency when the 30-nt context is known.")]
    public static OnTargetScoreResult calculate_on_target_doench2014(
        [Description("30-nt context: 4 nt upstream + 20 nt protospacer + 3 nt NGG PAM + 3 nt downstream.")] string context_30mer)
    {
        if (string.IsNullOrEmpty(context_30mer))
            throw new System.ArgumentException("Context 30-mer cannot be null or empty.", nameof(context_30mer));

        return new OnTargetScoreResult(CrisprDesigner.CalculateOnTargetDoench2014(context_30mer));
    }

    [McpServerTool(Name = "calculate_on_target_rule_set2", Title = "MolTools — On-Target Score (Doench 2016 Rule Set 2 / Azimuth)", ReadOnly = true), Description("Doench et al. 2016 \"Rule Set 2\" / Azimuth on-target efficacy score for an SpCas9 guide, conventionally in 0..1 (higher = predicted more active) — the \"Doench '16\" efficiency score reported by CRISPOR. Rule Set 2 is a trained gradient-boosted-tree model, reproduced here from Microsoft Research's Azimuth model. Input is the model's 30-nt context: 4 nt upstream + 20 nt protospacer + 3 nt PAM (must be NGG) + 3 nt downstream, A/C/G/T only. Pass amino_acid_cut_position and percent_peptide together to use Azimuth's gene-context (full) model instead of the sequence-only one. Call to rank guides by the published on-target activity model.")]
    public static OnTargetScoreResult calculate_on_target_rule_set2(
        [Description("30-nt context: 4 nt upstream + 20 nt protospacer + 3 nt NGG PAM + 3 nt downstream.")] string context_30mer,
        [Description("Optional amino-acid position of the cut site in the target protein (requires percent_peptide); enables the gene-context model.")] int? amino_acid_cut_position = null,
        [Description("Optional cut position as a percentage 0..100 along the coding sequence (requires amino_acid_cut_position).")] double? percent_peptide = null)
    {
        if (string.IsNullOrEmpty(context_30mer))
            throw new System.ArgumentException("Context 30-mer cannot be null or empty.", nameof(context_30mer));
        if (amino_acid_cut_position.HasValue != percent_peptide.HasValue)
            throw new System.ArgumentException(
                "amino_acid_cut_position and percent_peptide must be supplied together (gene-context model) or both omitted.",
                nameof(amino_acid_cut_position));

        double score = amino_acid_cut_position.HasValue
            ? CrisprDesigner.CalculateOnTargetRuleSet2(context_30mer, amino_acid_cut_position.Value, percent_peptide!.Value)
            : CrisprDesigner.CalculateOnTargetRuleSet2(context_30mer);
        return new OnTargetScoreResult(score);
    }

    #endregion

    #region ProbeDesigner

    [McpServerTool(Name = "design_probes", Title = "MolTools — Design Hybridization Probes", ReadOnly = true), Description("Designs hybridization probes by scanning the target for length-window candidates and ranking them with an additive penalty score (GC%, Tm, homopolymers, self-structure, simple repeats; returned sorted by score, descending). Tm is Primer3's seqtm (SantaLucia 1998 nearest-neighbour ≤ 36 nt, long_seq_tm above) at the parameters' conditions (default Primer3 probe conditions: 50 nM, 50 mM monovalent, no Mg/dNTP); probes ≤ 60 nt are screened with Primer3's ntthal self-dimer/hairpin Tm limit (47 °C). Use one of the ProbeParameters presets (Microarray | FISH | NorthernBlot | qPCR | SouthernBlot) or pass custom values; default = Microarray. Returns up to max_probes top-scoring probes; a target shorter than the minimum probe length yields an empty list.")]
    public static ProbesResult design_probes(
        [Description("Target DNA sequence.")] string target_sequence,
        [Description("Optional probe-design parameters (lengths, Tm range, GC range, max homopolymer, self-complementarity threshold). Defaults to Microarray when null.")] ProbeDesigner.ProbeParameters? parameters = null,
        [Description("Maximum probes to return (default 10).")] int max_probes = 10)
    {
        if (string.IsNullOrEmpty(target_sequence))
            throw new System.ArgumentException("Target sequence cannot be null or empty.", nameof(target_sequence));
        if (max_probes <= 0)
            throw new System.ArgumentException("Maximum probes must be positive.", nameof(max_probes));

        var probes = ProbeDesigner.DesignProbes(target_sequence, parameters, max_probes).ToList();
        return new ProbesResult(probes);
    }

    [McpServerTool(Name = "design_tiling_probes", Title = "MolTools — Design Tiling Probes", ReadOnly = true), Description("Generates fixed-length probes covering the entire target with a configurable overlap (step = probe_length − overlap). Sub-optimal candidates are still emitted (each with a 'Suboptimal probe' warning) so coverage is preserved. Returns the probe set plus covered-position count, mean Tm, and Tm range.")]
    public static ProbeDesigner.TilingProbeSet design_tiling_probes(
        [Description("Target DNA sequence.")] string target_sequence,
        [Description("Tiling probe length in bp (default 60).")] int probe_length = 60,
        [Description("Overlap between adjacent probes in bp (default 20).")] int overlap = 20,
        [Description("Optional probe-design parameters (Tm/GC bounds etc.).")] ProbeDesigner.ProbeParameters? parameters = null)
    {
        if (string.IsNullOrEmpty(target_sequence))
            throw new System.ArgumentException("Target sequence cannot be null or empty.", nameof(target_sequence));
        if (probe_length <= 0)
            throw new System.ArgumentException("Probe length must be positive.", nameof(probe_length));
        if (overlap < 0 || overlap >= probe_length)
            throw new System.ArgumentException("Overlap must be non-negative and less than the probe length.", nameof(overlap));

        return ProbeDesigner.DesignTilingProbes(target_sequence, probe_length, overlap, parameters);
    }

    [McpServerTool(Name = "design_antisense_probes", Title = "MolTools — Design Antisense Probes", ReadOnly = true), Description("Reverse-complements the supplied mRNA-sense sequence and runs the probe designer on it; every returned probe is tagged type=Antisense. Returns up to max_probes top-scoring antisense probes.")]
    public static ProbesResult design_antisense_probes(
        [Description("mRNA-sense sequence (will be reverse-complemented).")] string mrna_sequence,
        [Description("Optional probe-design parameters.")] ProbeDesigner.ProbeParameters? parameters = null,
        [Description("Maximum probes to return (default 5).")] int max_probes = 5)
    {
        if (string.IsNullOrEmpty(mrna_sequence))
            throw new System.ArgumentException("mRNA sequence cannot be null or empty.", nameof(mrna_sequence));
        if (max_probes <= 0)
            throw new System.ArgumentException("Maximum probes must be positive.", nameof(max_probes));

        var probes = ProbeDesigner.DesignAntisenseProbes(mrna_sequence, parameters, max_probes).ToList();
        return new ProbesResult(probes);
    }

    [McpServerTool(Name = "design_molecular_beacon", Title = "MolTools — Design Molecular Beacon", ReadOnly = true), Description("Designs a hairpin molecular-beacon probe: GC-rich complementary stems (stem5 = ⌊stem_length/2⌋ Gs + remaining Cs, stem3 = its reverse complement) flanking the best target-specific loop of probe_length bases, for real-time detection. The reported Tm is the loop (probe–target) Tm (Primer3 seqtm) and Start/End mark the loop in the target; warnings carry the ntthal stem-loop Tm. With detection_temperature T the loop Tm window is [T+7, T+10] °C and the stem-loop Tm is checked against T+7 °C (Tyagi & Kramer molecular-beacon rules). Returns probe=null when the target is shorter than probe_length.")]
    public static MolecularBeaconResult design_molecular_beacon(
        [Description("Target DNA sequence.")] string target_sequence,
        [Description("Loop (target-specific) length in bp (default 25).")] int probe_length = 25,
        [Description("Stem length in bp (default 5).")] int stem_length = 5,
        [Description("Optional detection (annealing) temperature in °C for the Tyagi & Kramer 7–10 °C rules.")] double? detection_temperature = null)
    {
        if (string.IsNullOrEmpty(target_sequence))
            throw new System.ArgumentException("Target sequence cannot be null or empty.", nameof(target_sequence));
        if (probe_length <= 0)
            throw new System.ArgumentException("Probe (loop) length must be positive.", nameof(probe_length));
        if (stem_length <= 0)
            throw new System.ArgumentException("Stem length must be positive.", nameof(stem_length));

        return new MolecularBeaconResult(
            ProbeDesigner.DesignMolecularBeacon(target_sequence, probe_length, stem_length, detection_temperature));
    }

    [McpServerTool(Name = "validate_probe", Title = "MolTools — Validate Probe Specificity", ReadOnly = true), Description("Validates a hybridization probe. (1) Ungapped k-mismatch (Hamming) scan of the reference sequences: off-target hit count (intended site included; > 1 hit is an issue) and a library uniqueness score (0 hits → 0.0, N hits → 1/N). (2) Self-structure: for ≤ 60-nt A/C/G/T probes Primer3's thermodynamic probe screen — ntthal self-dimer, 3′ self-dimer and hairpin Tm at 50 nM oligo / 50 mM monovalent / no Mg²⁺ (Primer3 probe conditions) must not exceed 47 °C (PRIMER_INTERNAL_MAX_*_TH); longer or non-ACGT probes use the sequence-only fold-back fraction (> self_complementarity_threshold) and inverted-repeat screens. (3) Optional non_target_sequences: Kane et al. (2000) cross-hybridization criteria on both strands — overall identity of the best local (BLAST-scored Smith–Waterman–Gotoh) alignment over the probe length > 75 % or a contiguous identical stretch > 15 nt; each site also reports its ntthal duplex Tm with the probe (primer3-py calc_heterodimer), optionally thresholded by max_duplex_tm (OligoArray 2.0). isValid = no issue recorded. Call to check whether a designed probe is specific and structure-free.")]
    public static ProbeDesigner.ProbeValidation validate_probe(
        [Description("Probe sequence to validate.")] string probe_sequence,
        [Description("Reference sequences to scan for off-target hits.")] string[] reference_sequences,
        [Description("Maximum allowed mismatches (default 3).")] int max_mismatches = 3,
        [Description("Fold-back-fraction limit of the sequence-only fallback screen (default 0.3).")] double self_complementarity_threshold = 0.3,
        [Description("Optional known non-target sequences for the Kane et al. (2000) cross-hybridization criteria (both strands).")] string[]? non_target_sequences = null,
        [Description("Kane identity threshold in [0,1]; a non-target strand with identity strictly above it is flagged (default 0.75).")] double max_non_target_identity = 0.75,
        [Description("Kane contiguous-identity threshold in nt; a longer identical stretch is flagged (default 15).")] int max_contiguous_match = 15,
        [Description("Optional OligoArray-style threshold (°C): a non-target site whose ntthal duplex Tm with the probe is above it is flagged (default none).")] double? max_duplex_tm = null)
    {
        if (probe_sequence is null)
            throw new System.ArgumentException("Probe sequence cannot be null.", nameof(probe_sequence));
        if (reference_sequences is null)
            throw new System.ArgumentException("Reference sequences cannot be null.", nameof(reference_sequences));
        if (max_mismatches < 0)
            throw new System.ArgumentException("Maximum mismatches cannot be negative.", nameof(max_mismatches));
        if (double.IsNaN(max_non_target_identity) || max_non_target_identity < 0 || max_non_target_identity > 1)
            throw new System.ArgumentException("Non-target identity threshold must be in [0, 1].", nameof(max_non_target_identity));
        if (max_contiguous_match < 0)
            throw new System.ArgumentException("Contiguous-match threshold cannot be negative.", nameof(max_contiguous_match));

        return ProbeDesigner.ValidateProbe(probe_sequence, reference_sequences, max_mismatches, self_complementarity_threshold,
            nonTargetSequences: non_target_sequences,
            maxNonTargetIdentity: max_non_target_identity,
            maxContiguousMatch: max_contiguous_match,
            maxDuplexTm: max_duplex_tm);
    }

    [McpServerTool(Name = "design_probes_primer3", Title = "MolTools — Primer3 Hybridization-Probe Picker", ReadOnly = true), Description("Picks hybridization probes exactly as Primer3 does for PRIMER_TASK=pick_hyb_probe_only (internal-oligo picker; verified against primer3-py design_primers): every A/C/G/T window of min_size..max_size within the G+C % window, poly-X ≤ max_poly_x, Primer3 seqtm Tm within [min_tm, max_tm] and ntthal self-dimer / 3′ self-dimer / hairpin Tm ≤ their limits, enumerated per 3′ end with Primer3's 5′-extension break; ranked by the Primer3 penalty |Tm − opt_tm| + |length − opt_size| (penalty ascending, then start descending, then length ascending). Defaults are Primer3's PRIMER_INTERNAL_* defaults (18/20/27 nt, Tm 57/60/63 °C, GC 20–80 %, poly-X 5, 47 °C structure limits, 50 mM monovalent, no Mg²⁺/dNTP, 50 nM). Returns up to num_return probes (start is 0-based).")]
    public static Primer3ProbesResult design_probes_primer3(
        [Description("Template DNA sequence (probes are picked on this strand; case-insensitive).")] string template,
        [Description("PRIMER_NUM_RETURN: maximum probes to return (default 5).")] int num_return = 5,
        [Description("PRIMER_INTERNAL_MIN_SIZE (default 18).")] int min_size = 18,
        [Description("PRIMER_INTERNAL_OPT_SIZE (default 20).")] int opt_size = 20,
        [Description("PRIMER_INTERNAL_MAX_SIZE (default 27; at most 36).")] int max_size = 27,
        [Description("PRIMER_INTERNAL_MIN_TM in °C (default 57).")] double min_tm = 57.0,
        [Description("PRIMER_INTERNAL_OPT_TM in °C (default 60).")] double opt_tm = 60.0,
        [Description("PRIMER_INTERNAL_MAX_TM in °C (default 63).")] double max_tm = 63.0,
        [Description("PRIMER_INTERNAL_MIN_GC in percent (default 20).")] double min_gc_percent = 20.0,
        [Description("PRIMER_INTERNAL_MAX_GC in percent (default 80).")] double max_gc_percent = 80.0,
        [Description("PRIMER_INTERNAL_MAX_POLY_X (default 5).")] int max_poly_x = 5,
        [Description("PRIMER_INTERNAL_MAX_SELF_ANY_TH in °C (default 47).")] double max_self_any_th = PrimerDesigner.Primer3MaxStructureTm,
        [Description("PRIMER_INTERNAL_MAX_SELF_END_TH in °C (default 47).")] double max_self_end_th = PrimerDesigner.Primer3MaxStructureTm,
        [Description("PRIMER_INTERNAL_MAX_HAIRPIN_TH in °C (default 47).")] double max_hairpin_th = PrimerDesigner.Primer3MaxStructureTm,
        [Description("PRIMER_INTERNAL_SALT_MONOVALENT in mM (default 50).")] double monovalent_mm = PrimerDesigner.Primer3InternalMonovalentMillimolar,
        [Description("PRIMER_INTERNAL_SALT_DIVALENT (Mg²⁺) in mM (default 0).")] double divalent_mm = PrimerDesigner.Primer3InternalDivalentMillimolar,
        [Description("PRIMER_INTERNAL_DNTP_CONC in mM (default 0).")] double dntp_mm = PrimerDesigner.Primer3InternalDntpMillimolar,
        [Description("PRIMER_INTERNAL_DNA_CONC in nM (default 50).")] double dna_conc_nm = PrimerDesigner.Primer3InternalDnaConcentrationNanomolar)
    {
        if (string.IsNullOrEmpty(template))
            throw new System.ArgumentException("Template sequence cannot be null or empty.", nameof(template));
        if (num_return < 0)
            throw new System.ArgumentException("num_return cannot be negative.", nameof(num_return));
        if (min_size < 1 || max_size < min_size || max_size > 36)
            throw new System.ArgumentException("Sizes must satisfy 1 ≤ min_size ≤ max_size ≤ 36.", nameof(max_size));

        var settings = new ProbeDesigner.Primer3ProbeSettings(
            min_size, opt_size, max_size, min_tm, opt_tm, max_tm, min_gc_percent, max_gc_percent, max_poly_x,
            max_self_any_th, max_self_end_th, max_hairpin_th, monovalent_mm, divalent_mm, dntp_mm, dna_conc_nm);
        return new Primer3ProbesResult(ProbeDesigner.DesignProbesPrimer3(template, settings, num_return));
    }

    [McpServerTool(Name = "analyze_oligo", Title = "MolTools — Oligonucleotide Property Analysis", ReadOnly = true), Description("Returns Tm, GC fraction, molecular weight (Da), and 260 nm extinction coefficient (M⁻¹·cm⁻¹) for a short oligonucleotide. Call when the user needs the basic physical properties of an oligo/primer/probe. Tm is Primer3's seqtm at the Primer3 hybridization-probe conditions (50 nM oligo, 50 mM monovalent, no Mg/dNTP; SantaLucia 1998 nearest-neighbour for ≤ 36 nt, long_seq_tm above) and is null when not computable (fewer than 2 bases or a non-ACGT base, e.g. RNA). Molecular weight is the single-stranded Biopython molecular_weight (RNA when the oligo has U and no T); ε260 is the mononucleotide sum. GC is returned as a fraction (0-1).")]
    public static OligoAnalysisResult analyze_oligo(
        [Description("Oligonucleotide sequence (non-empty; A/C/G/T/U, case-insensitive).")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        var (tm, gc, mw, eps) = ProbeDesigner.AnalyzeOligo(sequence);
        return new OligoAnalysisResult(double.IsNaN(tm) ? null : tm, gc, mw, eps);
    }

    [McpServerTool(Name = "oligo_extinction_coefficient", Title = "MolTools — Oligo Extinction Coefficient", ReadOnly = true), Description("Estimates an oligonucleotide's 260 nm molar extinction coefficient (M⁻¹·cm⁻¹). Default: sum of per-base contributions (A=15400, C=7400, G=11500, T=8700, U=9900; any other base = 10000). With nearest_neighbor=true: the nearest-neighbour model (Cantor, Warshaw & Shapiro 1970 DNA / Warshaw & Tinoco 1966 RNA table, ε = Σ ε(dinucleotides) − Σ ε(internal mononucleotides)). Call to estimate an oligo's ε₂₆₀ for concentration calculations.")]
    public static ExtinctionCoefficientResult oligo_extinction_coefficient(
        [Description("Oligonucleotide sequence.")] string sequence,
        [Description("Use the nearest-neighbour model (default false = mononucleotide sum).")] bool nearest_neighbor = false,
        [Description("For nearest_neighbor: true = DNA table (A/C/G/T), false = RNA table (A/C/G/U) (default true).")] bool is_dna = true)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        if (!nearest_neighbor)
            return new ExtinctionCoefficientResult(ProbeDesigner.CalculateExtinctionCoefficient(sequence));

        double eps = ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor(sequence, is_dna);
        if (double.IsNaN(eps))
            throw new System.ArgumentException(
                $"Nearest-neighbour ε260 needs only {(is_dna ? "A/C/G/T" : "A/C/G/U")} bases.", nameof(sequence));
        return new ExtinctionCoefficientResult(eps);
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
