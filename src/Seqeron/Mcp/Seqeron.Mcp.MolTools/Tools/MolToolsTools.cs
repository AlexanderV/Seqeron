using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;

namespace Seqeron.Mcp.MolTools.Tools;

[McpServerToolType]
public class MolToolsTools
{
    /// <summary>
    /// The wrappers' shared argument guard: a null or empty text is an <see cref="System.ArgumentException"/>
    /// (not <see cref="System.ArgumentNullException"/>) saying "<c>{label} cannot be null or empty.</c>". It replaces a
    /// three-line guard that was written out at every call (B08 round 5 duplication sweep).
    /// </summary>
    private static void RequireNotEmpty(string? value, string paramName, string label)
    {
        if (string.IsNullOrEmpty(value))
            throw new System.ArgumentException($"{label} cannot be null or empty.", paramName);
    }

    // Utility holder for static MCP tools; never instantiated (S1118).
    private MolToolsTools() { }

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

    [McpServerTool(Name = "longest_dinucleotide_repeat", Title = "MolTools — Longest Dinucleotide Repeat", ReadOnly = true), Description("Returns the number of repeat units in the longest dinucleotide tandem repeat (e.g. ATATAT = 3 units of AT; any 2-mer counts, so AAAA = 2 units of AA), case-insensitive. Sequences shorter than 4 nt return 0. Library screen (unsourced heuristic): Primer3 has no dinucleotide-repeat limit and no published definition of this screen or of the library default limit 4 was found; Primer3's sourced alternatives are PRIMER_MAX_POLY_X (longest_homopolymer) and repeat-library mispriming (design_primers mispriming_library). Call to flag microsatellite-like dinucleotide repeats in a primer/probe.")]
    public static DinucleotideRepeatResult longest_dinucleotide_repeat(
        [Description("Nucleotide sequence.")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        return new DinucleotideRepeatResult(PrimerDesigner.FindLongestDinucleotideRepeat(sequence));
    }

    [McpServerTool(Name = "hairpin_potential", Title = "MolTools — Hairpin Potential", ReadOnly = true), Description("Sequence-only library hairpin screen (heuristic, not Primer3's): true when an exact Watson-Crick self-complementary stem of at least min_stem_length bp is separated by a loop of at least min_loop_length nt (no wobble, mismatches or energies). The default loop minimum 3 is Primer3's (thal.c min_hrpn_loop); the default stem minimum 4 is an unsourced library threshold. It can disagree with Primer3 (e.g. AAAACCCTTTT is flagged although Primer3 calc_hairpin finds no structure). For Primer3's sourced hairpin screen (ntthal hairpin Tm vs PRIMER_MAX_HAIRPIN_TH 47 °C) use evaluate_primer (hairpinTh). Uses an O(n²) scan for short sequences and a suffix-tree scan for sequences ≥ 100 bp. Call for a quick exact-stem check of a primer/probe.")]
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

    [McpServerTool(Name = "get_enzyme", Title = "MolTools — Look Up Restriction Enzyme", ReadOnly = true), Description("Looks up a restriction enzyme by name (case-insensitive) among the 754 REBASE enzymes with a defined cut (Biopython 1.88 REBASE snapshot; the 39 common laboratory enzymes included, Type IIS such as BsaI and two-sided Type IIB such as BaeI too) and returns its recognition sequence, cut positions, organism, end chemistry (blunt / 5' / 3') and the overhang sequence it leaves — the top-strand bases between its two strand cuts, 5'→3', which is REBASE/Biopython ovhgseq (EcoRI \"AATT\", PstI \"TGCA\", empty for a blunt cutter, and an IUPAC code where the enzyme's own site is degenerate, e.g. SfiI \"NNN\"). Also lists the other enzymes of enzyme_set (default \"common\" = the 39 common enzymes; \"rebase\" = all 754) recognising the same sequence, split into equischizomers (same site AND same cut positions, so interchangeable in a digest — MspI/HpaII) and neoschizomers (same site, different cut, so different ends — MboI vs DpnI). recognitionLength is the full site string, so the interrupted palindrome SfiI (GGCCNNNNNGGCC) reports 13 although only 8 bases are specified. rebaseMethylationField is REBASE allenz.104 field <6> verbatim (EcoRI \"3(6)\") and cognateMethylationSites its typed parse (1-based position in the site or null for \"?\", complementary-strand flag, m4C/m5C/m6A): the base the enzyme's own cognate methyltransferase modifies — NOT whether Dam/Dcm/CpG methylation blocks the enzyme, which is not modelled. A two-sided (Type IIB) cutter also reports secondCutForward/secondCutReverse, its second cut pair in the same convention (BaeI 23/18); both are null for a one-cut enzyme. Returns enzyme=null when the name is not a REBASE enzyme with a defined cut (an unknown name, or one of the 334 REBASE enzymes whose cut position is undefined, e.g. Dpi3084I). Call to fetch a single enzyme's properties by name.")]
    public static EnzymeLookupResult get_enzyme(
        [Description("Enzyme name (e.g. EcoRI, BamHI, BsaI).")] string name,
        [Description("Which enzyme table the equischizomer/neoschizomer lists are drawn from: \"common\" (default; the 39 common laboratory enzymes) or \"rebase\" (all 754 REBASE enzymes with a defined cut). Case-insensitive. The looked-up enzyme itself is always resolved against all 754.")] string? enzyme_set = "common")
    {
        if (string.IsNullOrWhiteSpace(name))
            throw new System.ArgumentException("Enzyme name cannot be null or blank.", nameof(name));
        var set = ParseEnzymeSet(enzyme_set);

        var enzyme = RestrictionAnalyzer.GetEnzyme(name);
        if (enzyme is null)
            return new EnzymeLookupResult(null);

        return new EnzymeLookupResult(
            enzyme,
            RestrictionAnalyzer.GetOverhangSequence(enzyme),
            RestrictionAnalyzer.FindEquischizomers(enzyme.Name, set).Select(e => e.Name).ToList(),
            RestrictionAnalyzer.FindNeoschizomers(enzyme.Name, set).Select(e => e.Name).ToList());
    }

    /// <summary>
    /// Parses the <c>enzyme_set</c> argument of the listing tools: <c>"common"</c> (default, also for
    /// <see langword="null"/>) or <c>"rebase"</c>, case-insensitive (B08 D4 part 2).
    /// </summary>
    private static EnzymeSet ParseEnzymeSet(string? enzyme_set)
    {
        if (enzyme_set is null || string.Equals(enzyme_set, "common", System.StringComparison.OrdinalIgnoreCase))
            return EnzymeSet.Common;
        if (string.Equals(enzyme_set, "rebase", System.StringComparison.OrdinalIgnoreCase))
            return EnzymeSet.Rebase;
        throw new System.ArgumentException(
            $"Unknown enzyme set '{enzyme_set}'; expected \"common\" or \"rebase\".", nameof(enzyme_set));
    }

    [McpServerTool(Name = "enzymes_by_cut_length", Title = "MolTools — Enzymes by Recognition Length", ReadOnly = true), Description("Lists the restriction enzymes of enzyme_set (default \"common\" = the 39 common laboratory enzymes; \"rebase\" = all 754 REBASE enzymes with a defined cut) whose recognition sequence has exactly the specified length in base pairs (e.g. 4 for frequent cutters, 6 for typical cloning enzymes, 8 for rare cutters), or, when max_length is set, every enzyme whose recognition string is length..max_length long \u2014 so \"all 4-6 bp cutters\" is one call rather than three. The length measured is the full recognition string (Biopython .size), so an interrupted palindrome counts its spacer: SfiI (GGCCNNNNNGGCC) is returned by length 13, not by length 8, although only 8 of its bases are specified. Call to pick enzymes by how often they cut.")]
    public static EnzymeListResult enzymes_by_cut_length(
        [Description("Recognition-sequence length in bp (must be positive). With max_length set this is the inclusive lower bound of a range instead.")] int length,
        [Description("Which enzyme table to list from: \"common\" (default; the 39 common laboratory enzymes) or \"rebase\" (all 754 REBASE enzymes with a defined cut, the 39 included). Case-insensitive.")] string? enzyme_set = "common",
        [Description("Optional inclusive upper bound: with it set, every enzyme whose recognition string is length..max_length long is returned, so \"all 4-6 bp cutters\" is one call. Must be positive; a max_length below length returns an empty list, as the range is then empty. Omit for the exact-length listing.")] int? max_length = null)
    {
        if (length <= 0)
            throw new System.ArgumentException("Recognition-sequence length must be positive.", nameof(length));
        if (max_length is <= 0)
            throw new System.ArgumentException("max_length must be positive.", nameof(max_length));
        var set = ParseEnzymeSet(enzyme_set);

        // The ranged C# overload had no MCP expression, so "all 4-6 bp cutters" took three calls
        // (review 2026-09 B08 H28).
        return new EnzymeListResult(max_length is { } max
            ? RestrictionAnalyzer.GetEnzymesByCutLength(length, max, set).ToList()
            : RestrictionAnalyzer.GetEnzymesByCutLength(length, set).ToList());
    }

    [McpServerTool(Name = "blunt_cutters", Title = "MolTools — Blunt-End Restriction Enzymes", ReadOnly = true), Description("Lists the restriction enzymes of enzyme_set (default \"common\" = the 39 common laboratory enzymes; \"rebase\" = all 754 REBASE enzymes with a defined cut) that produce blunt ends (both strands cut at the same position, Biopython is_blunt). Any blunt end ligates to any other blunt end. Call when the user wants blunt-cutting enzymes for blunt-end cloning.")]
    public static EnzymeListResult blunt_cutters(
        [Description("Which enzyme table to list from: \"common\" (default; the 39 common laboratory enzymes) or \"rebase\" (all 754 REBASE enzymes with a defined cut, the 39 included). Case-insensitive.")] string? enzyme_set = "common")
    {
        var set = ParseEnzymeSet(enzyme_set);
        return new EnzymeListResult(RestrictionAnalyzer.GetBluntCutters(set).ToList());
    }

    [McpServerTool(Name = "sticky_cutters", Title = "MolTools — Sticky-End Restriction Enzymes", ReadOnly = true), Description("Lists the restriction enzymes of enzyme_set (default \"common\" = the 39 common laboratory enzymes; \"rebase\" = all 754 REBASE enzymes with a defined cut) that produce sticky (cohesive) ends — a staggered cut leaving a 5' or 3' single-stranded overhang (Biopython is_5overhang / is_3overhang); the complement of blunt_cutters. Use get_enzyme for the overhang sequence an individual enzyme leaves, and enzymes_compatible to test whether two of them can be ligated. Call when the user wants overhang-generating enzymes for directional / sticky-end cloning.")]
    public static EnzymeListResult sticky_cutters(
        [Description("Which enzyme table to list from: \"common\" (default; the 39 common laboratory enzymes) or \"rebase\" (all 754 REBASE enzymes with a defined cut, the 39 included). Case-insensitive.")] string? enzyme_set = "common")
    {
        var set = ParseEnzymeSet(enzyme_set);
        return new EnzymeListResult(RestrictionAnalyzer.GetStickyCutters(set).ToList());
    }

    [McpServerTool(Name = "find_restriction_sites", Title = "MolTools — Find Restriction Sites", ReadOnly = true), Description("Finds restriction sites for one or more named enzymes (any of the 754 REBASE enzymes with a defined cut, case-insensitive; an unknown name is an error) on both strands of a DNA sequence. IUPAC degenerate codes in recognition sequences are matched against ACGT input. position and cutPosition are 0-based forward-strand coordinates for both strands; cutPosition is the number of forward-strand bases 5' of the cut made on that record's own strand, so a palindromic enzyme reports a forward record (top-strand cut) and a reverse record (bottom-strand cut) at the same position, and recognizedSequence is the window as read on that strand. Overlapping sites are all reported. On a Linear molecule a site is reported only when both of its cuts fall on a bond that exists (so an enzyme cutting at the edge of its own site, e.g. Sau3AI ^GATC, is not reported at the very ends); set topology=Circular for a plasmid, where the site may span the origin and cut positions wrap modulo the length. Call to locate where specific enzymes cut a sequence. Optional methylation (none|dam|dcm|dam+dcm|emboss) makes methylated bases block recognition exactly as EMBOSS restrict -methylation does (default none: unchanged).")]
    public static RestrictionSiteListResult find_restriction_sites(
        [Description("DNA sequence to scan. IUPAC ambiguity codes are accepted anywhere in the sequence (N, R, Y, ...), so an N-containing contig or a soft-masked assembly can be scanned as it is: a site N matches any of the 15 IUPAC codes, every other site code matches only a concrete A/C/G/T base it expands to, and a character outside the 15 codes (U, a space, a FASTA header line) matches nothing \u2014 Biopython\u0027s own rule; such text is not rejected, but it still occupies coordinates and cannot be part of a site. Lower case is accepted; the returned fragment and site sequences are upper-case (the soft mask is not carried over; positions and lengths are unaffected).")] string sequence,
        [Description("Names of restriction enzymes to look for.")] string[] enzyme_names,
        [Description("Linear (default) or Circular: whether the molecule is a linear fragment or a closed circle (plasmid).")] MoleculeTopology topology = MoleculeTopology.Linear,
        [Description("Optional methylation-aware recognition (EMBOSS restrict -methylation semantics): \"none\" (default), \"dam\" (GmATC), \"dcm\" (mCmCAGG, mCmCTGG), \"dam+dcm\" or \"emboss\" (all five Emethylsites.dat entries incl. EcoKI). Case-insensitive. A recognition window containing a methylated base is NOT a site (a window rule, not an enzyme-specific sensitivity table; BamHI is blocked by Dam).")] string? methylation = null)
    {
        RequireNotEmpty(sequence, nameof(sequence), "Sequence");
        if (enzyme_names is null || enzyme_names.Length == 0)
            throw new System.ArgumentException("At least one enzyme name is required.", nameof(enzyme_names));

        // Raw-string entry point: an N-containing contig or a soft-masked assembly is a sequence
        // this scan handles, and DnaSequence would reject it (review 2026-09 B08 H27).
        var sites = RestrictionAnalyzer.FindSites(
            sequence, topology, MethylaseSets.FromName(methylation), enzyme_names).ToList();
        return new RestrictionSiteListResult(sites);
    }

    [McpServerTool(Name = "find_all_restriction_sites", Title = "MolTools — Find All Restriction Sites", ReadOnly = true), Description("Finds sites for EVERY restriction enzyme of enzyme_set (default \"common\" = the 39 common laboratory enzymes; \"rebase\" = all 754 REBASE enzymes with a defined cut) on both strands of a DNA sequence. Coordinates, the two records per palindromic site and the linear cut-feasibility rule are the same as find_restriction_sites; set topology=Circular for a plasmid. Call for a comprehensive restriction-site scan when the enzyme set is not known in advance. Optional methylation (none|dam|dcm|dam+dcm|emboss) makes methylated bases block recognition exactly as EMBOSS restrict -methylation does (default none: unchanged).")]
    public static RestrictionSiteListResult find_all_restriction_sites(
        [Description("DNA sequence to scan. IUPAC ambiguity codes are accepted anywhere in the sequence (N, R, Y, ...), so an N-containing contig or a soft-masked assembly can be scanned as it is: a site N matches any of the 15 IUPAC codes, every other site code matches only a concrete A/C/G/T base it expands to, and a character outside the 15 codes (U, a space, a FASTA header line) matches nothing \u2014 Biopython\u0027s own rule; such text is not rejected, but it still occupies coordinates and cannot be part of a site. Lower case is accepted; the returned fragment and site sequences are upper-case (the soft mask is not carried over; positions and lengths are unaffected).")] string sequence,
        [Description("Linear (default) or Circular: whether the molecule is a linear fragment or a closed circle (plasmid).")] MoleculeTopology topology = MoleculeTopology.Linear,
        [Description("Which enzyme table to list from: \"common\" (default; the 39 common laboratory enzymes) or \"rebase\" (all 754 REBASE enzymes with a defined cut, the 39 included). Case-insensitive.")] string? enzyme_set = "common",
        [Description("Optional methylation-aware recognition (EMBOSS restrict -methylation semantics): \"none\" (default), \"dam\" (GmATC), \"dcm\" (mCmCAGG, mCmCTGG), \"dam+dcm\" or \"emboss\" (all five Emethylsites.dat entries incl. EcoKI). Case-insensitive. A recognition window containing a methylated base is NOT a site (a window rule, not an enzyme-specific sensitivity table; BamHI is blocked by Dam).")] string? methylation = null)
    {
        RequireNotEmpty(sequence, nameof(sequence), "Sequence");
        var set = ParseEnzymeSet(enzyme_set);

        var sites = RestrictionAnalyzer.FindAllSites(sequence, topology, set, MethylaseSets.FromName(methylation)).ToList();
        return new RestrictionSiteListResult(sites);
    }

    [McpServerTool(Name = "restriction_digest", Title = "MolTools — Restriction Digest", ReadOnly = true), Description("Simulates a restriction digest with one or more named enzymes (any of the 754 REBASE enzymes with a defined cut, case-insensitive; an unknown name is an error) and yields the resulting fragments in ascending start position (each with its sequence, start, length, flanking enzymes and the blunt/5'/3' overhang of each cut end). Set topology=Circular for a plasmid: the site search then wraps past the origin, k cut positions give k fragments (vs k+1 linear), and the origin-spanning fragment is reported last; an enzyme that does not cut leaves the circle whole. Call to predict the fragment pattern of a digest. Optional methylation (none|dam|dcm|dam+dcm|emboss) makes methylated bases block recognition exactly as EMBOSS restrict -methylation does (default none: unchanged).")]
    public static DigestResult restriction_digest(
        [Description("DNA sequence to digest. IUPAC ambiguity codes are accepted anywhere in the sequence (N, R, Y, ...), so an N-containing contig or a soft-masked assembly can be digested as it is: a site N matches any of the 15 IUPAC codes, every other site code matches only a concrete A/C/G/T base it expands to, and a character outside the 15 codes (U, a space, a FASTA header line) matches nothing \u2014 Biopython\u0027s own rule; such text is not rejected, but it still occupies coordinates and cannot be part of a site. Lower case is accepted; the returned fragment and site sequences are upper-case (the soft mask is not carried over; positions and lengths are unaffected).")] string sequence,
        [Description("Names of restriction enzymes to use (must be non-empty).")] string[] enzyme_names,
        [Description("Linear (default) or Circular: whether the molecule is a linear fragment or a closed circle (plasmid).")] MoleculeTopology topology = MoleculeTopology.Linear,
        [Description("Optional methylation-aware recognition (EMBOSS restrict -methylation semantics): \"none\" (default), \"dam\" (GmATC), \"dcm\" (mCmCAGG, mCmCTGG), \"dam+dcm\" or \"emboss\" (all five Emethylsites.dat entries incl. EcoKI). Case-insensitive. A recognition window containing a methylated base is NOT a site (a window rule, not an enzyme-specific sensitivity table; BamHI is blocked by Dam).")] string? methylation = null)
    {
        RequireNotEmpty(sequence, nameof(sequence), "Sequence");
        if (enzyme_names is null || enzyme_names.Length == 0)
            throw new System.ArgumentException("At least one enzyme name is required.", nameof(enzyme_names));

        var fragments = RestrictionAnalyzer.Digest(
            sequence, topology, MethylaseSets.FromName(methylation), enzyme_names).ToList();
        return new DigestResult(fragments);
    }

    [McpServerTool(Name = "digest_summary", Title = "MolTools — Restriction Digest Summary", ReadOnly = true), Description("Aggregate statistics over a simulated restriction digest with one or more named enzymes (any of the 754 REBASE enzymes with a defined cut, case-insensitive; an unknown name is an error): total fragment count, fragment sizes (descending), largest/smallest fragment, average fragment size, enzymes used and the topology digested. Set topology=Circular for a plasmid, which yields k fragments for k cut positions instead of k+1. Call for a gel-like summary of a digest without the full fragment sequences. Optional methylation (none|dam|dcm|dam+dcm|emboss) makes methylated bases block recognition exactly as EMBOSS restrict -methylation does (default none: unchanged).")]
    public static DigestSummary digest_summary(
        [Description("DNA sequence to digest. IUPAC ambiguity codes are accepted anywhere in the sequence (N, R, Y, ...), so an N-containing contig or a soft-masked assembly can be digested as it is: a site N matches any of the 15 IUPAC codes, every other site code matches only a concrete A/C/G/T base it expands to, and a character outside the 15 codes (U, a space, a FASTA header line) matches nothing \u2014 Biopython\u0027s own rule; such text is not rejected, but it still occupies coordinates and cannot be part of a site. Lower case is accepted; the returned fragment and site sequences are upper-case (the soft mask is not carried over; positions and lengths are unaffected).")] string sequence,
        [Description("Names of restriction enzymes to use.")] string[] enzyme_names,
        [Description("Linear (default) or Circular: whether the molecule is a linear fragment or a closed circle (plasmid).")] MoleculeTopology topology = MoleculeTopology.Linear,
        [Description("Optional methylation-aware recognition (EMBOSS restrict -methylation semantics): \"none\" (default), \"dam\" (GmATC), \"dcm\" (mCmCAGG, mCmCTGG), \"dam+dcm\" or \"emboss\" (all five Emethylsites.dat entries incl. EcoKI). Case-insensitive. A recognition window containing a methylated base is NOT a site (a window rule, not an enzyme-specific sensitivity table; BamHI is blocked by Dam).")] string? methylation = null)
    {
        RequireNotEmpty(sequence, nameof(sequence), "Sequence");
        if (enzyme_names is null || enzyme_names.Length == 0)
            throw new System.ArgumentException("At least one enzyme name is required.", nameof(enzyme_names));

        return RestrictionAnalyzer.GetDigestSummary(
            sequence, topology, MethylaseSets.FromName(methylation), enzyme_names);
    }

    [McpServerTool(Name = "partial_digest", Title = "MolTools — Partial Restriction Digest", ReadOnly = true), Description("Simulates a PARTIAL restriction digest: every fragment that can arise when any subset of the cut sites is left uncut, each with missedCuts = the number of uncut sites inside it. Enzymes are named as for restriction_digest (any of the 754 REBASE enzymes with a defined cut); uses exactly its cut positions, enzyme lists and end overhangs; the missedCuts = 0 fragments come first and are the complete digest. Linear: all intervals between two of {0, cuts, length} — (k+1)(k+2)/2 fragments for k cuts. Circular: every cut-to-cut arc (wrapping past the origin; k-1 missed = the molecule linearised at one cut) plus the intact circle (missedCuts = k) — k²+1 fragments. Order: missedCuts, then start position. max_missed_cuts bounds the uncut sites per fragment; above 100000 fragments the call is refused unless max_missed_cuts is lowered. Which fragments CAN arise only — no yield/abundance. Call to predict extra bands from an incomplete digest. Optional methylation (none|dam|dcm|dam+dcm|emboss) makes methylated bases block recognition exactly as EMBOSS restrict -methylation does (default none: unchanged). A blocked site is neither cut nor counted as a missed cut.")]
    public static PartialDigestResult partial_digest(
        [Description("DNA sequence to digest. IUPAC ambiguity codes are accepted anywhere in the sequence (N, R, Y, ...), so an N-containing contig or a soft-masked assembly can be digested as it is: a site N matches any of the 15 IUPAC codes, every other site code matches only a concrete A/C/G/T base it expands to, and a character outside the 15 codes (U, a space, a FASTA header line) matches nothing \u2014 Biopython\u0027s own rule; such text is not rejected, but it still occupies coordinates and cannot be part of a site. Lower case is accepted; the returned fragment and site sequences are upper-case (the soft mask is not carried over; positions and lengths are unaffected).")] string sequence,
        [Description("Names of restriction enzymes to use (must be non-empty).")] string[] enzyme_names,
        [Description("Linear (default) or Circular: whether the molecule is a linear fragment or a closed circle (plasmid).")] MoleculeTopology topology = MoleculeTopology.Linear,
        [Description("Optional upper bound on the number of uncut sites inside one fragment (0 = complete digest); omit for no bound.")] int? max_missed_cuts = null,
        [Description("Optional methylation-aware recognition (EMBOSS restrict -methylation semantics): \"none\" (default), \"dam\" (GmATC), \"dcm\" (mCmCAGG, mCmCTGG), \"dam+dcm\" or \"emboss\" (all five Emethylsites.dat entries incl. EcoKI). Case-insensitive. A recognition window containing a methylated base is NOT a site (a window rule, not an enzyme-specific sensitivity table; BamHI is blocked by Dam).")] string? methylation = null)
    {
        RequireNotEmpty(sequence, nameof(sequence), "Sequence");
        if (enzyme_names is null || enzyme_names.Length == 0)
            throw new System.ArgumentException("At least one enzyme name is required.", nameof(enzyme_names));

        var fragments = RestrictionAnalyzer.PartialDigest(
            sequence, topology, max_missed_cuts, MethylaseSets.FromName(methylation), enzyme_names).ToList();
        return new PartialDigestResult(fragments);
    }

    [McpServerTool(Name = "restriction_map", Title = "MolTools — Restriction Map", ReadOnly = true), Description("Builds a restriction map of a DNA sequence: all forward+reverse sites, sites grouped by enzyme, the total number of CLEAVAGES, the unique cutters (enzymes that cleave the molecule exactly once), and the non-cutters. Enzyme names (any of the 754 REBASE enzymes with a defined cut) are resolved case-insensitively and de-duplicated, and every result is keyed by the table's own spelling. An empty enzyme_names list considers every enzyme of enzyme_set (default \"common\" = the 39 common laboratory enzymes; \"rebase\" = all 754), so nonCutters then lists every such enzyme without a site; enzyme_set applies only to an empty list; sitesByEnzyme and nonCutters always partition the considered set. Set topology=Circular for a plasmid. totalSites and uniqueCutters count cleavages, not recognition windows: every forward-strand record counts, plus every reverse-strand record whose backwards reading cuts other bonds. So a palindromic site with symmetric cuts (all 39 common enzymes) counts once, a non-palindromic enzyme\u0027s reverse-strand sites count on their own — a Type IIS enzyme such as BsaI is NOT a unique cutter just because it has one forward-strand site, since GGTCTC on the other strand is a second cleavage — and the five enzymes whose two strand readings of one window cut different bonds (HauII, SgeI, MspJI, LpnPI, SgrTI) count twice for that one window. Call to plan cloning around single-cutter enzymes. Optional methylation (none|dam|dcm|dam+dcm|emboss) makes methylated bases block recognition exactly as EMBOSS restrict -methylation does (default none: unchanged).")]
    public static RestrictionMap restriction_map(
        [Description("DNA sequence to map. IUPAC ambiguity codes are accepted anywhere in the sequence (N, R, Y, ...), so an N-containing contig or a soft-masked assembly can be mapped as it is: a site N matches any of the 15 IUPAC codes, every other site code matches only a concrete A/C/G/T base it expands to, and a character outside the 15 codes (U, a space, a FASTA header line) matches nothing \u2014 Biopython\u0027s own rule; such text is not rejected, but it still occupies coordinates and cannot be part of a site. Lower case is accepted; the returned fragment and site sequences are upper-case (the soft mask is not carried over; positions and lengths are unaffected).")] string sequence,
        [Description("Names of restriction enzymes to consider; empty to use every enzyme of enzyme_set.")] string[] enzyme_names,
        [Description("Linear (default) or Circular: whether the molecule is a linear fragment or a closed circle (plasmid).")] MoleculeTopology topology = MoleculeTopology.Linear,
        [Description("Enzyme table used when enzyme_names is empty: \"common\" (default; the 39 common laboratory enzymes) or \"rebase\" (all 754 REBASE enzymes with a defined cut). Case-insensitive; ignored when names are given.")] string? enzyme_set = "common",
        [Description("Optional methylation-aware recognition (EMBOSS restrict -methylation semantics): \"none\" (default), \"dam\" (GmATC), \"dcm\" (mCmCAGG, mCmCTGG), \"dam+dcm\" or \"emboss\" (all five Emethylsites.dat entries incl. EcoKI). Case-insensitive. A recognition window containing a methylated base is NOT a site (a window rule, not an enzyme-specific sensitivity table; BamHI is blocked by Dam).")] string? methylation = null)
    {
        RequireNotEmpty(sequence, nameof(sequence), "Sequence");
        var set = ParseEnzymeSet(enzyme_set);
        var names = enzyme_names ?? System.Array.Empty<string>();
        var methylases = MethylaseSets.FromName(methylation);

        // The full set goes through the record overload (same CreateMapCore, B08 F58); the common set
        // keeps the name overload so the default output is unchanged.
        if (names.Length == 0 && set == EnzymeSet.Rebase)
            return RestrictionAnalyzer.CreateMap(
                sequence, topology, RestrictionAnalyzer.RebaseEnzymes.Values, methylases);

        return RestrictionAnalyzer.CreateMap(sequence, topology, methylases, names);
    }

    [McpServerTool(Name = "compatible_enzymes", Title = "MolTools — Compatible Enzyme Pairs", ReadOnly = true), Description("Enumerates pairs of restriction enzymes whose ends can be ligated to each other: either both produce blunt ends, or both produce the same polarity of overhang (both 5' or both 3') and overhangs that can anneal — same length, and one the reverse complement of the other. compatibleEnd is \"blunt\" or the FIRST listed enzyme\u0027s own overhang window, and end2 the second\u0027s; they are equal unless the pair anneals only through the reverse complement (12 REBASE pairs, all GN against CN). A 5' overhang is never compatible with a 3' overhang even when the letters match (SacI 3' AGCT vs HindIII 5' AGCT are NOT compatible). Without enzyme_name every unordered pair of distinct enzymes of enzyme_set is listed (default \"common\" = the 39 laboratory enzymes, 64 pairs; \"rebase\" = all 754 REBASE enzymes with a defined cut, 27869 pairs, table order, e.g. EcoRI/MfeI share the 5' end AATT). The all-pairs response is bounded: it returns at most max_pairs pairs (default 1000, maximum 10000) starting at offset, with totalPairs (all pairs of the set) and truncated (true when more pairs follow the returned page) — page through a large set with offset. With enzyme_name set (any of the 754 REBASE enzymes with a defined cut), lists what that one enzyme can be ligated to among enzyme_set (default \"common\" = the 39; \"rebase\" = all 754), itself included when it is in that set — use this for an enzyme with a unique overhang (EcoRI, PstI, SfiI), which appears in no pair at all; this per-enzyme form is not paged and rejects max_pairs / offset. Call when the user needs enzymes that yield compatible (ligatable) ends for cloning.")]
    public static CompatibleEnzymesResult compatible_enzymes(
        [Description("Optional: restrict the result to pairs involving this enzyme (case-insensitive); omit or leave blank for every pair of enzyme_set.")] string? enzyme_name = null,
        [Description("Table to list: \"common\" (default; the 39 common laboratory enzymes) or \"rebase\" (all 754 REBASE enzymes with a defined cut). Case-insensitive.")] string? enzyme_set = "common",
        [Description("All-pairs listing only: maximum number of pairs returned (1-10000; default 1000, which holds the whole common set). Not allowed together with enzyme_name.")] int? max_pairs = null,
        [Description("All-pairs listing only: number of pairs to skip before the returned page (>= 0; default 0). Not allowed together with enzyme_name.")] int? offset = null)
    {
        var set = ParseEnzymeSet(enzyme_set);
        if (string.IsNullOrWhiteSpace(enzyme_name))
        {
            int limit = max_pairs ?? DefaultCompatiblePairsPage;
            int skip = offset ?? 0;
            if (limit < 1 || limit > MaxCompatiblePairsPage)
                throw new System.ArgumentOutOfRangeException(
                    nameof(max_pairs), limit, $"max_pairs must be between 1 and {MaxCompatiblePairsPage}.");
            if (skip < 0)
                throw new System.ArgumentOutOfRangeException(nameof(offset), skip, "offset cannot be negative.");

            var page = new List<EnzymeCompatibilityPair>();
            int total = 0;
            foreach (var t in RestrictionAnalyzer.FindCompatibleEnzymes(set))
            {
                if (total >= skip && page.Count < limit)
                {
                    page.Add(new EnzymeCompatibilityPair(
                        t.Enzyme1, t.Enzyme2, t.CompatibleEnd, OtherEnd(t.Enzyme2, t.CompatibleEnd)));
                }
                total++;
            }

            return new CompatibleEnzymesResult(page, total, Truncated: skip + page.Count < total);
        }

        if (max_pairs is not null || offset is not null)
            throw new System.ArgumentException(
                "max_pairs and offset page the all-pairs listing; they cannot be combined with enzyme_name.",
                max_pairs is not null ? nameof(max_pairs) : nameof(offset));

        var partners = RestrictionAnalyzer.FindCompatibleEnzymes(enzyme_name, set).ToList();
        // The enzyme the caller named, under its table spelling (" bamhi " and "CviKI-1" resolve; echoing the raw text would
        // report an enzyme1 that no other tool knows). FindCompatibleEnzymes has just refused an unknown name.
        string canonicalName = RestrictionAnalyzer.GetEnzyme(enzyme_name)!.Name;
        var pairs = partners
            .Select(t => new EnzymeCompatibilityPair(
                canonicalName, t.Enzyme, t.CompatibleEnd, OtherEnd(t.Enzyme, t.CompatibleEnd)))
            .ToList();
        return new CompatibleEnzymesResult(pairs);
    }

    private const int DefaultCompatiblePairsPage = 1000;
    private const int MaxCompatiblePairsPage = 10000;

    [McpServerTool(Name = "enzymes_compatible", Title = "MolTools — Enzyme End Compatibility", ReadOnly = true), Description("Returns whether two named restriction enzymes (any of the 754 REBASE enzymes with a defined cut, case-insensitive) produce ligatable (compatible) ends: both blunt, or the same overhang polarity (both 5' or both 3') with overhangs that can anneal — same length, and one the reverse complement of the other. compatibleEnd is \"blunt\" for a blunt pair, otherwise the FIRST enzyme\u0027s own overhang window — not necessarily a string both enzymes leave — and is null when they are incompatible; end2 is the second enzyme\u0027s window, reported so the asymmetry is visible. The two are the same whenever the windows intersect directly, but 12 REBASE pairs anneal only through the reverse complement and there they differ: all 12 are GN against CN (BsmI/BsrI, Bse1I/BsmI, BsrI/LmnI and so on), and e.g. BsrI never leaves CN. Same letters with opposite polarity are NOT compatible (SacI 3' AGCT vs HindIII 5' AGCT → false). An unknown enzyme name is an ERROR, as it is for find_restriction_sites, so a typo is never reported as \"not compatible\"; use get_enzyme first if the name may not exist. Call to check a specific pair for cloning compatibility.")]
    public static EnzymeCompatibilityResult enzymes_compatible(
        [Description("First enzyme name.")] string enzyme1_name,
        [Description("Second enzyme name.")] string enzyme2_name)
    {
        if (string.IsNullOrWhiteSpace(enzyme1_name))
            throw new System.ArgumentException("First enzyme name cannot be null or blank.", nameof(enzyme1_name));
        if (string.IsNullOrWhiteSpace(enzyme2_name))
            throw new System.ArgumentException("Second enzyme name cannot be null or blank.", nameof(enzyme2_name));

        bool compatible = RestrictionAnalyzer.AreCompatible(enzyme1_name, enzyme2_name);
        if (!compatible)
            return new EnzymeCompatibilityResult(false, null, null);

        // AreCompatible resolved the partner through GetEnzyme, which trims the name and maps a REBASE hyphen onto the
        // table's underscore ("CviKI-1" -> CviKI_1); the shared end is therefore looked up by the CANONICAL name, not
        // by the caller's text (the raw text made First() throw for " HpaII" and "CviKI-1"; B08 round 5 A5R-1).
        string partner = RestrictionAnalyzer.GetEnzyme(enzyme2_name)!.Name;
        // Search the full set: enzyme2 may be any REBASE enzyme (GetEnzyme resolves all 754), and
        // a common-set search would not contain it (B08 WP12b).
        string? compatibleEnd = RestrictionAnalyzer.FindCompatibleEnzymes(enzyme1_name, EnzymeSet.Rebase)
            .First(t => string.Equals(t.Enzyme, partner, System.StringComparison.OrdinalIgnoreCase))
            .CompatibleEnd;
        return new EnzymeCompatibilityResult(true, compatibleEnd, OtherEnd(partner, compatibleEnd));
    }

    /// <summary>
    /// The partner enzyme's own overhang window, so a caller sees both ends of a compatible pair.
    /// <c>compatibleEnd</c> is the first/queried enzyme's window, which is not necessarily a string
    /// the partner leaves: 12 REBASE pairs anneal only through the reverse complement and there the
    /// two differ (all GN against CN). <c>"blunt"</c> is shared by definition (review 2026-09 B08 H30).
    /// </summary>
    private static string? OtherEnd(string otherEnzymeName, string? compatibleEnd)
        => compatibleEnd == "blunt" ? "blunt" : RestrictionAnalyzer.GetOverhangSequence(otherEnzymeName);

    [McpServerTool(Name = "find_pattern_occurrences", Title = "MolTools \u2014 Find Pattern Occurrences", ReadOnly = true), Description("Enumerates every 0-based start index at which an IUPAC recognition pattern occurs in a sequence, on the GIVEN STRAND ONLY. This is the raw pattern scan underneath find_restriction_sites, exposed for a caller that holds a pattern rather than an enzyme: no enzyme, no cut arithmetic, no reverse-strand scan and no linear-molecule cut-feasibility filter. Overlapping occurrences are all reported, in ascending order. Matching is Biopython Bio.Restriction's site-regex rule: a pattern N matches any of the 15 IUPAC codes in the sequence (an N included), every other pattern code matches only a concrete A/C/G/T base it expands to \u2014 so \"GAANTC\" does NOT match GAATTC and \"GTYRAC\" does NOT match GTYRAC \u2014 and a character outside the 15 codes matches nothing. Both the sequence and the pattern may carry IUPAC codes and the scan is case-insensitive. To cover both strands, call again with the reverse complement of the pattern. A null or empty sequence yields no occurrences; an invalid pattern is an error. Call to locate a recognition pattern, a degenerate site or a cloning scar when no enzyme is involved.")]
    public static PatternOccurrencesResult find_pattern_occurrences(
        [Description("Sequence to scan. IUPAC ambiguity codes are accepted anywhere; null or empty yields no occurrences.")] string? sequence,
        [Description("Recognition pattern: 1 or more of the 15 IUPAC nucleotide codes A, C, G, T, N, R, Y, S, W, K, M, B, D, H, V (case-insensitive). At most 64 characters.")] string recognition_pattern,
        [Description("First index to consider; occurrences starting before it are skipped. Defaults to 0, must not be negative.")] int start_index = 0)
    {
        // The C# raw pattern scan had no MCP expression at all, so the documented IUPAC-tolerant
        // scan was unreachable over MCP (review 2026-09 B08 H27/H28).
        var positions = RestrictionAnalyzer
            .FindPatternOccurrences(sequence, recognition_pattern, start_index)
            .ToList();
        return new PatternOccurrencesResult(positions);
    }

    /// <summary>
    /// The caller-supplied enzyme of the two custom-enzyme tools: a recognition sequence of IUPAC codes, a top- and a
    /// bottom-strand cut offset and, for a two-sided (Type IIB) cutter, a second pair (both or neither).
    /// </summary>
    private static RestrictionEnzyme BuildCustomEnzyme(
        string? name, string recognition_sequence, int cut_position_forward, int cut_position_reverse,
        int? second_cut_forward, int? second_cut_reverse)
    {
        if ((second_cut_forward is null) != (second_cut_reverse is null))
        {
            throw new System.ArgumentException(
                "Pass both second_cut_forward and second_cut_reverse, or neither.",
                second_cut_forward is null ? nameof(second_cut_forward) : nameof(second_cut_reverse));
        }

        // The IEnumerable<RestrictionEnzyme> overloads had no MCP expression, so a caller-supplied
        // or cut-undefined REBASE enzyme could not be used at all (review 2026-09 B08 H28).
        return new RestrictionEnzyme(
            string.IsNullOrWhiteSpace(name) ? "CustomEnzyme" : name.Trim(),
            recognition_sequence,
            cut_position_forward,
            cut_position_reverse,
            Organism: "")
        {
            SecondCut = second_cut_forward is { } f && second_cut_reverse is { } r ? (f, r) : null,
        };
    }

    [McpServerTool(Name = "find_sites_custom_enzyme", Title = "MolTools \u2014 Find Sites of a Custom Enzyme", ReadOnly = true), Description("Finds the sites of a CALLER-SUPPLIED enzyme definition \u2014 one that is not in the 754-enzyme REBASE table, or one of the 334 REBASE enzymes whose cut position REBASE leaves undefined and which therefore cannot be resolved by name. Give the recognition sequence (IUPAC codes allowed, need not be palindromic) and the two cut offsets, counted from the first base of the site: cut_position_forward is the number of site bases 5' of the top-strand cut and cut_position_reverse the same for the bottom-strand cut, so EcoRI G^AATTC / GAATT^C is (1, 5) and the Type IIS BsaI GGTCTCN^NNNN is (7, 11). Offsets outside the site are allowed and are how a distal (Type IIS/IIB) cutter is expressed; a negative offset cuts upstream of the site. second_cut_forward/second_cut_reverse add the second cut pair of a two-sided (Type IIB) cutter such as BaeI (23, 18) and must leave the same stagger as the first pair. Coordinates, the two records per palindromic site, the linear cut-feasibility rule and the circular wrap are exactly as in find_restriction_sites. Call when the enzyme is not in the table \u2014 use find_restriction_sites for any named one.")]
    public static RestrictionSiteListResult find_sites_custom_enzyme(
        [Description("DNA sequence to scan. IUPAC ambiguity codes are accepted anywhere in the sequence.")] string sequence,
        [Description("Recognition sequence of the custom enzyme: 1 or more of the 15 IUPAC nucleotide codes; need not be palindromic. At most 64 bases (REBASE's longest is 15).")] string recognition_sequence,
        [Description("Top-strand cut offset: the number of recognition-sequence bases 5' of the cut. May be negative or exceed the site length for a distal cutter.")] int cut_position_forward,
        [Description("Bottom-strand cut offset in the same frame. For a palindromic enzyme cut_position_forward + cut_position_reverse equals the site length. The two cuts may be at most 64 bases apart (the stagger, i.e. the overhang length); a larger one is refused.")] int cut_position_reverse,
        [Description("Name reported on each site record; defaults to \"CustomEnzyme\". At most 64 characters.")] string? name = null,
        [Description("Linear (default) or Circular: whether the molecule is a linear fragment or a closed circle (plasmid).")] MoleculeTopology topology = MoleculeTopology.Linear,
        [Description("Optional second top-strand cut offset of a two-sided (Type IIB) cutter; pass both second offsets or neither.")] int? second_cut_forward = null,
        [Description("Optional second bottom-strand cut offset; must leave the same stagger as the first pair.")] int? second_cut_reverse = null,
        [Description("Optional methylation-aware recognition (EMBOSS restrict -methylation semantics): \"none\" (default), \"dam\", \"dcm\", \"dam+dcm\" or \"emboss\". A recognition window containing a methylated base is NOT a site.")] string? methylation = null)
    {
        RequireNotEmpty(sequence, nameof(sequence), "Sequence");
        var enzyme = BuildCustomEnzyme(
            name, recognition_sequence, cut_position_forward, cut_position_reverse, second_cut_forward, second_cut_reverse);

        var sites = RestrictionAnalyzer
            .FindSites(sequence, topology, [enzyme], MethylaseSets.FromName(methylation))
            .ToList();
        return new RestrictionSiteListResult(sites);
    }

    [McpServerTool(Name = "restriction_digest_custom_enzyme", Title = "MolTools \u2014 Digest with a Custom Enzyme", ReadOnly = true), Description("Simulates a restriction digest with a CALLER-SUPPLIED enzyme definition and yields the resulting fragments, each with its sequence, start, length, flanking enzymes and the blunt/5'/3' overhang of each cut end. The enzyme is given exactly as in find_sites_custom_enzyme (recognition sequence plus the two cut offsets, optionally a second pair for a two-sided Type IIB cutter). Fragment conventions, the circular k-fragments-for-k-cuts rule and the methylation mode are exactly as in restriction_digest. Call when the enzyme is not one of the 754 named REBASE enzymes \u2014 for a named one use restriction_digest.")]
    public static DigestResult restriction_digest_custom_enzyme(
        [Description("DNA sequence to digest. IUPAC ambiguity codes are accepted anywhere in the sequence.")] string sequence,
        [Description("Recognition sequence of the custom enzyme: 1 or more of the 15 IUPAC nucleotide codes; need not be palindromic. At most 64 bases (REBASE's longest is 15).")] string recognition_sequence,
        [Description("Top-strand cut offset: the number of recognition-sequence bases 5' of the cut. May be negative or exceed the site length for a distal cutter.")] int cut_position_forward,
        [Description("Bottom-strand cut offset in the same frame. The two cuts may be at most 64 bases apart (the stagger, i.e. the overhang length); a larger one is refused.")] int cut_position_reverse,
        [Description("Name reported on each fragment's flanking-enzyme fields; defaults to \"CustomEnzyme\". At most 64 characters.")] string? name = null,
        [Description("Linear (default) or Circular: whether the molecule is a linear fragment or a closed circle (plasmid).")] MoleculeTopology topology = MoleculeTopology.Linear,
        [Description("Optional second top-strand cut offset of a two-sided (Type IIB) cutter; pass both second offsets or neither.")] int? second_cut_forward = null,
        [Description("Optional second bottom-strand cut offset; must leave the same stagger as the first pair.")] int? second_cut_reverse = null,
        [Description("Optional methylation-aware recognition (EMBOSS restrict -methylation semantics): \"none\" (default), \"dam\", \"dcm\", \"dam+dcm\" or \"emboss\".")] string? methylation = null)
    {
        RequireNotEmpty(sequence, nameof(sequence), "Sequence");
        var enzyme = BuildCustomEnzyme(
            name, recognition_sequence, cut_position_forward, cut_position_reverse, second_cut_forward, second_cut_reverse);

        var fragments = RestrictionAnalyzer
            .Digest(sequence, topology, [enzyme], MethylaseSets.FromName(methylation))
            .ToList();
        return new DigestResult(fragments);
    }

    [McpServerTool(Name = "count_codons", Title = "MolTools — Count Codons", ReadOnly = true), Description("Counts occurrences of each codon in a coding sequence (frame 0, non-overlapping triplets). DNA or RNA, case-insensitive: U is read as T and keys are uppercase DNA codons. A triplet containing any other character (e.g. N) is skipped without shifting the reading frame; a trailing partial codon is ignored. Independent of the genetic code. Call to get the raw codon-frequency table of a gene.")]
    public static CodonCountsResult count_codons(
        [Description("Coding sequence, DNA or RNA (frame 0).")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));

        return new CodonCountsResult(CodonUsageAnalyzer.CountCodons(sequence));
    }

    [McpServerTool(Name = "rscu", Title = "MolTools — Relative Synonymous Codon Usage", ReadOnly = true), Description("Relative Synonymous Codon Usage (RSCU) per codon: observed count / count-expected-if-uniform among its synonymous codons. RSCU = 1 means no bias, > 1 over-represented, < 1 under-represented. Synonymous families (and the stop family) come from the NCBI genetic code genetic_code (default 1 = Standard; RSCU is genetic-code dependent, CodonW -code). Optional pseudocount c (default 0 = CodonW, absent amino acids report 0) is added to every codon count before the per-family normalisation, RSCU = n(x + c)/(sum x + n c), as codonbias RelativeSynonymousCodonUsage(pseudocount) and cubar est_rscu(pseudo_cnt); with c > 0 an absent amino acid is uniform (1). Optional family_level Subfamily groups codons as the R package cubar does by default (est_rscu level = \"subfam\": amino acid + first two codon bases, so Leu/Ser/Arg split into a 4-fold and a 2-fold sub-family; the stop codons by their first two bases too); default AminoAcid (Sharp & Li / CodonW families). Call to quantify codon bias per codon in a coding sequence.")]
    public static RscuResult rscu(
        [Description("Coding DNA sequence (frame 0).")] string sequence,
        [Description("NCBI translation table number defining the synonymous codon families and stop codons (default 1 = Standard; e.g. 2 = vertebrate mitochondrial, 11 = bacterial/plastid). Must be one of the NCBI gc.prt tables (1-6, 9-16, 21-33).")] int genetic_code = 1,
        [Description("Pseudocount added to every codon count before normalising within its synonymous family (finite, >= 0; default 0 = no pseudocount, CodonW; codonbias and cubar default to 1).")] double pseudocount = 0,
        [Description("Synonymous grouping: AminoAcid (default; one family per amino acid plus the stop family, Sharp & Li / CodonW) or Subfamily (cubar est_rscu level = \"subfam\", its default: amino acid + first two codon bases, e.g. Leu_TT {TTA, TTG} and Leu_CT {CTN}).")] RscuFamilyLevel family_level = RscuFamilyLevel.AminoAcid)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        var code = ResolveGeneticCode(genetic_code);
        if (!double.IsFinite(pseudocount) || pseudocount < 0)
            throw new System.ArgumentException($"pseudocount must be a finite non-negative number (got {pseudocount}).", nameof(pseudocount));

        if (!System.Enum.IsDefined(family_level))
            throw new System.ArgumentException(
                $"Unknown RSCU family level: {family_level}. Supported: AminoAcid, Subfamily.", nameof(family_level));

        return new RscuResult(CodonUsageAnalyzer.CalculateRscu(sequence, code, pseudocount, family_level));
    }

    [McpServerTool(Name = "codon_adaptation_index", Title = "MolTools — Codon Adaptation Index (RSCU)", ReadOnly = true), Description("Codon Adaptation Index (Sharp & Li 1987) using a caller-supplied reference RSCU table (typically derived from highly expressed genes), codons in the DNA alphabet. Output range 0..1; per-codon relative adaptiveness w = RSCU/max-synonymous-RSCU, and CAI is their geometric mean. Single-codon amino acids (Met/Trp) and stop codons are excluded; a codon with w < 0.0001 (absent from the reference) is scored with w = 0.01 (CodonW convention); non-ACGT(U) triplets are skipped. Stops, single-codon families and the synonymous families used to derive w follow the NCBI genetic code genetic_code (default 1 = Standard; CodonW -code), e.g. under table 2 ATA joins Met and TGA joins Trp while AGA/AGG are stops. Call to score how well a gene matches an organism's preferred codons given an RSCU reference (distinct from cai_from_organism_table, which takes a frequency table).")]
    public static CaiResult codon_adaptation_index(
        [Description("Coding DNA sequence (frame 0).")] string sequence,
        [Description("Reference RSCU table: codon (DNA alphabet) → RSCU value.")] Dictionary<string, double> reference_rscu,
        [Description("NCBI translation table number defining the synonymous codon families and stop codons (default 1 = Standard; e.g. 2 = vertebrate mitochondrial, 11 = bacterial/plastid). Must be one of the NCBI gc.prt tables (1-6, 9-16, 21-33).")] int genetic_code = 1)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (reference_rscu is null)
            throw new System.ArgumentException("Reference RSCU table cannot be null.", nameof(reference_rscu));
        var code = ResolveGeneticCode(genetic_code);

        return new CaiResult(CodonUsageAnalyzer.CalculateCai(sequence, reference_rscu, code));
    }

    [McpServerTool(Name = "effective_number_of_codons", Title = "MolTools — Effective Number of Codons (ENC)", ReadOnly = true), Description("Effective Number of Codons (Wright 1990 Nc; CodonW enc_out conventions), measuring how far a gene departs from uniform synonymous-codon usage under the NCBI genetic code genetic_code (default 1 = Standard; CodonW -code). Range 20..61 under the standard code (20 = extreme bias, 61 = no bias; overshoot re-adjusted to the code's sense-codon count, e.g. 60 for table 2). Amino acids seen once, or with every observed codon used once, are not estimable; a missing isoleucine class is replaced by the mean of the 2- and 4-fold classes. When Nc cannot be calculated (some other synonymous class has no estimable amino acid - gene too short or amino-acid usage too skewed; CodonW prints *****) returns enc = 0 and calculable = false; otherwise calculable = true. DNA or RNA, case-insensitive; non-ACGT(U) triplets are skipped. Optional method selects another estimator: Novembre2002 (background-corrected Nc', chi-square of the codon frequencies against those expected from the base composition, unweighted class means, no pseudocount) or SunYangXia2013 (pseudocount 1, F = (chi-square + 1)/k, count-weighted class means), both as in codonbias 0.5.0, where a class with no estimable amino acid is imputed with F = 1/k so the value is always calculable (calculable = false only without any sense codon); pseudocount / background_correction / background tune them. EncPrimeNc / EncPrimeNcPrime reproduce the Nc / Ncp columns of ENCprime (Novembre's reference program): only amino acids observed more than 5 times are estimated, a missing class is imputed with F = 1/k (3-fold: mean of the 2- and 4-fold classes), class means are floored at 1/k, no pseudocount and no cap; Nc' uses the gene's own A/C/G/T composition or background (calculable = false only for an empty sequence). codonw_cap re-adjusts Wright Nc to CodonW's fixed 61 instead of the code's sense-codon count. Call to summarise a gene's overall codon bias with a single number.")]
    public static EncResult effective_number_of_codons(
        [Description("Coding DNA sequence (frame 0).")] string sequence,
        [Description("NCBI translation table number defining the synonymous codon families and stop codons (default 1 = Standard; e.g. 2 = vertebrate mitochondrial, 11 = bacterial/plastid). Must be one of the NCBI gc.prt tables (1-6, 9-16, 21-33).")] int genetic_code = 1,
        [Description("Estimator: Wright1990 (default; CodonW enc_out), Novembre2002 (background-corrected Nc', codonbias bg_correction=True, robust=False, pseudocount=0, unweighted) , SunYangXia2013 (codonbias defaults: robust, pseudocount 1, weighted), EncPrimeNc or EncPrimeNcPrime (ENCprime program Nc / Ncp).")] EncMethod method = EncMethod.Wright1990,
        [Description("Novembre2002/SunYangXia2013 only: pseudocount added to every sense-codon count (>= 0; omitted = 0 for Novembre2002, 1 for SunYangXia2013).")] double? pseudocount = null,
        [Description("Novembre2002/SunYangXia2013 only: correct the expected codon frequencies for the background base composition (omitted = true for Novembre2002, false for SunYangXia2013).")] bool? background_correction = null,
        [Description("Novembre2002/SunYangXia2013 with background correction, or EncPrimeNcPrime: sequence whose A/C/G/T composition is the background (omitted = the gene itself; ENCprime: the SeqCount -n nucleotide-composition input).")] string? background = null,
        [Description("Wright1990 only: re-adjust Nc to CodonW's fixed 61 (codon_us.c enc_out) instead of the genetic code's sense-codon count (differs for codes with other than 61 sense codons, e.g. table 2 caps at 60, tables 3/9 at 62). Default false.")] bool codonw_cap = false)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        var code = ResolveGeneticCode(genetic_code);
        if (!System.Enum.IsDefined(method))
            throw new System.ArgumentException(
                $"Unknown ENC method: {method}. Supported: Wright1990, Novembre2002, SunYangXia2013, EncPrimeNc, EncPrimeNcPrime.", nameof(method));
        bool encPrime = method is EncMethod.EncPrimeNc or EncMethod.EncPrimeNcPrime;
        if (encPrime && (pseudocount is not null || background_correction is not null))
            throw new System.ArgumentException(
                "pseudocount and background_correction do not apply to the ENCprime methods (fixed program rules).", nameof(method));
        if (method == EncMethod.EncPrimeNc && background is not null)
            throw new System.ArgumentException("background applies to EncPrimeNcPrime only (ENCprime's Nc ignores the base composition).", nameof(background));
        if (method == EncMethod.Wright1990 && (pseudocount is not null || background_correction is not null || background is not null))
            throw new System.ArgumentException(
                "pseudocount, background_correction and background apply to the Novembre2002 and SunYangXia2013 methods only.", nameof(method));
        if (method != EncMethod.Wright1990 && codonw_cap)
            throw new System.ArgumentException("codonw_cap applies to the Wright1990 method only.", nameof(codonw_cap));
        if (pseudocount is double pc && (!double.IsFinite(pc) || pc < 0))
            throw new System.ArgumentException($"pseudocount must be a finite non-negative number (got {pc}).", nameof(pseudocount));
        if (background is not null && (background.Length == 0
            || !(method == EncMethod.EncPrimeNcPrime || (background_correction ?? method == EncMethod.Novembre2002))))
            throw new System.ArgumentException(
                "background must be a non-empty sequence and requires background correction.", nameof(background));

        var options = new EncOptions
        {
            Method = method,
            Pseudocount = pseudocount,
            BackgroundCorrection = background_correction,
            Background = background,
            CodonWCap = codonw_cap,
        };
        bool calculable = CodonUsageAnalyzer.TryCalculateEnc(sequence, code, options, out double enc);
        return new EncResult(enc, calculable);
    }

    [McpServerTool(Name = "codon_usage_statistics", Title = "MolTools — Codon-Usage Statistics", ReadOnly = true), Description("Aggregate codon-usage report for a coding sequence: per-codon counts, RSCU, Effective Number of Codons (ENC; 0 with isEncCalculable = false when not calculable, see effective_number_of_codons), total codons, GC% at codon positions 1/2/3, GC3s (synonymous third-position GC), and overall GC. RSCU, ENC and GC3s use the NCBI genetic code genetic_code (default 1 = Standard; CodonW -code). gc_convention selects GC1/2/3 and overall GC: EmbossCusp (default) counts every codon including stop codons (EMBOSS cusp) and does not depend on the genetic code; CodonW excludes the stop codons of genetic_code (CodonW gc_out; overall GC = CodonW GC index). Set include_fop_cbi = true to also report CodonW 1.4.4's Fop (frequency of optimal codons, Ikemura 1981) and CBI (codon bias index, Bennetzen & Hall 1982) against a CodonW built-in optimal-codon set (fop_cbi_organism, default Escherichia coli = CodonW default) or a user list optimal_codons; fop/cbi are null when not requested or not calculable. Call for a one-shot codon-usage summary of a gene.")]
    public static CodonUsageStatistics codon_usage_statistics(
        [Description("Coding DNA sequence (frame 0).")] string sequence,
        [Description("NCBI translation table number defining the synonymous codon families and stop codons (default 1 = Standard; e.g. 2 = vertebrate mitochondrial, 11 = bacterial/plastid). Must be one of the NCBI gc.prt tables (1-6, 9-16, 21-33).")] int genetic_code = 1,
        [Description("Convention for GC1/GC2/GC3 and overall GC: EmbossCusp (default; all codons, stop codons included, EMBOSS cusp) or CodonW (stop codons of genetic_code excluded, CodonW 1.4.4 gc_out).")] CodonGcConvention gc_convention = CodonGcConvention.EmbossCusp,
        [Description("false (default): fop and cbi are null (output otherwise unchanged). true: also compute CodonW 1.4.4 Fop = optimal/total and CBI = (optimal - expected)/(total - expected), expected = sum over amino acids of n_aa * (optimal codons of the aa)/(codons of the aa), over the amino acids with more than one codon in genetic_code that have at least one optimal codon (stops, Met/Trp excluded); null when not calculable (CodonW prints 0.000).")] bool include_fop_cbi = false,
        [Description("CodonW built-in optimal-codon set for Fop/CBI (CodonW -f_type): EscherichiaColi (default, Ikemura 1985), BacillusSubtilis, DictyosteliumDiscoideum, AspergillusNidulans, SaccharomycesCerevisiae, DrosophilaMelanogaster, CaenorhabditisElegans, NeurosporaCrassa (case/spaces ignored; 'E. coli' style accepted). Ignored when optimal_codons is given or include_fop_cbi is false.")] string fop_cbi_organism = "EscherichiaColi",
        [Description("Omitted/null (default): use fop_cbi_organism. A list of optimal codons (DNA or RNA, any case) instead, like a CodonW -fop_file/-cbi_file; codons of stops or single-codon amino acids are ignored.")] string[]? optimal_codons = null)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new System.ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        var code = ResolveGeneticCode(genetic_code);
        if (!System.Enum.IsDefined(gc_convention))
            throw new System.ArgumentException(
                $"Unknown GC convention: {gc_convention}. Supported: EmbossCusp, CodonW.", nameof(gc_convention));

        if (!include_fop_cbi)
            return CodonUsageAnalyzer.GetStatistics(sequence, code, gc_convention);

        IReadOnlySet<string> optimal;
        if (optimal_codons is not null)
        {
            optimal = optimal_codons.ToHashSet();
        }
        else if (CodonWOptimalCodonSets.TryGetByName(fop_cbi_organism, out var set))
        {
            optimal = set.OptimalCodons;
        }
        else
        {
            throw new System.ArgumentException(
                $"Unknown Fop/CBI organism: {fop_cbi_organism}. Supported: EscherichiaColi, BacillusSubtilis, DictyosteliumDiscoideum, AspergillusNidulans, SaccharomycesCerevisiae, DrosophilaMelanogaster, CaenorhabditisElegans, NeurosporaCrassa.",
                nameof(fop_cbi_organism));
        }

        try
        {
            return CodonUsageAnalyzer.GetStatistics(sequence, code, gc_convention, optimal);
        }
        catch (System.ArgumentException ex) when (ex.ParamName == "optimalCodons")
        {
            throw new System.ArgumentException(ex.Message.Split(" (Parameter")[0], nameof(optimal_codons));
        }
    }

    // Shared description of the optional genetic_code parameter of the CodonOptimizer tools (F34 pattern).
    private const string GeneticCodeParameterDescription =
        "NCBI translation table number defining the encoded amino acids, synonymous codon families and stop codons (default 1 = Standard, identical to omitting it; e.g. 2 = vertebrate mitochondrial, 11 = bacterial/plastid). Must be one of the NCBI gc.prt tables (1-6, 9-16, 21-33).";

    // Maps an NCBI translation table number to a GeneticCode; an unknown number is reported
    // with the list of supported tables (GeneticCode.SupportedTableNumbers).
    private static GeneticCode ResolveGeneticCode(int genetic_code) =>
        GeneticCode.GetByTableNumberOrThrow(genetic_code, nameof(genetic_code));

    [McpServerTool(Name = "optimize_codons", Title = "MolTools — Optimize Codons for Expression", ReadOnly = true), Description("Optimizes a coding sequence for expression in a target organism using one of five strategies (MaximizeCAI, BalancedOptimization (default), HarmonizeExpression, MinimizeSecondary, AvoidRareCodons — also accepted under its historical misspelling AvoidRareCodeons; strategy names are matched case-insensitively). Internally trims to whole codons and converts T→U; stop codons and single-codon amino acids (Met/Trp) are left unchanged. Returns the original/optimized RNA, translated protein, original/optimized CAI, GC fractions, the number of changed codons, and each codon change. HarmonizeExpression matches the target codon-usage profile as closely as integer codon counts allow (DNA Chisel match_codon_usage) and is deterministic. The GC pass constrains the global GC fraction by default; set gc_window to constrain every sliding window of that many nucleotides instead (DNA Chisel EnforceGCContent(mini, maxi, window=...)). MinimizeSecondary = BalancedOptimization (codon selection + GC pass) followed by the sliding-window structure pass of reduce_secondary_structure, scored by secondary_scoring (PairCount heuristic by default; Turner2004Mfe for the Turner 2004 free-energy model, ~9 s per kilobase). The encoded protein and the synonymous-codon families follow the NCBI genetic code genetic_code (default 1 = Standard; DNA Chisel EnforceTranslation(genetic_table=...)), e.g. under table 2 AGA/AGG are stops and never used for Arg, UGA is Trp and AUA is Met. Pass avoid_sites to additionally clear restriction recognition sequences from the optimized sequence (hard constraint applied last, DNA Chisel AvoidPattern).")]
    public static OptimizationResultDto optimize_codons(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Target organism: a preset id (EColiK12 | Yeast | Human) or an inline custom table (organismName + codonFrequencies in RNA alphabet).")] CodonUsageTableInput target_organism,
        [Description("Optimization strategy.")] CodonOptimizer.OptimizationStrategy strategy = CodonOptimizer.OptimizationStrategy.BalancedOptimization,
        [Description("Lower bound for target GC fraction (default 0.40).")] double gc_target_min = 0.40,
        [Description("Upper bound for target GC fraction (default 0.60).")] double gc_target_max = 0.60,
        [Description("Frequency threshold below which a codon is considered rare (default 0.15).")] double rare_codon_threshold = 0.15,
        [Description("Window scoring model of the MinimizeSecondary structure pass: PairCount (default, canonical-pair fraction heuristic) or Turner2004Mfe (Turner 2004 nearest-neighbour free energy; thermodynamic but much slower). Ignored by the other strategies.")] CodonOptimizer.SecondaryStructureScoring secondary_scoring = CodonOptimizer.SecondaryStructureScoring.PairCount,
        [Description("Omitted/null (default): the GC pass constrains the GLOBAL GC fraction only (DNA Chisel EnforceGCContent(mini, maxi)). A positive window size in nucleotides instead constrains EVERY sliding window of that size to [gc_target_min, gc_target_max] (DNA Chisel EnforceGCContent(mini, maxi, window=...)), resolved by a deterministic greedy synonymous-codon repair. Used by BalancedOptimization and MinimizeSecondary.")] int? gc_window = null,
        [Description(GeneticCodeParameterDescription)] int genetic_code = 1,
        [Description("Omitted/null (default): no restriction-site avoidance. A list of recognition sequences (DNA or RNA, IUPAC codes allowed, same form as remove_restriction_sites) instead runs the site-removal pass on the optimized sequence AFTER codon selection, the GC pass and the structure pass, so no later pass can re-introduce a site (DNA Chisel composes CodonOptimize with AvoidPattern(EnzymeSitePattern(...)) + EnforceTranslation in one problem and resolves the hard constraints last). Both strands are cleared; each removal changes exactly one codon; an occurrence no single synonymous substitution can remove is left in place. Each site is at most 64 characters long (the scan costs sequence length x site length); a longer site, or one with a character that is not an IUPAC nucleotide code, is refused whatever the sequence is (an empty one included). The pass costs time linear in the sequence length, except that with a GC target every candidate is also scored against the whole sequence. For the strategies with a GC target (BalancedOptimization, MinimizeSecondary) the pass is GC-constraint-aware: it keeps gc_target_min/gc_target_max, and every gc_window window, satisfied whenever a site-removing synonymous codon can (preferring, in order, GC-satisfying candidates, then codon frequency, then the leftmost codon, then NCBI codon order); if none can, the least-breach candidate is used and the site is still removed (site avoidance stays the hard constraint). Reported CAI, GC and codon changes always describe the final sequence.")] string[]? avoid_sites = null)
    {
        if (string.IsNullOrEmpty(coding_sequence))
            throw new System.ArgumentException("Coding sequence cannot be null or empty.", nameof(coding_sequence));
        if (gc_target_min < 0 || gc_target_min > 1)
            throw new System.ArgumentException("GC target minimum must be a fraction in [0, 1].", nameof(gc_target_min));
        if (gc_target_max < 0 || gc_target_max > 1)
            throw new System.ArgumentException("GC target maximum must be a fraction in [0, 1].", nameof(gc_target_max));
        if (gc_target_min > gc_target_max)
            throw new System.ArgumentException("GC target minimum must not exceed the maximum.", nameof(gc_target_min));
        if (!System.Enum.IsDefined(secondary_scoring))
            throw new System.ArgumentException(
                $"Unknown secondary-structure scoring: {secondary_scoring}. Supported: PairCount, Turner2004Mfe.", nameof(secondary_scoring));
        if (gc_window is <= 0)
            throw new System.ArgumentException("GC window must be a positive number of nucleotides.", nameof(gc_window));
        var code = ResolveGeneticCode(genetic_code);

        var table = ResolveCodonUsageTable(target_organism, code);
        var result = CodonOptimizer.OptimizeSequence(
            coding_sequence, table, strategy, gc_target_min, gc_target_max, rare_codon_threshold, secondary_scoring, gc_window, code,
            avoid_sites);
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

    [McpServerTool(Name = "cai_from_organism_table", Title = "MolTools — CAI from Codon-Usage Table", ReadOnly = true), Description("Computes the Codon Adaptation Index (Sharp & Li 1987) for a coding sequence against an organism codon-usage FREQUENCY table (distinct from codon_adaptation_index, which takes a reference RSCU dictionary). CAI is the geometric mean of per-codon relative adaptiveness w = f(codon)/max f(synonymous). Stop codons and single-codon amino acids (Met/Trp) are excluded (Sharp & Li; CodonW); a codon with w < 0.0001 (absent from the table while a synonym is present) is scored with w = 0.01 (CodonW); amino acids without frequency data are skipped. Stops, single-codon families and synonymous families follow the NCBI genetic code genetic_code (default 1 = Standard; CodonW -code). Same core as codon_adaptation_index. Call when scoring how well a gene matches an organism's preferred codons.")]
    public static CaiResult cai_from_organism_table(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Target organism: preset id (EColiK12 | Yeast | Human) or inline custom table.")] CodonUsageTableInput target_organism,
        [Description(GeneticCodeParameterDescription)] int genetic_code = 1)
    {
        if (string.IsNullOrEmpty(coding_sequence))
            throw new System.ArgumentException("Coding sequence cannot be null or empty.", nameof(coding_sequence));
        var code = ResolveGeneticCode(genetic_code);

        var table = ResolveCodonUsageTable(target_organism, code);
        return new CaiResult(CodonOptimizer.CalculateCAI(coding_sequence, table, geneticCode: code));
    }

    [McpServerTool(Name = "remove_restriction_sites", Title = "MolTools — Remove Restriction Sites", ReadOnly = true), Description("Synonymously rewrites codons to eliminate the listed restriction recognition sequences from a coding sequence while preserving the encoded protein (RNA-alphabet output). Site strings may be DNA or RNA and may contain IUPAC ambiguity codes; both strands are cleared (a non-palindromic site such as BsaI GGTCTC is also removed where its reverse complement occurs); each removal changes exactly one codon, choosing the highest-frequency synonymous codon in the supplied table; sites with no synonymous alternative are left in place. Each site is at most 64 characters long (the scan costs sequence length x site length); a longer site, or one with a character that is not an IUPAC nucleotide code, is refused whatever the coding sequence is. The time is linear in the sequence length. Synonymous codons are those of the NCBI genetic code genetic_code (default 1 = Standard), so the protein is preserved under that code. Call to make a gene compatible with a cloning strategy.")]
    public static OptimizedSequenceResult remove_restriction_sites(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Restriction recognition sequences to eliminate (DNA or RNA, IUPAC codes allowed, each at most 64 characters; empty entries are skipped).")] string[] restriction_sites,
        [Description("Target organism: preset id or inline custom table (used for synonymous-codon lookup).")] CodonUsageTableInput target_organism,
        [Description(GeneticCodeParameterDescription)] int genetic_code = 1)
    {
        if (string.IsNullOrEmpty(coding_sequence))
            throw new System.ArgumentException("Coding sequence cannot be null or empty.", nameof(coding_sequence));
        if (restriction_sites is null || restriction_sites.Length == 0)
            throw new System.ArgumentException("At least one restriction site is required.", nameof(restriction_sites));
        var code = ResolveGeneticCode(genetic_code);

        var table = ResolveCodonUsageTable(target_organism, code);
        return new OptimizedSequenceResult(CodonOptimizer.RemoveRestrictionSites(coding_sequence, restriction_sites, table, code));
    }

    [McpServerTool(Name = "reduce_secondary_structure", Title = "MolTools — Reduce mRNA Secondary Structure", ReadOnly = true), Description("Greedy synonymous-codon swap inside a sliding window that lowers the window structure score, reducing mRNA secondary structure while preserving the protein. Two scoring models: PairCount (default) is a heuristic local self-complementarity score (fraction of canonical pairs incl. G\u00b7U wobble; no stacking, loop penalty or nesting); Turner2004Mfe scores every window with the Turner 2004 nearest-neighbour minimum free energy (the energy model of the analysis-server minimum_free_energy tool / ViennaRNA RNAfold) and maximises it, i.e. makes it less negative \u2014 thermodynamic but ~150x slower (~9 s per kilobase vs ~0.06 s; window folds are memoised and the fold count is capped). Output is upper-case RNA trimmed to whole codons; sequences shorter than window_size are returned unchanged apart from that normalisation. Synonymous swaps use the families of the NCBI genetic code genetic_code (default 1 = Standard). Call to relax strong secondary structure in a coding sequence.")]
    public static OptimizedSequenceResult reduce_secondary_structure(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Target organism: preset id or inline custom table.")] CodonUsageTableInput target_organism,
        [Description("Sliding-window size in nucleotides (default 40).")] int window_size = 40,
        [Description("Window scoring model: PairCount (default, canonical-pair fraction heuristic) or Turner2004Mfe (Turner 2004 nearest-neighbour free energy; thermodynamic but much slower).")] CodonOptimizer.SecondaryStructureScoring scoring = CodonOptimizer.SecondaryStructureScoring.PairCount,
        [Description(GeneticCodeParameterDescription)] int genetic_code = 1)
    {
        if (string.IsNullOrEmpty(coding_sequence))
            throw new System.ArgumentException("Coding sequence cannot be null or empty.", nameof(coding_sequence));
        if (window_size <= 0)
            throw new System.ArgumentException("Window size must be positive.", nameof(window_size));
        if (!System.Enum.IsDefined(scoring))
            throw new System.ArgumentException(
                $"Unknown secondary-structure scoring: {scoring}. Supported: PairCount, Turner2004Mfe.", nameof(scoring));
        var code = ResolveGeneticCode(genetic_code);

        var table = ResolveCodonUsageTable(target_organism, code);
        return new OptimizedSequenceResult(
            CodonOptimizer.ReduceSecondaryStructure(coding_sequence, table, window_size, scoring: scoring, geneticCode: code));
    }

    [McpServerTool(Name = "find_rare_codons", Title = "MolTools — Find Rare Codons", ReadOnly = true), Description("Reports every codon in a coding sequence whose frequency in the target organism's codon-usage table is below the threshold (default 0.15; strict <, or <= with inclusive_threshold = true as in Sherlocc), with its 0-based position, codon (RNA), amino acid (translated with the NCBI genetic code genetic_code, default 1 = Standard; * for stop - e.g. UGA is * under table 1 and W under table 2; the rare-codon selection itself is a per-codon table lookup and does not depend on the code) and frequency. Frame-0 complete triplets only; ambiguous triplets (N, R, Y, …) are skipped without shifting the frame. Call to locate translation-slowing rare codons before optimization.")]
    public static RareCodonsResult find_rare_codons(
        [Description("Coding sequence (DNA or RNA).")] string coding_sequence,
        [Description("Target organism: preset id or inline custom table.")] CodonUsageTableInput target_organism,
        [Description("Frequency threshold (default 0.15).")] double threshold = 0.15,
        [Description("false (default): frequency < threshold; true: frequency <= threshold (Sherlocc 'equal or lower than the threshold', Chartier et al. 2012).")] bool inclusive_threshold = false,
        [Description(GeneticCodeParameterDescription)] int genetic_code = 1)
    {
        if (string.IsNullOrEmpty(coding_sequence))
            throw new System.ArgumentException("Coding sequence cannot be null or empty.", nameof(coding_sequence));
        var code = ResolveGeneticCode(genetic_code);

        var table = ResolveCodonUsageTable(target_organism, code);
        var rare = CodonOptimizer.FindRareCodons(coding_sequence, table, threshold, inclusive_threshold, code)
            .Select(t => new RareCodon(t.Position, t.Codon, t.AminoAcid, t.Frequency))
            .ToList();
        return new RareCodonsResult(rare);
    }

    [McpServerTool(Name = "codon_minmax_profile", Title = "MolTools — %MinMax Codon-Usage Profile", ReadOnly = true), Description("%MinMax codon-usage profile (Clarke & Clark 2008, PLoS ONE 3:e3412; term-for-term as the Clark lab calculateMinMax reference code) over a sliding window of window_size codons (default 18). One value per window start (0-based codon index 0 .. codonCount - window_size): positive = %Max (run of predominantly common synonymous codons, +100 = all most-common), negative = %Min (run of rare codons, -100 = all rarest), 0 = average usage. Synonymous families (including the stop family) are those of the NCBI genetic code genetic_code (default 1 = Standard); single-codon families (Met/Trp under table 1) contribute 0; an ambiguous triplet occupies its window position but contributes nothing. The reference takes per-thousand usage (Kazusa /1000 column): with usage_basis PerThousandWhenAvailable (default) the EColiK12, Yeast and Human presets use their Kazusa 316407 / 4932 / 9606 per-thousand usage and reproduce the reference; an inline table is used as given (supply per-thousand values for reference-equal results); usage_basis CodonFrequencies uses the per-amino-acid fractions on any preset, which weight every residue equally. Empty when the sequence has fewer than window_size complete codons. Call to locate slowly-translated (rare-codon) regions along a gene.")]
    public static MinMaxProfileResult codon_minmax_profile(
        [Description("Coding sequence (DNA or RNA, frame 0).")] string coding_sequence,
        [Description("Reference organism: preset id or inline custom table (per-thousand frequencies reproduce the reference values exactly).")] CodonUsageTableInput target_organism,
        [Description("Sliding-window width in codons (default 18, Clarke & Clark 2008); must be >= 1.")] int window_size = 18,
        [Description("Usage values: PerThousandWhenAvailable (default; preset per-thousand usage when present, as the reference) or CodonFrequencies (always the per-amino-acid fractions / inline values).")] CodonOptimizer.MinMaxUsageBasis usage_basis = CodonOptimizer.MinMaxUsageBasis.PerThousandWhenAvailable,
        [Description(GeneticCodeParameterDescription)] int genetic_code = 1)
    {
        if (string.IsNullOrEmpty(coding_sequence))
            throw new System.ArgumentException("Coding sequence cannot be null or empty.", nameof(coding_sequence));
        if (!System.Enum.IsDefined(usage_basis))
            throw new System.ArgumentException(
                $"Unknown usage basis: {usage_basis}. Supported: PerThousandWhenAvailable, CodonFrequencies.", nameof(usage_basis));
        var code = ResolveGeneticCode(genetic_code);

        var table = ResolveCodonUsageTable(target_organism, code);
        var windows = CodonOptimizer.CalculateMinMaxProfile(coding_sequence, table, window_size, usage_basis, code)
            .Select(w => new MinMaxWindowDto(w.WindowStartCodon, w.PercentMinMax))
            .ToList();
        return new MinMaxProfileResult(windows);
    }

    [McpServerTool(Name = "find_rare_codon_clusters", Title = "MolTools — Find Rare-Codon Clusters", ReadOnly = true), Description("Rare-codon clusters by the Sherlocc rule (Chartier et al. 2012, Bioinformatics 28:1438): a window of window_size codons (default 7) is a cluster when it holds at least min_rare_codons (default 4) rare codons; overlapping cluster windows are merged into maximal regions reported once. A codon is rare when its frequency in the organism table is equal to or lower than threshold (default 0.15; Sherlocc sherlocc.pl `$avrg <= $avrg_lim`); inclusive_threshold = false gives the strict < comparison of find_rare_codons' default; ambiguous triplets are never rare. Returns 0-based inclusive codon indices (start_codon, end_codon) and the rare-codon count per region. Call to find runs of rare codons (translational pause sites) rather than isolated rare codons.")]
    public static RareCodonClustersResult find_rare_codon_clusters(
        [Description("Coding sequence (DNA or RNA, frame 0).")] string coding_sequence,
        [Description("Target organism: preset id or inline custom table.")] CodonUsageTableInput target_organism,
        [Description("Per-codon rare-frequency cutoff (default 0.15).")] double threshold = 0.15,
        [Description("Cluster window width in codons (default 7, Sherlocc); must be >= 1.")] int window_size = 7,
        [Description("Minimum rare codons in a window to call a cluster (default 4, Sherlocc); must be >= 1.")] int min_rare_codons = 4,
        [Description("true (default, Sherlocc): frequency <= threshold is a pause; false: frequency < threshold.")] bool inclusive_threshold = true)
    {
        if (string.IsNullOrEmpty(coding_sequence))
            throw new System.ArgumentException("Coding sequence cannot be null or empty.", nameof(coding_sequence));

        var table = ResolveCodonUsageTable(target_organism);
        var clusters = CodonOptimizer.FindRareCodonClusters(coding_sequence, table, threshold, window_size, min_rare_codons, inclusive_threshold)
            .Select(c => new RareCodonClusterDto(c.StartCodon, c.EndCodon, c.RareCount))
            .ToList();
        return new RareCodonClustersResult(clusters);
    }

    [McpServerTool(Name = "find_conserved_rare_codon_clusters", Title = "MolTools — Conserved Rare-Codon Clusters (Sherlocc)", ReadOnly = true), Description("Sherlocc 1.0 conserved rare-codon cluster analysis of a protein-family CODON ALIGNMENT (Chartier, Gaudreault & Najmanovich 2012, Bioinformatics 28:1438; sherlocc.pl stage 3, github.com/mtthchrtr/sherlocc commit b3b4d15). Input: equal-length aligned in-frame coding sequences (DNA or RNA, any case; length a multiple of 3; a codon of only '-'/'.' is a gap; a partially gapped, ambiguous or unlisted triplet has no usage value). Every cell gets its codon's usage from the organism table - by default the per-thousand usage (Sherlocc's CUTG /1000 column: the EColiK12, Yeast and Human presets carry it (Kazusa 316407 / 4932 / 9606), an inline table's codonFrequencies values are taken AS per-thousand values; use_codon_frequencies = true scores the per-amino-acid fractions / inline values instead, threshold then on that scale). For each column the valued cells of columns i-3..i+3 of ALL sequences are averaged; the column is representative only when it has >= min_coverage x sequences x averaging_window valued cells (default 0.75, Sherlocc's gap/conservation filter), and a representative column with average <= threshold is a pause. Every cluster_window-column window (default 7) holding >= min_pauses pauses (default 4) yields a cluster from its first to its last pause; a cluster starting inside the previous one is merged (abutting ones are not). Returns column averages (null = not representative), valued-cell counts, pause columns, clusters (0-based inclusive start/end columns, pause count, Sherlocc mean usage with non-representative columns counted as 0, 1-based middle point), total cluster length, the fraction of the alignment in clusters and the mean cluster usage. Give either target_organism (one table for all rows) or target_organisms_per_sequence (one table per row, Sherlocc scores each member with its own organism's usage). Independent of the genetic code (per-codon table lookup). No P-values (not computed by sherlocc.pl). Call to find rare-codon clusters conserved across a family of homologous genes.")]
    public static ConservedRareCodonClustersResult find_conserved_rare_codon_clusters(
        [Description("Aligned in-frame coding sequences of the family (equal length, multiple of 3; gaps '-' or '.').")] string[] aligned_coding_sequences,
        [Description("Sherlocc threshold: a column is a pause when its window average is <= threshold (per thousand by default; Sherlocc README examples 14-16).")] double threshold,
        [Description("One codon-usage table for every sequence: preset id (EColiK12 | Yeast | Human, all carry Kazusa per-thousand usage) or inline custom table whose codonFrequencies values are used as per-thousand usage. Omit when target_organisms_per_sequence is given.")] CodonUsageTableInput? target_organism = null,
        [Description("Optional: one codon-usage table per aligned sequence, same order (instead of target_organism).")] CodonUsageTableInput[]? target_organisms_per_sequence = null,
        [Description("false (default): score per-thousand usage (Sherlocc); true: score the tables' codonFrequencies (per-amino-acid fractions for presets) - SherloccOptions.UseCodonFrequencies.")] bool use_codon_frequencies = false,
        [Description("Centred averaging window in columns (default 7 = columns i-3..i+3, Sherlocc); positive odd number.")] int averaging_window = 7,
        [Description("Cluster window width in columns (default 7, Sherlocc); must be >= 1.")] int cluster_window = 7,
        [Description("Minimum pause columns in a cluster window (default 4, Sherlocc); must be >= 1.")] int min_pauses = 4,
        [Description("Representativeness filter: minimum fraction of valued cells per averaging window (default 0.75, Sherlocc); in [0, 1].")] double min_coverage = 0.75)
    {
        if (aligned_coding_sequences is null || aligned_coding_sequences.Length == 0)
            throw new System.ArgumentException("At least one aligned coding sequence is required.", nameof(aligned_coding_sequences));
        if ((target_organism is null) == (target_organisms_per_sequence is null))
            throw new System.ArgumentException(
                "Supply exactly one of target_organism or target_organisms_per_sequence.", nameof(target_organism));

        var options = new CodonOptimizer.SherloccOptions
        {
            AveragingWindow = averaging_window,
            ClusterWindow = cluster_window,
            MinPauses = min_pauses,
            MinCoverage = min_coverage,
            UseCodonFrequencies = use_codon_frequencies,
        };

        var result = target_organism is not null
            ? CodonOptimizer.FindConservedRareCodonClusters(aligned_coding_sequences, ResolveSherloccTable(target_organism), threshold, options)
            : CodonOptimizer.FindConservedRareCodonClusters(
                aligned_coding_sequences, target_organisms_per_sequence!.Select(ResolveSherloccTable).ToArray(), threshold, options);

        return new ConservedRareCodonClustersResult(
            result.SequenceCount,
            result.AlignmentLength,
            result.ColumnAverages,
            result.ColumnUsageCounts,
            result.PauseColumns,
            result.Clusters.Select(c => new ConservedRareCodonClusterDto(c.StartColumn, c.EndColumn, c.PauseCount, c.MeanUsage, c.MiddlePoint)).ToList(),
            result.TotalClusterLength,
            result.ClusterFraction,
            result.MeanClusterUsage);
    }

    // Sherlocc tables: presets as-is (EColiK12/Yeast/Human carry PerThousand); an inline table's values are
    // the caller's per-thousand usage (Sherlocc files/cut), so they also populate PerThousand
    // (RNA-normalised keys of CreateCodonUsageTable).
    private static CodonOptimizer.CodonUsageTable ResolveSherloccTable(CodonUsageTableInput input)
    {
        var table = ResolveCodonUsageTable(input);
        return string.IsNullOrEmpty(input.Preset) ? table with { PerThousand = table.CodonFrequencies } : table;
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

    [McpServerTool(Name = "build_codon_table", Title = "MolTools — Build Codon-Usage Table", ReadOnly = true), Description("Derives a per-organism CodonUsageTable from a reference coding sequence by computing per-amino-acid relative codon frequencies (RNA alphabet, U not T). All 64 codons are returned: a codon absent from the reference set is counted as 0.5 (Sharp & Li 1987, as in Biopython CodonAdaptationIndex), so the relative adaptiveness derived from the table matches Biopython exactly. Frequencies are normalised within the synonymous families (stop family included) of the NCBI genetic code genetic_code (default 1 = Standard), which also gives the returned codon-to-amino-acid mapping. Call when the user wants a custom codon-usage table built from their own reference gene(s).")]
    public static CodonUsageTableDto build_codon_table(
        [Description("Reference coding sequence (DNA or RNA).")] string reference_sequence,
        [Description("Organism name to attach to the resulting table.")] string organism_name,
        [Description(GeneticCodeParameterDescription)] int genetic_code = 1)
    {
        if (string.IsNullOrEmpty(reference_sequence))
            throw new System.ArgumentException("Reference sequence cannot be null or empty.", nameof(reference_sequence));
        if (string.IsNullOrWhiteSpace(organism_name))
            throw new System.ArgumentException("Organism name cannot be null or blank.", nameof(organism_name));
        var code = ResolveGeneticCode(genetic_code);

        var table = CodonOptimizer.CreateCodonTableFromSequence(reference_sequence, organism_name, code);
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
    private static CodonOptimizer.CodonUsageTable ResolveCodonUsageTable(CodonUsageTableInput input, GeneticCode? code = null)
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
            input.CodonFrequencies,
            code);
    }

    [McpServerTool(Name = "crispr_system_info", Title = "MolTools — CRISPR System Metadata", ReadOnly = true), Description("Returns the metadata record (name, PAM sequence, guide length, PAM placement relative to target, description, cutOffset) for a known CRISPR nuclease system. cutOffset is the number of nucleotides between the PAM and the cut on the PAM-bearing strand that design_guide_rnas uses to place the cut (3 for the Cas9 systems, 18 for Cas12a, 13 for CasX - the midpoint of the 12-14 nt non-target-strand cleavage Liu 2019 report). An unknown system value is rejected with an ArgumentException. Call when the user needs the PAM/guide parameters of SpCas9, SaCas9, Cas12a, CasX, etc.")]
    public static CrisprSystem crispr_system_info(
        [Description("CRISPR system: SpCas9 | SpCas9_NAG | SpCas9_NGA | SaCas9 | SaCas9_NNGRRN | Cas12a | AsCas12a | LbCas12a | CasX.")] CrisprSystemType system_type)
    {
        return CrisprDesigner.GetSystem(system_type);
    }

    [McpServerTool(Name = "find_pam_sites", Title = "MolTools — Find CRISPR PAM Sites", ReadOnly = true), Description("Finds all PAM matches (forward + reverse strand) for the chosen CRISPR system. PAM matching honours IUPAC codes (e.g. NGG, NNGRRT, TTTV). Each site reports the PAM, the adjacent guide/target window, its position and strand; sites whose target window falls outside the sequence are skipped. Coordinates (position, targetStart) are always 0-based forward-strand; the PAM and guide sequences are read on the protospacer strand (CRISPOR convention). Call to enumerate targetable protospacers in a sequence.")]
    public static PamSitesResult find_pam_sites(
        [Description("DNA sequence to scan: A, C, G, T or N, in either case. Any other character (U, an IUPAC code other than N, a space, a trailing newline, a FASTA header line) is rejected with the character and its index.")] string sequence,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9)
    {
        RequireNotEmpty(sequence, nameof(sequence), "Sequence");

        var sites = CrisprDesigner.FindPamSites(sequence, system_type).ToList();
        return new PamSitesResult(sites);
    }

    [McpServerTool(Name = "design_guide_rnas", Title = "MolTools — Design CRISPR Guide RNAs", ReadOnly = true), Description("Generates and scores guide-RNA candidates whose Cas9/Cas12a cut site falls inside the requested region, ranked best-first. The cut site follows the CRISPOR convention: 3 bp 5\u0027 of the PAM on the PAM-bearing strand for Cas9, after the 18th protospacer base for Cas12a, 13 nt after the PAM for CasX (midpoint of Liu 2019\u0027s 12-14 nt non-target-strand cleavage; CRISPOR has no CasX cut) - on both strands. Candidates scoring below parameters.MinScore, or with a GC outside parameters.minGcContent..maxGcContent, are filtered out (the GC window is an independent filter exactly as CHOPCHOP's filterGCmin/filterGCmax, whatever scoreGcContent says; set 0..100 for none). The composite score is CHOPCHOP's own penalty sum shown as clamp(100 - penalty, 0, 100): 1 point per self-complementary 4-bp stem (selfComplementaryStems, SCORE['FOLDING']) plus CHOPCHOP's flat SCORE['CRISPR_BAD_GC'] = 300 (so score 0) for a GC outside parameters.minGcContent..maxGcContent (default 40-70, edges inclusive); parameters.scoreGcContent (default true) = false drops the GC score term, as CHOPCHOP's own default does (its separate GC filter defaults to 0..100, i.e. none). TTTT, seed GC and common cloning sites are reported in issues only, with no score effect (as CRISPOR warns). parameters.minScore is not from any tool (CHOPCHOP/CRISPOR have no score threshold; their filters are the GC window and maxSelfComplementaryStems). 20-nt NGG guides with 4 nt of 5\u0027 and 3 nt of 3\u0027 flanking context additionally report context30Mer and onTargetScore (the published Doench 2016 Rule Set 2 / Azimuth on-target efficacy score, 0..1), and parameters.ranking = OnTargetRuleSet2 ranks by it; grafMotif flags the Graf 2019 TT-/GCC- inefficiency motifs. Region indices are 0-based; region_end is inclusive and must satisfy 0 <= region_start <= region_end < sequence.Length. Call to enumerate high-quality guides targeting a locus.")]
    public static GuideRnasResult design_guide_rnas(
        [Description("DNA sequence containing the target region.")] string sequence,
        [Description("0-based start of the target region.")] int region_start,
        [Description("0-based inclusive end of the target region.")] int region_end,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9,
        [Description("Optional guide-RNA design parameters (minGcContent, maxGcContent, minScore, avoidPolyT, checkSelfComplementarity, scoreGcContent (default true; false = CHOPCHOP's own default of no GC score term), ranking, maxSelfComplementaryStems, selfComplementarityBackboneRegions). maxSelfComplementaryStems (default null = no filter, CHOPCHOP filterSelfCompMax = -1) drops guides with more CHOPCHOP 4-bp self-complementary stems; selfComplementarityBackboneRegions (default null = guide-internal stems only) lists sgRNA scaffold regions, on the guide's strand, the stems may also pair against (CHOPCHOP --backbone). Defaults are used when null.")] GuideRnaParameters? parameters = null)
    {
        RequireNotEmpty(sequence, nameof(sequence), "Sequence");
        if (region_start < 0 || region_start >= sequence.Length)
            throw new System.ArgumentException("Region start must be within the sequence.", nameof(region_start));
        if (region_end < region_start || region_end >= sequence.Length)
            throw new System.ArgumentException("Region end must satisfy region_start <= region_end < sequence.Length.", nameof(region_end));

        var guides = CrisprDesigner
            .DesignGuideRnas(new DnaSequence(sequence), region_start, region_end, system_type, parameters)
            .ToList();
        return new GuideRnasResult(guides);
    }

    [McpServerTool(Name = "evaluate_guide_rna", Title = "MolTools — Evaluate Guide RNA", ReadOnly = true), Description("Scores a single guide RNA against on-target quality heuristics: overall GC%, seed-region GC%, polyT (Pol III terminator) presence, self-complementarity, and common-restriction-site presence; returns a 0..100 score plus an issues list. The composite score is CHOPCHOP's own penalty sum shown as clamp(100 - penalty, 0, 100): 1 point per CHOPCHOP 4-bp self-complementary stem (selfComplementaryStems, SCORE['FOLDING']) plus CHOPCHOP's flat SCORE['CRISPR_BAD_GC'] = 300 (so score 0) for a GC outside minGcContent..maxGcContent (default 40-70, edges inclusive); parameters.scoreGcContent = false reproduces CHOPCHOP's default (GC term off). TTTT, seed GC and common cloning sites are reported in issues only (no score effect, as in CRISPOR). selfComplementarityScore is selfComplementaryStems divided by the len-4 windows CHOPCHOP examines (in [0, 1]; it replaced the former unsourced pair fraction) and is reported only. avoidPolyT controls the TTTT issue; checkSelfComplementarity switches off the stem penalty (CHOPCHOP -SC). For NGG systems the Graf 2019 TT-/GCC- inefficiency motif is reported in grafMotif and in issues (a warning only, no deduction, as in CRISPOR). The published Doench 2016 Rule Set 2 score needs 30 nt of genomic context and is therefore only reported by design_guide_rnas or calculate_on_target_rule_set2. Position is -1 for ad-hoc evaluation. The guide must have the system's guide length (system_type: 20 nt SpCas9, 21 nt SaCas9, 23 nt Cas12a, 20 nt CasX); parameters are validated (0 <= minGcContent <= maxGcContent <= 100, maxSelfComplementaryStems null or >= 0). Call to QC one guide sequence.")]
    public static GuideRnaCandidate evaluate_guide_rna(
        [Description("Guide RNA sequence.")] string guide_sequence,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9,
        [Description("Optional guide-RNA design parameters.")] GuideRnaParameters? parameters = null)
    {
        RequireNotEmpty(guide_sequence, nameof(guide_sequence), "Guide sequence");

        return CrisprDesigner.EvaluateGuideRna(guide_sequence, system_type, parameters);
    }

    [McpServerTool(Name = "find_off_targets", Title = "MolTools — CRISPR Off-Target Scan", ReadOnly = true), Description("Naïve genome scan: enumerates all PAM sites in the genome and reports those whose target differs from the guide by 1..max_mismatches (the 0-mismatch on-target is excluded). Per-hit offTargetScore is THIS LIBRARY'S OWN HEURISTIC, NOT A PUBLISHED MODEL (5 per mismatch in the PAM-proximal 12-nt seed + 2 per other mismatch; Hsu 2013 motivates only the direction, the weights 5/2 and the 12-nt width are unsourced); for published models over the same hits use crispr_mit_specificity_score (MIT/Hsu 2013) or crispr_cfd_specificity_score (Doench 2016 CFD). O(genome × guide) — recommend genome ≤ ~1 Mb. Guide length must match the system's guide length. Call to list a guide's candidate off-target sites.")]
    public static OffTargetsResult find_off_targets(
        [Description("Guide RNA sequence (length must match the system's guide length).")] string guide_sequence,
        [Description("Genome / reference sequence to scan.")] string genome,
        [Description("Maximum allowed mismatches (range 0..5; default 3).")] int max_mismatches = 3,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9)
    {
        RequireNotEmpty(guide_sequence, nameof(guide_sequence), "Guide sequence");
        RequireNotEmpty(genome, nameof(genome), "Genome");
        if (max_mismatches < 0 || max_mismatches > 5)
            throw new System.ArgumentException("Maximum mismatches must be in the range 0..5.", nameof(max_mismatches));

        var hits = CrisprDesigner
            .FindOffTargets(guide_sequence, new DnaSequence(genome), max_mismatches, system_type)
            .ToList();
        return new OffTargetsResult(hits);
    }

    [McpServerTool(Name = "crispr_specificity_score", Title = "MolTools — CRISPR Guide Specificity Score (library heuristic)", ReadOnly = true), Description("THIS LIBRARY'S OWN HEURISTIC, NOT A PUBLISHED MODEL: aggregates off-target hits (≤4 mismatches) for a guide RNA against a genome as max(0, 100 − Σ(5 per seed-region mismatch + 2 per PAM-distal mismatch)), giving 0..100 (100 = no off-targets). Hsu 2013 motivates only the seed weighting's direction; the weights 5/2, the 12-nt seed width, the linear sum and the subtraction from 100 are unsourced, and the score saturates at 0. For the published models call crispr_mit_specificity_score (MIT/Hsu 2013 guide score) or crispr_cfd_specificity_score (Doench 2016 CFD) instead — both reproduce CRISPOR numerically. Call this one only as a coarse triage number. Guide length must match the system's guide length. (C# API: CrisprDesigner.CalculateHeuristicSpecificityScore; the shorter CalculateSpecificityScore is a backward-compatible alias for it.)")]
    public static SpecificityResult crispr_specificity_score(
        [Description("Guide RNA sequence (length must match the system's guide length).")] string guide_sequence,
        [Description("Genome / reference sequence to scan.")] string genome,
        [Description("CRISPR system (default SpCas9).")] CrisprSystemType system_type = CrisprSystemType.SpCas9)
    {
        RequireNotEmpty(guide_sequence, nameof(guide_sequence), "Guide sequence");
        RequireNotEmpty(genome, nameof(genome), "Genome");

        // Delegates to the primary, honestly-named library method; the TOOL name is unchanged on
        // purpose (renaming it would be a breaking catalog change). B08 WP2 / D13.
        return new SpecificityResult(
            CrisprDesigner.CalculateHeuristicSpecificityScore(guide_sequence, new DnaSequence(genome), system_type));
    }

    [McpServerTool(Name = "crispr_mit_specificity_score", Title = "MolTools — CRISPR Guide Specificity (MIT/Hsu 2013)", ReadOnly = true), Description("The PUBLISHED MIT / Hsu 2013 aggregate guide specificity score ('Guide score' of crispr.mit.edu, CRISPOR's calcMitGuideScore): scans the genome for off-targets with 1..max_mismatches mismatches, scores each with the Hsu position-weight model, and returns 100/(100 + Σ hitScores) × 100 in 0..100 (100 = no off-targets). Unrounded (the reference rounds to an integer). The per-hit model is picked by system exactly as the reference's annotateOfftargets does: a SaCas9 PAM (NNGRRT/NNNRRT) with its 21-nt guide is scored with the SaCas9-specific trained model of Tycko et al. 2018 (the reference's calcSaHitScore), every other PAM-after-target system with the Hsu model. Cas12a/Cpf1-class and CasX (PAM-before-target) systems are REFUSED: the reference reports no MIT off-target score for them. Set include_alternative_pams=true to also count off-targets sitting on the reference's alternative PAMs (its offtargetPams table: NGG → NAG, NGA; NNGRRT → NNGRRN), each NGG-guide alternative hit contributing only 0.2 × its hit score exactly as annotateOfftargets does ('alternative PAMs represent only ~10% of all cleavage events'; the rule is gated on the main PAM being NGG, so SaCas9 NNGRRN hits enter unweighted); the default false scans the one system's PAM only, which reports a guide as MORE specific than the reference. Like CRISPOR's annotateOfftargets, exactly ONE kind of perfect match is dropped as 'the on-target' and every further perfect match is scored as a full off-target (hit score 100): by default (CRISPOR's no-interval 'auto-ontarget' mode) the first exact match in scan order - forward-strand PAMs by ascending position, then reverse-strand PAMs - so a guide with 2 / 5 / 50 perfect sites scores 50 / 20 / 2 (CFD aggregate: 99.0 / 96.2 / 67.1); pass on_target_start/on_target_end (0-based, end exclusive, covering the protospacer+PAM) to name the on-target locus instead (exact matches inside it are dropped, every other one is scored). A guide with more than max_occurrences alignments (windows of the guide's length, either strand, any PAM, within max_mismatches substitutions; default 60000 = CRISPOR's MAXOCC) is a repeat: the reference suppresses it and the score is 0. Call for a publication-grade specificity number instead of crispr_specificity_score.")]
    public static SpecificityResult crispr_mit_specificity_score(
        [Description("Guide RNA sequence (length must match the system's guide length).")] string guide_sequence,
        [Description("Genome / reference sequence to scan.")] string genome,
        [Description("Maximum allowed mismatches (range 0..5; default 4, the reference's default).")] int max_mismatches = 4,
        [Description("CRISPR system (default SpCas9); must be PAM-after-target (Cas9-class).")] CrisprSystemType system_type = CrisprSystemType.SpCas9,
        [Description("Also scan the reference's alternative PAMs for this system (default false = main PAM only): SpCas9 adds NAG and NGA (downweighted by 0.2), SaCas9 adds NNGRRN (unweighted, as the reference).")] bool include_alternative_pams = false,
        [Description("Only with include_alternative_pams: minimum unweighted hit score (0..100) an alternative-PAM hit needs to count — the reference's ALTPAMMINSCORE / --minAltPamScore. Default 1.0, the reference's own default; 0 = no filtering.")] double min_alternative_pam_hit_score = 1.0,
        [Description("Optional on-target interval start (0-based, inclusive, forward-strand coordinates of the genome string; give both on_target_start and on_target_end or neither). Exact (0-mismatch) matches whose protospacer+PAM footprint lies inside [on_target_start, on_target_end) are skipped as the on-target - CRISPOR's annotateOfftargets skip - and every other perfect match is scored as an off-target (hit score 100). Default (omitted): the first exact match in scan order is the on-target (CRISPOR's auto-ontarget mode) and every further perfect match is scored.")] int? on_target_start = null,
        [Description("Optional on-target interval end (0-based, exclusive); see on_target_start.")] int? on_target_end = null,
        [Description("Genomic-occurrence limit of a repeat guide (default 60000, CRISPOR's MAXOCC): a guide with more alignments than this - windows of the guide's length, either strand, any PAM, within max_mismatches substitutions - is suppressed and scores 0.")] int max_occurrences = CrisprDesigner.DefaultMaxOccurrences)
    {
        RequireNotEmpty(guide_sequence, nameof(guide_sequence), "Guide sequence");
        RequireNotEmpty(genome, nameof(genome), "Genome");
        if (max_mismatches < 0 || max_mismatches > 5)
            throw new System.ArgumentException("Maximum mismatches must be in the range 0..5.", nameof(max_mismatches));
        if (max_occurrences < 1)
            throw new System.ArgumentException("Maximum occurrences must be at least 1.", nameof(max_occurrences));

        return new SpecificityResult(include_alternative_pams
            ? CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(
                guide_sequence, new DnaSequence(genome), max_mismatches, system_type, min_alternative_pam_hit_score, on_target_start, on_target_end, max_occurrences)
            : CrisprDesigner.CalculateMitSpecificityScore(
                guide_sequence, new DnaSequence(genome), max_mismatches, system_type, on_target_start, on_target_end, max_occurrences));
    }

    [McpServerTool(Name = "crispr_cfd_specificity_score", Title = "MolTools — CRISPR Guide Specificity (Doench 2016 CFD)", ReadOnly = true), Description("The PUBLISHED aggregate CFD (Cutting Frequency Determination, Doench 2016) guide specificity score, following CRISPOR's annotateOfftargets (`calcCfdScore` on each hit, summed and aggregated by calcMitGuideScore): scans the genome for off-targets with 1..max_mismatches mismatches, scores each with the CFD mismatch/PAM matrices using that site's own PAM, and aggregates the raw per-hit CFD values. By default only the system's own PAM is scanned, which reports a guide as MORE specific than the reference; set include_alternative_pams=true to also count off-targets on the reference's alternative PAMs (NGG → NAG, NGA) exactly as its wider search does (an alternative-PAM hit must pass the MIT-hit-score cutoff min_alternative_pam_hit_score = ALTPAMMINSCORE, default 1.0; no 0.2 factor on the CFD path — the PAM penalty is in the CFD PAM table). aggregation=MitStyle (CRISPOR's default) returns 100/(100 + Σ CFD) × 100, so 100 = no off-targets; aggregation=ParkinsonNormalized returns the doCfdFix variant 100 / Σ CFD, which is not bounded by 100 and errors when the CFD sum is 0. CFD is defined only for a 20-nt protospacer, so SaCas9 (21 nt), Cas12a-class and CasX are REFUSED: the reference does not score them either (per-hit cfdScore = -1 for a SaCas9 PAM, 0.0 for a Cpf1 PAM). Like CRISPOR's annotateOfftargets, exactly ONE kind of perfect match is dropped as 'the on-target' and every further perfect match is scored as a full off-target (hit score 100): by default (CRISPOR's no-interval 'auto-ontarget' mode) the first exact match in scan order - forward-strand PAMs by ascending position, then reverse-strand PAMs - so a guide with 2 / 5 / 50 perfect sites scores 50 / 20 / 2 (CFD aggregate: 99.0 / 96.2 / 67.1); pass on_target_start/on_target_end (0-based, end exclusive, covering the protospacer+PAM) to name the on-target locus instead (exact matches inside it are dropped, every other one is scored). A guide with more than max_occurrences alignments (windows of the guide's length, either strand, any PAM, within max_mismatches substitutions; default 60000 = CRISPOR's MAXOCC) is a repeat: the reference suppresses it and the score is 0. Call for a publication-grade specificity number instead of crispr_specificity_score.")]
    public static SpecificityResult crispr_cfd_specificity_score(
        [Description("Guide RNA sequence (must be 20 nt for the supported systems).")] string guide_sequence,
        [Description("Genome / reference sequence to scan.")] string genome,
        [Description("Maximum allowed mismatches (range 0..5; default 4, the reference's default).")] int max_mismatches = 4,
        [Description("CRISPR system (default SpCas9); must be a 20-nt PAM-after-target system.")] CrisprSystemType system_type = CrisprSystemType.SpCas9,
        [Description("Which published aggregation to use (default MitStyle, the reference's default).")] CfdAggregation aggregation = CfdAggregation.MitStyle,
        [Description("Also scan the reference's alternative PAMs for this system (default false = main PAM only): SpCas9 adds NAG and NGA; NGA adds NGG. Every hit is scored by its own PAM with the CFD matrices, no 0.2 factor (that factor is applied by the reference to the MIT score only).")] bool include_alternative_pams = false,
        [Description("Only with include_alternative_pams: minimum MIT/Hsu hit score (0..100) an alternative-PAM hit needs to count — the reference's ALTPAMMINSCORE / --minAltPamScore, tested on the MIT hit score even for the CFD aggregate. Default 1.0, the reference's own default; 0 = no filtering.")] double min_alternative_pam_hit_score = 1.0,
        [Description("Optional on-target interval start (0-based, inclusive, forward-strand coordinates of the genome string; give both on_target_start and on_target_end or neither). Exact (0-mismatch) matches whose protospacer+PAM footprint lies inside [on_target_start, on_target_end) are skipped as the on-target - CRISPOR's annotateOfftargets skip - and every other perfect match is scored as an off-target (hit score 100). Default (omitted): the first exact match in scan order is the on-target (CRISPOR's auto-ontarget mode) and every further perfect match is scored.")] int? on_target_start = null,
        [Description("Optional on-target interval end (0-based, exclusive); see on_target_start.")] int? on_target_end = null,
        [Description("Genomic-occurrence limit of a repeat guide (default 60000, CRISPOR's MAXOCC): a guide with more alignments than this - windows of the guide's length, either strand, any PAM, within max_mismatches substitutions - is suppressed and scores 0 under either aggregation.")] int max_occurrences = CrisprDesigner.DefaultMaxOccurrences)
    {
        RequireNotEmpty(guide_sequence, nameof(guide_sequence), "Guide sequence");
        RequireNotEmpty(genome, nameof(genome), "Genome");
        if (max_mismatches < 0 || max_mismatches > 5)
            throw new System.ArgumentException("Maximum mismatches must be in the range 0..5.", nameof(max_mismatches));
        if (max_occurrences < 1)
            throw new System.ArgumentException("Maximum occurrences must be at least 1.", nameof(max_occurrences));

        return new SpecificityResult(include_alternative_pams
            ? CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(
                guide_sequence, new DnaSequence(genome), max_mismatches, system_type, aggregation, min_alternative_pam_hit_score, on_target_start, on_target_end, max_occurrences)
            : CrisprDesigner.CalculateCfdSpecificityScore(
                guide_sequence, new DnaSequence(genome), max_mismatches, system_type, aggregation, on_target_start, on_target_end, max_occurrences));
    }

    [McpServerTool(Name = "calculate_on_target_doench2014", Title = "MolTools — On-Target Score (Doench 2014 Rule Set 1)", ReadOnly = true), Description("Doench et al. 2014 \"Rule Set 1\" on-target efficacy score for an SpCas9 guide, returned on a 0..100 scale (higher = predicted more active). This is the published logistic linear model: intercept 0.59763615 + a GC-count term over the 20-nt protospacer using the paper's TWO coefficients (-0.2026259 for GC <= 10, -0.1665878 for GC > 10, applied as abs(10 - GC) x coefficient) + 70 position-specific single- and di-nucleotide weights, through a logistic sigmoid. Source: Doench et al. 2014, Nat Biotechnol 32:1262 (PMID 25184501); coefficients taken verbatim from the reference implementation doenchScore.py shipped in CRISPOR (Haeussler et al. 2016, Genome Biol 17:148). Input is the model's 30-nt context: 4 nt upstream + 20 nt protospacer + 3 nt PAM (must be NGG) + 3 nt downstream, A/C/G/T only; the NGG requirement is this library's explicit guard (the model's own feature table presumes NGG). Returned UNROUNDED: CRISPOR's own table reports int(100 x sigmoid), i.e. the floor of this value. Call to rank guides by predicted cutting efficiency when the 30-nt context is known.")]
    public static OnTargetScoreResult calculate_on_target_doench2014(
        [Description("30-nt context: 4 nt upstream + 20 nt protospacer + 3 nt NGG PAM + 3 nt downstream.")] string context_30mer)
    {
        RequireNotEmpty(context_30mer, nameof(context_30mer), "Context 30-mer");

        return new OnTargetScoreResult(CrisprDesigner.CalculateOnTargetDoench2014(context_30mer));
    }

    [McpServerTool(Name = "calculate_on_target_rule_set2", Title = "MolTools — On-Target Score (Doench 2016 Rule Set 2 / Azimuth)", ReadOnly = true), Description("Doench et al. 2016 \"Rule Set 2\" / Azimuth on-target efficacy score for an SpCas9 guide, conventionally in 0..1 (higher = predicted more active) — the \"Doench '16\" efficiency score reported by CRISPOR. Rule Set 2 is a trained gradient-boosted-tree model, reproduced here from Microsoft Research's Azimuth model. Input is the model's 30-nt context: 4 nt upstream + 20 nt protospacer + 3 nt PAM + 3 nt downstream, A/C/G/T only. The PAM is not audited, as in CRISPOR (pam_audit=False): the model was trained on NGG, so a NAG / NGA context is scored but by extrapolation. Pass amino_acid_cut_position and percent_peptide together to use Azimuth's gene-context (full) model instead of the sequence-only one; percent_peptide = -1 is Azimuth's \"no gene position\" sentinel and is answered with the sequence-only model, as in the reference. Call to rank guides by the published on-target activity model.")]
    public static OnTargetScoreResult calculate_on_target_rule_set2(
        [Description("30-nt context: 4 nt upstream + 20 nt protospacer + 3 nt PAM (trained on NGG; not audited) + 3 nt downstream.")] string context_30mer,
        [Description("Optional amino-acid position of the cut site in the target protein, >= 0 (requires percent_peptide); enables the gene-context model. Ignored when percent_peptide = -1.")] int? amino_acid_cut_position = null,
        [Description("Optional cut position as a percentage 0..100 along the coding sequence, or -1 for Azimuth's \"no gene position\" sentinel (which falls back to the sequence-only model); requires amino_acid_cut_position.")] double? percent_peptide = null)
    {
        RequireNotEmpty(context_30mer, nameof(context_30mer), "Context 30-mer");
        if (amino_acid_cut_position.HasValue != percent_peptide.HasValue)
            throw new System.ArgumentException(
                "amino_acid_cut_position and percent_peptide must be supplied together (gene-context model) or both omitted.",
                nameof(amino_acid_cut_position));

        double score = amino_acid_cut_position.HasValue
            ? CrisprDesigner.CalculateOnTargetRuleSet2(context_30mer, amino_acid_cut_position.Value, percent_peptide!.Value)
            : CrisprDesigner.CalculateOnTargetRuleSet2(context_30mer);
        return new OnTargetScoreResult(score);
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
        [Description("Absorbance at 260 nm (A260; finite).")] double absorbance260,
        [Description("Extinction coefficient ε in M⁻¹·cm⁻¹ (finite, positive).")] double extinction_coefficient,
        [Description("Path length in cm (default 1.0; finite, positive).")] double path_length = 1.0)
    {
        if (!double.IsFinite(absorbance260))
            throw new System.ArgumentException("Absorbance must be finite.", nameof(absorbance260));
        if (!double.IsFinite(extinction_coefficient))
            throw new System.ArgumentException("Extinction coefficient must be finite.", nameof(extinction_coefficient));
        if (!double.IsFinite(path_length))
            throw new System.ArgumentException("Path length must be finite.", nameof(path_length));
        if (extinction_coefficient <= 0)
            throw new System.ArgumentException("Extinction coefficient must be positive.", nameof(extinction_coefficient));
        if (path_length <= 0)
            throw new System.ArgumentException("Path length must be positive.", nameof(path_length));

        return new ConcentrationResult(
            ProbeDesigner.CalculateConcentration(absorbance260, extinction_coefficient, path_length));
    }
}
