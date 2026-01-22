using System;
using System.Collections.Generic;
using System.Linq;

namespace SuffixTree.Genomics;

/// <summary>
/// Analyzes restriction enzyme sites in DNA sequences.
/// Supports common restriction enzymes with palindromic recognition sequences.
/// </summary>
public static class RestrictionAnalyzer
{
    #region Built-in Enzyme Database

    private static readonly Dictionary<string, RestrictionEnzyme> _enzymes = new(StringComparer.OrdinalIgnoreCase)
    {
        // Type II restriction enzymes (most common in molecular biology)
        // Format: Name, Recognition sequence, Cut positions (forward, reverse from 5')

        // 6-cutter enzymes
        ["EcoRI"] = new("EcoRI", "GAATTC", 1, 5, "Escherichia coli"),
        ["BamHI"] = new("BamHI", "GGATCC", 1, 5, "Bacillus amyloliquefaciens"),
        ["HindIII"] = new("HindIII", "AAGCTT", 1, 5, "Haemophilus influenzae"),
        ["XhoI"] = new("XhoI", "CTCGAG", 1, 5, "Xanthomonas holcicola"),
        ["SalI"] = new("SalI", "GTCGAC", 1, 5, "Streptomyces albus"),
        ["PstI"] = new("PstI", "CTGCAG", 5, 1, "Providencia stuartii"),
        ["SphI"] = new("SphI", "GCATGC", 5, 1, "Streptomyces phaeochromogenes"),
        ["KpnI"] = new("KpnI", "GGTACC", 5, 1, "Klebsiella pneumoniae"),
        ["SacI"] = new("SacI", "GAGCTC", 5, 1, "Streptomyces achromogenes"),
        ["XbaI"] = new("XbaI", "TCTAGA", 1, 5, "Xanthomonas badrii"),
        ["NcoI"] = new("NcoI", "CCATGG", 1, 5, "Nocardia corallina"),
        ["NdeI"] = new("NdeI", "CATATG", 2, 4, "Neisseria denitrificans"),
        ["NheI"] = new("NheI", "GCTAGC", 1, 5, "Neisseria mucosa"),
        ["SpeI"] = new("SpeI", "ACTAGT", 1, 5, "Sphaerotilus species"),
        ["AvrII"] = new("AvrII", "CCTAGG", 1, 5, "Anabaena variabilis"),
        ["ClaI"] = new("ClaI", "ATCGAT", 2, 4, "Caryophanon latum"),
        ["AgeI"] = new("AgeI", "ACCGGT", 1, 5, "Agrobacterium gelatinovorum"),

        // 8-cutter enzymes (rare cutters)
        ["NotI"] = new("NotI", "GCGGCCGC", 2, 6, "Nocardia otitidis-caviarum"),
        ["SfiI"] = new("SfiI", "GGCCNNNNNGGCC", 8, 5, "Streptomyces fimbriatus"), // Note: N = any base
        ["PacI"] = new("PacI", "TTAATTAA", 5, 3, "Pseudomonas alcaligenes"),
        ["AscI"] = new("AscI", "GGCGCGCC", 2, 6, "Arthrobacter species"),
        ["FseI"] = new("FseI", "GGCCGGCC", 6, 2, "Frankia species"),
        ["SwaI"] = new("SwaI", "ATTTAAAT", 4, 4, "Staphylococcus warneri"),

        // 4-cutter enzymes (frequent cutters)
        ["MspI"] = new("MspI", "CCGG", 1, 3, "Moraxella species"),
        ["HpaII"] = new("HpaII", "CCGG", 1, 3, "Haemophilus parainfluenzae"),
        ["TaqI"] = new("TaqI", "TCGA", 1, 3, "Thermus aquaticus"),
        ["AluI"] = new("AluI", "AGCT", 2, 2, "Arthrobacter luteus"), // Blunt
        ["RsaI"] = new("RsaI", "GTAC", 2, 2, "Rhodopseudomonas sphaeroides"), // Blunt
        ["HaeIII"] = new("HaeIII", "GGCC", 2, 2, "Haemophilus aegyptius"), // Blunt
        ["DpnI"] = new("DpnI", "GATC", 2, 2, "Diplococcus pneumoniae"), // Blunt, methylation-dependent
        ["Sau3AI"] = new("Sau3AI", "GATC", 0, 4, "Staphylococcus aureus"),
        ["MboI"] = new("MboI", "GATC", 0, 4, "Moraxella bovis"),

        // Common cloning enzymes
        ["EcoRV"] = new("EcoRV", "GATATC", 3, 3, "Escherichia coli"), // Blunt
        ["SmaI"] = new("SmaI", "CCCGGG", 3, 3, "Serratia marcescens"), // Blunt
        ["HincII"] = new("HincII", "GTYRAC", 3, 3, "Haemophilus influenzae"), // Degenerate, blunt
        ["ScaI"] = new("ScaI", "AGTACT", 3, 3, "Streptomyces caespitosus"), // Blunt
        ["StuI"] = new("StuI", "AGGCCT", 3, 3, "Streptomyces tubercidicus"), // Blunt
        ["BglII"] = new("BglII", "AGATCT", 1, 5, "Bacillus globigii"),
        ["ApaI"] = new("ApaI", "GGGCCC", 5, 1, "Acetobacter pasteurianus"),
    };

    /// <summary>
    /// Gets all available restriction enzymes.
    /// </summary>
    public static IReadOnlyDictionary<string, RestrictionEnzyme> Enzymes => _enzymes;

    /// <summary>
    /// Gets a restriction enzyme by name.
    /// </summary>
    public static RestrictionEnzyme? GetEnzyme(string name) =>
        _enzymes.TryGetValue(name, out var enzyme) ? enzyme : null;

    /// <summary>
    /// Gets all enzymes with a specific cut length.
    /// </summary>
    public static IEnumerable<RestrictionEnzyme> GetEnzymesByCutLength(int length) =>
        _enzymes.Values.Where(e => e.RecognitionSequence.Length == length);

    /// <summary>
    /// Gets all enzymes that produce blunt ends.
    /// </summary>
    public static IEnumerable<RestrictionEnzyme> GetBluntCutters() =>
        _enzymes.Values.Where(e => e.IsBluntEnd);

    /// <summary>
    /// Gets all enzymes that produce sticky ends.
    /// </summary>
    public static IEnumerable<RestrictionEnzyme> GetStickyCutters() =>
        _enzymes.Values.Where(e => !e.IsBluntEnd);

    #endregion

    #region Restriction Site Finding

    /// <summary>
    /// Finds all restriction sites for a specific enzyme in a sequence.
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze.</param>
    /// <param name="enzymeName">Name of the restriction enzyme.</param>
    /// <returns>Collection of restriction sites found.</returns>
    public static IEnumerable<RestrictionSite> FindSites(DnaSequence sequence, string enzymeName)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (string.IsNullOrEmpty(enzymeName))
            throw new ArgumentNullException(nameof(enzymeName));

        var enzyme = GetEnzyme(enzymeName) ??
            throw new ArgumentException($"Unknown enzyme: {enzymeName}", nameof(enzymeName));

        return FindSitesCore(sequence.Sequence, enzyme);
    }

    /// <summary>
    /// Finds all restriction sites for a specific enzyme in a raw sequence string.
    /// </summary>
    public static IEnumerable<RestrictionSite> FindSites(string sequence, string enzymeName)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        var enzyme = GetEnzyme(enzymeName) ??
            throw new ArgumentException($"Unknown enzyme: {enzymeName}", nameof(enzymeName));

        foreach (var site in FindSitesCore(sequence.ToUpperInvariant(), enzyme))
            yield return site;
    }

    /// <summary>
    /// Finds all restriction sites for a custom enzyme.
    /// </summary>
    public static IEnumerable<RestrictionSite> FindSites(DnaSequence sequence, RestrictionEnzyme enzyme)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(enzyme);

        return FindSitesCore(sequence.Sequence, enzyme);
    }

    private static IEnumerable<RestrictionSite> FindSitesCore(string seq, RestrictionEnzyme enzyme)
    {
        string pattern = enzyme.RecognitionSequence;

        for (int i = 0; i <= seq.Length - pattern.Length; i++)
        {
            if (MatchesPattern(seq, i, pattern))
            {
                int cutPosition = i + enzyme.CutPositionForward;

                yield return new RestrictionSite(
                    Position: i,
                    Enzyme: enzyme,
                    IsForwardStrand: true,
                    CutPosition: cutPosition,
                    RecognizedSequence: seq.Substring(i, pattern.Length));
            }
        }

        // Check reverse strand (complement sequence)
        string revComp = DnaSequence.GetReverseComplementString(seq);
        for (int i = 0; i <= revComp.Length - pattern.Length; i++)
        {
            if (MatchesPattern(revComp, i, pattern))
            {
                // Convert back to forward strand coordinates
                int forwardPos = seq.Length - i - pattern.Length;
                int cutPosition = forwardPos + enzyme.CutPositionReverse;

                yield return new RestrictionSite(
                    Position: forwardPos,
                    Enzyme: enzyme,
                    IsForwardStrand: false,
                    CutPosition: cutPosition,
                    RecognizedSequence: seq.Substring(forwardPos, pattern.Length));
            }
        }
    }

    /// <summary>
    /// Finds all restriction sites for multiple enzymes.
    /// </summary>
    public static IEnumerable<RestrictionSite> FindSites(DnaSequence sequence, params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);

        foreach (var name in enzymeNames)
        {
            foreach (var site in FindSites(sequence, name))
                yield return site;
        }
    }

    /// <summary>
    /// Finds all sites for all known enzymes (comprehensive analysis).
    /// </summary>
    public static IEnumerable<RestrictionSite> FindAllSites(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);

        foreach (var enzyme in _enzymes.Values)
        {
            foreach (var site in FindSitesCore(sequence.Sequence, enzyme))
                yield return site;
        }
    }

    private static bool MatchesPattern(string seq, int position, string pattern)
    {
        for (int i = 0; i < pattern.Length; i++)
        {
            char seqChar = seq[position + i];
            char patChar = pattern[i];

            if (!MatchesIupac(seqChar, patChar))
                return false;
        }
        return true;
    }

    private static bool MatchesIupac(char nucleotide, char iupacCode) => iupacCode switch
    {
        'A' => nucleotide == 'A',
        'C' => nucleotide == 'C',
        'G' => nucleotide == 'G',
        'T' => nucleotide == 'T',
        'N' => nucleotide is 'A' or 'C' or 'G' or 'T',
        'R' => nucleotide is 'A' or 'G',
        'Y' => nucleotide is 'C' or 'T',
        'S' => nucleotide is 'G' or 'C',
        'W' => nucleotide is 'A' or 'T',
        'K' => nucleotide is 'G' or 'T',
        'M' => nucleotide is 'A' or 'C',
        'B' => nucleotide is 'C' or 'G' or 'T',
        'D' => nucleotide is 'A' or 'G' or 'T',
        'H' => nucleotide is 'A' or 'C' or 'T',
        'V' => nucleotide is 'A' or 'C' or 'G',
        _ => nucleotide == iupacCode
    };

    #endregion

    #region Digest Simulation

    /// <summary>
    /// Simulates a restriction digest and returns the resulting fragments.
    /// </summary>
    /// <param name="sequence">DNA sequence to digest.</param>
    /// <param name="enzymeNames">Names of restriction enzymes to use.</param>
    /// <returns>Collection of DNA fragments after digestion.</returns>
    public static IEnumerable<DigestFragment> Digest(DnaSequence sequence, params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (enzymeNames == null || enzymeNames.Length == 0)
            throw new ArgumentException("At least one enzyme is required", nameof(enzymeNames));

        // Find all cut sites
        var cutPositions = new SortedSet<int>();
        var sitesByPosition = new Dictionary<int, RestrictionSite>();

        foreach (var name in enzymeNames)
        {
            foreach (var site in FindSites(sequence, name))
            {
                if (site.IsForwardStrand) // Only count forward strand to avoid double counting
                {
                    cutPositions.Add(site.CutPosition);
                    sitesByPosition[site.CutPosition] = site;
                }
            }
        }

        if (cutPositions.Count == 0)
        {
            // No cuts - return entire sequence as one fragment
            yield return new DigestFragment(
                Sequence: sequence.Sequence,
                StartPosition: 0,
                Length: sequence.Length,
                LeftEnzyme: null,
                RightEnzyme: null,
                FragmentNumber: 1);
            yield break;
        }

        // Generate fragments
        var cuts = new List<int> { 0 };
        cuts.AddRange(cutPositions);
        cuts.Add(sequence.Length);

        for (int i = 0; i < cuts.Count - 1; i++)
        {
            int start = cuts[i];
            int end = cuts[i + 1];
            int length = end - start;

            if (length > 0)
            {
                string fragmentSeq = sequence.Sequence.Substring(start, length);

                yield return new DigestFragment(
                    Sequence: fragmentSeq,
                    StartPosition: start,
                    Length: length,
                    LeftEnzyme: i == 0 ? null : sitesByPosition.GetValueOrDefault(start)?.Enzyme.Name,
                    RightEnzyme: i == cuts.Count - 2 ? null : sitesByPosition.GetValueOrDefault(end)?.Enzyme.Name,
                    FragmentNumber: i + 1);
            }
        }
    }

    /// <summary>
    /// Gets a summary of digestion with fragment sizes sorted.
    /// </summary>
    public static DigestSummary GetDigestSummary(DnaSequence sequence, params string[] enzymeNames)
    {
        var fragments = Digest(sequence, enzymeNames).ToList();
        var sizes = fragments.Select(f => f.Length).OrderByDescending(x => x).ToList();

        return new DigestSummary(
            TotalFragments: fragments.Count,
            FragmentSizes: sizes,
            LargestFragment: sizes.FirstOrDefault(),
            SmallestFragment: sizes.LastOrDefault(),
            AverageFragmentSize: sizes.Count > 0 ? sizes.Average() : 0,
            EnzymesUsed: enzymeNames.ToList());
    }

    #endregion

    #region Restriction Map

    /// <summary>
    /// Creates a restriction map showing all enzyme sites.
    /// </summary>
    public static RestrictionMap CreateMap(DnaSequence sequence, params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);

        var sites = enzymeNames.Length > 0
            ? FindSites(sequence, enzymeNames).ToList()
            : FindAllSites(sequence).ToList();

        var byEnzyme = sites
            .GroupBy(s => s.Enzyme.Name)
            .ToDictionary(g => g.Key, g => g.Select(s => s.Position).OrderBy(p => p).ToList());

        var uniqueCutters = byEnzyme.Where(kv => kv.Value.Count == 1).Select(kv => kv.Key).ToList();
        var nonCutters = enzymeNames.Where(e => !byEnzyme.ContainsKey(e)).ToList();

        return new RestrictionMap(
            SequenceLength: sequence.Length,
            Sites: sites,
            SitesByEnzyme: byEnzyme,
            TotalSites: sites.Count(s => s.IsForwardStrand), // Count only forward strand
            UniqueCutters: uniqueCutters,
            NonCutters: nonCutters);
    }

    #endregion

    #region Compatibility Analysis

    /// <summary>
    /// Finds enzymes that produce compatible (ligatable) ends.
    /// </summary>
    public static IEnumerable<(string Enzyme1, string Enzyme2, string CompatibleEnd)> FindCompatibleEnzymes()
    {
        var enzymes = _enzymes.Values.ToList();

        for (int i = 0; i < enzymes.Count; i++)
        {
            for (int j = i + 1; j < enzymes.Count; j++)
            {
                var e1 = enzymes[i];
                var e2 = enzymes[j];

                if (AreCompatible(e1, e2, out string? compatibleEnd))
                {
                    yield return (e1.Name, e2.Name, compatibleEnd!);
                }
            }
        }
    }

    /// <summary>
    /// Checks if two enzymes produce compatible sticky ends.
    /// </summary>
    public static bool AreCompatible(string enzyme1Name, string enzyme2Name)
    {
        var e1 = GetEnzyme(enzyme1Name);
        var e2 = GetEnzyme(enzyme2Name);

        if (e1 == null || e2 == null) return false;

        return AreCompatible(e1, e2, out _);
    }

    private static bool AreCompatible(RestrictionEnzyme e1, RestrictionEnzyme e2, out string? compatibleEnd)
    {
        compatibleEnd = null;

        // Blunt ends are always compatible with other blunt ends
        if (e1.IsBluntEnd && e2.IsBluntEnd)
        {
            compatibleEnd = "blunt";
            return true;
        }

        // For sticky ends, check if overhangs are compatible
        string overhang1 = GetOverhang(e1);
        string overhang2 = GetOverhang(e2);

        if (overhang1 == overhang2 && !string.IsNullOrEmpty(overhang1))
        {
            compatibleEnd = overhang1;
            return true;
        }

        return false;
    }

    private static string GetOverhang(RestrictionEnzyme enzyme)
    {
        if (enzyme.IsBluntEnd) return "";

        int cutDiff = enzyme.CutPositionForward - enzyme.CutPositionReverse;

        // This is a simplified calculation - in reality would need actual sequence
        return enzyme.OverhangType switch
        {
            OverhangType.FivePrime => enzyme.RecognitionSequence.Substring(
                enzyme.CutPositionForward,
                Math.Abs(cutDiff)),
            OverhangType.ThreePrime => enzyme.RecognitionSequence.Substring(
                enzyme.CutPositionReverse,
                Math.Abs(cutDiff)),
            _ => ""
        };
    }

    #endregion
}

/// <summary>
/// Represents a restriction enzyme with its properties.
/// </summary>
public sealed record RestrictionEnzyme(
    string Name,
    string RecognitionSequence,
    int CutPositionForward,
    int CutPositionReverse,
    string Organism)
{
    /// <summary>
    /// Gets whether this enzyme produces blunt ends.
    /// </summary>
    public bool IsBluntEnd => CutPositionForward == CutPositionReverse;

    /// <summary>
    /// Gets the type of overhang produced.
    /// </summary>
    public OverhangType OverhangType => CutPositionForward == CutPositionReverse
        ? OverhangType.Blunt
        : CutPositionForward < CutPositionReverse
            ? OverhangType.FivePrime
            : OverhangType.ThreePrime;

    /// <summary>
    /// Gets the length of the recognition sequence.
    /// </summary>
    public int RecognitionLength => RecognitionSequence.Length;
}

/// <summary>
/// Type of overhang produced by restriction enzyme.
/// </summary>
public enum OverhangType
{
    /// <summary>Blunt end (no overhang)</summary>
    Blunt,
    /// <summary>5' overhang (sticky end)</summary>
    FivePrime,
    /// <summary>3' overhang (sticky end)</summary>
    ThreePrime
}

/// <summary>
/// Represents a restriction site in a sequence.
/// </summary>
public sealed record RestrictionSite(
    int Position,
    RestrictionEnzyme Enzyme,
    bool IsForwardStrand,
    int CutPosition,
    string RecognizedSequence);

/// <summary>
/// Represents a DNA fragment after restriction digestion.
/// </summary>
public sealed record DigestFragment(
    string Sequence,
    int StartPosition,
    int Length,
    string? LeftEnzyme,
    string? RightEnzyme,
    int FragmentNumber);

/// <summary>
/// Summary of a restriction digest.
/// </summary>
public sealed record DigestSummary(
    int TotalFragments,
    IReadOnlyList<int> FragmentSizes,
    int LargestFragment,
    int SmallestFragment,
    double AverageFragmentSize,
    IReadOnlyList<string> EnzymesUsed);

/// <summary>
/// A restriction map of a sequence.
/// </summary>
public sealed record RestrictionMap(
    int SequenceLength,
    IReadOnlyList<RestrictionSite> Sites,
    IReadOnlyDictionary<string, List<int>> SitesByEnzyme,
    int TotalSites,
    IReadOnlyList<string> UniqueCutters,
    IReadOnlyList<string> NonCutters);
