namespace Seqeron.Genomics.MolTools;

/// <summary>
/// A DNA methyltransferase recognition site, in EMBOSS <c>Emethylsites.dat</c> notation: each
/// methylated base is preceded by a lower-case <c>m</c> (<c>GmATC</c> = Dam, N6-methyladenine at
/// the A). The site is a plain recognition pattern otherwise (IUPAC codes allowed).
/// </summary>
/// <remarks>
/// This is the *methylase's* site, not an enzyme's sensitivity: a restriction enzyme is blocked
/// when its recognition window contains a methylated base (see
/// <see cref="RestrictionAnalyzer.FindSites(DnaSequence, string, MoleculeTopology, IReadOnlyList{Methylase})"/>).
/// EMBOSS upper-cases the file and reads every <c>M</c> as the marker; here only a lower-case
/// <c>m</c> is the marker, so the IUPAC code <c>M</c> (A or C) stays usable in a custom site.
/// </remarks>
/// <param name="Name">Methylase name (informational — EMBOSS ignores it too).</param>
/// <param name="Site">Site with <c>m</c> before every methylated base, e.g. <c>mCmCAGG</c>.</param>
public sealed record Methylase(string Name, string Site);

/// <summary>
/// Built-in methylase sets for the opt-in methylation-aware site finding / digest of
/// <see cref="RestrictionAnalyzer"/> (<c>methylases</c> parameter). Default behaviour (no set) is
/// unchanged.
/// </summary>
/// <remarks>
/// <see cref="Emboss"/> is the data of EMBOSS <c>restrict -methylation -mfile Emethylsites.dat</c>:
/// https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/data/Emethylsites.dat
/// (11 lines, fetched 2026-10-21; byte-identical to <c>/usr/share/EMBOSS/data/Emethylsites.dat</c> of
/// Ubuntu emboss 6.6.0+dfsg-12ubuntu2). <see cref="Dam"/> and <see cref="Dcm"/> are its subsets.
/// </remarks>
public static class MethylaseSets
{
    /// <summary>Dam: <c>GmATC</c>.</summary>
    public static IReadOnlyList<Methylase> Dam { get; } = [new("dam", "GmATC")];

    /// <summary>Dcm: <c>mCmCAGG</c> and <c>mCmCTGG</c>.</summary>
    public static IReadOnlyList<Methylase> Dcm { get; } = [new("dcm", "mCmCAGG"), new("dcm", "mCmCTGG")];

    /// <summary>Dam and Dcm.</summary>
    public static IReadOnlyList<Methylase> DamDcm { get; } = [.. Dam, .. Dcm];

    /// <summary>All five entries of EMBOSS <c>Emethylsites.dat</c> (Dam, Dcm x2, EcoKI x2), in file order.</summary>
    public static IReadOnlyList<Methylase> Emboss { get; } =
    [
        new("dam", "GmATC"),
        new("dcm", "mCmCAGG"),
        new("dcm", "mCmCTGG"),
        new("MEcoKI", "mAmACNNNNNNGTGC"),
        new("MEcoKI", "GCmACNNNNNNGTT"),
    ];

    /// <summary>
    /// Resolves a set name: <c>none</c>, <c>dam</c>, <c>dcm</c>, <c>dam+dcm</c> or <c>emboss</c>
    /// (case-insensitive); <c>none</c> and <see langword="null"/>/empty give <see langword="null"/>
    /// (methylation off).
    /// </summary>
    /// <exception cref="ArgumentException">Unknown set name.</exception>
    public static IReadOnlyList<Methylase>? FromName(string? name)
    {
        switch (name?.Trim().ToLowerInvariant())
        {
            case null or "" or "none": return null;
            case "dam": return Dam;
            case "dcm": return Dcm;
            case "dam+dcm": return DamDcm;
            case "emboss": return Emboss;
            default:
                throw new ArgumentException(
                    $"Unknown methylation set '{name}'; expected none, dam, dcm, dam+dcm or emboss.", nameof(name));
        }
    }
}

/// <summary>
/// EMBOSS <c>restrict -methylation</c> semantics (EMBOSS 6.6 <c>nucleus/embpat.c</c>:
/// <c>patRestrictReadMethyl</c>, <c>patRestrictMethylMod</c>, <c>embPatRestrictMatch</c>,
/// <c>embPatRestrictScan</c>): every methylated base found by a methylase site is replaced by
/// <c>N</c> — the base and its partner on the other strand — and an enzyme window containing such
/// a base is never a site.
/// </summary>
internal static class MethylationMask
{
    private readonly record struct Parsed(char[] Codes, bool[] Methylated);

    private static Parsed Parse(Methylase m, string paramName)
    {
        ArgumentNullException.ThrowIfNull(m);
        if (string.IsNullOrEmpty(m.Site))
            throw new ArgumentException($"Methylase '{m.Name}' has an empty site.", paramName);

        var codes = new List<char>();
        var marks = new List<bool>();
        bool pendingMark = false;
        bool sawMark = false;
        foreach (char raw in m.Site)
        {
            if (raw == 'm')
            {
                if (pendingMark)
                    throw new ArgumentException($"Methylase '{m.Name}' site '{m.Site}': two consecutive 'm' markers.", paramName);
                pendingMark = true;
                sawMark = true;
                continue;
            }

            char c = char.IsAscii(raw) ? char.ToUpperInvariant(raw) : raw;
            if (!IupacHelper.IsNucleotideCode(c))
                throw new ArgumentException(
                    $"Methylase '{m.Name}' site '{m.Site}' contains '{raw}', which is neither an IUPAC nucleotide code nor the 'm' marker.",
                    paramName);
            codes.Add(c);
            // EMBOSS writes N over the methylated bases AND over every N of the site (its Replace string
            // keeps 'N' for a pattern N): EcoKI's six-base spacer is masked too (verified with a
            // one-window "NN" probe enzyme under `restrict -methylation`).
            marks.Add(pendingMark || c == 'N');
            pendingMark = false;
        }

        if (pendingMark)
            throw new ArgumentException($"Methylase '{m.Name}' site '{m.Site}': 'm' marker without a following base.", paramName);
        if (!sawMark)
            throw new ArgumentException($"Methylase '{m.Name}' site '{m.Site}' marks no methylated base (no 'm').", paramName);

        return new Parsed([.. codes], [.. marks]);
    }

    /// <summary>
    /// The per-forward-coordinate mask (<see langword="true"/> = the base pair is methylated on at
    /// least one strand), or <see langword="null"/> when <paramref name="methylases"/> is null/empty.
    /// </summary>
    /// <remarks>
    /// Masked bases are the <c>m</c>-marked ones plus every <c>N</c> of the site (EMBOSS does the
    /// same). Same order of work as EMBOSS: every entry (last file entry first) is applied to the forward
    /// strand, then every entry to the reverse strand, each entry seeing the <c>N</c>s already written (so for palindromic
    /// <c>GmATC</c> only one A of a GATC is marked, as in EMBOSS). A circular molecule also finds
    /// methylase sites spanning the origin (EMBOSS masks the linear string before doubling it and
    /// does not).
    /// </remarks>
    internal static bool[]? Compute(string seq, IReadOnlyList<Methylase>? methylases, bool circular)
    {
        if (methylases is null || methylases.Count == 0)
            return null;

        // EMBOSS builds its list with ajListPush (adds at the FRONT) and pops from the front, so every
        // pass visits the entries in the REVERSE of file order (EcoKI, EcoKI, dcm, dcm, dam).
        var parsed = methylases.Select(m => Parse(m, nameof(methylases))).Reverse().ToArray();
        int n = seq.Length;
        var mask = new bool[n];

        Pass(seq, reversed: false);
        Pass(DnaSequence.GetReverseComplementString(seq), reversed: true);
        return mask;

        // Forward-strand coordinate of index k of the (possibly reverse-complemented) strand.
        int Forward(int k, bool reversed)
        {
            if (reversed)
                return n - 1 - k;
            return k;
        }

        void Pass(string s, bool reversed)
        {
            foreach (var (codes, methylated) in parsed)
            {
                int len = codes.Length;
                if (len > n)
                    continue;

                // Collect first, mark after — EMBOSS finds all hits of one entry before writing the Ns.
                var hits = new List<int>();
                int starts = circular ? n : n - len + 1;
                for (int i = 0; i < starts; i++)
                {
                    bool ok = true;
                    for (int j = 0; j < len && ok; j++)
                    {
                        int k = (i + j) % n;
                        char c = mask[Forward(k, reversed)] ? 'N' : s[k];
                        ok = SequenceBaseMatchesSiteCode(c, codes[j]);
                    }

                    if (ok)
                        hits.Add(i);
                }

                foreach (int i in hits)
                {
                    for (int j = 0; j < len; j++)
                    {
                        if (!methylated[j])
                            continue;
                        int k = (i + j) % n;
                        mask[Forward(k, reversed)] = true;
                    }
                }
            }
        }
    }

    /// <summary>Whether the window of <paramref name="len"/> bases at forward <paramref name="start"/>
    /// (wrapping modulo the mask length) contains a methylated base.</summary>
    internal static bool IsBlocked(bool[]? mask, int start, int len)
    {
        if (mask is null)
            return false;
        for (int j = 0; j < len; j++)
        {
            if (mask[(start + j) % mask.Length])
                return true;
        }

        return false;
    }

    private static bool SequenceBaseMatchesSiteCode(char sequenceBase, char siteCode) =>
        siteCode == 'N'
            ? IupacHelper.IsNucleotideCode(sequenceBase)
            : IupacHelper.MatchesIupac(sequenceBase, siteCode);
}

public static partial class RestrictionAnalyzer
{
    // Opt-in methylation-aware overloads (review 2026-09 B08 round 2 E2). A methylase set blocks every
    // enzyme window that contains a base the set methylates, on either strand — the rule of EMBOSS
    // `restrict -methylation`, NOT an enzyme-specific sensitivity table (so BamHI is blocked by Dam
    // here, as in EMBOSS, although it is insensitive in the bench literature). A null or empty set is
    // exactly the overloads without the parameter.

    /// <summary>
    /// <see cref="FindSites(DnaSequence, string, MoleculeTopology)"/> with methylated bases
    /// blocking recognition (EMBOSS <c>restrict -methylation</c> semantics).
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze.</param>
    /// <param name="enzymeName">Name of the restriction enzyme (case-insensitive).</param>
    /// <param name="topology">Whether the molecule is linear or circular.</param>
    /// <param name="methylases">Methylase sites (<see cref="MethylaseSets"/>); null or empty = methylation off.</param>
    /// <remarks>
    /// Every base found inside a methylase site at its <c>m</c> mark — and its partner on the other
    /// strand — is methylated; a window of the enzyme's recognition sequence that contains one
    /// methylated base is not reported (EMBOSS <c>embPatRestrictScan</c>: such a base reads as
    /// <c>N</c> and breaks the match). <see cref="RestrictionSite.RecognizedSequence"/> still shows the
    /// real bases. Methylase sites are searched on both strands; on a circular molecule also across
    /// the origin.
    /// </remarks>
    /// <exception cref="ArgumentException">A methylase site is malformed, or the enzyme is unknown.</exception>
    public static IEnumerable<RestrictionSite> FindSites(
        DnaSequence sequence, string enzymeName, MoleculeTopology topology,
        IReadOnlyList<Methylase>? methylases)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var enzyme = RequireEnzyme(enzymeName);
        return FindSitesCore(sequence.Sequence, enzyme, topology,
            MethylationMask.Compute(sequence.Sequence, methylases, topology == MoleculeTopology.Circular));
    }

    /// <summary>
    /// <see cref="FindSites(DnaSequence, RestrictionEnzyme, MoleculeTopology)"/> with methylated
    /// bases blocking recognition; see
    /// <see cref="FindSites(DnaSequence, string, MoleculeTopology, IReadOnlyList{Methylase})"/>.
    /// </summary>
    public static IEnumerable<RestrictionSite> FindSites(
        DnaSequence sequence, RestrictionEnzyme enzyme, MoleculeTopology topology,
        IReadOnlyList<Methylase>? methylases)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(enzyme);
        return FindSitesCore(sequence.Sequence, enzyme, topology,
            MethylationMask.Compute(sequence.Sequence, methylases, topology == MoleculeTopology.Circular));
    }

    /// <summary>
    /// <see cref="FindSites(DnaSequence, MoleculeTopology, string[])"/> with methylated bases blocking
    /// recognition; the per-enzyme results are concatenated in the order given.
    /// </summary>
    public static IEnumerable<RestrictionSite> FindSites(
        DnaSequence sequence, MoleculeTopology topology, IReadOnlyList<Methylase>? methylases,
        params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (enzymeNames == null || enzymeNames.Length == 0)
            throw new ArgumentException("At least one enzyme name is required", nameof(enzymeNames));
        return enzymeNames.SelectMany(name => FindSites(sequence, name, topology, methylases));
    }

    /// <summary>
    /// <see cref="Digest(DnaSequence, MoleculeTopology, string[])"/> with methylated bases blocking
    /// recognition; see <see cref="FindSites(DnaSequence, string, MoleculeTopology, IReadOnlyList{Methylase})"/>.
    /// </summary>
    public static IEnumerable<DigestFragment> Digest(
        DnaSequence sequence, MoleculeTopology topology, IReadOnlyList<Methylase>? methylases,
        params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (enzymeNames == null || enzymeNames.Length == 0)
            throw new ArgumentException("At least one enzyme is required", nameof(enzymeNames));
        return DigestWithMask(sequence, topology, ResolveEnzymes(enzymeNames), methylases);
    }

    /// <summary>
    /// <see cref="Digest(DnaSequence, MoleculeTopology, IEnumerable{RestrictionEnzyme})"/> with
    /// methylated bases blocking recognition.
    /// </summary>
    public static IEnumerable<DigestFragment> Digest(
        DnaSequence sequence, MoleculeTopology topology, IEnumerable<RestrictionEnzyme> enzymes,
        IReadOnlyList<Methylase>? methylases)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return DigestWithMask(sequence, topology, RequireEnzymeList(enzymes), methylases);
    }

    private static IEnumerable<DigestFragment> DigestWithMask(
        DnaSequence sequence, MoleculeTopology topology, IEnumerable<RestrictionEnzyme> enzymes,
        IReadOnlyList<Methylase>? methylases)
    {
        bool circular = topology == MoleculeTopology.Circular;
        var blocked = MethylationMask.Compute(sequence.Sequence, methylases, circular);
        return circular ? DigestCircular(sequence, enzymes, blocked) : DigestCore(sequence, enzymes, blocked);
    }

    /// <summary>
    /// <see cref="GetDigestSummary(DnaSequence, MoleculeTopology, string[])"/> with methylated
    /// bases blocking recognition.
    /// </summary>
    public static DigestSummary GetDigestSummary(
        DnaSequence sequence, MoleculeTopology topology, IReadOnlyList<Methylase>? methylases,
        params string[] enzymeNames)
    {
        var fragments = Digest(sequence, topology, methylases, enzymeNames).ToList();
        return Summarize(fragments, enzymeNames.ToList(), topology);
    }

    /// <summary>
    /// <see cref="GetDigestSummary(DnaSequence, MoleculeTopology, IEnumerable{RestrictionEnzyme})"/>
    /// with methylated bases blocking recognition.
    /// </summary>
    public static DigestSummary GetDigestSummary(
        DnaSequence sequence, MoleculeTopology topology, IEnumerable<RestrictionEnzyme> enzymes,
        IReadOnlyList<Methylase>? methylases)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var list = RequireEnzymeList(enzymes);
        var fragments = Digest(sequence, topology, list, methylases).ToList();
        return Summarize(fragments, list.Select(e => e.Name).ToList(), topology);
    }

    /// <summary>
    /// <see cref="CreateMap(DnaSequence, MoleculeTopology, string[])"/> with methylated bases
    /// blocking recognition.
    /// </summary>
    public static RestrictionMap CreateMap(
        DnaSequence sequence, MoleculeTopology topology, IReadOnlyList<Methylase>? methylases,
        params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(enzymeNames);
        var considered = (enzymeNames.Length > 0
                ? enzymeNames.Select(n => RequireEnzyme(n))
                : _enzymes.Values)
            .DistinctBy(e => e.Name, StringComparer.OrdinalIgnoreCase)
            .ToList();
        return CreateMapCore(sequence, topology, considered,
            MethylationMask.Compute(sequence.Sequence, methylases, topology == MoleculeTopology.Circular));
    }

    /// <summary>
    /// <see cref="CreateMap(DnaSequence, MoleculeTopology, IEnumerable{RestrictionEnzyme})"/> with
    /// methylated bases blocking recognition.
    /// </summary>
    public static RestrictionMap CreateMap(
        DnaSequence sequence, MoleculeTopology topology, IEnumerable<RestrictionEnzyme> enzymes,
        IReadOnlyList<Methylase>? methylases)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var considered = RequireEnzymeList(enzymes)
            .DistinctBy(e => e.Name, StringComparer.OrdinalIgnoreCase)
            .ToList();
        return CreateMapCore(sequence, topology, considered,
            MethylationMask.Compute(sequence.Sequence, methylases, topology == MoleculeTopology.Circular));
    }
}
