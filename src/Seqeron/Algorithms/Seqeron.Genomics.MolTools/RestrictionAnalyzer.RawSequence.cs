namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Raw-sequence-string entry points: the same searches, digests and maps as the
/// <see cref="DnaSequence"/> overloads, over a sequence this library will not reject for carrying
/// IUPAC ambiguity codes.
/// </summary>
/// <remarks>
/// <para>
/// <see cref="DnaSequence"/> accepts only A/C/G/T, so an <c>N</c>-containing contig or a
/// soft-masked assembly cannot be handed to the <see cref="DnaSequence"/> overloads at all — and
/// before review 2026-09 B08 H27 every restriction MCP tool built a <see cref="DnaSequence"/>, so
/// the IUPAC-tolerant scan that
/// <see cref="RestrictionAnalyzer.FindSites(string, string, MoleculeTopology)"/> and
/// <see cref="RestrictionAnalyzer.FindPatternOccurrences"/> document was unreachable over MCP.
/// These overloads are the reachable form of it for the whole family.
/// </para>
/// <para>
/// The matching rule is the raw-string family's, which is Biopython 1.88's: the sequence is
/// upper-cased with an ASCII-only fold, a site <c>N</c> matches any of the 15 IUPAC nucleotide
/// codes in the sequence (an <c>N</c> included), every other site code matches only a concrete
/// A/C/G/T base it expands to — so an ambiguity code in the sequence blocks a concrete site code
/// rather than widening it — and a character outside the 15 codes matches nothing. A
/// <see langword="null"/> or empty sequence yields the natural empty result, as it does on
/// <see cref="RestrictionAnalyzer.FindSites(string, string)"/>; an unknown enzyme name is still
/// rejected eagerly, before the sequence is looked at.
/// </para>
/// <para>
/// Coordinates, the two-records-per-palindromic-site rule, the linear cut-feasibility filter, the
/// circular wrap and the methylation mask are all exactly as documented on the corresponding
/// <see cref="DnaSequence"/> overload; these methods only widen the accepted alphabet. Each one
/// forwards to the same private core, so there is no second implementation.
/// </para>
/// </remarks>
public static partial class RestrictionAnalyzer
{
    #region Raw-sequence-string entry points (B08 H27)

    /// <summary>
    /// <see cref="FindAllSites(DnaSequence, MoleculeTopology, EnzymeSet, IReadOnlyList{Methylase})"/>
    /// over a raw sequence string that may carry IUPAC ambiguity codes.
    /// </summary>
    /// <param name="sequence">Sequence to scan; null or empty yields no sites.</param>
    /// <param name="topology">Whether the molecule is linear or circular.</param>
    /// <param name="set">Which enzyme table to scan with.</param>
    /// <param name="methylases">Methylases whose sites block recognition, or null for none.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="set"/> is not a defined value.</exception>
    /// <exception cref="ArgumentException">A methylase site is malformed.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<RestrictionSite> FindAllSites(
        string? sequence, MoleculeTopology topology, EnzymeSet set, IReadOnlyList<Methylase>? methylases)
    {
        RequireTopology(topology);
        var table = GetEnzymes(set);
        string seq = Normalize(sequence);
        // The mask is computed for an empty sequence too, so a malformed methylase is rejected at the call whatever the
        // sequence (the documented ArgumentException; B08 round 5 A5R-6).
        var blocked = MethylationMask.Compute(seq, methylases, topology == MoleculeTopology.Circular);
        if (seq.Length == 0)
            return [];

        return table.Values.SelectMany(e => FindSitesCore(seq, e, topology, blocked));
    }

    /// <summary>
    /// <see cref="FindSites(DnaSequence, MoleculeTopology, IEnumerable{RestrictionEnzyme}, IReadOnlyList{Methylase})"/>
    /// over a raw sequence string that may carry IUPAC ambiguity codes. Accepts caller-supplied
    /// enzyme definitions, including ones that are not in the built-in table.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="enzymes"/> is null or contains null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymes"/> is empty, an enzyme definition
    /// is invalid, or a methylase site is malformed.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<RestrictionSite> FindSites(
        string? sequence, MoleculeTopology topology, IEnumerable<RestrictionEnzyme> enzymes,
        IReadOnlyList<Methylase>? methylases)
    {
        RequireTopology(topology);
        var list = RequireEnzymeList(enzymes);
        string seq = Normalize(sequence);
        var blocked = MethylationMask.Compute(seq, methylases, topology == MoleculeTopology.Circular);
        if (seq.Length == 0)
            return [];

        return list.SelectMany(e => FindSitesCore(seq, e, topology, blocked));
    }

    /// <summary>
    /// <see cref="FindSites(DnaSequence, MoleculeTopology, IReadOnlyList{Methylase}, string[])"/>
    /// over a raw sequence string that may carry IUPAC ambiguity codes.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeNames"/> is null, or a name in
    /// it is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeNames"/> is empty, a name in it is
    /// empty or not a REBASE enzyme with a defined cut, or a methylase site is malformed.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<RestrictionSite> FindSites(
        string? sequence, MoleculeTopology topology, IReadOnlyList<Methylase>? methylases,
        params string[] enzymeNames)
        => FindSites(sequence, topology, RequireEnzymeNames(enzymeNames, nameof(enzymeNames)), methylases);

    /// <summary>
    /// <see cref="Digest(DnaSequence, MoleculeTopology, IReadOnlyList{Methylase}, string[])"/> over a
    /// raw sequence string that may carry IUPAC ambiguity codes.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeNames"/> is null, or a name in
    /// it is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeNames"/> is empty, a name in it is
    /// empty or not a REBASE enzyme with a defined cut, or a methylase site is malformed.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<DigestFragment> Digest(
        string? sequence, MoleculeTopology topology, IReadOnlyList<Methylase>? methylases,
        params string[] enzymeNames)
        => Digest(sequence, topology, RequireEnzymeNames(enzymeNames, nameof(enzymeNames)), methylases);

    /// <summary>
    /// <see cref="Digest(DnaSequence, MoleculeTopology, IEnumerable{RestrictionEnzyme}, IReadOnlyList{Methylase})"/>
    /// over a raw sequence string that may carry IUPAC ambiguity codes.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="enzymes"/> is null or contains null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymes"/> is empty, an enzyme definition
    /// is invalid, or a methylase site is malformed.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<DigestFragment> Digest(
        string? sequence, MoleculeTopology topology, IEnumerable<RestrictionEnzyme> enzymes,
        IReadOnlyList<Methylase>? methylases)
    {
        RequireTopology(topology);
        var list = RequireEnzymeList(enzymes);
        string seq = Normalize(sequence);
        var blocked = MethylationMask.Compute(seq, methylases, topology == MoleculeTopology.Circular);

        return topology == MoleculeTopology.Linear
            ? DigestCore(seq, list, blocked)
            : DigestCircular(seq, list, blocked);
    }

    /// <summary>
    /// <see cref="GetDigestSummary(DnaSequence, MoleculeTopology, IReadOnlyList{Methylase}, string[])"/>
    /// over a raw sequence string that may carry IUPAC ambiguity codes.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeNames"/> is null, or a name in
    /// it is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeNames"/> is empty, a name in it is
    /// empty or not a REBASE enzyme with a defined cut, or a methylase site is malformed.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static DigestSummary GetDigestSummary(
        string? sequence, MoleculeTopology topology, IReadOnlyList<Methylase>? methylases,
        params string[] enzymeNames)
    {
        RequireTopology(topology);
        var list = RequireEnzymeNames(enzymeNames, nameof(enzymeNames));
        var fragments = Digest(sequence, topology, list, methylases).ToList();
        return Summarize(fragments, enzymeNames.ToList(), topology);
    }

    /// <summary>
    /// <see cref="GetDigestSummary(DnaSequence, MoleculeTopology, IEnumerable{RestrictionEnzyme}, IReadOnlyList{Methylase})"/>
    /// over a raw sequence string that may carry IUPAC ambiguity codes.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="enzymes"/> is null or contains null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymes"/> is empty, an enzyme definition
    /// is invalid, or a methylase site is malformed.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static DigestSummary GetDigestSummary(
        string? sequence, MoleculeTopology topology, IEnumerable<RestrictionEnzyme> enzymes,
        IReadOnlyList<Methylase>? methylases)
    {
        RequireTopology(topology);
        var list = RequireEnzymeList(enzymes);
        var fragments = Digest(sequence, topology, list, methylases).ToList();
        return Summarize(fragments, list.Select(e => e.Name).ToList(), topology);
    }

    /// <summary>
    /// <see cref="PartialDigest(DnaSequence, MoleculeTopology, int?, IReadOnlyList{Methylase}, string[])"/>
    /// over a raw sequence string that may carry IUPAC ambiguity codes.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeNames"/> is null, or a name in
    /// it is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMissedCuts"/> is negative.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeNames"/> is empty, a name in it is
    /// empty or not a REBASE enzyme with a defined cut, the fragment count exceeds the limit, or a methylase
    /// site is malformed.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<PartialDigestFragment> PartialDigest(
        string? sequence, MoleculeTopology topology, int? maxMissedCuts,
        IReadOnlyList<Methylase>? methylases, params string[] enzymeNames)
        => PartialDigest(
            sequence, topology, maxMissedCuts,
            RequireEnzymeNames(enzymeNames, nameof(enzymeNames)), methylases);

    /// <summary>
    /// <see cref="PartialDigest(DnaSequence, MoleculeTopology, int?, IEnumerable{RestrictionEnzyme}, IReadOnlyList{Methylase})"/>
    /// over a raw sequence string that may carry IUPAC ambiguity codes.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="enzymes"/> is null or contains null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMissedCuts"/> is negative.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymes"/> is empty, an enzyme definition
    /// is invalid, the fragment count exceeds the limit, or a methylase site is malformed.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<PartialDigestFragment> PartialDigest(
        string? sequence, MoleculeTopology topology, int? maxMissedCuts,
        IEnumerable<RestrictionEnzyme> enzymes, IReadOnlyList<Methylase>? methylases)
    {
        RequireTopology(topology);
        var list = RequireEnzymeList(enzymes);
        string seq = Normalize(sequence);
        return PartialDigestCore(
            seq, topology, maxMissedCuts, list,
            MethylationMask.Compute(seq, methylases, topology != MoleculeTopology.Linear));
    }

    /// <summary>
    /// <see cref="CreateMap(DnaSequence, MoleculeTopology, IReadOnlyList{Methylase}, string[])"/> over
    /// a raw sequence string that may carry IUPAC ambiguity codes. An empty
    /// <paramref name="enzymeNames"/> considers every common enzyme, exactly as that overload does.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeNames"/> is null, or a name in
    /// it is null.</exception>
    /// <exception cref="ArgumentException">A name in <paramref name="enzymeNames"/> is empty or not a REBASE
    /// enzyme with a defined cut, or a methylase site is malformed.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static RestrictionMap CreateMap(
        string? sequence, MoleculeTopology topology, IReadOnlyList<Methylase>? methylases,
        params string[] enzymeNames)
    {
        RequireTopology(topology);
        ArgumentNullException.ThrowIfNull(enzymeNames);
        var considered = (enzymeNames.Length > 0
                ? enzymeNames.Select(n => RequireEnzyme(n, nameof(enzymeNames)))
                : Enzymes.Values)
            .DistinctBy(e => e.Name, StringComparer.OrdinalIgnoreCase)
            .ToList();

        return CreateMapRaw(sequence, topology, considered, methylases);
    }

    /// <summary>
    /// <see cref="CreateMap(DnaSequence, MoleculeTopology, IEnumerable{RestrictionEnzyme}, IReadOnlyList{Methylase})"/>
    /// over a raw sequence string that may carry IUPAC ambiguity codes.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="enzymes"/> is null or contains null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymes"/> is empty, an enzyme definition
    /// is invalid, or a methylase site is malformed.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static RestrictionMap CreateMap(
        string? sequence, MoleculeTopology topology, IEnumerable<RestrictionEnzyme> enzymes,
        IReadOnlyList<Methylase>? methylases)
    {
        RequireTopology(topology);
        var considered = RequireEnzymeList(enzymes)
            .DistinctBy(e => e.Name, StringComparer.OrdinalIgnoreCase)
            .ToList();

        return CreateMapRaw(sequence, topology, considered, methylases);
    }

    private static RestrictionMap CreateMapRaw(
        string? sequence, MoleculeTopology topology, List<RestrictionEnzyme> considered,
        IReadOnlyList<Methylase>? methylases)
    {
        string seq = Normalize(sequence);
        return CreateMapCore(
            seq, topology, considered,
            MethylationMask.Compute(seq, methylases, topology == MoleculeTopology.Circular));
    }

    /// <summary>
    /// A raw sequence string in the form the private cores expect: ASCII-upper-cased, with null
    /// treated as empty. The ASCII-only fold matters because <c>ToUpperInvariant</c> would map
    /// U+017F 'ſ' onto the IUPAC code 'S' and make a non-IUPAC character match (B08 F49).
    /// </summary>
    private static string Normalize(string? sequence)
        => string.IsNullOrEmpty(sequence) ? string.Empty : SequenceExtensions.ToUpperAscii(sequence);

    #endregion
}
