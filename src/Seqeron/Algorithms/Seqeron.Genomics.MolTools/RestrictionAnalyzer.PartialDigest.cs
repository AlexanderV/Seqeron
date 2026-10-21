namespace Seqeron.Genomics.MolTools;

// ───────────────────────────────────────────────────────────────────────────────────────────
// Partial digestion (review 2026-09 B08 D8 / F53). Kept in its own partial-class file so that
// the line anchors other documents hold into RestrictionAnalyzer.cs do not move.
//
// Definition used (combinatorial, no kinetics): a partial digest is a digest in which any
// subset of the cut sites may be left uncut, so every fragment is an interval between two cut
// positions (or a linear molecule's terminus) with zero or more uncut sites inside it. For a
// linear molecule with cut positions c_1 < … < c_k this is the set of intervals [a, b), a < b,
// over X = {0, c_1, …, c_k, L} — the "difference multiset ΔX" of the partial digest problem,
// which has C(k+2, 2) = (k+1)(k+2)/2 elements (Rosalind PDPL "Creating a Restriction Map",
// after Jones & Pevzner, An Introduction to Bioinformatics Algorithms). The cut positions are
// exactly those the complete Digest uses (CollectForwardCuts), and every fragment end carries
// the same enzyme list and overhang as the complete digest's end at that cut.
//
// Out of scope: the yield of each fragment. The relative abundance of partial products depends
// on enzyme amount, time and per-site cleavage rates; no source opened gives a model for the
// built-in enzymes, so only WHICH fragments can arise is computed, never HOW MUCH of each.
// ───────────────────────────────────────────────────────────────────────────────────────────
public static partial class RestrictionAnalyzer
{
    /// <summary>
    /// The largest number of fragments <see cref="PartialDigest(DnaSequence, MoleculeTopology, int?, string[])"/>
    /// will enumerate. A linear molecule with <c>k</c> cuts has <c>(k+1)(k+2)/2</c> partial-digest
    /// fragments and a circular one <c>k² + 1</c>, so without a missed-cut limit the count grows
    /// quadratically in the number of sites; a request above this bound throws and asks the caller
    /// to pass <c>maxMissedCuts</c>.
    /// </summary>
    public const int MaxPartialDigestFragments = 100_000;

    /// <summary>
    /// Simulates a <b>partial</b> digest of a <b>linear</b> molecule: every fragment that can arise
    /// when any subset of the cut sites is left uncut. Equivalent to
    /// <see cref="PartialDigest(DnaSequence, MoleculeTopology, int?, string[])"/> with
    /// <see cref="MoleculeTopology.Linear"/> and no missed-cut limit.
    /// </summary>
    /// <param name="sequence">DNA sequence to digest.</param>
    /// <param name="enzymeNames">Names of restriction enzymes to use.</param>
    public static IEnumerable<PartialDigestFragment> PartialDigest(DnaSequence sequence, params string[] enzymeNames)
        => PartialDigest(sequence, MoleculeTopology.Linear, null, enzymeNames);

    /// <summary>
    /// Simulates a <b>partial</b> digest of a molecule with the given topology, with no missed-cut
    /// limit. See <see cref="PartialDigest(DnaSequence, MoleculeTopology, int?, string[])"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence to digest.</param>
    /// <param name="topology">Whether the molecule is linear or circular.</param>
    /// <param name="enzymeNames">Names of restriction enzymes to use.</param>
    public static IEnumerable<PartialDigestFragment> PartialDigest(
        DnaSequence sequence, MoleculeTopology topology, params string[] enzymeNames)
        => PartialDigest(sequence, topology, null, enzymeNames);

    /// <summary>
    /// Simulates a <b>partial</b> digest: returns every fragment that can arise when any subset of
    /// the cut sites is left uncut, each tagged with the number of uncut sites it contains.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Cut positions.</b> The same distinct top-strand cut positions, with the same enzymes per
    /// cut, as <see cref="Digest(DnaSequence, MoleculeTopology, string[])"/> — partial digestion
    /// needs no data beyond the complete digest's cut list. Every fragment end that is a cut
    /// carries exactly the information the complete digest's end at that cut carries
    /// (<see cref="DigestFragment.LeftEnzyme"/>/<see cref="DigestFragment.LeftEnzymes"/>, the
    /// overhang type and sequence, and the same on the right); an end that is a linear molecule's
    /// original terminus reports <see langword="null"/> enzymes and overhang.
    /// </para>
    /// <para>
    /// <b>Linear.</b> With cut positions <c>c_1 &lt; … &lt; c_k</c> the fragments are all intervals
    /// <c>[a, b)</c>, <c>a &lt; b</c>, over <c>{0, c_1, …, c_k, L}</c> — <c>(k+1)(k+2)/2</c> of
    /// them, the difference multiset of the partial digest problem. A fragment's
    /// <see cref="PartialDigestFragment.MissedCuts"/> is the number of cut positions strictly
    /// inside it; the whole molecule is the fragment with all <c>k</c> sites missed.
    /// </para>
    /// <para>
    /// <b>Circular.</b> With <c>k ≥ 1</c> cut positions every fragment runs from one cut to a later
    /// cut, cyclically, wrapping past the origin when it must: for each of the <c>k</c> start cuts
    /// one fragment per <c>0 … k−1</c> missed cuts, where <c>k−1</c> missed cuts is the full-length
    /// molecule linearised at that start cut (the complete digest's one-cut result). Leaving all
    /// <c>k</c> sites uncut gives the intact circle, reported exactly as the complete digest
    /// reports a circle without a cut — <c>StartPosition 0</c>, full length, no ends — with
    /// <see cref="PartialDigestFragment.MissedCuts"/> <c>= k</c>. Total <c>k² + 1</c>. With
    /// <c>k = 0</c> the only fragment is that intact circle, with <c>MissedCuts = 0</c>
    /// (no site exists to be missed) — the same single fragment the complete digest returns.
    /// </para>
    /// <para>
    /// <b>Order and numbering.</b> Fragments are ordered by
    /// <see cref="PartialDigestFragment.MissedCuts"/> ascending, then
    /// <see cref="DigestFragment.StartPosition"/> ascending, and
    /// <see cref="DigestFragment.FragmentNumber"/> is the 1-based position in that order. The
    /// fragments with <c>MissedCuts == 0</c> therefore come first and are, member for member
    /// (numbering included), the complete digest's fragments; with <paramref name="maxMissedCuts"/>
    /// = <c>m</c> the result is exactly the prefix of the unlimited result with
    /// <c>MissedCuts ≤ m</c>.
    /// </para>
    /// <para>
    /// <b>Out of scope.</b> Only which fragments CAN arise is computed; the yield (relative
    /// abundance) of each partial product is not modelled.
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to digest.</param>
    /// <param name="topology">Whether the molecule is linear or circular.</param>
    /// <param name="maxMissedCuts">Optional upper bound on the number of uncut sites inside a
    /// fragment; <c>0</c> gives the complete digest; <see langword="null"/> means no bound.</param>
    /// <param name="enzymeNames">Names of restriction enzymes to use.</param>
    /// <returns>The partial-digest fragments, enumerated lazily.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentException">No enzyme was given, an enzyme name is unknown, or the
    /// fragment count would exceed <see cref="MaxPartialDigestFragments"/> (pass a smaller
    /// <paramref name="maxMissedCuts"/>).</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMissedCuts"/> is negative.</exception>
    public static IEnumerable<PartialDigestFragment> PartialDigest(
        DnaSequence sequence,
        MoleculeTopology topology,
        int? maxMissedCuts,
        params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (enzymeNames == null || enzymeNames.Length == 0)
            throw new ArgumentException("At least one enzyme is required", nameof(enzymeNames));

        return PartialDigestCore(sequence, topology, maxMissedCuts, ResolveEnzymes(enzymeNames));
    }

    /// <summary>
    /// Simulates a <b>partial</b> digest with caller-supplied enzyme definitions — including
    /// enzymes that are not in the built-in table, such as Type IIS cutters whose cut lies outside
    /// the recognition site. Same contract as
    /// <see cref="PartialDigest(DnaSequence, MoleculeTopology, int?, string[])"/>; the cut list is
    /// the one <see cref="Digest(DnaSequence, MoleculeTopology, IEnumerable{RestrictionEnzyme})"/>
    /// uses.
    /// </summary>
    /// <param name="sequence">DNA sequence to digest.</param>
    /// <param name="topology">Whether the molecule is linear or circular.</param>
    /// <param name="maxMissedCuts">Optional upper bound on the number of uncut sites inside a
    /// fragment; <see langword="null"/> means no bound.</param>
    /// <param name="enzymes">Enzyme definitions to digest with; at least one.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or
    /// <paramref name="enzymes"/> is null, or <paramref name="enzymes"/> contains null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymes"/> is empty, an enzyme's
    /// recognition sequence is invalid, or the fragment count would exceed
    /// <see cref="MaxPartialDigestFragments"/>.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMissedCuts"/> is negative.</exception>
    public static IEnumerable<PartialDigestFragment> PartialDigest(
        DnaSequence sequence,
        MoleculeTopology topology,
        int? maxMissedCuts,
        IEnumerable<RestrictionEnzyme> enzymes)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return PartialDigestCore(sequence, topology, maxMissedCuts, RequireEnzymeList(enzymes));
    }

    private static IEnumerable<PartialDigestFragment> PartialDigestCore(
        DnaSequence sequence, MoleculeTopology topology, int? maxMissedCuts, IEnumerable<RestrictionEnzyme> enzymes, bool[]? blocked = null)
    {
        if (maxMissedCuts < 0)
            throw new ArgumentOutOfRangeException(nameof(maxMissedCuts), maxMissedCuts,
                "The missed-cut limit cannot be negative.");

        bool circular = topology != MoleculeTopology.Linear;
        var (cuts, sitesByPosition) = CollectForwardCuts(
            sequence, enzymes, circular ? MoleculeTopology.Circular : MoleculeTopology.Linear, blocked);

        long count = CountPartialDigestFragments(cuts.Count, circular, maxMissedCuts);
        if (count > MaxPartialDigestFragments)
            throw new ArgumentException(
                $"A partial digest with {cuts.Count} cut positions would produce {count} fragments, above the limit of " +
                $"{MaxPartialDigestFragments}; pass a smaller maxMissedCuts.", nameof(maxMissedCuts));

        return circular
            ? PartialDigestCircular(sequence.Sequence, cuts, sitesByPosition, maxMissedCuts)
            : PartialDigestLinear(sequence.Sequence, cuts, sitesByPosition, maxMissedCuts);
    }

    /// <summary>
    /// The number of partial-digest fragments for <paramref name="k"/> cut positions: linear
    /// <c>Σ_{j=0..min(m,k)} (k+1−j)</c>, i.e. <c>(k+1)(k+2)/2</c> unbounded; circular
    /// <c>k·min(m+1, k)</c> plus the intact circle when <c>m ≥ k</c>, i.e. <c>k² + 1</c> unbounded
    /// (1 when <c>k = 0</c>).
    /// </summary>
    private static long CountPartialDigestFragments(int k, bool circular, int? maxMissedCuts)
    {
        long m = maxMissedCuts ?? long.MaxValue;
        if (!circular)
        {
            long top = Math.Min(m, k);
            long total = 0;
            for (long j = 0; j <= top; j++)
                total += k + 1 - j;
            return total;
        }

        if (k == 0)
            return 1;
        return (long)k * (Math.Min(m, k - 1) + 1) + (m >= k ? 1 : 0);
    }

    private static IEnumerable<PartialDigestFragment> PartialDigestLinear(
        string seq, List<int> cuts, Dictionary<int, List<RestrictionSite>> sitesByPosition, int? maxMissedCuts)
    {
        // Boundaries X = {0, c_1, …, c_k, L}; boundary 0 and L are the molecule's own termini.
        var boundaries = new List<int>(cuts.Count + 2) { 0 };
        boundaries.AddRange(cuts);
        boundaries.Add(seq.Length);
        int last = boundaries.Count - 1;
        int k = cuts.Count;
        int maxMissed = Math.Min(maxMissedCuts ?? k, k);

        int fragmentNumber = 0;
        for (int missed = 0; missed <= maxMissed; missed++)
        {
            for (int a = 0; a + missed + 1 <= last; a++)
            {
                int b = a + missed + 1;
                int start = boundaries[a];
                int length = boundaries[b] - start;
                var leftSites = a == 0 ? null : sitesByPosition[start];
                var rightSites = b == last ? null : sitesByPosition[boundaries[b]];
                yield return new PartialDigestFragment(
                    PartialFragment(seq, seq.Substring(start, length), start, length,
                        leftSites, rightSites, MoleculeTopology.Linear, ++fragmentNumber),
                    missed);
            }
        }
    }

    private static IEnumerable<PartialDigestFragment> PartialDigestCircular(
        string seq, List<int> cuts, Dictionary<int, List<RestrictionSite>> sitesByPosition, int? maxMissedCuts)
    {
        int n = seq.Length;
        int k = cuts.Count;
        int fragmentNumber = 0;

        if (k == 0)
        {
            // No site at all: the intact circle, as the complete digest reports it; nothing missed.
            yield return new PartialDigestFragment(UncutCircle(seq, 1), 0);
            yield break;
        }

        int maxMissed = Math.Min(maxMissedCuts ?? k, k - 1);
        for (int missed = 0; missed <= maxMissed; missed++)
        {
            for (int i = 0; i < k; i++)
            {
                int start = cuts[i];
                int end = cuts[(i + missed + 1) % k];
                // Distance travelled forward from start to end; 0 means the end is the start cut
                // itself (all other k−1 sites missed) — the full-length linearised molecule.
                int length = (end - start + n) % n;
                if (length == 0)
                    length = n;

                string fragmentSeq = start + length <= n
                    ? seq.Substring(start, length)
                    : string.Concat(seq.AsSpan(start, n - start), seq.AsSpan(0, length - (n - start)));

                yield return new PartialDigestFragment(
                    PartialFragment(seq, fragmentSeq, start, length,
                        sitesByPosition[start], sitesByPosition[end], MoleculeTopology.Circular, ++fragmentNumber),
                    missed);
            }
        }

        if (maxMissedCuts is null || maxMissedCuts >= k)
            yield return new PartialDigestFragment(UncutCircle(seq, fragmentNumber + 1), k);
    }

    /// <summary>
    /// Builds one fragment with the same end information the complete digest attaches to a cut
    /// end: the first site's enzyme as the representative name, the whole ordinal list, and the
    /// representative's overhang (<see cref="DownstreamEnd"/>/<see cref="UpstreamEnd"/>).
    /// </summary>
    private static DigestFragment PartialFragment(
        string seq, string fragmentSeq, int start, int length,
        List<RestrictionSite>? leftSites, List<RestrictionSite>? rightSites,
        MoleculeTopology topology, int fragmentNumber)
    {
        var leftSite = leftSites?[0];
        var rightSite = rightSites?[0];
        var leftEnd = DownstreamEnd(seq, leftSite, topology);
        var rightEnd = UpstreamEnd(seq, rightSite, topology);

        return new DigestFragment(
            Sequence: fragmentSeq,
            StartPosition: start,
            Length: length,
            LeftEnzyme: leftSite?.Enzyme.Name,
            RightEnzyme: rightSite?.Enzyme.Name,
            FragmentNumber: fragmentNumber,
            LeftOverhangType: leftEnd.Type,
            LeftOverhang: leftEnd.Sequence,
            RightOverhangType: rightEnd.Type,
            RightOverhang: rightEnd.Sequence)
        {
            LeftEnzymes = EnzymeNames(leftSites),
            RightEnzymes = EnzymeNames(rightSites)
        };
    }

    /// <summary>The intact circle, exactly as the complete circular digest reports a non-cutter.</summary>
    private static DigestFragment UncutCircle(string seq, int fragmentNumber) =>
        new(Sequence: seq, StartPosition: 0, Length: seq.Length,
            LeftEnzyme: null, RightEnzyme: null, FragmentNumber: fragmentNumber);
}

/// <summary>
/// One fragment of a partial digest
/// (<see cref="RestrictionAnalyzer.PartialDigest(DnaSequence, MoleculeTopology, int?, string[])"/>).
/// </summary>
/// <param name="Fragment">The fragment, with the same members — sequence, start, length, end
/// enzymes and overhangs — a complete-digest fragment carries.</param>
/// <param name="MissedCuts">How many cut positions lie strictly inside the fragment, i.e. were
/// left uncut to produce it; <c>0</c> for a complete-digest fragment. For the intact circle of a
/// circular molecule with <c>k</c> sites it is <c>k</c>.</param>
public sealed record PartialDigestFragment(DigestFragment Fragment, int MissedCuts);
