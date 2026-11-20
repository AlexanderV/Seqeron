using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Analysis;

/// <summary>A cluster of <see cref="GenomicAnalyzer.FindNucmerClusters"/> with the strand of the query it lies on.</summary>
/// <param name="IsReverseComplement">The cluster is on the reverse-complement strand of the query: the second start of every match is a 1-based position in the reverse-complemented query (what <c>mummer</c> / <c>mgaps</c> print).</param>
/// <param name="Cluster">The cluster (starts 1-based, as <c>mgaps</c> prints them).</param>
public sealed record StrandedNucmerCluster(bool IsReverseComplement, NucmerCluster Cluster);

public static partial class GenomicAnalyzer
{
    #region nucmer clusters (mummer | mgaps)

    /// <summary>
    /// The anchor-clustering stage of MUMmer 3 <c>nucmer</c>: the exact matches of <c>mummer</c> (<paramref name="mode"/>: <c>-mum</c> = <see cref="CommonRegionMode.MaximalUnique"/> by default,
    /// <c>-mumreference</c>, <c>-maxmatch</c>; <paramref name="minMatchLength"/> = <c>-l</c>, default 20; only A/C/G/T match, as <c>mummer -n</c>) are clustered by <c>mgaps</c>
    /// (<see cref="NucmerClusters"/>: <c>-c</c> minimum cluster length 65, <c>-g</c> maximum separation 90, <c>-D</c> diagonal difference 5, <c>-d</c> diagonal factor 0.12).
    /// The forward strand is clustered first, then the reverse-complement one.
    /// </summary>
    /// <remarks>
    /// This is the part of nucmer that decides which anchors belong together; the stages after it (<c>postnuc</c>: banded extension of the clusters and of the gaps between their matches,
    /// then <c>delta-filter</c>) are not part of this method. Coordinates are 1-based as in the programs' output. Verified against the installed
    /// <c>mummer | mgaps</c> pipeline of <c>nucmer</c> 3.23 on 700 related pairs (docs/Evidence/ALIGN-MULTI-001-Evidence.md).
    /// </remarks>
    /// <param name="reference">Reference sequence (the indexed one).</param>
    /// <param name="query">Query sequence.</param>
    /// <param name="minMatchLength">Minimum length of an exact match (<c>nucmer -l</c>).</param>
    /// <param name="mode">Which exact matches seed the clusters.</param>
    /// <param name="options">Clustering options (null = nucmer defaults).</param>
    /// <param name="strands">Strands of the query to cluster (default both, <c>nucmer</c>'s default).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <returns>The clusters of the forward strand, then of the reverse-complement strand, in the order <c>mgaps</c> prints them.</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static IReadOnlyList<StrandedNucmerCluster> FindNucmerClusters(
        string reference, string query, int minMatchLength = 20, CommonRegionMode mode = CommonRegionMode.MaximalUnique,
        NucmerClusterOptions? options = null, StrandSelection strands = StrandSelection.Both, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(reference);
        ArgumentNullException.ThrowIfNull(query);
        var regions = FindCommonRegions(reference, query, minMatchLength, mode, strands, cancellationToken);
        var result = new List<StrandedNucmerCluster>();
        foreach (bool reverse in new[] { false, true })
        {
            var matches = regions.Where(r => r.IsReverseComplement == reverse)
                .Select(r => ((long)r.PositionInFirst + 1, reverse ? (long)query.Length - r.PositionInSecond - r.Length + 1 : r.PositionInSecond + 1L, (long)r.Length))
                .ToList();
            foreach (var cluster in NucmerClusters.Cluster(matches, options))
                result.Add(new StrandedNucmerCluster(reverse, cluster));
        }

        return result;
    }

    #endregion
}
