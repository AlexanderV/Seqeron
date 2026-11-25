using System.ComponentModel;
using ModelContextProtocol.Server;
using global::Seqeron.Genomics.Analysis;

namespace Seqeron.Mcp.Analysis.Tools;

// EMBOSS program ports on the Analysis server (B09 A3-8 / A3-12 / A3-13); every tool delegates to GenomicAnalyzer.
public partial class AnalysisTools
{

    [McpServerTool(Name = "nucmer_clusters", Title = "Genome — MUMmer nucmer Anchor Clusters", ReadOnly = true)]
    [Description("The anchor-clustering stage of MUMmer 3.23 nucmer = `mummer -mum|-mumreference|-maxmatch -b -l L -n | mgaps -l c -s g -d D -f F` (700 related pairs, both strands, all three match modes equal; 758 further mgaps inputs equal): the exact matches (mode 'mum' = maximal unique matches in both sequences, 'mumreference' = unique in the reference, 'maxmatch' = all maximal matches; minMatchLength = -l; only A/C/G/T match) are clustered: matches are sorted by their start in the query, matches inside a repeat are removed and overlapping matches of one diagonal merged, neighbours whose diagonals differ by at most max(diagonalDifference, diagonalFraction * separation) within maxSeparation are joined, and each component is cut into chains of maximum score (match length - overlap - diagonal difference); a chain is reported when its length (sum of the match lengths, or the extent with useExtent) reaches minClusterLength. Coordinates are 1-based as printed by mummer / mgaps; on the reverse strand the query start is a position in the reverse-complemented query; a match's overlap is the number of bases cut off its start because it overlaps the previous match. This is the stage that decides which anchors belong together; the banded extension (postnuc) and delta-filter are not part of this tool. Every option value mgaps accepts is accepted and gives the program's clusters (A4-4, measured on the installed mgaps 3.23: negative or zero separation / cluster length / diagonal difference, NaN, infinite or huge diagonal fraction - the C conversion of the product to long is LONG_MIN, so no fraction is added). One cluster of more than 100 000 matches is refused (SequenceAligner.NucmerMaxClusterMatches: the chain dynamic program of mgaps is quadratic in the cluster).")]
    public static NucmerClustersResult NucmerClusters(
        [Description("Reference DNA sequence.")] string reference,
        [Description("Query DNA sequence.")] string query,
        [Description("Minimum exact match length (nucmer -l, default 20).")] int minMatchLength = 20,
        [Description("Match mode: 'mum' (default), 'mumreference' or 'maxmatch'.")] string mode = "mum",
        [Description("Also cluster the reverse-complement strand of the query (nucmer default true).")] bool bothStrands = true,
        [Description("Largest separation between matches of a cluster (nucmer -g, default 90).")] long maxSeparation = 90,
        [Description("Smallest reported cluster length (nucmer -c, default 65).")] long minClusterLength = 65,
        [Description("Fixed diagonal difference that joins matches (nucmer -D, default 5).")] int diagonalDifference = 5,
        [Description("Fraction of the separation added to the allowed diagonal difference (nucmer -d, default 0.12).")] double diagonalFraction = 0.12,
        [Description("Cluster length = extent in the reference instead of the sum of the match lengths (mgaps -e, default false).")] bool useExtent = false,
        CancellationToken cancellationToken = default)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(query));
        var matchMode = (mode ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "mum" => CommonRegionMode.MaximalUnique,
            "mumreference" => CommonRegionMode.MaximalUniqueReference,
            "maxmatch" => CommonRegionMode.MaximalExactMatches,
            _ => throw new ArgumentException("mode must be 'mum', 'mumreference' or 'maxmatch'", nameof(mode)),
        };
        var options = new Seqeron.Genomics.Alignment.NucmerClusterOptions
        {
            MaxSeparation = maxSeparation, MinClusterLength = minClusterLength, DiagonalDifference = diagonalDifference, DiagonalFraction = diagonalFraction, UseExtent = useExtent,
        };
        var clusters = GenomicAnalyzer.FindNucmerClusters(reference, query, minMatchLength, matchMode, options, bothStrands ? StrandSelection.Both : StrandSelection.Forward, cancellationToken)
            .Select(c => new NucmerClusterItem(c.IsReverseComplement ? "reverse" : "forward", c.Cluster.Matches.Select(m => new NucmerClusterMatchItem(m.Start1, m.Start2, m.Length, m.Overlap)).ToArray()))
            .ToArray();
        return new NucmerClustersResult(clusters);
    }
}
