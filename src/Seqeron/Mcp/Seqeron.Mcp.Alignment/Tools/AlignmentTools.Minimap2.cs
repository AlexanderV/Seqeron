using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A3-9c: the chaining step of minimap2 2.26; delegates to Minimap2Chaining.
public partial class AlignmentTools
{
    [McpServerTool(Name = "minimap2_chains", Title = "Pairwise — minimap2 Seed Chaining", ReadOnly = true)]
    [Description("The seeding and chaining stage of minimap2 2.26 (450 related pairs: the seed array of `--print-seeds` equal anchor by anchor, the chains of `--print-chains` equal, chain scores = the PAF s1 / cm tags; the chaining dynamic program mg_lchain_dp, its backtracking and compaction are ported line by line): the (k, w)-minimizers of the targets and of the query (mm_sketch, k <= 28, w < 256; symmetric k-mers skipped, A/C/G/T/U only) are indexed, every query minimizer that occurs in a target gives anchors on the forward or reverse strand (target / query position of the LAST base of the seed, sorted by minimap2's in-place radix sort), and the anchors are chained: the score of a step is min(min(dq, dr), span) minus floor(0.8 * 0.01 * k * |dr - dq| + 0.5 * log2(|dr - dq| + 1)); a predecessor must lie within maxGap on both sequences and within bandwidth of the diagonal, at most 5000 predecessors are scanned per anchor and the scan stops after 25 non-improving ones; chains are cut where the score drops by more than the bandwidth below its peak and kept when they reach minChainScore and minAnchors. Returns the chains (score, target, strand, anchors in increasing order) ordered by the target position of their first anchor. Not part of this tool: the base-level alignment minimap2 -c adds (ksw2 gap filling and extension), the secondary / primary selection, and the handling of frequent minimizers (a minimizer of the query that occurs more than mid_occ = max(10, 0.9998 quantile + 1) times in the targets is an error here, as is homopolymer compression).")]
    public static Minimap2ChainsResult Minimap2Chains(
        [Description("Target (reference) sequences (A/C/G/T/U; other symbols break k-mers), at least one.")] string[] targets,
        [Description("Query sequence.")] string query,
        [Description("k-mer size (minimap2 -k, 1-28, default 15).")] int k = 15,
        [Description("Minimizer window size (minimap2 -w, 1-255, default 10).")] int w = 10,
        [Description("Largest distance between two anchors of a chain on both sequences (minimap2 -g, default 5000).")] int maxGap = 5000,
        [Description("Chaining bandwidth (minimap2 -r, first value, default 500).")] int bandwidth = 500,
        [Description("Smallest chain score (minimap2 -m, default 40).")] int minChainScore = 40,
        [Description("Fewest anchors of a chain (minimap2 -n, default 3).")] int minAnchors = 3)
    {
        if (targets is null || targets.Length == 0)
            throw new ArgumentException("At least one target sequence is required", nameof(targets));
        foreach (var target in targets)
            RequireSequence(target, nameof(targets));
        RequireSequence(query, nameof(query));
        var options = new global::Seqeron.Genomics.Alignment.MinimapChainOptions
        {
            MaxTargetGap = maxGap, MaxQueryGap = maxGap, Bandwidth = bandwidth, MinScore = minChainScore, MinCount = minAnchors,
        };
        var anchors = global::Seqeron.Genomics.Alignment.Minimap2Chaining.CollectAnchors(targets, query, k, w);
        var chains = global::Seqeron.Genomics.Alignment.Minimap2Chaining.Chain(anchors, options with { KmerSize = k })
            .Select(c => new Minimap2ChainDto(c.Score, anchors[c.Anchors[0]].TargetId, anchors[c.Anchors[0]].Reverse ? "reverse" : "forward",
                c.Anchors.Select(i => anchors[i]).Select(a => new Minimap2AnchorDto(a.TargetId, a.Reverse ? "reverse" : "forward", a.TargetEnd, a.QueryEnd, a.Span)).ToArray()))
            .ToArray();
        return new Minimap2ChainsResult(chains, anchors.Count);
    }
}
