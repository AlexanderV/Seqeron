using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A3-9c: the chaining step of minimap2 2.26; delegates to Minimap2Chaining.
public partial class AlignmentTools
{
    [McpServerTool(Name = "minimap2_chains", Title = "Pairwise — minimap2 Seed Chaining", ReadOnly = true)]
    [Description("The seeding and chaining stage of minimap2 2.26 (450 related pairs: the seed array of `--print-seeds` equal anchor by anchor, the chains of `--print-chains` equal, chain scores = the PAF s1 / cm tags; the chaining dynamic program mg_lchain_dp, its backtracking and compaction are ported line by line): the (k, w)-minimizers of the targets and of the query (mm_sketch, k <= 28, w < 256; symmetric k-mers skipped, A/C/G/T/U only) are indexed, every query minimizer that occurs in a target gives anchors on the forward or reverse strand (target / query position of the LAST base of the seed, sorted by minimap2's in-place radix sort), and the anchors are chained: the score of a step is min(min(dq, dr), span) minus floor(0.8 * 0.01 * k * |dr - dq| + 0.5 * log2(|dr - dq| + 1)); a predecessor must lie within maxGap on both sequences and within bandwidth of the diagonal, at most 5000 predecessors are scanned per anchor and the scan stops after 25 non-improving ones; chains are cut where the score drops by more than the bandwidth below its peak and kept when they reach minChainScore and minAnchors. Returns the chains (score, target, strand, anchors in increasing order) ordered by the target position of their first anchor. Frequent minimizers (A4-7): a minimizer that occurs more than mid_occ times in the targets (mid_occ = max(minMidOcc, the (1 - midOccFraction) quantile of the occurrence counts + 1), capped by maxMidOcc; midOccFraction >= 1 is the limit itself, as minimap2 -f) is dropped except for the least frequent ones of a streak of frequent minimizers (mm_seed_select: one per occDist query bases, never above 4095 occurrences); query minimizers that occur more than mid_occ times in the query and in more than queryOccFraction of its minimizers are dropped first; repeatLength = the query bases covered by dropped minimizers (the PAF rl tag). Also ported: homopolymerCompressed (-H), sdustThreshold (-T, the sdust low-complexity masking of the query minimizers), allVsAll (-X: with targetNames / queryName, hits on the diagonal of a same-named, same-length target and targets that sort before the query are skipped; all chains are reported, no long join), forwardOnly / reverseOnly. mode = chains (default, as before: every chain after the dynamic program, ordered by the first anchor) or paf: the whole mapping up to the end of the chaining stage (mm_map_frag without -c): the long-join re-chaining of -r b,B (bandwidthLong > bandwidth: mg_lchain_rmq, a range-minimum-query AVL tree, when the first chains do not cover the query), rmq chaining (--rmq), the primary / secondary marking (mm_set_parent: maskLevel, maskLength, hardMaskLevel; mm_select_sub: primaryRatio -p, bestN -N, strand retention), the divergence estimate dv (mm_est_err) and the mapping quality, i.e. every field of the PAF lines (tp P / S, cm, s1, s2, mapq, dv, rl) of minimap2 without -c (6 000+ oracle cases of the installed program: seeds in six groups of 330 cases, PAF lines in four groups of 330). Not part of this tool: the base-level alignment minimap2 -c adds (ksw2 gap filling and extension, registered A3-9c2), splice mode and the short-read / paired-end modes. Options are validated like the program (A4-4): k outside 1-28 is refused (minimap2: \"-k and -w must be positive\", assertion k <= 28 in mm_sketch); every other value minimap2 accepts is accepted and gives the program's chains - a negative bandwidth leaves no chain, a negative gap limit means the bandwidth, minAnchors 0 / -1 and a negative minChainScore are no limit (measured on 30 option values, see option_boundary_oracle.json). A query / target pair whose minimizers have more than 8 388 608 anchors in all (SequenceAligner.MinimapMaxAnchors, e.g. poly-A or a short tandem repeat against itself) is refused with an ArgumentException before the anchors are stored; minimap2 itself dies at the allocation or drops the frequent minimizers (not ported).")]
    public static Minimap2ChainsResult Minimap2Chains(
        [Description("Target (reference) sequences (A/C/G/T/U; other symbols break k-mers), at least one.")] string[] targets,
        [Description("Query sequence.")] string query,
        [Description("k-mer size (minimap2 -k, 1-28, default 15).")] int k = 15,
        [Description("Minimizer window size (minimap2 -w, 1-255, default 10).")] int w = 10,
        [Description("Largest distance between two anchors of a chain on both sequences (minimap2 -g, default 5000).")] int maxGap = 5000,
        [Description("Chaining bandwidth (minimap2 -r, first value, default 500).")] int bandwidth = 500,
        [Description("Smallest chain score (minimap2 -m, default 40).")] int minChainScore = 40,
        [Description("Fewest anchors of a chain (minimap2 -n, default 3).")] int minAnchors = 3,
        [Description("chains (default: every chain after the dynamic program, ordered by the first anchor) | paf (the mapping up to the end of the chaining stage: long-join re-chaining, primary / secondary marking, dv, mapq; the hits are returned in the PAF order).")] string mode = "chains",
        [Description("Homopolymer-compressed k-mers (minimap2 -H).")] bool homopolymerCompressed = false,
        [Description("sdust threshold masking low-complexity query minimizers (minimap2 -T, hidden option; 0 = off, default).")] int sdustThreshold = 0,
        [Description("minimap2 -f: below 1 the fraction of repeated minimizers to filter (default 0.0002), otherwise the occurrence limit mid_occ itself.")] double midOccFraction = 0.0002,
        [Description("minimap2 -U lower bound of mid_occ (default 10).")] int minMidOcc = 10,
        [Description("minimap2 -U upper bound of mid_occ (default 1000000).")] int maxMidOcc = 1000000,
        [Description("minimap2 --q-occ-frac (default 0.01).")] double queryOccFraction = 0.01,
        [Description("minimap2 -e: occ_dist, one frequent minimizer per this many query bases is kept (default 500; 0 = none).")] int occDist = 500,
        [Description("All-vs-all mode (minimap2 -X): needs targetNames and queryName; all chains are reported without primary selection and without long join.")] bool allVsAll = false,
        [Description("Names of the targets (for allVsAll), one per target.")] string[]? targetNames = null,
        [Description("Name of the query (for allVsAll; in paf mode its hash orders chains of equal score).")] string? queryName = null,
        [Description("Only hits on the forward strand (minimap2 --for-only).")] bool forwardOnly = false,
        [Description("Only hits on the reverse strand (minimap2 --rev-only).")] bool reverseOnly = false,
        [Description("Long-join bandwidth (minimap2 -r second value, default 20000; the re-chaining runs only when it exceeds bandwidth, paf mode).")] int bandwidthLong = 20000,
        [Description("Stop the search for a predecessor after this many non-improving anchors (minimap2 --max-chain-skip, default 25).")] int maxChainSkip = 25,
        [Description("Largest number of predecessors examined per anchor (minimap2 --max-chain-iter, default 5000).")] int maxChainIterations = 5000,
        [Description("Gap penalty scale (minimap2 --chain-gap-scale, default 0.8).")] double chainGapScale = 0.8,
        [Description("Skip penalty scale (minimap2 --chain-skip-scale, default 0).")] double chainSkipScale = 0,
        [Description("paf mode: retain at most this many secondary hits (minimap2 -N, default 5).")] int bestN = 5,
        [Description("paf mode: minimum secondary-to-primary score ratio (minimap2 -p, default 0.8; 0 = no selection).")] double primaryRatio = 0.8,
        [Description("paf mode: minimap2 --mask-level (default 0.5).")] double maskLevel = 0.5,
        [Description("paf mode: minimap2 --mask-len (default none).")] int maskLength = int.MaxValue,
        [Description("paf mode: minimap2 --hard-mask-level.")] bool hardMaskLevel = false,
        [Description("paf mode: no long-join re-chaining (minimap2 --no-long-join).")] bool noLongJoin = false,
        [Description("paf mode: chain with the RMQ algorithm (minimap2 --rmq=yes).")] bool rmq = false,
        [Description("paf mode: return the secondary hits too (default true; false = minimap2 --secondary=no).")] bool outputSecondary = true,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        if (targets is null || targets.Length == 0)
            throw new ArgumentException("At least one target sequence is required", nameof(targets));
        foreach (var target in targets)
            RequireSequence(target, nameof(targets));
        RequireSequence(query, nameof(query));
        if (mode is null || (!mode.Equals("chains", StringComparison.OrdinalIgnoreCase) && !mode.Equals("paf", StringComparison.OrdinalIgnoreCase)))
            throw new ArgumentException("mode must be chains or paf", nameof(mode));
        if (allVsAll && (targetNames is null || targetNames.Length != targets.Length || queryName is null))
            throw new ArgumentException("allVsAll needs one target name per target and the query name", nameof(targetNames));
        // minimap2 main.c: -f x[,y]: x < 1 is the fraction, otherwise the occurrence limit itself
        var seed = new global::Seqeron.Genomics.Alignment.MinimapSeedOptions
        {
            KmerSize = k, WindowSize = w, HomopolymerCompressed = homopolymerCompressed, SdustThreshold = sdustThreshold,
            MidOccFraction = midOccFraction >= 1.0 ? 2e-4f : (float)midOccFraction, MidOcc = midOccFraction >= 1.0 ? (int)(midOccFraction + .499) : 0,
            MinMidOcc = minMidOcc, MaxMidOcc = maxMidOcc, QueryOccurrenceFraction = (float)queryOccFraction, OccDist = occDist,
            SkipSelf = allVsAll, SkipDual = allVsAll, ForwardOnly = forwardOnly, ReverseOnly = reverseOnly,
        };
        if (mode.Equals("paf", StringComparison.OrdinalIgnoreCase))
        {
            var mapOptions = new global::Seqeron.Genomics.Alignment.MinimapMapOptions
            {
                Seed = seed, MaxGap = maxGap, Bandwidth = bandwidth, BandwidthLong = bandwidthLong, MaxChainSkip = maxChainSkip, MaxChainIterations = maxChainIterations, MinCount = minAnchors,
                MinChainScore = minChainScore, ChainGapScale = (float)chainGapScale, ChainSkipScale = (float)chainSkipScale, NoLongJoin = noLongJoin || allVsAll, UseRmq = rmq, AllChains = allVsAll,
                MaskLevel = (float)maskLevel, MaskLength = maskLength, HardMaskLevel = hardMaskLevel, PrimaryRatio = (float)primaryRatio, BestN = bestN, OutputSecondary = outputSecondary,
            };
            var mapped = global::Seqeron.Genomics.Alignment.Minimap2Chaining.Map(targets, query, mapOptions, targetNames, queryName);
            var hits = mapped.Hits.Select(h => new Minimap2HitDto(h.Id, h.Parent, h.IsPrimary, h.TargetId, h.Reverse ? "reverse" : "forward", h.QueryStart, h.QueryEnd, h.TargetStart, h.TargetEnd, h.MatchLength,
                h.BlockLength, h.MapQ, h.Score, h.SecondaryScore, h.AnchorCount, h.Divergence >= 0f ? Math.Round(h.Divergence, 4) : null,
                h.Anchors.Select(a => new Minimap2AnchorDto(a.TargetId, a.Reverse ? "reverse" : "forward", a.TargetEnd, a.QueryEnd, a.Span)).ToArray())).ToArray();
            var asChains = mapped.Hits.Select(h => new Minimap2ChainDto(h.Score, h.TargetId, h.Reverse ? "reverse" : "forward",
                h.Anchors.Select(a => new Minimap2AnchorDto(a.TargetId, a.Reverse ? "reverse" : "forward", a.TargetEnd, a.QueryEnd, a.Span)).ToArray())).ToArray();
            return new Minimap2ChainsResult(asChains, hits.Sum(h => h.Anchors.Length), hits, mapped.RepeatLength, mapped.MidOcc);
        }

        var options = new global::Seqeron.Genomics.Alignment.MinimapChainOptions
        {
            MaxTargetGap = maxGap, MaxQueryGap = maxGap, Bandwidth = bandwidth, MinScore = minChainScore, MinCount = minAnchors, MaxChainSkip = maxChainSkip, MaxChainIterations = maxChainIterations,
            ChainGapScale = (float)chainGapScale, ChainSkipScale = (float)chainSkipScale,
        };
        var seeds = global::Seqeron.Genomics.Alignment.Minimap2Chaining.CollectSeeds(targets, query, seed, targetNames, queryName);
        var anchors = seeds.Anchors;
        var chains = global::Seqeron.Genomics.Alignment.Minimap2Chaining.Chain(anchors, options with { KmerSize = k })
            .Select(c => new Minimap2ChainDto(c.Score, anchors[c.Anchors[0]].TargetId, anchors[c.Anchors[0]].Reverse ? "reverse" : "forward",
                c.Anchors.Select(i => anchors[i]).Select(a => new Minimap2AnchorDto(a.TargetId, a.Reverse ? "reverse" : "forward", a.TargetEnd, a.QueryEnd, a.Span)).ToArray()))
            .ToArray();
        return new Minimap2ChainsResult(chains, anchors.Count, null, seeds.RepeatLength, seeds.MidOcc);
    }
}
