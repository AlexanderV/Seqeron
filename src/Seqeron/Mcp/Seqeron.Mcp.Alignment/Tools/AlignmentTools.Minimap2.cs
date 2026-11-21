using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A3-9c: the chaining step of minimap2 2.26; delegates to Minimap2Chaining.
public partial class AlignmentTools
{
    [McpServerTool(Name = "minimap2_chains", Title = "Pairwise — minimap2 Seed Chaining", ReadOnly = true)]
    [Description("The seeding and chaining stage of minimap2 2.26 (450 related pairs: the seed array of `--print-seeds` equal anchor by anchor, the chains of `--print-chains` equal, chain scores = the PAF s1 / cm tags; the chaining dynamic program mg_lchain_dp, its backtracking and compaction are ported line by line): the (k, w)-minimizers of the targets and of the query (mm_sketch, k <= 28, w < 256; symmetric k-mers skipped, A/C/G/T/U only) are indexed, every query minimizer that occurs in a target gives anchors on the forward or reverse strand (target / query position of the LAST base of the seed, sorted by minimap2's in-place radix sort), and the anchors are chained: the score of a step is min(min(dq, dr), span) minus floor(0.8 * 0.01 * k * |dr - dq| + 0.5 * log2(|dr - dq| + 1)); a predecessor must lie within maxGap on both sequences and within bandwidth of the diagonal, at most 5000 predecessors are scanned per anchor and the scan stops after 25 non-improving ones; chains are cut where the score drops by more than the bandwidth below its peak and kept when they reach minChainScore and minAnchors. Returns the chains (score, target, strand, anchors in increasing order) ordered by the target position of their first anchor. Frequent minimizers (A4-7): a minimizer that occurs more than mid_occ times in the targets (mid_occ = max(minMidOcc, the (1 - midOccFraction) quantile of the occurrence counts + 1), capped by maxMidOcc; midOccFraction >= 1 is the limit itself, as minimap2 -f) is dropped except for the least frequent ones of a streak of frequent minimizers (mm_seed_select: one per occDist query bases, never above 4095 occurrences); query minimizers that occur more than mid_occ times in the query and in more than queryOccFraction of its minimizers are dropped first; repeatLength = the query bases covered by dropped minimizers (the PAF rl tag). Also ported: homopolymerCompressed (-H), sdustThreshold (-T, the sdust low-complexity masking of the query minimizers), allVsAll (-X: with targetNames / queryName, hits on the diagonal of a same-named, same-length target and targets that sort before the query are skipped; all chains are reported, no long join), forwardOnly / reverseOnly. mode = chains (default, as before: every chain after the dynamic program, ordered by the first anchor) or paf: the whole mapping up to the end of the chaining stage (mm_map_frag without -c): the long-join re-chaining of -r b,B (bandwidthLong > bandwidth: mg_lchain_rmq, a range-minimum-query AVL tree, when the first chains do not cover the query), rmq chaining (--rmq), the primary / secondary marking (mm_set_parent: maskLevel, maskLength, hardMaskLevel; mm_select_sub: primaryRatio -p, bestN -N, strand retention), the divergence estimate dv (mm_est_err) and the mapping quality, i.e. every field of the PAF lines (tp P / S, cm, s1, s2, mapq, dv, rl) of minimap2 without -c (6 000+ oracle cases of the installed program: seeds in six groups of 330 cases, PAF lines in four groups of 330). Not part of this tool: the base-level alignment minimap2 -c adds (ksw2 gap filling and extension), the SAM output, splice mode and the short-read / paired-end modes: use minimap2_map for the whole program. Options are validated like the program (A4-4): k outside 1-28 is refused (minimap2: \"-k and -w must be positive\", assertion k <= 28 in mm_sketch); every other value minimap2 accepts is accepted and gives the program's chains - a negative bandwidth leaves no chain, a negative gap limit means the bandwidth, minAnchors 0 / -1 and a negative minChainScore are no limit (measured on 30 option values, see option_boundary_oracle.json). A query / target pair whose minimizers have more than 8 388 608 anchors in all (SequenceAligner.MinimapMaxAnchors, e.g. poly-A or a short tandem repeat against itself) is refused with an ArgumentException before the anchors are stored; minimap2 itself dies at the allocation or drops the frequent minimizers (not ported).")]
    public static Minimap2ChainsResult Minimap2Chains(
        [Description("Target (reference) sequences (A/C/G/T/U; other symbols break k-mers), at least one.")] string[] targets,
        [Description("Query sequence.")] string query,
        [Description("k-mer size (minimap2 -k, 1-28; omit = 15, or the preset's value).")] int? k = null,
        [Description("Minimizer window size (minimap2 -w, 1-255; omit = 10, or the preset's value).")] int? w = null,
        [Description("Largest distance between two anchors of a chain on both sequences (minimap2 -g; omit = 5000, or the preset's value).")] int? maxGap = null,
        [Description("Chaining bandwidth (minimap2 -r, first value; omit = 500, or the preset's value).")] int? bandwidth = null,
        [Description("Smallest chain score (minimap2 -m; omit = 40, or the preset's value).")] int? minChainScore = null,
        [Description("Fewest anchors of a chain (minimap2 -n; omit = 3, or the preset's value).")] int? minAnchors = null,
        [Description("chains (default: every chain after the dynamic program, ordered by the first anchor) | paf (the mapping up to the end of the chaining stage: long-join re-chaining, primary / secondary marking, dv, mapq; the hits are returned in the PAF order).")] string mode = "chains",
        [Description("Homopolymer-compressed k-mers (minimap2 -H; the flag only switches on: false = not given, a preset that sets it keeps it).")] bool homopolymerCompressed = false,
        [Description("sdust threshold masking low-complexity query minimizers (minimap2 -T, hidden option; omit = 0 = off).")] int? sdustThreshold = null,
        [Description("minimap2 -f: below 1 the fraction of repeated minimizers (omit = 0.0002), otherwise the occurrence limit mid_occ itself.")] double? midOccFraction = null,
        [Description("minimap2 -U lower bound of mid_occ (omit = 10, or the preset's value).")] int? minMidOcc = null,
        [Description("minimap2 -U upper bound of mid_occ (omit = 1000000, or the preset's value).")] int? maxMidOcc = null,
        [Description("minimap2 --q-occ-frac (omit = 0.01).")] double? queryOccFraction = null,
        [Description("minimap2 -e: occ_dist, one frequent minimizer per this many query bases is kept (omit = 500, or the preset's value; 0 = none).")] int? occDist = null,
        [Description("All-vs-all mode (minimap2 -X): needs targetNames and queryName; all chains are reported without primary selection and without long join.")] bool allVsAll = false,
        [Description("Names of the targets (for allVsAll), one per target.")] string[]? targetNames = null,
        [Description("Name of the query (for allVsAll; in paf mode its hash orders chains of equal score).")] string? queryName = null,
        [Description("Only hits on the forward strand (minimap2 --for-only).")] bool forwardOnly = false,
        [Description("Only hits on the reverse strand (minimap2 --rev-only).")] bool reverseOnly = false,
        [Description("Long-join bandwidth (minimap2 -r second value; omit = 20000, or the preset's value; the re-chaining runs only when it exceeds bandwidth, paf mode).")] int? bandwidthLong = null,
        [Description("Stop the search for a predecessor after this many non-improving anchors (minimap2 --max-chain-skip; omit = 25, or the preset's value).")] int? maxChainSkip = null,
        [Description("Largest number of predecessors examined per anchor (minimap2 --max-chain-iter; omit = 5000).")] int? maxChainIterations = null,
        [Description("Gap penalty scale (minimap2 --chain-gap-scale; omit = 0.8).")] double? chainGapScale = null,
        [Description("Skip penalty scale (minimap2 --chain-skip-scale; omit = 0).")] double? chainSkipScale = null,
        [Description("paf mode: retain at most this many secondary hits (minimap2 -N; omit = 5, or the preset's value).")] int? bestN = null,
        [Description("paf mode: minimum secondary-to-primary score ratio (minimap2 -p; omit = 0.8, or the preset's value; 0 = no selection).")] double? primaryRatio = null,
        [Description("paf mode: minimap2 --mask-level (omit = 0.5).")] double? maskLevel = null,
        [Description("paf mode: minimap2 --mask-len (omit = none).")] int? maskLength = null,
        [Description("paf mode: minimap2 --hard-mask-level (a flag: false = not given).")] bool hardMaskLevel = false,
        [Description("paf mode: no long-join re-chaining (minimap2 --no-long-join; a flag: false = not given, a preset that sets it keeps it).")] bool noLongJoin = false,
        [Description("paf mode: chain with the RMQ algorithm (minimap2 --rmq=yes | --rmq=no; omit = off, or the preset's value: --rmq=no switches the asm presets' RMQ off).")] bool? rmq = null,
        [Description("paf mode: return the secondary hits (minimap2 --secondary=yes | --secondary=no; omit = yes, or the preset's value: --secondary=yes re-enables them for -x sr).")] bool? outputSecondary = null,
        CancellationToken cancellationToken = default,
        [Description("paf mode: the minimap2 -x preset (A5-9): map-ont, map-pb / map10k, map-hifi / map-ccs, ava-ont, ava-pb, asm5, asm10, asm20, sr / short, splice / cdna, splice:hq (options.c mm_set_opt). The program applies the preset first and every explicit option after it, so does this tool: every parameter that is GIVEN overrides the preset, also a value equal to the default (preset asm5 with k = 15 is k = 15, not the preset's 19); omitted parameters take the preset's value. The base-level scoring of the preset is not used (no -c); the chains mode takes no preset.")] string? preset = null,
        [Description("paf mode: minimap2 -F, the maximum fragment length (effective with preset sr, which sets 800; omit = the preset's value).")] int? maxFragmentLength = null,
        [Description("paf mode: minimap2 -G / --max-intron-len: with a splice preset the longest gap on the reference, the bandwidth and the long bandwidth (ignored otherwise, as the program does).")] int? maxIntronLength = null)
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
        // minimap2 applies -x first and every option after it: the one canonical merge is MinimapPresets.Apply(preset, overrides) (A6-7: a given value overrides the preset also when it equals the default)
        var overrides = new global::Seqeron.Genomics.Alignment.MinimapOverrides
        {
            KmerSize = k, WindowSize = w, HomopolymerCompressed = homopolymerCompressed, SdustThreshold = sdustThreshold, MidOccFraction = midOccFraction, MinMidOcc = minMidOcc, MaxMidOcc = maxMidOcc,
            QueryOccurrenceFraction = queryOccFraction, OccDist = occDist, ForwardOnly = forwardOnly, ReverseOnly = reverseOnly, AllVsAll = allVsAll, MaxGap = maxGap, Bandwidth = bandwidth, BandwidthLong = bandwidthLong,
            MaxChainSkip = maxChainSkip, MaxChainIterations = maxChainIterations, MinCount = minAnchors, MinChainScore = minChainScore, ChainGapScale = chainGapScale, ChainSkipScale = chainSkipScale,
            NoLongJoin = noLongJoin, UseRmq = rmq, MaskLevel = maskLevel, MaskLength = maskLength, HardMaskLevel = hardMaskLevel, PrimaryRatio = primaryRatio, BestN = bestN, OutputSecondary = outputSecondary,
            MaxFragmentLength = maxFragmentLength, MaxIntronLength = maxIntronLength,
        };
        if (mode.Equals("paf", StringComparison.OrdinalIgnoreCase))
        {
            var mapOptions = global::Seqeron.Genomics.Alignment.MinimapPresets.Apply(preset, overrides);
            var mapped = global::Seqeron.Genomics.Alignment.Minimap2Chaining.Map(targets, query, mapOptions, targetNames, queryName);
            string[] pafNames = targetNames ?? [.. targets.Select((_, i) => "target" + (i + 1).ToString(System.Globalization.CultureInfo.InvariantCulture))];
            string pafText = global::Seqeron.Genomics.Alignment.Minimap2Chaining.FormatPaf(mapped, queryName ?? "query", query.Length, pafNames, [.. targets.Select(t => t.Length)]);
            var hits = mapped.Hits.Select(h => new Minimap2HitDto(h.Id, h.Parent, h.IsPrimary, h.TargetId, h.Reverse ? "reverse" : "forward", h.QueryStart, h.QueryEnd, h.TargetStart, h.TargetEnd, h.MatchLength,
                h.BlockLength, h.MapQ, h.Score, h.SecondaryScore, h.AnchorCount, h.Divergence >= 0f ? Math.Round(h.Divergence, 4) : null,
                h.Anchors.Select(a => new Minimap2AnchorDto(a.TargetId, a.Reverse ? "reverse" : "forward", a.TargetEnd, a.QueryEnd, a.Span)).ToArray())).ToArray();
            var asChains = mapped.Hits.Select(h => new Minimap2ChainDto(h.Score, h.TargetId, h.Reverse ? "reverse" : "forward",
                h.Anchors.Select(a => new Minimap2AnchorDto(a.TargetId, a.Reverse ? "reverse" : "forward", a.TargetEnd, a.QueryEnd, a.Span)).ToArray())).ToArray();
            return new Minimap2ChainsResult(asChains, hits.Sum(h => h.Anchors.Length), hits, mapped.RepeatLength, mapped.MidOcc, pafText);
        }

        var resolved = global::Seqeron.Genomics.Alignment.MinimapPresets.Apply(null, overrides);        // the chains mode takes no preset: the defaults and the given options
        var seed = resolved.Seed;
        var options = new global::Seqeron.Genomics.Alignment.MinimapChainOptions
        {
            MaxTargetGap = resolved.MaxGap, MaxQueryGap = resolved.MaxGap, Bandwidth = resolved.Bandwidth, MinScore = resolved.MinChainScore, MinCount = resolved.MinCount, MaxChainSkip = resolved.MaxChainSkip,
            MaxChainIterations = resolved.MaxChainIterations, ChainGapScale = resolved.ChainGapScale, ChainSkipScale = resolved.ChainSkipScale,
        };
        var seeds = global::Seqeron.Genomics.Alignment.Minimap2Chaining.CollectSeeds(targets, query, seed, targetNames, queryName);
        var anchors = seeds.Anchors;
        var chains = global::Seqeron.Genomics.Alignment.Minimap2Chaining.Chain(anchors, options with { KmerSize = seed.KmerSize })
            .Select(c => new Minimap2ChainDto(c.Score, anchors[c.Anchors[0]].TargetId, anchors[c.Anchors[0]].Reverse ? "reverse" : "forward",
                c.Anchors.Select(i => anchors[i]).Select(a => new Minimap2AnchorDto(a.TargetId, a.Reverse ? "reverse" : "forward", a.TargetEnd, a.QueryEnd, a.Span)).ToArray())).ToArray();
        return new Minimap2ChainsResult(chains, anchors.Count, null, seeds.RepeatLength, seeds.MidOcc);
    }

    [McpServerTool(Name = "minimap2_map", Title = "Pairwise — minimap2 Program (mapping, -c / -a alignment, paired reads)", ReadOnly = true)]
    [Description("minimap2 2.26 as a program: the whole standard output (PAF or SAM) of `minimap2 [options] ref.fa query.fa [query2.fa]` on the given sequences (every option of the mapping and alignment stages, parsed with the program's getopt rules, `-x` applied first and cumulatively, the option checks of mm_check_opt with the program's exit status and error text). reference = the FASTA text of the target sequences (named ref.fa on the command line), queries = the FASTA / FASTQ text of one query file (query.fa; fragment mode = interleaved reads with the same name, /1 /2 suffixes ignored) or of two (query1.fa, query2.fa: the paired reads of `-x sr`), options = the command-line words before the files, e.g. [\"-x\", \"sr\", \"-a\"] or [\"-c\", \"--cs\", \"--MD\", \"-N\", \"10\"]. Base-level alignment (`-c`, `-a`, `--cs`, `--MD`: align.c, ksw2 gap filling between the chained anchors, left / right extension, Z-drop, inversion rescue, spliced alignment with `-x splice` / `-u` / `--junc-bed`, the CIGAR `cg`, `cs` short / long, `MD`, `--eqx`, `-Y`, `-L`, `--secondary-seq`, tags NM ms AS nn ts de zd), the SAM writer (flags, MAPQ, SA, mate fields, @HD / @SQ / @RG / @PG header), the multi-segment paired-read mode (segment ids in the anchors, mm_select_sub_multi, mm_seg_gen, mm_pair with pe_bonus, proper pairs), the index in parts (`-I`) and the `--split-prefix` merge, `--alt` / `--alt-drop`, `--qstrand`, `-X` / `-P`, `--paf-no-hit`, `-o FILE` (the output is returned in files) and the debug-only options (`-t -K --print-seeds --print-chains --print-aln-seq --no-kalloc --cap-kalloc --bucket-bits`, no effect on the output) are ported from the 2.26 source and compared with the installed program as text: 2 700 base-level PAF / SAM cases, 2 100 paired cases, 2 403 option / error cases (exit status and standard error included), all equal; the ksw2 kernels ksw_extz2_sse / ksw_extd2_sse / ksw_exts2_sse / ksw_ll_i16 equal the C functions compiled with SSE4.1 on 2 000 random inputs. The program's own log lines (`[M::...]`) and `-d` (index dump) are not reproduced; -t has no effect; a k above 28 or a window of 256 or more is refused (the C program aborts).")]
    public static Minimap2MapResult Minimap2Map(
        [Description("The FASTA text of the target (reference) sequences.")] string reference,
        [Description("The FASTA / FASTQ text of one query file, or of two query files (paired reads, `-x sr`).")] string[] queries,
        [Description("The command-line options before the files (e.g. [\"-x\", \"map-ont\", \"-c\", \"--cs\"]); omit for `minimap2 ref.fa query.fa`.")] string[]? options = null,
        [Description("The text of a BED file of annotated junctions for `--junc-bed` (BED12 blocks or BED6 intervals), or null.")] string? junctionBed = null,
        [Description("The text of an ALT contig list for `--alt` (the first word of each line is a target name), or null.")] string? altList = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        if (string.IsNullOrWhiteSpace(reference))
            throw new ArgumentException("A reference is required", nameof(reference));
        if (queries is null || queries.Length == 0 || queries.Length > 3)
            throw new ArgumentException("One query file (or two for paired reads) is required", nameof(queries));
        foreach (var q in queries)
            if (string.IsNullOrEmpty(q))
                throw new ArgumentException("A query file is empty", nameof(queries));
        var files = new Dictionary<string, string> { ["ref.fa"] = reference };
        var args = new List<string>(options ?? []);
        if (args.Any(a => a is null))
            throw new ArgumentException("An option is null", nameof(options));
        if (junctionBed is not null)
        {
            files["junctions.bed"] = junctionBed;
            args.AddRange(["--junc-bed", "junctions.bed"]);
        }

        if (altList is not null)
        {
            files["alt.txt"] = altList;
            args.AddRange(["--alt", "alt.txt"]);
        }

        args.Add("ref.fa");
        for (int i = 0; i < queries.Length; i++)
        {
            string name = queries.Length == 1 ? "query.fa" : $"query{i + 1}.fa";
            files[name] = queries[i];
            args.Add(name);
        }

        try
        {
            var r = global::Seqeron.Genomics.Alignment.Minimap2Program.Run(args, files);
            return new Minimap2MapResult(r.StdOut, r.StdErr, r.ExitCode, r.Files.Count == 0 ? null : new Dictionary<string, string>(r.Files));
        }
        catch (ArgumentOutOfRangeException ex)
        {
            throw new ArgumentException(ex.Message, nameof(options), ex);
        }
    }
}
