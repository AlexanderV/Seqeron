namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The <c>-x</c> presets of minimap2 2.26 (<c>options.c</c> <c>mm_set_opt</c>, read at raw.githubusercontent.com/lh3/minimap2 v2.26) as <see cref="MinimapMapOptions"/>: <c>map-ont</c> (the defaults), <c>map-pb</c> / <c>map10k</c>, <c>map-hifi</c> /
/// <c>map-ccs</c>, <c>ava-ont</c>, <c>ava-pb</c>, <c>asm5</c> / <c>asm10</c> / <c>asm20</c>, <c>sr</c> / <c>short</c>, <c>splice</c> / <c>cdna</c> and <c>splice:hq</c>. The program applies <c>-x</c> BEFORE every other option, so an explicit
/// option overrides the preset: <c>MinimapPresets.Apply("asm5") with { MinChainScore = 30 }</c> is <c>minimap2 -x asm5 -m 30</c>.
/// </summary>
/// <remarks>
/// Explicit options of the command line that FOLLOW <c>-x</c> are <see cref="MinimapOverrides"/>: <c>MinimapPresets.Apply("asm5", new MinimapOverrides { KmerSize = 15 })</c> is <c>minimap2 -x asm5 -k 15</c>, also when the value equals the default of the
/// option (the MCP tool <c>minimap2_chains</c> and the oracle test build the options this way: audit A6-7). The chaining-stage members (<c>-k -w -H -f/occ -g -r -n -m -p -N -S ... --rmq -P -D --no-dual --no-long-join -F --sr</c>) are read by <see cref="Minimap2Chaining.Map"/>; the base-level members (<c>-A -B -O -E -z -s -u</c>, splice flags,
/// <c>pe_ori</c>, the fragment mode, the heap sort) are read by the base-level alignment (<see cref="MinimapMapOptions.Cigar"/>, audit A3-9c2, F187) and the multi-segment mode (<see cref="Minimap2Chaining.MapFragment"/>, A5-9b / A5-9c, F188). <see cref="ApplyOn"/> applies a preset to options
/// that already carry another preset (the program applies every <c>-x</c> in turn and a preset sets only the options it lists), which <see cref="Minimap2Program"/> uses for <c>-x splice -x asm5</c>.
/// </remarks>
public static class MinimapPresets
{
    /// <summary>The preset names accepted by <see cref="Apply"/>.</summary>
    public static readonly IReadOnlyList<string> Names =
        ["map-ont", "map-pb", "map10k", "map-hifi", "map-ccs", "ava-ont", "ava-pb", "asm5", "asm10", "asm20", "sr", "short", "splice", "splice:hq", "cdna"];

    /// <summary>The options <c>minimap2 -x <paramref name="preset"/></c> sets (null / empty = no preset: the program's defaults).</summary>
    /// <exception cref="ArgumentException">Unknown preset (the program: <c>unknown preset</c>).</exception>
    public static MinimapMapOptions Apply(string? preset) => ApplyOn(new MinimapMapOptions(), preset);

    /// <summary>
    /// <c>mm_set_opt(preset)</c> applied to the options as they are: like the program, a preset sets only the options it lists, so <c>-x splice -x asm5</c> keeps the splice flags of the first preset (and fails the option check <c>--rmq doesn't work with --splice</c>).
    /// </summary>
    /// <exception cref="ArgumentException">Unknown preset.</exception>
    public static MinimapMapOptions ApplyOn(MinimapMapOptions current, string? preset)
    {
        ArgumentNullException.ThrowIfNull(current);
        var o = current;
        if (string.IsNullOrEmpty(preset))
            return o;
        switch (preset)
        {
            case "map-ont":
                return o;
            case "ava-ont":
                return AllVsAll(o with { Seed = o.Seed with { KmerSize = 15, WindowSize = 5, HomopolymerCompressed = false, OccDist = 0, SkipSelf = true, SkipDual = true }, Bandwidth = 2000, BandwidthLong = 2000 });
            case "map10k" or "map-pb":
                return o with { Seed = o.Seed with { HomopolymerCompressed = true, KmerSize = 19 } };
            case "ava-pb":
                return AllVsAll(o with { Seed = o.Seed with { HomopolymerCompressed = true, KmerSize = 19, WindowSize = 5, OccDist = 0, SkipSelf = true, SkipDual = true }, BandwidthLong = o.Bandwidth });
            case "map-hifi" or "map-ccs":
                return o with
                {
                    Seed = o.Seed with { HomopolymerCompressed = false, KmerSize = 19, WindowSize = 19, OccDist = 500, MinMidOcc = 50, MaxMidOcc = 500 },
                    MaxGap = 10000,
                    Scoring = o.Scoring with { Match = 1, Mismatch = 4, GapOpen = 6, GapOpen2 = 26, GapExtend = 2, GapExtend2 = 1, MinDpMax = 200 },
                };
            case "asm5" or "asm10" or "asm20":
                {
                    var asm = o with
                    {
                        Seed = o.Seed with { HomopolymerCompressed = false, KmerSize = 19, WindowSize = 19, MinMidOcc = 50, MaxMidOcc = 500 },
                        Bandwidth = 1000, BandwidthLong = 100000, MaxGap = 10000, UseRmq = true, BestN = 50,
                    };
                    return preset switch
                    {
                        "asm5" => asm with { Scoring = asm.Scoring with { Match = 1, Mismatch = 19, GapOpen = 39, GapOpen2 = 81, GapExtend = 3, GapExtend2 = 1, ZDrop = 200, ZDropInverted = 200, MinDpMax = 200 } },
                        "asm10" => asm with { Scoring = asm.Scoring with { Match = 1, Mismatch = 9, GapOpen = 16, GapOpen2 = 41, GapExtend = 2, GapExtend2 = 1, ZDrop = 200, ZDropInverted = 200, MinDpMax = 200 } },
                        _ => asm with
                        {
                            Seed = asm.Seed with { WindowSize = 10 },
                            Scoring = asm.Scoring with { Match = 1, Mismatch = 4, GapOpen = 6, GapOpen2 = 26, GapExtend = 2, GapExtend2 = 1, ZDrop = 200, ZDropInverted = 200, MinDpMax = 200 },
                        },
                    };
                }
            case "short" or "sr":
                return o with
                {
                    Seed = o.Seed with { HomopolymerCompressed = false, KmerSize = 21, WindowSize = 11, MidOcc = 1000, MaxOcc = 5000, HeapSort = true },
                    ShortRead = true, FragmentMode = true, PairOrientation = 1, OutputSecondary = false, MaxFragmentLength = 800, MaxGap = 100, Bandwidth = 100, BandwidthLong = 100, PrimaryRatio = 0.5f, MinCount = 2, MinChainScore = 25, BestN = 20,
                    Scoring = o.Scoring with { Match = 2, Mismatch = 8, GapOpen = 12, GapExtend = 2, GapOpen2 = 24, GapExtend2 = 1, ZDrop = 100, ZDropInverted = 100, EndBonus = 10, MinDpMax = 40 },
                };
            case "splice" or "cdna" or "splice:hq":
                {
                    var splice = o with
                    {
                        Seed = o.Seed with { HomopolymerCompressed = false, KmerSize = 15, WindowSize = 5 },
                        Splice = true, SpliceForward = true, SpliceReverse = true, SpliceFlank = true, MaxGap = 2000, MaxGapRef = 200000, Bandwidth = 200000, BandwidthLong = 200000,
                        Scoring = o.Scoring with { Match = 1, Mismatch = 2, GapOpen = 2, GapExtend = 1, GapOpen2 = 32, GapExtend2 = 0, NonCanonical = 9, JunctionBonus = 9, ZDrop = 200, ZDropInverted = 100, MaxSmithWatermanMatrix = 0 },
                    };
                    return preset == "splice:hq" ? splice with { Scoring = splice.Scoring with { JunctionBonus = 5, Mismatch = 4, GapOpen = 6, GapOpen2 = 24 } } : splice;
                }
            default:
                throw new ArgumentException($"Unknown minimap2 preset '{preset}'. Supported: {string.Join(", ", Names)}.", nameof(preset));
        }
    }

    /// <summary>
    /// <c>minimap2 -x <paramref name="preset"/> [options]</c>: the preset (<see cref="Apply(string?)"/>) and then every option that is SET in <paramref name="overrides"/> (the program applies <c>-x</c> before every other option), so an option
    /// whose value equals the default still overrides the preset (<c>-x asm5 -k 15</c> is k = 15, preset k = 19). <c>-G</c> (<see cref="MinimapOverrides.MaxIntronLength"/>) is applied last, as <see cref="WithMaxIntronLength"/>.
    /// </summary>
    /// <exception cref="ArgumentException">Unknown preset.</exception>
    public static MinimapMapOptions Apply(string? preset, MinimapOverrides overrides)
    {
        ArgumentNullException.ThrowIfNull(overrides);
        var o = Apply(preset);
        var seed = o.Seed with
        {
            KmerSize = overrides.KmerSize ?? o.Seed.KmerSize,
            WindowSize = overrides.WindowSize ?? o.Seed.WindowSize,
            HomopolymerCompressed = overrides.HomopolymerCompressed || o.Seed.HomopolymerCompressed,                 // -H has no negation in minimap2 2.26 (ipt.flag |= MM_I_HPC)
            SdustThreshold = overrides.SdustThreshold ?? o.Seed.SdustThreshold,
            MinMidOcc = overrides.MinMidOcc ?? o.Seed.MinMidOcc,
            MaxMidOcc = overrides.MaxMidOcc ?? o.Seed.MaxMidOcc,
            QueryOccurrenceFraction = overrides.QueryOccurrenceFraction is { } q ? (float)q : o.Seed.QueryOccurrenceFraction,
            OccDist = overrides.OccDist ?? o.Seed.OccDist,
            ForwardOnly = overrides.ForwardOnly || o.Seed.ForwardOnly,
            ReverseOnly = overrides.ReverseOnly || o.Seed.ReverseOnly,
            SkipSelf = overrides.AllVsAll || o.Seed.SkipSelf,
            SkipDual = overrides.AllVsAll || o.Seed.SkipDual,
        };
        if (overrides.MidOccFraction is { } f)                                  // main.c -f: x < 1 is the fraction (mid_occ = 0), otherwise the limit itself
            seed = f < 1.0 ? seed with { MidOccFraction = (float)f, MidOcc = 0 } : seed with { MidOcc = (int)(f + .499) };
        o = o with
        {
            Seed = seed,
            MaxGap = overrides.MaxGap ?? o.MaxGap,
            Bandwidth = overrides.Bandwidth ?? o.Bandwidth,
            BandwidthLong = overrides.BandwidthLong ?? o.BandwidthLong,
            MaxChainSkip = overrides.MaxChainSkip ?? o.MaxChainSkip,
            MaxChainIterations = overrides.MaxChainIterations ?? o.MaxChainIterations,
            MinCount = overrides.MinCount ?? o.MinCount,
            MinChainScore = overrides.MinChainScore ?? o.MinChainScore,
            ChainGapScale = overrides.ChainGapScale is { } g ? (float)g : o.ChainGapScale,
            ChainSkipScale = overrides.ChainSkipScale is { } ss ? (float)ss : o.ChainSkipScale,
            NoLongJoin = overrides.NoLongJoin || overrides.AllVsAll || o.NoLongJoin,                                // --no-long-join / -X only switch ON
            AllChains = overrides.AllVsAll || o.AllChains,
            UseRmq = overrides.UseRmq ?? o.UseRmq,                                                                   // --rmq=yes|no
            MaskLevel = overrides.MaskLevel is { } ml ? (float)ml : o.MaskLevel,
            MaskLength = overrides.MaskLength ?? o.MaskLength,
            HardMaskLevel = overrides.HardMaskLevel || o.HardMaskLevel,
            PrimaryRatio = overrides.PrimaryRatio is { } pr ? (float)pr : o.PrimaryRatio,
            BestN = overrides.BestN ?? o.BestN,
            OutputSecondary = overrides.OutputSecondary ?? o.OutputSecondary,                                        // --secondary=yes|no
            MaxFragmentLength = overrides.MaxFragmentLength ?? o.MaxFragmentLength,
        };
        return overrides.MaxIntronLength is { } intron ? WithMaxIntronLength(o, intron) : o;
    }

    /// <summary>
    /// <c>-G NUM</c> / <c>--max-intron-len</c> (<c>mm_mapopt_max_intron_len</c>): only in splice mode, the longest gap on the reference, the bandwidth and the long bandwidth become NUM; without splice mode the program ignores the option.
    /// </summary>
    public static MinimapMapOptions WithMaxIntronLength(MinimapMapOptions options, int maxIntronLength)
    {
        ArgumentNullException.ThrowIfNull(options);
        return options.Splice && maxIntronLength > 0 ? options with { MaxGapRef = maxIntronLength, Bandwidth = maxIntronLength, BandwidthLong = maxIntronLength } : options;
    }

    // MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN with the chaining members the two ava presets share
    private static MinimapMapOptions AllVsAll(MinimapMapOptions o) => o with { AllChains = true, NoLongJoin = true, MinChainScore = 100, PrimaryRatio = 0f, MaxChainSkip = 25 };
}

/// <summary>
/// The options of a <c>minimap2</c> command line that follow <c>-x PRESET</c>, up to the end of the chaining stage: a member that is <c>null</c> was not given on the command line (the preset's value, or the default, stays); any given value overrides
/// the preset, also one that equals the default (<see cref="MinimapPresets.Apply(string?, MinimapOverrides)"/>, audit A6-7). Booleans: <c>-H</c>, <c>-X</c>, <c>--no-long-join</c>, <c>--hard-mask-level</c>, <c>--for-only</c> and <c>--rev-only</c>
/// only switch a flag ON in minimap2 2.26 (<c>false</c> = not given); <c>--rmq=yes|no</c> and <c>--secondary=yes|no</c> take a value (<c>null</c> = not given).
/// </summary>
public sealed record MinimapOverrides
{
    /// <summary><c>-k</c>.</summary>
    public int? KmerSize { get; init; }

    /// <summary><c>-w</c>.</summary>
    public int? WindowSize { get; init; }

    /// <summary><c>-H</c>.</summary>
    public bool HomopolymerCompressed { get; init; }

    /// <summary><c>-T</c>.</summary>
    public int? SdustThreshold { get; init; }

    /// <summary><c>-f</c> first value: below 1 the fraction of repeated minimizers, otherwise the occurrence limit itself.</summary>
    public double? MidOccFraction { get; init; }

    /// <summary><c>-U</c> first value.</summary>
    public int? MinMidOcc { get; init; }

    /// <summary><c>-U</c> second value.</summary>
    public int? MaxMidOcc { get; init; }

    /// <summary><c>--q-occ-frac</c>.</summary>
    public double? QueryOccurrenceFraction { get; init; }

    /// <summary><c>-e</c>.</summary>
    public int? OccDist { get; init; }

    /// <summary><c>--for-only</c>.</summary>
    public bool ForwardOnly { get; init; }

    /// <summary><c>--rev-only</c>.</summary>
    public bool ReverseOnly { get; init; }

    /// <summary><c>-X</c>: all-vs-all (skip self and dual, all chains, no long join).</summary>
    public bool AllVsAll { get; init; }

    /// <summary><c>-g</c>.</summary>
    public int? MaxGap { get; init; }

    /// <summary><c>-r</c> first value.</summary>
    public int? Bandwidth { get; init; }

    /// <summary><c>-r</c> second value.</summary>
    public int? BandwidthLong { get; init; }

    /// <summary><c>--max-chain-skip</c>.</summary>
    public int? MaxChainSkip { get; init; }

    /// <summary><c>--max-chain-iter</c>.</summary>
    public int? MaxChainIterations { get; init; }

    /// <summary><c>-n</c>.</summary>
    public int? MinCount { get; init; }

    /// <summary><c>-m</c>.</summary>
    public int? MinChainScore { get; init; }

    /// <summary><c>--chain-gap-scale</c>.</summary>
    public double? ChainGapScale { get; init; }

    /// <summary><c>--chain-skip-scale</c>.</summary>
    public double? ChainSkipScale { get; init; }

    /// <summary><c>--no-long-join</c>.</summary>
    public bool NoLongJoin { get; init; }

    /// <summary><c>--rmq=yes</c> / <c>--rmq=no</c>.</summary>
    public bool? UseRmq { get; init; }

    /// <summary><c>--mask-level</c> (<c>-M</c>).</summary>
    public double? MaskLevel { get; init; }

    /// <summary><c>--mask-len</c>.</summary>
    public int? MaskLength { get; init; }

    /// <summary><c>--hard-mask-level</c>.</summary>
    public bool HardMaskLevel { get; init; }

    /// <summary><c>-p</c>.</summary>
    public double? PrimaryRatio { get; init; }

    /// <summary><c>-N</c>.</summary>
    public int? BestN { get; init; }

    /// <summary><c>--secondary=yes</c> / <c>--secondary=no</c>.</summary>
    public bool? OutputSecondary { get; init; }

    /// <summary><c>-F</c>.</summary>
    public int? MaxFragmentLength { get; init; }

    /// <summary><c>-G</c> / <c>--max-intron-len</c> (applied last; only in splice mode, see <see cref="MinimapPresets.WithMaxIntronLength"/>).</summary>
    public int? MaxIntronLength { get; init; }
}
