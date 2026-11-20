namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The <c>-x</c> presets of minimap2 2.26 (<c>options.c</c> <c>mm_set_opt</c>, read at raw.githubusercontent.com/lh3/minimap2 v2.26) as <see cref="MinimapMapOptions"/>: <c>map-ont</c> (the defaults), <c>map-pb</c> / <c>map10k</c>, <c>map-hifi</c> /
/// <c>map-ccs</c>, <c>ava-ont</c>, <c>ava-pb</c>, <c>asm5</c> / <c>asm10</c> / <c>asm20</c>, <c>sr</c> / <c>short</c>, <c>splice</c> / <c>cdna</c> and <c>splice:hq</c>. The program applies <c>-x</c> BEFORE every other option, so an explicit
/// option overrides the preset: <c>MinimapPresets.Apply("asm5") with { MinChainScore = 30 }</c> is <c>minimap2 -x asm5 -m 30</c>.
/// </summary>
/// <remarks>
/// The chaining-stage members (<c>-k -w -H -f/occ -g -r -n -m -p -N -S ... --rmq -P -D --no-dual --no-long-join -F --sr</c>) are read by <see cref="Minimap2Chaining.Map"/>; the base-level members (<c>-A -B -O -E -z -s -u</c>, splice flags) are stored
/// in <see cref="MinimapMapOptions.Scoring"/> / <see cref="MinimapMapOptions.Splice"/> for the base-level stage (audit A3-9c2). <c>-x splice</c> additionally needs the cDNA chaining of <c>mg_lchain_dp</c> (<c>is_cdna</c>) and the
/// multi-segment fragment mode of <c>-x sr</c> with paired reads; both stay with the carried items (QUEUE.md / B09.md A5-9).
/// </remarks>
public static class MinimapPresets
{
    /// <summary>The preset names accepted by <see cref="Apply"/>.</summary>
    public static readonly IReadOnlyList<string> Names =
        ["map-ont", "map-pb", "map10k", "map-hifi", "map-ccs", "ava-ont", "ava-pb", "asm5", "asm10", "asm20", "sr", "short", "splice", "splice:hq", "cdna"];

    /// <summary>The options <c>minimap2 -x <paramref name="preset"/></c> sets (null / empty = no preset: the program's defaults).</summary>
    /// <exception cref="ArgumentException">Unknown preset (the program: <c>unknown preset</c>).</exception>
    public static MinimapMapOptions Apply(string? preset)
    {
        var o = new MinimapMapOptions();
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
                    Seed = o.Seed with { HomopolymerCompressed = false, KmerSize = 21, WindowSize = 11, MidOcc = 1000, MaxOcc = 5000 },
                    ShortRead = true, OutputSecondary = false, MaxFragmentLength = 800, MaxGap = 100, Bandwidth = 100, BandwidthLong = 100, PrimaryRatio = 0.5f, MinCount = 2, MinChainScore = 25, BestN = 20,
                    Scoring = o.Scoring with { Match = 2, Mismatch = 8, GapOpen = 12, GapExtend = 2, GapOpen2 = 24, GapExtend2 = 1, ZDrop = 100, ZDropInverted = 100, EndBonus = 10, MinDpMax = 40 },
                };
            case "splice" or "cdna" or "splice:hq":
                {
                    var splice = o with
                    {
                        Seed = o.Seed with { HomopolymerCompressed = false, KmerSize = 15, WindowSize = 5 },
                        Splice = true, MaxGap = 2000, MaxGapRef = 200000, Bandwidth = 200000, BandwidthLong = 200000,
                        Scoring = o.Scoring with { Match = 1, Mismatch = 2, GapOpen = 2, GapExtend = 1, GapOpen2 = 32, GapExtend2 = 0, NonCanonical = 9, JunctionBonus = 9, ZDrop = 200, ZDropInverted = 100, MaxSmithWatermanMatrix = 0 },
                    };
                    return preset == "splice:hq" ? splice with { Scoring = splice.Scoring with { JunctionBonus = 5, Mismatch = 4, GapOpen = 6, GapOpen2 = 24 } } : splice;
                }
            default:
                throw new ArgumentException($"Unknown minimap2 preset '{preset}'. Supported: {string.Join(", ", Names)}.", nameof(preset));
        }
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
