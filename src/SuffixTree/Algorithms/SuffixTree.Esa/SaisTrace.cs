using System.Diagnostics;

namespace SuffixTree.Esa;

/// <summary>One level of the SA-IS recursion: the text length n, its alphabet size, and the number m of LMS positions (the length of the next level).</summary>
internal readonly record struct SaisLevel(int Length, int AlphabetSize, int LmsCount);

/// <summary>
/// The exclusive phases of one SA-IS level (the time of the recursive call is not part of any phase of the calling level): symbol counting, type
/// classification, round 1 (LMS suffixes into bucket ends, induce L, induce S: sorts the LMS substrings), moving the sorted LMS positions, naming,
/// gathering the reduced string, the reduced-problem bookkeeping around the recursion, round 2 (sorted LMS suffixes into bucket ends, induce L, induce S).
/// </summary>
internal enum SaisPhase
{
    Count,
    Classify,
    Round1Place,
    Round1InduceL,
    Round1InduceS,
    CompactLms,
    Name,
    CompactNames,
    Reduced,
    Round2Place,
    Round2InduceL,
    Round2InduceS,
}

/// <summary>
/// Optional instrumentation of a suffix-array build: the levels in recursion order (level 0 is the input text) and stopwatch ticks per level and
/// <see cref="SaisPhase"/>. Used by the tests to assert the depth bound of ADR-002 section 2 item 5 and by the baseline harness (work package 1a.9)
/// for the time share of the phases; the kernels pay one null check per phase when no trace is passed.
/// </summary>
internal sealed class SaisTrace
{
    private readonly List<long[]> _ticks = new();

    /// <summary>The levels in recursion order.</summary>
    internal List<SaisLevel> Levels { get; } = new();

    /// <summary>Number of levels, i.e. 1 + the deepest recursion nesting.</summary>
    internal int Depth => Levels.Count;

    /// <summary>Stopwatch ticks of <paramref name="phase"/> at recursion <paramref name="level"/> (0 when the level was not reached).</summary>
    internal long Ticks(int level, SaisPhase phase) => level < _ticks.Count ? _ticks[level][(int)phase] : 0;

    /// <summary>Stopwatch ticks of <paramref name="phase"/> summed over all levels.</summary>
    internal long TotalTicks(SaisPhase phase) => _ticks.Sum(t => t[(int)phase]);

    /// <summary>Start timestamp for a phase; 0 without a trace.</summary>
    internal static long Begin(SaisTrace? trace) => trace is null ? 0 : Stopwatch.GetTimestamp();

    /// <summary>Adds the time since <paramref name="begin"/> to a phase of a level.</summary>
    internal void End(int level, SaisPhase phase, long begin)
    {
        long now = Stopwatch.GetTimestamp();
        while (_ticks.Count <= level) _ticks.Add(new long[Enum.GetValues<SaisPhase>().Length]);
        _ticks[level][(int)phase] += now - begin;
    }
}

/// <summary>Stopwatch ticks of the four steps of <see cref="LcpReference"/>: Phi, PLCP, escape counting, and writing the compact LCP.</summary>
internal sealed class LcpTrace
{
    /// <summary>Phi: <c>Phi[SA[i]] = SA[i-1]</c>.</summary>
    internal long PhiTicks { get; set; }

    /// <summary>PLCP: the character comparisons in text order.</summary>
    internal long PlcpTicks { get; set; }

    /// <summary>Counting the values at or above the escape thresholds and choosing the representation.</summary>
    internal long CountTicks { get; set; }

    /// <summary>Writing <c>LCP[r] = PLCP[SA[r]]</c> into the compact array.</summary>
    internal long OutputTicks { get; set; }
}
