namespace SuffixTree.Esa;

/// <summary>One level of the SA-IS recursion: the text length n, its alphabet size, and the number m of LMS positions (the length of the next level).</summary>
internal readonly record struct SaisLevel(int Length, int AlphabetSize, int LmsCount);

/// <summary>
/// Optional instrumentation of a suffix-array build: the levels in recursion order (level 0 is the input text). Used by the
/// tests to assert the depth bound of ADR-002 section 2 item 5; the kernels pay nothing when no trace is passed.
/// </summary>
internal sealed class SaisTrace
{
    /// <summary>The levels in recursion order.</summary>
    internal List<SaisLevel> Levels { get; } = new();

    /// <summary>Number of levels, i.e. 1 + the deepest recursion nesting.</summary>
    internal int Depth => Levels.Count;
}
