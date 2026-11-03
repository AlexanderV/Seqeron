namespace SuffixTree.Esa.Tests.Kernels;

/// <summary>Collects the events of an <see cref="ISaisProbe"/> into <see cref="SaisLevelState"/> records (the form the oracle produces).</summary>
internal sealed class SaisStateRecorder
{
    internal sealed class Partial
    {
        internal int[] Counts = [];
        internal int[] Heads = [];
        internal int[] Tails = [];
        internal int[] LmsSorted = [];
        internal int[] LmsPositions = [];
        internal int[] Reduced = [];
        internal int Names;
    }

    private readonly SortedDictionary<int, Partial> _levels = [];

    internal void Record(int level, Action<Partial> set)
    {
        if (!_levels.TryGetValue(level, out Partial? p)) _levels[level] = p = new Partial();
        set(p);
    }

    /// <summary>The recorded levels in recursion order.</summary>
    internal List<SaisLevelState> ToStates() =>
        [.. _levels.Values.Select(p => new SaisLevelState(p.Counts.Sum(), p.Counts.Length, p.Counts, p.Heads, p.Tails, p.LmsPositions, p.LmsSorted, p.Names, p.Reduced))];
}

/// <summary>The recording probe: copies every reported span.</summary>
internal readonly struct RecordingSaisProbe(SaisStateRecorder recorder) : ISaisProbe
{
    public static bool Enabled => true;

    public void Counts(int level, ReadOnlySpan<int> counts)
    {
        int[] copy = counts.ToArray();
        recorder.Record(level, p => p.Counts = copy);
    }

    public void Tails(int level, ReadOnlySpan<int> tails)
    {
        int[] copy = tails.ToArray();
        recorder.Record(level, p => p.Tails = copy);
    }

    public void Heads(int level, ReadOnlySpan<int> heads)
    {
        int[] copy = heads.ToArray();
        recorder.Record(level, p => p.Heads = copy);
    }

    public void LmsSorted(int level, ReadOnlySpan<int> positions)
    {
        int[] copy = positions.ToArray();
        recorder.Record(level, p => p.LmsSorted = copy);
    }

    public void LmsPositions(int level, ReadOnlySpan<int> positions)
    {
        int[] copy = positions.ToArray();
        recorder.Record(level, p => p.LmsPositions = copy);
    }

    public void Reduced(int level, int names, ReadOnlySpan<int> reduced)
    {
        int[] copy = reduced.ToArray();
        recorder.Record(level, p =>
        {
            p.Names = names;
            p.Reduced = copy;
        });
    }
}
