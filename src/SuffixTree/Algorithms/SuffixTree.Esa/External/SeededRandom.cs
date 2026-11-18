#pragma warning disable S2325, S2933, S1144, CA1823, IDE0052, CA1822, S1172, IDE0060, CA1801, S3251, S1186 // STUB: removed when the implementation lands
namespace SuffixTree.Esa.External;

/// <summary>STUB (work package 1b.7d, tests first).</summary>
internal sealed class SeededRandom(ulong seed)
{
    private ulong _state = seed;

    internal ulong Next() => throw new NotImplementedException();

    internal int NextInt(int bound) => throw new NotImplementedException();
}

/// <summary>STUB (work package 1b.7d, tests first).</summary>
internal sealed class WorkPool<T>
    where T : class
{
    internal Exception? Error => throw new NotImplementedException();

    internal void Push(T task) => throw new NotImplementedException();

    internal bool Take(out T task) => throw new NotImplementedException();

    internal void Done() => throw new NotImplementedException();

    internal void Fail(Exception e) => throw new NotImplementedException();
}
