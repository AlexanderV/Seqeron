using SuffixTree.Esa;

namespace SuffixTree.Esa.Tests.Tree;

/// <summary>Records the events of <see cref="LcpIntervalScan"/> (test support).</summary>
internal sealed class ScanLog
{
    internal List<(int Previous, int K)> LIndices { get; } = [];

    internal List<(int Ell, int Lb, int Rb, int First, bool LastChild)> Closes { get; } = [];
}

/// <summary>An <see cref="IIntervalSink"/> that appends every event to a <see cref="ScanLog"/>.</summary>
internal readonly struct LoggingSink(ScanLog log) : IIntervalSink
{
    public void LIndex(int previous, int k) => log.LIndices.Add((previous, k));

    public void Close(int ell, int lb, int rb, int first, bool lastChild) => log.Closes.Add((ell, lb, rb, first, lastChild));
}
