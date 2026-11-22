#pragma warning disable S2325 // tests-first stub: the members are instance members of the implementation that follows
#pragma warning disable IDE0060 // the parameters are used by the implementation that follows
namespace SuffixTree.Esa;

/// <summary>The points of the reader protocol at which a test can stop a thread or count (work package 1b.9g).</summary>
internal interface IReaderGateSeams
{
    /// <summary>The process-wide barrier of the releaser.</summary>
    void ProcessWideBarrier();

    /// <summary>A reader has stored its entry and has not yet read the flag.</summary>
    void AfterPublish();

    /// <summary>A reader has cleared its entry and has not yet read the flag.</summary>
    void AfterRetract();

    /// <summary>The barrier has returned; no cell has been read yet.</summary>
    void BeforeScan();

    /// <summary>The cells have been read; the release has not been claimed yet.</summary>
    void AfterScan(bool inside);
}

/// <summary>The production seams. Tests-first stub.</summary>
internal readonly struct ReaderGateSeams : IReaderGateSeams
{
    public void ProcessWideBarrier() => Interlocked.MemoryBarrierProcessWide();

    public void AfterPublish()
    {
    }

    public void AfterRetract()
    {
    }

    public void BeforeScan()
    {
    }

    public void AfterScan(bool inside)
    {
    }
}

/// <summary>The lifetime protocol of one file-backed tree (work package 1b.9g). Tests-first stub: the members exist so that the tests compile and fail.</summary>
/// <typeparam name="TSeams">The seams of the protocol.</typeparam>
internal sealed class ReaderGate<TSeams> : IDisposable
    where TSeams : struct, IReaderGateSeams
{
    internal ReaderGate(object owner, Action release, TSeams seams = default)
    {
    }

    internal long Id => throw new NotImplementedException();

    internal bool IsDisposed => throw new NotImplementedException();

    internal bool IsReleased => throw new NotImplementedException();

    internal ReaderCell Enter() => throw new NotImplementedException();

    internal void Leave(ReaderCell cell) => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
