using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>
/// The points of the reader protocol at which a test can stop a thread or count (work package 1b.9g). The production implementation, <see cref="ReaderGateSeams"/>, is an empty struct: the gate is generic over it so that the calls
/// are specialised and inlined away and the path of a call carries no test code.
/// </summary>
internal interface IReaderGateSeams
{
    /// <summary>The process-wide barrier of the releaser (production: <see cref="Interlocked.MemoryBarrierProcessWide"/>).</summary>
    void ProcessWideBarrier();

    /// <summary>A reader has stored its entry and has not yet read the flag.</summary>
    void AfterPublish();

    /// <summary>A reader has cleared its entry and has not yet read the flag.</summary>
    void AfterRetract();

    /// <summary>The barrier has returned; no cell has been read yet.</summary>
    void BeforeScan();

    /// <summary>The cells have been read (<paramref name="inside"/>: a reader was found); the release has not been claimed yet.</summary>
    void AfterScan(bool inside);
}

/// <summary>The production seams: the real barrier, and nothing at the other points.</summary>
internal readonly struct ReaderGateSeams : IReaderGateSeams
{
    public void ProcessWideBarrier() => Interlocked.MemoryBarrierProcessWide();

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void AfterPublish()
    {
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void AfterRetract()
    {
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void BeforeScan()
    {
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void AfterScan(bool inside)
    {
    }
}

/// <summary>
/// The lifetime protocol of one file-backed tree (work package 1b.9g, replacing the striped reader count of 1b.8): readers publish themselves in their thread-exclusive <see cref="ReaderCell"/> with plain stores, <see cref="Dispose"/> raises a flag and a
/// process-wide barrier makes every reader's earlier stores visible before the releaser reads the cells. A call that was in flight when the flag was raised completes and the last of them (or <see cref="Dispose"/>, when none is inside) runs
/// the release once; a call that starts later throws <see cref="ObjectDisposedException"/> without touching what is released. ADR-004 section 5a, "1b.9g", has the proof (x64 and arm64) and the contract.
/// </summary>
/// <typeparam name="TSeams">The seams of the protocol (<see cref="ReaderGateSeams"/> in production).</typeparam>
internal sealed class ReaderGate<TSeams> : IDisposable
    where TSeams : struct, IReaderGateSeams
{
    private readonly object _owner;
    private readonly Action _release;
    private readonly long _id;
    private readonly TSeams _seams;
    private int _disposed;
    private int _released;

    /// <param name="owner">The object that the <see cref="ObjectDisposedException"/> names.</param>
    /// <param name="release">What the last reader out (or <see cref="Dispose"/>) runs once: the unmapping and the unpinning.</param>
    /// <param name="seams">The seams (tests); the default is the production behaviour.</param>
    internal ReaderGate(object owner, Action release, TSeams seams = default)
    {
        ArgumentNullException.ThrowIfNull(owner);
        ArgumentNullException.ThrowIfNull(release);
        _owner = owner;
        _release = release;
        _seams = seams;
        _id = ReaderCell.NextId();
    }

    /// <summary>The id that this gate's readers carry in their cells.</summary>
    internal long Id => _id;

    /// <summary>Whether <see cref="Dispose"/> has run.</summary>
    internal bool IsDisposed => Volatile.Read(ref _disposed) != 0;

    /// <summary>Whether the release has been claimed (and runs or has run).</summary>
    internal bool IsReleased => Volatile.Read(ref _released) != 0;

    /// <summary>
    /// A call begins: the flag is read (a call that starts after <see cref="Dispose"/> throws without publishing anything), the entry is stored in the cell of the thread, and the flag is read again; a flag that is raised by then is
    /// the other half of the handshake with the releaser (it may not have seen the entry), so the entry is taken back and the call goes away. No fence: the barrier of the releaser is what orders the store before the second read.
    /// </summary>
    /// <returns>The cell of the calling thread, to be handed to <see cref="Leave"/>.</returns>
    /// <exception cref="ObjectDisposedException">The gate is disposed.</exception>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal ReaderCell Enter()
    {
        ReaderCell cell = ReaderCell.Current;
        if (Volatile.Read(ref _disposed) != 0) ThrowDisposed();
        cell.Publish(_id);
        _seams.AfterPublish();
        if (IsDisposedAfterStore()) EnterDisposed(cell);
        return cell;
    }

    /// <summary>The call ends: the entry is cleared (a release store) and the flag is read; a raised flag makes this thread try the release, because a releaser that ran before the entry was gone has left it to the last reader.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void Leave(ReaderCell cell)
    {
        cell.Retract();
        _seams.AfterRetract();
        if (IsDisposedAfterStore()) ReleaseWhenIdle();
    }

    /// <summary>Raises the flag and runs the release if no call is inside. Idempotent; does not wait for the calls that are inside.</summary>
    public void Dispose()
    {
        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0) return;
        ReleaseWhenIdle();
    }

    // The read of the flag that follows the store of the entry (or the clearing of it). A separate method that is not inlined, so that no compiler can move the store behind the read: a store may not move across a call of a method the compiler cannot see through.
    [MethodImpl(MethodImplOptions.NoInlining)]
    private bool IsDisposedAfterStore() => Volatile.Read(ref _disposed) != 0;

    [MethodImpl(MethodImplOptions.NoInlining)]
    private void EnterDisposed(ReaderCell cell)
    {
        cell.Retract();
        _seams.AfterRetract();
        ReleaseWhenIdle();
        ThrowDisposed();
    }

    [DoesNotReturn]
    private void ThrowDisposed() => throw new ObjectDisposedException(_owner.GetType().FullName);

    /// <summary>
    /// Runs the release once when the flag is up and no call is inside. Every caller has read the flag as raised. The process-wide barrier comes after that reading and before the scan of the cells: a reader that is not seen by the scan
    /// has its entry after the barrier in its own order and reads the flag as raised (ADR-004, lemma 1). The barrier is not run once the release is claimed, so a stream of late callers costs nothing.
    /// </summary>
    private void ReleaseWhenIdle()
    {
        if (Volatile.Read(ref _released) != 0) return;
        _seams.ProcessWideBarrier();
        Interlocked.MemoryBarrier();
        _seams.BeforeScan();
        bool inside = ReaderCell.AnyInside(_id);
        _seams.AfterScan(inside);
        if (inside) return;
        if (Interlocked.CompareExchange(ref _released, 1, 0) != 0) return;
        _release();
    }
}
