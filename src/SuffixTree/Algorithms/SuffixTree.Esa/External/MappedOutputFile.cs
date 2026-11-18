#pragma warning disable // STUB (tests first): removed by the implementation commit
namespace SuffixTree.Esa.External;

/// <summary>STUB (tests first): how the blocks of a new output file are reserved.</summary>
internal enum ReservationMode
{
    Sparse,
    Preallocate,
    Auto,
}

/// <summary>STUB (tests first): what the creation of an output file did about the disk space.</summary>
internal enum ReservationOutcome
{
    Sparse,
    Preallocated,
    FreeSpaceChecked,
    Unchecked,
}

/// <summary>STUB (tests first): the platform side of the reservation.</summary>
internal interface IReservationSource
{
    bool CanPreallocate { get; }

    FileStream CreatePreallocated(string path, long length);

    long AvailableFreeBytes(string directory);
}

/// <summary>STUB (tests first): the real source.</summary>
internal sealed class SystemReservationSource : IReservationSource
{
    internal static SystemReservationSource Instance { get; } = new();

    public bool CanPreallocate => throw new NotImplementedException();

    public FileStream CreatePreallocated(string path, long length) => throw new NotImplementedException();

    public long AvailableFreeBytes(string directory) => throw new NotImplementedException();
}

/// <summary>STUB (tests first): options of a new output file.</summary>
internal sealed class MappedFileOptions
{
    internal ReservationMode Reservation { get; init; } = ReservationMode.Auto;

    internal IReservationSource? Source { get; init; }
}

/// <summary>STUB (tests first): a file of known length that several writers fill through windows of one shared mapping.</summary>
internal sealed class MappedOutputFile : IDisposable
{
    internal const ReservationMode AutoMode = ReservationMode.Preallocate;

    private MappedOutputFile()
    {
    }

    internal string Path => throw new NotImplementedException();

    internal long Length => throw new NotImplementedException();

    internal ReservationOutcome Outcome => throw new NotImplementedException();

    internal static MappedOutputFile Create(string path, long length, MappedFileOptions? options = null) => throw new NotImplementedException();

    internal MappedWindowWriter CreateWriter(long start, long length, int windowBytes = MappedWindowWriter.DefaultWindowBytes, IoCounters? counters = null) => throw new NotImplementedException();

    internal MappedWindowWriter CreateWriter(int windowBytes = MappedWindowWriter.DefaultWindowBytes, IoCounters? counters = null) => throw new NotImplementedException();

    internal void Complete(bool fsync = true) => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
