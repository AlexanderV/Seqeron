#pragma warning disable S2325, S2933, S1144, CA1823, IDE0052, CA1822, S1172, IDE0060, CA1801, S3251, S1186, S3453, CS9113 // STUB: removed when the implementation lands
namespace SuffixTree.Esa.External;

/// <summary>STUB (work package 1b.7d, tests first). Where the leaf buckets write their slices of the suffix array (entry offsets, 4 bytes per entry).</summary>
internal interface ISaSliceSink
{
    /// <summary>A writer for one worker (not thread-safe, disposed by the worker).</summary>
    ISaSliceWriter CreateWriter();
}

/// <summary>STUB (work package 1b.7d, tests first).</summary>
internal interface ISaSliceWriter : IDisposable
{
    void Write(long entryOffset, ReadOnlySpan<uint> data);
}

/// <summary>STUB (work package 1b.7d, tests first). The suffix array file: a <see cref="MappedOutputFile"/> of <c>4 (n + 1)</c> bytes, one sliding-window writer per worker.</summary>
internal sealed class MappedSaSliceSink(MappedOutputFile file, int windowBytes, IoCounters? counters) : ISaSliceSink
{
    public ISaSliceWriter CreateWriter() => throw new NotImplementedException();
}
