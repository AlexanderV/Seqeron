#pragma warning disable // STUB (tests first): removed by the implementation commit
namespace SuffixTree.Esa.External;

/// <summary>STUB (tests first): options of the streaming packed-text writer.</summary>
internal sealed class PackedTextWriterOptions
{
    internal int BlockSymbols { get; init; } = 1 << 20;

    internal int WindowBytes { get; init; } = MappedWindowWriter.DefaultWindowBytes;

    internal int MinBits { get; init; } = PackedTextFormat.MinBits;

    internal long MaxSymbols { get; init; } = PackedTextFormat.MaxSymbols;

    internal IoCounters? Counters { get; init; }

    internal MappedFileOptions? File { get; init; }

    internal Action? BetweenPasses { get; init; }
}

/// <summary>STUB (tests first): the streaming PH0-style writer of the packed text.</summary>
internal static class PackedTextWriter
{
    internal const int BlockMultiple = 8;

    internal static long InputLength(string inputPath) => throw new NotImplementedException();

    internal static PackedTextInfo Scan(string inputPath, PackedTextWriterOptions? options = null) => throw new NotImplementedException();

    internal static void Pack(string inputPath, PackedTextInfo info, MappedWindowWriter target, PackedTextWriterOptions? options = null) => throw new NotImplementedException();

    internal static PackedTextInfo WriteFile(string inputPath, string packedPath, PackedTextWriterOptions? options = null) => throw new NotImplementedException();
}
