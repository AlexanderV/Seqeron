using SuffixTree.Esa.External;

namespace SuffixTree.Esa;

/// <summary>
/// The pattern descent over the sections of a container (work package 1b.9a): STUB of the red commit; the tests of <c>MappedSearchTests</c> are written against it and fail until the next commit.
/// </summary>
internal sealed unsafe class MappedSearch
{
    private readonly SesaFile _file;

    internal MappedSearch(SesaFile file)
    {
        ArgumentNullException.ThrowIfNull(file);
        _file = file;
    }

    internal bool TryMatch(ReadOnlySpan<byte> pattern, out int lo, out int hi) => throw new NotImplementedException($"1b.9a: the specialised descent is not written yet ({_file.Length} bytes, {pattern.Length} symbols)");
}
