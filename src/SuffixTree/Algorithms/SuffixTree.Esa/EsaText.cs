using SuffixTree.Esa.External;
using SuffixTree.Indexing;

namespace SuffixTree.Esa;

/// <summary>
/// The text of an <see cref="EsaIndex"/> as a sequence of UTF-16 code units, whichever way it is stored (work package 1b.7f): the caller's <see cref="string"/> (an index built in memory) or the packed TEXT section of a mapped
/// SESA container, whose dense symbols the alphabet maps back to code units. Reads of the packed text are checked against the length (a position that a damaged suffix array supplied is an <see cref="InvalidDataException"/>,
/// never a read past the section). Immutable and safe for concurrent readers.
/// </summary>
internal sealed class EsaText
{
    private const int ChunkSymbols = 4096;

    private readonly string? _string;
    private readonly MappedTextReader? _reader;
    private readonly char[]? _alphabet;
    private readonly int _length;

    private EsaText(string? text, MappedTextReader? reader, char[]? alphabet, int length)
    {
        _string = text;
        _reader = reader;
        _alphabet = alphabet;
        _length = length;
    }

    /// <summary>The text of a string (no copy).</summary>
    internal static EsaText FromString(string text) => new(text, null, null, text.Length);

    /// <summary>The text of a packed section: <paramref name="alphabet"/>[symbol] is the code unit of a dense symbol.</summary>
    internal static EsaText FromPacked(MappedTextReader reader, char[] alphabet, int length) => new(null, reader, alphabet, length);

    /// <summary>Number of code units.</summary>
    internal int Length => _length;

    /// <summary>True when the text is the packed section of a container.</summary>
    internal bool IsPacked => _reader is not null;

    /// <summary>The code unit at <paramref name="index"/>.</summary>
    /// <exception cref="IndexOutOfRangeException">An in-memory text: the index is outside <c>0..n-1</c>.</exception>
    /// <exception cref="InvalidDataException">A packed text: the index is outside <c>0..n-1</c> (a damaged table sent the read there) or the symbol has no letter.</exception>
    internal char this[int index]
    {
        get
        {
            if (_string is not null) return _string[index];
            if ((uint)index >= (uint)_length) ThrowOutside(index);
            return Letter(_reader!.Symbol(index));
        }
    }

#pragma warning disable S2325 // stubs of the tests-first step
    /// <summary>The text of a packed section of any length (work package S2.3).</summary>
    internal static EsaText FromPacked(MappedTextReader reader, char[] alphabet, long length) => throw new NotImplementedException();

    /// <summary>The code unit at a 64-bit position.</summary>
    internal char this[long index] => throw new NotImplementedException();

    /// <summary>The code units <c>[start, start + length)</c> as a new string (the 64-bit start).</summary>
    internal string Substring(long start, int length) => throw new NotImplementedException();

    /// <summary>The code units from the 64-bit position <paramref name="start"/> to the end as a new string.</summary>
    internal string Substring(long start) => throw new NotImplementedException();

    /// <summary>Copies <c>destination.Length</c> code units that start at the 64-bit position <paramref name="start"/>.</summary>
    internal void CopyTo(long start, Span<char> destination) => throw new NotImplementedException();

    /// <summary>The reader of the dense symbols of a packed text.</summary>
    internal MappedTextReader Reader => _reader ?? throw new InvalidOperationException("the text is a string");

    /// <summary>The letters of the alphabet of a packed text.</summary>
    internal char[] Alphabet => _alphabet ?? throw new InvalidOperationException("the text is a string");

    /// <summary>The code units <c>[start, start + length)</c> as a new string.</summary>
    internal string Substring(int start, int length)
    {
        if (_string is not null) return _string.Substring(start, length);
        CheckRange(start, length);
        return string.Create(length, (this, start), static (span, state) => state.Item1.CopyChars(state.start, span));
    }

    /// <summary>The code units from <paramref name="start"/> to the end as a new string.</summary>
    internal string Substring(int start)
    {
        if (_string is not null) return _string.Substring(start);
        if ((uint)start > (uint)_length) ThrowOutside(start);
        return Substring(start, _length - start);
    }

    /// <summary>The code units as memory: a slice of the string (no copy) or a new string's memory.</summary>
    internal ReadOnlyMemory<char> AsMemory(int start, int length) => _string is not null ? _string.AsMemory(start, length) : Substring(start, length).AsMemory();

    /// <summary>Copies <c>destination.Length</c> code units that start at <paramref name="start"/> (the range is checked).</summary>
    internal void CopyTo(int start, Span<char> destination)
    {
        if (_string is not null) _string.AsSpan(start, destination.Length).CopyTo(destination);
        else
        {
            CheckRange(start, destination.Length);
            CopyChars(start, destination);
        }
    }

    /// <summary>The text as <see cref="ITextSource"/> (the facade's <c>Text</c> and the navigator's).</summary>
    internal ITextSource AsTextSource() => _string is not null ? new StringTextSource(_string) : new PackedTextSource(this);

    /// <summary>The text as <see cref="IIndexText"/> (the long contract's).</summary>
    internal IIndexText AsIndexText() => _string is not null ? new EsaIndexText(_string) : new PackedIndexText(this);

    private void CopyChars(int start, Span<char> destination)
    {
        Span<byte> buffer = stackalloc byte[ChunkSymbols];
        for (int done = 0; done < destination.Length; done += ChunkSymbols)
        {
            int take = Math.Min(ChunkSymbols, destination.Length - done);
            _reader!.ReadSymbols(start + done, buffer[..take]);
            for (int i = 0; i < take; i++) destination[done + i] = Letter(buffer[i]);
        }
    }

    private char Letter(uint symbol)
    {
        char[] alphabet = _alphabet!;
        if (symbol >= (uint)alphabet.Length) throw new InvalidDataException($"the text holds the symbol {symbol}, but the alphabet has {alphabet.Length} letters: the file is damaged");
        return alphabet[symbol];
    }

    private void CheckRange(int start, int length)
    {
        if (start < 0 || length < 0 || start > _length - length) ThrowOutside(start);
    }

    private void ThrowOutside(long index) => throw new InvalidDataException($"the text is read at position {index}, but it has {_length} symbols: the file is damaged");

    private sealed class PackedTextSource(EsaText text) : ITextSource
    {
        public int Length => text.Length;

        public char this[int index] => text[index];

        public string Substring(int start, int length) => text.Substring(start, length);

        public ReadOnlySpan<char> Slice(int start, int length) => text.Substring(start, length).AsSpan();

        public IEnumerator<char> GetEnumerator()
        {
            for (int i = 0; i < text.Length; i++) yield return text[i];
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();

        public override string ToString() => text.Substring(0, text.Length);
    }

    private sealed class PackedIndexText(EsaText text) : IIndexText
    {
        public long Length => text.Length;

        public int SymbolAt(long position)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(position);
            return position >= text.Length ? -1 : text[(int)position];
        }

        public int CopyTo(long start, Span<char> destination)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(start);
            ArgumentOutOfRangeException.ThrowIfGreaterThan(start, text.Length);
            int count = (int)Math.Min(destination.Length, text.Length - start);
            text.CopyTo((int)start, destination[..count]);
            return count;
        }
    }
}
