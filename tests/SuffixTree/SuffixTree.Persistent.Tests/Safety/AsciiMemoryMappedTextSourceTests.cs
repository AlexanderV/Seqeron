using System;
using System.Collections;
using System.IO;
using System.IO.MemoryMappedFiles;
using System.Linq;
using System.Text;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests.Safety;

/// <summary>
/// Tests for <see cref="AsciiMemoryMappedTextSource"/> — an unsafe memory-mapped
/// ITextSource that reads 1-byte ASCII and widens to char. Covers constructor
/// validation, indexing, Substring, Slice, enumeration, and dispose semantics.
/// </summary>
[TestFixture]
[Category("Safety")]
public class AsciiMemoryMappedTextSourceTests
{
    private string _tempFile = string.Empty;
    private MemoryMappedFile? _mmf;
    private MemoryMappedViewAccessor? _accessor;

    [SetUp]
    public void SetUp()
    {
        _tempFile = Path.GetTempFileName();
    }

    [TearDown]
    public void TearDown()
    {
        _accessor?.Dispose();
        _mmf?.Dispose();

        if (File.Exists(_tempFile))
            File.Delete(_tempFile);
    }

    private AsciiMemoryMappedTextSource CreateSource(string text, long offset = 0)
    {
        byte[] ascii = Encoding.ASCII.GetBytes(text);
        File.WriteAllBytes(_tempFile, ascii);

        _mmf = MemoryMappedFile.CreateFromFile(_tempFile, FileMode.Open, null, ascii.Length,
            MemoryMappedFileAccess.Read);
        _accessor = _mmf.CreateViewAccessor(0, ascii.Length, MemoryMappedFileAccess.Read);

        return new AsciiMemoryMappedTextSource(_accessor, offset, text.Length - (int)offset);
    }

    private AsciiMemoryMappedTextSource CreateSourceWithOffset(string text, long offset, int length)
    {
        byte[] ascii = Encoding.ASCII.GetBytes(text);
        File.WriteAllBytes(_tempFile, ascii);

        _mmf = MemoryMappedFile.CreateFromFile(_tempFile, FileMode.Open, null, ascii.Length,
            MemoryMappedFileAccess.Read);
        _accessor = _mmf.CreateViewAccessor(0, ascii.Length, MemoryMappedFileAccess.Read);

        return new AsciiMemoryMappedTextSource(_accessor, offset, length);
    }

    #region Constructor Validation

    [Test]
    public void Constructor_NullAccessor_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new AsciiMemoryMappedTextSource(null!, 0, 10));
    }

    [Test]
    public void Constructor_NegativeOffset_ThrowsArgumentOutOfRangeException()
    {
        byte[] ascii = Encoding.ASCII.GetBytes("hello");
        File.WriteAllBytes(_tempFile, ascii);

        using var mmf = MemoryMappedFile.CreateFromFile(_tempFile, FileMode.Open, null, ascii.Length,
            MemoryMappedFileAccess.Read);
        using var accessor = mmf.CreateViewAccessor(0, ascii.Length, MemoryMappedFileAccess.Read);

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new AsciiMemoryMappedTextSource(accessor, -1, 5));
    }

    [Test]
    public void Constructor_NegativeLength_ThrowsArgumentOutOfRangeException()
    {
        byte[] ascii = Encoding.ASCII.GetBytes("hello");
        File.WriteAllBytes(_tempFile, ascii);

        using var mmf = MemoryMappedFile.CreateFromFile(_tempFile, FileMode.Open, null, ascii.Length,
            MemoryMappedFileAccess.Read);
        using var accessor = mmf.CreateViewAccessor(0, ascii.Length, MemoryMappedFileAccess.Read);

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new AsciiMemoryMappedTextSource(accessor, 0, -1));
    }

    [Test]
    public void Constructor_OffsetPlusLengthExceedsCapacity_Throws()
    {
        byte[] ascii = Encoding.ASCII.GetBytes("hello");
        File.WriteAllBytes(_tempFile, ascii);

        using var mmf = MemoryMappedFile.CreateFromFile(_tempFile, FileMode.Open, null, ascii.Length,
            MemoryMappedFileAccess.Read);
        using var accessor = mmf.CreateViewAccessor(0, ascii.Length, MemoryMappedFileAccess.Read);

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new AsciiMemoryMappedTextSource(accessor, 3, 10));
    }

    #endregion

    #region Length and Indexer

    [Test]
    public void Length_ReturnsCorrectValue()
    {
        using var source = CreateSource("GATTACA");
        Assert.That(source.Length, Is.EqualTo(7));
    }

    [Test]
    public void Indexer_ReturnsCorrectCharacters()
    {
        using var source = CreateSource("ACGT");

        Assert.Multiple(() =>
        {
            Assert.That(source[0], Is.EqualTo('A'));
            Assert.That(source[1], Is.EqualTo('C'));
            Assert.That(source[2], Is.EqualTo('G'));
            Assert.That(source[3], Is.EqualTo('T'));
        });
    }

    [Test]
    public void Indexer_OutOfRange_Throws()
    {
        using var source = CreateSource("ABC");
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = source[-1]);
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = source[3]);
    }

    #endregion

    #region Substring

    [Test]
    public void Substring_ReturnsCorrectContent()
    {
        using var source = CreateSource("GATTACA");

        Assert.Multiple(() =>
        {
            Assert.That(source.Substring(0, 4), Is.EqualTo("GATT"));
            Assert.That(source.Substring(3, 3), Is.EqualTo("TAC"));
            Assert.That(source.Substring(0, 7), Is.EqualTo("GATTACA"));
        });
    }

    [Test]
    public void Substring_InvalidRange_Throws()
    {
        using var source = CreateSource("ABC");
        Assert.Throws<ArgumentOutOfRangeException>(() => source.Substring(2, 5));
    }

    #endregion

    #region Slice

    [Test]
    public void Slice_ReturnsCorrectContent()
    {
        using var source = CreateSource("GATTACA");

        Assert.Multiple(() =>
        {
            Assert.That(source.Slice(0, 4).ToString(), Is.EqualTo("GATT"));
            Assert.That(source.Slice(4, 3).ToString(), Is.EqualTo("ACA"));
        });
    }

    [Test]
    public void Slice_InvalidRange_Throws()
    {
        using var source = CreateSource("ABC");
        Assert.Throws<ArgumentOutOfRangeException>(() => source.Slice(2, 5));
    }

    #endregion

    #region Enumeration

    [Test]
    public void GetEnumerator_Generic_IteratesAllChars()
    {
        using var source = CreateSource("ACGT");
        var chars = source.ToArray();
        Assert.That(new string(chars), Is.EqualTo("ACGT"));
    }

    [Test]
    public void GetEnumerator_NonGeneric_IteratesAllChars()
    {
        using var source = CreateSource("XYZ");
        var chars = new System.Collections.Generic.List<char>();
        IEnumerable enumerable = source;
        foreach (char c in enumerable)
            chars.Add(c);
        Assert.That(new string(chars.ToArray()), Is.EqualTo("XYZ"));
    }

    #endregion

    #region ToString

    [Test]
    public void ToString_ReturnsFullContent()
    {
        using var source = CreateSource("GATTACA");
        Assert.That(source.ToString(), Is.EqualTo("GATTACA"));
    }

    #endregion

    #region Offset Support

    [Test]
    public void Constructor_WithOffset_ReadsFromCorrectPosition()
    {
        // Write "HELLO" (5 bytes), create source starting at offset 2, length 3 → "LLO"
        using var source = CreateSourceWithOffset("HELLO", offset: 2, length: 3);

        Assert.Multiple(() =>
        {
            Assert.That(source.Length, Is.EqualTo(3));
            Assert.That(source[0], Is.EqualTo('L'));
            Assert.That(source[1], Is.EqualTo('L'));
            Assert.That(source[2], Is.EqualTo('O'));
            Assert.That(source.ToString(), Is.EqualTo("LLO"));
        });
    }

    #endregion

    #region Dispose Semantics

    [Test]
    public void Dispose_ThenAccess_ThrowsObjectDisposedException()
    {
        var source = CreateSource("ABC");
        source.Dispose();

        Assert.Multiple(() =>
        {
            Assert.Throws<ObjectDisposedException>(() => _ = source[0]);
            Assert.Throws<ObjectDisposedException>(() => source.Substring(0, 1));
            Assert.Throws<ObjectDisposedException>(() => source.Slice(0, 1));
            Assert.Throws<ObjectDisposedException>(() => source.ToString());
            Assert.Throws<ObjectDisposedException>(() => source.GetEnumerator());
        });
    }

    [Test]
    public void Dispose_MultipleTimes_DoesNotThrow()
    {
        var source = CreateSource("ABC");
        source.Dispose();
        Assert.DoesNotThrow(() => source.Dispose());
    }

    #endregion
}

