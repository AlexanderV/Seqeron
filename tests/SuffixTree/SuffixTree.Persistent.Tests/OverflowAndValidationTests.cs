using System;
using System.IO;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Tests for P1 (integer overflow in bounds), P7 (zero-capacity infinite loop),
/// P11 (file-ctor parameter validation).
/// Written RED-first: these tests expose existing bugs.
/// </summary>
[TestFixture]
public class OverflowAndValidationTests
{
    private string _tempFile = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _tempFile = Path.GetTempFileName();
    }

    [TearDown]
    public void TearDown()
    {
        try { if (File.Exists(_tempFile)) File.Delete(_tempFile); } catch { }
    }

    // ─── P1: Integer overflow in Substring / Slice bounds ────────────

    [Test]
    public void Substring_OverflowingStartPlusLength_ShouldThrow()
    {
        string text = "ABCD";
        File.WriteAllBytes(_tempFile, System.Text.Encoding.Unicode.GetBytes(text));

        using var source = new MemoryMappedTextSource(_tempFile, 0, text.Length);
        // start=1, length=int.MaxValue → start+length overflows to int.MinValue < _length → BUG: passes check
        Assert.Throws<IndexOutOfRangeException>(() => source.Substring(1, int.MaxValue));
    }

    [Test]
    public void Slice_OverflowingStartPlusLength_ShouldThrow()
    {
        string text = "ABCD";
        File.WriteAllBytes(_tempFile, System.Text.Encoding.Unicode.GetBytes(text));

        using var source = new MemoryMappedTextSource(_tempFile, 0, text.Length);
        Assert.Throws<IndexOutOfRangeException>(() => source.Slice(1, int.MaxValue));
    }

    [Test]
    public void Substring_NegativeLength_ShouldThrow()
    {
        string text = "ABCD";
        File.WriteAllBytes(_tempFile, System.Text.Encoding.Unicode.GetBytes(text));

        using var source = new MemoryMappedTextSource(_tempFile, 0, text.Length);
        Assert.Throws<IndexOutOfRangeException>(() => source.Substring(0, -1));
    }

    [Test]
    public void Slice_NegativeLength_ShouldThrow()
    {
        string text = "ABCD";
        File.WriteAllBytes(_tempFile, System.Text.Encoding.Unicode.GetBytes(text));

        using var source = new MemoryMappedTextSource(_tempFile, 0, text.Length);
        Assert.Throws<IndexOutOfRangeException>(() => source.Slice(0, -1));
    }

    // ─── P7: HeapStorageProvider(0) → EnsureCapacity infinite loop ──

    [Test]
    public void HeapProvider_ZeroCapacity_AllocateDoesNotHang()
    {
        var provider = new HeapStorageProvider(0);
        // Without fix: EnsureCapacity loops forever (0 * 2 = 0)
        long offset = provider.Allocate(1);
        Assert.That(offset, Is.EqualTo(0));
        provider.Dispose();
    }

    [Test]
    public void HeapProvider_ZeroCapacity_WriteAndRead()
    {
        var provider = new HeapStorageProvider(0);
        provider.Allocate(4);
        provider.WriteInt32(0, 42);
        Assert.That(provider.ReadInt32(0), Is.EqualTo(42));
        provider.Dispose();
    }
}
