using System;
using System.IO;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Tests for Q4 (accessor leak if AcquirePointer throws) and
/// Q5 (ToString after Dispose → AccessViolationException instead of ODE).
/// </summary>
[TestFixture]
public class TextSourceSafetyTests
{
    private string _tempFile = string.Empty;

    [SetUp]
    public void SetUp() => _tempFile = Path.GetTempFileName();

    [TearDown]
    public void TearDown()
    {
        try { if (File.Exists(_tempFile)) File.Delete(_tempFile); } catch { }
    }

    // ─── Q5: ToString after Dispose must throw ODE ──────────────────

    [Test]
    public void ToString_AfterDispose_ThrowsObjectDisposedException()
    {
        string text = "ABCD";
        File.WriteAllBytes(_tempFile, System.Text.Encoding.Unicode.GetBytes(text));

        var source = new MemoryMappedTextSource(_tempFile, 0, text.Length);
        // Verify it works before dispose
        Assert.That(source.ToString(), Is.EqualTo("ABCD"));

        source.Dispose();

        // After dispose, must throw ODE — not AccessViolationException
        Assert.Throws<ObjectDisposedException>(() => _ = source.ToString());
    }

    // ─── Q4: Constructor leak verification (structural test) ────────
    // We can't easily force AcquirePointer to fail, but we verify that
    // a successful construction properly disposes when Dispose is called
    // (this tests the happy-path resource lifecycle).

    [Test]
    public void Constructor_FromFile_DisposesCleanly()
    {
        string text = "Hello";
        File.WriteAllBytes(_tempFile, System.Text.Encoding.Unicode.GetBytes(text));

        var source = new MemoryMappedTextSource(_tempFile, 0, text.Length);
        Assert.That(source.Length, Is.EqualTo(5));
        source.Dispose();

        // Double dispose should not throw
        Assert.DoesNotThrow(() => source.Dispose());
    }
}
