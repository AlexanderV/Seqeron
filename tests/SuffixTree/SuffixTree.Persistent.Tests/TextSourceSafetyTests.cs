using System;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
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

    // ─── C17: Concurrent Dispose must not cause NullReferenceException ──

    [Test]
    public void ConcurrentDispose_NeverThrowsNullReferenceException()
    {
        // The TOCTOU fix snapshots _ptr before the disposed check so that
        // concurrent Dispose (which nulls _ptr) can never cause NRE.
        // Expected: readers see either valid data or ObjectDisposedException.
        const int Iterations = 200;
        for (int i = 0; i < Iterations; i++)
        {
            string text = "ABCDEFGH";
            File.WriteAllBytes(_tempFile, System.Text.Encoding.Unicode.GetBytes(text));

            var source = new MemoryMappedTextSource(_tempFile, 0, text.Length);
            var barrier = new Barrier(2);
            Exception? readerException = null;

            var reader = Task.Run(() =>
            {
                barrier.SignalAndWait();
                try
                {
                    // Rapid reads — some will race with Dispose
                    for (int j = 0; j < 50; j++)
                    {
                        _ = source[0];
                        _ = source.Substring(0, 2);
                        _ = source.ToString();
                    }
                }
                catch (ObjectDisposedException) { /* expected */ }
                catch (Exception ex) { Volatile.Write(ref readerException, ex); }
            });

            var disposer = Task.Run(() =>
            {
                barrier.SignalAndWait();
                source.Dispose();
            });

            Task.WaitAll(reader, disposer);

            var caught = Volatile.Read(ref readerException);
            Assert.That(caught, Is.Null,
                $"C17: Reader must see ODE or valid data, never {caught?.GetType().Name}: {caught?.Message}");
        }
    }
}
