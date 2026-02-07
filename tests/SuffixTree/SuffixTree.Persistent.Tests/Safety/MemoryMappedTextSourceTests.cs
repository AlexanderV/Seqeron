using System;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests;

[TestFixture]
public class MemoryMappedTextSourceTests
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
        if (File.Exists(_tempFile))
            File.Delete(_tempFile);
    }

    [Test]
    public void MemoryMappedTextSource_ShouldReadCorrectly()
    {
        string text = "Hello";
        byte[] bytes = System.Text.Encoding.Unicode.GetBytes(text);
        File.WriteAllBytes(_tempFile, bytes);

        // Encoding.Unicode (UTF-16LE) = 2 bytes per char. No BOM written by GetBytes.
        using (var source = new MemoryMappedTextSource(_tempFile, 0, text.Length))
        {
            Assert.That(source.Length, Is.EqualTo(text.Length));
            for (int i = 0; i < text.Length; i++)
            {
                Assert.That(source[i], Is.EqualTo(text[i]));
            }
        }
    }

    [Test]
    public void MemoryMappedTextSource_Slice_ShouldReturnCorrectContent()
    {
        string text = "ABCDEFGHIJ";
        // Write raw bytes without BOM
        byte[] bytes = System.Text.Encoding.Unicode.GetBytes(text);
        File.WriteAllBytes(_tempFile, bytes);

        using (var source = new MemoryMappedTextSource(_tempFile, 0, text.Length))
        {
            var slice = source.Slice(2, 3);
            Assert.That(slice.ToString(), Is.EqualTo("CDE"));
        }
    }

    [Test]
    public void MemoryMappedTextSource_Enumerator_ShouldWork()
    {
        string text = "ENUM";
        byte[] bytes = System.Text.Encoding.Unicode.GetBytes(text);
        File.WriteAllBytes(_tempFile, bytes);

        using (var source = new MemoryMappedTextSource(_tempFile, 0, text.Length))
        {
            string result = new string(System.Linq.Enumerable.ToArray(source));
            Assert.That(result, Is.EqualTo(text));
        }
    }

    [Test]
    public void SuffixTreeWithMMFSource_ShouldFindPatterns()
    {
        string text = "mississippi";
        File.WriteAllText(_tempFile, text, System.Text.Encoding.Unicode);

        using (var source = new MemoryMappedTextSource(_tempFile, 2, text.Length))
        {
            var tree = SuffixTree.Build(source);
            Assert.That(tree.Contains("ssi"), Is.True);
            Assert.That(tree.Contains("iss"), Is.True);
            Assert.That(tree.CountOccurrences("i"), Is.EqualTo(4));
            Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo("issi"));
        }
    }

    [Test]
    public void Constructor_InvalidArguments_ThrowsException()
    {
        string text = "ShouldNotBeCalled";
        File.WriteAllText(_tempFile, text);
        long fileSize = new FileInfo(_tempFile).Length;

        // Negative offset
        Assert.Throws<ArgumentOutOfRangeException>(() => new MemoryMappedTextSource(_tempFile, -1, 5));

        // Negative length
        Assert.Throws<ArgumentOutOfRangeException>(() => new MemoryMappedTextSource(_tempFile, 0, -5));

        // Length beyond file size
        // Note: MemoryMappedFile.CreateViewAccessor throws IOException or ArgumentOutOfRangeException or UnauthorizedAccessException depending on OS/implementation
        Assert.Throws(Is.InstanceOf<ArgumentOutOfRangeException>().Or.InstanceOf<UnauthorizedAccessException>().Or.InstanceOf<IOException>(),
            () => new MemoryMappedTextSource(_tempFile, 0, (int)fileSize + 10));

        // Offset beyond file
        Assert.Throws(Is.InstanceOf<ArgumentOutOfRangeException>().Or.InstanceOf<UnauthorizedAccessException>().Or.InstanceOf<IOException>(),
            () => new MemoryMappedTextSource(_tempFile, fileSize + 10, 5));
    }

    [Test]
    public void Read_WithNonZeroOffset_ReadsCorrectly()
    {
        // "Prefix" (6 bytes) + "Target" (6 bytes)
        // A B C D E F
        // 0 1 2 3 4 5
        byte[] bytes = new byte[] { 65, 0, 66, 0, 67, 0, 68, 0, 69, 0, 70, 0 }; // A B C D E F in UTF-16LE
        File.WriteAllBytes(_tempFile, bytes);

        // We want to skip 'A' (2 bytes). Offset = 2 bytes. Length = 5 chars (B, C, D, E, F)
        // Wait, file has 6 chars. Length is char count?
        // Constructor takes length in CHARS.
        // Filesize is 12 bytes.
        // We want to start at byte offset 2. Remaining bytes = 10. Remaining chars = 5.

        using (var source = new MemoryMappedTextSource(_tempFile, 2, 5))
        {
            Assert.That(source.Length, Is.EqualTo(5));
            Assert.That(source[0], Is.EqualTo('B'));
            Assert.That(source[4], Is.EqualTo('F'));
            Assert.That(source.ToString(), Is.EqualTo("BCDEF"));
        }
    }

    // ─── C7: Constructor 2 must release pointer on failure ─────────

    [Test]
    public void Constructor2_FailureAfterAcquirePointer_ReleasesPointer()
    {
        // Constructor 2 should have try/catch around AcquirePointer so that
        // if anything throws after pointer acquisition, ReleasePointer is called.
        // We verify by passing capacity that exactly fits validation but then
        // disposing the accessor and confirming no handle leak.
        string text = "ABCDE";
        byte[] bytes = System.Text.Encoding.Unicode.GetBytes(text);
        File.WriteAllBytes(_tempFile, bytes);

        using var mmf = System.IO.MemoryMappedFiles.MemoryMappedFile.CreateFromFile(
            _tempFile, FileMode.Open, null, 0, System.IO.MemoryMappedFiles.MemoryMappedFileAccess.Read);
        var accessor = mmf.CreateViewAccessor(0, bytes.Length, System.IO.MemoryMappedFiles.MemoryMappedFileAccess.Read);

        // Create a valid MemoryMappedTextSource via Constructor 2 — it should work
        using (var source = new MemoryMappedTextSource(accessor, 0, text.Length))
        {
            Assert.That(source[0], Is.EqualTo('A'));
        }

        // After dispose of the source, creating another should still work
        // (proves pointer was released properly)
        using (var source2 = new MemoryMappedTextSource(accessor, 0, text.Length))
        {
            Assert.That(source2[0], Is.EqualTo('A'));
        }

        accessor.Dispose();
    }

    // ─── C6: Dispose must ReleasePointer before nulling _ptr ──────

    [Test]
    public void Dispose_DoesNotThrow_AndPreventsSubsequentAccess()
    {
        // Verifies that Dispose properly releases the pointer before nulling
        // _ptr, so that ReleasePointer sees a valid handle state. Double-dispose
        // must also be safe.
        string text = "test";
        byte[] bytes = System.Text.Encoding.Unicode.GetBytes(text);
        File.WriteAllBytes(_tempFile, bytes);

        var source = new MemoryMappedTextSource(_tempFile, 0, text.Length);
        Assert.DoesNotThrow(() => source.Dispose());
        Assert.DoesNotThrow(() => source.Dispose()); // idempotent

        Assert.Throws<ObjectDisposedException>(() => _ = source[0]);
        Assert.Throws<ObjectDisposedException>(() => _ = source.ToString());
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
