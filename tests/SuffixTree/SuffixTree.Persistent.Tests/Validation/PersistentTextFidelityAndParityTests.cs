using System.Text;
namespace SuffixTree.Persistent.Tests.Validation;

/// <summary>
/// SUFFIXTREE-PERSISTENT (review 2026-09, B05): structure-level parity with the in-memory
/// <see cref="global::SuffixTree.SuffixTree"/> on edge-case texts across heap, MMF, reloaded-MMF
/// and hybrid (compact→large) builds, plus lossless text persistence/serialization (F16/F17).
/// Tied LRS/LCS representatives are implementation-specific (SUFFIXTREE-CORE) → lengths compared.
/// </summary>
[TestFixture]
[Category("Parity")]
public class PersistentTextFidelityAndParityTests
{
    private sealed class Collector : ISuffixTreeVisitor
    {
        public readonly List<string> Items = new();
        private readonly Stack<int> keys = new();
        public void VisitNode(int s, int e, int lc, int cc, int d) => Items.Add($"{s},{e},{lc},{cc},{d}");
        public void EnterBranch(int key) { keys.Push(key); Items.Add("K" + key); }
        public void ExitBranch() { keys.Pop(); Items.Add("X"); }
    }
    private static List<string> Topo(ISuffixTree t) { var c = new Collector(); t.Traverse(c); return c.Items; }

    private static IEnumerable<string> Texts()
    {
        yield return ""; yield return "a"; yield return "aaaaaaaaaa"; yield return "$"; yield return "a$b$"; yield return "\0\0a\0";
        yield return "￿￿a￿"; yield return "ab😀ab😀"; yield return "\uD800a\uD800"; yield return "é€ab€é";
        var rnd = new Random(42);
        for (int i = 0; i < 30; i++)
        {
            int n = rnd.Next(1, 300);
            var sb = new StringBuilder();
            for (int j = 0; j < n; j++) { sb.Append("ACGT"[rnd.Next(i % 2 == 0 ? 2 : 4)]); }
            yield return sb.ToString();
        }
    }

    private static void Compare(string text, ISuffixTree p, string label)
    {
        var r = global::SuffixTree.SuffixTree.Build(text);
        Assert.That(p.NodeCount, Is.EqualTo(r.NodeCount), label + " nodes");
        Assert.That(p.LeafCount, Is.EqualTo(r.LeafCount), label + " leaves");
        Assert.That(Topo(p), Is.EqualTo(Topo(r)), label + " topo");
        Assert.That(p.LongestRepeatedSubstring().Length, Is.EqualTo(r.LongestRepeatedSubstring().Length), label + " lrs");
        if (text.Length <= 400)
        {
            Assert.That(p.GetAllSuffixes(), Is.EqualTo(r.GetAllSuffixes()), label + " suffixes");
        }
        var q = text.Length > 2 ? string.Concat(text.AsSpan(1, Math.Min(text.Length - 1, 600)), "xy" + text[^Math.Min(text.Length, 600)..]) : "zz";
        Assert.That(p.LongestCommonSubstring(q).Length, Is.EqualTo(r.LongestCommonSubstring(q).Length), label + " lcs");
        Assert.That(p.FindExactMatchAnchors(q, 2).ToList(), Is.EqualTo(r.FindExactMatchAnchors(q, 2).ToList()), label + " anchors");
        // Long texts: probe the first and last 64 positions (covers chunk-boundary code units)
        foreach (int i in Enumerable.Range(0, text.Length).Where(i => text.Length <= 400 || i < 64 || i >= text.Length - 64))
            for (int len = 1; len <= 3 && i + len <= text.Length; len++)
            {
                var pat = text.Substring(i, len);
                Assert.That(p.FindAllOccurrences(pat).OrderBy(x => x), Is.EqualTo(r.FindAllOccurrences(pat).OrderBy(x => x)), label + " occ " + pat);
            }
    }

    [Test]
    public void Parity_Heap_Mmf_Hybrid_Reload()
    {
        int idx = 0;
        foreach (var t in Texts())
        {
            idx++;
            using (var h = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(t))) Compare(t, h, $"heap#{idx}");
            var path = Path.Combine(Path.GetTempPath(), $"probe_{Guid.NewGuid():N}.st");
            try
            {
                using (var m = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(t), path)) Compare(t, m, $"mmf#{idx}");
                using (var l = PersistentSuffixTreeFactory.LoadPersistent(path)) { Compare(t, l, $"load#{idx}"); Assert.That(l.Text.ToString(), Is.EqualTo(t)); }
            }
            finally { File.Delete(path); }
            if (t.Length > 3)
                using (var hy = PersistentSuffixTreeFactory.CreateCore(new StringTextSource(t), null, 88 + 24 * 3)) Compare(t, hy, $"hybrid#{idx}");
        }
    }

    // ── F16: empty text must round-trip through an MMF file (text region starts at EOF) ──
    [Test]
    public void MmfLoad_EmptyText_Loads()
    {
        var path = Path.Combine(Path.GetTempPath(), $"st_empty_{Guid.NewGuid():N}.st");
        try
        {
            using (PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(string.Empty), path)) { }
            using var loaded = PersistentSuffixTreeFactory.LoadPersistent(path);
            Assert.Multiple(() =>
            {
                Assert.That(loaded.IsEmpty, Is.True);
                Assert.That(loaded.Text.Length, Is.EqualTo(0));
                Assert.That(loaded.LeafCount, Is.EqualTo(0));
                Assert.That(loaded.Contains("a"), Is.False);
            });
        }
        finally { File.Delete(path); }
    }

    // ── F17: persisted text must be the exact code units the tree was built over ──
    private static IEnumerable<TestCaseData> NonWellFormedTexts()
    {
        yield return new TestCaseData("a\uD800b\uDC00a\uD800").SetArgDisplayNames("LoneSurrogates");
        yield return new TestCaseData(new string('A', 4095) + "\uD83D\uDE00AC\uD83D\uDE00")
            .SetArgDisplayNames("PairStraddles4096ChunkBoundary");
        yield return new TestCaseData(new string('G', 8191) + "\uD800\uDFFF").SetArgDisplayNames("PairAtSecondChunkEnd");
    }

    [TestCaseSource(nameof(NonWellFormedTexts))]
    public void Build_TextPersistedVerbatim_HeapReloadAndMmfReload(string text)
    {
        var storage = new HeapStorageProvider();
        new PersistentSuffixTreeBuilder(storage).Build(new StringTextSource(text));
        using (var heapLoaded = PersistentSuffixTree.Load(storage))
        {
            Assert.That(heapLoaded.Text.ToString(), Is.EqualTo(text), "heap reload text");
            Compare(text, heapLoaded, "heap reload");
        }

        var path = Path.Combine(Path.GetTempPath(), $"st_utf16_{Guid.NewGuid():N}.st");
        try
        {
            using (PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text), path)) { }
            using var loaded = PersistentSuffixTreeFactory.LoadPersistent(path);
            Assert.That(loaded.Text.ToString(), Is.EqualTo(text), "mmf reload text");
            Compare(text, loaded, "mmf reload");
        }
        finally { File.Delete(path); }
    }

    [TestCaseSource(nameof(NonWellFormedTexts))]
    public void Serializer_RoundTrip_PreservesCodeUnits(string text)
    {
        var reference = global::SuffixTree.SuffixTree.Build(text);
        using var ms = new MemoryStream();
        SuffixTreeSerializer.Export(reference, ms);
        ms.Position = 0;
        var imported = SuffixTreeSerializer.Import(ms, new HeapStorageProvider());
        Assert.That(imported.Text.ToString(), Is.EqualTo(text));
        Assert.That(SuffixTreeSerializer.CalculateLogicalHash(imported),
            Is.EqualTo(SuffixTreeSerializer.CalculateLogicalHash(reference)));
    }

    [Test]
    public void Export_WellFormedText_PayloadIsUtf16LE_FormatUnchanged()
    {
        const string text = "h\u00E9llo \uD83D\uDE00";
        using var ms = new MemoryStream();
        SuffixTreeSerializer.Export(global::SuffixTree.SuffixTree.Build(text), ms);
        byte[] all = ms.ToArray();
        byte[] expected = Encoding.Unicode.GetBytes(text);
        // magic(8) + version(4) + 7-bit length(1) = 13
        Assert.That(all.AsSpan(13, expected.Length).ToArray(), Is.EqualTo(expected));
        Assert.That(all.Length, Is.EqualTo(13 + expected.Length + 4 + 4 + 32));
    }

    private static byte[] ExportBytes(string text)
    {
        using var ms = new MemoryStream();
        SuffixTreeSerializer.Export(global::SuffixTree.SuffixTree.Build(text), ms);
        return ms.ToArray();
    }

    [Test]
    public void Import_CorruptedOrTruncatedTrailer_ThrowsInvalidDataException()
    {
        byte[] full = ExportBytes("banana"); // 13 + 12 text + 4 nodeCount + 4 hashLen + 32 hash
        int trailer = 13 + 12;

        byte[] noHashLen = full.AsSpan(0, trailer + 4).ToArray();
        byte[] shortHash = full.AsSpan(0, full.Length - 5).ToArray();
        byte[] badHashLen = (byte[])full.Clone();
        BitConverter.GetBytes(-1).CopyTo(badHashLen, trailer + 4);
        byte[] negLen = full.AsSpan(0, 12).ToArray().Concat(new byte[] { 0xFF, 0xFF, 0xFF, 0xFF, 0x0F }).ToArray();

        Assert.Multiple(() =>
        {
            foreach (var bytes in new[] { noHashLen, shortHash, badHashLen, negLen })
                Assert.Throws<InvalidDataException>(() =>
                    SuffixTreeSerializer.Import(new MemoryStream(bytes), new HeapStorageProvider()));
        });
    }
}

/// <summary>
/// F18: text payload allocation must not wrap for &gt; 2 GiB (UTF-16 text &gt; 2^30 − 1 chars).
/// Verified with a recording provider — no real memory is allocated.
/// </summary>
[TestFixture]
[Category("Safety")]
public class TextPayloadAllocationOverflowTests
{
    private sealed class RecordingBumpProvider : IStorageProvider
    {
        public readonly List<int> Requests = new();
        private long _position = 1000;
        public long Size => _position;
        public long Allocate(int size)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(size);
            Requests.Add(size);
            long o = _position;
            _position += size;
            return o;
        }
        public void EnsureCapacity(long capacity) { }
        public int ReadInt32(long offset) => 0;
        public void WriteInt32(long offset, int value) { }
        public uint ReadUInt32(long offset) => 0;
        public void WriteUInt32(long offset, uint value) { }
        public long ReadInt64(long offset) => 0;
        public void WriteInt64(long offset, long value) { }
        public char ReadChar(long offset) => '\0';
        public void WriteChar(long offset, char value) { }
        public void ReadBytes(long offset, byte[] buffer, int start, int count) { }
        public void WriteBytes(long offset, byte[] buffer, int start, int count) { }
        public void Dispose() { }
    }

    [Test]
    public void AllocateContiguous_Utf16TextOver2GiB_DoesNotWrap()
    {
        const int textLength = (1 << 30) + 5;          // UTF-16 chars
        long bytes = (long)textLength * sizeof(char);   // 2^31 + 10 bytes > int.MaxValue
        Assert.That(unchecked((int)bytes), Is.LessThan(0), "precondition: old (int) cast wrapped negative");

        var provider = new RecordingBumpProvider();
        long offset = PersistentSuffixTreeBuilder.AllocateContiguous(provider, bytes);

        Assert.Multiple(() =>
        {
            Assert.That(offset, Is.EqualTo(1000));
            Assert.That(provider.Requests.Sum(r => (long)r), Is.EqualTo(bytes));
            Assert.That(provider.Requests, Is.EqualTo(new[] { int.MaxValue, 11 }));
            Assert.That(provider.Size, Is.EqualTo(1000 + bytes));
        });
    }

    [TestCase(0L)]
    [TestCase(24L)]
    public void AllocateContiguous_SmallSizes_SingleAllocation(long bytes)
    {
        var provider = new RecordingBumpProvider();
        Assert.That(PersistentSuffixTreeBuilder.AllocateContiguous(provider, bytes), Is.EqualTo(1000));
        Assert.That(provider.Requests, Is.EqualTo(new[] { (int)bytes }));
    }
}

/// <summary>F18b: file-backed UTF-16 text source ≥ 2^30 chars (sparse file, no real I/O).</summary>
[TestFixture]
[Category("Safety")]
public class MemoryMappedTextSourceLargeLengthTests
{
    [Test]
    public void FileCtor_LengthOver2Pow30Chars_DoesNotOverflow()
    {
        const int length = (1 << 30) + 5;
        var path = Path.Combine(Path.GetTempPath(), $"mmts_large_{Guid.NewGuid():N}.bin");
        try
        {
            using (var fs = new FileStream(path, FileMode.CreateNew)) fs.SetLength((long)length * sizeof(char));
            using var src = new MemoryMappedTextSource(path, 0, length);
            Assert.Multiple(() =>
            {
                Assert.That(src.Length, Is.EqualTo(length));
                Assert.That(src[length - 1], Is.EqualTo('\0'));
            });
        }
        finally { File.Delete(path); }
    }
}
