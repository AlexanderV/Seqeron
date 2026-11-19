using System.Buffers.Binary;
using System.Diagnostics;
using System.Globalization;
using System.IO.MemoryMappedFiles;
using Microsoft.Win32.SafeHandles;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// Gates G2, G3 and G5 of work package 1b.7f at their real sizes (outside the default run: <c>[Explicit]</c>, driven by <c>analysis/em/run_1b7f_g2.sh</c> and <c>run_1b7f_g3.sh</c>, which set the <c>SESA_*</c> variables and run one variant per
/// process, cold cache where the protocol says so). Variants (<c>SESA_VARIANT</c>): <c>open</c> (the open of <c>SESA_FILE</c> at <c>SESA_LEVEL</c> / <c>SESA_ACCESS</c>: wall time and RssAnon before and at the peak, the report of the
/// validation; <c>SESA_CHUNK</c> sets the chunk of the structure level), <c>identity</c> (sampled counts and locate sets of 1,000 patterns of the text and 1,000 patterns that the oracle finds absent against an independent oracle that does a binary search over the suffix array and the input text, no engine),
/// <c>corrupt</c> (a damaged section with its checksum recomputed must be caught by the structure level), <c>throughput</c> (count, locate and a child walk, in memory against the file at <c>Checked</c> and <c>Trusted</c>, warm).
/// </summary>
[TestFixture]
[Category("Cost")]
[Explicit("gates G2 / G3 / G5 of 1b.7f: run by analysis/em/run_1b7f_g2.sh and run_1b7f_g3.sh, which set the SESA_* variables")]
public class SesaReaderCostTests
{
    private static string Env(string name) => Environment.GetEnvironmentVariable(name) ?? throw new InvalidOperationException($"{name} is not set");

    private static string EnvOr(string name, string fallback) => Environment.GetEnvironmentVariable(name) ?? fallback;

    private static string F(double v) => v.ToString("0.#####", CultureInfo.InvariantCulture);

    private static void Row(string kind, string label, params object[] fields)
    {
        string line = string.Join('\t', new object[] { kind, label }.Concat(fields).Select(x => x is double d ? F(d) : Convert.ToString(x, CultureInfo.InvariantCulture)));
        File.AppendAllText(Env("SESA_ROWS"), line + "\n");
    }

    private static long RssAnonKb()
    {
        foreach (string line in File.ReadLines("/proc/self/status"))
        {
            if (line.StartsWith("RssAnon:", StringComparison.Ordinal)) return long.Parse(line.Split(' ', StringSplitOptions.RemoveEmptyEntries)[1], CultureInfo.InvariantCulture);
        }
        return -1;
    }

    [Test]
    [CancelAfter(7_200_000)]
    public void Run()
    {
        switch (Env("SESA_VARIANT"))
        {
            case "open": OpenVariant(); break;
            case "identity": Identity(); break;
            case "corrupt": Corrupt(); break;
            case "throughput": Throughput(); break;
            default: throw new InvalidOperationException("SESA_VARIANT must be open, identity, corrupt or throughput");
        }
    }

    // ---------------------------------------------------------------- open

    private static void OpenVariant()
    {
        string file = Env("SESA_FILE"), label = Env("SESA_LABEL");
        var level = Enum.Parse<EsaValidationLevel>(EnvOr("SESA_LEVEL", "Header"));
        var access = Enum.Parse<EsaAccessMode>(EnvOr("SESA_ACCESS", "Checked"));
        bool sha = EnvOr("SESA_SHA", "0") == "1";
        if (EnvOr("SESA_CHUNK", "") is { Length: > 0 } chunk) SesaVerifierTuning.ChunkSymbols = int.Parse(chunk, CultureInfo.InvariantCulture); // G5: the time and the memory of the structure level against the size of a chunk of the Phi pass
        long before = RssAnonKb();
        long max = before;
        using var stop = new CancellationTokenSource();
        var sampler = new Thread(() =>
        {
            while (!stop.IsCancellationRequested)
            {
                long v = RssAnonKb();
                if (v > max) max = v;
                Thread.Sleep(10);
            }
        }) { IsBackground = true };
        sampler.Start();
        long t0 = Stopwatch.GetTimestamp();
        using EsaSuffixTree tree = SesaTestFile.Open(file, level, access, sha);
        double seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        stop.Cancel();
        sampler.Join();
        long after = RssAnonKb();
        max = Math.Max(max, after);
        SesaOpenReport r = tree.File!.Report;
        Row("open", label, level, access, seconds, r.HeaderSeconds, r.SectionsSeconds, r.StructureSeconds, r.HardwareCrc ? "hw" : "table", r.SectionBytesChecksummed, r.BitmapBytes, r.BufferBytes, before, max, after, tree.File.Header.TextLength);
        if (EnvOr("SESA_WARM", "0") == "1")
        {
            // the same open again in the same process: the page cache is warm and the code is compiled (what the first call paid for the JIT is the difference)
            long t1 = Stopwatch.GetTimestamp();
            using EsaSuffixTree again = SesaTestFile.Open(file, level, access, sha);
            Row("open-warm", label, level, access, (Stopwatch.GetTimestamp() - t1) / (double)Stopwatch.Frequency);
        }
        Assert.That(tree.File.Header.TextLength, Is.GreaterThan(0));
    }

    // ---------------------------------------------------------------- identity at 10^8: an oracle that is not the engine

    /// <summary>
    /// The oracle: the suffix array is read from the container by a header parser of its own (offsets from the table, nothing of the product), the text is the INPUT file mapped as it is (bytes), and a pattern is looked up by two
    /// binary searches that compare it with the suffixes at the rank in the middle.
    /// </summary>
    private sealed unsafe class BinarySearchOracle : IDisposable
    {
        private readonly MemoryMappedFile _containerMap, _textMap;
        private readonly MemoryMappedViewAccessor _containerView, _textView;
        private readonly byte* _container;
        private readonly byte* _text;
        private readonly int* _sa;

        internal BinarySearchOracle(string container, string input)
        {
            _containerMap = MemoryMappedFile.CreateFromFile(container, FileMode.Open, null, 0, MemoryMappedFileAccess.Read);
            _containerView = _containerMap.CreateViewAccessor(0, 0, MemoryMappedFileAccess.Read);
            byte* p = null;
            _containerView.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
            _container = p + _containerView.PointerOffset;
            _textMap = MemoryMappedFile.CreateFromFile(input, FileMode.Open, null, 0, MemoryMappedFileAccess.Read);
            _textView = _textMap.CreateViewAccessor(0, 0, MemoryMappedFileAccess.Read);
            byte* q = null;
            _textView.SafeMemoryMappedViewHandle.AcquirePointer(ref q);
            _text = q + _textView.PointerOffset;
            N = BinaryPrimitives.ReadInt64LittleEndian(new ReadOnlySpan<byte>(_container + 16, 8));
            int count = BinaryPrimitives.ReadInt32LittleEndian(new ReadOnlySpan<byte>(_container + 48, 4));
            long saOffset = -1;
            for (int i = 0; i < count; i++)
            {
                var e = new ReadOnlySpan<byte>(_container + 128 + (64 * i), 64);
                if (BinaryPrimitives.ReadInt32LittleEndian(e) == 3) saOffset = BinaryPrimitives.ReadInt64LittleEndian(e[8..]);
            }
            _sa = (int*)(_container + saOffset);
            using (SafeFileHandle h = File.OpenHandle(input))
            {
                if (RandomAccess.GetLength(h) != N) throw new InvalidOperationException("the input text and the container differ in length"); // through the handle: the input may be a symbolic link, whose FileInfo length is the link's
            }
        }

        internal long N { get; }

        internal long TextByte(long i) => _text[i];

        /// <summary>Compares the pattern with the suffix of rank <paramref name="rank"/> over the pattern's length: negative when the suffix is smaller, 0 when the pattern is a prefix of it.</summary>
        private int Compare(string pattern, int rank)
        {
            long start = _sa[rank];
            for (int k = 0; k < pattern.Length; k++)
            {
                if (start + k >= N) return -1;
                int d = _text[start + k] - pattern[k];
                if (d != 0) return d < 0 ? -1 : 1;
            }
            return 0;
        }

        /// <summary>The ranks <c>[lo, hi)</c> of the suffixes that start with the pattern.</summary>
        internal (long Lo, long Hi) Range(string pattern)
        {
            long lo = 1, hi = N + 1;
            while (lo < hi)
            {
                long mid = lo + ((hi - lo) / 2);
                if (Compare(pattern, (int)mid) < 0) lo = mid + 1;
                else hi = mid;
            }
            long first = lo;
            hi = N + 1;
            while (lo < hi)
            {
                long mid = lo + ((hi - lo) / 2);
                if (Compare(pattern, (int)mid) <= 0) lo = mid + 1;
                else hi = mid;
            }
            return (first, lo);
        }

        internal int Sa(long rank) => _sa[rank];

        public void Dispose()
        {
            _containerView.SafeMemoryMappedViewHandle.ReleasePointer();
            _containerView.Dispose();
            _containerMap.Dispose();
            _textView.SafeMemoryMappedViewHandle.ReleasePointer();
            _textView.Dispose();
            _textMap.Dispose();
        }
    }

    private static unsafe void Identity()
    {
        string file = Env("SESA_FILE"), input = Env("SESA_INPUT"), label = Env("SESA_LABEL");
        var level = Enum.Parse<EsaValidationLevel>(EnvOr("SESA_LEVEL", "Header"));
        using EsaSuffixTree tree = SesaTestFile.Open(file, level, EsaAccessMode.Checked, sha: false);
        using var oracle = new BinarySearchOracle(file, input);
        long n = oracle.N;
        var rng = new SplitMix64(0x1B7F2);
        var patterns = new List<(string Text, bool Present)>();
        for (int i = 0; i < 1_000; i++)
        {
            int len = 1 + rng.NextInt(40);
            long at = (long)(rng.Next() % (ulong)(n - len));
            var chars = new char[len];
            for (int k = 0; k < len; k++) chars[k] = (char)oracle.TextByte(at + k);
            patterns.Add((new string(chars), true));
        }
        var alphabet = new List<char>();
        for (long i = 0; i < Math.Min(n, 100_000); i++) if (!alphabet.Contains((char)oracle.TextByte(i))) alphabet.Add((char)oracle.TextByte(i));
        // 1,000 patterns that the ORACLE says are absent: random strings over the alphabet of the text (a short one is usually present and is not kept)
        int kept = 0, tries = 0;
        while (kept < 1_000)
        {
            Assert.That(++tries, Is.LessThan(1_000_000), "too few absent random patterns");
            int len = 1 + rng.NextInt(40);
            var chars = new char[len];
            for (int k = 0; k < len; k++) chars[k] = alphabet[rng.NextInt(alphabet.Count)];
            string candidate = new(chars);
            (long clo, long chi) = oracle.Range(candidate);
            if (chi != clo) continue;
            patterns.Add((candidate, false));
            kept++;
        }
        int countMismatch = 0, locateMismatch = 0, absent = 0, locatedFully = 0, locatedPrefix = 0;
        long t0 = Stopwatch.GetTimestamp();
        var positions = new long[1_000];
        foreach ((string text, bool present) in patterns)
        {
            (long lo, long hi) = oracle.Range(text);
            long expected = hi - lo;
            if (present && expected == 0) Assert.Fail($"a pattern taken from the text has no occurrence by the oracle: '{text}'");
            if (expected == 0) absent++;
            long count = tree.CountOccurrences(text);
            if (count != expected)
            {
                countMismatch++;
                continue;
            }
            if (expected == 0) continue;
            if (expected <= 100_000)
            {
                long[] want = [.. Enumerable.Range(0, (int)expected).Select(k => (long)oracle.Sa(lo + k)).Order()];
                long[] got = [.. tree.FindAllOccurrences(text).Select(x => (long)x).Order()];
                if (!want.SequenceEqual(got)) locateMismatch++;
                locatedFully++;
            }
            else
            {
                // a pattern of tens of millions of occurrences: count equal, and the first 1,000 locate results in suffix-array order equal the oracle's slice of the array
                IndexLocus locus = ((ISuffixIndex)tree).Root;
                foreach (char c in text)
                {
                    Assert.That(((ISuffixIndex)tree).TryExtendRight(locus, c, out IndexLocus next), Is.True);
                    locus = next;
                }
                int got = ((ISuffixIndex)tree).Locate(locus, 0, positions);
                bool same = got == positions.Length;
                for (int k = 0; same && k < got; k++) same = positions[k] == oracle.Sa(lo + k);
                if (!same) locateMismatch++;
                locatedPrefix++;
            }
        }
        double seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        Row("identity", label, level, patterns.Count, absent, locatedFully, locatedPrefix, countMismatch, locateMismatch, seconds);
        Assert.That(countMismatch, Is.Zero, "counts that differ from the oracle");
        Assert.That(locateMismatch, Is.Zero, "locate sets that differ from the oracle");
        Assert.That(absent, Is.EqualTo(1_000), "the oracle found the 1,000 absent patterns absent");
    }

    // ---------------------------------------------------------------- a damaged body whose checksum was recomputed

    private static unsafe void Corrupt()
    {
        string file = Env("SESA_FILE"), label = Env("SESA_LABEL");
        byte[] page = new byte[4096];
        using (var fs = new FileStream(file, FileMode.Open, FileAccess.Read)) fs.ReadExactly(page);
        int count = BinaryPrimitives.ReadInt32LittleEndian(page.AsSpan(48));
        long saOffset = 0, childOffset = 0;
        for (int i = 0; i < count; i++)
        {
            ReadOnlySpan<byte> e = page.AsSpan(128 + (64 * i), 64);
            int kind = BinaryPrimitives.ReadInt32LittleEndian(e);
            if (kind == 3) saOffset = BinaryPrimitives.ReadInt64LittleEndian(e[8..]);
            if (kind == 6) childOffset = BinaryPrimitives.ReadInt64LittleEndian(e[8..]);
        }
        long n = BinaryPrimitives.ReadInt64LittleEndian(page.AsSpan(16));
        long rank = n / 2;
        var cases = new (string Name, long Offset, int Length, Action<Span<byte>> Edit)[]
        {
            ("sa-adjacent-swap", saOffset + (4 * rank), 8, b => { Span<byte> t = stackalloc byte[4]; b[..4].CopyTo(t); b.Slice(4, 4).CopyTo(b); t.CopyTo(b.Slice(4, 4)); }),
            ("child-cell", childOffset + (4 * (n / 3)), 4, b => b[0] ^= 0x01),
        };
        foreach ((string name, long offset, int length, Action<Span<byte>> edit) in cases)
        {
            var saved = new byte[length];
            using (var fs = new FileStream(file, FileMode.Open, FileAccess.ReadWrite))
            {
                fs.Position = offset;
                fs.ReadExactly(saved);
                var changed = (byte[])saved.Clone();
                edit(changed);
                if (changed.AsSpan().SequenceEqual(saved)) Assert.Inconclusive($"{name}: the edit changed nothing");
                fs.Position = offset;
                fs.Write(changed);
            }
            ResealWholeFile(file);
            long t0 = Stopwatch.GetTimestamp();
            InvalidDataException? caught = null;
            try
            {
                using EsaSuffixTree tree = SesaTestFile.Open(file, EsaValidationLevel.Structure);
            }
            catch (InvalidDataException e)
            {
                caught = e;
            }
            double seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
            // restore before asserting, so that a failure leaves the container as it was
            using (var fs = new FileStream(file, FileMode.Open, FileAccess.ReadWrite))
            {
                fs.Position = offset;
                fs.Write(saved);
            }
            ResealWholeFile(file);
            Row("corrupt", label, name, caught is null ? "NOT CAUGHT" : "caught", seconds, caught?.Message.Replace('\t', ' ').Replace('\n', ' ') ?? string.Empty);
            Assert.That(caught, Is.Not.Null, $"{name}: the structure level did not see the damage");
            using EsaSuffixTree sums = SesaTestFile.Open(file, EsaValidationLevel.Sections);
        }
    }

    /// <summary>Recomputes the CRC-32C of every section and the header CRC of a container file in place (one pass over the mapping).</summary>
    private static unsafe void ResealWholeFile(string path)
    {
        using var map = MemoryMappedFile.CreateFromFile(path, FileMode.Open, null, 0, MemoryMappedFileAccess.ReadWrite);
        using MemoryMappedViewAccessor view = map.CreateViewAccessor(0, 0, MemoryMappedFileAccess.ReadWrite);
        byte* p = null;
        view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
        try
        {
            byte* data = p + view.PointerOffset;
            int count = BinaryPrimitives.ReadInt32LittleEndian(new ReadOnlySpan<byte>(data + 48, 4));
            for (int i = 0; i < count; i++)
            {
                byte* e = data + 128 + (64 * i);
                int flags = BinaryPrimitives.ReadInt32LittleEndian(new ReadOnlySpan<byte>(e + 4, 4));
                if ((flags & 2) == 0) continue;
                long offset = BinaryPrimitives.ReadInt64LittleEndian(new ReadOnlySpan<byte>(e + 8, 8));
                long length = BinaryPrimitives.ReadInt64LittleEndian(new ReadOnlySpan<byte>(e + 16, 8));
                BinaryPrimitives.WriteUInt32LittleEndian(new Span<byte>(e + 56, 4), Crc32C.Append(0, data + offset, length));
            }
            uint crc = Crc32C.Append(Crc32C.Append(0, data, 124), data + 128, 64L * count);
            BinaryPrimitives.WriteUInt32LittleEndian(new Span<byte>(data + 124, 4), crc);
            view.Flush();
        }
        finally
        {
            view.SafeMemoryMappedViewHandle.ReleasePointer();
        }
    }

    // ---------------------------------------------------------------- throughput

    private static void Throughput()
    {
        string file = Env("SESA_FILE"), input = Env("SESA_INPUT"), label = Env("SESA_LABEL");
        int repeats = int.Parse(EnvOr("SESA_REPEATS", "5"), CultureInfo.InvariantCulture);
        string text = File.ReadAllText(input, System.Text.Encoding.Latin1);
        // SESA_TREES selects the trees of a process (memory, checked, trusted; default all three in this order): the protocol of G3 runs one tree per process, in alternating orders, so that the order inside a process is not a bias
        string[] wanted = EnvOr("SESA_TREES", "memory,checked,trusted").Split(',', StringSplitOptions.RemoveEmptyEntries);
        using EsaSuffixTree? memory = wanted.Contains("memory") ? EsaSuffixTree.Build(text) : null;
        using EsaSuffixTree? checkedTree = wanted.Contains("checked") ? SesaTestFile.Open(file, EsaValidationLevel.Header, EsaAccessMode.Checked) : null;
        using EsaSuffixTree? trusted = wanted.Contains("trusted") ? SesaTestFile.Open(file, EsaValidationLevel.Header, EsaAccessMode.Trusted) : null;
        var rng = new SplitMix64(1234);
        string[] twenty = [.. Enumerable.Range(0, 200_000).Select(_ => text.Substring(rng.NextInt(text.Length - 20), 20))];
        string[] eight = [.. Enumerable.Range(0, 1_000).Select(_ => text.Substring(rng.NextInt(text.Length - 8), 8))];
        string[] absent = [.. Enumerable.Range(0, 100_000).Select(_ => new string(Enumerable.Range(0, 24).Select(_ => "ACGT"[rng.NextInt(4)]).ToArray()))];
        var trees = new (string Name, EsaSuffixTree? Tree)[] { ("memory", memory), ("checked", checkedTree), ("trusted", trusted) };
        foreach ((string name, EsaSuffixTree? tree) in trees)
        {
            if (tree is null) continue;
            Warm(tree, twenty, eight, absent);
            var count20 = new List<double>();
            var locate8 = new List<double>();
            var count24 = new List<double>();
            var walk = new List<double>();
            for (int r = 0; r < repeats; r++)
            {
                long t = Stopwatch.GetTimestamp();
                long sum = 0;
                foreach (string p in twenty) sum += tree.CountOccurrences(p);
                count20.Add(twenty.Length / Seconds(t));
                t = Stopwatch.GetTimestamp();
                long located = 0;
                foreach (string p in eight) located += tree.FindAllOccurrences(p).Count;
                locate8.Add(eight.Length / Seconds(t));
                t = Stopwatch.GetTimestamp();
                foreach (string p in absent) sum += tree.CountOccurrences(p);
                count24.Add(absent.Length / Seconds(t));
                t = Stopwatch.GetTimestamp();
                long nodes = ChildWalk(tree, 9);
                walk.Add(nodes / Seconds(t));
                GC.KeepAlive(sum + located + nodes);
            }
            Row("throughput", label, name, Median(count20), Median(locate8), Median(count24), Median(walk), repeats);
        }
    }

    private static double Seconds(long t0) => (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;

    private static double Median(List<double> v)
    {
        double[] s = [.. v.Order()];
        return s[s.Length / 2];
    }

    private static void Warm(EsaSuffixTree tree, string[] a, string[] b, string[] c)
    {
        foreach (string p in a.Take(20_000)) _ = tree.CountOccurrences(p);
        foreach (string p in b) _ = tree.FindAllOccurrences(p).Count;
        foreach (string p in c.Take(20_000)) _ = tree.CountOccurrences(p);
        _ = ChildWalk(tree, 9);
    }

    /// <summary>Visits every node of the top levels (the contract's child walk): the number of nodes visited.</summary>
    private static long ChildWalk(ISuffixTreeNavigation tree, int depth)
    {
        long visited = 0;
        var buffer = new IndexLocus[8];
        void Go(IndexLocus node, int d)
        {
            visited++;
            if (d == 0) return;
            int count = tree.GetChildren(node, buffer);
            var mine = new IndexLocus[count];
            Array.Copy(buffer, mine, Math.Min(count, buffer.Length));
            for (int i = 0; i < Math.Min(count, mine.Length); i++) Go(mine[i], d - 1);
        }
        Go(tree.Root, depth);
        return visited;
    }
}
