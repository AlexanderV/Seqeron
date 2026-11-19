using System.Buffers.Binary;
using System.Runtime.InteropServices;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The <c>Checked</c> access mode (work package 1b.7f, gate G1; ADR-003 section 6): a body that passes the header and was damaged afterwards is an <see cref="InvalidDataException"/> or a plain answer, never a read outside the
/// mapping, a hang or another exception. The file is mapped between two inaccessible pages (an anonymous mapping and <c>mprotect</c>), so that a read one byte outside the container ends the test host.
/// </summary>
[TestFixture]
public class SesaReaderHostileTests
{
    /// <summary>The image in memory with an inaccessible page before and after the pages that hold it; <paramref name="flush"/> puts its last byte right against the page behind, else its first byte right behind the page before.</summary>
    private sealed unsafe class Guarded : IDisposable
    {
        private const int ProtNone = 0, ProtRead = 1, ProtWrite = 2, MapPrivate = 2, MapAnonymous = 0x20;

        [DllImport("libc", SetLastError = true)]
        private static extern nint mmap(nint addr, nuint length, int prot, int flags, int fd, nint offset);

        [DllImport("libc", SetLastError = true)]
        private static extern int mprotect(nint addr, nuint length, int prot);

        [DllImport("libc", SetLastError = true)]
        private static extern int munmap(nint addr, nuint length);

        private readonly nint _base;
        private readonly nuint _size;

        internal Guarded(byte[] image, bool flush)
        {
            const int page = 4096;
            int pages = (image.Length + page - 1) / page;
            _size = (nuint)((pages + 2) * page);
            _base = mmap(0, _size, ProtRead | ProtWrite, MapPrivate | MapAnonymous, -1, 0);
            if (_base == -1) throw new InvalidOperationException("mmap failed");
            byte* data = (byte*)_base + page + (flush ? (pages * page) - image.Length : 0);
            image.CopyTo(new Span<byte>(data, image.Length));
            if (mprotect(_base, page, ProtNone) != 0 || mprotect(_base + (nint)((pages + 1) * page), page, ProtNone) != 0) throw new InvalidOperationException("mprotect failed");
            Data = data;
            Length = image.Length;
        }

        internal byte* Data { get; }

        internal int Length { get; }

        public void Dispose() => munmap(_base, _size);
    }

    private sealed class Counting : ISuffixTreeVisitor
    {
        internal long Nodes;

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) => Nodes++;

        public void EnterBranch(int key)
        {
        }

        public void ExitBranch()
        {
        }
    }

    /// <summary>Runs every kind of operation on the tree; any exception other than <see cref="InvalidDataException"/> (or the plain rejection of a locus that a damaged table produced) is a failure.</summary>
    private static void Battery(EsaSuffixTree tree, string text, string label)
    {
        int ok = 0, bad = 0;
        void Step(string what, Action action)
        {
            if (Environment.GetEnvironmentVariable("ESA_STEP_LOG") is { } logPath) File.AppendAllText(logPath, $"{label}: {what}\n");
            try
            {
                action();
                ok++;
            }
            catch (InvalidDataException)
            {
                bad++;
            }
            catch (ArgumentException e) when (e.GetType() == typeof(ArgumentException) && what.StartsWith("nav", StringComparison.Ordinal))
            {
                bad++; // a locus that a damaged table produced fails the O(1) checks of the next call
            }
            catch (Exception e)
            {
                Assert.Fail($"{label}: {what} threw {e.GetType().Name}: {e.Message}");
            }
        }

        var rng = new SplitMix64((ulong)text.Length);
        var patterns = new List<string> { "a", "ab", "q", text.Length > 3 ? text[1..4] : "a" };
        for (int i = 0; i < 12 && text.Length > 8; i++) patterns.Add(text.Substring(rng.NextInt(text.Length - 6), 1 + rng.NextInt(6)));
        foreach (string p in patterns)
        {
            Step("count", () => tree.CountOccurrences(p));
            Step("contains", () => tree.Contains(p));
            Step("locate", () => tree.FindAllOccurrences(p));
        }
        Step("lrs", () => tree.LongestRepeatedSubstring());
        Step("lrs memory", () => tree.LongestRepeatedSubstringMemory());
        Step("all lrs", () => tree.FindAllLongestRepeatedSubstrings());
        Step("distinct", () => tree.CountDistinctSubstrings());
        Step("distinct by length", () => tree.CountDistinctSubstringsByLength(9));
        Step("suffixes", () => tree.GetAllSuffixes());
        Step("enumerate suffixes", () => tree.EnumerateSuffixes().Take(5).ToList());
        Step("traverse", () => tree.Traverse(new Counting()));
        Step("print", () => tree.PrintTree());
        Step("lcs", () => tree.LongestCommonSubstring(text.Length > 20 ? text[3..20] : "ab"));
        Step("mem", () => tree.FindMaximalExactMatches(text.Length > 30 ? text[2..30] + "q" + text[5..15] : "abab", 3));
        Step("mum", () => tree.FindMaximalUniqueMatches(text.Length > 30 ? text[2..30] : "abab", 3));
        Step("anchors", () => tree.FindExactMatchAnchors(text.Length > 30 ? text[1..30] : "abab", 3));

        ISuffixTreeNavigation nav = tree;
        long rank = Math.Min(3, tree.Text.Length);
        ISuffixArrayAccess sa = tree;
        Step("nav suffix", () => sa.SuffixAt(rank));
        Step("nav lcp", () => sa.LcpAt(rank));
        Step("nav rank of", () => sa.RankOf(rank));
        Step("nav lcp pair", () => sa.Lcp(1, Math.Max(1, tree.Text.Length)));
        Step("nav lca", () => sa.Lca(1, Math.Max(1, tree.Text.Length)));
        Step("nav walk", () => Walk(nav, nav.Root, 0, 6));
        Step("nav extend left", () => nav.TryExtendLeft(nav.Root, 'a', out _));
        Assert.That(ok + bad, Is.GreaterThan(0), label);
    }

    private static void Walk(ISuffixTreeNavigation nav, IndexLocus node, int depth, int maxDepth)
    {
        if (depth >= maxDepth) return;
        var children = new IndexLocus[64];
        int count = nav.GetChildren(node, children);
        Assert.That(count, Is.GreaterThanOrEqualTo(0));
        for (int i = 0; i < Math.Min(count, children.Length); i++)
        {
            _ = nav.TryParent(children[i], out _);
            _ = nav.TrySuffixLink(children[i], out _);
            Walk(nav, children[i], depth + 1, maxDepth);
        }
    }

    private static void RunWithWatchdog(Action action, string label, int seconds = 30)
    {
        Exception? failure = null;
        var thread = new Thread(() =>
        {
            try
            {
                action();
            }
            catch (Exception e)
            {
                failure = e;
            }
        }) { IsBackground = true };
        thread.Start();
        Assert.That(thread.Join(TimeSpan.FromSeconds(seconds)), Is.True, $"{label}: hangs");
        if (failure is not null) throw new AssertionException($"{label}: {failure.GetType().Name}: {failure.Message}", failure);
    }

    /// <summary>The sealed random damage of one section: a share of its words is replaced by random or extreme values.</summary>
    private static byte[] Hostile(byte[] clean, SplitMix64 rng, int kind)
    {
        byte[] image = (byte[])clean.Clone();
        SesaOracle o = SesaOracle.FromBytes(image);
        OracleSection s = o.Get(kind);
        int words = (int)(s.Length / 4);
        if (words == 0) return image;
        int hits = 1 + rng.NextInt(Math.Max(1, words / 5));
        long n = o.TextLen;
        uint[] extreme = [0, 1, (uint)n, (uint)n + 1, (uint)n + 2, 0x7FFFFFFF, 0x80000000, 0xFFFFFFFF, 0xFFFFFFFE, 0x40000000, 0x3FFFFFFF, 0xC0000001];
        for (int h = 0; h < hits; h++)
        {
            int w = rng.NextInt(words);
            uint v = rng.NextInt(3) switch { 0 => extreme[rng.NextInt(extreme.Length)], 1 => (uint)rng.Next(), _ => (uint)rng.NextInt((int)Math.Max(2, n + 2)) };
            BinaryPrimitives.WriteUInt32LittleEndian(image.AsSpan((int)s.Offset + (4 * w)), v);
        }
        SesaTestFile.Reseal(image);
        return image;
    }

    [Test]
    public void ADamagedBody_InCheckedMode_NeverLeavesTheMapping_NeverHangs_AndFailsOnlyAsInvalidData()
    {
        if (!OperatingSystem.IsLinux()) Assert.Ignore("the guard pages need Linux mmap / mprotect");
        string text = StringFamilies.Random(150, 3, 5) + "abababababababababababab" + new string('a', 60);
        text = new string(text.Select(c => (char)('a' + ((c - 'a') % 3))).ToArray());
        LcpChildCase c = LcpChildCase.Of(text);
        byte[] clean = SesaTestFile.Image(c);
        var rng = new SplitMix64(0xBADC0DE);
        int[] kinds = [SesaFormat.KindSa, SesaFormat.KindChild, SesaFormat.KindLcp, SesaFormat.KindPsi];
        var options = new EsaOpenOptions { Validation = EsaValidationLevel.Header, Access = EsaAccessMode.Checked };
        for (int round = 0; round < 240; round++)
        {
            int kind = kinds[round % kinds.Length];
            byte[] image = Hostile(clean, rng, kind);
            using var guarded = new Guarded(image, flush: round % 2 == 0);
            int r = round;
            RunWithWatchdog(() =>
            {
                unsafe
                {
                    using SesaFile file = SesaFile.FromMemory(guarded.Data, guarded.Length, options);
                    using EsaSuffixTree tree = EsaSuffixTree.FromFile(file);
                    Battery(tree, text, $"round {r} (kind {kind}, {(r % 2 == 0 ? "flush" : "leading")})");
                }
            }, $"round {round} kind {kind}");
        }
    }

    [Test]
    public void ADamagedOverflowList_InCheckedMode_IsInvalidData_NotAWrongAnswerOrACrash()
    {
        if (!OperatingSystem.IsLinux()) Assert.Ignore("the guard pages need Linux mmap / mprotect");
        string random = StringFamilies.Random(1_500, 4, 21);
        string text = random + random[..300];
        LcpChildCase c = LcpChildCase.Of(text);
        byte[] clean = SesaTestFile.Image(c);
        var rng = new SplitMix64(77);
        var options = new EsaOpenOptions { Validation = EsaValidationLevel.Header, Access = EsaAccessMode.Checked };
        for (int round = 0; round < 60; round++)
        {
            byte[] image = Hostile(clean, rng, SesaFormat.KindLcpOverflow);
            using var guarded = new Guarded(image, flush: round % 2 == 0);
            int r = round;
            RunWithWatchdog(() =>
            {
                unsafe
                {
                    using SesaFile file = SesaFile.FromMemory(guarded.Data, guarded.Length, options);
                    using EsaSuffixTree tree = EsaSuffixTree.FromFile(file);
                    Battery(tree, text, $"overflow round {r}");
                }
            }, $"overflow round {round}");
        }
    }

    [Test]
    public void ACleanImage_InGuardedMemory_AnswersLikeTheMemoryTree_AtTheEdgesOfTheMapping()
    {
        if (!OperatingSystem.IsLinux()) Assert.Ignore("the guard pages need Linux mmap / mprotect");
        string text = StringFamilies.Random(500, 4, 11);
        LcpChildCase c = LcpChildCase.Of(text);
        byte[] image = SesaTestFile.Image(c);
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        foreach (bool flush in new[] { true, false })
        {
            using var guarded = new Guarded(image, flush);
            unsafe
            {
                using SesaFile file = SesaFile.FromMemory(guarded.Data, guarded.Length, new EsaOpenOptions { Validation = EsaValidationLevel.Structure });
                using EsaSuffixTree tree = EsaSuffixTree.FromFile(file);
                for (int i = 0; i < 300; i++)
                {
                    string p = text.Substring(i, 1 + (i % 12));
                    Assert.That(tree.CountOccurrences(p), Is.EqualTo(memory.CountOccurrences(p)));
                }
                Assert.That(tree.GetAllSuffixes(), Is.EqualTo(memory.GetAllSuffixes()));
                Assert.That(tree.FindMaximalExactMatches(text[100..300], 6), Is.EqualTo(memory.FindMaximalExactMatches(text[100..300], 6)));
            }
        }
    }

    [Test]
    public void AChildTableThatLoops_InCheckedMode_IsBounded()
    {
        // every l-index chain of the root points back at itself: the walk must end with InvalidData, not loop
        string text = "abcabcabcabcabcabcabc";
        LcpChildCase c = LcpChildCase.Of(text);
        byte[] image = SesaTestFile.Image(c);
        SesaOracle o = SesaOracle.FromBytes(image);
        OracleSection child = o.Get(SesaFormat.KindChild);
        for (int cell = 0; cell <= c.N; cell++)
            BinaryPrimitives.WriteUInt32LittleEndian(image.AsSpan((int)child.Offset + (4 * cell)), (uint)Math.Max(0, cell - 1));
        SesaTestFile.Reseal(image);
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image);
        RunWithWatchdog(() => Battery(tree, text, "looping child table"), "looping child table");
    }

    [Test]
    public void ALyingLcp_WithATooSmallChildLcp_DoesNotLoopThePatternDescent()
    {
        string text = "abababababababababab";
        LcpChildCase c = LcpChildCase.Of(text);
        byte[] image = SesaTestFile.Image(c);
        SesaOracle o = SesaOracle.FromBytes(image);
        OracleSection lcp = o.Get(SesaFormat.KindLcp);
        for (int r = 2; r <= c.N; r++) image[lcp.Offset + r] = 1;
        SesaTestFile.Reseal(image);
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image);
        RunWithWatchdog(() => Battery(tree, text, "constant LCP"), "constant LCP");
    }

    [Test]
    public void AnEngineCreatedFromAFileThatWasDisposedUnderIt_IsRefusedAtTheBoundary()
    {
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, "abracadabra");
        SesaFile file = SesaFile.Open(path, new EsaOpenOptions());
        EsaSuffixTree tree = EsaSuffixTree.FromFile(file);
        tree.Dispose();
        Assert.That(() => tree.CountOccurrences("a"), Throws.TypeOf<ObjectDisposedException>());
    }
}
