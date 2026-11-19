using System.Runtime.InteropServices;
using System.Security.Cryptography;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The oracle and the plumbing of the PH6 tests: the naive BWT and Psi, a suffix array written at an offset of a container-like file, and one run of the production stages on it.</summary>
internal static class Ph6Run
{
    internal sealed record Outcome(BwtResult Bwt, PsiStats Psi, byte[] BwtBytes, uint[] PsiEntries, byte[] Container, int[] Sa, string BwtPath, long SaOffset, long PsiOffset);

    /// <summary>The suffix array of dense symbols from the definition (the empty suffix first).</summary>
    internal static int[] Sa(byte[] s)
    {
        int n = s.Length;
        int[] sa = [.. Enumerable.Range(0, n + 1)];
        Array.Sort(sa, (a, b) => a == b ? 0 : s.AsSpan(a).SequenceCompareTo(s.AsSpan(b)));
        return sa;
    }

    /// <summary><c>BWT[j] = T[SA[j] - 1]</c>, the byte 0 at the row of the suffix 0.</summary>
    internal static byte[] Bwt(byte[] s, int[] sa) => [.. sa.Select(p => p == 0 ? (byte)0 : s[p - 1])];

    /// <summary>A container-like file: <paramref name="size"/> bytes of 0xA5, the suffix array at <paramref name="saOffset"/>.</summary>
    internal static string WriteContainer(TempDir dir, int[] sa, long saOffset, long size, string name = "container.bin")
    {
        string path = dir.File(name);
        byte[] bytes = new byte[size];
        Array.Fill(bytes, (byte)0xA5);
        MemoryMarshal.AsBytes(sa.AsSpan()).CopyTo(bytes.AsSpan((int)saOffset));
        File.WriteAllBytes(path, bytes);
        return path;
    }

    /// <summary>The container as an open <see cref="MappedOutputFile"/> (created sparse, filled with 0xA5, the suffix array at <paramref name="saOffset"/>).</summary>
    internal static MappedOutputFile CreateContainer(string path, int[] sa, long saOffset, long size)
    {
        MappedOutputFile file = MappedOutputFile.Create(path, size, new MappedFileOptions { Reservation = ReservationMode.Sparse });
        using (MappedWindowWriter w = file.CreateWriter(512)) w.Append(new byte[size].Select(_ => (byte)0xA5).ToArray());
        using (MappedWindowWriter w = file.CreateWriter(saOffset, 4L * sa.Length, 512)) w.Append(MemoryMarshal.AsBytes(sa.AsSpan()));
        return file;
    }

    /// <summary>Runs the BWT pass and the Psi stream on <paramref name="s"/> with the suffix array inside a larger file; the container keeps 0xA5 outside the Psi range.</summary>
    internal static Outcome Run(TempDir dir, byte[] s, int sigma, BwtOptions? bwt = null, PsiOptions? psi = null, long saOffset = 64, long? psiOffsetOrNull = null, int[]? sa = null, long[]? counts = null, string tag = "")
    {
        int n = s.Length;
        long psiOffset = psiOffsetOrNull ?? (((saOffset + (4L * (n + 1)) + 4 + 63) / 64) * 64);
        int bits = TextOracle.BitsOf(Math.Max(1, sigma));
        sa ??= Sa(s);
        string textPath = TextOracle.WriteFile(dir, $"text{tag}.pk", s, bits);
        long size = Math.Max(saOffset + (4L * (n + 1)), psiOffset + (4L * (n + 1))) + 100;
        // the PSI range and the SA range must not overlap
        if (psiOffset < saOffset + (4L * (n + 1)) && saOffset < psiOffset + (4L * (n + 1))) throw new ArgumentException("ranges overlap");
        string container = dir.File($"container{tag}.bin");
        string bwtPath = dir.File($"bwt{tag}.bin");
        using MappedTextReader reader = MappedTextReader.Open(textPath, n, bits);
        BwtResult b;
        byte[] bwtBytes;
        PsiStats stats;
        using (MappedOutputFile file = CreateContainer(container, sa, saOffset, size))
        {
            file.Complete(false);
            b = BwtPass.Build(container, saOffset, reader, sigma, bwtPath, bwt);
            bwtBytes = File.ReadAllBytes(bwtPath);
            counts ??= b.Counts;
            stats = PsiStream.Build(bwtPath, n, b.PrimaryRow, counts, file, psiOffset, psi);
            file.Complete(false);
        }
        Assert.That(ProcMaps.Contains(container) ?? false, Is.False, "the container is not mapped any more: every window and writer was released");
        Assert.That(ProcMaps.Contains(bwtPath) ?? false, Is.False, "the BWT file is not mapped any more");
        byte[] all = File.ReadAllBytes(container);
        uint[] entries = MemoryMarshal.Cast<byte, uint>(all.AsSpan((int)psiOffset, 4 * (n + 1))).ToArray();
        return new Outcome(b, stats, bwtBytes, entries, all, sa, bwtPath, saOffset, psiOffset);
    }

    /// <summary>Compares everything with the definitions.</summary>
    internal static void Compare(byte[] s, int sigma, Outcome o, string label)
    {
        long saOffset = o.SaOffset, psiOffset = o.PsiOffset;
        int n = s.Length;
        Assert.That(o.Bwt.N, Is.EqualTo(n), $"{label}: n");
        Assert.That(o.BwtBytes, Is.EqualTo(Bwt(s, o.Sa)), $"{label}: BWT");
        Assert.That(o.Bwt.PrimaryRow, Is.EqualTo(Array.IndexOf(o.Sa, 0)), $"{label}: the row of the suffix 0");
        var counts = new long[sigma];
        foreach (byte x in s) counts[x]++;
        Assert.That(o.Bwt.Counts, Is.EqualTo(counts), $"{label}: counts");
        Assert.That(o.PsiEntries, Is.EqualTo(SesaRun.NaivePsi(o.Sa)), $"{label}: PSI");
        for (int i = 0; i < o.Container.Length; i++)
        {
            bool inSa = i >= saOffset && i < saOffset + (4L * (n + 1));
            bool inPsi = i >= psiOffset && i < psiOffset + (4L * (n + 1));
            if (!inPsi && !inSa && o.Container[i] != 0xA5) Assert.Fail($"{label}: byte {i} outside the ranges was changed");
        }
    }
}

/// <summary>PH6 of ADR-004 section 3.3: the BWT pass and the Psi distribution against the definitions (gate G1 of work package 1b.7e).</summary>
[TestFixture]
public sealed class Ph6Tests
{
    private static readonly BwtOptions TinyBwt = new() { ReadBytes = 8, WindowBytes = 4 };

    private static PsiOptions TinyPsi(int v) => new() { WindowBytes = 4 * (1 + (v % 3)), SymbolsPerPass = 1 + (v % 4), ReadBytes = 1 + (v % 5), ArenaGuards = ArenaOptions.Guarded(0xCD) };

    private static BwtOptions TinyBwtOptions(int v) => new() { ReadBytes = 4 * (1 + (v % 4)), WindowBytes = 1 + (v % 7), ArenaGuards = ArenaOptions.Guarded(0xCD) };

    [Test]
    public void TheDefaults_AreWindowsOfAtLeastOneMiB_AndTheCapsOfTheOwner()
    {
        Assert.That(new BwtOptions().WindowBytes, Is.GreaterThanOrEqualTo(1 << 20), "the BWT writer");
        Assert.That(new PsiOptions().WindowBytes, Is.GreaterThanOrEqualTo(1 << 20), "the Psi cursors (E9: the view accessor flushes every window)");
        Assert.That(new SesaPh6Tuning().WindowBytes, Is.Null);
        Assert.That(new PsiOptions().SymbolsPerPass, Is.EqualTo(64));
        Assert.That(new PsiOptions().WindowCapBytes, Is.EqualTo(160L << 20));
        Assert.That(new BwtOptions().ReadBytes, Is.EqualTo(1 << 20));
        Assert.That(new PsiOptions().ReadBytes, Is.EqualTo(1 << 20));
    }

    [Test]
    public void TheDefaultWindows_OfARealRun_AreNotSmallerThanOneMiB()
    {
        using var dir = new TempDir();
        byte[] s = TextOracle.Random(3_000_000, 4, 5);
        Ph6Run.Outcome o;
        var counters = new IoCounters();
        // a real-size run needs a real-size suffix array: the reversed identity is a valid permutation for the plumbing (the stage does not check the order)
        int[] id = [.. Enumerable.Range(0, s.Length + 1).Select(i => s.Length - i)];
        o = Ph6Run.Run(dir, s, 4, new BwtOptions { Counters = counters }, new PsiOptions { Counters = counters }, sa: id);
        Assert.That(counters.MappedMaxWindowBytes, Is.EqualTo(1 << 20), "the largest window mapped by the BWT writer and the four Psi cursors");
        Assert.That(o.Psi.Windows, Is.GreaterThanOrEqualTo(4));
    }

    [Test]
    public void EveryStringOfSmallAlphabets_GivesTheDefinitionsBwtAndPsi()
    {
        using var dir = new TempDir();
        int v = 0;
        foreach ((int sigma, int maxLength) in new[] { (1, 8), (2, 8), (3, 5), (4, 4) })
        {
            foreach (byte[] s in TextOracle.AllTexts(sigma, maxLength))
            {
                v++;
                long saOffset = new long[] { 0, 64, 12, 4096 }[v % 4];
                long psiOffset = saOffset + (4L * (s.Length + 1)) + new long[] { 0, 8, 100 }[v % 3];
                Ph6Run.Outcome o = Ph6Run.Run(dir, s, sigma, TinyBwtOptions(v), TinyPsi(v), saOffset, psiOffset, tag: "x");
                Ph6Run.Compare(s, sigma, o, $"sigma {sigma}, text {string.Join(',', s)}, variant {v}");
            }
        }
        Assert.That(v, Is.GreaterThan(900));
    }

    [Test]
    public void TheFamilies_GiveTheDefinitionsBwtAndPsi_WithDefaultAndWithTinyOptions()
    {
        using var dir = new TempDir();
        int v = 0;
        foreach (int n in new[] { 600, 3000 })
        {
            foreach ((string name, string text) in LcpChildRun.Families(n))
            {
                v++;
                char[] letters = [.. text.Distinct().OrderBy(c => c)];
                byte[] s = [.. text.Select(c => (byte)Array.BinarySearch(letters, c))];
                int sigma = letters.Length;
                Ph6Run.Outcome o = Ph6Run.Run(dir, s, sigma, tag: "d");
                Ph6Run.Compare(s, sigma, o, $"{name} {n} defaults");
                o = Ph6Run.Run(dir, s, sigma, TinyBwtOptions(v), TinyPsi(v), 4096, 4096 + (4L * (s.Length + 1)) + 64, tag: "t");
                Ph6Run.Compare(s, sigma, o, $"{name} {n} tiny {v}");
            }
        }
    }

    [Test]
    public void ARandomTextOver200Symbols_TakesFourPassesAtSixtyFour_AndGivesTheSameBytesAtEveryPassWidth()
    {
        using var dir = new TempDir();
        byte[] s = TextOracle.Random(6000, 200, 9);
        Assert.That(s.Distinct().Count(), Is.GreaterThan(190));
        int sigma = s.Max() + 1;
        string reference = "";
        foreach (int perPass in new[] { 64, 63, 200, 1, 7, 128 })
        {
            Ph6Run.Outcome o = Ph6Run.Run(dir, s, sigma, psi: new PsiOptions { SymbolsPerPass = perPass, WindowBytes = 40 }, tag: "p");
            Ph6Run.Compare(s, sigma, o, $"{perPass} symbols per pass");
            Assert.That(o.Psi.Passes, Is.EqualTo((sigma + perPass - 1) / perPass), $"passes at {perPass}");
            Assert.That(o.Psi.SymbolsPerPass, Is.EqualTo(Math.Min(perPass, sigma)));
            Assert.That(o.Psi.BwtBytesRead, Is.EqualTo(o.Psi.Passes * (long)(s.Length + 1)), "one read of the BWT per pass");
            string sha = Convert.ToHexString(SHA256.HashData(MemoryMarshal.AsBytes(o.PsiEntries.AsSpan())));
            if (reference.Length == 0) reference = sha;
            Assert.That(sha, Is.EqualTo(reference), $"the bytes do not depend on the passes ({perPass})");
        }
    }

    [Test]
    public void ABytesAlphabetOf256Symbols_IsHandled()
    {
        using var dir = new TempDir();
        byte[] s = [.. Enumerable.Range(0, 256).Select(i => (byte)i), .. TextOracle.Random(2000, 256, 3)];
        Ph6Run.Outcome o = Ph6Run.Run(dir, s, 256);
        Ph6Run.Compare(s, 256, o, "256 symbols");
        Assert.That(o.Psi.Passes, Is.EqualTo(4));
    }

    [Test]
    public void PsiComposedNPlusOneTimesFromThePrimaryRow_VisitsEveryRowOnce()
    {
        using var dir = new TempDir();
        foreach ((string name, string text) in LcpChildRun.Families(2000))
        {
            char[] letters = [.. text.Distinct().OrderBy(c => c)];
            byte[] s = [.. text.Select(c => (byte)Array.BinarySearch(letters, c))];
            Ph6Run.Outcome o = Ph6Run.Run(dir, s, letters.Length, tag: "c");
            int n = s.Length;
            var seen = new bool[n + 1];
            long row = o.Bwt.PrimaryRow;
            for (int step = 0; step <= n; step++)
            {
                Assert.That(seen[row], Is.False, $"{name}: row {row} visited twice at step {step}");
                seen[row] = true;
                row = o.PsiEntries[row];
            }
            Assert.That(seen.All(x => x), Is.True, $"{name}: not every row visited");
            Assert.That(row, Is.Zero, $"{name}: the suffix n is row 0, and PSI[0] = 0 by convention");
            Assert.That(o.PsiEntries[0], Is.Zero, $"{name}: PSI[0]");
        }
    }

    [Test]
    public void TheEmptyText_AndTheTextsOfOneAndTwoSymbols()
    {
        using var dir = new TempDir();
        foreach (byte[] s in new[] { Array.Empty<byte>(), [0], [0, 0], [0, 1], [1, 0], [2, 0, 1] })
        {
            int sigma = s.Length == 0 ? 0 : s.Max() + 1;
            Ph6Run.Outcome o = Ph6Run.Run(dir, s, sigma, tag: "e");
            Ph6Run.Compare(s, sigma, o, $"text {string.Join(',', s)}");
        }
    }

    [Test]
    public void ThePrimaryRow_CanBeRow1_AnInnerRow_AndTheLastRow()
    {
        using var dir = new TempDir();
        var rows = new HashSet<long>();
        foreach (byte[] s in TextOracle.AllTexts(2, 7))
        {
            if (s.Length < 3) continue;
            Ph6Run.Outcome o = Ph6Run.Run(dir, s, 2, tag: "r");
            rows.Add(o.Bwt.PrimaryRow == 1 ? 1 : o.Bwt.PrimaryRow == s.Length ? 2 : 3);
        }
        Assert.That(rows, Is.EquivalentTo(new long[] { 1, 2, 3 }));
    }

    [Test]
    public void TheCountedIo_IsOneReadOfTheSuffixArray_OneWriteOfTheBwt_OneReadOfTheBwtPerPass_AndOneWriteOfPsi()
    {
        using var dir = new TempDir();
        byte[] s = TextOracle.Random(5000, 4, 21);
        int n = s.Length;
        var bwtCounters = new IoCounters();
        var psiCounters = new IoCounters();
        Ph6Run.Outcome o = Ph6Run.Run(dir, s, 4, new BwtOptions { Counters = bwtCounters, ReadBytes = 1000, WindowBytes = 512 }, new PsiOptions { Counters = psiCounters, WindowBytes = 256 });
        Assert.That(bwtCounters.ReadBytes, Is.EqualTo(4L * (n + 1)), "the suffix array is read once");
        Assert.That(bwtCounters.MappedWriteBytes, Is.EqualTo(n + 1), "the BWT is written once, one byte per row");
        Assert.That(bwtCounters.MappedMaxWindowBytes, Is.EqualTo(512));
        Assert.That(psiCounters.ReadBytes, Is.EqualTo(n + 1L), "one pass reads the BWT once");
        Assert.That(psiCounters.MappedWriteBytes, Is.EqualTo(4L * (n + 1)), "PSI is written once");
        Assert.That(psiCounters.MappedMaxWindowBytes, Is.LessThanOrEqualTo(256));
        Assert.That((o.Psi.Passes, o.Psi.PsiBytesWritten, o.Psi.BwtBytesRead), Is.EqualTo((1, 4L * (n + 1), n + 1L)));
        Assert.That(o.Psi.MappedBoundBytes, Is.EqualTo(4 * 64L).Or.EqualTo(4 * 256L), "sigma windows of 256 bytes at most");
    }

    [Test]
    public void ABwtFileThatIsNotWritten_IsRemovedWhenTheSuffixArrayIsDamaged()
    {
        using var dir = new TempDir();
        byte[] s = [0, 1, 0, 1, 1, 0];
        int[] good = Ph6Run.Sa(s);
        int n = s.Length;
        (string Label, Action<int[]> Damage, bool Primary)[] damages =
        [
            ("SA[0] is not n", sa => sa[0] = 1, false),
            ("an entry above n", sa => sa[3] = n + 1, false),
            ("the suffix 0 twice", sa => sa[Array.IndexOf(sa, 0) == 2 ? 3 : 2] = 0, false),
            ("no suffix 0", sa => sa[Array.IndexOf(sa, 0)] = 3, false),
        ];
        foreach ((string label, Action<int[]> damage, _) in damages)
        {
            int[] sa = (int[])good.Clone();
            damage(sa);
            string textPath = TextOracle.WriteFile(dir, "dmg.pk", s, 1);
            string container = Ph6Run.WriteContainer(dir, sa, 64, 64 + (4L * (n + 1)) + 100, "dmg.bin");
            string bwt = dir.File("dmg.bwt");
            using MappedTextReader reader = MappedTextReader.Open(textPath, n, 1);
            Assert.Throws<InvalidDataException>(() => BwtPass.Build(container, 64, reader, 2, bwt), label);
            Assert.That(File.Exists(bwt), Is.False, $"{label}: the BWT file is left behind");
            ProcFds.AssertClosed(label, [container, bwt], []);
            Assert.That(ProcMaps.Contains(bwt) ?? false, Is.False, $"{label}: the BWT file is still mapped");
        }
    }

    [Test]
    public void ASuffixArrayRangeThatDoesNotFitTheFile_AndASymbolOutsideTheAlphabet_AreRefused()
    {
        using var dir = new TempDir();
        byte[] s = [0, 1, 2, 1, 0, 2, 2];
        int n = s.Length;
        int[] sa = Ph6Run.Sa(s);
        string textPath = TextOracle.WriteFile(dir, "t.pk", s, 2);
        string container = Ph6Run.WriteContainer(dir, sa, 64, 64 + (4L * (n + 1)), "c.bin");
        string bwt = dir.File("o.bwt");
        using MappedTextReader reader = MappedTextReader.Open(textPath, n, 2);
        Assert.Throws<InvalidDataException>(() => BwtPass.Build(container, 68, reader, 3, bwt), "the range ends 4 bytes after the file");
        Assert.Throws<InvalidDataException>(() => BwtPass.Build(container, 64, reader, 2, bwt), "the symbol 2 is outside an alphabet of 2");
        Assert.Throws<ArgumentOutOfRangeException>(() => BwtPass.Build(container, -4, reader, 3, bwt));
        Assert.Throws<ArgumentOutOfRangeException>(() => BwtPass.Build(container, 64, reader, 257, bwt));
        Assert.Throws<ArgumentOutOfRangeException>(() => BwtPass.Build(container, 64, reader, -1, bwt));
        Assert.Throws<ArgumentOutOfRangeException>(() => BwtPass.Build(container, 64, reader, 3, bwt, new BwtOptions { ReadBytes = 3 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => BwtPass.Build(container, 64, reader, 3, bwt, new BwtOptions { WindowBytes = 0 }));
        Assert.That(File.Exists(bwt), Is.False);
    }

    [Test]
    public void AFullDiskAtTheCreationOfTheBwt_IsAnIoException_AndLeavesNothing()
    {
        using var dir = new TempDir();
        byte[] s = TextOracle.Random(300, 4, 5);
        int[] sa = Ph6Run.Sa(s);
        string textPath = TextOracle.WriteFile(dir, "t.pk", s, 2);
        string container = Ph6Run.WriteContainer(dir, sa, 64, 64 + (4L * 301) + 100, "c.bin");
        string bwt = dir.File("o.bwt");
        using MappedTextReader reader = MappedTextReader.Open(textPath, s.Length, 2);
        var full = new FakeReservationSource { CanPreallocate = true, PreallocateThrows = new IOException("disk full") };
        var ex = Assert.Throws<IOException>(() => BwtPass.Build(container, 64, reader, 4, bwt, new BwtOptions { Files = new MappedFileOptions { Source = full } }));
        Assert.That(ex!.Message, Does.Contain("disk full"));
        Assert.That(File.Exists(bwt), Is.False);
        Assert.That(full.PreallocateCalls, Is.EqualTo(1));
        Assert.Throws<DirectoryNotFoundException>(() => BwtPass.Build(container, 64, reader, 4, Path.Combine(dir.Path, "missing", "o.bwt")));
    }

    [Test]
    public void ABwtThatDoesNotMatchTheCounts_ALengthThatIsWrong_AndBadArguments_AreRefusedByThePsiStream()
    {
        using var dir = new TempDir();
        byte[] s = [0, 1, 1, 0, 2, 2, 1];
        int n = s.Length;
        Ph6Run.Outcome o = Ph6Run.Run(dir, s, 3, tag: "ok");
        long size = o.Container.Length;
        long[] counts = [.. o.Bwt.Counts];
        using MappedOutputFile file = Ph6Run.CreateContainer(dir.File("container-ok2.bin"), o.Sa, 64, size);
        Assert.Throws<InvalidDataException>(() => PsiStream.Build(o.BwtPath, n, o.Bwt.PrimaryRow, [counts[0] + 1, counts[1] - 1, counts[2]], file, 192), "a symbol occurs more often than counted");
        Assert.Throws<InvalidDataException>(() => PsiStream.Build(o.BwtPath, n, o.Bwt.PrimaryRow, [counts[0] - 1, counts[1] + 1, counts[2]], file, 192), "a symbol occurs less often than counted");
        Assert.Throws<ArgumentException>(() => PsiStream.Build(o.BwtPath, n, o.Bwt.PrimaryRow, [counts[0] + 1, counts[1], counts[2]], file, 192), "the counts sum to more than n");
        Assert.Throws<ArgumentException>(() => PsiStream.Build(o.BwtPath, n, o.Bwt.PrimaryRow, [counts[0] - 1, counts[1], counts[2]], file, 192), "the counts sum to less than n");
        Assert.Throws<ArgumentException>(() => PsiStream.Build(o.BwtPath, n, o.Bwt.PrimaryRow, [-1, counts[1] + 1, counts[2] + counts[0]], file, 192), "a negative count");
        Assert.Throws<ArgumentOutOfRangeException>(() => PsiStream.Build(o.BwtPath, n, n + 1, counts, file, 192), "the primary row is outside the rows");
        Assert.Throws<ArgumentOutOfRangeException>(() => PsiStream.Build(o.BwtPath, n, -1, counts, file, 192));
        Assert.Throws<ArgumentOutOfRangeException>(() => PsiStream.Build(o.BwtPath, n, o.Bwt.PrimaryRow, counts, file, size - 10), "the range is outside the file");
        Assert.Throws<InvalidDataException>(() => PsiStream.Build(o.BwtPath, n + 1, o.Bwt.PrimaryRow, [.. counts.Take(2), counts[2] + 1], file, 192), "the BWT has n + 1 rows, not n + 2");
        Assert.Throws<ArgumentOutOfRangeException>(() => PsiStream.Build(o.BwtPath, n, o.Bwt.PrimaryRow, counts, file, 192, new PsiOptions { WindowBytes = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => PsiStream.Build(o.BwtPath, n, o.Bwt.PrimaryRow, counts, file, 192, new PsiOptions { SymbolsPerPass = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => PsiStream.Build(o.BwtPath, n, o.Bwt.PrimaryRow, counts, file, 192, new PsiOptions { ReadBytes = 0 }));
        Assert.Throws<ArgumentException>(() => PsiStream.Build(o.BwtPath, n, o.Bwt.PrimaryRow, [.. Enumerable.Repeat(0L, 257).Select((_, i) => i == 0 ? (long)n : 0L)], file, 192), "a BWT is one byte per row: at most 256 symbols");
    }

    [Test]
    public void TheWindowsOfThePsiStream_AreChargedToTheCap_AndGivenBack()
    {
        using var dir = new TempDir();
        byte[] s = TextOracle.Random(500, 4, 5);
        var budget = new ArenaBudget(long.MaxValue);
        Ph6Run.Outcome o = Ph6Run.Run(dir, s, 4, psi: new PsiOptions { Budget = budget, WindowBytes = 64 }, tag: "b");
        Assert.That(budget.LiveBytes, Is.Zero, "everything is given back");
        Assert.That(budget.PeakBytes, Is.GreaterThanOrEqualTo(o.Psi.MappedBoundBytes), "the windows are on the account");
        Assert.That(o.Psi.MappedBoundBytes, Is.EqualTo(4 * 64L));
        string container = dir.File("container-b2.bin");
        using MappedOutputFile file = Ph6Run.CreateContainer(container, o.Sa, 64, o.Container.Length);
        var small = new ArenaBudget(100);
        Assert.Throws<ArenaBudgetExceededException>(() => PsiStream.Build(o.BwtPath, s.Length, o.Bwt.PrimaryRow, o.Bwt.Counts, file, 192, new PsiOptions { Budget = small, WindowBytes = 64 }));
        Assert.That(small.LiveBytes, Is.Zero);
        Assert.Throws<ArenaBudgetExceededException>(() => PsiStream.Build(o.BwtPath, s.Length, o.Bwt.PrimaryRow, o.Bwt.Counts, file, 192, new PsiOptions { WindowBytes = 64, WindowCapBytes = 100 }), "the cap of the windows itself");
        var bwtSmall = new ArenaBudget(100);
        string textPath = TextOracle.WriteFile(dir, "tb.pk", s, 2);
        using MappedTextReader reader = MappedTextReader.Open(textPath, s.Length, 2);
        Assert.Throws<ArenaBudgetExceededException>(() => BwtPass.Build(container, 64, reader, 4, dir.File("b2.bwt"), new BwtOptions { Budget = bwtSmall }));
        Assert.That(bwtSmall.LiveBytes, Is.Zero);
        Assert.That(File.Exists(dir.File("b2.bwt")), Is.False);
    }

    [Test]
    public void APsiRangeIsFullyWritten_EvenWhereTheFileHeldSomethingElse()
    {
        using var dir = new TempDir();
        byte[] s = TextOracle.Random(100, 3, 5);
        Ph6Run.Outcome o = Ph6Run.Run(dir, s, 3, tag: "w");
        Assert.That(o.PsiEntries, Is.EqualTo(SesaRun.NaivePsi(o.Sa)));
        Assert.That(o.PsiEntries.Contains(0xA5A5A5A5u), Is.False, "no entry of the range kept the 0xA5 fill of the file");
    }
}
