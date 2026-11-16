using System.IO.MemoryMappedFiles;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph6;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.Ph6;

/// <summary>
/// Gate G1 of E9 (ADR-004 section 5, row E9): PSI by sigma-way distribution of the BWT is byte-identical to <c>ISA[SA[i] + 1]</c> of the in-RAM tables of the prototype P1 (<c>LinkTables</c>), for alphabets from 1 to 256 symbols
/// (single pass and several passes), windows from one entry up, the row of the suffix 0 anywhere, both conventions of entry 0, and it refuses a BWT that does not match its counts.
/// </summary>
[TestFixture]
public class PsiStreamTests
{
    private static readonly int[] Sigmas = [1, 2, 4, 5, 20, 64, 65, 255, 256];
    private static readonly int[] Windows = [1, 4, 6, 12, 64, 4096, 65536];
    private static readonly int[] Groups = [1, 3, 7, 64, 256];

    /// <summary>The BWT and its meta from the definition (not from <see cref="BwtPass"/>), so that this class tests the Psi pass alone.</summary>
    internal static (string Bwt, BwtMeta Meta) WriteBwt(Ph6Dir dir, Ph6Case c, string name = "x.bwt")
    {
        string path = dir.File(name);
        File.WriteAllBytes(path, c.ExpectedBwt());
        var meta = new BwtMeta(c.N, c.Primary, c.Counts());
        File.WriteAllText(BwtMeta.MetaPath(path), meta.Text());
        return (path, meta);
    }

    private static void CheckIdentity(Ph6Case c, int window, int group, bool cyclic, string label)
    {
        using var dir = new Ph6Dir();
        (string bwt, BwtMeta meta) = WriteBwt(dir, c);
        string psi = dir.File("x.psi");
        PsiStats st = PsiStream.Build(bwt, meta, psi, new PsiOptions { WindowBytes = window, SymbolsPerPass = group, CyclicEntryZero = cyclic, ReadBytes = 1 + (window % 7) });
        Assert.That(File.ReadAllBytes(psi), Is.EqualTo(c.ExpectedPsiBytes(cyclic)), $"{label}: window {window} group {group} cyclic {cyclic}");
        Assert.That(st.Passes, Is.EqualTo((c.Sigma + group - 1) / group), $"{label}: passes");
    }

    [Test]
    public void RandomTexts_AllAlphabets_RotatingWindowsGroupsAndConventions()
    {
        int v = 0;
        foreach (int sigma in Sigmas)
            for (int rep = 0; rep < 12; rep++, v++)
            {
                int n = sigma + 1 + ((rep * 41) % 240);
                Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(n, sigma, (ulong)(1000 + v)));
                Assert.That(c.Sigma, Is.EqualTo(sigma));
                CheckIdentity(c, Windows[v % Windows.Length], Groups[(v / 2) % Groups.Length], cyclic: v % 3 == 2, $"sigma {sigma} n {n}");
            }
    }

    [Test]
    public void OneText_EveryWindowGroupAndConvention()
    {
        foreach (int sigma in new[] { 5, 65, 256 })
        {
            Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(700, sigma, 77));
            foreach (int w in Windows)
                foreach (int g in Groups)
                    foreach (bool cyclic in new[] { false, true })
                        CheckIdentity(c, w, g, cyclic, $"sigma {sigma}");
        }
    }

    [Test]
    public void DefaultOptions_MultiPassAlphabets()
    {
        foreach (int sigma in new[] { 64, 65, 129, 255, 256 })
        {
            Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(3000, sigma, (ulong)sigma));
            using var dir = new Ph6Dir();
            (string bwt, BwtMeta meta) = WriteBwt(dir, c);
            string psi = dir.File("x.psi");
            PsiStats st = PsiStream.Build(bwt, meta, psi);
            Assert.That(File.ReadAllBytes(psi), Is.EqualTo(c.ExpectedPsiBytes(false)), $"sigma {sigma}");
            Assert.That(st.Passes, Is.EqualTo(sigma <= 64 ? 1 : (sigma + 63) / 64), $"sigma {sigma}: passes");
            Assert.That(st.SymbolsPerPass, Is.EqualTo(Math.Min(64, sigma)), $"sigma {sigma}: open cursors");
        }
    }

    [Test]
    public void TheRowOfTheSuffixZero_AtTheFirstMiddleAndLastRow()
    {
        // suffix 0 is the smallest non-empty suffix (row 1), a middle one, the largest (row n)
        string[] texts = ["abcdefg", "dcbaxyz", "gfedcba", "bab", "ba", "ab", "aab", "abab" + "b"];
        bool first = false, inner = false, last = false;
        foreach (string text in texts)
        {
            Ph6Case c = Ph6Case.OfNaive(text);
            first |= c.Primary == 1;
            last |= c.Primary == c.N;
            inner |= c.Primary > 1 && c.Primary < c.N;
            foreach (bool cyclic in new[] { false, true })
                foreach (int w in new[] { 1, 64 })
                    CheckIdentity(c, w, 64, cyclic, $"'{text}' primary row {c.Primary}");
        }
        Assert.That(first && inner && last, Is.True, "the cases cover row 1, an inner row and the last row");
        Assert.That(Ph6Case.OfNaive("abcdefg").Primary, Is.EqualTo(1));
        Assert.That(Ph6Case.OfNaive("gfedcba").Primary, Is.EqualTo(7));
    }

    [Test]
    public void EntryZero_IsZero_OrTheRowOfTheSuffixZero_AndNothingElseChanges()
    {
        Ph6Case c = Ph6Case.OfNaive("mississippi");
        using var dir = new Ph6Dir();
        (string bwt, BwtMeta meta) = WriteBwt(dir, c);
        string zero = dir.File("zero.psi"), cyc = dir.File("cyc.psi");
        PsiStream.Build(bwt, meta, zero, new PsiOptions { CyclicEntryZero = false });
        PsiStream.Build(bwt, meta, cyc, new PsiOptions { CyclicEntryZero = true });
        int[] z = PsiFiles.ReadInts(zero), y = PsiFiles.ReadInts(cyc);
        Assert.That(z[0], Is.EqualTo(0));
        Assert.That(y[0], Is.EqualTo(c.Primary));
        Assert.That(z[1..], Is.EqualTo(y[1..]));
        Assert.That(z[1..], Is.EqualTo(c.PsiRam[1..]));
        // the cyclic Psi is one cycle through all n + 1 rows
        int at = 0, steps = 0;
        do
        {
            at = y[at];
            steps++;
        }
        while (at != 0);
        Assert.That(steps, Is.EqualTo(c.N + 1));
    }

    [Test]
    public void ExhaustiveTinyTexts()
    {
        int v = 0;
        foreach ((int sigma, int len) in new[] { (2, 9), (3, 6), (4, 4) })
            foreach (string s in ExhaustiveStrings.All(sigma, len))
            {
                if (s.Length == 0) continue;
                Ph6Case c = Ph6Case.OfNaive(s);
                CheckIdentity(c, Windows[v % Windows.Length], Groups[v % Groups.Length], cyclic: v % 2 == 1, $"'{s}'");
                v++;
            }
    }

    [Test]
    public void AdversarialFamilies()
    {
        foreach (int n in new[] { 1, 2, 3, 5, 8, 17, 64, 100, 257, 400 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 11))
            {
                if (text.Distinct().Count() > 256) continue; // BWT symbols are bytes
                Ph6Case c = Ph6Case.Of(text);
                CheckIdentity(c, 12, 5, cyclic: n % 2 == 1, $"{name} n {n}");
                CheckIdentity(c, 65536, 64, cyclic: false, $"{name} n {n}");
            }
    }

    [Test]
    public void LongPeriodicFamilies_AllSameFibonacciThueMorse()
    {
        foreach (string text in new[] { StringFamilies.AllSame(30000), StringFamilies.Fibonacci(30000), StringFamilies.ThueMorse(30000) })
        {
            Ph6Case c = Ph6Case.Of(text);
            CheckIdentity(c, 4096, 64, false, "long family");
            CheckIdentity(c, 8, 1, true, "long family");
        }
    }

    [Test]
    public void EmptyAndSingleSymbolTexts()
    {
        // n = 0: one row (the empty suffix), no symbol; n = 1: two rows
        using (var dir = new Ph6Dir())
        {
            string bwt = dir.File("e.bwt");
            File.WriteAllBytes(bwt, [0]);
            var meta = new BwtMeta(0, 0, []);
            string psi = dir.File("e.psi");
            PsiStats st = PsiStream.Build(bwt, meta, psi);
            Assert.That(File.ReadAllBytes(psi), Is.EqualTo(new byte[4]));
            Assert.That(st.Passes, Is.EqualTo(0));
        }
        CheckIdentity(Ph6Case.OfNaive("a"), 4, 64, false, "'a'");
        CheckIdentity(Ph6Case.OfNaive("a"), 4, 64, true, "'a'");
    }

    [Test]
    public void ExistingLongerFile_IsReplaced()
    {
        Ph6Case c = Ph6Case.OfNaive("banana");
        using var dir = new Ph6Dir();
        (string bwt, BwtMeta meta) = WriteBwt(dir, c);
        string psi = dir.File("x.psi");
        File.WriteAllBytes(psi, Enumerable.Repeat((byte)0xAB, 1000).ToArray());
        PsiStream.Build(bwt, meta, psi);
        Assert.That(File.ReadAllBytes(psi), Is.EqualTo(c.ExpectedPsiBytes(false)));
    }

    [Test]
    public void Statistics_PassesWindowsAndBytes()
    {
        Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(1000, 20, 5));
        long[] counts = c.Counts();
        foreach (var (window, group) in new[] { (4, 64), (16, 7), (400, 20), (65536, 3) })
        {
            using var dir = new Ph6Dir();
            (string bwt, BwtMeta meta) = WriteBwt(dir, c);
            long r0 = IoCounters.ReadBytes, m0 = MappedIo.WriteBytes, w0 = MappedIo.Windows;
            PsiStats st = PsiStream.Build(bwt, meta, dir.File("x.psi"), new PsiOptions { WindowBytes = window, SymbolsPerPass = group, ReadBytes = 100 });
            int passes = (20 + group - 1) / group;
            int entries = Math.Max(1, window / 4);
            long expectedWindows = counts.Sum(k => (k + entries - 1) / entries);
            Assert.That(st.Passes, Is.EqualTo(passes), "passes");
            Assert.That(st.SymbolsPerPass, Is.EqualTo(Math.Min(group, 20)), "open cursors");
            Assert.That(st.BwtBytesRead, Is.EqualTo(passes * (c.N + 1L)), "BWT bytes read");
            Assert.That(IoCounters.ReadBytes - r0, Is.EqualTo(passes * (c.N + 1L)), "counted explicit reads");
            Assert.That(st.PsiBytesWritten, Is.EqualTo(4 * (c.N + 1L)), "PSI bytes written");
            Assert.That(MappedIo.WriteBytes - m0, Is.EqualTo(4 * (c.N + 1L)), "counted mapped writes (entry 0 included)");
            Assert.That(st.Windows, Is.EqualTo(expectedWindows), "windows of the cursors");
            Assert.That(MappedIo.Windows - w0, Is.EqualTo(expectedWindows + 1), "windows counted: the cursors' and the one of entry 0");
            Assert.That(st.MappedBoundBytes, Is.EqualTo(Math.Min(group, 20) * 4L * entries), "mapped working set: open cursors x window");
        }
    }

    [Test]
    public void BwtThatDoesNotMatchItsCounts_IsRefused()
    {
        Ph6Case c = Ph6Case.OfNaive("abracadabra");
        long[] counts = c.Counts();
        using var dir = new Ph6Dir();
        (string bwt, BwtMeta meta) = WriteBwt(dir, c);
        string psi = dir.File("x.psi");

        // one symbol counted less than it occurs (its cursor overflows) and another counted more (it would underflow)
        var fewer = (long[])counts.Clone();
        fewer[0]--;
        fewer[1]++;
        InvalidDataException over = Assert.Throws<InvalidDataException>(() => PsiStream.Build(bwt, meta with { Counts = fewer }, psi), "overflow of a cursor");
        Assert.That(over.Message, Does.Contain("more often"), "the overflow is refused when the cursor is asked for an entry beyond its range");
        var more = (long[])counts.Clone();
        more[0]++;
        more[1]--;
        Assert.Throws<InvalidDataException>(() => PsiStream.Build(bwt, meta with { Counts = more }, psi), "underflow of a cursor, the overflow of another symbol comes first");
        // one symbol per pass: the symbol that is short of its count is finished before the other one overflows, and it is the end-of-pass check that refuses it
        InvalidDataException under = Assert.Throws<InvalidDataException>(() => PsiStream.Build(bwt, meta with { Counts = more }, psi, new PsiOptions { SymbolsPerPass = 1 }), "underflow of a cursor");
        Assert.That(under.Message, Does.Contain("times in the BWT"));
        // windows of one entry: every window a cursor opens is filled completely, the cursor is short of its count by whole windows
        foreach (int window in new[] { 4, 8 })
        {
            // the symbol with the raised count ends its pass with every window full and the range not exhausted: the end-of-pass check names it (the symbol with the lowered count would overflow only in a later pass)
            InvalidDataException shortOfRange = Assert.Throws<InvalidDataException>(() => PsiStream.Build(bwt, meta with { Counts = more }, psi, new PsiOptions { WindowBytes = window, SymbolsPerPass = 1 }), $"underflow with windows of {window} bytes");
            Assert.That(shortOfRange.Message, Does.Contain("times in the BWT"), $"windows of {window} bytes");
            Assert.Throws<InvalidDataException>(() => PsiStream.Build(bwt, meta with { Counts = fewer }, psi, new PsiOptions { WindowBytes = window }), $"overflow with windows of {window} bytes");
        }
    }

    [Test]
    public void BwtSymbolOutsideTheAlphabet_IsRefused()
    {
        Ph6Case c = Ph6Case.OfNaive("abracadabra");
        using var dir = new Ph6Dir();
        (string bwt, BwtMeta meta) = WriteBwt(dir, c);
        byte[] bytes = File.ReadAllBytes(bwt);
        int row = Enumerable.Range(0, bytes.Length).First(j => j != c.Primary);
        bytes[row] = (byte)c.Sigma;
        File.WriteAllBytes(bwt, bytes);
        Assert.Throws<InvalidDataException>(() => PsiStream.Build(bwt, meta, dir.File("x.psi")));
    }

    [Test]
    public void InconsistentArguments_AreRefused()
    {
        Ph6Case c = Ph6Case.OfNaive("abracadabra");
        using var dir = new Ph6Dir();
        (string bwt, BwtMeta meta) = WriteBwt(dir, c);
        string psi = dir.File("x.psi");
        Assert.Throws<ArgumentException>(() => PsiStream.Build(bwt, meta with { N = meta.N + 1 }, psi), "counts sum to less than n");
        Assert.That(c.Primary, Is.LessThan(c.N - 1), "the primary row stays valid when n is lowered");
        Assert.Throws<ArgumentException>(() => PsiStream.Build(bwt, meta with { N = meta.N - 1 }, psi), "counts sum to more than n");
        Assert.Throws<ArgumentOutOfRangeException>(() => PsiStream.Build(bwt, meta with { PrimaryRow = meta.N + 1 }, psi), "primary row beyond the last row");
        Assert.Throws<ArgumentOutOfRangeException>(() => PsiStream.Build(bwt, meta with { PrimaryRow = -1 }, psi), "negative primary row");
        Assert.Throws<ArgumentOutOfRangeException>(() => PsiStream.Build(bwt, meta, psi, new PsiOptions { WindowBytes = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => PsiStream.Build(bwt, meta, psi, new PsiOptions { SymbolsPerPass = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => PsiStream.Build(bwt, meta, psi, new PsiOptions { ReadBytes = 0 }));
        File.AppendAllText(bwt, "x");
        Assert.Throws<InvalidDataException>(() => PsiStream.Build(bwt, meta, psi), "BWT file longer than n + 1");
        File.WriteAllBytes(bwt, new byte[3]);
        Assert.Throws<InvalidDataException>(() => PsiStream.Build(bwt, meta, psi), "BWT file shorter than n + 1");
    }

    [Test]
    public void MappedWorkingSet_IsBoundedByTheCap()
    {
        Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(600, 256, 3));
        using var dir = new Ph6Dir();
        (string bwt, BwtMeta meta) = WriteBwt(dir, c);
        string psi = dir.File("x.psi");
        // 256 cursors x 1 MiB = 256 MiB: over the 160 MiB cap of ADR-004 unless the pass is split into groups
        Assert.Throws<ArenaBudgetExceededException>(() => PsiStream.Build(bwt, meta, psi, new PsiOptions { WindowBytes = 1 << 20, SymbolsPerPass = 256 }));
        // the same alphabet in groups of 64 x 1 MiB = 64 MiB fits; a window larger than any cursor needs is not mapped larger than the cursor
        PsiStats st = PsiStream.Build(bwt, meta, psi, new PsiOptions { WindowBytes = 1 << 20, SymbolsPerPass = 64 });
        Assert.That(st.MappedBoundBytes, Is.EqualTo(64L << 20));
        Assert.That(File.ReadAllBytes(psi), Is.EqualTo(c.ExpectedPsiBytes(false)));
        // an explicit small cap and a shared budget are honoured
        var budget = new ArenaBudget(10_000);
        Assert.Throws<ArenaBudgetExceededException>(() => PsiStream.Build(bwt, meta, psi, new PsiOptions { WindowBytes = 1024, SymbolsPerPass = 64, Budget = budget, WindowCapBytes = 1L << 30 }), "read buffer + windows over the budget");
        Assert.That(budget.LiveBytes, Is.EqualTo(0), "everything the failed call reserved was released");
    }

    [Test]
    public void ArenaBudget_IsReleasedAfterABuild()
    {
        Ph6Case c = Ph6Case.OfNaive("mississippi");
        using var dir = new Ph6Dir();
        (string bwt, BwtMeta meta) = WriteBwt(dir, c);
        var budget = new ArenaBudget(160L << 20);
        PsiStats st = PsiStream.Build(bwt, meta, dir.File("x.psi"), new PsiOptions { Budget = budget });
        Assert.That(budget.LiveBytes, Is.EqualTo(0));
        // the read buffer (as many bytes as there are rows, the default read size is larger) plus the reservation of the windows
        Assert.That(budget.PeakBytes, Is.EqualTo(c.N + 1L + st.MappedBoundBytes));
    }
}
