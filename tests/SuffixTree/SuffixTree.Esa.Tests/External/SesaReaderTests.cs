using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.Facade;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Persistent;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The reader of the SESA container (work package 1b.7f, gate G1): a tree that reads the mapped sections in place is the tree of the text, on every member of <see cref="ISuffixTree"/> and of the navigation contract, against the
/// definitions, the Ukkonen tree, the persistent v6 tree and the in-memory ESA. The containers come from the end-to-end builder and from the test writer of <see cref="SesaTestFile"/>.
/// </summary>
[TestFixture]
public class SesaReaderTests
{
    private static readonly string[] Named =
    [
        "", "a", "aa", "ab", "ba", "aba", "abab", "banana", "mississippi", "abracadabra", "aaaaaaaaaa", "abababababab", "abcabxabcd", "xabxac", "GATTACAGATTACA", "ACGTACGTTACGT",
        "a$b$", "the quick brown fox jumps over the lazy dog", "aabaabaaabaabaaab", "café ÿres \u0001\u0000x",
    ];

    private static void Compare(string text, ulong seed, string context, int maxPatternLength = 6, EsaAccessMode access = EsaAccessMode.Checked)
    {
        using EsaSuffixTree tree = SesaTestFile.FileTree(text, psi: seed % 2 == 0, access: access);
        using var v6 = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        FacadeChecks.CheckEverything(tree, UkkonenTree.Build(text), v6, text, seed, $"[file] {context}", maxPatternLength);
        FacadeChecks.CheckAgainstTree(tree, memory, "esa-memory", text, $"[file] {context}");
    }

    [Test]
    public void NamedTexts_EveryMember_EqualTheOracles()
    {
        ulong seed = 1;
        foreach (string text in Named) Compare(text, seed++, "named");
    }

    [Test]
    public void NamedTexts_Trusted_GiveTheSameAnswers()
    {
        ulong seed = 100;
        foreach (string text in Named) Compare(text, seed++, "named, trusted", 5, EsaAccessMode.Trusted);
    }

    [Test]
    public void AdversarialFamilies_EveryMember_N150()
    {
        foreach ((string name, string text) in StringFamilies.Adversarial(150, 3))
            if (SesaTestFile.IsByteText(text)) Compare(text, 7, name, 4);
    }

    [Test]
    public void AdversarialFamilies_SmallSizes_EveryMember()
    {
        foreach (int n in new[] { 1, 2, 3, 4, 5, 8, 13 })
            foreach ((string name, string text) in StringFamilies.Adversarial(n, 11))
                if (SesaTestFile.IsByteText(text)) Compare(text, 5, $"{name} n={n}", 4);
    }

    [Test]
    public void RandomStrings_EveryMember()
    {
        var rng = new SplitMix64(2025);
        for (int i = 0; i < 40; i++)
        {
            int n = rng.NextInt(120);
            int sigma = new[] { 1, 2, 3, 4, 4, 20 }[rng.NextInt(6)];
            Compare(StringFamilies.Random(n, sigma, rng.Next()), rng.Next(), $"random {i} sigma {sigma} n {n}", 5);
        }
    }

    [Test]
    public void ARepeatThatNeedsTheOverflowList_IsReadBack()
    {
        // 4,000 random symbols and a copy of 320 of them: the LCP values of the copy reach 255 and more, so the compact LCP has an overflow list
        string random = StringFamilies.Random(3_700, 4, 17);
        string text = random + random[..320];
        LcpChildCase c = LcpChildCase.Of(text);
        Assert.That(c.Lcp.OverflowCount, Is.GreaterThan(0), "the text must have an overflow list");
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, SesaTestFile.Image(c), EsaValidationLevel.Structure);
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        FacadeChecks.CheckAgainstTree(tree, memory, "esa-memory", text, "overflow");
        Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo(memory.LongestRepeatedSubstring()));
        Assert.That(tree.CountDistinctSubstrings(), Is.EqualTo(memory.CountDistinctSubstrings()));
        foreach (string p in new[] { random[..300], random[10..400], text[^330..], "q" })
        {
            Assert.That(tree.CountOccurrences(p), Is.EqualTo(memory.CountOccurrences(p)), p.Length.ToString());
            Assert.That(tree.FindAllOccurrences(p), Is.EquivalentTo(memory.FindAllOccurrences(p)), p.Length.ToString());
        }
    }

    [Test]
    public void TheEndToEndBuilder_ContainersAreReadBack_AtEveryLevel()
    {
        var rng = new SplitMix64(99);
        for (int v = 0; v < 24; v++)
        {
            int n = 1 + rng.NextInt(200);
            int sigma = new[] { 1, 2, 3, 4, 5, 9 }[rng.NextInt(6)];
            string text = StringFamilies.Random(n, sigma, rng.Next());
            bool psi = v % 3 != 0;
            using var dir = new TempDir();
            SesaRun.Outcome o = SesaRun.Build(dir, text, SesaRun.Tiny(v, 1 + (v % 3), psi));
            foreach (EsaValidationLevel level in Enum.GetValues<EsaValidationLevel>())
            {
                using EsaSuffixTree tree = SesaTestFile.Open(o.Target, level, v % 2 == 0 ? EsaAccessMode.Checked : EsaAccessMode.Trusted, sha: level == EsaValidationLevel.Structure);
                using var v6 = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));
                FacadeChecks.CheckEverything(tree, UkkonenTree.Build(text), v6, text, (ulong)v, $"[builder v={v} {level}] sigma {sigma} n {n}", 4);
            }
        }
    }

    [Test]
    public void TheDefaultOptions_AreHeaderAndChecked()
    {
        var options = new EsaOpenOptions();
        Assert.That((options.Validation, options.Access, options.VerifyTextSha256), Is.EqualTo((EsaValidationLevel.Header, EsaAccessMode.Checked, false)));
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, "banana");
        using EsaSuffixTree tree = EsaSuffixTreeFactory.Open(path);
        Assert.That(tree.File, Is.Not.Null);
        Assert.That((tree.File!.Report.Level, tree.File.Report.Access), Is.EqualTo((EsaValidationLevel.Header, EsaAccessMode.Checked)));
    }

    [TestCase(EsaValidationLevel.Header, 0, 0)]
    [TestCase(EsaValidationLevel.Sections, 1, 0)]
    [TestCase(EsaValidationLevel.Structure, 1, 1)]
    public void TheReport_SaysWhichLevelRan(EsaValidationLevel level, int sections, int structure)
    {
        using var dir = new TempDir();
        string text = StringFamilies.Random(500, 4, 3);
        string path = SesaTestFile.Write(dir, text);
        using EsaSuffixTree tree = SesaTestFile.Open(path, level);
        SesaOpenReport r = tree.File!.Report;
        Assert.That(r.Level, Is.EqualTo(level));
        Assert.That(r.HeaderSeconds, Is.GreaterThanOrEqualTo(0));
        Assert.That(r.SectionBytesChecksummed > 0, Is.EqualTo(sections == 1), "the checksum pass runs from the Sections level");
        Assert.That(r.BitmapBytes > 0, Is.EqualTo(structure == 1), "the permutation bitmap exists at the Structure level only");
        if (structure == 1) Assert.That(r.BitmapBytes, Is.EqualTo((text.Length + 1 + 7) / 8), "n + 1 bits");
    }

    [Test]
    public void TheChecksumPass_ReadsEveryPayloadByte()
    {
        using var dir = new TempDir();
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(StringFamilies.Random(800, 4, 5)));
        SesaOracle o = SesaOracle.FromBytes(image);
        string path = dir.File("a.sesa");
        File.WriteAllBytes(path, image);
        using EsaSuffixTree tree = SesaTestFile.Open(path, EsaValidationLevel.Sections);
        Assert.That(tree.File!.Report.SectionBytesChecksummed, Is.EqualTo(o.Sections.Where(s => s.HasCrc).Sum(s => s.Length)));
    }

    [Test]
    public void TheTablePath_OfTheChecksum_GivesTheSameVerdict()
    {
        using var dir = new TempDir();
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(StringFamilies.Random(800, 4, 6)));
        string path = dir.File("a.sesa");
        File.WriteAllBytes(path, image);
        using (SesaFile f = SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Sections }, hardwareCrc: false))
        {
            Assert.That(f.Report.HardwareCrc, Is.False, "the table path was forced");
        }
        byte[] damaged = (byte[])image.Clone();
        damaged[SesaTestFile.Section(damaged, 3).Offset + 20] ^= 0x10;
        File.WriteAllBytes(path, damaged);
        Assert.That(() => SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Sections }, hardwareCrc: false), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Sections }, hardwareCrc: true), Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void Dispose_ReleasesTheMapping_IsIdempotent_AndEveryMemberThrowsAfterwards()
    {
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, "mississippi");
        EsaSuffixTree tree = SesaTestFile.Open(path);
        if (ProcMaps.Contains(path) is { } mapped) Assert.That(mapped, Is.True, "the file is mapped while the tree is open");
        Assert.That(tree.CountOccurrences("ss"), Is.EqualTo(2));
        tree.Dispose();
        tree.Dispose();
        if (ProcMaps.Contains(path) is { } still) Assert.That(still, Is.False, "the mapping is released by Dispose");
        Assert.That(() => tree.CountOccurrences("ss"), Throws.TypeOf<ObjectDisposedException>());
        Assert.That(() => tree.LongestRepeatedSubstring(), Throws.TypeOf<ObjectDisposedException>());
        Assert.That(() => tree.Traverse(new NullVisitor()), Throws.TypeOf<ObjectDisposedException>());
        File.Delete(path);
    }

    [Test]
    public void OneMappingOfTheWholeFile_IsMadeAndNothingElse()
    {
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, StringFamilies.Random(300, 4, 8));
        if (ProcMaps.Contains(path) is null) Assert.Ignore("no /proc/self/maps on this platform");
        using EsaSuffixTree tree = SesaTestFile.Open(path, EsaValidationLevel.Structure);
        string full = Path.GetFullPath(path);
        int lines = File.ReadLines("/proc/self/maps").Count(l => l.EndsWith(full, StringComparison.Ordinal));
        Assert.That(lines, Is.EqualTo(1), "one view of the whole file (ADR-003 section 9)");
    }

    [Test]
    public void TheTableIsExposed_InFileOrder_WithTheHeaderFacts()
    {
        using var dir = new TempDir();
        string text = StringFamilies.Random(300, 4, 9);
        LcpChildCase c = LcpChildCase.Of(text);
        string path = dir.File("a.sesa");
        File.WriteAllBytes(path, SesaTestFile.Image(c));
        using SesaFile f = SesaFile.Open(path, new EsaOpenOptions());
        Assert.That(f.Header.TextLength, Is.EqualTo(300));
        Assert.That(f.Header.NodeCount, Is.EqualTo(c.Statistics.NodeCount));
        Assert.That(f.Header.LrsDepth, Is.EqualTo(c.Statistics.LrsDepth));
        Assert.That(f.Header.Sigma, Is.EqualTo(c.Sigma));
        Assert.That(f.Header.PosWidth, Is.EqualTo(4));
        Assert.That(f.Sections.Select(s => s.Kind), Is.EqualTo(new[] { 1, 2, 13, 3, 6, 12, 11, 4 }));
        Assert.That(f.Sections.Select(s => s.Offset), Is.Ordered.Ascending);
    }

    [Test]
    public void AnEmptyText_Opens_AndIsEmpty()
    {
        using EsaSuffixTree tree = SesaTestFile.FileTree(string.Empty);
        Assert.That(tree.IsEmpty, Is.True);
        Assert.That(tree.NodeCount, Is.EqualTo(1));
        Assert.That(tree.CountOccurrences("a"), Is.Zero);
        Assert.That(tree.LongestRepeatedSubstring(), Is.Empty);
        Assert.That(tree.GetAllSuffixes(), Is.Empty);
    }

    [Test]
    public void TheSameFile_CanBeOpenedTwice_AndReadFromManyThreads()
    {
        using var dir = new TempDir();
        string text = StringFamilies.Random(2_000, 4, 10);
        string path = SesaTestFile.Write(dir, text);
        using EsaSuffixTree a = SesaTestFile.Open(path);
        using EsaSuffixTree b = SesaTestFile.Open(path, EsaValidationLevel.Structure);
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        var patterns = Enumerable.Range(0, 200).Select(i => text.Substring(i * 7, 3 + (i % 9))).ToArray();
        var expected = patterns.Select(p => memory.CountOccurrences(p)).ToArray();
        Parallel.For(0, 8, new ParallelOptions { MaxDegreeOfParallelism = 8 }, t =>
        {
            for (int round = 0; round < 20; round++)
            {
                for (int i = 0; i < patterns.Length; i++)
                {
                    Assert.That(a.CountOccurrences(patterns[i]), Is.EqualTo(expected[i]));
                    Assert.That(b.FindAllOccurrences(patterns[i]).Count, Is.EqualTo(expected[i]));
                }
                _ = a.FindMaximalExactMatches(text.Substring(t * 50, 120), 8);
            }
        });
    }

    private sealed class NullVisitor : ISuffixTreeVisitor
    {
        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
        }

        public void EnterBranch(int key)
        {
        }

        public void ExitBranch()
        {
        }
    }
}
