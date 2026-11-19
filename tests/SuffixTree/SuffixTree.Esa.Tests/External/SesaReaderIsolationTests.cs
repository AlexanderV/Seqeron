using System.Buffers.Binary;
using System.Diagnostics;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The tests that the mutation audit of work package 1b.7f asked for (round 1): a rule of the reader that another check of the same pass also catches (the permutation and the order and the Phi pass see the same damaged array;
/// the descent, the pyramid and the safe views see the same damaged table) survives every test that only asks for "an <see cref="InvalidDataException"/>". Each case here damages ONE thing, seals the file again so that no checksum
/// sees it, and pins the message of the rule that is meant to see it, so that a reader without that rule fails the test even when a later check would still refuse the file with another message.
/// </summary>
[TestFixture]
public class SesaReaderIsolationTests
{
    private const string Mississippi = "mississippi"; // SA 11 10 7 4 1 0 9 8 6 3 5 2, LCP 0 0 1 1 4 0 0 1 0 2 1 3

    private static byte[] Image(string text, bool psi = true, bool counts = true, bool buildInfo = true) => SesaTestFile.Image(LcpChildCase.Of(text), psi, counts, buildInfo);

    private static void SetI32(byte[] image, long at, int value) => BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan((int)at), value);

    private static void SetI64(byte[] image, long at, long value) => BinaryPrimitives.WriteInt64LittleEndian(image.AsSpan((int)at), value);

    /// <summary>Opens the image at the level and returns the message of the <see cref="InvalidDataException"/> that refuses it.</summary>
    private static string Refusal(byte[] image, EsaValidationLevel level, bool sha = false)
    {
        using var dir = new TempDir();
        string path = dir.File("i.sesa");
        File.WriteAllBytes(path, image);
        InvalidDataException? e = Assert.Throws<InvalidDataException>(() =>
        {
            using EsaSuffixTree tree = SesaTestFile.Open(path, level, EsaAccessMode.Checked, sha);
        });
        return e!.Message;
    }

    private static void Opens(byte[] image, EsaValidationLevel level)
    {
        using var dir = new TempDir();
        string path = dir.File("o.sesa");
        File.WriteAllBytes(path, image);
        using EsaSuffixTree tree = SesaTestFile.Open(path, level, EsaAccessMode.Checked, sha: true);
        Assert.That(tree, Is.Not.Null);
    }

    private static T WithTimeout<T>(Func<T> work, int seconds = 20)
    {
        Task<T> task = Task.Run(work);
        Assert.That(task.Wait(TimeSpan.FromSeconds(seconds)), Is.True, "the call does not return");
        return task.Result;
    }

    /// <summary>The image of a text with one thing damaged and every checksum recomputed.</summary>
    private static byte[] Sealed(string text, Action<byte[], SesaOracle> damage, bool psi = true, bool counts = true, bool buildInfo = true)
    {
        byte[] image = Image(text, psi, counts, buildInfo);
        damage(image, SesaOracle.FromBytes(image));
        SesaTestFile.Reseal(image);
        return image;
    }

    // ---------------------------------------------------------------- the numbers of the header page: each rule on its own (the Header level reads no payload)

    [TestCase(0)]
    [TestCase(1)]
    [TestCase(2)]
    [TestCase(3)]
    [TestCase(5)]
    [TestCase(6)]
    [TestCase(7)]
    [TestCase(9)]
    [TestCase(16)]
    [TestCase(255)]
    [TestCase(-1)]
    public void APositionWidthThatIsNeitherFourNorEight_IsCorrupt_AtHeader(int width)
    {
        byte[] image = Image(Mississippi);
        SetI32(image, 72, width);
        SesaTestFile.ResealHeader(image);
        Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain("POS_WIDTH"));
    }

    [TestCase("ab", 3)]
    [TestCase("abc", 4)]
    [TestCase(Mississippi, 12)]
    [TestCase(Mississippi, 65_537)]
    [TestCase(Mississippi, -1)]
    public void ASigmaAboveTheTextLength_OrBelowZero_IsCorrupt_AtHeader(string text, int sigma)
    {
        byte[] image = Image(text);
        SetI32(image, 76, sigma);
        SesaTestFile.ResealHeader(image);
        Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain("SIGMA is"));
    }

    [TestCase(24, 2L)]
    [TestCase(24, 0L)]
    [TestCase(56, 1L)]
    [TestCase(64, 0L)]
    [TestCase(64, 5L)]
    public void TheStatisticsOfAnEmptyText_AreFixed_AtHeader(int at, long value)
    {
        byte[] image = Image(string.Empty);
        SetI64(image, at, value);
        SesaTestFile.ResealHeader(image);
        Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain("empty text"));
    }

    [TestCase(12L)]
    [TestCase(24L)]
    [TestCase(-3L)]
    public void ANodeCountOutsideTheTreesBounds_IsCorrupt_AtHeader(long nodeCount)
    {
        byte[] image = Image(Mississippi);
        SetI64(image, 24, nodeCount);
        SesaTestFile.ResealHeader(image);
        Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain("NODE_COUNT"));
    }

    [TestCase(11L)]
    [TestCase(12L)]
    [TestCase(-1L)]
    public void ALongestRepeatAsLongAsTheText_IsCorrupt_AtHeader(long depth)
    {
        byte[] image = Image(Mississippi);
        SetI64(image, 56, depth);
        SesaTestFile.ResealHeader(image);
        Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain("LRS_DEPTH"));
    }

    [TestCase(11L)]
    [TestCase(12L)]
    [TestCase(-2L)]
    public void ARepeatRankAtTheTextLength_IsCorrupt_AtHeader(long rank)
    {
        byte[] image = Image(Mississippi);
        SetI64(image, 64, rank);
        SesaTestFile.ResealHeader(image);
        Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain("LRS_RANK"));
    }

    [Test]
    public void ARepeatRankWithoutARepeat_IsCorrupt_AtHeader()
    {
        byte[] image = Image(Mississippi);
        SetI64(image, 56, 0);
        SetI64(image, 64, 3);
        SesaTestFile.ResealHeader(image);
        Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain("LRS_RANK"));
    }

    // ---------------------------------------------------------------- the structure pass: each rule on its own (the file is sealed, so only the pass can see it)

    [Test]
    public void TheSuffixArray_ThatDoesNotStartWithTheEmptySuffix_IsRefusedByItsFirstRule()
    {
        // the n other entries are still a permutation of 0..n-1: nothing but the rule on SA[0] sees that the empty suffix is missing
        byte[] image = Sealed(Mississippi, (img, o) => SetI32(img, o.Get(SesaFormat.KindSa).Offset, 3));
        Assert.That(Refusal(image, EsaValidationLevel.Structure), Does.Contain("SA[0] is"));
    }

    [Test]
    public void ThePositionN_AtARankAboveZero_IsRefusedAsThePositionAtTheEnd()
    {
        byte[] image = Sealed(Mississippi, (img, o) => SetI32(img, o.Get(SesaFormat.KindSa).Offset + 4, 11));
        Assert.That(Refusal(image, EsaValidationLevel.Structure), Does.Contain("empty suffix at rank 0 only"));
    }

    [Test]
    public void ARepeatedPosition_IsRefusedAsNotAPermutation()
    {
        byte[] image = Sealed(Mississippi, (img, o) => SetI32(img, o.Get(SesaFormat.KindSa).Offset + (4 * 3), 7));
        Assert.That(Refusal(image, EsaValidationLevel.Structure), Does.Contain("occurs twice"));
    }

    [Test]
    public void AnLcpLongerThanTheShorterSuffix_IsRefusedWhereItIsRead()
    {
        byte[] image = Sealed(Mississippi, (img, o) => img[o.Get(SesaFormat.KindLcp).Offset + 3] = 200);
        Assert.That(Refusal(image, EsaValidationLevel.Structure), Does.Contain("longer than the shorter of the two suffixes"));
    }

    [TestCase(0)]
    [TestCase(1)]
    public void TheFirstTwoLcpEntries_AreZero(int rank)
    {
        byte[] image = Sealed(Mississippi, (img, o) => img[o.Get(SesaFormat.KindLcp).Offset + rank] = 1);
        Assert.That(Refusal(image, EsaValidationLevel.Structure), Does.Contain("LCP[0] and LCP[1] must be 0"));
    }

    [Test]
    public void TwoAdjacentSuffixes_InTheWrongOrder_AreRefusedByTheOrderRule()
    {
        byte[] image = Sealed(Mississippi, (img, o) =>
        {
            long sa = o.Get(SesaFormat.KindSa).Offset;
            int a = BinaryPrimitives.ReadInt32LittleEndian(img.AsSpan((int)sa + (4 * 4)));
            int b = BinaryPrimitives.ReadInt32LittleEndian(img.AsSpan((int)sa + (4 * 5)));
            SetI32(img, sa + (4 * 4), b);
            SetI32(img, sa + (4 * 5), a);
        });
        Assert.That(Refusal(image, EsaValidationLevel.Structure), Does.Contain("not in strictly ascending order"));
    }

    [Test]
    public void AnLcpThatIsTooShort_IsRefusedByTheOrderRule_WhichMeasuresTheCommonPrefixItself()
    {
        // LCP[4] is 4 (issippi against ississippi); with 3 the suffixes are still in order, and only the measured common prefix differs from the stored one
        byte[] image = Sealed(Mississippi, (img, o) => img[o.Get(SesaFormat.KindLcp).Offset + 4] = 3);
        Assert.That(Refusal(image, EsaValidationLevel.Structure), Does.Contain("not in strictly ascending order after their 3 common symbols"));
    }

    [Test]
    public void ALetterOfTheAlphabetThatDoesNotOccur_IsRefused_WhenNoCountsAreStored()
    {
        // COUNTS is optional; without it only the recount of the text sees a letter that is listed and never used
        LcpChildCase c = LcpChildCase.Of(Mississippi);
        byte[] image = SesaTestFile.Image(c, psi: true, counts: false);
        byte[] symbols = [.. c.Symbols.Select(s => s == 2 ? (byte)0 : s)];
        OracleSection text = SesaTestFile.Section(image, SesaFormat.KindText);
        TextOracle.Pack(symbols, c.Bits).CopyTo(image.AsSpan((int)text.Offset));
        SesaTestFile.Reseal(image);
        Assert.That(Refusal(image, EsaValidationLevel.Structure), Does.Contain("does not occur in the text"));
    }

    [Test]
    public void TheCounts_ThatDisagreeWithTheText_AreRefusedByTheRecount()
    {
        byte[] image = Sealed(Mississippi, (img, o) =>
        {
            OracleSection counts = o.Get(SesaFormat.KindCounts);
            BinaryPrimitives.WriteUInt64LittleEndian(img.AsSpan((int)counts.Offset), BinaryPrimitives.ReadUInt64LittleEndian(img.AsSpan((int)counts.Offset)) + 1);
            BinaryPrimitives.WriteUInt64LittleEndian(img.AsSpan((int)counts.Offset + 24), BinaryPrimitives.ReadUInt64LittleEndian(img.AsSpan((int)counts.Offset + 24)) - 1); // i 4 -> 5, s 4 -> 3: the sum is still n and no count is 0
        });
        Assert.That(Refusal(image, EsaValidationLevel.Structure), Does.Contain("COUNTS says"));
    }

    [Test]
    public void ABuildInfoThatIsNotUtf8_OrHasTheWrongNumberOfLines_IsRefused()
    {
        byte[] bad = Sealed(Mississippi, (img, o) => img[o.Get(SesaFormat.KindBuildInfo).Offset] = 0xFF);
        Assert.That(Refusal(bad, EsaValidationLevel.Structure), Does.Contain("not UTF-8"));
        byte[] lines = Sealed(Mississippi, (img, o) =>
        {
            OracleSection info = o.Get(SesaFormat.KindBuildInfo);
            int newline = Array.IndexOf(img, (byte)'\n', (int)info.Offset);
            img[newline] = (byte)' ';
        });
        Assert.That(Refusal(lines, EsaValidationLevel.Structure), Does.Contain("lines but its entry says"));
    }

    [Test]
    public void TheArraysOfAnEmptyText_AreCheckedByTheStructureLevel()
    {
        string[] expected = ["suffix array of an empty text", "PSI array of an empty text", "LCP[0] must be 0", "child table of an empty text"];
        int[] kinds = [SesaFormat.KindSa, SesaFormat.KindPsi, SesaFormat.KindLcp, SesaFormat.KindChild];
        for (int i = 0; i < kinds.Length; i++)
        {
            int kind = kinds[i];
            byte[] image = Sealed(string.Empty, (img, o) => img[o.Get(kind).Offset] = 1);
            Opens(Image(string.Empty), EsaValidationLevel.Structure);
            Assert.That(Refusal(image, EsaValidationLevel.Structure), Does.Contain(expected[i]), $"section {kind}");
        }
    }

    [TestCase(int.MaxValue - 3L, false)]
    [TestCase(int.MaxValue - 2L, false)]
    [TestCase(int.MaxValue - 1L, true)]
    [TestCase((long)int.MaxValue, true)]
    [TestCase(uint.MaxValue - 3L, true)]
    public void TheStructureLevel_IsNotSupported_ForATextThatTheThirtyTwoBitPassCannotIndex(long symbols, bool refused)
    {
        // a container of 2^31 symbols needs 20 GB of files: the rule is tested on its own, the pass calls it before it touches a payload
        if (refused) Assert.That(() => SesaVerifier.RequireSupportedLength(symbols), Throws.TypeOf<NotSupportedException>().With.Message.Contains("64-bit"));
        else Assert.That(() => SesaVerifier.RequireSupportedLength(symbols), Throws.Nothing);
    }

    // ---------------------------------------------------------------- the section checksums

    [Test]
    public void ASectionWithoutTheCrcFlag_IsNotChecksummed_AndTheFileOpens()
    {
        byte[] image = Image(Mississippi);
        OracleSection info = SesaTestFile.Section(image, SesaFormat.KindBuildInfo);
        SetI32(image, 128 + (64L * info.Index) + 4, 0); // flags: optional, no CRC
        BinaryPrimitives.WriteUInt32LittleEndian(image.AsSpan(128 + (64 * info.Index) + 56), 0);
        SesaTestFile.ResealHeader(image);
        Opens(image, EsaValidationLevel.Sections);
        Opens(image, EsaValidationLevel.Structure);
    }

    [TestCase(true)]
    [TestCase(false)]
    public void ASectionLongerThanOnePieceOfTheChecksum_IsChained_NotRestarted(bool hardware)
    {
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, StringFamilies.Random(900, 4, 12));
        long before = SesaVerifierTuning.CrcPieceBytes;
        try
        {
            SesaVerifierTuning.CrcPieceBytes = 64;
            using SesaFile file = SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Sections }, hardware);
            Assert.That(file.Report.SectionBytesChecksummed, Is.GreaterThan(64 * 10), "sections of many pieces were checksummed");
            Assert.That(file.Report.HardwareCrc, Is.EqualTo(hardware));
        }
        finally
        {
            SesaVerifierTuning.CrcPieceBytes = before;
        }
    }

    // ---------------------------------------------------------------- the mapping: nothing is left behind by a failure

    [Test]
    public void AFailedOpen_LeavesNothingMapped_AtEveryLevel()
    {
        using var dir = new TempDir();
        string probe = SesaTestFile.Write(dir, StringFamilies.Random(200, 4, 13), "probe.sesa");
        if (ProcMaps.Contains(probe) is null) Assert.Ignore("no /proc/self/maps on this platform");
        byte[] clean = Image(StringFamilies.Random(200, 4, 14));
        byte[] magic = (byte[])clean.Clone();
        magic[0] ^= 1;
        byte[] body = (byte[])clean.Clone();
        body[SesaTestFile.Section(body, SesaFormat.KindSa).Offset + 12] ^= 0x10;
        byte[] swapped = Sealed(StringFamilies.Random(200, 4, 14), (img, o) =>
        {
            long sa = o.Get(SesaFormat.KindSa).Offset;
            int a = BinaryPrimitives.ReadInt32LittleEndian(img.AsSpan((int)sa + 40));
            int b = BinaryPrimitives.ReadInt32LittleEndian(img.AsSpan((int)sa + 44));
            SetI32(img, sa + 40, b);
            SetI32(img, sa + 44, a);
        });
        foreach ((string name, byte[] image, EsaValidationLevel level) in new[] { ("magic", magic, EsaValidationLevel.Header), ("crc", body, EsaValidationLevel.Sections), ("structure", swapped, EsaValidationLevel.Structure) })
        {
            string path = dir.File(name + ".sesa");
            File.WriteAllBytes(path, image);
            for (int i = 0; i < 5; i++) Assert.Throws<InvalidDataException>(() => SesaFile.Open(path, new EsaOpenOptions { Validation = level }).Dispose(), name);
            Assert.That(ProcMaps.Contains(path), Is.False, $"a failed open at {level} ({name}) left the file mapped");
        }
    }

    [Test]
    public void AFileThatCannotSeek_IsRefusedByTheReaderAndByTheDispatcher_BeforeAnythingIsMapped()
    {
        if (!OperatingSystem.IsLinux()) Assert.Ignore("a named pipe");
        using var dir = new TempDir();
        string fifo = dir.File("p.fifo");
        using (Process? mk = Process.Start(new ProcessStartInfo("mkfifo", fifo) { UseShellExecute = false }))
        {
            if (mk is null) Assert.Ignore("no mkfifo");
            mk.WaitForExit();
            if (mk.ExitCode != 0) Assert.Ignore("mkfifo failed");
        }
        foreach ((string name, Action open) in new (string, Action)[] { ("reader", () => SesaFile.Open(fifo, new EsaOpenOptions()).Dispose()), ("dispatcher", () => ((IDisposable)SuffixIndexFactory.Open(fifo)).Dispose()) })
        {
            // opening a pipe for reading waits for a writer: give it one
            Task writer = Task.Run(() =>
            {
                using var w = new FileStream(fifo, FileMode.Open, FileAccess.Write);
            });
            Exception? caught = WithTimeout(() => Assert.Catch(() => open()));
            Assert.That(caught, Is.TypeOf<InvalidDataException>(), name);
            Assert.That(caught!.Message, Does.Contain("not a regular file"), name);
            Assert.That(writer.Wait(TimeSpan.FromSeconds(20)), Is.True, name);
        }
    }

    [Test]
    public void ATextReaderThatOutlivesItsFile_IsDisposedWithIt_AndNeverReadsTheUnmappedMemory()
    {
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, StringFamilies.Random(500, 4, 15));
        SesaFile file = SesaFile.Open(path, new EsaOpenOptions());
        MappedTextReader reader = file.TextReader;
        Assert.That(reader.Window(0), Is.Not.Zero, "a window of the live text");
        file.Dispose();
        Assert.That(() => reader.Window(0), Throws.TypeOf<ObjectDisposedException>(), "the reader of a disposed file must not read its unmapped pages");
    }

    // ---------------------------------------------------------------- the descent in checked mode: a damaged cell or depth on the path of a pattern

    [Test]
    public void AChildAsDeepAsItsParent_IsRefusedByTheDescent_TheWalkAndTheNavigators()
    {
        // 's' = [8, 11] has depth 1; its first child [8, 9] ('si', depth 2, first l-index 9) is made as deep as its parent by LCP[9] = 1. The sibling chain of 's' starts at the cell 8 and never looks at rank 9, so the damage is seen
        // only where the depth of the child is read: in the descent, in the walk, in the contract's navigation and in the generic navigator
        byte[] image = Sealed(Mississippi, (img, o) => img[o.Get(SesaFormat.KindLcp).Offset + 9] = 1);
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image);
        Exception? count = WithTimeout(() => Assert.Catch(() => tree.CountOccurrences("si")), 10);
        Assert.That(count, Is.TypeOf<InvalidDataException>(), "the descent");
        var nav = (ISuffixTreeNavigation)tree;
        Exception? contract = WithTimeout(() => Assert.Catch(() => Descend(nav, nav.Root, 0)), 10);
        Assert.That(contract, Is.TypeOf<InvalidDataException>(), "the contract's navigation");
        Exception? printed = WithTimeout(() => Assert.Catch(() => tree.PrintTree()), 10);
        Assert.That(printed, Is.TypeOf<InvalidDataException>(), "the generic navigator (PrintTree)");
        AssertWalkStopsBeforeTheDamagedNode(tree);
    }

    [Test]
    public void AChildDeeperThanTheText_IsRefusedByTheDescent_AndTheWalk()
    {
        byte[] image = Sealed(Mississippi, (img, o) => img[o.Get(SesaFormat.KindLcp).Offset + 9] = 200);
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image);
        Exception? count = WithTimeout(() => Assert.Catch(() => tree.CountOccurrences("si")), 10);
        Assert.That(count, Is.TypeOf<InvalidDataException>(), "the descent");
        AssertWalkStopsBeforeTheDamagedNode(tree);
    }

    /// <summary>
    /// The walk over a damaged table is refused, and it is refused right after the damaged branch was entered: no node is reported whose depth or links the walk has not proved (the occurrences of the nodes come from the whole LCP array, so what was
    /// reported before is compared with the clean walk only where the damage is in the child table, see <see cref="TheWalk_ReportsNothing_OfAnIntervalWhoseFirstLIndexLiesOutsideIt"/>).
    /// </summary>
    private static void AssertWalkStopsBeforeTheDamagedNode(EsaSuffixTree damagedTree)
    {
        var damaged = new RecordingVisitor();
        Exception? walk = WithTimeout(() => Assert.Catch(() => damagedTree.Traverse(damaged)), 10);
        Assert.That(walk, Is.TypeOf<InvalidDataException>(), "the walk");
        Assert.That(damaged.Events, Is.Not.Empty);
        Assert.That(damaged.Events[^1], Does.StartWith("enter"), "the last thing the walk reported is the branch that it then refused");
    }

    [Test]
    public void ALeafWhoseSuffixIsNoLongerThanTheDepthOfItsParent_IsRefusedByTheNavigator()
    {
        // SA[11] = 9 puts the suffix "pi" (a leaf of depth 3 with the terminator) under the node "ssi" of depth 3
        byte[] image = Sealed(Mississippi, (img, o) => SetI32(img, o.Get(SesaFormat.KindSa).Offset + (4 * 11), 9));
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image);
        var nav = (ISuffixTreeNavigation)tree;
        Exception? e = WithTimeout(() => Assert.Catch(() => Descend(nav, nav.Root, 0)), 10);
        Assert.That(e, Is.TypeOf<InvalidDataException>(), "the contract's navigation");
        Exception? printed = WithTimeout(() => Assert.Catch(() => tree.PrintTree()), 10); // the generic navigator (EsaNavigator) makes its own children
        Assert.That(printed, Is.TypeOf<InvalidDataException>(), "PrintTree");
    }

    private static void Descend(ISuffixTreeNavigation nav, IndexLocus node, int depth)
    {
        if (depth > 8) return;
        var children = new IndexLocus[64];
        int count = nav.GetChildren(node, children);
        for (int i = 0; i < Math.Min(count, children.Length); i++) Descend(nav, children[i], depth + 1);
    }

    [Test]
    public void ALargeSection_MustStartAtAMultipleOf4096_AtHeader()
    {
        // a section of 1 MiB or more is aligned to 4,096 (ADR-003 section 3): the SA of 270,000 symbols is 1.08 MB
        byte[] image = Image(StringFamilies.Random(270_000, 4, 33), psi: false, counts: true, buildInfo: false);
        OracleSection sa = SesaTestFile.Section(image, SesaFormat.KindSa);
        Assert.That(sa.Length, Is.GreaterThanOrEqualTo(1 << 20));
        Assert.That(sa.Offset % 4096, Is.Zero);
        SetI64(image, 128 + (64L * sa.Index) + 8, sa.Offset + 64);
        SesaTestFile.ResealHeader(image);
        Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain("not a multiple of 4096"));
    }

    [Test]
    public void ARootWhoseFirstLIndexIsOutsideTheText_IsRefused_BeforeTheWalkReportsAnything()
    {
        byte[] image = Sealed(Mississippi, (img, o) => SetI32(img, o.Get(SesaFormat.KindChild).Offset, 11 + 3));
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image);
        var visitor = new CountingVisitor();
        Exception? e = WithTimeout(() => Assert.Catch(() => tree.Traverse(visitor)), 10);
        Assert.That(e, Is.TypeOf<InvalidDataException>());
        Assert.That(visitor.Calls, Is.Zero, "a damaged root is found before the root node is reported");
    }

    [Test]
    public void AnAlphabetWithARepeatedLetter_OrALetterAbove255_IsCorrupt_AtHeader()
    {
        byte[] repeated = Image(Mississippi);
        OracleSection alphabet = SesaTestFile.Section(repeated, SesaFormat.KindAlphabet);
        BinaryPrimitives.WriteUInt16LittleEndian(repeated.AsSpan((int)alphabet.Offset + 2), BinaryPrimitives.ReadUInt16LittleEndian(repeated.AsSpan((int)alphabet.Offset)));
        SesaTestFile.Reseal(repeated);
        Assert.That(Refusal(repeated, EsaValidationLevel.Header), Does.Contain("not strictly ascending"));
        byte[] wide = Image(Mississippi);
        BinaryPrimitives.WriteUInt16LittleEndian(wide.AsSpan((int)alphabet.Offset + 6), 300);
        SesaTestFile.Reseal(wide);
        Assert.That(Refusal(wide, EsaValidationLevel.Header), Does.Contain("below 256"));
    }

    [Test]
    public void ASectionThatEndsPastTheFile_IsCorrupt_AtHeader_EvenByOneByte()
    {
        byte[] longer = Image(Mississippi);
        OracleSection sa = SesaTestFile.Section(longer, SesaFormat.KindSa);
        SetI64(longer, 128 + (64L * sa.Index) + 16, longer.Length + 1);
        SesaTestFile.ResealHeader(longer);
        Assert.That(Refusal(longer, EsaValidationLevel.Header), Does.Contain("lies outside the file"));
        byte[] shifted = Image(Mississippi);
        OracleSection last = SesaOracle.FromBytes(shifted).Sections.OrderBy(x => x.Offset).Last();
        SetI64(shifted, 128 + (64L * last.Index) + 8, last.Offset + 1); // the last section now ends one byte past the end of the file
        SesaTestFile.ResealHeader(shifted);
        Assert.That(Refusal(shifted, EsaValidationLevel.Header), Does.Contain("lies outside the file"));
        byte[] trailing = [.. Image(Mississippi), .. new byte[64]];
        SetI64(trailing, 32, trailing.Length); // FILE_SIZE says what the file has, but the last section ends 64 bytes before
        SesaTestFile.ResealHeader(trailing);
        Assert.That(Refusal(trailing, EsaValidationLevel.Header), Does.Contain("the last section ends at"));
    }

    [Test]
    public void ASectionKindThatOccursTwice_IsCorrupt_AtHeader()
    {
        byte[] image = Image(Mississippi);
        OracleSection info = SesaTestFile.Section(image, SesaFormat.KindBuildInfo);
        BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan(128 + (64 * info.Index)), SesaFormat.KindSa);
        SesaTestFile.ResealHeader(image);
        Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain("occurs twice"));
    }

    [Test]
    public void Locate_WithASkipAndASmallBuffer_EqualsTheIndexInMemory_OnBothWidthsOfPosition()
    {
        string text = StringFamilies.Random(400, 3, 16);
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        using EsaSuffixTree file = SesaTestFile.FileTree(text);
        foreach (string pattern in new[] { "a", "ab", "b", "abc", "ca" })
        {
            Assert.That(memory.Index.TryMatch(pattern, out int mlo, out int mhi), Is.True);
            Assert.That(file.Index.TryMatch(pattern, out int flo, out int fhi), Is.True);
            Assert.That((flo, fhi), Is.EqualTo((mlo, mhi)));
            int count = mhi - mlo + 1;
            foreach (long skip in new long[] { 0, 1, 2, count / 2, Math.Max(0, count - 1), count })
            {
                foreach (int size in new[] { 1, 3, 7, count + 2 })
                {
                    var wantInt = new int[size];
                    var gotInt = new int[size];
                    int wn = memory.Index.Locate(mlo, mhi, skip, wantInt);
                    int gn = file.Index.Locate(flo, fhi, skip, gotInt);
                    Assert.That(gn, Is.EqualTo(wn), $"'{pattern}' skip {skip} size {size}");
                    Assert.That(gotInt.AsSpan(0, gn).ToArray(), Is.EqualTo(wantInt.AsSpan(0, wn).ToArray()), $"'{pattern}' skip {skip} size {size}");
                    var wantLong = new long[size];
                    var gotLong = new long[size];
                    int wl = memory.Index.Locate(mlo, mhi, skip, wantLong);
                    int gl = file.Index.Locate(flo, fhi, skip, gotLong);
                    Assert.That(gotLong.AsSpan(0, gl).ToArray(), Is.EqualTo(wantLong.AsSpan(0, wl).ToArray()), $"'{pattern}' skip {skip} size {size} (64-bit)");
                }
            }
        }
    }

    [Test]
    public void TheLcpEntry_CountsAtMostOneEscapePerRank_AndNoneFor32BitEntries_AtHeader()
    {
        string random = StringFamilies.Random(1_500, 4, 21);
        byte[] u8 = Image(random + random[..300]); // U8 with an overflow list: rows = 1,801
        OracleSection lcp = SesaTestFile.Section(u8, SesaFormat.KindLcp);
        SetI64(u8, 128 + (64L * lcp.Index) + 40, 1_802);
        SesaTestFile.ResealHeader(u8);
        Assert.That(Refusal(u8, EsaValidationLevel.Header), Does.Contain("LCP PARAM0"));
        byte[] u32 = Image(new string('A', 200_000)); // 32-bit entries: no overflow list
        OracleSection wide = SesaTestFile.Section(u32, SesaFormat.KindLcp);
        Assert.That(wide.Codec, Is.EqualTo(3), "32-bit LCP entries");
        SetI64(u32, 128 + (64L * wide.Index) + 40, 1);
        SesaTestFile.ResealHeader(u32);
        Assert.That(Refusal(u32, EsaValidationLevel.Header), Does.Contain("LCP PARAM0"));
    }

    [Test]
    public void TheTextEntry_MustDescribeTheTextOfTheHeader_AtHeader()
    {
        string text = StringFamilies.Random(600, 4, 18);
        (int Field, long Value, string Message)[] cases =
        [
            (24, 601, "does not describe"),      // COUNT
            (32, 1, "does not describe"),        // ELEM_WIDTH
            (40, 3, "does not describe"),        // PARAM0 (bits per symbol)
            (48, 1, "does not describe"),        // PARAM1
        ];
        foreach ((int field, long value, string message) in cases)
        {
            byte[] image = Image(text);
            OracleSection t = SesaTestFile.Section(image, SesaFormat.KindText);
            if (field == 32) SetI32(image, 128 + (64L * t.Index) + field, (int)value);
            else SetI64(image, 128 + (64L * t.Index) + field, value);
            SesaTestFile.ResealHeader(image);
            Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain(message), $"field at +{field}");
        }
        byte[] longer = Image(text);
        OracleSection entry = SesaTestFile.Section(longer, SesaFormat.KindText);
        SetI64(longer, 128 + (64L * entry.Index) + 16, entry.Length + 8);
        SesaTestFile.ResealHeader(longer);
        Assert.That(Refusal(longer, EsaValidationLevel.Header), Does.Contain("TEXT has LENGTH"));
        byte[] wide = Image(text);
        SetI32(wide, 76, 300);
        SesaTestFile.ResealHeader(wide);
        Assert.That(Refusal(wide, EsaValidationLevel.Header), Does.Contain("SIGMA is 300, more than"));
    }

    [Test]
    public void TheOverflowEntry_AndTheOverflowList_AgreeAboutWhetherThereAreEscapes_AtHeader()
    {
        string random = StringFamilies.Random(1_500, 4, 21);
        byte[] image = Image(random + random[..300]);
        OracleSection lcp = SesaTestFile.Section(image, SesaFormat.KindLcp);
        SetI64(image, 128 + (64L * lcp.Index) + 40, 0);
        SesaTestFile.ResealHeader(image);
        Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain("no overflow entries"));
    }

    [Test]
    public void EveryArrayEntry_MustHaveItsCountAndItsWidth_AtHeader()
    {
        string random = StringFamilies.Random(1_500, 4, 21);
        string text = random + random[..300]; // every kind of section: an 8-bit LCP with an overflow list, PSI, COUNTS
        foreach ((int kind, string name) in new[] { (SesaFormat.KindAlphabet, "ALPHABET"), (SesaFormat.KindSa, "SA"), (SesaFormat.KindChild, "CHILD"), (SesaFormat.KindLcp, "LCP"), (SesaFormat.KindLcpOverflow, "LCP_OVERFLOW"), (SesaFormat.KindPsi, "PSI"), (SesaFormat.KindCounts, "COUNTS") })
        {
            foreach (int field in new[] { 24, 32 }) // COUNT, ELEM_WIDTH
            {
                byte[] image = Image(text);
                OracleSection s = SesaTestFile.Section(image, kind);
                long at = 128 + (64L * s.Index) + field;
                if (field == 24) SetI64(image, at, s.Count + 1);
                else SetI32(image, at, s.ElemWidth + 1);
                SesaTestFile.ResealHeader(image);
                Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain($"the {name} entry has COUNT"), $"{name}, field at +{field}");
            }
        }
    }

    [Test]
    public void ADuplicatedPosition_IsRefusedWhenTheInverseIsBuilt_InBothModes()
    {
        byte[] image = Sealed(Mississippi, (img, o) => SetI32(img, o.Get(SesaFormat.KindSa).Offset + (4 * 3), 7));
        foreach (EsaAccessMode mode in new[] { EsaAccessMode.Checked, EsaAccessMode.Trusted })
        {
            using var dir = new TempDir();
            using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Header, mode);
            Exception? e = WithTimeout(() => Assert.Catch(() => ((ISuffixArrayAccess)tree).RankOf(2)), 10);
            Assert.That(e, Is.TypeOf<InvalidDataException>(), mode.ToString());
            Assert.That(e!.Message, Does.Contain("not a permutation"), mode.ToString());
        }
    }

    [Test]
    public void APositionAsLargeAsTheCount_IsRefusedWhenTheInverseIsBuilt()
    {
        // n + 1 = 12 is one past the text and the empty suffix: not a position of the inverse
        byte[] image = Sealed(Mississippi, (img, o) => SetI32(img, o.Get(SesaFormat.KindSa).Offset + (4 * 3), 12));
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Header, EsaAccessMode.Trusted);
        Exception? e = WithTimeout(() => Assert.Catch(() => ((ISuffixArrayAccess)tree).RankOf(2)), 10);
        Assert.That(e, Is.TypeOf<InvalidDataException>());
        Assert.That(e!.Message, Does.Contain("not a permutation"));
    }

    [Test]
    public void TheInverse_ReadAtAPositionOutsideTheText_IsInvalidDataForAFile_AndAnArgumentErrorInMemory()
    {
        using EsaSuffixTree file = SesaTestFile.FileTree(Mississippi);
        using EsaSuffixTree memory = EsaSuffixTree.Build(Mississippi);
        Assert.That(() => file.Index.Navigation.Isa(100), Throws.TypeOf<InvalidDataException>().With.Message.Contains("inverse suffix array"));
        Assert.That(() => memory.Index.Navigation.Isa(100), Throws.InstanceOf<ArgumentException>());
    }

    [Test]
    public void CountsThatSumToLessThanTheText_AreCorrupt_AtHeader()
    {
        byte[] image = Sealed(Mississippi, (img, o) =>
        {
            long at = o.Get(SesaFormat.KindCounts).Offset;
            BinaryPrimitives.WriteUInt64LittleEndian(img.AsSpan((int)at), BinaryPrimitives.ReadUInt64LittleEndian(img.AsSpan((int)at)) - 1);
        });
        Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain("not to the 11 symbols"));
    }

    [Test]
    public void ACountOfZero_IsCorrupt_AtHeader()
    {
        byte[] image = Sealed(Mississippi, (img, o) =>
        {
            OracleSection counts = o.Get(SesaFormat.KindCounts);
            ulong first = BinaryPrimitives.ReadUInt64LittleEndian(img.AsSpan((int)counts.Offset));
            BinaryPrimitives.WriteUInt64LittleEndian(img.AsSpan((int)counts.Offset), 0);
            BinaryPrimitives.WriteUInt64LittleEndian(img.AsSpan((int)counts.Offset + 8), BinaryPrimitives.ReadUInt64LittleEndian(img.AsSpan((int)counts.Offset + 8)) + first);
        });
        Assert.That(Refusal(image, EsaValidationLevel.Header), Does.Contain("does not occur"));
        byte[] more = Sealed(Mississippi, (img, o) => BinaryPrimitives.WriteUInt64LittleEndian(img.AsSpan((int)o.Get(SesaFormat.KindCounts).Offset), 12));
        Assert.That(Refusal(more, EsaValidationLevel.Header), Does.Contain("sums to more than"));
    }

    [Test]
    public void TheWalk_ReportsNothing_OfAnIntervalWhoseFirstLIndexLiesOutsideIt()
    {
        // the cell 4 holds the first l-index of 'i' = [1, 4]; 9 is in the text but not in (1, 4], and the lcp at 9 is deeper than the root: only the rule on the first l-index stops the walk before that interval is reported
        var clean = new RecordingVisitor();
        using (EsaSuffixTree reference = SesaTestFile.FileTree(Mississippi)) reference.Traverse(clean);
        byte[] image = Sealed(Mississippi, (img, o) => SetI32(img, o.Get(SesaFormat.KindChild).Offset + (4 * 4), 9));
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image);
        var damaged = new RecordingVisitor();
        Exception? e = WithTimeout(() => Assert.Catch(() => tree.Traverse(damaged)), 10);
        Assert.That(e, Is.TypeOf<InvalidDataException>());
        Assert.That(damaged.Events.Count, Is.LessThan(clean.Events.Count));
        Assert.That(damaged.Events, Is.EqualTo(clean.Events.Take(damaged.Events.Count)), "what was reported before the refusal is what the clean walk reports first");
    }

    private sealed class RecordingVisitor : ISuffixTreeVisitor
    {
        internal List<string> Events { get; } = [];

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) => Events.Add($"node {startIndex} {endIndex} {leafCount} {childCount} {depth}");

        public void EnterBranch(int key) => Events.Add($"enter {key}");

        public void ExitBranch() => Events.Add("exit");
    }

    private sealed class CountingVisitor : ISuffixTreeVisitor
    {
        internal int Calls { get; private set; }

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) => Calls++;

        public void EnterBranch(int key) => Calls++;

        public void ExitBranch() => Calls++;
    }

    [Test]
    public void AFirstLIndexOutsideItsInterval_IsRefusedByTheDescent_WithTheTablesMessage()
    {
        // the interval [1, 4] of 'i': its first l-index is in the cell 4 and must be in (1, 4]; 10 is a rank of the text whose lcp (1) is a depth that the descent accepts for the child, so only the rule on the first l-index sees it when the descent reads it again for the next symbol
        byte[] image = Sealed(Mississippi, (img, o) => SetI32(img, o.Get(SesaFormat.KindChild).Offset + (4 * 4), 10));
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image);
        Exception? e = WithTimeout(() => Assert.Catch(() => tree.CountOccurrences("is")), 10);
        Assert.That(e, Is.TypeOf<InvalidDataException>());
        Assert.That(e!.Message, Does.Contain("inconsistent"));
    }

    [Test]
    public void AnEdgeLongerThanTheComparisonBuffer_IsComparedToItsLastSymbol()
    {
        string x = StringFamilies.RandomDna(700, 7);
        string text = x + x; // a repeat of 700 symbols: one edge of the tree holds nearly all of it
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        using EsaSuffixTree file = SesaTestFile.FileTree(text);
        Assert.That(file.CountOccurrences(x[..650]), Is.EqualTo(2));
        for (int at = 0; at < 650; at++)
        {
            char[] p = x[..650].ToCharArray();
            p[at] = p[at] == 'A' ? 'C' : 'A';
            string pattern = new(p);
            Assert.That(file.CountOccurrences(pattern), Is.EqualTo(memory.CountOccurrences(pattern)), $"symbol {at} of the pattern changed");
        }
    }
}
