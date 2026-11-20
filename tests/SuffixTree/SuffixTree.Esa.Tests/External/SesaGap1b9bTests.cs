using System.Buffers.Binary;
using System.Text;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The test gaps of the orchestrator's verifications of 1b.7e and 1b.7f that work package 1b.9b closes (each one is damage or a hand-over that only one direction of was tested): <c>B102</c> the guard option of the builder reaches
/// every pass, including the BWT pass; <c>B118</c> a BUILD_INFO that would overrun its room is refused; <c>VE7</c> a non-zero byte in the first eight bytes of a padding gap between sections is refused at <c>Structure</c>;
/// <c>VE9</c> a PSI primary row stored below the real one is refused; <c>VE12</c> an <c>LRS_DEPTH</c> stored above the real one is refused at <c>Structure</c>.
/// </summary>
[TestFixture]
public class SesaGap1b9bTests
{
    // ---------------------------------------------------------------- B102: the guards of the builder reach every pass

    [Test]
    public void TheGuardOption_ReachesThePassOfTheSuffixArray_TheLcpAndChildPass_TheBwtPass_AndThePsiPass()
    {
        using var dir = new TempDir();
        var points = new List<(string Point, long Count)>();
        SesaBuildOptions guarded = SesaRun.Tiny(3, 2) with { Hook = p => points.Add((p, NativeArena.GuardedCreated)) };
        SesaRun.Build(dir, "mississippimississippi" + "abracadabra" + "banana", guarded, "-g");
        long At(string point) => points.Single(x => x.Point == point).Count;
        Assert.That(At(SesaHooks.SaBuilt), Is.GreaterThan(At(SesaHooks.Packed)), "the suffix-array pass creates guarded arenas");
        Assert.That(At(SesaHooks.LcpChildBuilt), Is.GreaterThan(At(SesaHooks.SaBuilt)), "the LCP and child pass creates guarded arenas");
        Assert.That(At(SesaHooks.BwtBuilt), Is.GreaterThan(At(SesaHooks.LcpChildBuilt)), "the BWT pass creates guarded arenas (B102)");
        Assert.That(At(SesaHooks.PsiBuilt), Is.GreaterThan(At(SesaHooks.BwtBuilt)), "the PSI pass creates guarded arenas");
    }

    [Test]
    public void WithoutTheGuardOption_NoPassCreatesAGuardedArena()
    {
        // the control of the test above: the counter moves because of the option and of nothing else
        using var dir = new TempDir();
        long before = NativeArena.GuardedCreated;
        SesaBuildOptions plain = SesaRun.Tiny(3, 2) with { ArenaGuards = null };
        SesaRun.Build(dir, "mississippimississippi" + "abracadabra" + "banana", plain, "-p");
        Assert.That(NativeArena.GuardedCreated, Is.EqualTo(before));
    }

    // ---------------------------------------------------------------- B118: BUILD_INFO may not overrun its room

    [Test]
    public void ABuildInfoThatWouldOverrunItsRoom_IsRefused_OneByteOverTheCapacityIsEnough()
    {
        const long capacity = SesaLayout.BuildInfoCapacity;
        Assert.That(capacity, Is.EqualTo(4096));
        SesaBuilder.CheckBuildInfoFits(0, capacity);
        SesaBuilder.CheckBuildInfoFits((int)capacity - 1, capacity);
        SesaBuilder.CheckBuildInfoFits((int)capacity, capacity);
        foreach (int over in new[] { 1, 2, 100, 4096 })
        {
            int length = (int)capacity + over;
            Assert.That(() => SesaBuilder.CheckBuildInfoFits(length, capacity), Throws.TypeOf<InvalidOperationException>().With.Message.Contains("BUILD_INFO"), $"{length} bytes in a room of {capacity}");
        }
    }

    [Test]
    public void TheBuildInfoOfABuild_FitsItsRoom_WithMuchToSpare()
    {
        using var dir = new TempDir();
        SesaRun.Outcome o = SesaRun.Build(dir, "mississippi", new SesaBuildOptions { Fsync = false });
        OracleSection info = o.Oracle.Get(SesaFormat.KindBuildInfo);
        Assert.That(info.Length, Is.LessThan(SesaLayout.BuildInfoCapacity / 4), "the text is a few hundred bytes in a room of 4,096");
        Assert.That(Encoding.UTF8.GetString(o.Oracle.Payload(info)), Does.StartWith("format=SESA 1.0\n"));
    }

    // ---------------------------------------------------------------- VE7: the first eight bytes of a padding gap

    private static IEnumerable<(string Label, byte[] Image)> Containers()
    {
        string random = StringFamilies.Random(1_500, 4, 21);
        yield return ("DNA 150, packed child, PSI", SesaTestFile.Image(LcpChildCase.Of(StringFamilies.Random(150, 4, 2))));
        yield return ("DNA 150, tagged child, no PSI", SesaTestFile.Image(LcpChildCase.Of(StringFamilies.Random(150, 4, 3)), psi: false, codec: ChildCodec.Tagged));
        yield return ("U8 with overflow, tagged child, PSI", SesaTestFile.Image(LcpChildCase.Of(random + random[..300]), codec: ChildCodec.Tagged));
        yield return ("mississippi, packed child, PSI", SesaTestFile.Image(LcpChildCase.Of("mississippi")));
    }

    [Test]
    public void ANonZeroByteInTheFirstEightBytesOfAGap_IsRefusedAtStructure_ForEveryGapOfEveryContainer()
    {
        using var dir = new TempDir();
        int tested = 0, gaps = 0;
        foreach ((string label, byte[] image) in Containers())
        {
            SesaOracle o = SesaOracle.FromBytes(image);
            long end = SesaFormat.PageBytes;
            foreach (OracleSection s in o.Sections.OrderBy(x => x.Offset).ThenBy(x => x.Length))
            {
                if (s.Offset > end)
                {
                    gaps++;
                    for (long b = end; b < Math.Min(end + 8, s.Offset); b++)
                    {
                        byte[] copy = (byte[])image.Clone();
                        copy[b] = 0x5A;
                        long at = b;
                        Assert.That(() => SesaTestFile.OpenImage(dir, copy, EsaValidationLevel.Structure, name: "g.sesa").Dispose(),
                            Throws.TypeOf<InvalidDataException>().With.Message.Contains("are not zero"), $"{label}: gap byte {at} (the gap starts at {end})");
                        tested++;
                    }
                    // the byte behind the eighth is also refused (the whole gap is checked)
                    if (s.Offset - end > 8)
                    {
                        byte[] copy = (byte[])image.Clone();
                        copy[s.Offset - 1] = 0x5A;
                        Assert.That(() => SesaTestFile.OpenImage(dir, copy, EsaValidationLevel.Structure, name: "g.sesa").Dispose(), Throws.TypeOf<InvalidDataException>().With.Message.Contains("are not zero"), $"{label}: the last byte of the gap before {s.Offset}");
                    }
                }
                end = Math.Max(end, s.Offset + s.Length);
            }
        }
        Assert.That(gaps, Is.GreaterThanOrEqualTo(10), "the containers have gaps");
        Assert.That(tested, Is.GreaterThanOrEqualTo(60));
    }

    [Test]
    public void TheCleanContainers_OfTheGapTest_OpenAtStructure()
    {
        using var dir = new TempDir();
        foreach ((string label, byte[] image) in Containers())
        {
            using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Structure, sha: true, name: "c.sesa");
            Assert.That(tree.Text.Length, Is.GreaterThan(0), label);
        }
    }

    // ---------------------------------------------------------------- VE9: the PSI primary row, in both directions

    private static void AssertOnlyStructureSees(string name, byte[] image)
    {
        using var dir = new TempDir();
        string path = dir.File("s.sesa");
        File.WriteAllBytes(path, image);
        using (EsaSuffixTree header = SesaTestFile.Open(path, EsaValidationLevel.Header)) Assert.That(header, Is.Not.Null, name);
        using (EsaSuffixTree sections = SesaTestFile.Open(path, EsaValidationLevel.Sections)) Assert.That(sections, Is.Not.Null, name);
        Assert.That(() => SesaTestFile.Open(path, EsaValidationLevel.Structure, EsaAccessMode.Checked, sha: true).Dispose(), Throws.TypeOf<InvalidDataException>(), $"{name}: Structure");
    }

    [TestCase("GATTACAGATTACAGGTA")]
    [TestCase("mississippi")]
    [TestCase("aaaaaaaaaaaa")]
    public void ThePsiPrimaryRow_StoredAtAnyRowButTheRealOne_IsRefusedAtStructure(string text)
    {
        LcpChildCase c = LcpChildCase.Of(text);
        int real = Array.IndexOf(c.Sa, 0);
        Assert.That(real, Is.GreaterThan(0));
        int tested = 0;
        for (long stored = -2; stored <= c.N + 2; stored++)
        {
            if (stored == real) continue;
            byte[] image = SesaTestFile.Image(c);
            OracleSection psi = SesaTestFile.Section(image, SesaFormat.KindPsi);
            BinaryPrimitives.WriteInt64LittleEndian(image.AsSpan(128 + (64 * psi.Index) + 40), stored);
            SesaTestFile.Reseal(image);
            long value = stored;
            if (stored is < 0 or > int.MaxValue || stored >= c.N + 1)
            {
                // outside 0..n the Header level may already refuse the table entry; Structure must in any case
                using var dir = new TempDir();
                Assert.That(() => SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Structure, sha: true).Dispose(), Throws.TypeOf<InvalidDataException>(), $"'{text}': primary row {value}");
            }
            else
            {
                AssertOnlyStructureSees($"'{text}': primary row {value} instead of {real}", image);
            }
            tested++;
        }
        Assert.That(tested, Is.GreaterThan(c.N));
    }

    // ---------------------------------------------------------------- VE12: LRS_DEPTH, in both directions

    [TestCase("GATTACAGATTACAGGTAGATTACA")]
    [TestCase("mississippi")]
    [TestCase("abracadabra")]
    public void TheLrsDepthStoredAboveOrBelowTheRealOne_IsRefusedAtStructure(string text)
    {
        LcpChildCase c = LcpChildCase.Of(text);
        long real = c.Statistics.LrsDepth;
        Assert.That(real, Is.GreaterThan(1));
        foreach (long delta in new long[] { 1, 2, 3, 7, 100, -1, -2 })
        {
            byte[] image = SesaTestFile.Image(c);
            BinaryPrimitives.WriteInt64LittleEndian(image.AsSpan(56), real + delta);
            SesaTestFile.Reseal(image);
            using var dir = new TempDir();
            string path = dir.File("l.sesa");
            File.WriteAllBytes(path, image);
            Assert.That(() => SesaTestFile.Open(path, EsaValidationLevel.Structure, EsaAccessMode.Checked, sha: true).Dispose(), Throws.TypeOf<InvalidDataException>(), $"'{text}': LRS_DEPTH {real} + {delta}");
        }
    }

    [Test]
    public void TheLrsDepthStoredAboveTheRealOne_IsRefusedByTheComparisonAtStructure_NotByTheHeaderRangeChecks()
    {
        // one above the real depth is inside every range the Header level can check (at most n), so only the comparison with the lcp-interval tree can see it
        string text = "GATTACAGATTACAGGTAGATTACA";
        LcpChildCase c = LcpChildCase.Of(text);
        byte[] image = SesaTestFile.Image(c);
        BinaryPrimitives.WriteInt64LittleEndian(image.AsSpan(56), c.Statistics.LrsDepth + 1);
        SesaTestFile.Reseal(image);
        using var dir = new TempDir();
        string path = dir.File("l.sesa");
        File.WriteAllBytes(path, image);
        using (EsaSuffixTree header = SesaTestFile.Open(path, EsaValidationLevel.Header)) Assert.That(header, Is.Not.Null);
        InvalidDataException? ex = Assert.Throws<InvalidDataException>(() => SesaTestFile.Open(path, EsaValidationLevel.Structure, EsaAccessMode.Checked, sha: true).Dispose());
        Assert.That(ex!.Message, Does.Contain("LRS_DEPTH"));
    }
}
