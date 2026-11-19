using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The placement of the sections (ADR-004 section 4.3): sizes known after PH0 first, the LCP behind them, the alignment rule of ADR-003, the length at which the file is created.</summary>
[TestFixture]
public sealed class SesaLayoutTests
{
    private static PackedTextInfo Info(long n, int sigma, int bits) => new(n, bits, [.. Enumerable.Range(0, sigma).Select(i => (byte)i)], [.. Enumerable.Range(0, sigma).Select(i => i == 0 ? n - Math.Max(0, sigma - 1) : 1L)], new string('0', 64));

    private static long Align(long offset, long length)
    {
        long a = length >= (1 << 20) ? 4096 : 64;
        return (offset + a - 1) / a * a;
    }

    public static IEnumerable<TestCaseData> Grid()
    {
        foreach (long n in new long[] { 0, 1, 2, 7, 100, 4095, 1_000_000, 100_000_000, 3_000_000_000 })
        {
            foreach ((int sigma, int bits) in new[] { (1, 1), (4, 2), (5, 3), (20, 5), (256, 8) })
            {
                if (n < sigma && n != 0) continue;
                if (n == 0 && sigma != 1) continue;
                foreach (bool psi in new[] { true, false }) yield return new TestCaseData(n, sigma, bits, psi);
            }
        }
    }

    [TestCaseSource(nameof(Grid))]
    public void ThePlan_PlacesEverySectionAtTheAlignmentOfADR003_InTheOrderOfADR004(long n, int sigma, int bits, bool psi)
    {
        PackedTextInfo info = Info(n, sigma, bits);
        SesaPlan p = SesaLayout.Plan(info, psi);
        long textLength = ((n * bits) + 7) / 8 + 64;
        Assert.That((p.N, p.Bits, p.Sigma, p.Psi), Is.EqualTo((n, bits, sigma, psi)));
        Assert.That(p.TextOffset, Is.EqualTo(4096));
        Assert.That(p.TextLength, Is.EqualTo(textLength));
        Assert.That(p.AlphabetOffset, Is.EqualTo(Align(p.TextOffset + textLength, 2L * sigma)));
        Assert.That(p.AlphabetLength, Is.EqualTo(2L * sigma));
        Assert.That(p.CountsOffset, Is.EqualTo(Align(p.AlphabetOffset + p.AlphabetLength, 8L * sigma)));
        Assert.That(p.CountsLength, Is.EqualTo(8L * sigma));
        Assert.That(p.SaOffset, Is.EqualTo(Align(p.CountsOffset + p.CountsLength, 4 * (n + 1))));
        Assert.That(p.SaLength, Is.EqualTo(4 * (n + 1)));
        Assert.That(p.ChildOffset, Is.EqualTo(Align(p.SaOffset + p.SaLength, 4 * (n + 1))));
        Assert.That(p.ChildLength, Is.EqualTo(4 * (n + 1)));
        long end = p.ChildOffset + p.ChildLength;
        if (psi)
        {
            Assert.That(p.PsiOffset, Is.EqualTo(Align(end, 4 * (n + 1))));
            Assert.That(p.PsiLength, Is.EqualTo(4 * (n + 1)));
            end = p.PsiOffset + p.PsiLength;
        }
        else
        {
            Assert.That((p.PsiOffset, p.PsiLength), Is.EqualTo((0L, 0L)));
        }
        Assert.That(p.BuildInfoCapacity, Is.EqualTo(SesaLayout.BuildInfoCapacity));
        Assert.That(p.BuildInfoOffset, Is.EqualTo(Align(end, SesaLayout.BuildInfoCapacity)));
        Assert.That(p.LcpOffset, Is.EqualTo(Align(p.BuildInfoOffset + p.BuildInfoCapacity, n + 1)));
        Assert.That(p.ReservedLength, Is.EqualTo(p.LcpOffset + n + 1), "the file is created with the LCP at its narrowest width");
        foreach (long offset in new[] { p.TextOffset, p.AlphabetOffset, p.CountsOffset, p.SaOffset, p.ChildOffset, p.BuildInfoOffset, p.LcpOffset }) Assert.That(offset % 64, Is.Zero);
        foreach ((long offset, long length) in new[] { (p.TextOffset, p.TextLength), (p.SaOffset, p.SaLength), (p.ChildOffset, p.ChildLength), (p.PsiOffset, p.PsiLength), (p.LcpOffset, n + 1) })
        {
            if (length >= (1 << 20)) Assert.That(offset % 4096, Is.Zero, $"offset {offset} of a section of {length} bytes");
        }
    }

    [Test]
    public void TheSmallestText_FitsInOneHeaderPageAndAFewHundredBytes()
    {
        SesaPlan p = SesaLayout.Plan(new PackedTextInfo(0, 1, [], [], new string('0', 64)), false);
        Assert.That((p.TextOffset, p.TextLength, p.AlphabetLength, p.CountsLength, p.SaLength, p.ChildLength), Is.EqualTo((4096L, 64L, 0L, 0L, 4L, 4L)));
        Assert.That(p.ReservedLength, Is.LessThan(4096 + 64 + 64 + 64 + 64 + 4096 + 64 + 64));
    }

    [TestCase(0, 0L)]
    [TestCase(0, 1L)]
    [TestCase(0, 12345L)]
    [TestCase(1, 0L)]
    [TestCase(1, 77L)]
    [TestCase(2, 0L)]
    public void TheTail_FollowsTheLcpAtItsWidth_AndTheOverflowListBehindIt(int kindId, long overflow)
    {
        LcpKind kind = kindId switch { 0 => LcpKind.U8, 1 => LcpKind.U16, _ => LcpKind.U32 };
        foreach (long n in new long[] { 0, 5, 1_000_000, 3_000_000_000 })
        {
            if (overflow > n + 1) continue;
            SesaPlan p = SesaLayout.Plan(Info(n, 4, 2), true);
            SesaTail t = SesaLayout.PlaceTail(p, kind, overflow);
            int w = kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 4 };
            Assert.That((t.Kind, t.LcpWidth, t.LcpLength, t.OverflowCount), Is.EqualTo((kind, w, (n + 1) * w, overflow)));
            long lcpEnd = p.LcpOffset + t.LcpLength;
            if (overflow == 0)
            {
                Assert.That((t.OverflowOffset, t.OverflowLength), Is.EqualTo((0L, 0L)));
                Assert.That(t.FinalLength, Is.EqualTo(lcpEnd));
            }
            else
            {
                Assert.That(t.OverflowLength, Is.EqualTo(8 * overflow));
                Assert.That(t.OverflowOffset, Is.EqualTo(Align(lcpEnd, 8 * overflow)));
                Assert.That(t.FinalLength, Is.EqualTo(t.OverflowOffset + t.OverflowLength));
            }
            Assert.That(t.FinalLength, Is.GreaterThanOrEqualTo(p.ReservedLength), "no text can undercut the reserved length");
            if (kind == LcpKind.U8 && overflow == 0) Assert.That(t.FinalLength, Is.EqualTo(p.ReservedLength), "the common case does not grow");
        }
    }

    [Test]
    public void AnOverflowListOfAnLcpWithoutEscapes_IsRefused()
    {
        SesaPlan p = SesaLayout.Plan(Info(10, 4, 2), false);
        Assert.Throws<ArgumentException>(() => SesaLayout.PlaceTail(p, LcpKind.U32, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => SesaLayout.PlaceTail(p, LcpKind.U8, -1));
    }

    [Test]
    public void ThePlan_NeverPutsTwoSectionsOnTheSameBytes()
    {
        foreach (bool psi in new[] { true, false })
        {
            SesaPlan p = SesaLayout.Plan(Info(1_000_000, 4, 2), psi);
            (long Offset, long Length)[] parts = [(p.TextOffset, p.TextLength), (p.AlphabetOffset, p.AlphabetLength), (p.CountsOffset, p.CountsLength), (p.SaOffset, p.SaLength), (p.ChildOffset, p.ChildLength), (p.PsiOffset, p.PsiLength), (p.BuildInfoOffset, p.BuildInfoCapacity), (p.LcpOffset, 1_000_001)];
            long end = 4096;
            foreach ((long offset, long length) in parts.Where(x => x.Length > 0))
            {
                Assert.That(offset, Is.GreaterThanOrEqualTo(end));
                end = offset + length;
            }
        }
    }
}
