using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.Views;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// The generator-backed views of work package S2.3 at small n (default set): the views over <see cref="PosGen{TWidth}"/> equal the views of the in-memory index of S1 on the closed forms of S2.1 (<c>a^n</c>, <c>(a^k b)*</c>), so that the
/// generator and the views are checked against each other and against S1 before the same code runs at n above 2^31 in <see cref="GeneratedViewLargeTests"/>. The oracle of every row is the S1 index (<see cref="LcpChildCase"/>), which the closed
/// forms were themselves checked against in S2.1.
/// </summary>
[TestFixture]
public class GeneratedViewSelfTests
{
    /// <summary>The shapes: a^n and (a^k b)* for a few k and n, up to a few thousand symbols (the S1 index builds in microseconds).</summary>
    internal static IEnumerable<ClosedFormText> Shapes()
    {
        foreach (int n in new[] { 1, 2, 3, 7, 40, 301, 700 }) yield return new PowerOfA(n);
        foreach (int k in new[] { 1, 2, 3, 5 })
        {
            foreach (int n in new[] { 4, 9, 33, 250, 600 }) yield return new PeriodicAkB(n, k);
        }
    }

    private static string Name(ClosedFormText f) => f.ToString();

    [Test]
    public void TheSuffixArrayView_OverTheGenerator_EqualsTheS1SuffixArray([ValueSource(nameof(Shapes))] ClosedFormText f)
    {
        LcpChildCase c = LcpChildCase.Of(f.Materialize());
        var view = new SaView<PosGen<W5>>(new PosGen<W5>(f, GenSection.Sa));
        Assert.That(view.Length, Is.EqualTo(c.N + 1L), Name(f));
        for (int r = 0; r <= c.N; r++)
        {
            Assert.That(view[r], Is.EqualTo(c.Sa[r]), $"{Name(f)} rank {r}");
            Assert.That(view.At(r), Is.EqualTo(c.Sa[r]), $"{Name(f)} At({r})");
        }
        Assert.That(view.ToArray(0, c.N + 1), Is.EqualTo(c.Sa.Select(v => (long)v)));
        var ints = new int[Math.Min(5, c.N + 1)];
        view.CopyTo(c.N + 1 - ints.Length, ints);
        Assert.That(ints, Is.EqualTo(c.Sa.TakeLast(ints.Length)));
        Assert.That(() => view.At(c.N + 1L), Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void TheLcpArray_OverTheGenerator_EqualsTheS1Lcp_InEveryKind([ValueSource(nameof(Shapes))] ClosedFormText f)
    {
        LcpChildCase c = LcpChildCase.Of(f.Materialize());
        int[] dense = c.Lcp.ToArray();
        long length = c.N + 1L;
        CompactLcp u64 = CompactLcp.FromReaders(LcpKind.U64, new PosGen<W8>(f, GenSection.LcpBase), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);
        CompactLcp u32 = CompactLcp.FromReaders(LcpKind.U32, new PosGen<W4>(f, GenSection.LcpBase), new PosGen<W4>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);
        Assert.That(u64.Length, Is.EqualTo(length));
        Assert.That(u32.Length, Is.EqualTo(length));
        for (int r = 0; r < length; r++)
        {
            Assert.That(u64[(long)r], Is.EqualTo(dense[r]), $"{Name(f)} U64 rank {r}");
            Assert.That(u32[(long)r], Is.EqualTo(dense[r]), $"{Name(f)} U32 rank {r}");
            Assert.That(u64[r], Is.EqualTo(dense[r]), $"{Name(f)} U64 int rank {r}");
        }
        Assert.That(() => u64[length], Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void TheOverflowListOfPowerOfA_IsSearchedByRank([Values(1, 255, 256, 300, 700, 5_000)] int n, [Values(255, 65_535)] int escape)
    {
        var f = new PowerOfA(n);
        LcpKind kind = escape == 255 ? LcpKind.U8 : LcpKind.U16;
        long pairs = Math.Max(0, n - escape);
        CompactLcp lcp = kind == LcpKind.U8
            ? CompactLcp.FromReaders(kind, new PosGen<W1>(f, GenSection.LcpBase, escape), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, escape), isChecked: true)
            : CompactLcp.FromReaders(kind, new PosGen<W2>(f, GenSection.LcpBase, escape), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, escape), isChecked: true);
        Assert.That(lcp.OverflowCount, Is.EqualTo(pairs));
        for (long r = 0; r <= n; r++) Assert.That(lcp[r], Is.EqualTo(f.Lcp(r)), $"a^{n} rank {r} escape {escape}");
    }

    [Test]
    public void TheChildCells_OverTheGenerator_EqualTheS1Layouts_OnEveryPairOfRanks([ValueSource(nameof(Shapes))] ClosedFormText f, [Values("packed", "tagged")] string codec)
    {
        bool tagged = codec == "tagged";
        LcpChildCase c = LcpChildCase.Of(f.Materialize());
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U64, new PosGen<W8>(f, GenSection.LcpBase), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);
        ChildTable expected = c.Index.WithLayout(tagged ? ChildLayoutKind.Tagged : ChildLayoutKind.Packed).Child;
        ChildCells<PosGen<W4>> cells = new(new PosGen<W4>(f, tagged ? GenSection.ChildTagged : GenSection.ChildPacked), lcp, tagged);
        Assert.That(cells.Length, Is.EqualTo(c.N + 1L));
        for (int lb = 0; lb <= c.N; lb++)
        {
            for (int rb = lb; rb <= c.N; rb++)
            {
                Assert.That(cells.FirstLIndex(lb, rb, true), Is.EqualTo(expected.FirstLIndex(lb, rb, true)), $"{Name(f)} FirstLIndex({lb}, {rb}, last)");
                Assert.That(cells.FirstLIndex(lb, rb, false), Is.EqualTo(expected.FirstLIndex(lb, rb, false)), $"{Name(f)} FirstLIndex({lb}, {rb}, not last)");
                Assert.That(cells.NextLIndex(lb, rb), Is.EqualTo(expected.NextLIndex(lb, rb)), $"{Name(f)} NextLIndex({lb}, {rb})");
            }
        }
    }

    [Test]
    public void TheTaggedCellOfTheClosedForm_IsTheTaggedCellOfS1([ValueSource(nameof(Shapes))] ClosedFormText f)
    {
        LcpChildCase c = LcpChildCase.Of(f.Materialize());
        uint[] tagged = c.Tagged;
        for (int r = 0; r <= c.N; r++) Assert.That(PosGen<W4>.Tagged(f, r), Is.EqualTo(tagged[r]), $"{Name(f)} cell {r}");
    }

    [Test]
    public void TheSyntheticLcp_IsALcpArrayForThePyramid()
    {
        var f = new PowerOfA(2_000);
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U32, new PosGen<W4>(f, GenSection.Synthetic), new PosGen<W4>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: false);
        using LcpPyramid pyramid = LcpPyramid.Build(lcp, 8);
        for (long r = 0; r <= 2_000; r += 13)
        {
            Assert.That(pyramid.ValueAt(r), Is.EqualTo((long)(PosGen<W4>.Hash(r) % 50)));
            long expected = r;
            while (expected >= 0 && (long)(PosGen<W4>.Hash(expected) % 50) >= 5) expected--;
            Assert.That(pyramid.PrevLess(r, 5), Is.EqualTo(expected));
        }
    }

    [Test]
    public void TheGeneratorRefusesAReadOutsideItsSection()
    {
        var f = new PowerOfA(10);
        var sa = new PosGen<W5>(f, GenSection.Sa);
        Assert.That(sa.Count, Is.EqualTo(11));
        Assert.That(sa.Get(0), Is.EqualTo(10UL));
        Assert.That(sa.Get(10), Is.EqualTo(0UL));
        Assert.That(() => sa.Get(11), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => sa.Get(-1), Throws.TypeOf<InvalidDataException>());
        Assert.That(PosGen<W5>.ElementBytes, Is.EqualTo(5));
        Assert.That(PosGen<W1>.ElementBytes, Is.EqualTo(1));
        Assert.That(new PosGen<W8>(f, GenSection.LcpPairsOfPowerOfA, 4).Count, Is.EqualTo(12));
    }
}
