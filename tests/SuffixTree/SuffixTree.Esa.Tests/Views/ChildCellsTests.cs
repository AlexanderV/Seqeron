using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.Fuzz;

namespace SuffixTree.Esa.Tests.Views;

/// <summary>
/// The child cells of work package S2.3 over a position reader and over a mapped section of any stored width: <see cref="ChildCells{TPos}"/> (the generic form) and <see cref="MappedChild"/> (the form that a container holds): the packed
/// codec over 4- and 5-byte cells, the tagged codec over 4-byte cells only, 64-bit ranks, the bounds of Checked mode, and the equality with the layouts of the in-memory index of S1 on every pair of ranks of small texts.
/// </summary>
[TestFixture]
public unsafe class ChildCellsTests
{
    private static readonly string[] Texts = ["banana", "mississippi", "aaaaaaaaaaaa", "abracadabra", "ACGTACGTTTGAC", "ab", "a", "abababababab"];

    /// <summary>The three forms that exist: 4-byte packed cells, 4-byte tagged cells, 5-byte packed cells (a tagged codec at width 5 does not exist, ADR-005 section 6.1).</summary>
    private static (int Width, bool Tagged) Form(string form) => form switch { "4" => (4, false), "4 tagged" => (4, true), _ => (5, false) };

    private static CompactLcp LcpOf(ViewImage image) =>
        CompactLcp.FromMappedWide(image.Kind, image.Count, image.LcpBasePointer, image.PairsPointer, image.Pairs.Count, image.PosWidth, isChecked: true);

    private static IEnumerable<(int Lb, int Rb, bool Last)> Intervals(int n)
    {
        for (int lb = 0; lb <= n; lb++)
        {
            for (int rb = 0; rb <= n; rb++)
            {
                yield return (lb, rb, true);
                yield return (lb, rb, false);
            }
        }
    }

    // ---------------------------------------------------------------------------------------------- equality with the layouts of S1

    [Test]
    public void TheGenericCells_EqualTheLayoutsOfTheInMemoryIndex_OnEveryPairOfRanks([ValueSource(nameof(Texts))] string text, [Values("4", "4 tagged", "5")] string form, [Values(true, false)] bool isChecked)
    {
        (int width, bool tagged) = Form(form);
        using ViewImage image = ViewImage.Of(text, width, LcpKind.U8, tagged);
        CompactLcp lcp = LcpOf(image);
        ChildTable expected = image.Case.Index.WithLayout(tagged ? ChildLayoutKind.Tagged : ChildLayoutKind.Packed).Child;
        int n = image.Case.N;
        if (width == 4 && isChecked) Compare(new ChildCells<Pos32<CheckedAccess>>(new Pos32<CheckedAccess>(image.ChildPointer, image.Count), lcp, tagged), expected, n);
        else if (width == 4) Compare(new ChildCells<Pos32<TrustedAccess>>(new Pos32<TrustedAccess>(image.ChildPointer, image.Count), lcp, tagged), expected, n);
        else if (isChecked) Compare(new ChildCells<Pos40<CheckedAccess>>(new Pos40<CheckedAccess>(image.ChildPointer, image.Count), lcp, false), expected, n);
        else Compare(new ChildCells<Pos40<TrustedAccess>>(new Pos40<TrustedAccess>(image.ChildPointer, image.Count), lcp, false), expected, n);
    }

    private static void Compare<TPos>(ChildCells<TPos> cells, ChildTable expected, int n)
        where TPos : struct, IPosReader
    {
        Assert.That(cells.Length, Is.EqualTo(n + 1L));
        foreach ((int lb, int rb, bool last) in Intervals(n))
        {
            Assert.That(cells.FirstLIndex(lb, rb, last), Is.EqualTo(expected.FirstLIndex(lb, rb, last)), $"FirstLIndex({lb}, {rb}, {last})");
            if (lb <= rb) Assert.That(cells.NextLIndex(lb, rb), Is.EqualTo(expected.NextLIndex(lb, rb)), $"NextLIndex({lb}, {rb})");
        }
    }

    [Test]
    public void TheMappedChild_AtEveryForm_EqualsTheLayoutsOfTheInMemoryIndex([ValueSource(nameof(Texts))] string text, [Values("4", "4 tagged", "5")] string form, [Values(true, false)] bool isChecked)
    {
        (int width, bool tagged) = Form(form);
        using ViewImage image = ViewImage.Of(text, width, LcpKind.U16, tagged);
        CompactLcp lcp = LcpOf(image);
        ChildTable expected = image.Case.Index.WithLayout(tagged ? ChildLayoutKind.Tagged : ChildLayoutKind.Packed).Child;
        MappedChild mapped = MappedChild.FromMappedWide(image.ChildPointer, image.Count, width, lcp, isChecked, tagged);
        Assert.That(mapped.Length, Is.EqualTo(image.Count));
        Assert.That(mapped.IsTagged, Is.EqualTo(tagged));
        int n = image.Case.N;
        foreach ((int lb, int rb, bool last) in Intervals(n))
        {
            Assert.That(mapped.FirstLIndex((long)lb, rb, last), Is.EqualTo(expected.FirstLIndex(lb, rb, last)), $"FirstLIndex({lb}, {rb}, {last})");
            Assert.That(mapped.FirstLIndex(lb, rb, last), Is.EqualTo(expected.FirstLIndex(lb, rb, last)), $"int FirstLIndex({lb}, {rb}, {last})");
            if (lb > rb) continue;
            Assert.That(mapped.NextLIndex((long)lb, rb), Is.EqualTo(expected.NextLIndex(lb, rb)), $"NextLIndex({lb}, {rb})");
            Assert.That(mapped.NextLIndex(lb, rb), Is.EqualTo(expected.NextLIndex(lb, rb)), $"int NextLIndex({lb}, {rb})");
        }
        ChildTable table = ChildTable.FromMapped(mapped);
        Assert.That(table.FirstLIndex(0L, n, true), Is.EqualTo(expected.FirstLIndex(0, n, true)));
        Assert.That(table.NextLIndex(0L, n), Is.EqualTo(expected.NextLIndex(0, n)));
        Assert.That(table.FirstLIndex(0, n, true), Is.EqualTo(expected.FirstLIndex(0, n, true)));
    }

    [Test]
    public void TheWideFormAtWidth4_AgreesWithTheLegacyMappedChild([ValueSource(nameof(Texts))] string text, [Values(true, false)] bool tagged)
    {
        using ViewImage image = ViewImage.Of(text, 4, LcpKind.U8, tagged);
        CompactLcp lcp = LcpOf(image);
        var legacy = new MappedChild((int*)image.ChildPointer, (int)image.Count, lcp, isChecked: true, tagged: tagged);
        MappedChild wide = MappedChild.FromMappedWide(image.ChildPointer, image.Count, 4, lcp, true, tagged);
        int n = image.Case.N;
        foreach ((int lb, int rb, bool last) in Intervals(n))
        {
            Assert.That(wide.FirstLIndex((long)lb, rb, last), Is.EqualTo(legacy.FirstLIndex(lb, rb, last)));
            if (lb <= rb) Assert.That(wide.NextLIndex((long)lb, rb), Is.EqualTo(legacy.NextLIndex(lb, rb)));
        }
        Assert.That(wide.PayloadBytes, Is.EqualTo(legacy.PayloadBytes));
        Assert.That((long)wide.Pointer, Is.EqualTo((long)legacy.Pointer));
    }

    // ---------------------------------------------------------------------------------------------- bounds, damage, shapes

    [Test]
    public void Checked_ACellOutsideTheSection_IsInvalidData([Values("4", "4 tagged", "5")] string form)
    {
        (int width, bool tagged) = Form(form);
        using ViewImage image = ViewImage.Of("mississippi", width, LcpKind.U8, tagged);
        CompactLcp lcp = LcpOf(image);
        long n = image.Count;
        MappedChild mapped = MappedChild.FromMappedWide(image.ChildPointer, n, width, lcp, isChecked: true, tagged);
        foreach (long cell in new[] { n, n + 1, -1, 1L << 33, long.MaxValue, long.MinValue })
        {
            Assert.That(() => mapped.FirstLIndex(cell, 0, true), Throws.TypeOf<InvalidDataException>(), $"lb {cell}");
            Assert.That(() => mapped.FirstLIndex(0, cell, false), Throws.TypeOf<InvalidDataException>(), $"rb {cell}");
        }
        Assert.That(() => mapped.NextLIndex(n, n + 2), Throws.TypeOf<InvalidDataException>().Or.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => mapped.FirstLIndex((int)n, 0, true), Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void ATaggedTableOverFiveByteCells_IsRefused([Values(true, false)] bool isChecked)
    {
        using ViewImage image = ViewImage.Of("banana", 5, LcpKind.U8);
        CompactLcp lcp = LcpOf(image);
        Assert.That(() => MappedChild.FromMappedWide(image.ChildPointer, image.Count, 5, lcp, isChecked, tagged: true), Throws.TypeOf<ArgumentException>());
        Assert.That(() => MappedChild.FromMappedWide(image.ChildPointer, image.Count, 3, lcp, isChecked, tagged: false), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => MappedChild.FromMappedWide(image.ChildPointer, -1, 4, lcp, isChecked, tagged: false), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => MappedChild.FromMappedWide(image.ChildPointer, image.Count, 4, null!, isChecked, tagged: false), Throws.TypeOf<ArgumentNullException>());
        using var four = ViewImage.Of("banana", 4, LcpKind.U8);
        Assert.That(() => new ChildCells<Pos40<CheckedAccess>>(new Pos40<CheckedAccess>(four.ChildPointer, 7), lcp, tagged: true), Throws.TypeOf<ArgumentException>());
        Assert.That(() => new ChildCells<PosInt>(new PosInt([1, 2]), lcp, tagged: true), Throws.Nothing, "a 4-byte reader may be tagged");
    }

    [Test]
    public void ThePackedValueOfADamagedCell_IsNotUsedAsARank_WhenItLiesOutsideTheText([Values(4, 5)] int width)
    {
        // a packed cell whose value is 2^33 (width 5) or 2^31 + 7 (width 4): NextLIndex follows it into the LCP array, which refuses the rank
        using ViewImage image = ViewImage.Of("mississippi", width, LcpKind.U8);
        CompactLcp lcp = LcpOf(image);
        ulong bad = width == 5 ? 1UL << 33 : (1UL << 31) + 7;
        ulong[] cells = [.. image.Cells];
        // the cell of rank 1: every rank k below rb whose next entry differs reads cell k; use the smallest such k
        int k = 1;
        while (k < image.Count - 2 && image.Lcp[k + 1] == image.Lcp[k]) k++;
        cells[k] = bad;
        using var guarded = new GuardedImage(ViewImage.Encode(cells, width, width == 5 ? 3 : 0), flush: true);
        MappedChild mapped = MappedChild.FromMappedWide(guarded.Data, image.Count, width, lcp, isChecked: true, tagged: false);
        Assert.That(mapped.FirstLIndex((long)k, 0, true), Is.EqualTo((long)bad), "a first l-index is returned as it is: the caller validates it");
        Assert.That(() => mapped.NextLIndex((long)k, image.Count - 1), Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void TheTaggedCell_KeepsTheThirtyBitsOfTheValue_AndTheKindOfNext()
    {
        // cell 0: next with value 2^30 - 1; cell 1: up with the same value; cell 2: down; cell 3: empty. Kind in bits 30-31.
        uint[] cells = [(1u << 30) | ((1u << 30) - 1), (2u << 30) | ((1u << 30) - 1), (3u << 30) | 12345, 0];
        using var guarded = new GuardedImage(ViewImage.Encode(cells.Select(c => (ulong)c), 4, 0), flush: true);
        LcpChildCase c = LcpChildCase.Of("banana");
        using ViewImage image = ViewImage.Of(c, 4, LcpKind.U8, tagged: true);
        CompactLcp lcp = LcpOf(image);
        MappedChild mapped = MappedChild.FromMappedWide(guarded.Data, 4, 4, lcp, isChecked: true, tagged: true);
        Assert.That(mapped.FirstLIndex(0L, 0, true), Is.EqualTo((1L << 30) - 1));
        Assert.That(mapped.FirstLIndex(0L, 1, false), Is.EqualTo((1L << 30) - 1));
        Assert.That(mapped.FirstLIndex(2L, 0, true), Is.EqualTo(12345));
        Assert.That(mapped.NextLIndex(0L, 3), Is.EqualTo((1L << 30) - 1));
        Assert.That(mapped.NextLIndex(1L, 3), Is.EqualTo(-1), "an up cell is not a next link");
        Assert.That(mapped.NextLIndex(2L, 3), Is.EqualTo(-1), "a down cell is not a next link");
        Assert.That(mapped.NextLIndex(3L, 3), Is.EqualTo(-1));
    }

    [Test]
    public void ThePackedNextLIndex_AtTheRightBound_IsMinusOne_WithoutReading([Values(4, 5)] int width)
    {
        using ViewImage image = ViewImage.Of("banana", width, LcpKind.U8);
        CompactLcp lcp = LcpOf(image);
        MappedChild mapped = MappedChild.FromMappedWide(image.ChildPointer, image.Count, width, lcp, true, false);
        Assert.That(mapped.NextLIndex(5L, 5), Is.EqualTo(-1));
        Assert.That(mapped.NextLIndex(6L, 5), Is.EqualTo(-1));
        Assert.That(mapped.NextLIndex(long.MaxValue, 5), Is.EqualTo(-1), "k >= rb answers before any read");
    }

    [Test]
    public void ThePayloadBytes_OfAWideTable_CountThePadding()
    {
        using ViewImage four = ViewImage.Of("banana", 4, LcpKind.U8);
        using ViewImage five = ViewImage.Of("banana", 5, LcpKind.U8);
        Assert.That(MappedChild.FromMappedWide(four.ChildPointer, four.Count, 4, LcpOf(four), true, false).PayloadBytes, Is.EqualTo(4 * four.Count));
        Assert.That(MappedChild.FromMappedWide(five.ChildPointer, five.Count, 5, LcpOf(five), true, false).PayloadBytes, Is.EqualTo((5 * five.Count) + 3));
    }

    [Test]
    public void AnIntegerMemberOfAWideTable_NeverNarrows()
    {
        // a width-5 table of 5 cells whose cell 0 is 2^33 + 1: the long member reads it, the 32-bit member says it does not fit
        ulong[] cells = [(1UL << 33) + 1, 4, 0, 0, 0];
        using var guarded = new GuardedImage(ViewImage.Encode(cells, 5, 3), flush: true);
        using ViewImage image = ViewImage.Of("abab", 5, LcpKind.U8);
        MappedChild mapped = MappedChild.FromMappedWide(guarded.Data, 5, 5, LcpOf(image), isChecked: true, tagged: false);
        Assert.That(mapped.FirstLIndex(0L, 0, true), Is.EqualTo((1L << 33) + 1));
        Assert.That(() => mapped.FirstLIndex(0, 0, true), Throws.TypeOf<InvalidDataException>());
        Assert.That(mapped.FirstLIndex(1, 1, true), Is.EqualTo(4));
        Assert.That(() => mapped.Pointer, Throws.TypeOf<NotSupportedException>());
    }

    [Test]
    public void TheManagedTable_TakesTheLongOverloads_AndRefusesWhatIsNotAnIntegerRank()
    {
        LcpChildCase c = LcpChildCase.Of("mississippi");
        foreach (ChildLayoutKind kind in new[] { ChildLayoutKind.Reference3, ChildLayoutKind.Packed, ChildLayoutKind.Tagged })
        {
            ChildTable table = c.Index.WithLayout(kind).Child;
            for (int lb = 0; lb <= c.N; lb++)
            {
                Assert.That(table.FirstLIndex((long)lb, c.N, true), Is.EqualTo(table.FirstLIndex(lb, c.N, true)));
                Assert.That(table.NextLIndex((long)lb, c.N), Is.EqualTo(table.NextLIndex(lb, c.N)));
            }
            Assert.That(() => table.FirstLIndex(1L << 33, 0, true), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => table.NextLIndex(0, 1L << 33), Throws.TypeOf<ArgumentOutOfRangeException>());
        }
    }
}
