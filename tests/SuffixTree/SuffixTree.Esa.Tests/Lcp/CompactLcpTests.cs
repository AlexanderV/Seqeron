using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Lcp;

/// <summary>1a.4: compact LCP storage (U8 + overflow, U16 + overflow, U32), the selection rule, the cursor reader.</summary>
[TestFixture]
public class CompactLcpTests
{
    private static readonly int[] BoundaryValues =
        { 0, 1, 2, 253, 254, 255, 256, 257, 65533, 65534, 65535, 65536, 65537, 1_000_000, int.MaxValue - 1, int.MaxValue };

    private static void AssertRepresents(CompactLcp lcp, int[] expected, LcpKind kind, string context)
    {
        Assert.That(lcp.Kind, Is.EqualTo(kind), context);
        Assert.That(lcp.Length, Is.EqualTo(expected.Length), context);
        Assert.That(lcp.ToArray(), Is.EqualTo(expected), context);
        for (int i = 0; i < expected.Length; i++)
            if (lcp[i] != expected[i]) Assert.Fail($"{context}: random access at {i}: {lcp[i]} != {expected[i]}");
        var cursor = lcp.CreateCursor();
        for (int i = 0; i < expected.Length; i++)
            if (cursor.Get(i) != expected[i]) Assert.Fail($"{context}: cursor at {i}: {cursor.Get(i)} != {expected[i]}");
        long escapes = kind switch
        {
            LcpKind.U8 => expected.Count(v => v >= 255),
            LcpKind.U16 => expected.Count(v => v >= 65535),
            _ => 0,
        };
        Assert.That(lcp.OverflowCount, Is.EqualTo(escapes), context);
    }

    [TestCase(0, TestName = "ForcedKind_RoundTripsBoundaryValues_U8")]
    [TestCase(1, TestName = "ForcedKind_RoundTripsBoundaryValues_U16")]
    [TestCase(2, TestName = "ForcedKind_RoundTripsBoundaryValues_U32")]
    public void ForcedKind_RoundTripsBoundaryValues(int kindValue)
    {
        var kind = (LcpKind)kindValue; // 0 = U8, 1 = U16, 2 = U32 (the enum is internal, NUnit needs a public signature)
        // every boundary value alone, in pairs, and shuffled with many repetitions
        foreach (int v in BoundaryValues)
            AssertRepresents(CompactLcp.FromDense(new[] { 0, v, 0 }, kind), new[] { 0, v, 0 }, kind, $"single {v}");

        var all = BoundaryValues.Concat(BoundaryValues.Reverse()).ToArray();
        AssertRepresents(CompactLcp.FromDense(all, kind), all, kind, "all boundary values");

        var rng = new SplitMix64(91);
        var shuffled = new int[3000];
        for (int i = 0; i < shuffled.Length; i++)
            shuffled[i] = rng.NextInt(3) == 0 ? BoundaryValues[rng.NextInt(BoundaryValues.Length)] : rng.NextInt(300);
        AssertRepresents(CompactLcp.FromDense(shuffled, kind), shuffled, kind, "shuffled");
    }

    [TestCase(0, TestName = "ForcedKind_EmptyAndTinyArrays_U8")]
    [TestCase(1, TestName = "ForcedKind_EmptyAndTinyArrays_U16")]
    [TestCase(2, TestName = "ForcedKind_EmptyAndTinyArrays_U32")]
    public void ForcedKind_EmptyAndTinyArrays(int kindValue)
    {
        var kind = (LcpKind)kindValue;
        AssertRepresents(CompactLcp.FromDense(new[] { 0 }, kind), new[] { 0 }, kind, "n = 0");
        AssertRepresents(CompactLcp.FromDense(new[] { 0, 0 }, kind), new[] { 0, 0 }, kind, "n = 1");
        Assert.Throws<ArgumentException>(() => CompactLcp.FromDense(ReadOnlySpan<int>.Empty, kind));
        Assert.Throws<ArgumentOutOfRangeException>(() => CompactLcp.FromDense(new[] { 0, -1 }, kind));
    }

    [Test]
    public void Cursor_ReadsForwardInAmortisedConstantTime_AndAlsoBackwardsAndRandomly()
    {
        var rng = new SplitMix64(5);
        var values = new int[5000];
        for (int i = 0; i < values.Length; i++) values[i] = rng.NextInt(4) == 0 ? 255 + rng.NextInt(100_000) : rng.NextInt(255);
        values[0] = 0;
        foreach (var kind in new[] { LcpKind.U8, LcpKind.U16, LcpKind.U32 })
        {
            var lcp = CompactLcp.FromDense(values, kind);
            var backwards = lcp.CreateCursor();
            for (int i = values.Length - 1; i >= 0; i--) Assert.That(backwards.Get(i), Is.EqualTo(values[i]), $"{kind} backwards {i}");
            var random = lcp.CreateCursor();
            for (int k = 0; k < 20_000; k++)
            {
                int i = rng.NextInt(values.Length);
                if (random.Get(i) != values[i]) Assert.Fail($"{kind} random access {i}");
            }
            var again = lcp.CreateCursor();
            Assert.That(again.Get(100), Is.EqualTo(values[100]));
            Assert.That(again.Get(100), Is.EqualTo(values[100]));
            Assert.Throws<ArgumentOutOfRangeException>(() => again.Get(values.Length));
            Assert.Throws<ArgumentOutOfRangeException>(() => _ = lcp[-1]);
            Assert.Throws<ArgumentOutOfRangeException>(() => _ = lcp[values.Length]);
        }
    }

    [Test]
    public void SelectKind_U8_WhenTheOverflowListIsAtMostHalfTheBaseArray()
    {
        // ADR-002 section 4: U8 + overflow iff escapes * 8 <= n / 2, i.e. 16 * escapes <= n
        Assert.That(CompactLcp.SelectKind(0, 0, 0), Is.EqualTo(LcpKind.U8));
        Assert.That(CompactLcp.SelectKind(1600, 100, 0), Is.EqualTo(LcpKind.U8));
        Assert.That(CompactLcp.SelectKind(1600, 101, 0), Is.EqualTo(LcpKind.U16));
        Assert.That(CompactLcp.SelectKind(1599, 100, 0), Is.EqualTo(LcpKind.U16));
        Assert.That(CompactLcp.SelectKind(1601, 100, 0), Is.EqualTo(LcpKind.U8));
        Assert.That(CompactLcp.SelectKind(15, 1, 0), Is.EqualTo(LcpKind.U16));
        Assert.That(CompactLcp.SelectKind(16, 1, 0), Is.EqualTo(LcpKind.U8));
        Assert.That(CompactLcp.SelectKind(int.MaxValue - 1, int.MaxValue / 16, 0), Is.EqualTo(LcpKind.U8));
    }

    [Test]
    public void SelectKind_U16_WhenItsOverflowListIsAtMostHalfTheBaseArray_ElseU32()
    {
        // the same shape one level up (base array 2n bytes, overflow entries 8 bytes): 8 * escapes16 <= n
        Assert.That(CompactLcp.SelectKind(800, 800, 100), Is.EqualTo(LcpKind.U16));
        Assert.That(CompactLcp.SelectKind(800, 800, 101), Is.EqualTo(LcpKind.U32));
        Assert.That(CompactLcp.SelectKind(799, 799, 100), Is.EqualTo(LcpKind.U32));
        Assert.That(CompactLcp.SelectKind(70_000, 69_745, 4465), Is.EqualTo(LcpKind.U16)); // a^70000
        Assert.That(CompactLcp.SelectKind(100, 100, 0), Is.EqualTo(LcpKind.U16));
        Assert.Throws<ArgumentException>(() => CompactLcp.SelectKind(100, 5, 6)); // values >= 65535 are also >= 255
        Assert.Throws<ArgumentOutOfRangeException>(() => CompactLcp.SelectKind(-1, 0, 0));
    }

    [Test]
    public void FromDense_AppliesTheSelectionRule_AtTheBoundaries()
    {
        // n = 1600 (1601 entries): 100 values of 255 -> U8 with 100 escapes; 101 -> U16 (no value >= 65535)
        foreach (var (escapes, kind) in new[] { (0, LcpKind.U8), (100, LcpKind.U8), (101, LcpKind.U16), (1600, LcpKind.U16) })
        {
            var values = new int[1601];
            for (int i = 1; i <= escapes; i++) values[i] = 255;
            AssertRepresents(CompactLcp.FromDense(values), values, kind, $"{escapes} escapes");
        }

        // values of 254 are not escapes
        var near = Enumerable.Repeat(254, 1601).ToArray();
        near[0] = 0;
        AssertRepresents(CompactLcp.FromDense(near), near, LcpKind.U8, "all 254");

        // n = 800: 100 values >= 65535 -> U16 (U8 impossible: 100 * 16 > 800); 101 -> U32
        foreach (var (escapes, kind) in new[] { (100, LcpKind.U16), (101, LcpKind.U32) })
        {
            var values = new int[801];
            for (int i = 1; i <= escapes; i++) values[i] = 65535;
            AssertRepresents(CompactLcp.FromDense(values), values, kind, $"{escapes} escapes >= 65535");
        }
    }

    [Test]
    public void BytesPerEntry_MatchTheKind()
    {
        Assert.That(CompactLcp.FromDense(new int[101], LcpKind.U8).PayloadBytes, Is.EqualTo(101));
        Assert.That(CompactLcp.FromDense(new int[101], LcpKind.U16).PayloadBytes, Is.EqualTo(202));
        Assert.That(CompactLcp.FromDense(new int[101], LcpKind.U32).PayloadBytes, Is.EqualTo(404));
        var withOverflow = new int[101];
        withOverflow[5] = 300;
        withOverflow[6] = 301;
        Assert.That(CompactLcp.FromDense(withOverflow, LcpKind.U8).PayloadBytes, Is.EqualTo(101 + (2 * 8)));
    }
}
