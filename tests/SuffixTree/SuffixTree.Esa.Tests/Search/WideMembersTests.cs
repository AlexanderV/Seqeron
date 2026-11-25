namespace SuffixTree.Esa.Tests.Search;

/// <summary>
/// Small members of the long search side that the audit of work package S2.4 (round 1) found without a test of their own: <see cref="EsaNode64"/>, the 64-bit <c>TryMatch</c> of an index in memory, the narrowing of the count of distinct
/// substrings to 64 bits.
/// </summary>
[TestFixture]
public class WideMembersTests
{
    [Test]
    public void EsaNode64_ShapeMembers_AreTheThoseOfTheThirtyTwoBitNode()
    {
        var inner = new EsaNode64(3, 9, 7, 2, 40, true);
        Assert.That(inner.IsLeaf, Is.False);
        Assert.That(inner.IsNull, Is.False);
        Assert.That(inner.EdgeLength, Is.EqualTo(5));
        Assert.That(inner.RankCount, Is.EqualTo(7));
        var leaf = new EsaNode64(5, 5, 12, 4, 1, false);
        Assert.That(leaf.IsLeaf, Is.True);
        Assert.That(leaf.RankCount, Is.EqualTo(1));
        Assert.That(leaf.EdgeLength, Is.EqualTo(8));
        Assert.That(new EsaNode64(4, 5, 12, 4, 1, false).IsLeaf, Is.False, "an interval of two ranks is not a leaf");
        Assert.That(new EsaNode64(4, 4, 12, 4, 1, false).IsLeaf, Is.True);
        Assert.That(EsaNode64.Null.IsNull, Is.True);
        Assert.That(new EsaNode64(0, 0, 0, 0, 0, true).IsNull, Is.False, "rank 0 is a rank");
        Assert.That(new EsaNode64(-1, 3, 0, 0, 0, true).IsNull, Is.True, "a negative left bound is the null node");
    }

    [Test]
    public void EsaNode64_RoundTripsTheThirtyTwoBitNode_WithEveryFieldAndTheLastChildFlag([Values(false, true)] bool lastChild)
    {
        var node = new EsaNode(3, 9, 7, 2, 40, lastChild);
        EsaNode64 wide = EsaNode64.Of(node);
        Assert.That((wide.Lb, wide.Rb, wide.Depth, wide.ParentDepth, wide.Pos, wide.LastChild), Is.EqualTo((3L, 9L, 7L, 2L, 40L, lastChild)));
        Assert.That(wide.ToNode("the test"), Is.EqualTo(node));
    }

    [Test]
    public void EsaNode64_ToNode_RefusesEveryFieldAboveTheIntegers_NamingTheMember()
    {
        const long big = (long)int.MaxValue + 1;
        EsaNode64[] nodes =
        [
            new(big, big, 1, 0, 0, true), new(0, big, 1, 0, 0, true), new(0, 1, big, 0, 0, true), new(0, 1, 1, big, 0, true), new(0, 1, 1, 0, big, true),
        ];
        foreach (EsaNode64 node in nodes)
        {
            NotSupportedException e = Assert.Throws<NotSupportedException>(() => node.ToNode("EsaNode64Test"))!;
            Assert.That(e.Message, Does.Contain("EsaNode64Test"));
        }
        Assert.That(new EsaNode64(int.MaxValue, int.MaxValue, int.MaxValue, int.MaxValue, int.MaxValue, true).ToNode("x").Pos, Is.EqualTo(int.MaxValue), "int.MaxValue itself fits");
    }

    [TestCase("banana")]
    [TestCase("mississippi")]
    [TestCase("aaaaaaaaaaaa")]
    [TestCase("abracadabraabracadabra")]
    public void TryMatchAt64Bits_OfAnIndexInMemory_IsTheIntervalOfTheThirtyTwoBitMember_ForEverySubstringAndItsNeighbours(string text)
    {
        EsaIndex index = EsaIndex.Build(text);
        var patterns = new HashSet<string>(StringComparer.Ordinal);
        for (int i = 0; i < text.Length; i++)
        {
            for (int len = 1; len <= text.Length - i; len++)
            {
                string p = text.Substring(i, len);
                patterns.Add(p);
                patterns.Add(p + "z");
                patterns.Add(p + text[0]);
            }
        }
        foreach (string pattern in patterns)
        {
            bool narrow = index.TryMatch(pattern.AsSpan(), out int lo, out int hi);
            bool wide = index.TryMatch(pattern.AsSpan(), out long wlo, out long whi);
            Assert.That(wide, Is.EqualTo(narrow), pattern);
            if (narrow) Assert.That((wlo, whi), Is.EqualTo(((long)lo, (long)hi)), pattern);
            Assert.That(index.LongCount(pattern), Is.EqualTo(narrow ? hi - lo + 1L : 0L), pattern);
        }
    }

    [Test]
    public void TheCountOfDistinctSubstrings_IsNarrowedTo64Bits_AndRefusesWhatDoesNotFit()
    {
        Assert.That(EsaIndex.ToCount(0), Is.EqualTo(0));
        Assert.That(EsaIndex.ToCount(12345), Is.EqualTo(12345));
        Assert.That(EsaIndex.ToCount(long.MaxValue), Is.EqualTo(long.MaxValue));
        NotSupportedException e = Assert.Throws<NotSupportedException>(() => EsaIndex.ToCount((Int128)long.MaxValue + 1))!;
        Assert.That(e.Message, Does.Contain("more than 64 bits"));
        Assert.That(() => EsaIndex.ToCount(Int128.MaxValue), Throws.TypeOf<NotSupportedException>());
    }
}
