using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Contract;

/// <summary>The value types and enums of the index contract (work package 1b.0).</summary>
[TestFixture]
public class IndexContractTypesTests
{
    [Test]
    public void IndexLocus_ExposesItsPartsAndComparesComponentWise()
    {
        var a = new IndexLocus(3, 9, 4, 2);
        Assert.That((a.Lo, a.Hi, a.Length, a.Aux), Is.EqualTo((3L, 9L, 4L, 2L)));
        Assert.That(a, Is.EqualTo(new IndexLocus(3, 9, 4, 2)));
        Assert.That(a == new IndexLocus(3, 9, 4, 2), Is.True);
        Assert.That(a != new IndexLocus(3, 9, 4, 2), Is.False);
        Assert.That(a.Equals((object)new IndexLocus(3, 9, 4, 2)), Is.True);
        Assert.That(a.GetHashCode(), Is.EqualTo(new IndexLocus(3, 9, 4, 2).GetHashCode()));
        foreach (var other in new[] { new IndexLocus(4, 9, 4, 2), new IndexLocus(3, 8, 4, 2), new IndexLocus(3, 9, 5, 2), new IndexLocus(3, 9, 4, 3) })
        {
            Assert.That(a, Is.Not.EqualTo(other), other.ToString());
            Assert.That(a == other, Is.False, other.ToString());
            Assert.That(a != other, Is.True, other.ToString());
        }
        Assert.That(a.Equals(null), Is.False);
        Assert.That(a.Equals("x"), Is.False);
    }

    [Test]
    public void IndexLocus_DefaultIsTheAllZeroLocus()
    {
        Assert.That(default(IndexLocus), Is.EqualTo(new IndexLocus(0, 0, 0, 0)));
        Assert.That(new IndexLocus(1, 2, 3, 4).ToString(), Is.EqualTo("[1..2] depth 3 aux 4"));
    }

    [Test]
    public void IndexLocus_ValuesAre64Bit()
    {
        var big = new IndexLocus(long.MaxValue, long.MaxValue - 1, 1L << 40, -1);
        Assert.That(big.Lo, Is.EqualTo(long.MaxValue));
        Assert.That(big.Length, Is.EqualTo(1L << 40));
        Assert.That(big.Aux, Is.EqualTo(-1L));
    }

    [Test]
    public void Enums_HaveTheValuesOfADR001()
    {
        Assert.That(Enum.GetNames<IndexOperation>(), Is.EqualTo(new[]
        {
            "Count", "ExtendRight", "ExtendLeft", "LazyLocate", "TextAccess", "Parent", "Child", "SuffixLink", "StringDepth",
            "SuffixArrayAccess", "InverseSuffixArrayAccess", "LcpAccess", "Lcp", "Lca",
        }));
        Assert.That((int)IndexSupport.Unsupported, Is.EqualTo(0));
        Assert.That((int)IndexSupport.Degraded, Is.EqualTo(1));
        Assert.That((int)IndexSupport.Native, Is.EqualTo(2));
    }
}
