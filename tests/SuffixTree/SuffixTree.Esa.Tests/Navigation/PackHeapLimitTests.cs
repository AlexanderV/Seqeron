using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Navigation;

/// <summary>
/// Work package S2.5: the boundary of <see cref="EsaNavigationPack.PackHeapLimit"/> at both sides, with no allocation (the orchestrator found at S2.4 that a boundary-only mutant survives when no test passes the guard at exactly the limit):
/// a pack made over sources whose text has exactly the limit symbols passes the guard of the inverse suffix array and ends at the next line (it has no index to build it from: a different refusal that does not name the limit), a pack of one symbol
/// more is refused by name, and the policy of the 32-bit navigator (<see cref="EsaNavigationPack.InPlaceFrom"/>) is a constant between the two limits.
/// </summary>
[TestFixture]
public class PackHeapLimitTests
{
    /// <summary>The ranges of a text that is never read.</summary>
    private sealed class NoRanges : ILcpRanges
    {
        public long Length => 0;

        public int Fanout => 32;

        public long LevelBytes => 0;

        public long ValueAt(long rank) => throw new NotSupportedException();

        public long Min(long lo, long hi) => throw new NotSupportedException();

        public long PrevLess(long pos, long bound) => throw new NotSupportedException();

        public long NextLess(long pos, long bound) => throw new NotSupportedException();

        public void Borders(long a, long b, long t, out long left, out long right, out long leftValue, out long rightValue) => throw new NotSupportedException();

        public void Dispose()
        {
        }
    }

    private static EsaNavigationPack PackOf(long n) => EsaNavigationPack.FromSources(n, PosSource.Of(new Large.PowerReaders.Sa(n)), null, new NoRanges(), mapped: true);

    [TestCase(EsaNavigationPack.PackHeapLimit - 1, false)]
    [TestCase(EsaNavigationPack.PackHeapLimit, false)]
    [TestCase(EsaNavigationPack.PackHeapLimit + 1, true)]
    [TestCase(EsaNavigationPack.PackHeapLimit + 2, true)]
    public void TheGuardOfTheInverseSuffixArray_RefusesByName_AboveTheLimitAndNotAtIt(long n, bool refusedByName)
    {
        EsaNavigationPack pack = PackOf(n);
        NotSupportedException e = Assert.Throws<NotSupportedException>(() => pack.Isa(1, "TheMember"))!;
        if (refusedByName) Assert.That(e.Message, Does.Contain("TheMember").And.Contain("PackHeapLimit = 268435456").And.Contain("ISA"));
        else Assert.That(e.Message, Does.Contain("made over sources").And.Not.Contain("PackHeapLimit"), "at the limit the guard passes: the refusal is the pack's own (no index to build the array from), and no array was allocated");
    }

    [Test]
    public void ThePackOfAnIndexInMemory_IsNotLimited_AndTheLimitsAreOrdered()
    {
        Assert.That(EsaNavigationPack.PackHeapLimit, Is.EqualTo(1L << 28), "the 2^28 of ADR-005 section 3.3");
        Assert.That(EsaNavigationPack.InPlaceFrom, Is.LessThan(EsaNavigationPack.PackHeapLimit), "PSI is read in place from a size at which the heap array is still allowed");
        Assert.That(EsaNavigationPack.InPlaceFrom, Is.GreaterThan(10_000_000), "the texts of the S1 gate rows (10^7 symbols) keep the pack of S1");
        EsaIndex memory = EsaIndex.Build("mississippi");
        Assert.That(memory.Navigation.Isa(3, "RankOf"), Is.EqualTo(memory.Navigation.IsaArray[3]));
    }

    [Test]
    public void ALinkOverAPackWithoutPsiOrIsa_IsRefusedByName_NotOutOfMemory()
    {
        EsaNavigationPack pack = PackOf(EsaNavigationPack.PackHeapLimit + 1);
        var node = new EsaNode64(5, 9, 4, 1, 3, true);
        NotSupportedException e = Assert.Throws<NotSupportedException>(() => pack.SuffixLink(node))!;
        Assert.That(e.Message, Does.Contain("the suffix link").And.Contain("PackHeapLimit"));
    }
}
