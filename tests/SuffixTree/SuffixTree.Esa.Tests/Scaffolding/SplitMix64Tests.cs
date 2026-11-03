using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Scaffolding;

[TestFixture]
public class SplitMix64Tests
{
    [Test]
    public void Seed0_MatchesPublishedReferenceVector()
    {
        var rng = new SplitMix64(0);
        Assert.That(rng.Next(), Is.EqualTo(0xE220A8397B1DCDAFUL));
        Assert.That(rng.Next(), Is.EqualTo(0x6E789E6AA1B965F4UL));
        Assert.That(rng.Next(), Is.EqualTo(0x06C45D188009454FUL));
    }

    [Test]
    public void Seed1234567_MatchesPublishedReferenceVector()
    {
        var rng = new SplitMix64(1234567);
        Assert.That(rng.Next(), Is.EqualTo(6457827717110365317UL));
        Assert.That(rng.Next(), Is.EqualTo(3203168211198807973UL));
        Assert.That(rng.Next(), Is.EqualTo(9817491932198370423UL));
        Assert.That(rng.Next(), Is.EqualTo(4593380528125082431UL));
        Assert.That(rng.Next(), Is.EqualTo(16408922859458223821UL));
    }

    [Test]
    public void SameSeed_SameSequence_DifferentSeed_DifferentSequence()
    {
        var a = new SplitMix64(42);
        var b = new SplitMix64(42);
        var c = new SplitMix64(43);
        bool differs = false;
        for (int i = 0; i < 100; i++)
        {
            ulong x = a.Next();
            Assert.That(b.Next(), Is.EqualTo(x));
            differs |= c.Next() != x;
        }
        Assert.That(differs, Is.True);
    }

    [Test]
    public void NextInt_IsModuloOfNext_AndWithinBound()
    {
        var a = new SplitMix64(7);
        var b = new SplitMix64(7);
        foreach (int bound in new[] { 1, 2, 3, 4, 20, 255, 256, 1000, 65536, int.MaxValue })
        {
            int v = a.NextInt(bound);
            Assert.That(v, Is.InRange(0, bound - 1));
            Assert.That((ulong)v, Is.EqualTo(b.Next() % (ulong)bound));
        }
    }

    [Test]
    public void NextInt_NonPositiveBound_Throws()
    {
        var rng = new SplitMix64(1);
        Assert.Throws<ArgumentOutOfRangeException>(() => rng.NextInt(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => rng.NextInt(-5));
    }
}
