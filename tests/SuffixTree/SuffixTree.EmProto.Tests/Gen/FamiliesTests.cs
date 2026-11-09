using SuffixTree.EmProto.Gen;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.Gen;

[TestFixture]
public class FamiliesTests
{
    private static string Decode(SymbolText t) => new([.. t.Symbols.Select(s => t.Letters[s])]);

    [TestCase(0)]
    [TestCase(1)]
    [TestCase(2)]
    [TestCase(3)]
    [TestCase(100)]
    [TestCase(1000)]
    [TestCase(12345)]
    public void AdversarialFamilies_EqualTheEsaTestGenerators(int n)
    {
        Assert.That(Decode(Families.AllSame(n)), Is.EqualTo(StringFamilies.AllSame(n)));
        Assert.That(Decode(Families.Fibonacci(n)), Is.EqualTo(StringFamilies.Fibonacci(n)));
        Assert.That(Decode(Families.ThueMorse(n)), Is.EqualTo(StringFamilies.ThueMorse(n)));
        Assert.That(Decode(Families.RandomDna(n, 5)), Is.EqualTo(StringFamilies.RandomDna(n, 5)));
    }

    [Test]
    public void Tandem_IsTheUnitRepeated()
    {
        SymbolText t = Families.Tandem(171, 20, 3);
        Assert.That(t.Symbols, Has.Length.EqualTo(171 * 20));
        for (int i = 171; i < t.Symbols.Length; i++) Assert.That(t.Symbols[i], Is.EqualTo(t.Symbols[i - 171]));
        Assert.That(t.Sigma, Is.EqualTo(4));
    }

    [Test]
    public void NRun_HasOneRunOfNInsideAFiveLetterAlphabet()
    {
        SymbolText t = Families.NRun(10_000, 3_000, 1_000, 9);
        Assert.That(t.Letters, Is.EqualTo("ACGNT"));
        for (int i = 3_000; i < 4_000; i++) Assert.That(t.Symbols[i], Is.EqualTo(3));
        Assert.That(t.Symbols.Count(s => s == 3), Is.EqualTo(1_000), "N occurs only inside the run");
        Assert.That(t.Symbols.Max(), Is.EqualTo(4));
    }

    [Test]
    public void GenomeLike_IsDeterministic_HasAllSymbolsAndRepeatsAndSatellites()
    {
        SymbolText a = Families.GenomeLike(2_000_000, 11);
        SymbolText b = Families.GenomeLike(2_000_000, 11);
        SymbolText c = Families.GenomeLike(2_000_000, 12);
        Assert.That(a.Symbols, Is.EqualTo(b.Symbols));
        Assert.That(a.Symbols, Is.Not.EqualTo(c.Symbols));
        Assert.That(a.Letters, Is.EqualTo("ACGNT"));
        Assert.That(a.Symbols, Has.Length.EqualTo(2_000_000));
        long nCount = a.Symbols.LongCount(s => s == 3);
        Assert.That(nCount / 2_000_000.0, Is.InRange(0.005, 0.02), "about 1 % N");
        // interspersed repeats and satellites make 16-mers recur far more often than in random DNA
        var seen = new Dictionary<long, int>();
        long repeated = 0;
        for (int i = 0; i + 16 <= a.Symbols.Length; i++)
        {
            long key = 0;
            for (int k = 0; k < 16; k++) key = (key * 5) + a.Symbols[i + k];
            repeated += seen.TryGetValue(key, out int c0) && c0 > 0 ? 1 : 0;
            seen[key] = seen.GetValueOrDefault(key) + 1;
        }
        Assert.That(repeated, Is.GreaterThan(50_000), "random DNA would give about 500 recurring 16-mers here");
    }

    [Test]
    public void ToAscii_AndFromFile_RoundTripThroughAFile()
    {
        SymbolText t = Families.NRun(5_000, 1_000, 300, 4);
        string path = Path.Combine(Path.GetTempPath(), "emproto-fam-" + Guid.NewGuid().ToString("N") + ".txt");
        try
        {
            File.WriteAllBytes(path, Families.ToAscii(t));
            SymbolText back = Families.FromFile(path, "x");
            Assert.That(back.Symbols, Is.EqualTo(t.Symbols));
            Assert.That(back.Letters, Is.EqualTo(t.Letters));
        }
        finally
        {
            File.Delete(path);
        }
    }
}
