using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Scaffolding;

[TestFixture]
public class StringFamiliesTests
{
    [Test]
    public void AllSame_HasLengthAndSingleSymbol()
    {
        Assert.That(StringFamilies.AllSame(0), Is.EqualTo(""));
        Assert.That(StringFamilies.AllSame(5), Is.EqualTo("aaaaa"));
        Assert.That(StringFamilies.AllSame(3, 'z'), Is.EqualTo("zzz"));
    }

    [Test]
    public void Fibonacci_StartsWithTheStandardWord_AndTruncates()
    {
        Assert.That(StringFamilies.Fibonacci(8), Is.EqualTo("abaababa"));
        Assert.That(StringFamilies.Fibonacci(13), Is.EqualTo("abaababaabaab"));
        Assert.That(StringFamilies.Fibonacci(0), Is.EqualTo(""));
        Assert.That(StringFamilies.Fibonacci(1), Is.EqualTo("a"));
        Assert.That(StringFamilies.Fibonacci(1000), Has.Length.EqualTo(1000));
    }

    [Test]
    public void ThueMorse_IsPopcountParity()
    {
        Assert.That(StringFamilies.ThueMorse(16), Is.EqualTo("abbabaabbaababba"));
        string t = StringFamilies.ThueMorse(1000);
        for (int i = 0; i < t.Length; i++)
            Assert.That(t[i], Is.EqualTo((char)('a' + (int.PopCount(i) & 1))));
    }

    [TestCase(2, 1)]
    [TestCase(2, 2)]
    [TestCase(2, 3)]
    [TestCase(2, 8)]
    [TestCase(2, 10)]
    [TestCase(4, 1)]
    [TestCase(4, 2)]
    [TestCase(4, 4)]
    [TestCase(4, 5)]
    public void DeBruijn_EveryCyclicWindowAppearsExactlyOnce(int k, int order)
    {
        string s = StringFamilies.DeBruijn(k, order);
        int expected = (int)Math.Pow(k, order);
        Assert.That(s, Has.Length.EqualTo(expected));
        var seen = new HashSet<string>();
        string cyc = string.Concat(s, s.AsSpan(0, order - 1));
        for (int i = 0; i < s.Length; i++)
            Assert.That(seen.Add(cyc.Substring(i, order)), Is.True, $"window at {i} repeats");
        Assert.That(s.All(c => c >= 'a' && c < 'a' + k), Is.True);
    }

    [Test]
    public void DeBruijn_Truncation_IsPrefixOfTheFullSequence()
    {
        string full = StringFamilies.DeBruijn(2, 6);
        Assert.That(StringFamilies.DeBruijn(2, 6, 40), Is.EqualTo(full.Substring(0, 40)));
        Assert.Throws<ArgumentOutOfRangeException>(() => StringFamilies.DeBruijn(2, 6, 65));
    }

    [Test]
    public void AkB_IsPowerOfAThenB()
    {
        Assert.That(StringFamilies.AkB(11, 3), Is.EqualTo("aaabaaabaaa"));
        Assert.That(StringFamilies.AkB(4, 0), Is.EqualTo("bbbb"));
    }

    [Test]
    public void AbRepeat_Alternates()
    {
        Assert.That(StringFamilies.AbRepeat(7), Is.EqualTo("abababa"));
    }

    [TestCase(2)]
    [TestCase(16)]
    [TestCase(256)]
    [TestCase(65536)]
    public void Runs_AlternateLettersAndRespectTheMaximum(int maxRun)
    {
        string s = StringFamilies.Runs(5000, maxRun, 99);
        Assert.That(s, Has.Length.EqualTo(5000));
        int i = 0;
        char expected = 'a';
        while (i < s.Length)
        {
            int j = i;
            while (j < s.Length && s[j] == s[i]) j++;
            Assert.That(s[i], Is.EqualTo(expected));
            if (j < s.Length) Assert.That(j - i, Is.InRange(1, maxRun));
            expected = expected == 'a' ? 'b' : 'a';
            i = j;
        }
        Assert.That(StringFamilies.Runs(5000, maxRun, 99), Is.EqualTo(s));
    }

    [TestCase(2)]
    [TestCase(3)]
    [TestCase(5)]
    [TestCase(7)]
    [TestCase(64)]
    [TestCase(65)]
    public void Periodic_HasThePeriod(int p)
    {
        foreach (string s in new[] { StringFamilies.Periodic(1000, p, 5), StringFamilies.PeriodicAaB(1000, p) })
        {
            Assert.That(s, Has.Length.EqualTo(1000));
            for (int i = 0; i + p < s.Length; i++)
                Assert.That(s[i + p], Is.EqualTo(s[i]));
        }
        Assert.That(StringFamilies.PeriodicAaB(2 * p, p), Is.EqualTo(new string('a', p - 1) + "b" + new string('a', p - 1) + "b"));
    }

    [TestCase(1)]
    [TestCase(2)]
    [TestCase(3)]
    [TestCase(4)]
    [TestCase(20)]
    [TestCase(256)]
    [TestCase(1000)]
    public void Random_UsesOnlyTheRequestedAlphabet_AndIsDeterministic(int sigma)
    {
        string s = StringFamilies.Random(3000, sigma, 11);
        var alphabet = StringFamilies.Alphabet(sigma);
        Assert.That(alphabet, Has.Length.EqualTo(sigma));
        Assert.That(alphabet.Distinct().Count(), Is.EqualTo(sigma));
        Assert.That(s, Has.Length.EqualTo(3000));
        Assert.That(s.All(c => Array.IndexOf(alphabet, c) >= 0), Is.True);
        Assert.That(StringFamilies.Random(3000, sigma, 11), Is.EqualTo(s));
        if (sigma > 1)
            Assert.That(StringFamilies.Random(3000, sigma, 12), Is.Not.EqualTo(s));
        if (sigma <= 20)
            Assert.That(s.Distinct().Count(), Is.EqualTo(sigma));
    }

    [Test]
    public void Alphabet_IsStrictlyIncreasing()
    {
        foreach (int sigma in new[] { 1, 2, 4, 20, 26, 27, 256, 1000, 65536 })
        {
            var a = StringFamilies.Alphabet(sigma);
            for (int i = 1; i < a.Length; i++) Assert.That(a[i], Is.GreaterThan(a[i - 1]));
        }
    }

    [Test]
    public void RandomDna_RandomProtein_UseTheirLetters()
    {
        Assert.That(StringFamilies.RandomDna(2000, 3).All(c => "ACGT".Contains(c)), Is.True);
        string p = StringFamilies.RandomProtein(5000, 3);
        Assert.That(p.All(c => "ACDEFGHIKLMNPQRSTVWY".Contains(c)), Is.True);
        Assert.That(p.Distinct().Count(), Is.EqualTo(20));
    }

    [Test]
    public void Skewed_HasOneRareSymbol()
    {
        string s = StringFamilies.Skewed(20000, 8);
        int rare = s.Count(c => c == 'z');
        Assert.That(rare, Is.InRange(1, 500));
        Assert.That(s.Where(c => c != 'z').All(c => c is 'a' or 'b' or 'c'), Is.True);
    }

    [Test]
    public void AllDistinct_IsAPermutationOfTheFirstNCodeUnits()
    {
        foreach (int n in new[] { 0, 1, 2, 300, 65536 })
        {
            string s = StringFamilies.AllDistinct(n, 4);
            Assert.That(s, Has.Length.EqualTo(n));
            Assert.That(s.Distinct().Count(), Is.EqualTo(n));
            Assert.That(s.All(c => c < n), Is.True);
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => StringFamilies.AllDistinct(65537, 1));
    }

    [Test]
    public void Utf16Families_ContainTheClaimedUnits()
    {
        string bmp = StringFamilies.Bmp(500, 1);
        Assert.That(bmp.All(c => c > 255 && !char.IsSurrogate(c)), Is.True);

        string pairs = StringFamilies.SurrogatePairs(500, 1);
        Assert.That(pairs, Has.Length.EqualTo(500));
        for (int i = 0; i < pairs.Length; i += 2)
        {
            if (i + 1 < pairs.Length) Assert.That(char.IsSurrogatePair(pairs[i], pairs[i + 1]), Is.True);
        }

        string lone = StringFamilies.LoneSurrogates(500, 1);
        int loneCount = 0;
        for (int i = 0; i < lone.Length; i++)
        {
            bool validPair = (char.IsHighSurrogate(lone[i]) && i + 1 < lone.Length && char.IsLowSurrogate(lone[i + 1]))
                          || (char.IsLowSurrogate(lone[i]) && i > 0 && char.IsHighSurrogate(lone[i - 1]));
            if (char.IsSurrogate(lone[i]) && !validPair) loneCount++;
        }
        Assert.That(loneCount, Is.GreaterThan(0));

        string nul = StringFamilies.NulAndFfff(500, 1);
        Assert.That(nul, Does.Contain('\0'));
        Assert.That(nul, Does.Contain('￿'));
    }

    [Test]
    public void Adversarial_ListsUniquelyNamedFamilies_OfExactLength()
    {
        foreach (int n in new[] { 0, 1, 2, 3, 17, 300 })
        {
            var all = StringFamilies.Adversarial(n, 1).ToList();
            Assert.That(all.Select(x => x.Name).Distinct().Count(), Is.EqualTo(all.Count));
            Assert.That(all.Count, Is.GreaterThanOrEqualTo(25));
            foreach (var (name, text) in all)
                Assert.That(text, Has.Length.EqualTo(n), name);
        }
    }

    [Test]
    public void Adversarial_IsDeterministic()
    {
        var a = StringFamilies.Adversarial(200, 7).ToList();
        var b = StringFamilies.Adversarial(200, 7).ToList();
        Assert.That(b, Is.EqualTo(a));
    }
}
