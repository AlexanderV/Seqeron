using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.4, the small-n form of the generator-backed matrix and of the container reader's new members (the default set; the same checks at 2^31 + 10 .. 2^33 are <see cref="WideEngineLargeTests"/> and
/// <see cref="WideEngineSparseLargeTests"/>): the cheap readers of <c>a^n</c> are the closed form, the descent over them counts every run, the statistics of a container are the numbers of the in-memory index at 64 bits, and the
/// sections that the long descent reads (<see cref="SesaFile.WideSections"/>) answer what the 32-bit descent answers on a container of S1.
/// </summary>
[TestFixture]
public class WideEngineSelfTests
{
    [Test]
    public void ThePowerReaders_AreTheClosedFormOfPowerOfA([Values(1L, 2L, 5L, 64L, 1000L)] long n)
    {
        var f = new PowerOfA(n);
        var sa = new PowerReaders.Sa(n);
        var child = new PowerReaders.Child(n);
        var lcp = new PowerReaders.Lcp(n);
        var text = new PowerReaders.Text(n);
        for (long r = 0; r <= n; r++)
        {
            Assert.That((long)sa.Get(r), Is.EqualTo(f.SuffixAt(r)), $"SA[{r}]");
            Assert.That((long)child.Get(r), Is.EqualTo(f.PackedCell(r)), $"CHILD[{r}]");
            Assert.That((long)lcp.Get(r), Is.EqualTo(f.Lcp(r)), $"LCP[{r}]");
        }
        foreach (long outside in new[] { n + 1, -1 })
        {
            Assert.That(() => sa.Get(outside), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => child.Get(outside), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => lcp.Get(outside), Throws.TypeOf<InvalidDataException>());
        }
        Assert.That(text.Symbol(n - 1), Is.EqualTo(0));
        Assert.That(() => text.Symbol(n), Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    [Test]
    public void TheCountOfARun_OverThePowerReaders_IsTheClosedForm_AtEveryK([Values(1L, 2L, 3L, 100L, 1200L)] long n)
    {
        var text = new PowerReaders.Text(n);
        var cells = new WideCells<CheckedAccess, PowerReaders.Child, PowerReaders.Lcp>(new PowerReaders.Child(n), new PowerReaders.Lcp(n), GenDescent.Lcp64(new PowerOfA(n)), tagged: false, n);
        for (long k = 1; k <= n + 2; k++)
        {
            bool found = WideSearch.Run<CheckedAccess, PowerReaders.Sa, WideCells<CheckedAccess, PowerReaders.Child, PowerReaders.Lcp>, PowerReaders.Text, RunPattern>(new PowerReaders.Sa(n), cells, text, n, new RunPattern(k), out long lo, out long hi);
            Assert.That(found, Is.EqualTo(k <= n), $"a^{k} in a^{n}");
            if (found) Assert.That((lo, hi, hi - lo + 1), Is.EqualTo((k, n, n - k + 1)), $"a^{k} in a^{n}");
            else Assert.That((lo, hi), Is.EqualTo((0L, -1L)));
        }
    }

    // ------------------------------------------------------------------------------------------ the container reader

    [TestCase("")]
    [TestCase("a")]
    [TestCase("mississippi")]
    [TestCase("abracadabra abracadabra")]
    public void TheStatisticsOfAContainer_AreThoseOfTheInMemoryIndex_At64Bits(string text)
    {
        EsaIndex memory = EsaIndex.Build(text);
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(text));
        string path = Path.Combine(Path.GetTempPath(), "esa-wide-" + Guid.NewGuid().ToString("N") + ".sesa");
        File.WriteAllBytes(path, image);
        try
        {
            using SesaFile file = SesaFile.Open(path, new EsaOpenOptions());
            EsaStatistics st = file.Statistics;
            EsaStatistics expected = memory.Statistics;
            Assert.That((st.Length, st.NodeCount, st.InternalCount, st.LeafCount, st.MaxDepth, st.LrsDepth, st.LrsRank), Is.EqualTo((expected.Length, expected.NodeCount, expected.InternalCount, expected.LeafCount, expected.MaxDepth, expected.LrsDepth, expected.LrsRank == -1 ? 0 : expected.LrsRank)));
            Assert.That(file.IsNarrow, Is.True, "a text of 2^31 - 1 symbols and more is the only wide container at POS_WIDTH 4");
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Test]
    public void TheWideSectionsOfAContainerOfS1_AnswerWhatTheThirtyTwoBitDescentAnswers([Values(false, true)] bool tagged, [Values(EsaAccessMode.Checked, EsaAccessMode.Trusted)] EsaAccessMode mode)
    {
        var rng = new SplitMix64(0x51DE);
        var chars = new char[600];
        for (int i = 0; i < chars.Length; i++) chars[i] = (char)('A' + rng.NextInt(4));
        string text = new string(chars) + new string('A', 300) + new string(chars, 0, 200);
        EsaIndex memory = EsaIndex.Build(text);
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(text), codec: tagged ? ChildCodec.Tagged : ChildCodec.Packed);
        string path = Path.Combine(Path.GetTempPath(), "esa-wide-" + Guid.NewGuid().ToString("N") + ".sesa");
        File.WriteAllBytes(path, image);
        try
        {
            using SesaFile file = SesaFile.Open(path, new EsaOpenOptions { Access = mode });
            var narrow = new MappedSearch(file);
            Assert.That(narrow.IsWide, Is.False);
            WideSearcher wide = WideSearcher.Create(file.WideSections());
            char[] alphabet = file.TextView.Alphabet;
            for (int i = 0; i < 400; i++)
            {
                int len = 1 + rng.NextInt(40);
                string pattern = i % 3 == 0 ? new string([.. Enumerable.Range(0, len).Select(_ => (char)('A' + rng.NextInt(4)))]) : text.Substring(rng.NextInt(text.Length - len), len);
                byte[] dense = [.. pattern.Select(c => (byte)Array.IndexOf(alphabet, c))];
                bool expected = memory.TryMatch(pattern.AsSpan(), out int elo, out int ehi);
                Assert.That(narrow.TryMatch(dense, out int nlo, out int nhi), Is.EqualTo(expected));
                Assert.That(wide.TryMatch(dense, out long wlo, out long whi), Is.EqualTo(expected));
                Assert.That((wlo, whi), Is.EqualTo(((long)elo, (long)ehi)), pattern);
                Assert.That((nlo, nhi), Is.EqualTo((elo, ehi)), pattern);
                Assert.That(narrow.TryMatch(dense, out long l64, out long h64), Is.EqualTo(expected));
                Assert.That((l64, h64), Is.EqualTo(((long)elo, (long)ehi)));
            }
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Test]
    public void TheIndexOfAContainerOfS1_ServesBothWidthsOfTheInterface_ThroughTheThirtyTwoBitLoop()
    {
        const string text = "mississippimississippi";
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(text));
        string path = Path.Combine(Path.GetTempPath(), "esa-wide-" + Guid.NewGuid().ToString("N") + ".sesa");
        File.WriteAllBytes(path, image);
        try
        {
            using SesaFile file = SesaFile.Open(path, new EsaOpenOptions());
            EsaIndex index = EsaIndex.FromFile(file);
            EsaIndex memory = EsaIndex.Build(text);
            Assert.That(index.TextLength, Is.EqualTo(text.Length));
            Assert.That(index.Length, Is.EqualTo(text.Length));
            foreach (string p in new[] { "i", "ss", "issi", "ppi", "mississippi", "z", "ssippim" })
            {
                Assert.That(index.LongCount(p), Is.EqualTo((long)memory.CountOccurrences(p)), p);
                Assert.That(index.CountOccurrences(p), Is.EqualTo(memory.CountOccurrences(p)), p);
                Assert.That(index.Contains(p), Is.EqualTo(memory.Contains(p)), p);
                Assert.That(index.FindAllOccurrences(p), Is.EqualTo(memory.FindAllOccurrences(p)), p);
            }
        }
        finally
        {
            File.Delete(path);
        }
    }
}
