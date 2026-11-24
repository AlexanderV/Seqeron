using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.Fuzz;
using SuffixTree.Esa.Tests.Large;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Esa.Tests.Views;

namespace SuffixTree.Esa.Tests.Search;

/// <summary>
/// The long descent (<see cref="WideSearch"/>, <see cref="WideSearcher"/>; work package S2.4) at both stored widths, in both access modes, over every kind of LCP array and both CHILD codecs: on a valid image it answers what the memory
/// index answers, interval for interval (and with the interval of an absent pattern the same as the 32-bit loop gives); on a damaged one it ends as the 32-bit descent and the readable references do; over the generated views of a
/// closed form (<see cref="PosGen{TWidth}"/>) it finds the factors of the text and the runs <c>a^k</c> at every length. Written before the loop existed (the stub said <see cref="NotImplementedException"/>).
/// </summary>
[TestFixture]
public unsafe class WideDescentTests
{
    private static readonly EsaAccessMode[] Modes = [EsaAccessMode.Checked, EsaAccessMode.Trusted];

    /// <summary>The layouts of the matrix: the stored width, the kind of the LCP array, and the CHILD codec (the tagged one exists over 4 bytes only). 12 rows.</summary>
    internal static readonly (int Width, LcpKind Kind, bool Tagged)[] Layouts = WideImage.Layouts;

    private static byte[] Dense(char[] alphabet, string pattern)
    {
        var dense = new byte[pattern.Length];
        for (int i = 0; i < dense.Length; i++)
        {
            int s = Array.IndexOf(alphabet, pattern[i]);
            if (s < 0) return [];
            dense[i] = (byte)s;
        }
        return dense;
    }

    private static string RandomText(SplitMix64 rng, int n, int sigma)
    {
        var chars = new char[n];
        for (int i = 0; i < n; i++) chars[i] = (char)('A' + rng.NextInt(sigma));
        return new string(chars);
    }

    private static string RandomPattern(SplitMix64 rng, string text, int sigma, int length, bool fromText)
    {
        if (fromText && text.Length >= length) return text.Substring(rng.NextInt(text.Length - length + 1), length);
        var chars = new char[length];
        for (int i = 0; i < length; i++) chars[i] = (char)('A' + rng.NextInt(sigma));
        return new string(chars);
    }

    private static IEnumerable<string> AllSubstrings(string text, int maxLength)
    {
        for (int i = 0; i < text.Length; i++)
        {
            for (int len = 1; len <= Math.Min(maxLength, text.Length - i); len++) yield return text.Substring(i, len);
        }
    }

    /// <summary>Every layout and mode against the memory index for the given patterns.</summary>
    private static void CompareAll(string text, string label, IEnumerable<string> patterns)
    {
        EsaIndex memory = EsaIndex.Build(text);
        LcpChildCase c = LcpChildCase.Of(text);
        char[] alphabet = [.. text.Distinct().Order()];
        string[] list = [.. patterns];
        foreach ((int width, LcpKind kind, bool tagged) in Layouts)
        {
            using var image = new WideImage(c, width, kind, tagged);
            foreach (EsaAccessMode mode in Modes)
            {
                WideSearcher searcher = WideSearcher.Create(image.Sections(mode == EsaAccessMode.Checked));
                foreach (string pattern in list)
                {
                    bool expected = memory.TryMatch(pattern.AsSpan(), out int elo, out int ehi);
                    byte[] dense = Dense(alphabet, pattern);
                    if (dense.Length != pattern.Length)
                    {
                        Assert.That(expected, Is.False, $"{label}: a pattern with a letter outside the alphabet cannot occur");
                        continue;
                    }
                    bool found = searcher.TryMatch(dense, out long lo, out long hi);
                    string where = $"{label}: width {width} {kind}{(tagged ? " tagged" : string.Empty)} {mode}, pattern '{Show(pattern)}'";
                    Assert.That(found, Is.EqualTo(expected), where);
                    Assert.That((lo, hi), Is.EqualTo(((long)elo, (long)ehi)), where);
                }
            }
        }
    }

    private static string Show(string s) => s.Length <= 60 ? s : s[..60] + "...";

    [TestCase("banana")]
    [TestCase("mississippi")]
    [TestCase("abracadabra")]
    [TestCase("aaaaaaaaaaaa")]
    [TestCase("abababababababab")]
    [TestCase("a")]
    [TestCase("ab")]
    [TestCase("aa")]
    [TestCase("abcabcabcabcabcabcx")]
    public void EverySubstring_AndItsNeighbours_OfANamedText_EqualsTheMemoryIndex_AtEveryLayout(string text)
    {
        var patterns = new List<string>(AllSubstrings(text, 16));
        patterns.AddRange(patterns.Take(200).Select(p => p + text[0]));      // a pattern that usually does not occur
        patterns.AddRange(patterns.Take(200).Select(p => text[^1] + p));
        patterns.Add(text + text);                                          // longer than the text
        patterns.Add(text + "Z");                                           // a letter outside the alphabet
        CompareAll(text, text, patterns);
    }

    [TestCase(2)]
    [TestCase(4)]
    [TestCase(5)]
    [TestCase(20)]
    [TestCase(100)]
    public void RandomTexts_OfSeveralAlphabets_EqualTheMemoryIndex_AtEveryLayout(int sigma)
    {
        var rng = new SplitMix64(0x5EED1000UL + (ulong)sigma);
        foreach (int n in new[] { 7, 64, 900 })
        {
            string text = RandomText(rng, n, sigma);
            var patterns = new List<string>();
            for (int i = 0; i < 200; i++)
            {
                int len = 1 + rng.NextInt(Math.Min(70, n + 3));
                patterns.Add(RandomPattern(rng, text, sigma, len, fromText: i % 3 != 0));
            }
            CompareAll(text, $"sigma {sigma} n {n}", patterns);
        }
    }

    [Test]
    public void EveryKindOfTheLcpArray_WithItsEscapes_EqualsTheMemoryIndex()
    {
        var rng = new SplitMix64(0x1CC4);
        string repeat = RandomText(rng, 700, 4);
        var texts = new Dictionary<string, string>
        {
            ["a block twice (an 8-bit LCP array with escapes)"] = RandomText(rng, 600, 4) + repeat + RandomText(rng, 600, 4) + repeat,
            ["a^2000 (escapes at every rank above 254)"] = new string('A', 2000),
            ["a^70000 (a 16-bit escape at 65,535)"] = new string('A', 70_000),
        };
        foreach ((string name, string text) in texts)
        {
            var patterns = new List<string>();
            for (int i = 0; i < 100; i++) patterns.Add(RandomPattern(rng, text, 4, 1 + rng.NextInt(900), fromText: i % 4 != 0));
            patterns.AddRange(AllSubstrings(text.Length > 400 ? text[..30] : text, 30));
            if (text.Length == 70_000) patterns.AddRange(new[] { 254, 255, 256, 65_534, 65_535, 65_536, 69_999, 70_000, 70_001 }.Select(len => new string('A', len)));
            else if (text.Length == 2000) patterns.AddRange(new[] { 253, 254, 255, 256, 1999, 2000, 2001 }.Select(len => new string('A', len)));
            CompareAll(text, name, patterns);
        }
    }

    [Test]
    public void ThePatternOfASpan_IsBoundedByItsSymbols()
    {
        var pattern = new SpanPattern([3, 1, 2]);
        Assert.That(pattern.Length, Is.EqualTo(3));
        Assert.That((pattern[0], pattern[1], pattern[2]), Is.EqualTo((3, 1, 2)));
        foreach (long outside in new long[] { 3, 4, -1, long.MaxValue })
        {
            try
            {
                _ = pattern[outside];
                Assert.Fail($"the symbol {outside} of a pattern of 3 symbols was read");
            }
            catch (ArgumentOutOfRangeException)
            {
                // the refusal
            }
        }
    }

    // ---------------------------------------------------------------- the instantiations

    [Test]
    public void TheSwitchOfTheMatrix_MakesOneSearcherPerAccessWidthAndKind_16Instantiations()
    {
        LcpChildCase c = LcpChildCase.Of("mississippi");
        var seen = new HashSet<string>();
        foreach ((int width, LcpKind kind, bool tagged) in Layouts)
        {
            using var image = new WideImage(c, width, kind, tagged);
            foreach (bool isChecked in new[] { true, false })
            {
                Type searcher = WideSearcher.Create(image.Sections(isChecked)).GetType();
                Type[] args = searcher.GenericTypeArguments;
                Assert.That(args, Has.Length.EqualTo(3));
                Assert.That(args[0], Is.EqualTo(isChecked ? typeof(CheckedAccess) : typeof(TrustedAccess)));
                Assert.That(args[1].GetGenericTypeDefinition(), Is.EqualTo(width == 5 ? typeof(Pos40<>) : typeof(Pos32<>)));
                Type lcp = args[2].GetGenericTypeDefinition();
                Assert.That(lcp, Is.EqualTo(kind switch { LcpKind.U8 => typeof(Pos8<>), LcpKind.U16 => typeof(Pos16<>), LcpKind.U32 => typeof(Pos32<>), _ => typeof(Pos64<>) }));
                seen.Add(searcher.ToString());
            }
        }
        Assert.That(seen, Has.Count.EqualTo(2 * 2 * 4), "access x width x LCP kind; the tagged codec is a flag, not a type");
    }

    [Test]
    public void TheSwitchOfTheMatrix_RefusesAWidthThatIsNot4Or5_AndATaggedTableOver5ByteCells()
    {
        LcpChildCase c = LcpChildCase.Of("banana");
        using var image = new WideImage(c, 4, LcpKind.U32, tagged: false);
        WideSections ok = image.Sections(true);
        WideSections bad6 = new(true, 6, ok.LcpKind, false, ok.Sa, ok.Child, ok.LcpBase, ok.Lcp, ok.Text, ok.Bits, ok.N);
        Assert.That(() => WideSearcher.Create(bad6), Throws.TypeOf<ArgumentOutOfRangeException>());
        WideSections taggedWide = new(true, 5, ok.LcpKind, true, ok.Sa, ok.Child, ok.LcpBase, ok.Lcp, ok.Text, ok.Bits, ok.N);
        Assert.That(() => WideSearcher.Create(taggedWide), Throws.TypeOf<ArgumentException>());
    }

    // ---------------------------------------------------------------- damaged images

    private static string Outcome(Func<(bool Found, long Lo, long Hi)> run)
    {
        try
        {
            (bool found, long lo, long hi) = run();
            return found ? $"found [{lo}..{hi}]" : "absent";
        }
        catch (InvalidDataException)
        {
            return "InvalidDataException";
        }
    }

    private static void Damage(byte[] image, SplitMix64 rng, int kind)
    {
        SesaOracle oracle = SesaOracle.FromBytes(image);
        OracleSection section = oracle.Get(kind);
        int bytes = 1 + rng.NextInt(4);
        for (int i = 0; i < bytes; i++)
        {
            long at = section.Offset + rng.NextInt((int)section.Length);
            image[at] = rng.NextInt(3) switch { 0 => 0, 1 => 0xFF, _ => (byte)rng.NextInt(256) };
        }
    }

    private static WideSections SectionsOf(SesaFile file) =>
        new(file.IsChecked, file.Header.PosWidth, file.Parsed.LcpKind, file.Parsed.Child.Codec == SesaFormat.ChildCodecTagged, file.Data + file.Parsed.Sa.Offset, file.Data + file.Parsed.Child.Offset, file.Data + file.Parsed.Lcp.Offset, file.Lcp, file.Data + file.Parsed.Text.Offset, file.Parsed.Bits, file.Header.TextLength);

    [Test]
    public void ADamagedContainer_InCheckedMode_EndsAsTheThirtyTwoBitDescentAndTheReferencesDo_AndNeverReadsOutsideTheMapping()
    {
        Assume.That(OperatingSystem.IsLinux(), "the guard pages need mmap and mprotect");
        var rng = new SplitMix64(0xDA6A6E77);
        string text = RandomText(rng, 700, 4) + new string('A', 300) + RandomText(rng, 400, 4);
        LcpChildCase c = LcpChildCase.Of(text);
        byte[][] cleans = [SesaTestFile.Image(c, codec: ChildCodec.Packed), SesaTestFile.Image(c, codec: ChildCodec.Tagged)];
        int[] kinds = [SesaFormat.KindSa, SesaFormat.KindChild, SesaFormat.KindLcp];
        var options = new EsaOpenOptions { Access = EsaAccessMode.Checked };
        int refused = 0, answered = 0;
        for (int round = 0; round < 600; round++)
        {
            byte[] image = (byte[])cleans[round % 2].Clone();
            int damaged = 1 + rng.NextInt(2);
            for (int d = 0; d < damaged; d++) Damage(image, rng, kinds[rng.NextInt(kinds.Length)]);
            using var guarded = new SesaReaderHostileTests.Guarded(image, flush: round % 2 == 0);
            SesaFile file;
            try
            {
                file = SesaFile.FromMemory(guarded.Data, guarded.Length, options);
            }
            catch (InvalidDataException)
            {
                continue; // the header level refused it
            }
            using (file)
            {
                var narrow = new MappedSearch(file);
                WideSearcher wide = WideSearcher.Create(SectionsOf(file));
                char[] alphabet = file.TextView.Alphabet;
                for (int i = 0; i < 40; i++)
                {
                    string pattern = RandomPattern(rng, text, 4, 1 + rng.NextInt(60), fromText: i % 3 != 0);
                    byte[] dense = Dense(alphabet, pattern);
                    string expected = Outcome(() => (EsaIndex.DescendMapped(file.Child.Mapped, file.Sa, file.Lcp, file.TextReader, dense, true, out int lo, out int hi), lo, hi));
                    string thirtyTwo = Outcome(() => (narrow.TryMatch(dense, out int lo, out int hi), lo, hi));
                    string actual = Outcome(() => (wide.TryMatch(dense, out long lo, out long hi), lo, hi));
                    Assert.That(thirtyTwo, Is.EqualTo(expected), $"round {round}, pattern '{pattern}': the 32-bit loop against its reference");
                    Assert.That(actual, Is.EqualTo(expected), $"round {round}, pattern '{pattern}'");
                    if (actual == "InvalidDataException") refused++;
                    else answered++;
                }
            }
        }
        Assert.That(refused, Is.GreaterThan(0), "the damage must reach the descent");
        Assert.That(answered, Is.GreaterThan(0));
    }

    [TestCase(false)]
    [TestCase(true)]
    public void ATableThatSendsTheWalkBackwardsOrOutsideItsInterval_InCheckedMode_IsInvalidData_AsInThe32BitDescent(bool tagged)
    {
        ChildCodec codec = tagged ? ChildCodec.Tagged : ChildCodec.Packed;
        string text = "mississippimississippi";
        byte[] clean = SesaTestFile.Image(LcpChildCase.Of(text), codec: codec);
        OracleSection child = SesaOracle.FromBytes(clean).Get(SesaFormat.KindChild);
        int cells = (int)child.Count;
        var options = new EsaOpenOptions { Access = EsaAccessMode.Checked };
        int refused = 0;
        for (int cell = 0; cell < cells; cell++)
        {
            foreach (int value in new[] { 0, cell, cells - 1, cells, cells + 5, int.MaxValue, -1, (1 << 30) | 3, (2 << 30) | cell, (3 << 30) | (cells - 1) })
            {
                byte[] image = (byte[])clean.Clone();
                BitConverter.TryWriteBytes(image.AsSpan((int)child.Offset + (4 * cell)), value);
                using var guarded = new SesaReaderHostileTests.Guarded(image, flush: cell % 2 == 0);
                using SesaFile file = SesaFile.FromMemory(guarded.Data, guarded.Length, options);
                WideSearcher wide = WideSearcher.Create(SectionsOf(file));
                foreach (string pattern in new[] { "i", "ss", "issi", "ppi", "mississippi", "ssissippimi" })
                {
                    byte[] dense = Dense(file.TextView.Alphabet, pattern);
                    string expected = Outcome(() => (EsaIndex.DescendMapped(file.Child.Mapped, file.Sa, file.Lcp, file.TextReader, dense, true, out int lo, out int hi), lo, hi));
                    string actual = Outcome(() => (wide.TryMatch(dense, out long lo, out long hi), lo, hi));
                    if (value < 0 && !tagged)
                    {
                        // a packed cell above 2^31 is a negative number to the 32-bit descent (it reads "no boundary" and goes on); to the long descent it is a rank that no text of this size has: damage, refused. Stricter, never different otherwise.
                        Assert.That(actual, Is.EqualTo(expected).Or.EqualTo("InvalidDataException"), $"{codec} cell {cell} = {value}, pattern '{pattern}'");
                    }
                    else
                    {
                        Assert.That(actual, Is.EqualTo(expected), $"{codec} cell {cell} = {value}, pattern '{pattern}'");
                    }
                    if (actual == "InvalidDataException") refused++;
                }
            }
        }
        Assert.That(refused, Is.GreaterThan(0));
    }

    // ---------------------------------------------------------------- the generated views of a closed form

    public static IEnumerable<ClosedFormText> SmallForms()
    {
        foreach (long n in new long[] { 1, 2, 3, 17, 300, 1200 })
        {
            yield return new PowerOfA(n);
            yield return new PeriodicAkB(n, 1);
            yield return new PeriodicAkB(n, 3);
            yield return new PeriodicAkB(n, 11);
        }
    }

    [Test]
    public void TheFactorsOfAClosedFormText_AreFoundInTheGeneratedViews_AtBothWidths([ValueSource(nameof(SmallForms))] ClosedFormText f)
    {
        long n = f.Length;
        foreach (long start in new[] { 0L, 1, 2, n / 3, n / 2, Math.Max(0, n - 7) })
        {
            foreach (long length in new long[] { 1, 2, 3, 5, 17, 100, 600, n })
            {
                if (start >= n || start + length > n) continue;
                RankRange expected = f.RangeOfFactor(start, length);
                var pat = new FactorPattern(f, start, length);
                bool found5 = GenDescent.Width5(f, pat, out long lo5, out long hi5);
                bool found4 = GenDescent.Width4(f, pat, out long lo4, out long hi4);
                Assert.That(found5, Is.True, $"{f} factor [{start}, +{length}) over 5-byte views");
                Assert.That(new RankRange(lo5, hi5).SameAs(expected), Is.True, $"{f} factor [{start}, +{length}): {lo5}..{hi5}, the closed form says {expected}");
                Assert.That((found4, lo4, hi4), Is.EqualTo((found5, lo5, hi5)), $"{f} factor [{start}, +{length}): the 4-byte views agree");
            }
        }
    }

    [Test]
    public void AbsentPatterns_OfAClosedFormText_AreAbsent_InTheGeneratedViews([ValueSource(nameof(SmallForms))] ClosedFormText f)
    {
        // a pattern that is not a factor: the text with one symbol changed, a letter that the text lacks, a run longer than the text
        long n = f.Length;
        for (long start = 0; start < Math.Min(n, 5); start++)
        {
            long length = Math.Min(n - start, 40);
            var chars = new char[length];
            for (long i = 0; i < length; i++) chars[i] = f.CharAt(start + i);
            chars[length - 1] = chars[length - 1] == 'a' ? 'b' : 'a';
            string mutated = new(chars);
            RankRange expected = f.RangeOf(mutated);
            byte[] dense = [.. mutated.Select(ch => (byte)(ch - 'a'))];
            bool found = GenDescent.Width5(f, new DensePattern(dense), out long lo, out long hi);
            Assert.That(found, Is.EqualTo(!expected.IsEmpty), $"{f} '{mutated}'");
            if (found) Assert.That(new RankRange(lo, hi).SameAs(expected), Is.True, $"{f} '{mutated}'");
            else Assert.That((lo, hi), Is.EqualTo((0L, -1L)), "an absent pattern has the empty interval");
        }
        Assert.That(GenDescent.Width5(f, new RunPattern(n + 1), out _, out _), Is.False, "a run longer than the text");
    }

    private readonly struct DensePattern(byte[] symbols) : IWidePattern
    {
        public long Length => symbols.Length;

        public int this[long index] => symbols[index];
    }

    [Test]
    public void TheRuns_OfPowerOfA_AreCountedAtEveryLength_ByTheLongDescent([Values(1L, 2L, 3L, 300L, 1500L)] long n)
    {
        var f = new PowerOfA(n);
        long steps = Math.Min(n + 3, 1600);
        for (long k = 1; k <= steps; k++)
        {
            bool found = GenDescent.Width5(f, new RunPattern(k), out long lo, out long hi);
            RankRange expected = f.RangeOfRun(k);
            Assert.That(found, Is.EqualTo(!expected.IsEmpty), $"a^{k} in a^{n}");
            if (found) Assert.That((lo, hi, hi - lo + 1), Is.EqualTo((expected.First, expected.Last, n - k + 1)), $"a^{k} in a^{n}");
        }
    }

    [Test]
    public void TheRuns_OfPowerOfA_AreCountedOverTheEscapedLcpKinds([Values(255, 65_535)] int escape, [Values(300L, 1500L, 70_000L)] long n)
    {
        var f = new PowerOfA(n);
        foreach (long k in new long[] { 1, 2, escape - 1, escape, escape + 1, escape + 2, n - 1, n, n + 1 })
        {
            if (k < 1) continue;
            bool found = GenDescent.EscapedPowerOfA(f, escape, new RunPattern(k), out long lo, out long hi);
            RankRange expected = f.RangeOfRun(k);
            Assert.That(found, Is.EqualTo(!expected.IsEmpty), $"a^{k} in a^{n}, escape {escape}");
            if (found) Assert.That((lo, hi), Is.EqualTo((expected.First, expected.Last)), $"a^{k} in a^{n}, escape {escape}");
        }
    }

    [Test]
    public void ATableThatIsReadOutsideItsSection_IsInvalidData_InTheGeneratedViews()
    {
        // the generated readers are always checked: a rank that a cell sends the walk to, outside 0..n, must end as InvalidDataException
        var f = new PowerOfA(50);
        var broken = new BrokenCells(f);
        Assert.That(() => RunBroken(f, broken), Throws.TypeOf<InvalidDataException>());
    }

    private static bool RunBroken(PowerOfA f, BrokenCells cells) =>
        WideSearch.Run<CheckedAccess, PosGen<W5>, BrokenCells, GenSymbols, RunPattern>(new PosGen<W5>(f, GenSection.Sa), cells, new GenSymbols(f), f.Length, new RunPattern(10), out _, out _);

    /// <summary>Cells that send the walk backwards: the first l-index of every interval is its left bound.</summary>
    private readonly struct BrokenCells(PowerOfA f) : IWideCells
    {
        public long First(long lb, long rb, bool lastChild) => lb;

        public long Next(long k, long rb) => -1;

        public long Depth(long first) => f.Lcp(first);
    }
}
