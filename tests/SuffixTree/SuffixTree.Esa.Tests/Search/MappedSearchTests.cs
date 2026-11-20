using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Search;

/// <summary>
/// The pattern descent over a container (<see cref="MappedSearch"/>, work package 1b.9a): on a valid container it returns what the in-memory index and the readable reference
/// (<see cref="EsaIndex.DescendMapped"/>) return, interval for interval, in both access modes, for every width of the packed text (1 to 8 bits), every representation of the LCP array and every pattern length, and on a damaged one
/// it ends as the reference does (the same answer or an <see cref="InvalidDataException"/>) without a read outside the mapping.
/// </summary>
[TestFixture]
public class MappedSearchTests
{
    private static readonly EsaAccessMode[] Modes = [EsaAccessMode.Checked, EsaAccessMode.Trusted];

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

    /// <summary>The three answers for one pattern: the memory index (the oracle), the reference descent and the specialised one.</summary>
    private static void Compare(EsaIndex memory, SesaFile file, MappedSearch search, string pattern, string label)
    {
        bool expected = memory.TryMatch(pattern.AsSpan(), out int elo, out int ehi);
        byte[] dense = Dense(file.TextView.Alphabet, pattern);
        if (dense.Length != pattern.Length)
        {
            Assert.That(expected, Is.False, $"{label}: a pattern with a letter outside the alphabet cannot occur");
            return;
        }
        bool reference = EsaIndex.DescendMapped(file.Child.Mapped, file.Sa, file.Lcp, file.TextReader, dense, file.IsChecked, out int rlo, out int rhi);
        bool found = search.TryMatch(dense, out int lo, out int hi);
        Assert.That(reference, Is.EqualTo(expected), $"{label}: the reference disagrees with the memory index for '{Show(pattern)}'");
        Assert.That(found, Is.EqualTo(expected), $"{label}: found, pattern '{Show(pattern)}'");
        if (expected)
        {
            Assert.That((lo, hi), Is.EqualTo((elo, ehi)), $"{label}: interval, pattern '{Show(pattern)}'");
            Assert.That((rlo, rhi), Is.EqualTo((elo, ehi)), $"{label}: reference interval, pattern '{Show(pattern)}'");
        }
    }

    private static string Show(string s) => s.Length <= 60 ? s : s[..60] + "...";

    private static void CompareAll(string text, string label, IEnumerable<string> patterns)
    {
        EsaIndex memory = EsaIndex.Build(text);
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(text));
        foreach (EsaAccessMode mode in Modes)
        {
            string path = Path.Combine(Path.GetTempPath(), "esa-ms-" + Guid.NewGuid().ToString("N") + ".sesa");
            File.WriteAllBytes(path, image);
            try
            {
                using SesaFile file = SesaFile.Open(path, new EsaOpenOptions { Access = mode });
                var search = new MappedSearch(file);
                foreach (string p in patterns) Compare(memory, file, search, p, $"{label} {mode}");
            }
            finally
            {
                File.Delete(path);
            }
        }
    }

    private static IEnumerable<string> AllSubstrings(string text, int maxLength)
    {
        for (int i = 0; i < text.Length; i++)
        {
            for (int len = 1; len <= Math.Min(maxLength, text.Length - i); len++) yield return text.Substring(i, len);
        }
    }

    [TestCase("banana")]
    [TestCase("mississippi")]
    [TestCase("abracadabra")]
    [TestCase("aaaaaaaaaaaa")]
    [TestCase("abababababababab")]
    [TestCase("a")]
    [TestCase("ab")]
    [TestCase("abcabcabcabcabcabcx")]
    public void EverySubstring_AndItsNeighbours_OfANamedText_EqualsTheMemoryIndexAndTheReference(string text)
    {
        var patterns = new List<string>(AllSubstrings(text, 16));
        patterns.AddRange(patterns.Take(200).Select(p => p + text[0]));      // a pattern that usually does not occur
        patterns.AddRange(patterns.Take(200).Select(p => text[^1] + p));
        patterns.Add(text + text);                                          // longer than the text
        patterns.Add(text + "Z");                                           // a letter outside the alphabet
        CompareAll(text, text, patterns);
    }

    [TestCase(1)]
    [TestCase(2)]
    [TestCase(3)]
    [TestCase(4)]
    [TestCase(5)]
    [TestCase(8)]
    [TestCase(9)]
    [TestCase(16)]
    [TestCase(17)]
    [TestCase(32)]
    [TestCase(33)]
    [TestCase(64)]
    [TestCase(65)]
    [TestCase(128)]
    [TestCase(129)]
    [TestCase(200)]
    public void RandomTexts_OfEveryWidthOfThePackedText_EqualTheMemoryIndexAndTheReference(int sigma)
    {
        var rng = new SplitMix64(0x5EED0000UL + (ulong)sigma);
        foreach (int n in new[] { 7, 64, 1500 })
        {
            // a code unit of 200 letters starts at 'A' = 65 and ends at 264: keep the letters below 256 by drawing them from a shifted range for the widest alphabets
            string text = sigma <= 190 ? RandomText(rng, n, sigma) : new string([.. RandomText(rng, n, sigma).Select(c => (char)(c - 20))]);
            var patterns = new List<string>();
            for (int i = 0; i < 300; i++)
            {
                int len = 1 + rng.NextInt(Math.Min(70, n + 3));
                patterns.Add(RandomPattern(rng, text, sigma, len, fromText: i % 3 != 0));
            }
            // every length around the stretch of one 64-bit window (64 / b symbols) and its multiples
            int bits = Math.Max(1, (int)Math.Ceiling(Math.Log2(Math.Max(2, new HashSet<char>(text).Count))));
            for (int per = 64 / bits, len = Math.Max(1, per - 1); len <= Math.Min(3 * per + 1, n); len++) patterns.Add(RandomPattern(rng, text, sigma, len, fromText: true));
            CompareAll(text, $"sigma {sigma} n {n}", patterns.Where(p => p.All(c => c < 256)));
        }
    }

    [Test]
    public void EveryPatternLength_UpToTwoWindowsOfTheWidestAndTheNarrowestText_EqualsTheMemoryIndex()
    {
        var rng = new SplitMix64(0xA11CE);
        foreach (int sigma in new[] { 2, 4, 5, 20 })
        {
            string text = RandomText(rng, 20_000, sigma);
            var patterns = new List<string>();
            for (int len = 1; len <= 140; len++)
            {
                patterns.Add(RandomPattern(rng, text, sigma, len, fromText: true));
                patterns.Add(RandomPattern(rng, text, sigma, len, fromText: true));
                string mutated = RandomPattern(rng, text, sigma, len, fromText: true);
                int at = rng.NextInt(len);
                patterns.Add(mutated[..at] + (char)('A' + ((mutated[at] - 'A' + 1) % sigma)) + mutated[(at + 1)..]);
            }
            CompareAll(text, $"sigma {sigma} lengths", patterns);
        }
    }

    [Test]
    public void APatternOfThousandsOfSymbols_EqualsTheMemoryIndex()
    {
        var rng = new SplitMix64(0xB16);
        string text = RandomText(rng, 6000, 4);
        var patterns = new List<string> { text, text[..5000], text[100..5900], text[10..] + "A", text + text };
        CompareAll(text, "long patterns", patterns);
    }

    [Test]
    public void EveryRepresentationOfTheLcpArray_EqualsTheMemoryIndex()
    {
        var rng = new SplitMix64(0x1CC);
        string repeat = RandomText(rng, 1000, 4);
        var texts = new Dictionary<string, string>
        {
            ["random (8-bit LCP)"] = RandomText(rng, 20_000, 4),
            ["a block twice (8-bit LCP and an overflow list)"] = RandomText(rng, 9000, 4) + repeat + RandomText(rng, 9000, 4) + repeat,
            ["a^2000 (16-bit LCP)"] = new string('A', 2000),
            ["a^200000 (32-bit LCP)"] = new string('A', 200_000),
        };
        var kinds = new HashSet<LcpKind>();
        bool overflow = false;
        foreach ((string name, string text) in texts)
        {
            using EsaSuffixTree probe = SesaTestFile.FileTree(text);
            kinds.Add(probe.File!.Lcp.Kind);
            overflow |= probe.File!.Lcp.Kind != LcpKind.U32 && probe.File!.Lcp.OverflowCount > 0;
            var patterns = new List<string>();
            for (int i = 0; i < 300; i++) patterns.Add(RandomPattern(rng, text, 4, 1 + rng.NextInt(1200), fromText: i % 4 != 0));
            patterns.AddRange(AllSubstrings(text.Length > 400 ? text[..40] : text, 40));
            patterns.Add(text[..Math.Min(text.Length, 1999)]);
            CompareAll(text, name, patterns);
        }
        Assert.That(kinds, Is.EquivalentTo(new[] { LcpKind.U8, LcpKind.U16, LcpKind.U32 }), "the texts must use every representation of the LCP array");
        Assert.That(overflow, Is.True, "one text must use an overflow list");
    }

    // ---------------------------------------------------------------- damaged files

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

    private static string Outcome(Func<(bool Found, int Lo, int Hi)> run)
    {
        try
        {
            (bool found, int lo, int hi) = run();
            return found ? $"found [{lo}..{hi}]" : "absent";
        }
        catch (InvalidDataException)
        {
            return "InvalidDataException";
        }
    }

    [Test]
    public void ADamagedContainer_InCheckedMode_EndsAsTheReferenceDoes_AndNeverReadsOutsideTheMapping()
    {
        Assume.That(OperatingSystem.IsLinux(), "the guard pages need mmap and mprotect");
        var rng = new SplitMix64(0xDA6A6E);
        string text = RandomText(rng, 700, 4) + new string('A', 300) + RandomText(rng, 400, 4);
        byte[] clean = SesaTestFile.Image(LcpChildCase.Of(text));
        int[] kinds = [SesaFormat.KindSa, SesaFormat.KindChild, SesaFormat.KindLcp];
        var options = new EsaOpenOptions { Access = EsaAccessMode.Checked };
        int refused = 0, answered = 0;
        for (int round = 0; round < 400; round++)
        {
            byte[] image = (byte[])clean.Clone();
            int damaged = 1 + rng.NextInt(2);
            for (int d = 0; d < damaged; d++) Damage(image, rng, kinds[rng.NextInt(kinds.Length)]);
            using var guarded = new SesaReaderHostileTests.Guarded(image, flush: round % 2 == 0);
            unsafe
            {
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
                    var search = new MappedSearch(file);
                    char[] alphabet = file.TextView.Alphabet;
                    for (int i = 0; i < 40; i++)
                    {
                        string pattern = RandomPattern(rng, text, 4, 1 + rng.NextInt(60), fromText: i % 3 != 0);
                        byte[] dense = Dense(alphabet, pattern);
                        string expected = Outcome(() => (EsaIndex.DescendMapped(file.Child.Mapped, file.Sa, file.Lcp, file.TextReader, dense, true, out int lo, out int hi), lo, hi));
                        string actual = Outcome(() => (search.TryMatch(dense, out int lo, out int hi), lo, hi));
                        Assert.That(actual, Is.EqualTo(expected), $"round {round}, pattern '{pattern}'");
                        if (actual == "InvalidDataException") refused++;
                        else answered++;
                    }
                }
            }
        }
        Assert.That(refused, Is.GreaterThan(0), "the damage must reach the descent");
        Assert.That(answered, Is.GreaterThan(0));
    }

    [Test]
    public void ATableThatSendsTheWalkBackwardsOrOutsideItsInterval_InCheckedMode_IsInvalidData()
    {
        string text = "mississippimississippi";
        byte[] clean = SesaTestFile.Image(LcpChildCase.Of(text));
        SesaOracle oracle = SesaOracle.FromBytes(clean);
        OracleSection child = oracle.Get(SesaFormat.KindChild);
        int cells = (int)child.Count;
        var options = new EsaOpenOptions { Access = EsaAccessMode.Checked };
        int refused = 0;
        // every cell of the child table set in turn to a value outside its interval, to itself and to the end of the text: no outcome but a plain answer or InvalidDataException, and the first kind at least once
        for (int cell = 0; cell < cells; cell++)
        {
            foreach (int value in new[] { 0, cell, cells - 1, cells, cells + 5, int.MaxValue, -1 })
            {
                byte[] image = (byte[])clean.Clone();
                BitConverter.TryWriteBytes(image.AsSpan((int)child.Offset + (4 * cell)), value);
                using var guarded = new SesaReaderHostileTests.Guarded(image, flush: cell % 2 == 0);
                unsafe
                {
                    using SesaFile file = SesaFile.FromMemory(guarded.Data, guarded.Length, options);
                    var search = new MappedSearch(file);
                    foreach (string pattern in new[] { "i", "ss", "issi", "ppi", "mississippi", "ssissippimi" })
                    {
                        byte[] dense = Dense(file.TextView.Alphabet, pattern);
                        string expected = Outcome(() => (EsaIndex.DescendMapped(file.Child.Mapped, file.Sa, file.Lcp, file.TextReader, dense, true, out int lo, out int hi), lo, hi));
                        string actual = Outcome(() => (search.TryMatch(dense, out int lo, out int hi), lo, hi));
                        Assert.That(actual, Is.EqualTo(expected), $"cell {cell} = {value}, pattern '{pattern}'");
                        if (actual == "InvalidDataException") refused++;
                    }
                }
            }
        }
        Assert.That(refused, Is.GreaterThan(0));
    }

    [Test]
    public void ALongerPatternThanTheSuffixOfAChild_IsAnAbsentPattern_NotAReadPastTheText()
    {
        Assume.That(OperatingSystem.IsLinux(), "the guard pages need mmap and mprotect");
        // a text whose last suffixes are shorter than the patterns that start like them: the stretch to compare would run past the end of the text
        string text = "abababababababababab";
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(text));
        foreach (bool flush in new[] { true, false })
        {
            foreach (EsaAccessMode mode in Modes)
            {
                using var guarded = new SesaReaderHostileTests.Guarded(image, flush);
                unsafe
                {
                    using SesaFile file = SesaFile.FromMemory(guarded.Data, guarded.Length, new EsaOpenOptions { Access = mode });
                    var search = new MappedSearch(file);
                    for (int len = 1; len <= 40; len++)
                    {
                        string pattern = string.Concat(Enumerable.Repeat("ab", len / 2 + 1))[..len] + "b";
                        byte[] dense = Dense(file.TextView.Alphabet, pattern);
                        bool found = search.TryMatch(dense, out _, out _);
                        Assert.That(found, Is.EqualTo(text.Contains(pattern, StringComparison.Ordinal)), $"{mode} flush {flush} '{pattern}'");
                    }
                }
            }
        }
    }
}
