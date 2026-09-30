using System.Text;

namespace SuffixTree.Persistent.Tests.Parity;

/// <summary>
/// Seeded-random property / parity tests for review 2026-09 B05 (SUFFIXTREE-CORE F14,
/// SUFFIXTREE-PERSISTENT F16/F17) and B01's CountDistinctSubstrings, checking both the
/// in-memory <see cref="global::SuffixTree.SuffixTree"/> and <see cref="PersistentSuffixTree"/>
/// against brute-force oracles:
///   • FindAllLongestCommonSubstrings: Substring length = brute-force LCS length; PositionsInText /
///     PositionsInOther = ALL occurrences of that substring, ascending and distinct (F14); identical
///     output from both trees.
///   • Text fidelity (F17): random UTF-16 including lone high/low surrogates, U+0000/U+FFFF and pairs
///     straddling the 4096-char chunk boundary round-trips verbatim through heap reload, MMF reload
///     and Export/Import; the reloaded tree is structurally equal to the in-memory tree.
///   • CountDistinctSubstrings(ByLength) = |{ distinct substrings }| by brute force (code units).
/// Fixed seeds, bounded sizes.
/// </summary>
[TestFixture]
[Category("Parity")]
public class B05SuffixTreePropertyParityTests
{
    private static string RandomString(Random rng, int length, string alphabet)
    {
        var c = new char[length];
        for (int i = 0; i < length; i++)
            c[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(c);
    }

    /// <summary>Random UTF-16 code units biased towards surrogates and boundary values.</summary>
    private static string RandomUtf16(Random rng, int length)
    {
        var sb = new StringBuilder(length + 1);
        while (sb.Length < length)
        {
            switch (rng.Next(8))
            {
                case 0: sb.Append((char)rng.Next(0xD800, 0xDC00)); break;                // lone high
                case 1: sb.Append((char)rng.Next(0xDC00, 0xE000)); break;                // lone low
                case 2: sb.Append((char)rng.Next(0xD800, 0xDC00)).Append((char)rng.Next(0xDC00, 0xE000)); break; // pair
                case 3: sb.Append(rng.Next(2) == 0 ? '\0' : '￿'); break;
                case 4: sb.Append((char)rng.Next(0xE000, 0x10000)); break;
                default: sb.Append("ACGT"[rng.Next(4)]); break;
            }
        }
        return sb.ToString();
    }

    private static List<int> AllOccurrences(string text, string pattern)
    {
        var result = new List<int>();
        if (pattern.Length == 0) return result;
        for (int i = 0; i + pattern.Length <= text.Length; i++)
            if (string.CompareOrdinal(text, i, pattern, 0, pattern.Length) == 0)
                result.Add(i);
        return result;
    }

    private static int BruteForceLcsLength(string a, string b)
    {
        int best = 0;
        var prev = new int[b.Length + 1];
        var curr = new int[b.Length + 1];
        for (int i = 1; i <= a.Length; i++)
        {
            for (int j = 1; j <= b.Length; j++)
            {
                curr[j] = a[i - 1] == b[j - 1] ? prev[j - 1] + 1 : 0;
                if (curr[j] > best) best = curr[j];
            }
            (prev, curr) = (curr, prev);
            Array.Clear(curr);
        }
        return best;
    }

    [Test]
    public void FindAllLcs_BruteForceOccurrences_AscendingDistinct_SameInBothTrees()
    {
        for (int seed = 0; seed < 120; seed++)
        {
            var rng = new Random(31_000 + seed);
            string alphabet = (seed % 3) switch { 0 => "ab", 1 => "ACGT", _ => "a𐀀b" };
            string text = RandomString(rng, rng.Next(1, 60), alphabet);
            string other = rng.Next(4) == 0
                ? RandomString(rng, rng.Next(0, 10), "xyz")
                : RandomString(rng, rng.Next(0, 40), alphabet);

            var memory = global::SuffixTree.SuffixTree.Build(text).FindAllLongestCommonSubstrings(other);
            using var persistent = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));
            var pers = persistent.FindAllLongestCommonSubstrings(other);

            string label = $"seed={seed} text={text} other={other}";
            int lcsLen = BruteForceLcsLength(text, other);
            Assert.That(memory.Substring.Length, Is.EqualTo(lcsLen), label + " LCS length");
            if (lcsLen == 0)
            {
                Assert.That(memory.PositionsInText, Is.Empty, label);
                Assert.That(memory.PositionsInOther, Is.Empty, label);
            }
            else
            {
                Assert.That(memory.PositionsInText, Is.EqualTo(AllOccurrences(text, memory.Substring)), label + " text positions");
                Assert.That(memory.PositionsInOther, Is.EqualTo(AllOccurrences(other, memory.Substring)), label + " other positions");
            }
            Assert.That(memory.PositionsInText, Is.Ordered.Ascending.And.Unique, label);
            Assert.That(memory.PositionsInOther, Is.Ordered.Ascending.And.Unique, label);

            Assert.That(pers.Substring, Is.EqualTo(memory.Substring), label + " persistent substring");
            Assert.That(pers.PositionsInText, Is.EqualTo(memory.PositionsInText), label + " persistent text positions");
            Assert.That(pers.PositionsInOther, Is.EqualTo(memory.PositionsInOther), label + " persistent other positions");
        }
    }

    [Test]
    public void RandomUtf16Text_RoundTripsVerbatim_HeapMmfAndSerializer()
    {
        var rng = new Random(37_000);
        for (int iter = 0; iter < 12; iter++)
        {
            int length = iter switch { 0 => 0, 1 => 1, 2 => 4095, 3 => 4097, _ => rng.Next(2, iter < 8 ? 200 : 5000) };
            string text = RandomUtf16(rng, length);
            if (text.Length > length && length > 0) text = text[..length]; // may cut a pair → lone high surrogate
            var reference = global::SuffixTree.SuffixTree.Build(text);
            string label = $"iter={iter} len={text.Length}";

            var storage = new HeapStorageProvider();
            new PersistentSuffixTreeBuilder(storage).Build(new StringTextSource(text));
            using (var heap = PersistentSuffixTree.Load(storage))
            {
                Assert.That(heap.Text.ToString(), Is.EqualTo(text), label + " heap text");
                AssertSameShape(reference, heap, label + " heap");
            }

            var path = Path.Combine(Path.GetTempPath(), $"b05_utf16_{Guid.NewGuid():N}.st");
            try
            {
                using (PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text), path)) { }
                using var loaded = PersistentSuffixTreeFactory.LoadPersistent(path);
                Assert.That(loaded.Text.ToString(), Is.EqualTo(text), label + " mmf text");
                AssertSameShape(reference, loaded, label + " mmf");
            }
            finally { File.Delete(path); }

            using var ms = new MemoryStream();
            SuffixTreeSerializer.Export(reference, ms);
            ms.Position = 0;
            var imported = SuffixTreeSerializer.Import(ms, new HeapStorageProvider());
            try
            {
                Assert.That(imported.Text.ToString(), Is.EqualTo(text), label + " import text");
                Assert.That(SuffixTreeSerializer.CalculateLogicalHash(imported),
                    Is.EqualTo(SuffixTreeSerializer.CalculateLogicalHash(reference)), label + " hash");
            }
            finally { (imported as IDisposable)?.Dispose(); }
        }
    }

    /// <summary>
    /// Regression (found by the heavy tier): Pass 1 of the persistent builder collected at most 256
    /// children per node into a fixed buffer and silently dropped the rest, so any text whose root
    /// (or an inner node) branches on more than 256 distinct UTF-16 code units lost suffixes
    /// (5000-symbol alphabet, n = 4095 → 262 leaves instead of 4095). A suffix tree has exactly n
    /// leaves (Gusfield §5.2) and must equal the in-memory tree, for heap, MMF and hybrid builds.
    /// </summary>
    [TestCase(257)]
    [TestCase(300)]
    [TestCase(5000)]
    public void NodeWithMoreThan256Children_KeepsAllChildren_HeapMmfHybrid(int alphabetSize)
    {
        var rng = new Random(43_000 + alphabetSize);
        string text = new(Enumerable.Range(0, 3000).Select(_ => (char)(0x4000 + rng.Next(alphabetSize))).ToArray());
        var reference = global::SuffixTree.SuffixTree.Build(text);
        Assert.That(reference.LeafCount, Is.EqualTo(text.Length));

        using (var heap = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text)))
            AssertSameShape(reference, heap, "heap");
        using (var hybrid = PersistentSuffixTreeFactory.CreateCore(new StringTextSource(text), null, 88 + 24 * 3))
            AssertSameShape(reference, hybrid, "hybrid");
        var path = Path.Combine(Path.GetTempPath(), $"b05_wide_{Guid.NewGuid():N}.st");
        try
        {
            using (var mmf = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text), path))
                AssertSameShape(reference, mmf, "mmf");
            using var loaded = PersistentSuffixTreeFactory.LoadPersistent(path);
            AssertSameShape(reference, loaded, "mmf reload");
        }
        finally { File.Delete(path); }
    }

    private static void AssertSameShape(ISuffixTree reference, ISuffixTree actual, string label)
    {
        Assert.That(actual.NodeCount, Is.EqualTo(reference.NodeCount), label + " nodes");
        Assert.That(actual.LeafCount, Is.EqualTo(reference.LeafCount), label + " leaves");
        Assert.That(actual.CountDistinctSubstrings(), Is.EqualTo(reference.CountDistinctSubstrings()), label + " distinct");
        string text = reference.Text.ToString()!;
        for (int i = 0; i < Math.Min(text.Length, 20); i++)
        {
            int start = (i * 7919) % text.Length;
            string pat = text.Substring(start, Math.Min(3, text.Length - start));
            Assert.That(actual.FindAllOccurrences(pat).OrderBy(x => x),
                Is.EqualTo(reference.FindAllOccurrences(pat).OrderBy(x => x)), label + " occ");
        }
    }

    [Test]
    public void CountDistinctSubstrings_EqualsBruteForce_BothTrees()
    {
        for (int seed = 0; seed < 150; seed++)
        {
            var rng = new Random(41_000 + seed);
            string alphabet = (seed % 4) switch { 0 => "a", 1 => "ab", 2 => "ACGT", _ => "x𐀀￿" };
            string text = RandomString(rng, rng.Next(0, 45), alphabet);

            var distinct = new HashSet<string>(StringComparer.Ordinal);
            for (int i = 0; i < text.Length; i++)
                for (int len = 1; i + len <= text.Length; len++)
                    distinct.Add(text.Substring(i, len));
            var byLength = new long[text.Length + 1];
            byLength[0] = 1;
            foreach (var s in distinct) byLength[s.Length]++;

            var memory = global::SuffixTree.SuffixTree.Build(text);
            using var persistent = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));
            string label = $"seed={seed} text={text}";

            Assert.That(memory.CountDistinctSubstrings(), Is.EqualTo(distinct.Count), label + " memory");
            Assert.That(persistent.CountDistinctSubstrings(), Is.EqualTo(distinct.Count), label + " persistent");
            if (text.Length == 0) continue; // ByLength requires maxLength ≥ 1
            Assert.That(memory.CountDistinctSubstringsByLength(text.Length), Is.EqualTo(byLength), label + " memory by length");
            Assert.That(persistent.CountDistinctSubstringsByLength(text.Length), Is.EqualTo(byLength), label + " persistent by length");
        }
    }
}
