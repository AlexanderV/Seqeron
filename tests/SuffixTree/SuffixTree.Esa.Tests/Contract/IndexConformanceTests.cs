using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Contract;

/// <summary>
/// Engine conformance suite of the index contract (ADR-001 section 6, work package 1b.0): every rule of the contract that an
/// <see cref="ISuffixIndex"/> must obey, checked against definitions and brute force. Parameterised over <see cref="IndexEngineFactories.All"/>.
/// </summary>
[TestFixtureSource(typeof(IndexEngineFactories), nameof(IndexEngineFactories.All))]
public class IndexConformanceTests(IndexEngineFactory factory)
{
    private static readonly string[] SmallTexts =
    [
        "", "a", "aa", "ab", "ba", "aba", "banana", "mississippi", "abracadabra", "aaaaaaaa", "abababab", "abaababaabaab",
        "\uD800a\uD800", "\0￿\0", "xyzxyz",
    ];

    private ISuffixIndex Build(string text) => factory.Create(text);

    private static IEnumerable<string> AllSubstrings(string text, int maxLength)
    {
        var seen = new HashSet<string>(StringComparer.Ordinal);
        for (int i = 0; i < text.Length; i++)
            for (int l = 1; l <= Math.Min(maxLength, text.Length - i); l++)
                if (seen.Add(text.Substring(i, l))) yield return text.Substring(i, l);
    }

    private static IEnumerable<string> AbsentPatterns(string text)
    {
        yield return "q";
        yield return text + "q";
        yield return "q" + text;
        if (text.Length > 0)
        {
            yield return text + text[0];
            yield return new string(text[^1], text.Length + 1);
            yield return string.Concat(text.AsSpan(0, text.Length - 1), "\u0001");
        }
    }

    private static bool TryDescend(ISuffixIndex index, string pattern, out IndexLocus locus)
    {
        locus = index.Root;
        foreach (char c in pattern)
        {
            if (!index.TryExtendRight(locus, c, out var next)) return false;
            locus = next;
        }
        return true;
    }

    [Test]
    public void Metadata_EngineNameAndLengths()
    {
        foreach (string text in SmallTexts)
        {
            var index = Build(text);
            Assert.That(index.Engine, Is.Not.Null.And.Not.Empty);
            Assert.That(index.TextLength, Is.EqualTo(text.Length), text);
            Assert.That(index.Text.Length, Is.EqualTo(text.Length), text);
        }
    }

    [Test]
    public void Text_SymbolAtReturnsTheCodeUnitsAndTheTerminatorBeyondTheEnd()
    {
        foreach (string text in SmallTexts)
        {
            var index = Build(text);
            for (int i = 0; i < text.Length; i++)
                Assert.That(index.Text.SymbolAt(i), Is.EqualTo((int)text[i]), $"'{text}' at {i}");
            Assert.That(index.Text.SymbolAt(text.Length), Is.EqualTo(-1), text);
            Assert.That(index.Text.SymbolAt(text.Length + 7), Is.EqualTo(-1), text);
            Assert.That(index.Text.SymbolAt(long.MaxValue), Is.EqualTo(-1), text);
            Assert.That(() => index.Text.SymbolAt(-1), Throws.TypeOf<ArgumentOutOfRangeException>(), text);
        }
    }

    [Test]
    public void Text_CopyToClipsToTheTextAndRejectsBadStarts()
    {
        const string text = "mississippi";
        var index = Build(text);
        var buffer = new char[20];
        Assert.That(index.Text.CopyTo(0, buffer), Is.EqualTo(text.Length));
        Assert.That(new string(buffer, 0, text.Length), Is.EqualTo(text));
        Assert.That(index.Text.CopyTo(4, buffer.AsSpan(0, 3)), Is.EqualTo(3));
        Assert.That(new string(buffer, 0, 3), Is.EqualTo("iss"));
        Assert.That(index.Text.CopyTo(9, buffer), Is.EqualTo(2));
        Assert.That(new string(buffer, 0, 2), Is.EqualTo("pi"));
        Assert.That(index.Text.CopyTo(text.Length, buffer), Is.EqualTo(0));
        Assert.That(index.Text.CopyTo(3, Span<char>.Empty), Is.EqualTo(0));
        Assert.That(() => index.Text.CopyTo(text.Length + 1, new char[1]), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.Text.CopyTo(-1, new char[1]), Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    [Test]
    public void Root_CountsTheWholeText()
    {
        foreach (string text in SmallTexts)
        {
            var index = Build(text);
            Assert.That(index.Root.Length, Is.EqualTo(0), text);
            Assert.That(index.Count(index.Root), Is.EqualTo(text.Length), text);
            Assert.That(index.Root, Is.EqualTo(Build(text).Root), text);
        }
    }

    [Test]
    public void ExtendRight_CountEqualsBruteForce_ForEverySubstringAndAbsentPatterns()
    {
        foreach (string text in SmallTexts)
        {
            var index = Build(text);
            foreach (string pattern in AllSubstrings(text, 8))
            {
                Assert.That(TryDescend(index, pattern, out var locus), Is.True, $"'{text}' / '{pattern}'");
                Assert.That(locus.Length, Is.EqualTo(pattern.Length), $"'{text}' / '{pattern}'");
                Assert.That(index.Count(locus), Is.EqualTo(BruteForce.Occurrences(text, pattern).Count), $"'{text}' / '{pattern}'");
            }
            foreach (string pattern in AbsentPatterns(text))
            {
                bool expected = BruteForce.Occurrences(text, pattern).Count > 0;
                Assert.That(TryDescend(index, pattern, out _), Is.EqualTo(expected), $"'{text}' / '{pattern}'");
            }
        }
    }

    [Test]
    public void ExtendRight_PatternOverloadEqualsSymbolBySymbol()
    {
        foreach (string text in SmallTexts)
        {
            var index = Build(text);
            foreach (string pattern in AllSubstrings(text, 8).Concat(AbsentPatterns(text)))
            {
                bool viaSymbols = TryDescend(index, pattern, out var expected);
                bool viaPattern = index.TryExtendRight(index.Root, pattern.AsSpan(), out var actual);
                Assert.That(viaPattern, Is.EqualTo(viaSymbols), $"'{text}' / '{pattern}'");
                Assert.That(actual, Is.EqualTo(viaSymbols ? expected : default(IndexLocus)), $"'{text}' / '{pattern}'");
            }
        }
    }

    [Test]
    public void ExtendRight_FromAnInnerLocusEqualsTheConcatenation()
    {
        const string text = "abracadabraabracadabra";
        var index = Build(text);
        foreach (string head in AllSubstrings(text, 4))
        {
            Assert.That(TryDescend(index, head, out var headLocus), Is.True);
            foreach (string tail in AllSubstrings(text, 4))
            {
                bool whole = TryDescend(index, head + tail, out var wholeLocus);
                bool piecewise = index.TryExtendRight(headLocus, tail.AsSpan(), out var pieceLocus);
                Assert.That(piecewise, Is.EqualTo(whole), $"{head}+{tail}");
                if (whole)
                {
                    Assert.That(index.Count(pieceLocus), Is.EqualTo(index.Count(wholeLocus)), $"{head}+{tail}");
                    Assert.That(pieceLocus.Length, Is.EqualTo(head.Length + tail.Length), $"{head}+{tail}");
                }
            }
        }
    }

    [Test]
    public void ExtendRight_EmptyPatternReturnsTheLocusItself()
    {
        var index = Build("banana");
        Assert.That(index.TryExtendRight(index.Root, ReadOnlySpan<char>.Empty, out var rootAgain), Is.True);
        Assert.That(rootAgain, Is.EqualTo(index.Root));
        Assert.That(TryDescend(index, "an", out var an), Is.True);
        Assert.That(index.TryExtendRight(an, ReadOnlySpan<char>.Empty, out var anAgain), Is.True);
        Assert.That(anAgain, Is.EqualTo(an));
    }

    [Test]
    public void ExtendRight_AFailedExtensionLeavesTheResultDefault()
    {
        var index = Build("banana");
        Assert.That(index.TryExtendRight(index.Root, 'z', out var result), Is.False);
        Assert.That(result, Is.EqualTo(default(IndexLocus)));
        Assert.That(index.TryExtendRight(index.Root, "bananas".AsSpan(), out result), Is.False);
        Assert.That(result, Is.EqualTo(default(IndexLocus)));
    }

    [Test]
    public void ExtendRight_SymbolsOutsideTheCodeUnitRangeFail()
    {
        // the text holds U+0000 and U+FFFF, the code units that an engine would reach by truncating -1, 65536, int.MinValue or int.MaxValue to 16 bits
        var index = Build("ab\0\uFFFFab\0\uFFFF");
        foreach (int symbol in new[] { -1, -2, 65536, 65537, int.MaxValue, int.MinValue })
        {
            Assert.That(index.TryExtendRight(index.Root, symbol, out var result), Is.False, $"symbol {symbol}");
            Assert.That(result, Is.EqualTo(default(IndexLocus)), $"symbol {symbol}");
        }
        Assert.That(index.TryExtendRight(index.Root, 0, out var zero), Is.True);
        Assert.That(index.Count(zero), Is.EqualTo(2));
        Assert.That(index.TryExtendRight(index.Root, 65535, out var ffff), Is.True);
        Assert.That(index.Count(ffff), Is.EqualTo(2));
    }

    [Test]
    public void EmptyText_HasNoExtensionAndNothingToLocate()
    {
        var index = Build("");
        Assert.That(index.Count(index.Root), Is.EqualTo(0));
        Assert.That(index.TryExtendRight(index.Root, 'a', out _), Is.False);
        Assert.That(index.Locate(index.Root, 0, new long[4]), Is.EqualTo(0));
    }

    [Test]
    public void Locate_PagesCoverEveryPositionExactlyOnce()
    {
        foreach (string text in new[] { "banana", "mississippi", "abababab", "aaaaaaaa", "abaababaabaab" })
        {
            var index = Build(text);
            var loci = new List<(string Pattern, IndexLocus Locus)> { ("", index.Root) };
            foreach (string pattern in AllSubstrings(text, 5))
            {
                Assert.That(TryDescend(index, pattern, out var locus), Is.True);
                loci.Add((pattern, locus));
            }

            foreach (var (pattern, locus) in loci)
            {
                var expected = pattern.Length == 0 ? Enumerable.Range(0, text.Length).ToList() : BruteForce.Occurrences(text, pattern);
                foreach (int pageSize in new[] { 1, 2, 3, 1000 })
                {
                    var collected = new List<long>();
                    var page = new long[pageSize];
                    long skip = 0;
                    while (true)
                    {
                        int written = index.Locate(locus, skip, page);
                        if (written == 0) break;
                        collected.AddRange(page.Take(written));
                        skip += written;
                    }
                    Assert.That(collected.Count, Is.EqualTo(expected.Count), $"'{text}' / '{pattern}' page {pageSize}");
                    Assert.That(collected.Distinct().Count(), Is.EqualTo(collected.Count), $"'{text}' / '{pattern}' page {pageSize}");
                    Assert.That(collected.OrderBy(x => x), Is.EqualTo(expected.Select(x => (long)x)), $"'{text}' / '{pattern}' page {pageSize}");
                }
            }
        }
    }

    [Test]
    public void Locate_SkipPastTheEndAndEmptyDestination()
    {
        var index = Build("banana");
        Assert.That(TryDescend(index, "an", out var an), Is.True);
        Assert.That(index.Locate(an, 2, new long[4]), Is.EqualTo(0));
        Assert.That(index.Locate(an, 1_000_000, new long[4]), Is.EqualTo(0));
        Assert.That(index.Locate(an, 0, Span<long>.Empty), Is.EqualTo(0));
        Assert.That(index.Locate(an, 1, new long[4]), Is.EqualTo(1));
        Assert.That(() => index.Locate(an, -1, new long[4]), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.Locate(an, -1, Span<long>.Empty), Throws.TypeOf<ArgumentOutOfRangeException>(), "a negative skip is rejected even when nothing could be written");
        Assert.That(() => index.Locate(an, long.MinValue, new long[4]), Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    [Test]
    public void Locate_PartialPageWritesOnlyTheRequestedPrefix()
    {
        var index = Build("aaaaaaaa");
        Assert.That(TryDescend(index, "a", out var a), Is.True);
        var page = new long[8];
        Array.Fill(page, -7L);
        Assert.That(index.Locate(a, 0, page.AsSpan(0, 3)), Is.EqualTo(3));
        Assert.That(page.Skip(3), Is.All.EqualTo(-7L));
        Assert.That(page.Take(3), Is.All.InRange(0L, 7L));
    }

    [Test]
    public void Support_ReportsEveryOperationAndTheCoreIsNeverUnsupported()
    {
        var index = Build("banana");
        foreach (IndexOperation op in Enum.GetValues<IndexOperation>())
        {
            IndexSupport support = index.Support(op);
            Assert.That(Enum.IsDefined(support), Is.True, op.ToString());
            if (support == IndexSupport.Degraded)
                Assert.That(index.SupportNote(op), Is.Not.Null.And.Not.Empty, $"{op} is Degraded and needs a note");
        }
        foreach (var core in new[] { IndexOperation.Count, IndexOperation.ExtendRight, IndexOperation.LazyLocate, IndexOperation.TextAccess })
            Assert.That(index.Support(core), Is.Not.EqualTo(IndexSupport.Unsupported), core.ToString());
        Assert.That(() => index.Support((IndexOperation)999), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => index.SupportNote((IndexOperation)999), Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    [Test]
    public void ExtendLeft_IsEitherAnsweredCorrectlyOrThrowsNotSupported()
    {
        var index = Build("banana");
        Assert.That(TryDescend(index, "an", out var an), Is.True);
        if (index.Support(IndexOperation.ExtendLeft) == IndexSupport.Unsupported)
        {
            Assert.That(() => index.TryExtendLeft(an, 'b', out _), Throws.TypeOf<NotSupportedException>());
            return;
        }
        Assert.That(index.TryExtendLeft(an, 'b', out var ban), Is.True);
        Assert.That(index.Count(ban), Is.EqualTo(1));
        Assert.That(index.TryExtendLeft(an, 'a', out var aan), Is.False);
        Assert.That(aan, Is.EqualTo(default(IndexLocus)));
    }

    [Test]
    public void SameQueryTwiceGivesEqualLoci()
    {
        var index = Build("mississippi");
        Assert.That(TryDescend(index, "issi", out var first), Is.True);
        Assert.That(TryDescend(index, "issi", out var second), Is.True);
        Assert.That(second, Is.EqualTo(first));
        Assert.That(second.GetHashCode(), Is.EqualTo(first.GetHashCode()));
        Assert.That(TryDescend(index, "ssi", out var other), Is.True);
        Assert.That(other, Is.Not.EqualTo(first));
    }

    [Test]
    public void ConcurrentReadersAgree()
    {
        string text = StringFamilies.RandomDna(400, 17);
        var index = Build(text);
        var patterns = AllSubstrings(text, 6).Take(300).ToArray();
        long[] expected = patterns.Select(p => (long)BruteForce.Occurrences(text, p).Count).ToArray();
        var errors = new System.Collections.Concurrent.ConcurrentBag<string>();
        Parallel.For(0, 4, new ParallelOptions { MaxDegreeOfParallelism = 4 }, _ =>
        {
            for (int round = 0; round < 5; round++)
                for (int k = 0; k < patterns.Length; k++)
                {
                    if (!TryDescend(index, patterns[k], out var locus) || index.Count(locus) != expected[k])
                        errors.Add(patterns[k]);
                }
        });
        Assert.That(errors, Is.Empty);
    }

    [Test]
    public void RandomAndAdversarialTexts_CountAndLocateEqualBruteForce()
    {
        var rng = new SplitMix64(1_000_001);
        for (int round = 0; round < 12; round++)
        {
            int n = 1 + rng.NextInt(120);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, i) => (i + round) % 6 == 0))
            {
                var index = Build(text);
                var patterns = AllSubstrings(text, 5).Take(40).Concat(new[] { "zz", text + "a" });
                foreach (string pattern in patterns)
                {
                    var expected = BruteForce.Occurrences(text, pattern);
                    bool found = TryDescend(index, pattern, out var locus);
                    Assert.That(found, Is.EqualTo(expected.Count > 0), $"{name} / '{pattern}'");
                    if (!found) continue;
                    Assert.That(index.Count(locus), Is.EqualTo(expected.Count), $"{name} / '{pattern}'");
                    var all = new long[expected.Count];
                    Assert.That(index.Locate(locus, 0, all), Is.EqualTo(expected.Count), $"{name} / '{pattern}'");
                    Assert.That(all.OrderBy(x => x), Is.EqualTo(expected.Select(x => (long)x)), $"{name} / '{pattern}'");
                }
            }
        }
    }
}
