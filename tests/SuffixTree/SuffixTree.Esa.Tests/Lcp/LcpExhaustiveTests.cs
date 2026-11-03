using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Lcp;

/// <summary>
/// Heavy LCP checks. Category <c>Exhaustive</c>, <c>[Explicit]</c>:
/// <c>dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter Category=Exhaustive</c>.
/// </summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter Category=Exhaustive (see PROGRESS-S1.md, Log)")]
public class LcpExhaustiveTests
{
    private static CompactLcp Build(string text, int[] sa)
    {
        var p = TextPreparation.Prepare(text);
        return p.Width == SymbolWidth.Byte ? LcpReference.Build(p.ByteReader(), sa) : LcpReference.Build(p.UShortReader(), sa);
    }

    [Test]
    public void AllBinaryStrings_UpTo20_PipelineLcpEqualsKasai()
    {
        long count = 0;
        foreach (string s in ExhaustiveStrings.All(2, 20))
        {
            int[] sa = SaisReference.BuildSuffixArray(s);
            if (!Build(s, sa).ToArray().AsSpan().SequenceEqual(KasaiLcp.Build(s, sa)))
                Assert.Fail($"PLCP LCP differs from Kasai on '{s}'");
            count++;
        }
        Assert.That(count, Is.EqualTo((1L << 21) - 1));
    }

    [Test]
    public void TextsOf1e6_PipelineLcpEqualsKasai_AndPicksTheExpectedKind()
    {
        const int n = 1_000_000;
        // Fibonacci and Thue-Morse words are so repetitive that most LCP values exceed 65,535: the rule picks U32 for them too
        foreach (var (name, text, kind) in new[]
        {
            ("a^n", StringFamilies.AllSame(n), LcpKind.U32),
            ("fibonacci", StringFamilies.Fibonacci(n), LcpKind.U32),
            ("random-dna", StringFamilies.RandomDna(n, 31), LcpKind.U8),
            ("abab", StringFamilies.AbRepeat(n), LcpKind.U32),
            ("thue-morse", StringFamilies.ThueMorse(n), LcpKind.U32),
        })
        {
            int[] sa = SaisReference.BuildSuffixArray(text);
            CompactLcp lcp = Build(text, sa);
            Assert.That(lcp.ToArray().AsSpan().SequenceEqual(KasaiLcp.Build(text, sa)), Is.True, name);
            TestContext.Out.WriteLine($"{name}: kind {lcp.Kind}, escapes {lcp.OverflowCount}, payload {lcp.PayloadBytes / (double)n:F2} B/char");
            Assert.That(lcp.Kind, Is.EqualTo(kind), name);
        }
    }
}
