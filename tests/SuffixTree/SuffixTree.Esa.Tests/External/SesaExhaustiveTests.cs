using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// Heavy checks of work package 1b.7e (gate G1). Category <c>Exhaustive</c>, <c>[Explicit]</c>:
/// <c>TMPDIR=/dev/shm dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter FullyQualifiedName~SesaExhaustiveTests</c>. Every string over sigma 1 up to length 10, sigma 2 up to 10,
/// sigma 3 up to 6 and sigma 4 up to 5 goes through the whole container writer with rotating tiny parameters (every stage), with and without PSI, and every section is compared with the in-RAM index; the families go
/// through it at 10^5 and 10^6 symbols at the default geometry and two scaled ones; the four widths of the LCP (U8 with and without an overflow list, U16 with one, U32) are written as their sections.
/// </summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter FullyQualifiedName~SesaExhaustiveTests (see PROGRESS-S1.md, Log 1b.7e)")]
public class SesaExhaustiveTests
{
    private static readonly (int Sigma, int Max)[] Tiers = [(1, 10), (2, 10), (3, 6), (4, 5)];

    [Test]
    [CancelAfter(3_600_000)]
    public void EveryString_ThroughTheContainerWriter_EqualsTheInMemoryIndex()
    {
        int v = 0;
        foreach ((int sigma, int max) in Tiers)
        {
            foreach (byte[] s in TextOracle.AllTexts(sigma, max))
            {
                v++;
                string text = new([.. s.Select(b => (char)('A' + b))]);
                bool psi = v % 3 != 0;
                using var dir = new TempDir();
                SesaRun.Outcome o = SesaRun.Build(dir, text, SesaRun.Tiny(v, 1 + (v % 4), psi));
                SesaRun.Compare(o, $"'{text}' variant {v}", psi);
            }
        }
        Assert.That(v, Is.GreaterThan(4000));
    }

    [Test]
    [CancelAfter(3_600_000)]
    public void TheFamilies_At10To5And10To6_EqualTheInMemoryIndex_AtThreeGeometries()
    {
        int v = 0;
        foreach (int n in new[] { 100_000, 1_000_000 })
        {
            foreach ((string name, string text) in LcpChildRun.Families(n))
            {
                LcpChildCase expected = LcpChildCase.Of(text);
                foreach (long scale in new long[] { 0, 3_000_000_000, 300_000_000 })
                {
                    v++;
                    using var dir = new TempDir();
                    SesaBuildOptions options = new SesaBuildOptions { Threads = 1 + (v % 4), ScaleN = scale, Fsync = false };
                    string input = SesaRun.WriteInput(dir, text);
                    Directory.CreateDirectory(dir.Sub("out"));
                    string target = Path.Combine(dir.Sub("out"), "x.sesa");
                    SesaBuildResult result = SesaBuilder.Build(input, target, options);
                    SesaRun.Compare(new SesaRun.Outcome(result, input, target, dir.Sub("out"), SesaOracle.Read(target), expected), $"{name} {n} scale {scale}", true);
                }
            }
        }
    }

    [Test]
    [CancelAfter(1_800_000)]
    public void TheWidthsOfTheLcp_U16WithOverflow_AndU32_AreWrittenAsTheirSections()
    {
        foreach ((string label, string text, LcpKind kind) in new[]
        {
            ("U16 with an overflow list", new string('A', 70_000) + StringFamilies.RandomDna(70_000, 4), LcpKind.U16),
            ("U32", new string('A', 200_000), LcpKind.U32),
        })
        {
            using var dir = new TempDir();
            SesaRun.Outcome o = SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false });
            SesaRun.Compare(o, label, true);
            Assert.That(o.Result.LcpKind, Is.EqualTo(kind), label);
            if (kind == LcpKind.U16) Assert.That(o.Result.OverflowCount, Is.GreaterThan(0), label);
        }
    }
}
