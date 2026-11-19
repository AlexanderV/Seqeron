using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.Facade;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// Heavy checks of the reader of work package 1b.7f (gate G1). Category <c>Exhaustive</c>, <c>[Explicit]</c>:
/// <c>TMPDIR=/dev/shm dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter FullyQualifiedName~SesaReaderExhaustiveTests</c>. Every string over sigma 1 up to length 10, sigma 2 up to 10, sigma 3 up to 6
/// and sigma 4 up to 5 is written as a container, opened at the strongest level and compared with the in-memory tree on every member; for every suffix array of them every adjacent pair is swapped in turn and the damage (sealed again)
/// is caught; the families at 10^5 and 10^6 symbols go through the end-to-end builder and the reader, at three geometries, at every level.
/// </summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter FullyQualifiedName~SesaReaderExhaustiveTests (see PROGRESS-S1.md, Log 1b.7f)")]
public class SesaReaderExhaustiveTests
{
    private static readonly (int Sigma, int Max)[] Tiers = [(1, 10), (2, 10), (3, 6), (4, 5)];

    [Test]
    [CancelAfter(3_600_000)]
    public void EveryString_ReadThroughTheFile_EqualsTheInMemoryTree_AtTheStructureLevel()
    {
        int strings = 0;
        foreach ((int sigma, int max) in Tiers)
        {
            foreach (byte[] s in TextOracle.AllTexts(sigma, max))
            {
                strings++;
                string text = new([.. s.Select(b => (char)('A' + b))]);
                using EsaSuffixTree tree = SesaTestFile.FileTree(text, psi: strings % 3 != 0, access: strings % 2 == 0 ? EsaAccessMode.Checked : EsaAccessMode.Trusted, level: EsaValidationLevel.Structure);
                using EsaSuffixTree memory = EsaSuffixTree.Build(text);
                FacadeChecks.CheckSearch(tree, text, FacadeChecks.Patterns(text, 4), $"'{text}'");
                FacadeChecks.CheckAnalysis(tree, text, $"'{text}'");
                FacadeChecks.CheckAgainstTree(tree, memory, "esa-memory", text, $"'{text}'");
            }
        }
        Assert.That(strings, Is.GreaterThan(4000));
    }

    [Test]
    [CancelAfter(3_600_000)]
    public void EveryAdjacentSwapOfTheSuffixArray_OfEveryShortString_IsCaughtByTheStructureLevel()
    {
        int swaps = 0;
        foreach ((int sigma, int max) in new[] { (1, 8), (2, 8), (3, 5), (4, 4) })
        {
            foreach (byte[] s in TextOracle.AllTexts(sigma, max))
            {
                if (s.Length < 2) continue;
                string text = new([.. s.Select(b => (char)('A' + b))]);
                LcpChildCase c = LcpChildCase.Of(text);
                byte[] clean = SesaTestFile.Image(c);
                OracleSection sa = SesaOracle.FromBytes(clean).Get(SesaFormat.KindSa);
                for (int rank = 1; rank < c.N; rank++)
                {
                    byte[] image = (byte[])clean.Clone();
                    Span<byte> a = image.AsSpan((int)sa.Offset + (4 * rank), 4), b = image.AsSpan((int)sa.Offset + (4 * rank) + 4, 4);
                    byte[] t = a.ToArray();
                    b.CopyTo(a);
                    t.CopyTo(b);
                    SesaTestFile.Reseal(image);
                    using var dir = new TempDir();
                    string path = dir.File("s.sesa");
                    File.WriteAllBytes(path, image);
                    using (SesaTestFile.Open(path, EsaValidationLevel.Sections)) swaps++;
                    Assert.That(() => SesaTestFile.Open(path, EsaValidationLevel.Structure).Dispose(), Throws.TypeOf<InvalidDataException>(), $"'{text}' swap at rank {rank}");
                }
            }
        }
        Assert.That(swaps, Is.GreaterThan(5000));
    }

    [Test]
    [CancelAfter(3_600_000)]
    public void TheFamilies_At10To5And10To6_ThroughTheBuilderAndTheReader_EqualTheInMemoryTree()
    {
        int v = 0;
        foreach (int n in new[] { 100_000, 1_000_000 })
        {
            foreach ((string name, string text) in LcpChildRun.Families(n))
            {
                if (!SesaTestFile.IsByteText(text)) continue;
                using var memory = EsaSuffixTree.Build(text);
                foreach (long scale in new long[] { 0, 3_000_000_000, 300_000_000 })
                {
                    v++;
                    using var dir = new TempDir();
                    string input = SesaRun.WriteInput(dir, text);
                    Directory.CreateDirectory(dir.Sub("out"));
                    string target = Path.Combine(dir.Sub("out"), "x.sesa");
                    SesaBuilder.Build(input, target, new SesaBuildOptions { Threads = 1 + (v % 4), ScaleN = scale, Fsync = false });
                    EsaValidationLevel level = (EsaValidationLevel)(v % 3);
                    using EsaSuffixTree tree = SesaTestFile.Open(target, level, v % 2 == 0 ? EsaAccessMode.Checked : EsaAccessMode.Trusted, sha: level == EsaValidationLevel.Structure);
                    FacadeChecks.CheckAgainstTree(tree, memory, "esa-memory", text, $"{name} {n} scale {scale} {level}");
                    Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo(memory.LongestRepeatedSubstring()), $"{name} {n}");
                    Assert.That(tree.CountDistinctSubstrings(), Is.EqualTo(memory.CountDistinctSubstrings()), $"{name} {n}");
                    var rng = new SplitMix64((ulong)v);
                    for (int i = 0; i < 200; i++)
                    {
                        string p = text.Substring(rng.NextInt(Math.Max(1, text.Length - 30)), 1 + rng.NextInt(Math.Min(30, text.Length)));
                        Assert.That(tree.CountOccurrences(p), Is.EqualTo(memory.CountOccurrences(p)), $"{name} {n}: '{p}'");
                    }
                }
            }
        }
    }
}
