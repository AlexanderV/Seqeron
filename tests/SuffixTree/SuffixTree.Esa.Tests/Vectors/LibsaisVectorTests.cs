using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Vectors;

/// <summary>
/// Work package 1a.6: our suffix array and LCP (<see cref="SaisReference"/>, <see cref="LcpReference"/>) against the libsais v2.9.1 output
/// recorded in <c>analysis/sa_vectors.tsv</c> (hashes only; libsais and the inputs are outside the repository). Category <c>Vectors</c>,
/// <c>[Explicit]</c>: not part of the default run, run with
/// <c>SEQERON_BENCH_DIR=&lt;dir&gt; dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter Category=Vectors</c>
/// where <c>&lt;dir&gt;</c> holds <c>data/dna_1m.txt, dna_10m.txt, rep_10m.txt, prot_1m.txt</c> (regenerate with
/// <c>python3 docs/refactoring/suffix-tree-upgrade/analysis/bench/gen/gen_inputs.py &lt;dir&gt;</c>, outside the repository).
/// A missing or modified input is a FAILURE, never a skip.
/// </summary>
[TestFixture]
[Category("Vectors")]
[Explicit("Heavy: 10 Mbp reference builds; run with --filter Category=Vectors (see PROGRESS-S1.md, Log 1a.6)")]
public class LibsaisVectorTests
{
    private const string BenchDirVariable = "SEQERON_BENCH_DIR";

    public static IEnumerable<TestCaseData> Corpora()
    {
        foreach (var e in VectorManifest.Load().Entries)
            yield return new TestCaseData(e.Corpus).SetName("Vector_" + e.Corpus);
    }

    [TestCaseSource(nameof(Corpora))]
    public void SaAndLcp_EqualLibsais(string corpus)
    {
        var entry = VectorManifest.Load()[corpus];
        string text = entry.Kind == VectorManifest.CorpusKind.File ? ReadInput(corpus) : VectorManifest.GenerateText(corpus, entry.N);
        Assert.That(text.Length, Is.EqualTo(entry.N), "input length differs from the manifest");
        Assert.That(text.All(c => c <= 255) ? "u8" : "u16", Is.EqualTo(entry.Width), "code-unit width differs from the manifest");

        var sw = System.Diagnostics.Stopwatch.StartNew();
        var (sa, lcp) = EsaPipeline.Build(text);
        double seconds = sw.Elapsed.TotalSeconds;
        string saHash = EsaPipeline.HashSa(sa);
        string lcpHash = EsaPipeline.HashLcp(lcp, entry.N);
        TestContext.Out.WriteLine($"{corpus}: n = {entry.N}, SaisReference + LcpReference {seconds:F2} s, LCP kind {lcp.Kind}, {lcp.OverflowCount} escapes");

        if (saHash != entry.SaSha256)
        {
            // Which side is wrong? The verifier is independent of both constructions.
            string? verdict = SuffixArrayVerifier.Check(text, sa);
            Assert.Fail($"{corpus}: SA hash {saHash} != libsais {entry.SaSha256}; our SA {(verdict is null ? "passes the independent verifier (suspect the libsais side, the manifest or the input)" : "FAILS the verifier: " + verdict)}");
        }
        if (lcpHash != entry.LcpSha256)
        {
            int[] kasai = KasaiLcp.Build(text, sa);
            bool oursIsKasai = lcp.ToArray().AsSpan().SequenceEqual(kasai);
            Assert.Fail($"{corpus}: LCP hash {lcpHash} != libsais {entry.LcpSha256}; our LCP {(oursIsKasai ? "equals Kasai (suspect the libsais side, the manifest or the LCP[0] convention)" : "differs from Kasai")}");
        }
    }

    private static string ReadInput(string corpus)
    {
        string? dir = Environment.GetEnvironmentVariable(BenchDirVariable);
        if (string.IsNullOrWhiteSpace(dir))
            Assert.Fail($"{BenchDirVariable} is not set: it must name a directory with data/{corpus}.txt (regenerate with analysis/bench/gen/gen_inputs.py <dir>, outside the repository)");
        string relative = $"data/{corpus}.txt";
        try
        {
            if (!CorpusManifest.LoadBenchManifest().Verify(relative, dir!))
                Assert.Fail($"{Path.Combine(dir!, relative)} does not match analysis/bench/MANIFEST.sha256 (wrong or modified input)");
        }
        catch (FileNotFoundException ex)
        {
            Assert.Fail($"input missing: {ex.FileName} (set {BenchDirVariable}; regenerate with analysis/bench/gen/gen_inputs.py <dir>)");
        }
        return System.Text.Encoding.Latin1.GetString(File.ReadAllBytes(Path.Combine(dir!, relative)));
    }
}
